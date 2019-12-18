import os
from contextlib import closing

import pandas as pd

from dask import distributed
from distributed import Client, LocalCluster

from sgains.config import Config

from sgains.pipelines.varbin_10x_pipeline import Varbin10xPipeline
import pytest


def relative_to_this_test_folder(path):
    return os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        path
    )


@pytest.fixture(scope='session')
def dataset10x():

    def builder(dataset):
        dataset_dir = os.path.join(
            relative_to_this_test_folder('data'),
            dataset
        )
        bins_dir = relative_to_this_test_folder('data/R100_B10k')
        config = Config.default()
        config.mapping_10x.data_10x_dir = dataset_dir
        config.bins.bins_dir = bins_dir
        config.bins.bins_boundaries = "hg19_R100_B10k_bins_boundaries.txt"
        config.mappable_regions.chrom_sizes = \
            relative_to_this_test_folder('data/chrom_sizes.yml')
        print(config)

        return config
    return builder


@pytest.mark.parametrize('dataset_name', [
    '1000000',
    '29799993',
    '57896185',
])
def test_varbin_10x_config(dataset_name, dataset10x):
    config = dataset10x(dataset_name)

    assert os.path.exists(config.bins_boundaries_filename())

    pipeline = Varbin10xPipeline(config)

    assert pipeline is not None


def test_split_bins(dataset10x):
    pipeline = Varbin10xPipeline(dataset10x('29799993'))

    regions = pipeline.split_bins(100, bins_region=(0, 200))
    # print(regions)
    assert len(regions) == 2

    regions = pipeline.split_bins(1000)
    # print(regions)
    assert len(regions) == 24

    regions = pipeline.split_bins(500)
    # print(regions)
    assert len(regions) == 33

    regions = pipeline.split_bins(50)
    # print(regions)
    assert len(regions) == 212

    regions = pipeline.split_bins(10)
    # print(regions)
    assert len(regions) == 1011

    regions = pipeline.split_bins(1)
    # print(regions)
    assert len(regions) == 10000


@pytest.mark.parametrize('bins_step', [
    5,
    1,
])
def test_process_region_reads(dataset10x, bins_step):

    pipeline = Varbin10xPipeline(dataset10x('29799993'))
    regions = pipeline.split_bins(bins_step, bins_region=(0, 100))

    for region in regions[:2]:
        cells_reads = pipeline.process_region_reads(region)
        # print(cells_reads)
        assert len(cells_reads) == 4
        print(region, {k: len(v) for k, v in cells_reads.items()})


@pytest.fixture(scope='session')
def dask_client():
    dask_cluster = LocalCluster(
        n_workers=10, threads_per_worker=2,
        dashboard_address=':28787')
    with closing(dask_cluster) as cluster:
        dask_client = Client(cluster)
        print(dask_client)
        with closing(dask_client) as client:
            yield client


@pytest.mark.parametrize('bins_step,bins_region', [
    [1, (0, 20)],
    [10, (0, 20)]
])
def test_process_region_reads_delayed(
        dask_client, dataset10x, bins_step, bins_region):

    pipeline = Varbin10xPipeline(dataset10x('29799993'))
    # regions = pipeline.split_bins(bins_step, bins_region=bins_region)
    # assert len(regions) == 2

    delayed_reads = pipeline.process_reads(
        dask_client, bins_step=bins_step, bins_region=bins_region)

    cells_reads = pipeline.merge_reads(delayed_reads)

    assert len(cells_reads) == 4
    print(set(cells_reads.keys()))

    print({k: len(v) for k, v in cells_reads.items()})
    assert len(set(cells_reads[271])) == 3997
    assert len(set(cells_reads[18])) == 2551


def test_process_region_reads_delayed_varbin_cell_reads(
        dask_client, dataset10x):
    pipeline = Varbin10xPipeline(dataset10x('29799993'))
    delayed_02reads = pipeline.process_reads(
        dask_client, bins_step=10, bins_region=(0, 20))
    delayed_20reads = pipeline.process_reads(
        dask_client, bins_step=1, bins_region=(0, 20))

    reads02 = pipeline.merge_reads(delayed_02reads)
    reads20 = pipeline.merge_reads(delayed_20reads)

    c18_02 = reads02[18]
    c18_20 = reads20[18]

    assert len(c18_02) == 3109
    assert len(c18_20) == 3111

    print(len(set(c18_02)))
    print(len(set(c18_20)))

    assert len(set(c18_02)) == len(set(c18_20))

    df1 = pipeline.varbin_cell_reads(c18_02)
    df2 = pipeline.varbin_cell_reads(c18_20)
    df3 = pipeline.varbin_cell_reads(set(c18_20))

    print(df1.head())
    print(df2.head())
    print(df3.head())

    pd.testing.assert_frame_equal(df1, df2)
    pd.testing.assert_frame_equal(df1, df3)
