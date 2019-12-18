import os
from contextlib import closing

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


@pytest.mark.parametrize('bins_step,bins_region', [
    [1, (0, 2)],
    [10, (0, 20)]
])
def test_process_region_reads_delayed(dataset10x, bins_step, bins_region):

    pipeline = Varbin10xPipeline(dataset10x('29799993'))
    regions = pipeline.split_bins(bins_step, bins_region=bins_region)
    assert len(regions) == 2
    # for region in regions:
    #     cells_reads = pipeline.process_region_reads(region)
    #     # print(cells_reads)
    #     assert len(cells_reads) == 4
    #     print(region, {k: len(v) for k, v in cells_reads.items()})

    dask_cluster = LocalCluster(
        n_workers=10, threads_per_worker=2,
        dashboard_address=':28787')
    with closing(dask_cluster) as cluster:
        dask_client = Client(cluster)
        print(dask_client)
        with closing(dask_client) as client:
            delayed_result = pipeline.process_reads(
                client, bins_step=bins_step, bins_region=bins_region)
            distributed.wait(delayed_result)
            print(delayed_result)
            result = [r.result() for r in delayed_result]
    # print(result)
    assert len(result) == 2
    assert all([len(r) == 4 for r in result])
