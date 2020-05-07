import os
import tempfile

from contextlib import closing

import pandas as pd

from dask.distributed import Client, LocalCluster, Queue, worker_client

from sgains.configuration.parser import Config

from sgains.pipelines.varbin_10x_pipeline import Varbin10xPipeline
import pytest


def relative_to_this_test_folder(path):
    return os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        path
    )


# def dataset_builder(dataset):
#     dataset_dir = os.path.join(
#         relative_to_this_test_folder('data'),
#         dataset
#     )
#     bins_dir = relative_to_this_test_folder('data/R100_B10k')
#     config = Config.default()
#     config.mapping_10x.data_10x_dir = dataset_dir
#     config.bins.bins_dir = bins_dir
#     config.bins.bins_boundaries = "hg19_R100_B10k_bins_boundaries.txt"
#     config.mappable_regions.chrom_sizes = \
#         relative_to_this_test_folder('data/chrom_sizes.yml')
#     print(config)

#     return config


@pytest.fixture(scope='session')
def dataset10x(tests_config):

    def dataset_builder(dataset):
        dataset_dir = os.path.join(
            relative_to_this_test_folder('data'),
            dataset
        )
        bins_dir = relative_to_this_test_folder('data/R100_B10k')
        data = tests_config.to_box()

        data.data_10x.data_10x_dir = dataset_dir
        data.data_10x.data_10x_cell_summary = os.path.join(
            dataset_dir,
            f"selected_test_cells_chr1_1_{dataset}_per_"
            f"cell_summary_metrics.csv"
        )
        data.data_10x.data_10x_bam = os.path.join(
            dataset_dir,
            f"selected_test_cells_chr1_1_{dataset}_possorted_bam.bam"
        )
        data.data_10x.data_10x_bai = os.path.join(
            dataset_dir,
            f"selected_test_cells_chr1_1_{dataset}_possorted_bam.bam.bai"
        )
        data.bins.bins_dir = bins_dir
        data.bins.bins_file = "hg19_R100_B10k_bins_boundaries.txt"

        data.varbin.varbin_dir = tempfile.mkdtemp("varbin_dir", "sgains")
        # data.mappable_regions.chrom_sizes = \
        #     relative_to_this_test_folder('data/chrom_sizes.yml')

        return Config.from_dict(
            data.to_dict(), relative_to_this_test_folder("."))

    return dataset_builder


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

    for index, region in enumerate(regions[:2]):
        data = pipeline.process_region_reads(region, index)
        # print(cells_reads)
        # cells_reads = Varbin10xPipeline.decompress_reads(data)
        # print(cells_reads.head())

        # assert len(cells_reads.cell_id.unique()) == 4
        # # print(region, {k: len(v) for k, v in cells_reads.items()})


@pytest.fixture(scope='function')
def dask_client():
    dask_cluster = LocalCluster(
        n_workers=15, threads_per_worker=2,
        dashboard_address=':28787')
    with closing(dask_cluster) as cluster:
        dask_client = Client(cluster)
        print("start:", dask_client)
        with closing(dask_client) as client:
            yield client
        print("done:", dask_client)


# @pytest.mark.parametrize('bins_step,bins_region', [
#     [1, (0, 20)],
#     [10, (0, 20)]
# ])
# def test_process_region_reads_delayed(
#         dask_client, dataset10x, bins_step, bins_region):

#     pipeline = Varbin10xPipeline(dataset10x('29799993'))
#     # regions = pipeline.split_bins(bins_step, bins_region=bins_region)
#     # assert len(regions) == 2

#     delayed_reads = pipeline.process_reads(
#         dask_client, bins_step=bins_step, bins_region=bins_region)

#     cells_reads = pipeline.merge_reads(dask_client, delayed_reads)

#     assert len(cells_reads) == 4
#     print(set(cells_reads.keys()))

#     print({k: len(v) for k, v in cells_reads.items()})
#     assert len(set(cells_reads[271])) == 3997
#     assert len(set(cells_reads[18])) == 2551


# def test_process_region_reads_delayed_varbin_cell_reads(
#         dask_client, dataset10x):
#     pipeline = Varbin10xPipeline(dataset10x('29799993'))
#     delayed_02reads = pipeline.process_reads(
#         dask_client, bins_step=10, bins_region=(0, 20))
#     delayed_20reads = pipeline.process_reads(
#         dask_client, bins_step=1, bins_region=(0, 20))

#     reads02 = pipeline.merge_reads(dask_client, delayed_02reads)
#     reads20 = pipeline.merge_reads(dask_client, delayed_20reads)

#     c18_02 = reads02[18]
#     c18_20 = reads20[18]

#     assert len(c18_02) == 3109
#     assert len(c18_20) == 3111

#     print(len(set(c18_02)))
#     print(len(set(c18_20)))

#     assert len(set(c18_02)) == len(set(c18_20))

#     df1 = pipeline.varbin_cell_reads(c18_02)
#     df2 = pipeline.varbin_cell_reads(c18_20)
#     df3 = pipeline.varbin_cell_reads(set(c18_20))

#     print(df1.head())
#     print(df2.head())
#     print(df3.head())

#     pd.testing.assert_frame_equal(df1, df2)
#     pd.testing.assert_frame_equal(df1, df3)


def test_varbin10x_run(dask_client, dataset10x):
    pipeline = Varbin10xPipeline(dataset10x('29799993'))
    pipeline.run(dask_client)


@pytest.mark.skip
def test_varbin10x_all(dask_client):
    dataset_dir = \
        "/home/lubo/Work/single-cell/data-single-cell/"\
        "10xGenomics/datasets/bj_mkn45_10pct"
    bins_dir = relative_to_this_test_folder('data/R100_B10k')

    config = Config.default()
    config.parallel = 15
    config.mapping_10x.data_10x_dir = dataset_dir
    config.bins.bins_dir = bins_dir
    config.bins.bins_boundaries = "hg19_R100_B10k_bins_boundaries.txt"
    config.mappable_regions.chrom_sizes = \
        relative_to_this_test_folder('data/chrom_sizes.yml')
    config.varbin.varbin_dir = \
        "/home/lubo/Work/single-cell/data-single-cell/"\
        "10xGenomics/datasets/process_bj_mkn45_10pct/bwa_varbin_10x"
    print(config)

    pipeline = Varbin10xPipeline(config)

    assert pipeline is not None

    pipeline.run(dask_client)


# def lazy_task(task_id, how_long=2):
#     print(f"task [{task_id}] started")
#     import time
#     time.sleep(how_long)
#     print(f"task [{task_id}] done")
#     return task_id


# def task_producer(total, task_queue):
#     print("inside task producer")
#     with worker_client() as client:

#         for task_id in range(total):
#             task = client.submit(lazy_task, task_id)
#             print(f"task {task_id} submitted")
#             task_queue.put(task)


# def task_consumer(task_queue):
#     collected_tasks = []
#     with worker_client() as client:
#         while task_queue.qsize() > 0:
#             task = task_queue.get()
#             result = client.gather(task)
#             print(f"task result {result} collected")
#             collected_tasks.append(result)
#     return collected_tasks


# def test_delayed_tasks(dask_client):
#     import time
#     task_queue = Queue(maxsize=20)
#     producer = dask_client.submit(task_producer, 100, task_queue)
#     print("producer started...")
#     while task_queue.qsize() == 0:
#         print("task queue empty. sleeping...")
#         time.sleep(1)

#     consumer = dask_client.submit(task_consumer, task_queue)
#     result = dask_client.gather(consumer)
#     print(result)
#     result = dask_client.gather(producer)
#     print(result)


# def test_compress_decompress_reads():
#     reads = [
#         Varbin10xPipeline.Read("1", "1", 1),
#         Varbin10xPipeline.Read("2", "1", 1),
#         Varbin10xPipeline.Read("1", "1", 2),
#         Varbin10xPipeline.Read("2", "1", 2),
#     ]
#     data = Varbin10xPipeline.compress_reads(reads)
#     assert data is not None

#     df = Varbin10xPipeline.decompress_reads(data)
#     assert df is not None
#     assert len(df) == 4
#     print(df)
