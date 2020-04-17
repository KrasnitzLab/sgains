'''
Created on Aug 2, 2017

@author: lubo
'''
from distributed import Client, LocalCluster
from contextlib import closing
from dask_jobqueue import SGECluster


class Command(object):

    def __init__(self, config):
        self.config = config

    def create_local_cluster(self):
        workers = self.config.parallel
        threads_per_worker = 1
        print("workers=", workers, " threads_per_worker=", threads_per_worker)
        cluster = LocalCluster(
            n_workers=workers, threads_per_worker=threads_per_worker,
            dashboard_address=':28787')
        return cluster

    def create_sge_cluster(self):
        workers = self.config.parallel
        queue = self.config.sge_options.queue
        queue = ",".join([q.strip() for q in queue.split(',')])
        memory = self.config.sge_options.memory
        processes = int(self.config.sge_options.processes)
        cores = int(self.config.sge_options.cores)
        resource_spec = self.config.sge_options.resource_spec
        job_extra = self.config.sge_options.job_extra
        print(
            "SGE:", "queue=", queue, "memory=", memory,
            "processes=", processes, "cores=", cores,
            "resource_spec=", resource_spec,
            "job_extra=", job_extra)
        cluster = SGECluster(
            queue=queue,
            processes=processes,
            memory=memory,
            cores=cores,
            resource_spec=resource_spec,
            name="sgains-tools",
            job_extra=job_extra,
            walltime='08:00:00',
            dashboard_address=':28787',
        )
        cluster.adapt(minimum=workers, maximum=workers)
        print("SGE cluster dashboard link:", cluster.dashboard_link)
        print(cluster)
        print(cluster.job_script())
        # print(cluster.job_file())
        print("SGE cluster dashboard link:", cluster.dashboard_link)

        return cluster

    def create_dask_cluster(self):
        if self.config.sge:
            return self.create_sge_cluster()
        else:
            return self.create_local_cluster()

    def create_dask_client(self, dask_cluster):
        client = Client(dask_cluster)
        return client

    def run_pipeline(self, pipeline):
        dask_cluster = self.create_dask_cluster()
        with closing(dask_cluster) as cluster:
            dask_client = self.create_dask_client(cluster)
            print(dask_client)
            with closing(dask_client) as client:

                pipeline.run(dask_client=client)

