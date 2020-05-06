

class CompositePipeline(object):

    def __init__(self, config, pipelines):
        self.config = config
        self.pipelines = pipelines
    
    def run(self, dask_client):
        for pipeline in self.pipelines:
            print(pipeline)
            pipeline.run(dask_client=dask_client)
