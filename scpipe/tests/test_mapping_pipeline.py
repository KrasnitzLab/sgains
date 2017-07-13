'''
Created on Jul 13, 2017

@author: lubo
'''
from mapping_pipeline import MappingPipeline


def test_mapping_pipeline_simple(tests_config):

    assert tests_config is not None
    assert tests_config.mapping_fastq_filenames()

    fastq_filenames = tests_config.mapping_fastq_filenames()
    print(fastq_filenames)

    pipeline = MappingPipeline(tests_config)
    assert pipeline is not None

    pipeline.run()
