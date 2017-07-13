'''
Created on Jul 13, 2017

@author: lubo
'''
from mapping_pipeline import MappingPipeline
from config import Config


# def test_mapping_pipeline_simple(tests_config):
#
#     assert tests_config is not None
#     assert tests_config.mapping_fastq_filenames()
#
#     fastq_filenames = tests_config.mapping_fastq_filenames()
#     print(fastq_filenames)
#
#     pipeline = MappingPipeline(tests_config)
#     assert pipeline is not None
#
#     pipeline.run()


def test_bin_count_simple(hg, tests_config):
    varbin_filenames = tests_config.varbin_data_filenames()
    print(varbin_filenames)

    for filename in varbin_filenames[3:]:
        print(filename)
        df = hg.bin_count(filename)
        cellname = tests_config.cellname(filename)
        outfile = tests_config.varbin_work_filename(cellname)
        df.to_csv(outfile, index=False, sep='\t')
        break
