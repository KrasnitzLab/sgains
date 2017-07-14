'''
Created on Jul 9, 2017

@author: lubo
'''
import subprocess
from config import Config


class Rpipeline(object):

    def __init__(self, config):
        self.config = config

    def run(self):
        bin_boundaries_filename = self.config.bin_boundaries_filename()
        varbin_filenames = self.config.varbin_work_filenames()
        results_dirname = self.config.results_work_dirname()
        study_name = self.config.results.study_name

        print(bin_boundaries_filename)
        print(varbin_filenames)

        subprocess.check_call(
            [
                'Rscript', 'scripts/pipeline.R',
                study_name,
                results_dirname,
                bin_boundaries_filename,
                *varbin_filenames
            ],
            shell=False)


if __name__ == "__main__":
    config = Config.load("scpipe_tests.yml")
    print(config)

    pipeline = Rpipeline(config)

    pipeline.run()
