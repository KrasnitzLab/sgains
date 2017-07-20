'''
Created on Jul 9, 2017

@author: lubo
'''
import subprocess
from config import Config, NonEmptyWorkDirectory
import os
from termcolor import colored


class Rpipeline(object):

    def __init__(self, config):
        self.config = config

    def run(self):
        bin_boundaries_filename = self.config.bins_boundaries_filename()
        varbin_filenames = self.config.varbin_work_filenames()
        results_dirname = self.config.segment_work_dirname()
        study_name = self.config.segment.study_name

        print(bin_boundaries_filename)
        print(varbin_filenames)
        print(results_dirname)

        if os.path.exists(results_dirname) and \
                len(os.listdir(results_dirname)) > 0 and \
                not self.config.force:
            print(colored(
                "results directory {} is not empty; "
                "use --force to overwrite".format(results_dirname), "red"))
            raise NonEmptyWorkDirectory(results_dirname)

        if not self.config.dry_run:
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
