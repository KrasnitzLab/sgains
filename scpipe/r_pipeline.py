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
        cells_filenames = self.config.cells_filenames()
        print(bin_boundaries_filename)
        print(cells_filenames)

        subprocess.check_call(
            [
                'Rscript', 'scripts/pipeline.R',
                bin_boundaries_filename,
                *cells_filenames
            ],
            shell=False)


if __name__ == "__main__":
    config = Config.load("scpipe_tests.yml")
    print(config)

    pipeline = Rpipeline(config)

    pipeline.run()
