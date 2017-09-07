'''
Created on Jul 31, 2017

@author: lubo
'''
from hg19 import HumanGenome19
from termcolor import colored
import os


class VarbinPipeline(object):

    def __init__(self, config):
        self.config = config

    def run(self):

        hg = HumanGenome19(self.config)

        varbin_filenames = self.config.varbin_data_filenames()
        print(colored(
            "processing files: {}".format(varbin_filenames),
            "green"))

        for filename in varbin_filenames:
            cellname = self.config.cellname(filename)
            outfile = self.config.varbin_work_filename(cellname)
            print(colored(
                "processing cell {}; reading from {}; writing to {}".format(
                    cellname, filename, outfile),
                "green"))

            if os.path.exists(outfile) and not self.config.force:
                print(
                    colored(
                        "output file {} exists; add --force to overwrite"
                        .format(
                            outfile
                        ),
                        "red")
                )
            else:
                if not self.config.dry_run:
                    df = hg.bin_count(filename)
                    df.to_csv(outfile, index=False, sep='\t')
