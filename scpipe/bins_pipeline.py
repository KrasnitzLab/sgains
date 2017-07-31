'''
Created on Jul 31, 2017

@author: lubo
'''
from hg19 import HumanGenome19
from termcolor import colored
import os


class BinsPipeline(object):

    def __init__(self, config):
        self.config = config
        assert self.config.genome.version == 'hg19'
        self.hg = HumanGenome19(self.config)

    def run(self):
        outfile = self.config.bins_boundaries_filename()
        print(colored(
            "going to compute bin boundaries from mappable regions: {} "
            "into bins boundaries file {}".format(
                self.config.mappable_regions_filename(),
                self.config.bins_boundaries_filename()
            ),
            "green"
        ))
        if os.path.exists(outfile) and not self.config.force:
            print(colored(
                "output file {} already exists; "
                "use --force to overwrite".format(outfile),
                "red"))
            raise ValueError("output file already exists")

        regions_df = self.hg.load_mappable_regions()

        if self.config.dry_run:
            return

        bins_df = self.hg.calc_bins_boundaries(
            self.hg.CHROMS,
            regions_df
        )
        df = self.hg.calc_bins_gc_content(self.hg.CHROMS, bins_df)

        df.to_csv(outfile, sep='\t', index=False)
