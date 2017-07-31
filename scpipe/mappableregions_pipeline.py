'''
Created on Jul 31, 2017

@author: lubo
'''
from hg19 import HumanGenome19
import asyncio
from termcolor import colored
import os


class MappableRegionsPipeline(object):

    def __init__(self, config):
        self.config = config
        assert self.config.genome.version == 'hg19'
        self.hg = HumanGenome19(self.config)

    def run(self):
        outfilename = self.config.mappable_regions_filename()
        print(colored(
            "going to generate mappable regions with length {} "
            "from genome {} into {}".format(
                self.config.mappable_regions.length,
                self.config.genome.work_dir,
                outfilename
            ),
            "green"))

        if os.path.exists(outfilename) and not self.config.force:
            print(colored(
                "output file {} already exists; use --force to overwrite",
                "red"))
            raise ValueError("output file already exists")

        if self.config.dry_run:
            return

        if not os.path.exists(self.config.mappable_regions.work_dir):
            os.makedirs(self.config.mappable_regions.work_dir)

        event_loop = asyncio.get_event_loop()
        try:
            with open(self.config.mappable_regions_filename(), "w") as outfile:
                event_loop.run_until_complete(
                    self.hg.async_generate_mappable_regions(
                        self.hg.CHROMS,
                        self.config.mappable_regions.length,
                        outfile=outfile,
                        bowtie_opts=self.config.mappable_regions.bowtie_opts
                    )
                )
        finally:
            event_loop.close()
