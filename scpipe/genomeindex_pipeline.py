'''
Created on Jul 31, 2017

@author: lubo
'''
from hg19 import HumanGenome19
import os
import shutil
from termcolor import colored


class GenomeIndexPipeline(object):

    def __init__(self, config):
        self.config = config
        assert self.config.genome.version == 'hg19'
        self.hg = HumanGenome19(self.config)

    def copy_chromes_files(self):
        self.config.check_nonempty_workdir(self.config.genome.work_dir)

        for chrom in self.hg.CHROMS:
            src = os.path.join(
                self.config.genome.pristine,
                "{}.fa".format(chrom)
            )
            dst = os.path.join(
                self.config.genome.work_dir,
                "{}.fa".format(chrom)
            )
            print(colored(
                "copying chromosome {} from {} into "
                "working directory {}".format(
                    chrom, src, dst),
                "green"))
            if not self.config.dry_run:
                shutil.copy(src, dst)

    def mask_pars(self):
        pass

    def concatenate_all_chroms(self):
        pass

    def build_bowtie_index(self):
        pass

    def run(self):

        assert self.hg is not None
