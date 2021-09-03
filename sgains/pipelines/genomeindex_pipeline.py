'''
Created on Jul 31, 2017

@author: lubo
'''
from sgains.genome import Genome
import os
import shutil
from termcolor import colored
import subprocess


class GenomeIndexPipeline(object):

    def __init__(self, config):
        self.config = config
        # assert self.config.genome.version == 'hg19'
        self.genome = Genome(self.config)
        assert self.genome.aligner is not None

    def copy_chromes_files(self):
        self.config.check_nonempty_workdir(self.config.genome.genome_dir)

        for chrom in self.genome.version.CHROMS_ALL:
            if chrom == 'chrY':
                continue
            src = os.path.join(
                self.config.genome.genome_pristine_dir,
                "{}.fa".format(chrom)
            )
            dst = os.path.join(
                self.config.genome.genome_dir,
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
        dst = self.genome.chrom_filename('chrY')
        print(colored(
            "masking pseudoautosomal regions in chrY",
            "green")
        )
        if os.path.exists(dst) and not self.config.force:
            print(colored(
                "destination file for masked chrY already exists",
                "red"
            ))
            raise ValueError("dst file already exists")
        if not self.config.dry_run:
            masked = self.genome.mask_chrY_pars()
            self.genome.save_chrom(masked, 'chrY')

    def concatenate_all_chroms(self):
        dirname = self.config.genome.genome_dir
        dst = os.path.join(
            dirname,
            'genome.fa'
        )
        if os.path.exists(dst) and not self.config.force:
            print(colored(
                "destination genome file already exists"
                "use --force to overwrite", "red"))
            raise ValueError("destination file exists... use --force")

        if not self.config.dry_run:
            with open(dst, 'wb') as output:
                for chrom in self.genome.version.CHROMS_ALL:
                    src = self.genome.chrom_filename(chrom, pristine=False)
                    print(colored(
                        "appending {} to {}".format(src, dst),
                        "green"))
                    with open(src, 'rb') as src:
                        if not self.config.dry_run:
                            shutil.copyfileobj(src, output, 1024 * 1024 * 10)

    def build_aligner_index(self):
        print(colored(
            f"building genome index of {self.genome.sequence_filename} "
            f"into {self.genome.index_prefix}",
            "green"))
        command = " ".join(self.genome.aligner.build_index_command(
            self.genome.sequence_filename,
            self.genome.index_prefix
        ))
        print(colored(
            f"going to execute aligner genome index build: {command}",
            "green"))

        test_filename = self.genome.aligner.genome_index_filenames[0]
        print(colored(f"checking for index file: {test_filename}", "green"))
        if os.path.exists(test_filename) and not self.config.force:
            print(colored(
                "output genome index {} already exists".format(test_filename),
                "red"))
            raise ValueError("destination file already exists")

        if not self.config.dry_run:
            subprocess.check_call(command, shell=True)

    def run(self, **kwargs):
        self.copy_chromes_files()
        self.mask_pars()
        self.concatenate_all_chroms()
        self.build_aligner_index()
