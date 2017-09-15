'''
Created on Jul 31, 2017

@author: lubo
'''
from collections import defaultdict
from hg19 import HumanGenome19
from termcolor import colored
import os
import pandas as pd
import pysam
import multiprocessing


class VarbinPipeline(object):

    def __init__(self, config):
        self.config = config
        self.hg = HumanGenome19(config)

    def varbin(self, filename):

        assert os.path.exists(filename), os.path.abspath(filename)

        infile = pysam.AlignmentFile(filename, 'rb')  # @UndefinedVariable
        bins_df = self.hg.bins_boundaries()
        assert bins_df is not None

        chrom_sizes = self.hg.chrom_sizes()
        chroms = set(self.hg.CHROMS)

        count = 0
        dups = 0
        total_reads = 0

        prev_pos = 0
        bin_counts = defaultdict(int)

        for seg in infile:
            total_reads += 1
            if seg.is_unmapped:
                continue
            chrom = seg.reference_name
            if chrom not in chroms:
                continue

            abspos = chrom_sizes[chrom].abspos + seg.reference_start
            if prev_pos == abspos:
                dups += 1
                continue
            count += 1
            index = bins_df['bin.start.abspos'].searchsorted(
                abspos, side='right')
            assert len(index) == 1

            index = index[0] - 1
            bin_counts[index] += 1
            prev_pos = abspos
        print(total_reads, dups, count)

        number_of_reads_per_bin = float(count) / len(bins_df)
        result = []
        for index, row in bins_df.iterrows():
            bin_count = bin_counts[index]
            ratio = float(bin_count) / number_of_reads_per_bin
            result.append(
                [
                    row['bin.chrom'],
                    row['bin.start'],
                    row['bin.start.abspos'],
                    bin_count,
                    ratio
                ]
            )
        df = pd.DataFrame.from_records(
            result,
            columns=[
                'chrom',
                'chrompos',
                'abspos',
                'bincount',
                'ratio',
            ])
        df.sort_values(by=['abspos'], inplace=True)
        return df

    def run_once(self, varbin_filename):
        cellname = self.config.cellname(varbin_filename)
        outfile = self.config.varbin_work_filename(cellname)
        print(colored(
            "processing cell {}; reading from {}; writing to {}".format(
                cellname, varbin_filename, outfile),
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
                df = self.varbin(varbin_filename)
                df.to_csv(outfile, index=False, sep='\t')

    def run(self):
        varbin_filenames = self.config.varbin_data_filenames()
        print(colored(
            "processing files: {}".format(varbin_filenames),
            "green"))

        pool = multiprocessing.Pool(processes=self.config.parallel)
        pool.map(self.run_once, varbin_filenames)
