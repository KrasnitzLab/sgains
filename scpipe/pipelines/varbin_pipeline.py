'''
Created on Jul 31, 2017

@author: lubo
'''
from collections import defaultdict
from hg19 import HumanGenome19
from termcolor import colored
import os
import pandas as pd
import numpy as np
import pysam
import multiprocessing
import traceback


class VarbinPipeline(object):

    def __init__(self, config):
        self.config = config
        self.hg = HumanGenome19(config)

    def find_bin_index(self, abspos, bins):
        index = np.searchsorted(
            abspos, bins, side='right')

        index = index - 1
        return index

    def find_bin_index_binsearch(self, bins, abspos):
        index_up = len(bins)
        index_down = 0
        index_mid = int((index_up - index_down) / 2.0)

        while True:
            if abspos >= int(bins[index_mid]):
                index_down = index_mid + 0
                index_mid = int((index_up - index_down) / 2.0) + index_mid
            else:
                index_up = index_mid + 0
                index_mid = int((index_up - index_down) / 2.0) + index_down

            if index_up - index_down < 2:
                break

        return index_down

    def varbin(self, filename):
        try:
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

            bins = bins_df['bin.start.abspos'].values

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
                index = self.find_bin_index_binsearch(bins, abspos)

                bin_counts[index] += 1
                prev_pos = abspos
        except Exception:
            traceback.print_exc()

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

    def run_once(self, mapping_filename):
        cellname = self.config.cellname(mapping_filename)
        outfile = self.config.varbin_filename(cellname)
        print(colored(
            "processing cell {}; reading from {}; writing to {}".format(
                cellname, mapping_filename, outfile),
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
                df = self.varbin(mapping_filename)
                df.to_csv(outfile, index=False, sep='\t')

    def run(self):
        mapping_filenames = self.config.mapping_filenames()
        print(colored(
            "processing files: {}".format(mapping_filenames),
            "green"))

        pool = multiprocessing.Pool(processes=self.config.parallel)
        pool.map(self.run_once, mapping_filenames)
