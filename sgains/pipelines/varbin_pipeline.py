'''
Created on Jul 31, 2017

@author: lubo
'''
from collections import defaultdict
from sgains.genome import Genome
from termcolor import colored
import os
import pandas as pd
import numpy as np
import pysam
import traceback
from dask import distributed


class VarbinPipeline(object):

    def __init__(self, config):
        self.config = config
        self.hg = Genome(config)

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

            infile = pysam.AlignmentFile(filename, 'rb')
            bins_df = self.hg.bins_boundaries()
            assert bins_df is not None
            chrom_sizes = self.hg.chrom_sizes()
            chroms = set(self.hg.version.CHROMS)

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
                assert seg.cigarstring == f'{seg.reference_length}M', \
                    (seg, seg.cigarstring)

                abspos = chrom_sizes[chrom].abspos + seg.reference_start
                if prev_pos == abspos:
                    dups += 1
                    continue
                count += 1
                index = self.find_bin_index_binsearch(bins, abspos)

                bin_counts[index] += 1
                prev_pos = abspos

            result = []
            for index, row in bins_df.iterrows():
                bin_count = bin_counts[index]
                result.append(
                    [
                        row['bin.chrom'],
                        row['bin.start'],
                        row['bin.start.abspos'],
                        bin_count,
                    ]
                )
            df = pd.DataFrame.from_records(
                result,
                columns=[
                    'chrom',
                    'chrompos',
                    'abspos',
                    'bincount',
                ])

            df.sort_values(by=['abspos'], inplace=True)
            total_count = df.bincount.sum()
            total_reads_per_bin = float(total_count)/len(bins_df)
            df['ratio'] = df.bincount / total_reads_per_bin

            return df
        except Exception as ex:
            traceback.print_exc()
            raise ex
        return None

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

    def run(self, dask_client):
        mapping_filenames = self.config.mapping_all_filenames()
        print(colored(
            "processing files: {}".format(mapping_filenames),
            "green"))

        if self.config.dry_run:
            return

        assert dask_client

        delayed_tasks = dask_client.map(self.run_once, mapping_filenames)
        distributed.wait(delayed_tasks)

        # for fut in delayed_tasks:
        #     print("fut done:", fut.done())
        #     print("fut exception:", fut.exception())
        #     print("fut traceback:", fut.traceback())
        #     if fut.traceback() is not None:
        #         traceback.print_tb(fut.traceback())
        #     print(fut.result())

