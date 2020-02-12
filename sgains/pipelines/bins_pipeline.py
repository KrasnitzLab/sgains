'''
Created on Jul 31, 2017

@author: lubo
'''
from dask import distributed

import os
import pandas as pd
from termcolor import colored

from sgains.genome import Genome
from sgains.utils import BinParams, MappableBin


class BinsPipeline(object):

    def __init__(self, config):
        self.config = config
        self.hg = Genome(self.config)

    def calc_bins_gc_content(self, chroms, bins_df):

        result = []
        for chrom in chroms:
            chrom_df = bins_df[bins_df['bin.chrom'] == chrom]
            gc_df = chrom_df.copy()
            gc_df.reset_index(inplace=True, drop=True)

            gc_series = pd.Series(index=gc_df.index)
            chrom_seq = self.hg.load_chrom(chrom)

            for index, row in gc_df.iterrows():
                start = row['bin.start']
                end = row['bin.end']
                seq = chrom_seq.seq[start:end]
                counts = [seq.count(x) for x in ['G', 'C', 'A', 'T']]
                gc = float(sum(counts[0:2])) / sum(counts)
                gc_series.iloc[index] = gc

            gc_df['gc.content'] = gc_series
            result.append(gc_df)
        assert len(result) > 0
        if len(result) == 1:
            return result[0]
        df = pd.concat(result)
        return df

    def bins_boundaries_generator(self, chroms, mappable_regions_df):
        chrom_sizes = self.hg.chrom_sizes()
        chrom_bins = self.hg.calc_chrom_bins()

        # if mappable_regions_df is None:
        #     mappable_regions_df = self.load_mappable_regions()

        for chrom in chroms:
            chrom_df = mappable_regions_df[mappable_regions_df.chrom == chrom]
            chrom_df = chrom_df.sort_values(
                by=['chrom', 'start_pos', 'end_pos'])

            params = BinParams.build(
                chrom_size=chrom_sizes[chrom],
                chrom_bin=chrom_bins[chrom])
            mappable_bin = None
            current_excess = 0
            bins_count = params.bins_count

            for _index, row in chrom_df.iterrows():
                if mappable_bin is None:
                    mappable_bin = MappableBin.from_start(params, start_pos=0)
                    current_excess = mappable_bin.adapt_excess(current_excess)
                if not mappable_bin.check_extend(row):
                    next_bin = mappable_bin.split_extend(row)

                    bins_count -= 1
                    if bins_count == 0:
                        # last bin a chromosome
                        mappable_bin.end_pos = chrom_sizes[chrom].size
                    yield mappable_bin
                    if next_bin.is_overfill():
                        current_excess, mappable_bins = \
                            next_bin.overfill_split(current_excess)

                        assert len(mappable_bins) > 1
                        for mb in mappable_bins[:-1]:
                            bins_count -= 1
                            yield mb
                        mappable_bin = mappable_bins[-1]
                    else:
                        mappable_bin = next_bin
                        current_excess = \
                            mappable_bin.adapt_excess(current_excess)

            mappable_bin = None

    def calc_bins_boundaries(self, chroms=None, regions_df=None):
        if chroms is None:
            chroms = self.hg.version.CHROMS
        bin_rows = []
        for mbin in self.bins_boundaries_generator(chroms, regions_df):
            bin_rows.append(
                (
                    mbin.chrom,
                    mbin.start_pos,
                    mbin.start_abspos,
                    mbin.end_pos,
                    mbin.end_pos - mbin.start_pos,
                    mbin.bin_size
                )
            )

        df = pd.DataFrame.from_records(
            bin_rows,
            columns=[
                'bin.chrom',
                'bin.start',
                'bin.start.abspos',
                'bin.end',
                'bin.length',
                'mappable.positions'
            ])
        df.sort_values(by=['bin.start.abspos'], inplace=True)
        return df

    def run_once(self, chrom):
        print(colored(
            f"started calculating bins for chromosome {chrom}",
            "green"
        ))
        regions_df = self.hg.load_mappable_regions()
        bins_df = self.calc_bins_boundaries(
            [chrom],
            regions_df
        )
        df = self.calc_bins_gc_content([chrom], bins_df)
        outfilename = self.config.bins_boundaries_filename(chrom)
        print(colored(
            f"saving bins for chromosome {chrom} into {outfilename}",
            "green"
        ))
        df.to_csv(outfilename, sep='\t', index=False)
        return outfilename

    def concatenate_all_chroms(self):
        outfilename = self.config.bins_boundaries_filename()
        if os.path.exists(outfilename) and not self.config.force:
            print(colored(
                "destination bins boundaries file already exists"
                "use --force to overwrite", "red"))
            raise ValueError("destination file exists... use --force")

        if self.config.dry_run:
            return

        dataframes = []
        for chrom in self.hg.version.CHROMS:
            srcfilename = self.config.bins_boundaries_filename(chrom)
            df = pd.read_csv(srcfilename, sep='\t')
            dataframes.append(df)
        outdf = pd.concat(dataframes, ignore_index=True)
        outdf.sort_values(
            by=['bin.start.abspos', 'bin.start', 'bin.end'], inplace=True)

        outdf.to_csv(outfilename, sep='\t', index=False)

    def run(self, dask_client):
        outfilename = self.config.bins_boundaries_filename()
        os.makedirs(os.path.dirname(outfilename), exist_ok=True)

        print(colored(
            "going to compute bin boundaries from mappable regions: {} "
            "into bins boundaries file {}".format(
                self.config.mappable_regions_filename(),
                outfilename
            ),
            "green"
        ))
        if os.path.exists(outfilename) and not self.config.force:
            print(colored(
                "output file {} already exists; "
                "use --force to overwrite".format(outfilename),
                "red"))
            raise ValueError("output file already exists")

        if self.config.dry_run:
            return

        # pool = multiprocessing.Pool(processes=self.config.parallel)
        # pool.map(self.run_once, self.hg.version.CHROMS)
        delayed_tasks = dask_client.map(
                self.run_once, self.hg.version.CHROMS)
        print(len(delayed_tasks), delayed_tasks)
        print(dask_client.scheduler_info())

        distributed.wait(delayed_tasks)
        for task in delayed_tasks:
            outfile = task.result()
            print(outfile, os.path.exists(outfile))

        self.concatenate_all_chroms()
