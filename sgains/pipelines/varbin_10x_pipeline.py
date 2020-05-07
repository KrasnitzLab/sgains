import os
import time
import glob
import shutil

from io import BytesIO
from collections import defaultdict, namedtuple

import pandas as pd
import numpy as np

from dask.distributed import Queue, worker_client, wait

from termcolor import colored

import pysam

from sgains.genome import Genome
from sgains.pipelines.extract_10x_pipeline import Base10xPipeline


class Varbin10xPipeline(Base10xPipeline):

    def __init__(self, config):
        super(Varbin10xPipeline, self).__init__(config)

        self.bins_df = self.genome.bins_boundaries()
        self.chrom2contig, self.contig2chrom = self._chrom2contig_mapping()
        self.chrom_sizes = self.genome.chrom_sizes()

    def _chrom2contig_mapping(self):
        with pysam.AlignmentFile(self.bam_filename, 'rb') as samfile:
            assert samfile.check_index(), \
                (self.bam_filename, self.bai_filename)
            sam_stats = samfile.get_index_statistics()
        chroms = set(self.bins_df['bin.chrom'].values)
        chrom2contig = {}
        contig2chrom = {}
        for stat in sam_stats:
            contig = stat.contig
            chrom_name = contig
            if chrom_name not in chroms:
                chrom_name = "chr{}".format(contig)
            if chrom_name not in chroms:
                continue
            chrom2contig[chrom_name] = contig
            contig2chrom[contig] = chrom_name
        return chrom2contig, contig2chrom

    Region = namedtuple('Region', ['chrom', 'contig', 'start', 'end'])

    def split_bins(self, bins_step, bins_region=None):
        total_bins = len(self.bins_df)
        regions = []
        bin_start = 0
        bin_end = total_bins
        if bins_region is not None:
            bin_start, bin_end = bins_region
            bin_end = min(bin_end, total_bins)

        index = bin_start
        while index < bin_end:
            start = index
            end = index + bins_step - 1
            if end >= total_bins:
                end = total_bins - 1

            chrom_start = self.bins_df.iloc[start, 0]
            chrom_end = self.bins_df.iloc[end, 0]
            while chrom_end != chrom_start:
                end -= 1
                assert end >= start
                chrom_end = self.bins_df.iloc[end, 0]

            pos_start = self.bins_df.iloc[start, 1]
            pos_end = self.bins_df.iloc[end, 3]

            regions.append((chrom_start, pos_start, pos_end))
            index = end + 1

        return [
            self.Region(chrom, self.chrom2contig[chrom], start, end)
            for chrom, start, end in regions
        ]

    # @staticmethod
    # def compress_reads(reads):
    #     if not reads:
    #         return None
    #     df = pd.DataFrame(reads, columns=['cell_id', 'chrom', 'pos'])
    #     buff = BytesIO()
    #     df.to_feather(buff)

    #     return str(buff.getbuffer(), encoding='latin1')

    @staticmethod
    def _cell_name(cell_id):
        cell_name = f"C{cell_id:0>8}"
        return cell_name

    def _cell_reads_dirname(self, cell_id):
        cell_name = self._cell_name(cell_id)
        dirname = os.path.join(
            self.config.varbin.varbin_dir,
            cell_name)
        os.makedirs(dirname, exist_ok=True)
        return dirname

    def _cell_region_filename(self, cell_id, region_index):
        cell_name = self._cell_name(cell_id)
        region_name = f"{region_index:0>8}"
        filename = os.path.join(
            self.config.varbin.varbin_dir,
            cell_name,
            f"{cell_name}_{region_name}{self.config.varbin.varbin_suffix}"
        )
        return filename

    def store_reads(self, reads, region_index):
        if not reads:
            return None
        df = pd.DataFrame(reads, columns=['cell_id', 'chrom', 'pos'])

        for cell_id, group_df in df.groupby(by="cell_id"):
            cell_region_filename = self._cell_region_filename(
                cell_id, region_index)

            cell_dirname = os.path.dirname(cell_region_filename)
            os.makedirs(cell_dirname, exist_ok=True)
            group_df.to_csv(cell_region_filename, index=False, sep="\t")
        return "done"

    def load_reads(self, cell_id):
        cell_name = self._cell_name(cell_id)

        pattern = f"{cell_name}_*{self.config.varbin.varbin_suffix}"
        pattern = os.path.join(
            self.config.varbin.varbin_dir,
            cell_name,
            pattern
        )
        print(colored(
            "merging reads cell files {} ".format(
                pattern
            ),
            "green"))
        filenames = glob.glob(pattern)
        filenames = sorted(filenames)

        dataframes = []
        for filename in filenames:
            df = pd.read_csv(filename, sep="\t")
            dataframes.append(df)
        if len(dataframes) == 0:
            return None
        elif len(dataframes) == 1:
            return dataframes[0]
        else:
            result_df = pd.concat(dataframes, ignore_index=True)
            result_df = result_df.sort_values(by=["cell_id", "chrom", "pos"])
            return result_df

    # @staticmethod
    # def decompress_reads(data):
    #     if data is None:
    #         return None
    #     buff = BytesIO(bytes(data, encoding='latin1'))
    #     df = pd.read_feather(buff)
    #     return df

    Read = namedtuple('Read', ['cell_id', 'chrom', 'pos'])

    def process_region_reads(self, region, region_index):
        print(f"started region {region}")
        with pysam.AlignmentFile(self.bam_filename, 'rb') as samfile:
            assert samfile.check_index(), \
                (self.bam_filename, self.bai_filename)
            cells_reads = []
            mapped = 0
            for rec in samfile.fetch(region.contig, region.start, region.end):
                if not rec.has_tag('CB'):
                    continue
                mapped += 1
                barcode = rec.get_tag('CB')
                if barcode not in self.barcodes:
                    continue
                cell_id = self.barcodes[barcode]
                contig = rec.reference_name
                if contig not in self.contig2chrom:
                    continue
                read = self.Read(
                    cell_id, self.contig2chrom[contig], rec.reference_start)
                cells_reads.append(read)
            print(f"done region {region}; reads processed {len(cells_reads)}")
            return self.store_reads(cells_reads, region_index)

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

    def varbin_cell_reads(self, reads_df):
        assert self.bins_df is not None
        count = 0
        dups = 0
        total_reads = 0

        prev_pos = 0
        bin_counts = defaultdict(int)

        bins = self.bins_df['bin.start.abspos'].values

        for _, read in reads_df.iterrows():
            total_reads += 1

            abspos = self.chrom_sizes[read.chrom].abspos + read.pos
            if prev_pos == abspos:
                dups += 1
                continue
            count += 1
            index = self.find_bin_index_binsearch(bins, abspos)

            bin_counts[index] += 1
            prev_pos = abspos

        number_of_reads_per_bin = float(count) / len(self.bins_df)
        result = []
        for index, row in self.bins_df.iterrows():
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

    def process_cell_reads(self, cell_id):
        reads_df = self.load_reads(cell_id)
        if reads_df is None:
            print(f"data for cell {cell_id} not found... skipping...")
            return

        df = self.varbin_cell_reads(reads_df)
        outfile = self.config.varbin_filename(f"C{cell_id:0>5}")
        df.to_csv(outfile, sep='\t', index=False)

        cell_dirname = self._cell_reads_dirname(cell_id)
        shutil.rmtree(cell_dirname)

        return cell_id

    def run(self, dask_client, bins_step=20, bins_region=None, outdir='.'):

        regions = self.split_bins(bins_step=bins_step, bins_region=bins_region)

        delayed_tasks = dask_client.map(
            lambda region_tuple: 
            self.process_region_reads(region_tuple[1], region_tuple[0]),
            list(enumerate(regions))
        )
        wait(delayed_tasks)

        delayed_tasks = dask_client.map(
            lambda cell_id: self.process_cell_reads(cell_id),
            list(self.barcodes.values())
        )
        wait(delayed_tasks)
