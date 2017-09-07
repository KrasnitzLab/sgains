'''
Created on Jun 10, 2017

@author: lubo
'''
from collections import defaultdict
import io
import os

from Bio import SeqIO  # @UnresolvedImport
from Bio.SeqRecord import SeqRecord  # @UnresolvedImport
from box import Box

import pandas as pd
from utils import MappableState, Mapping, MappableRegion, \
    MappableBin, BinParams, LOG
import pysam
from termcolor import colored


class HumanGenome19(object):
    VERSION = "hg19"

    CHROMS = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
        "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
        "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
        "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
    ]

    def __init__(self, config):
        assert config.genome.version == self.VERSION
        self.config = config

        assert os.path.exists(config.genome.data_dir)
        if not os.path.exists(config.genome.work_dir):
            os.makedirs(config.genome.work_dir)

        assert os.path.exists(config.genome.work_dir)
        self._chrom_sizes = None
        self._chrom_bins = None
        self._chrom_mappable_positions_count = None

    def load_chrom(self, chrom, pristine=False):
        infile = self.config.chrom_filename(chrom, pristine)
        assert os.path.exists(infile), infile
        seq_record = SeqIO.read(infile, 'fasta')
        seq_record.seq = seq_record.seq.upper()
        return seq_record

    def save_chrom(self, record, chrom):
        outfile = self.config.chrom_filename(chrom)
        SeqIO.write([record], outfile, 'fasta')

    def calc_chrom_sizes(self):
        result = Box(default_box=True)
        abspos = 0
        for chrom in self.CHROMS:
            record = self.load_chrom(chrom)
            size = len(record)
            result[chrom].size = size
            result[chrom].abspos = abspos
            abspos += size
        return result

    def chrom_sizes(self):
        if self._chrom_sizes is None:
            filename = self.config.chrom_sizes_filename()
            if not os.path.exists(filename):
                result = self.calc_chrom_sizes()
                result.to_yaml(filename)

            with open(filename, 'r') as infile:
                result = Box.from_yaml(infile)
                self._chrom_sizes = result
        return self._chrom_sizes

    def bins_boundaries(self):
        bins_boundaries_filename = self.config.bins_boundaries_filename()
        print(os.path.abspath(bins_boundaries_filename))

        if os.path.exists(bins_boundaries_filename):
            df = pd.read_csv(bins_boundaries_filename, sep='\t')
            return df.sort_values(by=['bin.start.abspos'])

        return None



    CHRY_PARs = [
        (10000, 2649520),
        (59034049, 59363566),
    ]

    def mask_chrY_pars(self):
        chr_y = self.load_chrom("chrY", pristine=True)

        masked = chr_y.seq.tomutable()
        for par in self.CHRY_PARs:
            start, end = par
            masked[start:end] = 'N' * (end - start)

        rec = SeqRecord(masked, id=chr_y.id, description=chr_y.description)
        return rec

    def calc_chrom_mappable_positions_count(self):
        filename = self.config.mappable_regions_filename()
        assert os.path.exists(filename)

        result = defaultdict(lambda: 0)
        with open(filename, 'r') as infile:
            for line in infile.readlines():
                row = [r.strip() for r in line.strip().split('\t')]
                result[row[0]] += int(row[2]) - int(row[1])
        return Box(result)

    def chrom_mappable_positions_count(self):
        if self._chrom_mappable_positions_count is None:
            filename = self.config.mappable_positions_count_filename()
            if not os.path.exists(filename):
                result = self.calc_chrom_mappable_positions_count()
                result.to_yaml(filename)
            with open(filename, 'r') as infile:
                result = Box.from_yaml(infile)
                self._chrom_mappable_positions_count = result
        return self._chrom_mappable_positions_count

    def total_mappable_positions_count(self):
        counts = self.chrom_mappable_positions_count()
        total = 0
        for chrom in self.CHROMS:
            total += counts[chrom]
        return total

    def chrom_bins(self):
        if self._chrom_bins is None:
            self._chrom_bins = self.calc_chrom_bins()
        return self._chrom_bins

    def calc_chrom_bins(self):
        bins_count = self.config.bins.bins_count
        chrom_mappable_positions_count = self.chrom_mappable_positions_count()
        total_mappable_positions_count = self.total_mappable_positions_count()

        chrom_bins = Box(default_box=True)
        bins_count_used = 0
        for chrom in self.CHROMS:
            chrom_mappable_positions = chrom_mappable_positions_count[chrom]
            cbc = float(bins_count) * \
                float(chrom_mappable_positions) / \
                float(total_mappable_positions_count)
            chrom_bins[chrom].chrom = chrom
            chrom_bins[chrom].bins_count = int(cbc)
            chrom_bins[chrom].remainder = \
                cbc - chrom_bins[chrom].bins_count

            bins_count_used += chrom_bins[chrom].bins_count

        remainder_bins_count = bins_count - bins_count_used
        chroms_remainders_order = [
            c for (_, c) in
            sorted([
                (-chrom_bins[chrom].remainder, chrom)
                for chrom in chrom_bins
            ])
        ]

        for r in range(remainder_bins_count):
            chrom = chroms_remainders_order[r]
            chrom_bins[chrom].bins_count += 1

        for chrom in chroms_remainders_order:
            chrom_mappable_positions = chrom_mappable_positions_count[chrom]
            chrom_bins_count = chrom_bins[chrom].bins_count

            bin_size = float(chrom_mappable_positions) / \
                float(chrom_bins_count)
            chrom_bins[chrom].bin_size = bin_size

        return chrom_bins

    def load_mappable_regions(self):
        df = pd.read_csv(
            self.config.mappable_regions_filename(),
            names=['chrom', 'start_pos', 'end_pos'],
            sep='\t')
        df.sort_values(by=['chrom', 'start_pos', 'end_pos'], inplace=True)

        return df

    @staticmethod
    def to_fasta_string(rec):
        out_handle = io.StringIO()
        out_handle.write(">{}\n".format(rec.id))
        out_handle.write(str(rec.seq).upper())
        out_handle.write("\n")
        return out_handle.getvalue().encode('utf-8')

    @staticmethod
    def write_fasta_read(outfile, rec):
        outfile.write(HumanGenome19.to_fasta_string(rec))
