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
from termcolor import colored


class HumanGenome19(object):
    VERSION = "hg19"

    CHROMS = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
        "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
        "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
        "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
    ]
    CHROMS_ALL = [
        'chr1',
        'chr1_gl000191_random',
        'chr1_gl000192_random',
        'chr2',
        'chr3',
        # 'chr4_ctg9_hap1',
        'chr4',
        'chr4_gl000193_random',
        'chr4_gl000194_random',
        'chr5',
        # 'chr6_apd_hap1',
        # 'chr6_cox_hap2',
        # 'chr6_dbb_hap3',
        'chr6',
        # 'chr6_mann_hap4',
        # 'chr6_mcf_hap5',
        # 'chr6_qbl_hap6',
        # 'chr6_ssto_hap7',
        'chr7',
        'chr7_gl000195_random',
        'chr8',
        'chr8_gl000196_random',
        'chr8_gl000197_random',
        'chr9',
        'chr9_gl000198_random',
        'chr9_gl000199_random',
        'chr9_gl000200_random',
        'chr9_gl000201_random',
        'chr10',
        'chr11',
        'chr11_gl000202_random',
        'chr12',
        'chr13',
        'chr14',
        'chr15',
        'chr16',
        # 'chr17_ctg5_hap1',
        'chr17',
        'chr17_gl000203_random',
        'chr17_gl000204_random',
        'chr17_gl000205_random',
        'chr17_gl000206_random',
        'chr18',
        'chr18_gl000207_random',
        'chr19',
        'chr19_gl000208_random',
        'chr19_gl000209_random',
        'chr20',
        'chr21',
        'chr21_gl000210_random',
        'chr22',
        'chrM',
        'chrUn_gl000211',
        'chrUn_gl000212',
        'chrUn_gl000213',
        'chrUn_gl000214',
        'chrUn_gl000215',
        'chrUn_gl000216',
        'chrUn_gl000217',
        'chrUn_gl000218',
        'chrUn_gl000219',
        'chrUn_gl000220',
        'chrUn_gl000221',
        'chrUn_gl000222',
        'chrUn_gl000223',
        'chrUn_gl000224',
        'chrUn_gl000225',
        'chrUn_gl000226',
        'chrUn_gl000227',
        'chrUn_gl000228',
        'chrUn_gl000229',
        'chrUn_gl000230',
        'chrUn_gl000231',
        'chrUn_gl000232',
        'chrUn_gl000233',
        'chrUn_gl000234',
        'chrUn_gl000235',
        'chrUn_gl000236',
        'chrUn_gl000237',
        'chrUn_gl000238',
        'chrUn_gl000239',
        'chrUn_gl000240',
        'chrUn_gl000241',
        'chrUn_gl000242',
        'chrUn_gl000243',
        'chrUn_gl000244',
        'chrUn_gl000245',
        'chrUn_gl000246',
        'chrUn_gl000247',
        'chrUn_gl000248',
        'chrUn_gl000249',
        'chrX',
        'chrY',
    ]

    CHRY_PARs = [
        (10000, 2649520),
        (59034049, 59363566),
    ]

    def __init__(self, config):
        assert config.genome.version == self.VERSION
        self.config = config

        # assert os.path.exists(config.genome.data_dir)
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

    def load_mappable_regions(self):
        df = pd.read_csv(
            self.config.mappable_regions_filename(),
            names=['chrom', 'start_pos', 'end_pos'],
            sep='\t')
        df.sort_values(by=['chrom', 'start_pos', 'end_pos'], inplace=True)
        assert len(df) > 0

        return df

    def bins_boundaries(self):
        bins_boundaries_filename = self.config.bins_boundaries_filename()
        if os.path.exists(bins_boundaries_filename):
            df = pd.read_csv(bins_boundaries_filename, sep='\t')
            return df.sort_values(by=['bin.start.abspos'])

        print(colored("bins boundaries file not found: {}".format(
            bins_boundaries_filename), "red"))

        return None

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
                directory = os.path.dirname(filename)
                if not os.path.exists(directory):
                    os.makedirs(directory)
                result.to_yaml(filename)
                self._chrom_mappable_positions_count = result
            else:
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
