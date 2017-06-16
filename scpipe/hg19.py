'''
Created on Jun 10, 2017

@author: lubo
'''
import os
from Bio import SeqIO  # @UnresolvedImport
from Bio.SeqRecord import SeqRecord  # @UnresolvedImport
import pysam
from collections import defaultdict


class HumanGenome19(object):
    VERSION = "hg19"

    CHROMS = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
        "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
        "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
        "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
    ]

    OUT = 0
    IN = 1

    def __init__(self, config):
        assert config.genome.version == self.VERSION
        self.config = config

        assert os.path.exists(config.genome.src)
        if not os.path.exists(config.genome.dst):
            os.makedirs(config.genome.dst)

        assert os.path.exists(config.genome.dst)

    def load_chrom(self, chrom, src=None):
        if src is None:
            src = self.config.genome.src
        infile = os.path.join(
            src,
            "{}.fa".format(chrom)
        )
        assert os.path.exists(infile)
        seq_record = SeqIO.read(infile, 'fasta')
        return seq_record

    def save_chrom(self, record, chrom, dst=None):
        if dst is None:
            dst = self.config.genome.dst
        outfile = os.path.join(
            dst,
            "{}.fa".format(chrom)
        )
        SeqIO.write([record], outfile, 'fasta')

    def mask_pseudoautosomal_chrY(self):
        chr_y = self.load_chrom("chrY")

        masked = chr_y.seq.tomutable()
        masked[10000:2649520] = 'N' * (2649520 - 10000)
        masked[59034049:59363566] = 'N' * (59363566 - 59034049)

        rec = SeqRecord(masked, id=chr_y.id, description=chr_y.description)
        return rec

    def generate_reads(self, chroms, read_length, src=None):
        for chrom in chroms:
            seq_record = self.load_chrom(chrom, src=src)
            for i in range(len(seq_record) - read_length + 1):
                out_record = SeqRecord(
                    seq_record.seq[i: i + read_length],
                    id="{}.{}".format(chrom, i),
                    description="generated_read"
                )
                yield out_record

    def mappable_generator(self, source):
        infile = pysam.AlignmentFile(source, 'r')  # @UndefinedVariable

        class MappableRegion(object):
            def __init__(self, mapping):
                self.chrom = mapping.reference_name
                self.start = mapping.reference_start + 1
                self.end = self.start + 1

            def extend(self, mapping):
                self.end = mapping.reference_start + 2

            def __repr__(self):
                return "{}\t{}\t{}".format(
                    self.chrom, self.start, self.end)

        prev = None
        state = self.OUT

        for mapping in infile.fetch():
            if state == self.OUT:
                if mapping.flag == 0:
                    prev = MappableRegion(mapping)
                    state = self.IN
            else:
                if mapping.flag == 0:
                    if mapping.reference_name == prev.chrom:
                        prev.extend(mapping)
                    else:
                        yield prev
                        prev = MappableRegion(mapping)
                else:
                    yield prev
                    state = self.OUT

        if state == self.IN:
            yield prev

    def count_chrom_mappable_regions(self, filename):
        result = defaultdict(lambda: 0)
        with open(filename, 'r') as infile:
            for line in infile.readlines():
                row = [r.strip() for r in line.strip().split('\t')]
                result[row[0]] += int(row[2]) - int(row[1])
        return result
