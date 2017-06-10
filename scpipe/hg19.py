'''
Created on Jun 10, 2017

@author: lubo
'''
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


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

        print(chr_y.id)
        print(len(chr_y))

        masked = chr_y.seq.tomutable()

        masked[10000:2649520] = 'N' * (2649520 - 10000)
        masked[59034049:59363566] = 'N' * (59363566 - 59034049)

        rec = SeqRecord(masked, id=chr_y.id, description=chr_y.description)
        return rec
