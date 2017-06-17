'''
Created on Jun 10, 2017

@author: lubo
'''
import os
from Bio import SeqIO  # @UnresolvedImport
from Bio.SeqRecord import SeqRecord  # @UnresolvedImport
import pysam
from collections import defaultdict
from box import Box
from subprocess import Popen, PIPE
import shlex
from multiprocessing import Process
import io


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


class MappableState(object):
    OUT = 0
    IN = 1


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

        masked = chr_y.seq.tomutable()
        masked[10000:2649520] = 'N' * (2649520 - 10000)
        masked[59034049:59363566] = 'N' * (59363566 - 59034049)

        rec = SeqRecord(masked, id=chr_y.id, description=chr_y.description)
        return rec

    def generate_reads(self, chroms, read_length, src=None):
        try:
            for chrom in chroms:
                seq_record = self.load_chrom(chrom, src=src)
                for i in range(len(seq_record) - read_length + 1):
                    out_record = SeqRecord(
                        seq_record.seq[i: i + read_length],
                        id="{}.{}".format(chrom, i),
                        description="generated_read"
                    )
                    yield out_record
        finally:
            print("closing generate reads...")

    def mappable_regions_generator(self, mappings_generator):
        try:
            prev = None
            state = MappableState.OUT

            for mapping in mappings_generator:
                if state == MappableState.OUT:
                    if mapping.flag == 0:
                        prev = MappableRegion(mapping)
                        state = MappableState.IN
                else:
                    if mapping.flag == 0:
                        if mapping.reference_name == prev.chrom:
                            prev.extend(mapping)
                        else:
                            yield prev
                            prev = MappableRegion(mapping)
                    else:
                        yield prev
                        state = MappableState.OUT

            if state == MappableState.IN:
                yield prev
        finally:
            print("closing mappable_regions_generator...")
            mappings_generator.close()

    @staticmethod
    def write_fasta_read(outfile, rec):
        out_handle = io.StringIO()
        out_handle.write(">{}\n".format(rec.id))
        out_handle.write(str(rec.seq).upper())
        out_handle.write("\n")
        outfile.write(out_handle.getvalue().encode('utf-8'))

    def mappings_generator(self, chroms, read_length):
        genomeindex = self.config.abspath(self.config.genome.index)
        bowtie_command = "bowtie -S -t -v 0 -m 1 -f {} -".format(
            genomeindex
        )
        args = shlex.split(bowtie_command)

        with Popen(args, stdin=PIPE, stdout=PIPE) as bowtie:

            def push_genome_reads():
                for rec in self.generate_reads(chroms, read_length):
                    self.write_fasta_read(bowtie.stdin, rec)
                bowtie.stdin.close()

            genome_reads = Process(target=push_genome_reads)
            genome_reads.start()

            infile = pysam.AlignmentFile(  # @UndefinedVariable
                bowtie.stdout, 'r')

            try:
                for mapping in infile.fetch():
                    yield mapping
            finally:
                print("mappings generator cleanup")
                genome_reads.terminate()
                bowtie.stdin.close()
                bowtie.stdout.close()

    def mappable_regions_filename(self):
        filename = os.path.join(
            self.config.bins.cache_dir,
            self.config.bins.mappable_regions
        )
        return self.config.abspath(filename)

    def mappable_positions_count_filename(self):
        filename = os.path.join(
            self.config.bins.cache_dir,
            self.config.bins.mappable_positions_count
        )
        return self.config.abspath(filename)

    def chrom_sizes_filename(self):
        filename = os.path.join(
            self.config.bins.cache_dir,
            self.config.bins.chrom_sizes
        )
        filename = self.config.abspath(filename)
        return filename

    def calc_chrom_mappable_positions_count(self):
        filename = self.mappable_regions_filename()
        assert os.path.exists(filename)

        result = defaultdict(lambda: 0)
        with open(filename, 'r') as infile:
            for line in infile.readlines():
                row = [r.strip() for r in line.strip().split('\t')]
                result[row[0]] += int(row[2]) - int(row[1])
        return Box(result)

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
        filename = self.chrom_sizes_filename()
        if not os.path.exists(filename):
            result = self.calc_chrom_sizes()
            result.to_yaml(filename)

        with open(filename, 'r') as infile:
            result = Box.from_yaml(infile)
            return result

    def chrom_mappable_positions_count(self):
        filename = self.mappable_positions_count_filename()
        if not os.path.exists(filename):
            result = self.calc_chrom_mappable_positions_count()
            result.to_yaml(filename)
        with open(filename, 'r') as infile:
            result = Box.from_yaml(infile)
            return result
