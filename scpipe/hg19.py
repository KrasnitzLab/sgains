'''
Created on Jun 10, 2017

@author: lubo
'''
import asyncio
from collections import defaultdict
import io
import os
import sys

from Bio import SeqIO  # @UnresolvedImport
from Bio.SeqRecord import SeqRecord  # @UnresolvedImport
from box import Box

import pandas as pd


class Mapping(object):
    def __init__(self, name=None, flag=None, chrom=None, start=None):
        self.flag = flag
        self.chrom = chrom
        self.start = start
        self.name = name

    @staticmethod
    def parse_sam(line):
        row = line.rstrip().split("\t")
        name = row[0]
        flag = int(row[1])
        chrom = row[2]
        start = int(row[3])
        return Mapping(name=name, flag=flag, chrom=chrom, start=start)


class MappableRegion(object):
    def __init__(self, mapping=None):
        self.flag = mapping.flag
        self.chrom = mapping.chrom
        self.start = mapping.start

        self.end = self.start + 1

    def extend(self, start):
        self.end = start + 1

    def __repr__(self):
        return "{}\t{}\t{}".format(
            self.chrom, self.start, self.end)


class MappableState(object):
    OUT = 0
    IN = 1


class MappableBin(object):

    def __init__(self, prev=None, start_pos=None,
                 chrom=None, bin_size=None, chrom_abspos=0):
        if prev is None:
            self.chrom = chrom
            self.expected_size = bin_size
            self.chrom_abspos = chrom_abspos
            self.start_pos = 0
            self.end_pos = 0
            self.mappable_possitions = self.end_pos - self.start_pos
        else:
            assert prev is not None
            assert start_pos is not None

            self.chrom = prev.chrom
            self.chrom_abspos = prev.chrom_abspos
            self.start_pos = start_pos
            self.end_pos = prev.end_pos
            self.mappable_possitions = self.end_pos - self.start_pos
            self.expected_size = prev.expected_size

    def check_extend(self, region):
        assert region['start_pos'] >= self.end_pos
        assert region['end_pos'] > region['start_pos']

        region_size = region['end_pos'] - region['start_pos']
        if region_size < self.missing_mappable_positions():
            self.end_pos = region['end_pos']
            self.mappable_possitions += region_size
            return True

    def split_extend(self, region):
        assert region['start_pos'] >= self.end_pos
        assert region['end_pos'] > region['start_pos']

        region_size = region['end_pos'] - region['start_pos']
        missing_mappable_positions = self.missing_mappable_positions()
        assert region_size >= missing_mappable_positions

        self.end_pos += missing_mappable_positions
        self.mappable_possitions += missing_mappable_positions

        next_bin = MappableBin(prev=self, start_pos=self.end_pos)
        next_bin.end_pos += region_size - missing_mappable_positions

        return next_bin

    def excess(self):
        return self.size - self.expected_size

    def missing_mappable_positions(self):
        return self.expected_size - self.mappable_possitions

    def is_full(self):
        return self.mappable_possitions >= self.expected_size

    def __repr__(self):
        return "{}\t{}\t{}\t{}".format(
            self.chrom,
            self.start_pos,
            self.end_pos,
            self.mappable_possitions
        )


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
            pass

    async def async_mappable_regions_generator(self, infile):
        prev = None
        state = MappableState.OUT

        while True:
            line = await infile.readline()
            if not line:
                break
            line = line.decode()
            if line[0] == '@':
                # comment
                continue

            mapping = Mapping.parse_sam(line)

            if state == MappableState.OUT:
                if mapping.flag == 0:
                    prev = MappableRegion(mapping)
                    state = MappableState.IN
            else:
                if mapping.flag == 0:
                    if mapping.chrom == prev.chrom:
                        prev.extend(mapping.start)
                    else:
                        yield prev
                        prev = MappableRegion(mapping)
                else:
                    yield prev
                    state = MappableState.OUT

        if state == MappableState.IN:
            yield prev

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

    @staticmethod
    async def async_write_fasta(outfile, rec):
        out = HumanGenome19.to_fasta_string(rec)
        outfile.write(out)
        await outfile.drain()

    async def async_start_bowtie(self):
        genomeindex = self.config.abspath(self.config.genome.index)
        create = asyncio.create_subprocess_exec(
            'bowtie', '-S', '-t', '-v', '0', '-m', '1',
            '-f', genomeindex, '-',
            stdin=asyncio.subprocess.PIPE,
            stdout=asyncio.subprocess.PIPE,
        )
        proc = await create
        return proc

    async def async_write_reads_generator(self, out, reads_generator):
        for rec in reads_generator:
            await self.async_write_fasta(out, rec)
        out.close()

    async def async_generate_mappable_regions(
            self, chroms, read_length, outfile=None):

        bowtie = await self.async_start_bowtie()
        reads_generator = self.generate_reads(chroms, read_length)
        writer = asyncio.Task(
            self.async_write_reads_generator(bowtie.stdin, reads_generator)
        )
        if outfile is None:
            outfile = sys.stdout
        async for mapping in self.async_mappable_regions_generator(
                bowtie.stdout):
            outfile.write(str(mapping))
            outfile.write('\n')
        await bowtie.wait()
        await writer

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

    def mappable_regions(self):
        filename = self.mappable_regions_filename()
        assert os.path.exists(filename)
        result = []
        with open(filename, 'r') as infile:
            for line in infile.readlines():
                row = [r.strip() for r in line.strip().split('\t')]
                result.append([row[0], int(row[1]), int(row[2])])

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

    def total_mappable_positions_count(self):
        counts = self.chrom_mappable_positions_count()
        total = 0
        for chrom in self.CHROMS:
            total += counts[chrom]
        return total

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

    def calc_bin_boundaries(self):
        chrom_sizes = self.chrom_sizes()
        chrom_mappable_positions_count = self.chrom_mappable_positions_count()
        total_mappable_positions_count = self.total_mappable_positions_count()

        chrom_bins = self.calc_chrom_bins()
        print(chrom_bins)

        df = pd.read_csv(
            self.mappable_regions_filename(),
            header=0,
            names=['chrom', 'start_pos', 'end_pos'],
            sep='\t')
        df.sort_values(by=['chrom', 'start_pos', 'end_pos'])

        current_excess = 0

        for chrom in self.CHROMS:
            mappable_bins = []
            chrom_df = df[df.chrom == chrom]
            bins_count = chrom_bins[chrom].bins_count
            bin_size = int(chrom_bins[chrom].bin_size)
            bin_size_excess = chrom_bins[chrom].bin_size - bin_size
            current_excess += bin_size_excess
            if current_excess >= 1.0:
                bin_size += 1
                current_excess -= 1.0

            print(bin_size)
            print(chrom_df.head())

            mappable_bin = None
            for index, row in chrom_df.iterrows():
                if mappable_bin is None:
                    mappable_bin = MappableBin(
                        chrom=chrom,
                        bin_size=bin_size,
                        chrom_abspos=chrom_sizes[chrom].abspos)
                if not mappable_bin.check_extend(row):
                    next_bin = mappable_bin.split_extend(row)
                    mappable_bins.append(mappable_bin)
                    print(mappable_bin)
                    if len(mappable_bins) >= 10:
                        break
                    mappable_bin = next_bin

            mappable_bin = None

            break
