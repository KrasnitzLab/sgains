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
from utils import MappableState, Mapping, MappableRegion, \
    MappableBin, BinParams


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
        self._chrom_sizes = None
        self._chrom_bins = None
        self._chrom_mappable_positions_count = None

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

    CHRY_PAR1 = [
        (10000, 2649520),
        (59034049, 59363566)
    ]

    def mask_chrY_pars(self):
        chr_y = self.load_chrom("chrY")

        masked = chr_y.seq.tomutable()
        for par in self.CHRY_PAR1:
            start, end = par
            masked[start:end] = 'N' * (end - start)

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

    async def async_mappings_generator(self, reads_generator, bowtie):
        writer = asyncio.Task(
            self.async_write_reads_generator(bowtie.stdin, reads_generator)
        )

        while True:
            line = await bowtie.stdout.readline()
            if not line:
                break
            yield line.decode()

        await bowtie.wait()
        await writer

    async def async_generate_mappings(
            self, chroms, read_length, outfile=None):
        if outfile is None:
            outfile = sys.stdout

        bowtie = await self.async_start_bowtie()
        reads_generator = self.generate_reads(chroms, read_length)
        async for mappings in self.async_mappings_generator(
                reads_generator, bowtie):
            outfile.write(mappings)

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

    def calc_chrom_mappable_positions_count(self):
        filename = self.config.mappable_regions_filename()
        assert os.path.exists(filename)

        result = defaultdict(lambda: 0)
        with open(filename, 'r') as infile:
            for line in infile.readlines():
                row = [r.strip() for r in line.strip().split('\t')]
                result[row[0]] += int(row[2]) - int(row[1])
        return Box(result)

#     def mappable_regions(self):
#         filename = self.mappable_regions_filename()
#         assert os.path.exists(filename)
#         result = []
#         with open(filename, 'r') as infile:
#             for line in infile.readlines():
#                 row = [r.strip() for r in line.strip().split('\t')]
#                 result.append([row[0], int(row[1]), int(row[2])])

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

    def chrom_bins(self):
        if self._chrom_bins is None:
            self._chrom_bins = self.calc_chrom_bins()
        return self._chrom_bins

    def load_mappable_regions(self):
        df = pd.read_csv(
            self.config.mappable_regions_filename(),
            names=['chrom', 'start_pos', 'end_pos'],
            sep='\t')
        df.sort_values(by=['chrom', 'start_pos', 'end_pos'], inplace=True)

        return df

    def calc_bin_boundaries(self, chroms, mappable_regions_df=None):
        chrom_sizes = self.chrom_sizes()
        chrom_bins = self.calc_chrom_bins()

        if mappable_regions_df is None:
            mappable_regions_df = self.load_mappable_regions()

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
