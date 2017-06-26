'''
Created on Jun 23, 2017

@author: lubo
'''


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


class BinParams:
    def __init__(self):
        self.chrom_bin = None
        self.chrom_size = None

        self.chrom = None
        self.chrom_abspos = None

        self.bin_count = None
        self.bin_size = None
        self.bin_size_excess = None

    @staticmethod
    def build(*, chrom_size=None, chrom_bin=None):
        params = BinParams()
        params.chrom_bin = chrom_bin
        params.chrom_size = chrom_size
        params.chrom_abspos = chrom_size.abspos

        params.chrom = chrom_bin.chrom
        params.bins_count = chrom_bin.bins_count
        params.bin_size = int(chrom_bin.bin_size)
        params.bin_size_excess = chrom_bin.bin_size - params.bin_size

        return params


class MappableBin(object):

    @staticmethod
    def from_start(params, start_pos):
        res = MappableBin()
        res.params = params
        res.start_pos = start_pos
        res.end_pos = start_pos
        res.current_size = 0
        res.bin_size = res.params.bin_size

        return res

    @staticmethod
    def from_prev(prev_bin, start_pos):
        res = MappableBin()
        res.params = prev_bin.params
        res.start_pos = start_pos
        res.end_pos = prev_bin.end_pos
        res.current_size = res.end_pos - res.start_pos
        res.bin_size = res.params.bin_size

        return res

    @property
    def chrom(self):
        return self.params.chrom

    @property
    def chrom_abspos(self):
        return self.params.chrom_size.abspos

    def __init__(self):
        self.start_pos = 0
        self.end_pos = 0
        self.current_size = 0
        self.bin_size = 0

    def check_extend(self, region):
        assert region['start_pos'] >= self.end_pos, (self, region)
        assert region['end_pos'] > region['start_pos']

        region_size = region['end_pos'] - region['start_pos']
        if region_size < self.missing_mappable_positions():
            self.end_pos = region['end_pos']
            self.current_size += region_size
            return True
        return False

    def split_extend(self, region):
        assert region['start_pos'] >= self.end_pos
        assert region['end_pos'] > region['start_pos']

        region_size = region['end_pos'] - region['start_pos']
        missing_mappable_positions = self.missing_mappable_positions()

        assert region_size >= missing_mappable_positions
        # assert region_size - missing_mappable_positions <= self.bin_size

        self.end_pos = region['start_pos'] + missing_mappable_positions
        self.current_size += missing_mappable_positions

        next_bin = MappableBin.from_prev(self, start_pos=self.end_pos)
        next_bin.end_pos = region['end_pos']
        next_bin.current_size = region_size - missing_mappable_positions
        return next_bin

    def overfill_split(self, current_excess):
        assert self.is_overfill()

        mappable_bins = []
        while self.current_size > 0:
            # current_excess += self.params.bin_size_excess
            mb = MappableBin.from_prev(self, start_pos=self.start_pos)
            current_excess = mb.adapt_excess(current_excess)

            current_size = self.current_size
            if current_size > mb.bin_size:
                current_size = mb.bin_size
            mb.current_size = current_size
            mb.end_pos = mb.start_pos + current_size

            mappable_bins.append(mb)

            self.current_size -= mb.bin_size
            self.start_pos = mb.end_pos
        return current_excess, mappable_bins

    def adapt_excess(self, excess):
        excess += self.params.bin_size_excess
        # print("start_pos:", self.start_pos, "excess: ", excess)
        if excess >= 1:
            self.bin_size += 1
            excess -= 1
        # print("after excess: ", excess, "bin_size: ", self.bin_size)
        return excess

    def missing_mappable_positions(self):
        return self.bin_size - self.current_size

    def is_full(self):
        return self.current_size == self.bin_size

    def is_overfill(self):
        return self.current_size > self.bin_size

    def __repr__(self):
        return "{}\t{}\t{}\t{}\t{}".format(
            self.chrom,
            self.start_pos,
            self.end_pos,
            self.current_size,
            self.bin_size
        )
