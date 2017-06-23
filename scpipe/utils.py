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


class MappableBin(object):

    def __init__(self, prev=None, start_pos=None,
                 chrom=None, bin_size=None, chrom_abspos=0):
        if prev is None:
            self.chrom = chrom
            self.expected_size = bin_size
            self.chrom_abspos = chrom_abspos
            self.start_pos = 0
            self.end_pos = 0
            self.mappable_positions = 0
        else:
            assert prev is not None
            assert start_pos is not None

            self.chrom = prev.chrom
            self.chrom_abspos = prev.chrom_abspos
            self.start_pos = start_pos
            self.end_pos = prev.end_pos
            self.mappable_positions = self.end_pos - self.start_pos
            self.expected_size = prev.expected_size

    def check_extend(self, region):
        assert region['start_pos'] >= self.end_pos
        assert region['end_pos'] > region['start_pos']

        region_size = region['end_pos'] - region['start_pos']
        if region_size < self.missing_mappable_positions():
            print("check_extend1: ", self, 
                  region['start_pos'], region['end_pos'], region_size)
            self.end_pos = region['end_pos']
            self.mappable_positions += region_size
            print("check_extend2: ", self)
            return True
        return False

    def split_extend(self, region):
        assert region['start_pos'] >= self.end_pos
        assert region['end_pos'] > region['start_pos']

        region_size = region['end_pos'] - region['start_pos']
        missing_mappable_positions = self.missing_mappable_positions()
        print("region size: ", region_size,
              "missing_mappable_positions:", missing_mappable_positions)

        assert region_size >= missing_mappable_positions

        self.end_pos = region['start_pos'] + missing_mappable_positions
        self.mappable_positions += missing_mappable_positions

        next_bin = MappableBin(prev=self, start_pos=self.end_pos)
        next_bin.end_pos = region['end_pos']
        next_bin.mappable_positions = region_size - missing_mappable_positions

#         assert next_bin.mappable_positions == \
#             region_size - missing_mappable_positions

        return next_bin

    def adapt_excess(self, excess):
        if excess >= 1:
            self.expected_size += 1
            excess -= 1
        return excess

    def missing_mappable_positions(self):
        return self.expected_size - self.mappable_positions

    def is_full(self):
        return self.mappable_positions >= self.expected_size

    def __repr__(self):
        return "{}\t{}\t{}\t{}".format(
            self.chrom,
            self.start_pos,
            self.end_pos,
            self.mappable_positions
        )
