'''
Created on Jun 23, 2017

@author: lubo
'''
from utils import MappableBin
import pytest

CHROM1_BINSIZE = 276115
CHROM1_BIN1_END = 932957
CHROM1_BIN5_END = 1822513


@pytest.fixture
def chr1_bin1():
    bin1 = MappableBin(start_pos=0, chrom='chr1', bin_size=CHROM1_BINSIZE)
    bin1.end_pos = 932112
    bin1.mappable_positions = 275309

    return bin1


@pytest.fixture
def chr1_bin5():
    bin1 = MappableBin(
        start_pos=1822513, chrom='chr1', bin_size=CHROM1_BINSIZE)
    bin1.end_pos = 1817060
    bin1.mappable_positions = 270674

    return bin1


def test_mappable_bin_split_extend_chr1_bin1(chr1_bin1):
    bin1 = chr1_bin1

    assert bin1.missing_mappable_positions() == \
        (CHROM1_BINSIZE - bin1.mappable_positions)

    region = {
        'start_pos': 932151,
        'end_pos': 943123,
    }
    assert not bin1.check_extend(region)

    bin2 = bin1.split_extend(region)

    assert bin1.missing_mappable_positions() == 0
    assert bin1.mappable_positions == CHROM1_BINSIZE

    assert bin1.end_pos == CHROM1_BIN1_END
    assert bin2.start_pos == CHROM1_BIN1_END
    assert bin2.expected_size == CHROM1_BINSIZE
    assert bin2.end_pos == region['end_pos']


def test_mappable_bin_split_extend_chr1_bin5(chr1_bin5):
    bin1 = chr1_bin5

    assert bin1.missing_mappable_positions() == \
        (CHROM1_BINSIZE - bin1.mappable_positions)

    region = {
        'start_pos': 1817072,
        'end_pos': 1828506,
    }

    assert not bin1.check_extend(region)

    bin2 = bin1.split_extend(region)

    assert bin1.missing_mappable_positions() == 0
    assert bin1.end_pos == CHROM1_BIN5_END

    assert bin2.start_pos == CHROM1_BIN5_END
    assert bin2.expected_size == CHROM1_BINSIZE
    assert bin2.end_pos == region['end_pos']


def test_adapt_excess(chr1_bin1):
    excess = 0.5
    assert chr1_bin1.adapt_excess(excess) == 0.5
    assert chr1_bin1.expected_size == CHROM1_BINSIZE

    excess = 1.5
    assert chr1_bin1.adapt_excess(excess) == 0.5
    assert chr1_bin1.expected_size == CHROM1_BINSIZE + 1


def test_check_extend():
    bin1 = MappableBin(start_pos=0, chrom='chr1', bin_size=CHROM1_BINSIZE)
    bin1.end_pos = 931936
    bin1.mappable_positions = 275191

    region = {
        'start_pos': 931955,
        'end_pos': 931981,
    }

    assert bin1.check_extend(region)
    assert bin1.mappable_positions == 275217
    assert bin1.end_pos == 931981
