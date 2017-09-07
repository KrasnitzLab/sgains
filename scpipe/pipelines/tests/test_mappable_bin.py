'''
Created on Jun 26, 2017

@author: lubo
'''
from utils import MappableBin, BinParams
import pytest


@pytest.fixture
def chr1_bin(hg):
    chrom_bin = hg.chrom_bins()['chr1']
    return chrom_bin


@pytest.fixture
def chr1_size(hg):
    chrom_size = hg.chrom_sizes()['chr1']
    return chrom_size


@pytest.fixture
def chr1_params(chr1_bin, chr1_size):
    params = BinParams.build(chrom_size=chr1_size, chrom_bin=chr1_bin)
    assert params is not None
    return params


def test_mappable_bin_from_region(chr1_bin, chr1_params):

    mb = MappableBin.from_start(chr1_params, 0)

    assert mb.chrom == 'chr1'
    assert mb.start_pos == 0
    assert mb.end_pos == 0
    assert mb.current_size == 0
    assert mb.bin_size == chr1_params.bin_size


def test_check_extend(chr1_bin, chr1_params):
    mb = MappableBin.from_start(chr1_params, 0)
    region = {
        'start_pos': 10,
        'end_pos': 20,
    }

    assert mb.check_extend(region)
    assert mb.start_pos == 0
    assert mb.end_pos == 20
    assert mb.current_size == 10


def test_check_extend_overflow(chr1_params):
    chr1_params.bin_size = 10000
    mb = MappableBin.from_start(chr1_params, 0)
    region = {
        'start_pos': 10,
        'end_pos': 10011,
    }
    assert not mb.check_extend(region)

    assert mb.start_pos == 0
    assert mb.end_pos == 0
    assert mb.current_size == 0
    assert mb.bin_size == 10000


def test_split_extend(chr1_params):
    chr1_params.bin_size = 10000
    mb = MappableBin.from_start(chr1_params, 0)
    region = {
        'start_pos': 10,
        'end_pos': 10011,
    }
    next_bin = mb.split_extend(region)

    assert mb.bin_size == 10000
    assert mb.start_pos == 0
    assert mb.end_pos == 10010
    assert mb.current_size == 10000

    assert next_bin.bin_size == 10000
    assert next_bin.start_pos == 10010
    assert next_bin.end_pos == 10011
    assert next_bin.current_size == 1


def test_split_extend_overfill(chr1_params):
    chr1_params.bin_size = 10000
    mb = MappableBin.from_start(chr1_params, 0)
    region = {
        'start_pos': 10,
        'end_pos': 20011,
    }
    next_bin = mb.split_extend(region)

    assert mb.is_full()
    assert not mb.is_overfill()

    assert not next_bin.is_full()
    assert next_bin.is_overfill()


def test_overfill_split(chr1_params):
    chr1_params.bin_size = 10000
    chr1_params.bin_size_excess = 0.1

    mb = MappableBin.from_start(chr1_params, 0)
    mb.end_pos = 20001
    mb.current_size = 20001

    current_excess = 0.1

    current_excess, mappable_bins = \
        mb.overfill_split(current_excess)

    assert len(mappable_bins) == 3
    assert current_excess == 0.4

    assert all([mb.is_full() for mb in mappable_bins[0:2]])
    assert all([mb.current_size == 10000 for mb in mappable_bins[0:2]])

    last_mb = mappable_bins[-1]
    assert last_mb.current_size == 1
    assert last_mb.bin_size == 10000
    assert last_mb.start_pos == 20000
    assert last_mb.end_pos == 20001
