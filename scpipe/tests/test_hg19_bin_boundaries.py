'''
Created on Jun 23, 2017

@author: lubo
'''
import pytest
import numpy as np

from hg19 import HumanGenome19


@pytest.mark.parametrize("chromosome", HumanGenome19.CHROMS)
def test_bin_boundaries_generator(hg, bin_boundaries, chromosome):

    mappable_regions_df = hg.load_mappable_regions()

    chrom_df = bin_boundaries[bin_boundaries['bin.chrom'] == chromosome]
    for index, mappable_bin in enumerate(hg.bin_boundaries_generator(
            [chromosome], mappable_regions_df)):
        fixture_bin = chrom_df.iloc[index, :]

        assert mappable_bin.chrom == fixture_bin['bin.chrom']
        assert mappable_bin.chrom == fixture_bin['bin.chrom']
        assert mappable_bin.current_size == \
            fixture_bin['mappable.positions'], \
            (index, mappable_bin, fixture_bin)
        assert mappable_bin.start_pos == fixture_bin['bin.start']
        assert mappable_bin.end_pos == \
            fixture_bin['bin.end'], \
            (index, mappable_bin, fixture_bin)


@pytest.mark.parametrize("chromosome", HumanGenome19.CHROMS)
def test_bin_boundaries(hg, bin_boundaries, chromosome):
    fixture_df = bin_boundaries[bin_boundaries['bin.chrom'] == chromosome]
    fixture_df = fixture_df.reset_index()

    regions_df = hg.load_mappable_regions()
    bins_df = hg.calc_bin_boundaries([chromosome], regions_df)

    df = bins_df[bins_df['bin.chrom'] == chromosome]
    df = df.reset_index()

    assert np.all(df.columns == fixture_df.columns)
    assert np.all(df['bin.start'] == fixture_df['bin.start'])
    assert np.all(df['bin.start.abspos'] == fixture_df['bin.start.abspos'])
    assert np.all(df['bin.end'] == fixture_df['bin.end'])
    assert np.all(df['bin.length'] == fixture_df['bin.length'])
    assert np.all(df['mappable.positions'] == fixture_df['mappable.positions'])
