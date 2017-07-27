'''
Created on Jun 29, 2017

@author: lubo
'''
import numpy as np
import pytest
from hg19 import HumanGenome19


@pytest.mark.parametrize("chromosome", HumanGenome19.CHROMS)
def test_bin_gc_content(hg, chromosome, bin_boundaries, gc_bin_boundaries):

    fixture_df = gc_bin_boundaries[
        gc_bin_boundaries['bin.chrom'] == chromosome]
    assert fixture_df is not None
    fixture_df = fixture_df.reset_index()

    df = hg.calc_bins_gc_content([chromosome], bin_boundaries)
    assert df is not None

    assert np.all(np.abs(df['gc.content'] -
                         fixture_df['gc.content']) < 1E-5)


@pytest.mark.parametrize("chromosome", HumanGenome19.CHROMS)
def test_bin_boundaries_full(hg, gc_bin_boundaries, chromosome):
    fixture_df = gc_bin_boundaries[
        gc_bin_boundaries['bin.chrom'] == chromosome
    ]
    fixture_df = fixture_df.reset_index()

    regions_df = hg.load_mappable_regions()
    bins_df = hg.calc_bins_boundaries([chromosome], regions_df)
    df = hg.calc_bins_gc_content([chromosome], bins_df)

    # assert np.all(df.columns == fixture_df.columns)
    assert np.all(df['bin.start'] == fixture_df['bin.start'])
    assert np.all(df['bin.start.abspos'] == fixture_df['bin.start.abspos'])
    assert np.all(df['bin.end'] == fixture_df['bin.end'])
    assert np.all(df['bin.length'] == fixture_df['bin.length'])
    assert np.all(df['mappable.positions'] == fixture_df['mappable.positions'])

    assert np.all(np.abs(df['gc.content'] -
                         fixture_df['gc.content']) < 1E-5)
