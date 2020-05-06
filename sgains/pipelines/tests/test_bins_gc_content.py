'''
Created on Jun 29, 2017

@author: lubo
'''
import numpy as np
import pandas as pd
import pytest
from sgains.pipelines.bins_pipeline import BinsPipeline


# pytestmark = pytest.mark.xfail


@pytest.mark.parametrize("chromosome", ['chr1'])
def test_bin_gc_content(
        tests_config, chromosome, bin_boundaries):

    pipeline = BinsPipeline(tests_config)

    fixture_df = bin_boundaries[
        bin_boundaries['bin.chrom'] == chromosome]
    print(len(fixture_df))
    assert fixture_df is not None
    fixture_df = fixture_df.reset_index()
    print("len(fixture_df) 2:", len(fixture_df))

    df = pipeline.calc_bins_gc_content([chromosome], bin_boundaries)
    assert df is not None
    print("len(fixture_df) 3:", len(fixture_df), "len(df):", len(df))

    print(fixture_df.head())
    print(df.head())

    pd.testing.assert_series_equal(
        df['gc.content'], fixture_df['gc.content'])


@pytest.mark.parametrize("chromosome", ['chr1'])
def test_bin_boundaries_full(tests_config, bin_boundaries, chromosome):
    fixture_df = bin_boundaries[
        bin_boundaries['bin.chrom'] == chromosome
    ]
    fixture_df = fixture_df.reset_index()

    pipeline = BinsPipeline(tests_config)

    regions_df = pipeline.genome.load_mappable_regions()

    bins_df = pipeline.calc_bins_boundaries([chromosome], regions_df)
    df = pipeline.calc_bins_gc_content([chromosome], bins_df)

    # assert np.all(df.columns == fixture_df.columns)
    assert np.all(df['bin.start'] == fixture_df['bin.start'])
    assert np.all(df['bin.start.abspos'] == fixture_df['bin.start.abspos'])
    assert np.all(df['bin.end'] == fixture_df['bin.end'])
    assert np.all(df['bin.length'] == fixture_df['bin.length'])
    assert np.all(df['mappable.positions'] == fixture_df['mappable.positions'])

    pd.testing.assert_series_equal(
        df['gc.content'], fixture_df['gc.content'])
