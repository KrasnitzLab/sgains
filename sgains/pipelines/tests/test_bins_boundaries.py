'''
Created on Jun 23, 2017

@author: lubo
'''
import pytest
import numpy as np

from sgains.pipelines.bins_pipeline import BinsPipeline


# pytestmark = pytest.mark.xfail


@pytest.mark.parametrize("chromosome", ['chr1'])
def test_bins_boundaries_generator(tests_config, bin_boundaries, chromosome):

    pipeline = BinsPipeline(tests_config)

    mappable_regions_df = pipeline.hg.load_mappable_regions()

    chrom_df = bin_boundaries[bin_boundaries['bin.chrom'] == chromosome]

    for index, mappable_bin in enumerate(pipeline.bins_boundaries_generator(
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


@pytest.mark.parametrize("chromosome", ['chr1'])  # HumanGenome19.CHROMS)
def test_bins_boundaries(tests_config, bin_boundaries, chromosome):
    pipeline = BinsPipeline(tests_config)

    fixture_df = bin_boundaries[bin_boundaries['bin.chrom'] == chromosome]
    fixture_df = fixture_df[[
        'bin.chrom', 'bin.start', 'bin.start.abspos',
        'bin.end', 'bin.length', 'mappable.positions']]
    fixture_df = fixture_df.reset_index()

    regions_df = pipeline.hg.load_mappable_regions()
    bins_df = pipeline.calc_bins_boundaries([chromosome], regions_df)

    df = bins_df[bins_df['bin.chrom'] == chromosome]
    df = df.reset_index()

    print(fixture_df.head())
    print(df.head())
    
    assert np.all(df.columns == fixture_df.columns)
    assert np.all(df['bin.start'] == fixture_df['bin.start'])
    assert np.all(df['bin.start.abspos'] == fixture_df['bin.start.abspos'])
    assert np.all(df['bin.end'] == fixture_df['bin.end'])
    assert np.all(df['bin.length'] == fixture_df['bin.length'])
    assert np.all(df['mappable.positions'] == fixture_df['mappable.positions'])
