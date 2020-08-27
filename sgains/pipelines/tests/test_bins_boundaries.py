'''
Created on Jun 23, 2017

@author: lubo
'''
import pytest
import numpy as np
import pandas as pd

from sgains.configuration.parser import Config
from sgains.pipelines.bins_pipeline import BinsPipeline


# pytestmark = pytest.mark.xfail


@pytest.mark.parametrize("chromosome", ['chr1'])
def test_bins_boundaries_generator(tests_config, bin_boundaries, chromosome):

    pipeline = BinsPipeline(tests_config)

    mappable_regions_df = pipeline.load_mappable_regions()

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


@pytest.mark.parametrize("chromosome", ['chr1'])
def test_bins_boundaries(tests_config, bin_boundaries, chromosome):
    pipeline = BinsPipeline(tests_config)

    fixture_df = bin_boundaries[bin_boundaries['bin.chrom'] == chromosome]
    fixture_df = fixture_df[[
        'bin.chrom', 'bin.start', 'bin.start.abspos',
        'bin.end', 'bin.length', 'mappable.positions']]
    fixture_df = fixture_df.reset_index()

    regions_df = pipeline.load_mappable_regions()
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


def test_bins_boundaries_bwa_R100():
    config = Config.parse(
        "/home/lubo/Work/single-cell/data/demo_data/bwa_hg19/"
        "sgains_bwa_10x_bj_mkn45_10pct_varbin_10x.yml")
    assert config is not None

    print(config.config.bins)
    bins_filename = config.bins_boundaries_filename(chrom="chr2")
    print(bins_filename)

    mappable_regions_filename = config.mappable_regions_filename(chrom="chr2")
    print(mappable_regions_filename)

    regions_df = pd.read_csv(
        mappable_regions_filename,
        names=['chrom', 'start_pos', 'end_pos'],
        sep='\t')
    regions_df = regions_df.sort_values(
        by=['chrom', 'start_pos', 'end_pos'])

    print(regions_df.head())
    print(regions_df.dtypes)

    pipeline = BinsPipeline(config)
    bins_df = pipeline.calc_bins_boundaries(["chr2"], regions_df)
    assert bins_df is not None

    print(bins_df.head())
