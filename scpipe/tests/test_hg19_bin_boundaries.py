'''
Created on Jun 23, 2017

@author: lubo
'''
import pandas as pd
import os


def test_calc_bin_boundaries(hg):
    bins_boundaries_fixture = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "data/bin.boundaries.bowtie.txt"
    )
    df = pd.read_csv(bins_boundaries_fixture, sep='\t')
    assert df is not None

    for index, mappable_bin in enumerate(hg.calc_bin_boundaries()):
        print(index, mappable_bin)

        fixture_bin = df.iloc[index, :]

        assert mappable_bin.chrom == fixture_bin['chrom']
        assert mappable_bin.mappable_positions == fixture_bin['mappable.positions']
        assert mappable_bin.start_pos == fixture_bin['bin.start.chrompos']
        assert mappable_bin.end_pos == fixture_bin['bin.end.chrompos']
