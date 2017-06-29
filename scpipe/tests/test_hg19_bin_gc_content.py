'''
Created on Jun 29, 2017

@author: lubo
'''
import os
import pandas as pd
import pytest


@pytest.fixture(scope='session')
def bin_boundaries():
    bins_boundaries_fixture = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "data/bin.boundaries.bowtie.txt"
    )
    df = pd.read_csv(bins_boundaries_fixture, sep='\t')
    return df


def test_bin_gc_content(hg, bin_boundaries):
    hg.bins_gc_content(['chr1'], bin_boundaries)
