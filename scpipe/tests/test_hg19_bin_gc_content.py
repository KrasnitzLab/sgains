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

    df = hg.bins_gc_content([chromosome], bin_boundaries)
    assert df is not None

    assert np.all(np.abs(df['gc.content'] -
                         fixture_df['gc.content']) < 1E-5)
