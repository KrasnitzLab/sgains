'''
Created on Jul 6, 2017

@author: lubo
'''
from hg19 import HumanGenome19
import pytest


@pytest.mark.parametrize("chromosome", ['chr1'])  # HumanGenome19.CHROMS)
def test_bin_counts_simple(hg, varbin_counts, chromosome):

    print(varbin_counts.head())

    hg.bin_count("CJA1247.rmdup.bam")
