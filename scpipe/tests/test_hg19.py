'''
Created on Jun 10, 2017

@author: lubo
'''
from hg19 import HumanGenome19
from config import load_config
import pytest


@pytest.fixture(scope='session')
def hg(request):
    config = load_config("scpipe.yml")
    return HumanGenome19(config)


def test_hg19_simple(hg):

    assert hg is not None

    chr_y = hg.load_chrom("chrY")
    assert chr_y.id == "chrY"
    assert len(chr_y) == 59373566


def test_mask_chrY(hg):
    rec = hg.mask_pseudoautosomal_chrY()
    assert rec.id == "chrY"
    assert len(rec) == 59373566
    hg.save_chrom(rec, "chrY.psr.test")
