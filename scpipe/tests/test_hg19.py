'''
Created on Jun 10, 2017

@author: lubo
'''
from hg19 import HumanGenome19
import pytest
from config import Config


@pytest.fixture(scope='session')
def hg(request):
    config = Config.load("scpipe.yml")
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


def test_psr_chrY(hg):
    chr_y = hg.load_chrom("chrY.psr", hg.config.genome.dst)
    chr_ty = hg.load_chrom("chrY.psr.test", hg.config.genome.dst)

    assert chr_y.seq == chr_ty.seq


def test_generate_reads(hg):

    for num, rec in enumerate(hg.generate_reads(['chr1'], 100)):
        if num >= 100:
            break
        print(rec.id, len(rec))


def test_count_chrom_mappable_positions(hg):
    result = hg.count_chrom_mappable_positions()

    assert result['chr1'] == 217026582
    assert result['chr10'] == 126701210


def test_chrom_sizes(hg):
    result = hg.chrom_sizes()

    assert result.chr1.size == 249250621
    assert result.chr1.abspos == 0

    assert result.chr10.size == 135534747
    assert result.chr10.abspos == 1680373143
