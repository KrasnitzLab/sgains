'''
Created on Jun 10, 2017

@author: lubo
'''
import pytest

# pytestmark = pytest.mark.xfail


def test_hg19_simple(hg):

    assert hg is not None

    chr_y = hg.load_chrom("chrY")
    assert chr_y.id == "chrY"
    assert len(chr_y) == 59373566


def test_mask_chrY(hg):
    rec = hg.mask_chrY_pars()
    assert rec.id == "chrY"
    assert len(rec) == 59373566
    hg.save_chrom(rec, "chrY.psr.test")


def test_count_chrom_mappable_positions(hg):
    result = hg.calc_chrom_mappable_positions_count()

    assert result['chr1'] == 217027151  # 217026582
    assert result['chr10'] == 126701497  # 126701210


def test_chrom_sizes(hg):
    result = hg.chrom_sizes()

    assert result.chr1.size == 249250621
    assert result.chr1.abspos == 0

    assert result.chr10.size == 135534747
    assert result.chr10.abspos == 1680373143


def test_chrom_mappable_positions_count(hg):
    result = hg.chrom_mappable_positions_count()

    assert result['chr1'] == 217027151  # 217026582
    assert result['chr10'] == 126701497  # 126701210


def test_mappable_positions_total_count(hg):
    total = hg.total_mappable_positions_count()

    assert total == 2764542640  # 2758913563


def test_calc_chrom_bins(hg):
    chrom_bins = hg.calc_chrom_bins()

    assert chrom_bins['chr8'].bins_count == 505
    assert chrom_bins['chr8'].bin_size == pytest.approx(276248.417, rel=1e-3)


def test_chr_bins(hg):
    chrom_bins = hg.chrom_bins()
    assert chrom_bins is not None
    chrom_bin = chrom_bins['chr1']

    assert chrom_bin.chrom == 'chr1'
