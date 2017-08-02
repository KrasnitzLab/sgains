'''
Created on Jun 10, 2017

@author: lubo
'''
import os

from config import Config


def test_load_config():
    conf = Config.load("scpipe_tests.yml")
    assert conf is not None

    assert conf.genome is not None
    assert os.path.basename(conf.genome.version) == "hg19"


def test_bins_config():
    conf = Config.load("scpipe_tests.yml")
    assert conf is not None

    assert conf.mappable_regions is not None
    assert conf.mappable_regions.length == 100

    assert conf.bins is not None
    assert conf.bins.bins_count == 10000
    assert conf.bins.work_dir == "data/R100_B10k"
