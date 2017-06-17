'''
Created on Jun 10, 2017

@author: lubo
'''
import os

from config import Config


def test_load_config():
    conf = Config.load("scpipe.yml")
    assert conf is not None

    assert conf.genome is not None
    assert os.path.basename(conf.genome.src) == "hg19_safe"
    assert os.path.basename(conf.genome.dst) == "hg19"

    assert os.path.isabs(conf.abspath(conf.genome.src))
    assert os.path.isabs(conf.abspath(conf.genome.src))


def test_bins_config():
    conf = Config.load("scpipe.yml")
    assert conf is not None

    assert conf.bins is not None
    assert conf.bins.reads_length == 50
    assert conf.bins.bins_count == 10000
    assert conf.bins.cache_dir == "data/10k.50"
