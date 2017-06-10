'''
Created on Jun 10, 2017

@author: lubo
'''
from config import load_config


def test_load_config():
    config = load_config("scpipe.yml")
    assert config is not None

    assert config.genome is not None
    assert config.genome.src == "data/hg19_safe"
    assert config.genome.dst == "data/hg19"
