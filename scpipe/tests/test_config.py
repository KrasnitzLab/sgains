'''
Created on Jun 10, 2017

@author: lubo
'''
import os

from config import load_config


def test_load_config():
    config = load_config("scpipe.yml")
    assert config is not None

    assert config.genome is not None
    assert os.path.basename(config.genome.src) == "hg19_safe"
    assert os.path.basename(config.genome.dst) == "hg19"
