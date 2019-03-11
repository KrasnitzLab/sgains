'''
Created on Jun 10, 2017

@author: lubo
'''
import os

from sgains.config import Config
import pickle


def test_load_config():
    conf = Config.load("tests/data/scpipe_tests.yml")
    assert conf is not None

    assert conf.genome is not None
    assert os.path.basename(conf.genome.version) == "hg19"


def test_bins_config():
    conf = Config.load("tests/data/scpipe_tests.yml")
    assert conf is not None

    assert conf.mappable_regions is not None
    assert conf.mappable_regions.length == 100

    assert conf.bins is not None
    assert conf.bins.bins_count == 10000
    assert conf.bins.bins_dir == "test_data/R100_B10k"


def test_pickle_depickle():
    conf = Config.load("tests/data/scpipe_tests.yml")
    assert conf is not None

    pkl = pickle.dumps(conf)
    assert pkl is not None

    res = pickle.loads(pkl)
    assert res is not None

    assert isinstance(res, Config)
    assert conf.to_dict() == res.to_dict()
