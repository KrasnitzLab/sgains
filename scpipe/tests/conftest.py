'''
Created on Jun 22, 2017

@author: lubo
'''
import pytest
from config import Config
from hg19 import HumanGenome19


@pytest.fixture(scope='session')
def hg():
    config = Config.load("scpipe_tests.yml")
    return HumanGenome19(config)
