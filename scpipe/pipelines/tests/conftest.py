'''
Created on Jun 22, 2017

@author: lubo
'''
import os

import pytest

from config import Config
from hg19 import HumanGenome19
import pandas as pd


@pytest.fixture(scope='session')
def tests_config():
    config = Config.load("tests/data/scpipe_tests.yml")
    return config


@pytest.fixture(scope='session')
def hg():
    config = Config.load("tests/data/scpipe_tests.yml")
    return HumanGenome19(config)


@pytest.fixture(scope='session')
def bin_boundaries():
    bins_boundaries_fixture = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "data/bin.boundaries.bowtie.txt"
    )
    df = pd.read_csv(bins_boundaries_fixture, sep='\t')
    return df


@pytest.fixture(scope='session')
def gc_bin_boundaries():
    gc_bins_boundaries_fixture = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "data/varbin.gc.content.bowtie.txt"
    )
    df = pd.read_csv(gc_bins_boundaries_fixture, sep='\t')
    return df


@pytest.fixture(scope='session')
def varbin_counts():
    fixture_filename = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "data/varbin.txt"
    )
    df = pd.read_csv(fixture_filename, sep='\t')
    return df


@pytest.fixture(scope='session')
def varbin0918():
    fixture_filename = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "data/CJA0918.varbin.txt"
    )
    df = pd.read_csv(fixture_filename, sep='\t')
    return df
