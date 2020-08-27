import os
import pytest

import pandas as pd

from sgains.configuration.parser import Config
from sgains.aligners import Hisat2, BWA
from sgains.genome import Genome


@pytest.fixture(scope='session')
def tests_config():
    sgains_data = os.environ.get('SGAINS_DATA')
    assert sgains_data is not None
    assert os.path.exists(sgains_data)
    assert os.path.isdir(sgains_data)

    # work_dir = os.path.join(sgains_data, 'test_data/hg19')
    # pristine_dir = os.path.join(sgains_data, 'test_data/hg19_pristine')

    config = Config.parse("sgains/tests/data/scpipe_tests.yml")
    # config.genome.work_dir = work_dir
    # config.genome.data_dir = pristine_dir

    return config



@pytest.fixture(scope='session')
def tests_genome(tests_config):
    genome = Genome(tests_config)
    assert genome is not None
    assert genome.version.VERSION == 'hg19'
    return genome


@pytest.fixture(scope='session')
def hisat2(tests_config, tests_genome):
    assert tests_config.genome.genome_version == 'hg19'
    return Hisat2(tests_config, tests_genome.version)


@pytest.fixture(scope='session')
def bwa(tests_config, tests_genome):
    assert tests_config.genome.genome_version == 'hg19'
    return BWA(tests_config, tests_genome.version)


@pytest.fixture(scope='session')
def hg(tests_config):
    return Genome(tests_config)


@pytest.fixture(scope='session')
def bin_boundaries(tests_config):
    bins_boundaries_fixture = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "test_data/R100_B10k/hg19_R100_B10k_bins_boundaries.txt"
    )
    df = pd.read_csv(
        bins_boundaries_fixture, sep='\t')
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
        "pipelines/tests/data/CJA0918.varbin.lubo.txt"
    )
    df = pd.read_csv(fixture_filename, sep='\t')
    return df
