import os
import pytest

from sgains.config import Config
from sgains.aligners import Hisat2
from sgains.genome import Genome


@pytest.fixture
def tests_config():
    sgains_data = os.environ.get('SGAINS_DATA')
    assert sgains_data is not None
    assert os.path.exists(sgains_data)
    assert os.path.isdir(sgains_data)

    work_dir = os.path.join(sgains_data, 'test_data/hg19')
    pristine_dir = os.path.join(sgains_data, 'test_data/hg19_pristine')

    config = Config.load("sgains/tests/data/scpipe_tests.yml")
    config.genome.work_dir = work_dir
    config.genome.data_dir = pristine_dir

    return config


@pytest.fixture
def tests_genome(tests_config):
    genome = Genome(tests_config)
    assert genome is not None
    assert genome.version.VERSION == 'hg19'
    return genome


@pytest.fixture
def hisat2(tests_config, tests_genome):
    assert tests_config.genome.version == 'hg19'
    return Hisat2(tests_config, tests_genome.version)
