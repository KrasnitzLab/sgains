'''
Created on Aug 2, 2017

@author: lubo
'''
import pytest
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from commands.genomeindex import GenomeIndex
from config import Config


@pytest.fixture
def tests_config():
    config = Config.load("scpipe/tests/data/scpipe_tests.yml")
    return config


@pytest.fixture
def argparser():
    parser = ArgumentParser(
        description="test program description",
        formatter_class=RawDescriptionHelpFormatter)

    return parser


@pytest.fixture
def argsubparser(argparser):
    subparser = argparser.add_subparsers(
        title="subcommands"
    )
    return subparser


def test_genomeindex_simple(argparser, argsubparser, tests_config):
    assert argparser is not None

    options = GenomeIndex(tests_config, argparser, argsubparser)
    assert options is not None

    options.add_options()

    argv = ["genomeindex", "--work-dir", "proba"]

    args = argparser.parse_args(argv)
    args.func(args)

    assert tests_config.genome.work_dir == "proba"
