'''
Created on Aug 2, 2017

@author: lubo
'''
import pytest
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from commands.genomeindex import GenomeIndexCommand
from config import Config
from cli_commands import parser_common_options


@pytest.fixture
def tests_config():
    config = Config.load("scpipe/tests/data/scpipe_tests.yml")
    return config


@pytest.fixture
def argparser():
    parser = ArgumentParser(
        description="test program description",
        formatter_class=RawDescriptionHelpFormatter)
    parser_common_options(parser)
    return parser


@pytest.fixture
def argsubparser(argparser):
    subparser = argparser.add_subparsers(
        title="subcommands"
    )
    return subparser


@pytest.fixture
def genomeindex_command(tests_config, argparser, argsubparser):
    command = GenomeIndexCommand(tests_config, argparser, argsubparser)
    assert command is not None

    return command


def test_genomeindex_long(argparser, tests_config, genomeindex_command):
    genomeindex_command.add_options()

    argv = [
        "--dry-run", "--force",
        "genomeindex",
        "--work-dir", "proba",
        "--genome-index", "probaindex",
        "--data-dir", "data/hg19_safe/",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    assert tests_config.force
    assert tests_config.dry_run

    assert tests_config.genome.work_dir == "proba"
    assert tests_config.genome.index == "probaindex"
    assert tests_config.genome.data_dir == "data/hg19_safe/"


def test_genomeindex_short(argparser, tests_config, genomeindex_command):
    genomeindex_command.add_options()

    argv = [
        "-n", "-F",
        "genomeindex",
        "-o", "proba",
        "-G", "probaindex",
        "-i", "data/hg19_safe/",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    assert tests_config.force
    assert tests_config.dry_run

    assert tests_config.genome.work_dir == "proba"
    assert tests_config.genome.index == "probaindex"
    assert tests_config.genome.data_dir == "data/hg19_safe/"
