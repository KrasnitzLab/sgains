'''
Created on Aug 2, 2017

@author: lubo
'''
import pytest
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from config import Config
from cli_commands import parser_common_options
from commands.mappable_regions_command import MappableRegionsCommand


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
def mappable_regions_command(tests_config, argparser, argsubparser):
    command = MappableRegionsCommand(tests_config, argparser, argsubparser)
    assert command is not None

    return command


def test_mappable_regions_long(
        argparser, tests_config, mappable_regions_command):
    mappable_regions_command.add_options()

    argv = [
        "--dry-run", "--force",
        "mappable-regions",
        "--mappable-dir", "data/proba",
        "--genome-index", "probaindex",
        "--genome-dir", "data/hg19_safe/",
        "--read-length", "200",
        "--bowtie-opts", "-1 -2 -3",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    assert tests_config.force
    assert tests_config.dry_run

    assert tests_config.genome.work_dir == "data/hg19_safe/"
    assert tests_config.genome.index == "probaindex"

    assert tests_config.mappable_regions.length == 200
    assert tests_config.mappable_regions.work_dir == "data/proba"
    assert tests_config.mappable_regions.bowtie_opts == "-1 -2 -3"


def test_mappable_regions_short(
        argparser, tests_config, mappable_regions_command):
    mappable_regions_command.add_options()

    argv = [
        "-n", "-F",
        "mappable-regions",
        "-m", "data/proba",
        "-G", "probaindex",
        "--genome-dir", "data/hg19_safe/",
        "-l", "200",
        "--bowtie-opts", "-1 -2 -3",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    assert tests_config.force
    assert tests_config.dry_run

    assert tests_config.genome.work_dir == "data/hg19_safe/"
    assert tests_config.genome.index == "probaindex"

    assert tests_config.mappable_regions.length == 200
    assert tests_config.mappable_regions.work_dir == "data/proba"
    assert tests_config.mappable_regions.bowtie_opts == "-1 -2 -3"
