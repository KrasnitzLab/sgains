'''
Created on Aug 3, 2017

@author: lubo
'''
import pytest
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from config import Config
from cli_commands import parser_common_options
from commands.mappable_regions_command import MappableRegionsCommand
from commands.genomeindex_command import GenomeIndexCommand
from commands.bins_command import BinsCommand
from commands.prepare_command import PrepareCommand


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


@pytest.fixture
def genomeindex_command(tests_config, argparser, argsubparser):
    command = GenomeIndexCommand(tests_config, argparser, argsubparser)
    assert command is not None

    return command


@pytest.fixture
def bins_command(tests_config, argparser, argsubparser):
    command = BinsCommand(tests_config, argparser, argsubparser)
    assert command is not None

    return command


@pytest.fixture
def prepare_command(tests_config, argparser, argsubparser):
    command = PrepareCommand(tests_config, argparser, argsubparser)
    assert command is not None

    return command
