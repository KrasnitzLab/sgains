'''
Created on Aug 3, 2017

@author: lubo
'''
import pytest
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from config import Config
from commands.mappable_regions_command import MappableRegionsCommand
from commands.genomeindex_command import GenomeIndexCommand
from commands.bins_command import BinsCommand
from commands.prepare_command import PrepareCommand
from commands.mapping_command import MappingCommand
from commands.varbin_command import VarbinCommand
from commands.sc_clust_command import SCclustCommand
from commands.process_command import ProcessCommand
from commands.common import OptionsBase


@pytest.fixture
def tests_config():
    config = Config.load("tests/data/scpipe_tests.yml")
    return config


@pytest.fixture
def argparser():
    parser = ArgumentParser(
        description="test program description",
        formatter_class=RawDescriptionHelpFormatter)
    OptionsBase.common_options(parser)
    return parser


@pytest.fixture
def argsubparser(argparser):
    subparser = argparser.add_subparsers(
        title="subcommands"
    )
    return subparser


@pytest.fixture
def mappable_regions_command(argparser, argsubparser):
    command = MappableRegionsCommand(argparser, argsubparser)
    assert command is not None

    return command


@pytest.fixture
def genomeindex_command(argparser, argsubparser):
    command = GenomeIndexCommand(argparser, argsubparser)
    assert command is not None

    return command


@pytest.fixture
def bins_command(argparser, argsubparser):
    command = BinsCommand(argparser, argsubparser)
    assert command is not None

    return command


@pytest.fixture
def prepare_command(argparser, argsubparser):
    command = PrepareCommand(argparser, argsubparser)
    assert command is not None

    return command


@pytest.fixture
def mapping_command(argparser, argsubparser):
    command = MappingCommand(argparser, argsubparser)
    assert command is not None

    return command


@pytest.fixture
def varbin_command(argparser, argsubparser):
    command = VarbinCommand(argparser, argsubparser)
    assert command is not None

    return command


@pytest.fixture
def sc_clust_command(argparser, argsubparser):
    command = SCclustCommand(argparser, argsubparser)
    assert command is not None

    return command


@pytest.fixture
def process_command(argparser, argsubparser):
    command = ProcessCommand(argparser, argsubparser)
    assert command is not None

    return command
