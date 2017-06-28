'''
Created on Jun 28, 2017

@author: lubo
'''
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import pytest
from common_arguments import Arguments
import os


@pytest.fixture
def parser():
    argparser = ArgumentParser(
        description="test program description",
        formatter_class=RawDescriptionHelpFormatter)
    result = Arguments.scpipe_arguments(argparser)
    return result


def test_parser_simple(parser):
    assert parser is not None

    argv = ["-c", "scpipe_tests.yml"]
    config = parser.process_scpipe_agrments(argv)

    assert config is not None
    assert os.path.basename(config.filename) == "scpipe_tests.yml"
