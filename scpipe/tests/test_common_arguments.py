'''
Created on Jun 28, 2017

@author: lubo
'''
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import pytest
from common_arguments import Parser
import os
from hg19 import HumanGenome19


@pytest.fixture
def parser():
    argparser = ArgumentParser(
        description="test program description",
        formatter_class=RawDescriptionHelpFormatter)
    result = Parser.from_argument_parser(argparser)
    return result


def test_parser_simple(parser):
    assert parser is not None

    argv = ["-c", "scpipe_tests.yml"]
    config = parser.parse_arguments(argv)

    assert config is not None
    assert os.path.basename(config.filename) == "scpipe_tests.yml"


def test_parser_default(parser):
    argv = []
    config = parser.parse_arguments(argv)

    assert config is not None
    assert config.filename is None
    assert config.dirname is not None

    assert config.genome.version == "hg19"


def test_parser_genome_dir_overwrite(parser):
    argv = ["-c", "scpipe_tests.yml", "--genome-dir", "/test"]
    config = parser.parse_arguments(argv)

    assert config.genome.work_dir == "/test"


def test_parser_chroms_default(parser):
    argv = ["-c", "scpipe_tests.yml"]
    config = parser.parse_arguments(argv)

    assert set(config.chroms) == set(HumanGenome19.CHROMS)


def test_parser_chroms(parser):
    argv = ["-c", "scpipe_tests.yml", "-C", "chrM,chrY"]
    config = parser.parse_arguments(argv)

    assert set(config.chroms) == {'chrM', 'chrY'}


def test_default_reads_filenames(parser):
    argv = ["-c", "scpipe_tests.yml", "-C", "chrM,chrY"]
    config = parser.parse_arguments(argv)

    assert config.reads.mappable_regions == "mappable_regions.tsv"
    assert config.reads.mappable_positions_count == \
        "mappable_positions_count.yml"
    assert config.reads.chrom_sizes == "chrom_sizes.yml"


def test_default_bins_filenames(parser):
    argv = ["-c", "scpipe_tests.yml", "-C", "chrM,chrY"]
    config = parser.parse_arguments(argv)

    assert config.bins.bin_boundaries == "bin_boundaries.tsv"


def test_reads_work_dir(parser):
    argv = ["-c", "scpipe_tests.yml", "--reads-dir", "alaBalaPortocala"]
    config = parser.parse_arguments(argv)
    assert config.reads.work_dir == "alaBalaPortocala"


def test_reads_length(parser):
    argv = ["-c", "scpipe_tests.yml", "-l", "42"]
    config = parser.parse_arguments(argv)
    assert config.reads.length == 42


def test_bins_work_dir(parser):
    argv = ["-c", "scpipe_tests.yml", "--bins-dir", "alaBalaPortocala"]
    config = parser.parse_arguments(argv)
    assert config.bins.work_dir == "alaBalaPortocala"


def test_bins_length(parser):
    argv = ["-c", "scpipe_tests.yml", "-b", "10042"]
    config = parser.parse_arguments(argv)
    assert config.bins.bins_count == 10042


def test_output_abs(parser):
    argv = ["-c", "scpipe_tests.yml", "-o", "/out.txt"]
    config = parser.parse_arguments(argv)
    assert config.output == "/out.txt"


def test_output_stdout(parser):
    argv = ["-c", "scpipe_tests.yml", "-o", "-"]
    config = parser.parse_arguments(argv)
    assert config.output == "-"


def test_output_rel(parser):
    argv = ["-c", "scpipe_tests.yml", "-o", "out.txt"]
    config = parser.parse_arguments(argv)
    assert os.path.basename(config.output) == "out.txt"


def test_output_none(parser):
    argv = ["-c", "scpipe_tests.yml"]
    config = parser.parse_arguments(argv)
    assert config.output is None
