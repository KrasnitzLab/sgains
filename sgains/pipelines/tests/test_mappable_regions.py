import pytest
import os
import io
import numpy as np
import pandas as pd

from sgains.pipelines.mappableregions_pipeline import MappableRegionsPipeline


def test_test_config(tests_config):
    print(tests_config.genome.genome_dir)


def test_bowtie_command(tests_config):
    pipeline = MappableRegionsPipeline(tests_config)
    command = pipeline.genome.aligner.build_mappable_regions_command()
    print(command)
    assert command[0] == 'hisat2'


def test_reads_generator(tests_config):

    pipeline = MappableRegionsPipeline(tests_config)
    reads_generator = pipeline.generate_reads(['chrM'], 100)
    # ouptut_writer = pipeline.mappable_regions_writer()
    count = 0
    for _read in reads_generator:
        # print(read)
        count += 1
    print(count)
    assert 16472 == count


def test_generate_mappable_regions(tests_config):
    pipeline = MappableRegionsPipeline(tests_config)
    pipeline.generate_mappable_regions(['chrM'], 100)


@pytest.mark.parametrize("chrom", ['chrM'])
def test_async_mappable_regions_50(tests_config, chrom):

    filename = os.path.join(
        'sgains/tests/data',
        "{}.50mer.mappable.regions.txt.gz".format(chrom)
    )
    gold_df = pd.read_csv(
        filename, sep='\t', compression='gzip',
        header=None,
        names=['chrom', 'start', 'end'])
    with io.StringIO() as outfile:
        pipeline = MappableRegionsPipeline(tests_config)
        pipeline.generate_mappable_regions(
            [chrom], 50, outfile=outfile)
        outfile.flush()

        infile = io.StringIO(outfile.getvalue())
        df = pd.read_csv(
            infile, sep='\t', header=None,
            names=['chrom', 'start', 'end'])
        print(df.head())
    assert np.all(df.columns == gold_df.columns)
    assert len(df) == len(gold_df)

    assert np.all(df.start == gold_df.start)
    assert np.all(df.end == gold_df.end)


def test_generate_mappable_regions_hisat(tests_config, hisat2):

    assert hisat2 is not None
    pipeline = MappableRegionsPipeline(tests_config, hisat2)

    pipeline.generate_mappable_regions(['chrM'], 100)


def test_generate_mappable_regions_bwa(tests_config, bwa):

    assert bwa is not None
    pipeline = MappableRegionsPipeline(tests_config, bwa)

    pipeline.generate_mappable_regions(['chrM'], 100)

@pytest.mark.parametrize("chrom", ['chrM'])
def test_async_mappable_regions_50_hisat(tests_config, chrom, hisat2):

    filename = os.path.join(
        'sgains/tests/data',
        "{}.50mer.mappable.regions.txt.gz".format(chrom)
    )
    gold_df = pd.read_csv(
        filename, sep='\t', compression='gzip',
        header=None,
        names=['chrom', 'start', 'end'])
    with io.StringIO() as outfile:
        pipeline = MappableRegionsPipeline(tests_config, hisat2)
        pipeline.generate_mappable_regions(
            [chrom], 50, outfile=outfile)
        outfile.flush()

        infile = io.StringIO(outfile.getvalue())
        df = pd.read_csv(
            infile, sep='\t', header=None,
            names=['chrom', 'start', 'end'])
        print(df.head())
    assert np.all(df.columns == gold_df.columns)
    assert len(df) == len(gold_df)

    assert np.all(df.start == gold_df.start)
    assert np.all(df.end == gold_df.end)
