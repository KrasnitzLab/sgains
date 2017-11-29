'''
Created on Jun 22, 2017

@author: lubo
'''
import asyncio
import pytest
import logging
import sys
from pipelines.mappableregions_pipeline import MappableRegionsPipeline
from utils import Mapping
import io
import pandas as pd
import os
import numpy as np


logging.basicConfig(
    level=logging.DEBUG,
    format='%(levelname)7s: %(message)s',
    stream=sys.stderr,
)
LOG = logging.getLogger('')

pytestmark = pytest.mark.asyncio


async def test_async_mappings_generator(tests_config, event_loop):
    # Enable debugging
    event_loop.set_debug(True)
    pipeline = MappableRegionsPipeline(tests_config)

    bowtie = await pipeline.async_start_bowtie()

    reads_generator = pipeline.generate_reads(['chrM'], 100)
    writer = asyncio.Task(
        pipeline.async_write_reads_generator(bowtie.stdin, reads_generator)
    )

    async for mapping in pipeline.async_mappable_regions_generator(
            bowtie.stdout):
        print(mapping)

    await bowtie.wait()
    await writer


async def test_async_generate_mappable_regions(tests_config, event_loop):
    # Enable debugging
    event_loop.set_debug(True)

    pipeline = MappableRegionsPipeline(tests_config)

    await pipeline.async_generate_mappable_regions(['chrM'], 100)


async def test_generate_reads(tests_config):
    pipeline = MappableRegionsPipeline(tests_config)
    generator = pipeline.generate_reads(['chr1'], 100)

    for num, rec in enumerate(generator):
        print(rec.id, len(rec))
        if num >= 10:
            break
    generator.close()


# @pytest.mark.parametrize("chrom", ['chrM', 'chr21'])
@pytest.mark.parametrize("chrom", ['chrM'])
async def test_bowtie_mappings(tests_config, event_loop, chrom):
    event_loop.set_debug(True)
    pipeline = MappableRegionsPipeline(tests_config)

    bowtie = await pipeline.async_start_bowtie()

    reads_generator = pipeline.generate_reads([chrom], 100)
    writer = asyncio.Task(
        pipeline.async_write_reads_generator(bowtie.stdin, reads_generator)
    )

    infile = bowtie.stdout

    while True:
        line = await infile.readline()
        if not line:
            break
        line = line.decode()
        if line[0] == '@':
            # comment
            continue
        mapping = Mapping.parse_sam(line)
        if mapping.flag == 0:
            # print(mapping)
            chromosome, pos = mapping.name.split('.')
            assert int(pos) == mapping.start
            assert chrom == chromosome

    await writer


# @pytest.mark.parametrize("chrom", ['chrM', 'chr21'])
@pytest.mark.parametrize("chrom", ['chrM'])
async def test_async_mappable_regions_50(tests_config, event_loop, chrom):
    event_loop.set_debug(True)
    pipeline = MappableRegionsPipeline(tests_config)

    filename = os.path.join(
        'scpipe/tests/data',
        "{}.50mer.mappable.regions.txt.gz".format(chrom)
    )
    gold_df = pd.read_csv(
        filename, sep='\t', compression='gzip',
        header=None,
        names=['chrom', 'start', 'end'])
    with io.StringIO() as outfile:
        await pipeline.async_generate_mappable_regions(
            [chrom], 50, outfile=outfile)
        infile = io.StringIO(outfile.getvalue())
        df = pd.read_csv(
            infile, sep='\t', header=None,
            names=['chrom', 'start', 'end'])
        print(df.head())

    assert np.all(df.columns == gold_df.columns)
    assert np.all(df.start == gold_df.start)
    assert np.all(df.end == gold_df.end)
