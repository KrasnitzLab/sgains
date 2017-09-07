'''
Created on Jun 22, 2017

@author: lubo
'''
import asyncio
import pytest
import logging
import sys
from pipelines.mappableregions_pipeline import MappableRegionsPipeline

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
