'''
Created on Jun 22, 2017

@author: lubo
'''
import asyncio
import pytest
import logging
import sys

logging.basicConfig(
    level=logging.DEBUG,
    format='%(levelname)7s: %(message)s',
    stream=sys.stderr,
)
LOG = logging.getLogger('')

pytestmark = pytest.mark.asyncio


async def test_async_mappings_generator(hg, event_loop):
    # Enable debugging
    # event_loop.set_debug(True)

    bowtie = await hg.async_start_bowtie()

    reads_generator = hg.generate_reads(['chrM'], 100)
    writer = asyncio.Task(
        hg.async_write_reads_generator(bowtie.stdin, reads_generator)
    )

    async for mapping in hg.async_mappable_regions_generator(bowtie.stdout):
        print(mapping)

    await bowtie.wait()
    await writer


async def test_async_generate_mappable_regions(hg, event_loop):
    # Enable debugging
    event_loop.set_debug(True)
    await hg.async_generate_mappable_regions(['chrM'], 100)
