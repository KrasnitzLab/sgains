'''
Created on Jun 27, 2017

@author: lubo
'''
from hg19 import HumanGenome19
from config import Config
import asyncio
import logging
import sys


logging.basicConfig(
    level=logging.DEBUG,
    format='%(levelname)7s: %(message)s',
    stream=sys.stderr,
)
LOG = logging.getLogger('')


async def main():
    config = Config.load('scpipe_tests.yml')
    hg = HumanGenome19(config)
    assert hg is not None

    reads_generator = hg.generate_reads(['chrM'], 100)
    bowtie = await hg.async_start_bowtie()

    writer = asyncio.Task(
        hg.async_write_reads_generator(bowtie.stdin, reads_generator)
    )

    while True:
        line = await bowtie.stdout.readline()
        sys.stdout.write(line.decode())
        if not line:
            break

    await bowtie.wait()
    await writer


if __name__ == "__main__":
    event_loop = asyncio.get_event_loop()
    # Enable debugging
    event_loop.set_debug(True)

    try:
        event_loop.run_until_complete(main())
    finally:
        event_loop.close()
