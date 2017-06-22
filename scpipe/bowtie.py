'''
Created on Jun 22, 2017

@author: lubo
'''
import asyncio
import functools
import logging
import sys
from hg19 import HumanGenome19, MappableRegion, MappableState
from config import Config


logging.basicConfig(
    level=logging.DEBUG,
    format='%(levelname)7s: %(message)s',
    stream=sys.stderr,
)
LOG = logging.getLogger('')


class BowtieProtocol(asyncio.SubprocessProtocol):
    FD_NAMES = ['stdin', 'stdout', 'stderr']

    OK = 0

    def __init__(self, done_future):
        self.done = done_future
        self.buffer = bytearray()
        self.transport = None
        super(BowtieProtocol, self).__init__()

    def connection_made(self, transport):
        print('process started {}'.format(transport.get_pid()))
        self.transport = transport

    def pipe_data_received(self, fd, data):
        print('read {} bytes from {}'.format(
            len(data),
            self.FD_NAMES[fd]))

    def process_exited(self):
        print('process exited')
        return_code = self.transport.get_returncode()
        print('return code {}'.format(return_code))
        if return_code == self.OK:
            print('return code is OK')
            results = ['OK']
        else:
            results = []
        self.done.set_result((return_code, results))


async def run_bowtie(loop):
    print('in run_bowtie...')
    cmd_done = asyncio.Future(loop=loop)

    factory = functools.partial(BowtieProtocol, cmd_done)
    proc = loop.subprocess_exec(
        factory,
        'bowtie', '-S', '-t', '-v', '0', '-m', '1', '-f', 'data/hg19/genomeindex', '-',
        stdin=asyncio.subprocess.PIPE,
        stderr=None,
    )
    try:
        print('launching process')
        transport, protocol = await proc
        print('waiting for process to complete')
        print(dir(transport))
        print(dir(protocol))

        loop.connect_write_pipe()

        await cmd_done
    finally:
        transport.close()

    return cmd_done.result()


async def start_bowtie():
    create = asyncio.create_subprocess_exec(
        'bowtie', '-S', '-t', '-v', '0', '-m', '1',
        '-f', 'data/hg19/genomeindex', '-',
        stdin=asyncio.subprocess.PIPE,
        stdout=asyncio.subprocess.PIPE,
    )
    print('launching process')
    proc = await create
    print('pid {}'.format(proc.pid))
    return proc


def close_bowtie(bowtie):
    print("close_bowtie called...")
    bowtie.stdin.close()


async def bowtie_write_reads(bowtie, reads_generator):
    for rec in reads_generator:
        await HumanGenome19.async_write_fasta(bowtie.stdin, rec)

    print("close_bowtie called...")
    bowtie.stdin.close()
    return "OK"


def parse_mapping(line):
    row = line.rstrip().split("\t")
    name = row[0]
    flag = int(row[1])
    chrom = row[2]
    start = int(row[3])
    end = start + 1
    return MappableRegion(flag=flag, chrom=chrom, start=start)


async def bowtie_mappings_reader(bowtie):
    prev = None
    state = MappableState.OUT

    while True:
        line = await bowtie.stdout.readline()
        if not line:
            break
        line = line.decode()
        if line[0] == '@':
            # comment
            continue
        mapping = parse_mapping(line)

        if state == MappableState.OUT:
            if mapping.flag == 0:
                prev = mapping
                state = MappableState.IN
        else:
            if mapping.flag == 0:
                if mapping.chrom == prev.chrom:
                    prev.extend(mapping.start)
                else:
                    print(prev)
                    prev = mapping
            else:
                print(prev)
                state = MappableState.OUT

    if state == MappableState.IN:
        print(prev)


async def bowtie_simple_reader(bowtie):
    print("bowtie_simple_reader started...")
    while True:
        line = await bowtie.stdout.readline()
        if not line:
            break
        line = line.decode()
        if line[0] == '@':
            print("comment")
            continue
        row = line.split('\t')
        if row[4] == '0':
            print("ROW:", row)

#         if row[4] == '0':
#             print("OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOWWWWWWWWWWWWWWWWW!!")
    return


async def bowtie_experiments(loop):
    config = Config.load("scpipe10k100.yml")
    hg = HumanGenome19(config)
    bowtie = await start_bowtie()

    reads_generator = hg.generate_reads(['chrM'], 100)

    reader = asyncio.Task(
        bowtie_mappings_reader(bowtie))

    writer = asyncio.Task(
        bowtie_write_reads(bowtie, reads_generator))

    await bowtie.wait()

    await writer
    await reader

    return bowtie


if __name__ == '__main__':

    event_loop = asyncio.get_event_loop()
    # Enable debugging
    # event_loop.set_debug(True)

    try:
        event_loop.run_until_complete(bowtie_experiments(event_loop))
    finally:
        event_loop.close()
