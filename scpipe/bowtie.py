'''
Created on Jun 22, 2017

@author: lubo
'''
import asyncio
import functools


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


if __name__ == '__main__':

    event_loop = asyncio.get_event_loop()
    try:
        return_code, results = event_loop.run_until_complete(
            run_bowtie(event_loop)
        )
    finally:
        event_loop.close()

    if return_code:
        print('error exit {}'.format(return_code))
    else:
        print('\nbowtie exited OK:')
