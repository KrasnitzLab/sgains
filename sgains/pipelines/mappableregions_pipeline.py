'''
Created on Jul 31, 2017

@author: lubo
'''
import traceback
import distributed
import queue

from subprocess import Popen, PIPE
from threading import Thread

from sgains.genome import Genome
from termcolor import colored
import os
from Bio.SeqRecord import SeqRecord  # @UnresolvedImport

from sgains.utils import MappableState, Mapping, MappableRegion, LOG
import sys
import shutil


class BowtieInputGenerateThread(Thread):

    def __init__(self, control_queue, bowtie_input, input_function_generator):
        super(BowtieInputGenerateThread, self).__init__()

        self.control_queue = control_queue
        self.bowtie_input = bowtie_input
        self.input_function_generator = input_function_generator

    def run(self):
        for rec in self.input_function_generator:
            out = Genome.to_fasta_string(rec)
            self.bowtie_input.write(out)
            self.bowtie_input.flush()
        self.bowtie_input.close()


class BowtieOutputProcessThread(Thread):

    def __init__(self, control_queue, bowtie_output, output_process_function):
        super(BowtieOutputProcessThread, self).__init__()

        self.control_queue = control_queue
        self.bowtie_output = bowtie_output
        self.output_process_function = output_process_function

    def run(self):
        prev = None
        state = MappableState.OUT

        while True:
            line = self.bowtie_output.readline()
            if not line:
                break

            line = line.decode()
            if line[0] == '@':
                # comment
                continue

            mapping = Mapping.parse_sam(line)

            if state == MappableState.OUT:
                if mapping.flag == 0:
                    prev = MappableRegion(mapping)
                    state = MappableState.IN
            else:
                if mapping.flag == 0:
                    if mapping.chrom == prev.chrom:
                        prev.extend(mapping.start)
                    else:
                        self.output_process_function(prev)
                        prev = MappableRegion(mapping)
                else:
                    self.output_process_function(prev)
                    state = MappableState.OUT

        if state == MappableState.IN:
            self.output_process_function(prev)

        self.bowtie_output.close()
        self.control_queue.put("out_done")


class MappableRegionsPipeline(object):

    def __init__(self, config):
        self.config = config
        self.hg = Genome(self.config)

    def mappable_regions_check(self, chroms, mappable_regions_df):

        for chrom in chroms:
            chrom_df = mappable_regions_df[mappable_regions_df.chrom == chrom]
            chrom_df = chrom_df.sort_values(
                by=['chrom', 'start_pos', 'end_pos'])
            start_pos_count = len(chrom_df.start_pos.unique())
            if start_pos_count < len(chrom_df):
                LOG.error(
                    "chrom {} has duplicate mappable regions".format(chrom))

    def generate_reads(self, chroms, read_length):
        try:
            for chrom in chroms:
                seq_record = self.hg.load_chrom(chrom)
                for i in range(len(seq_record) - read_length + 1):
                    out_record = SeqRecord(
                        seq_record.seq[i: i + read_length],
                        id="{}.{}".format(chrom, i + 1),
                        description="generated_read"
                    )
                    yield out_record
        finally:
            pass

    def bowtie_command(self, bowtie_opts=""):
        genomeindex = self.config.genome_index_filename()
        if bowtie_opts:
            command = [
                'bowtie', '-S', '-t', '-v', '0', '-m', '1',
                *bowtie_opts.split(' '),
                '-f', genomeindex, '-',
            ]
        else:
            command = [
                'bowtie', '-S', '-t', '-v', '0', '-m', '1',
                '-f', genomeindex, '-',
            ]
        print(colored(
            "going to execute bowtie: {}".format(" ".join(command)),
            "green"
        ))
        return command

    def generate_mappable_regions(
            self, chroms, read_length,
            outfile=None, bowtie_opts=""):

        if outfile is None:
            outfile = sys.stdout

        bowtie_command = self.bowtie_command(bowtie_opts=bowtie_opts)
        reads_generator = self.generate_reads(chroms, read_length)

        def bowtie_output_process_function(line):
            outfile.write(str(line))
            outfile.write("\n")

        with Popen(bowtie_command, stdout=PIPE, stdin=PIPE) as proc:

            control_queue = queue.Queue()
            input_thread = BowtieInputGenerateThread(
                control_queue,
                proc.stdin,
                reads_generator)
            output_thread = BowtieOutputProcessThread(
                control_queue,
                proc.stdout,
                bowtie_output_process_function)

            input_thread.start()
            output_thread.start()

            while True:
                msg = None
                try:
                    msg = control_queue.get()
                except queue.Empty:
                    print("timeout - queue empty")
                    msg = None
                if msg == 'out_done':
                    print("output done")
                    break
                if msg == 'in_done':
                    print('input done')
            input_thread.join()
            output_thread.join()

    def run_once(self, chrom):

        outfilename = self.config.mappable_regions_filename(chrom)
        with open(outfilename, "w") as outfile:
            self.generate_mappable_regions(
                [chrom], read_length=50,
                outfile=outfile)

    def concatenate_all_chroms(self):
        dst = self.config.mappable_regions_filename()
        if os.path.exists(dst) and not self.config.force:
            print(colored(
                "destination mappable regions file already exists"
                "use --force to overwrite", "red"))
            raise ValueError("destination file exists... use --force")

        if not self.config.dry_run:
            with open(dst, 'wb') as output:
                for chrom in self.hg.version.CHROMS:
                    src = self.config.mappable_regions_filename(chrom)
                    print(colored(
                        "appending {} to {}".format(src, dst),
                        "green"))
                    with open(src, 'rb') as src:
                        if not self.config.dry_run:
                            shutil.copyfileobj(src, output, 1024 * 1024 * 10)

    def run(self, dask_client):
        outfilename = self.config.mappable_regions_filename()
        print(colored(
            "going to generate mappable regions with length {} "
            "from genome {} into {}".format(
                self.config.mappable_regions.length,
                self.config.genome.work_dir,
                outfilename
            ),
            "green"))

        if os.path.exists(outfilename) and not self.config.force:
            print(colored(
                "output file {} already exists; "
                "use --force to overwrite".format(outfilename),
                "red"))
            raise ValueError("output file already exists")

        if not self.config.genome_index_filename_exists():
            print(colored(
                "genome index file {} not found".format(
                    self.config.genome_index_filename()),
                "red"))
            raise ValueError("genome index file not found")

        if self.config.dry_run:
            return

        if not os.path.exists(self.config.mappable_regions.work_dir):
            os.makedirs(self.config.mappable_regions.work_dir)

        assert dask_client

        delayed_tasks = dask_client.map(
                self.run_once, self.hg.version.CHROMS[:1])
        print(delayed_tasks)
        print(dask_client.scheduler_info())

        distributed.wait(delayed_tasks)

        for fut in delayed_tasks:
            print("fut done:", fut.done())
            print("fut exception:", fut.exception())
            print("fut traceback:", fut.traceback())
            if fut.traceback() is not None:
                traceback.print_tb(fut.traceback())
            if fut.exception() is None:
                print(fut.result())

        self.concatenate_all_chroms()
