'''
Created on Jul 13, 2017

@author: lubo
'''
import os
import glob
import subprocess

from dask import distributed

from sgains.configuration.parser import Config
from sgains.genome import Genome

from termcolor import colored
import functools


class MappingPipeline(object):

    def __init__(self, config):
        self.config = config
        self.genome = Genome(config)

    @staticmethod
    def cells(filenames):
        return [
            Config.cellname(filename) for filename in filenames
        ]

    def reads_filenames(self):
        pattern = os.path.join(
            self.config.reads.reads_dir,
            "*{}".format(self.config.reads.reads_suffix)
        )
        filenames = glob.glob(pattern)
        return filenames


    @staticmethod
    def execute_once(dry_run, pipeline):
        mapping_command, index_command = pipeline
        print(colored(mapping_command, "green"))
        if not dry_run:
            res = subprocess.check_call(
                mapping_command,
                # stdout=subprocess.DEVNULL,
                # stderr=subprocess.DEVNULL,
                shell=True)
            assert res == 0
        print(colored(index_command, "green"))
        if not dry_run:
            res = subprocess.check_call(
                index_command,
                # stdout=subprocess.DEVNULL,
                # stderr=subprocess.DEVNULL,
                shell=True)
            assert res == 0

    def run(self, dask_client):
        fastq_filenames = self.reads_filenames()
        assert fastq_filenames

        work_dirname = self.config.mapping.mapping_dir
        os.makedirs(work_dirname, exist_ok=True)

        self.config.check_nonempty_workdir(
            work_dirname
        )

        commands = []
        mapping_opts = self.config.mapping.mapping_aligner_options.split(' ')

        for fastq_filename in fastq_filenames:
            mapping_pipeline = self.genome.aligner\
                .build_mapping_pipeline(
                    fastq_filename,
                    num_lines=0, options=mapping_opts)
            mapping_command = ' | '.join(mapping_pipeline)
            index_command = self.genome.aligner\
                .build_samtools_index_command(fastq_filename)
            index_command = " ".join(index_command)

            commands.append((mapping_command, index_command))

        for command in commands:
            print(colored(
                f"going to execute {command}",
                "green"))

        if self.config.dry_run:
            return

        assert dask_client

        delayed_tasks = dask_client.map(
                functools.partial(
                    MappingPipeline.execute_once, self.config.dry_run),
                commands)

        distributed.wait(delayed_tasks)

        # for fut in delayed_tasks:
        #     print("fut done:", fut.done())
        #     print("fut exception:", fut.exception())
        #     print("fut traceback:", fut.traceback())
        #     if fut.traceback() is not None:
        #         traceback.print_tb(fut.traceback())
        #     print(fut.result())
