'''
Created on Jul 13, 2017

@author: lubo
'''
import os
import subprocess

from dask import distributed

from sgains.config import Config
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

    # @staticmethod
    # def unarchive_stage(filename):
    #     _base, ext = os.path.splitext(filename)
    #     if ext == '.gz':
    #         return ['gunzip', '-c', filename]
    #     elif ext == '.fastq':
    #         return ['cat', filename]
    #     assert False, 'unexptected archive extention: {}, {}'.format(
    #         ext, filename)

    # @staticmethod
    # def head_stage(_filename, lines=1000000):
    #     return [
    #         'head',
    #         '-n',
    #         '{}'.format(lines),
    #         '-'
    #     ]

    # def mapping_stage(self, filename):
    #     cellname = Config.cellname(filename)
    #     reportfile = os.path.join(
    #         self.config.mapping_dirname(),
    #         "{}.aligner_report.log".format(cellname)
    #     )
    #     reportfile = self.config.abspath(reportfile)
    #     mapping_opts = self.config.mapping.mapping_opts.split(' ')
    #     command = self.genome.aligner.build_mapping_command(mapping_opts)
    #     command.extend([
    #         '2>',
    #         reportfile,
    #     ])
    #     filter_command = self.genome.aligner.build_postmapping_filter()
    #     if filter_command:
    #         command.extend(['|', *filter_command])

    #     return command

    # @staticmethod
    # def samtools_view_stage(_filename):
    #     return [
    #         'samtools',
    #         'view',
    #         '-bu',
    #         '-o',
    #         '-',
    #         '-'
    #     ]

    # @staticmethod
    # def samtools_sort_stage(_filename):
    #     return [
    #         'samtools',
    #         'sort',
    #         '-o',
    #         '-',
    #         '-'
    #     ]

    # @staticmethod
    # def samtools_rmdup_stage(_filename):
    #     return [
    #         'samtools',
    #         'rmdup',
    #         '-s',
    #         '-',
    #         '-'
    #     ]

    # def samtools_view_store_stage(self, filename):
    #     cellname = Config.cellname((filename))
    #     outfile = os.path.join(
    #         self.config.mapping_dirname(),
    #         "{}.rmdup.bam".format(cellname)
    #     )
    #     outfile = self.config.abspath(outfile)
    #     return [
    #         'samtools',
    #         'view',
    #         '-b',
    #         '-o',
    #         outfile,
    #         '-',
    #     ]

    # @staticmethod
    # def samtools_view_remove_unmappted_stage(_filename):
    #     return [
    #         'samtools',
    #         'view',
    #         '-bu',
    #         '-F',
    #         '4',
    #         '-o',
    #         '-',
    #         '-'
    #     ]

    # def samtools_index_bam(self, filename):
    #     cellname = Config.cellname((filename))
    #     outfile = os.path.join(
    #         self.config.mapping_work_dirname(),
    #         "{}.rmdup.bam".format(cellname)
    #     )
    #     outfile = self.config.abspath(outfile)
    #     return [
    #         'samtools'
    #         'index',
    #         outfile,
    #     ]

    @staticmethod
    def execute_once(dry_run, pipeline):
        print(os.environ)
        mapping_command, index_command = pipeline
        print(colored(mapping_command, "green"))
        if not dry_run:
            res = subprocess.check_call(
                mapping_command,
                # stdout=subprocess.DEVNULL,
                # stderr=subprocess.DEVNULL,
                shell=True)
            assert res == 0
            # print(res)
            # print(res.cmd)
            # print(res.stdout)
            # print(res.stderr)
        print(colored(index_command, "green"))
        if not dry_run:
            res = subprocess.check_call(
                index_command,
                # stdout=subprocess.DEVNULL,
                # stderr=subprocess.DEVNULL,
                shell=True)
            assert res == 0
            # print(res)

    def run(self, dask_client):
        fastq_filenames = self.config.mapping_reads_filenames()
        assert fastq_filenames

        work_dirname = self.config.abspath(self.config.mapping.mapping_dir)
        self.config.check_nonempty_workdir(
            work_dirname
        )

        commands = []
        mapping_opts = self.config.mapping.mapping_opts.split(' ')

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
