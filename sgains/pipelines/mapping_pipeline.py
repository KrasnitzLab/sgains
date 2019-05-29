'''
Created on Jul 13, 2017

@author: lubo
'''
import os
import subprocess
from sgains.config import Config
from termcolor import colored
import multiprocessing
import functools


class MappingPipeline(object):

    def __init__(self, config):
        self.config = config

    @staticmethod
    def cells(filenames):
        return [
            Config.cellname(filename) for filename in filenames
        ]

    @staticmethod
    def unarchive_stage(filename):
        _base, ext = os.path.splitext(filename)
        if ext == '.gz':
            return ['gunzip', '-c', filename]
        elif ext == '.fastq':
            return ['cat', filename]
        assert False, 'unexptected archive extention: {}, {}'.format(
            ext, filename)

    @staticmethod
    def head_stage(_filename, lines=1000000):
        return [
            'head',
            '-n',
            '{}'.format(lines),
            '-'
        ]

    def bowtie_stage(self, filename):
        cellname = Config.cellname((filename))
        reportfile = os.path.join(
            self.config.mapping_dirname(),
            "{}.bowtie_report.log".format(cellname)
        )
        reportfile = self.config.abspath(reportfile)

        bowtie_opts = self.config.mapping.mapping_bowtie_opts.split(' ')
        return [
            'bowtie',
            *bowtie_opts,
            self.config.genome_index_filename(),
            '-',
            '2>',
            reportfile,
        ]

    @staticmethod
    def samtools_view_stage(_filename):
        return [
            'samtools',
            'view',
            '-bu',
            '-o',
            '-',
            '-'
        ]

    @staticmethod
    def samtools_sort_stage(_filename):
        return [
            'samtools',
            'sort',
            '-o',
            '-',
            '-'
        ]

    @staticmethod
    def samtools_rmdup_stage(filename):
        return [
            'samtools',
            'rmdup',
            '-s',
            '-',
            '-'
        ]

    def samtools_view_store_stage(self, filename):
        cellname = Config.cellname((filename))
        outfile = os.path.join(
            self.config.mapping_dirname(),
            "{}.rmdup.bam".format(cellname)
        )
        outfile = self.config.abspath(outfile)
        return [
            'samtools',
            'view',
            '-b',
            '-o',
            outfile,
            '-',
        ]

    @staticmethod
    def samtools_view_remove_unmappted_stage(_filename):
        return [
            'samtools',
            'view',
            '-bu',
            '-F',
            '4',
            '-o',
            '-',
            '-'
        ]

    def samtools_index_bam(self, filename):
        cellname = Config.cellname((filename))
        outfile = os.path.join(
            self.config.mapping_work_dirname(),
            "{}.rmdup.bam".format(cellname)
        )
        outfile = self.config.abspath(outfile)
        return [
            'samtools'
            'index',
            outfile,
        ]

    @staticmethod
    def execute_once(dry_run, pipeline):
        print(colored(' '.join(pipeline), "green"))
        if not dry_run:
            subprocess.check_call(
                ' '.join(pipeline),
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                shell=True)

    def run(self, dask_client):
        fastq_filenames = self.config.mapping_reads_filenames()
        assert fastq_filenames

        work_dirname = self.config.abspath(self.config.mapping.mapping_dir)
        self.config.check_nonempty_workdir(
            work_dirname
        )

        commands = []
        for filename in fastq_filenames:
            pipeline = [
                *self.unarchive_stage(filename),
                '|',
                #                 *self.head_stage(filename, lines=40000),
                #                 '|',
                *self.bowtie_stage(filename),
                '|',
                *self.samtools_view_stage(filename),
                '|',
                *self.samtools_view_remove_unmappted_stage(filename),
                '|',
                *self.samtools_sort_stage(filename),
                '|',
                *self.samtools_rmdup_stage(filename),
                '|',
                *self.samtools_view_store_stage(filename),
            ]
            commands.append(pipeline)

        # pool = multiprocessing.Pool(processes=self.config.parallel)
        # pool.map(functools.partial(
        #     MappingPipeline.execute_once, self.config.dry_run), commands)
        assert dask_client

        delayed_tasks = [
            dask_client.submit(
                functools.partial(
                    MappingPipeline.execute_once, self.config.dry_run),
                cmd) for cmd in commands
        ]
        fut = dask_client.compute(delayed_tasks)
        dask_client.gather(fut)
