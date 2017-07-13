'''
Created on Jul 13, 2017

@author: lubo
'''
import os
import subprocess


class MappingPipeline(object):

    def __init__(self, config):
        self.config = config

    @staticmethod
    def cells(filenames):
        return [
            MappingPipeline.cellname(filename) for filename in filenames
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
    def cellname(filename):
        return os.path.basename(filename).split(os.extsep, 1)[0]

    @staticmethod
    def head_stage(_filename, lines=1000000):
        return [
            'head',
            '-n',
            '{}'.format(lines),
            '-'
        ]

    def bowtie_stage(self, _filename):
        bowtie_ops = self.config.mapping.bowtie_ops.split(' ')
        return [
            'bowtie',
            *bowtie_ops,
            '-m', '1', '--best', '--strata', '--solexa-quals',
            self.config.genome_index_filename(),
            '-'
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
        cellname = self.cellname((filename))
        outfile = os.path.join(
            self.config.mapping_work_dirname(),
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

    def run(self):
        fastq_filenames = self.config.mapping_fastq_filenames()
        assert fastq_filenames

        for filename in fastq_filenames:
            pipeline = [
                *self.unarchive_stage(filename),
                '|',
                *self.head_stage(filename, lines=40000),
                '|',
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
            print(' '.join(pipeline))
            subprocess.check_call(' '.join(pipeline), shell=True)
