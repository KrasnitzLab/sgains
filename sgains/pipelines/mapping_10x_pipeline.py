'''
Created on Jul 13, 2017

@author: lubo
'''
import os
import gzip
import glob
import subprocess
from collections import defaultdict

import pysam

import pandas as pd
from dask import distributed
from termcolor import colored

from sgains.config import NonEmptyWorkDirectory
from sgains.pipelines.mapping_pipeline import MappingPipeline
from sgains.genome import Genome


class Mapping10xPipeline(object):

    def __init__(self, config):
        self.config = config
        self.summary_filename = config.build_data_10x_summary()
        self.fastq_dirname = config.build_mapping_10x_fastqdir()
        self.bam_filename = config.build_data_10x_bam()
        self.bai_filename = config.build_data_10x_bai()
        self.mapping_dirname = config.build_mapping_10x_dir()

        assert os.path.exists(self.summary_filename), self.summary_filename
        assert os.path.exists(self.bam_filename), self.bam_filename
        assert os.path.exists(self.bai_filename), self.bai_filename

        self.summary_df = pd.read_csv(self.summary_filename, sep=',')
        self.barcodes = {
            k: v for (k, v) in
            self.summary_df[['barcode', 'cell_id']].to_records(index=False)
        }
        self.genome = Genome(self.config)
        assert self.genome is not None

    def _build_segment_regions(self, segment_len=25_000_000):
        
        with pysam.AlignmentFile(self.bam_filename, 'rb') as samfile:
            assert samfile.check_index(), \
                (self.bam_filename, self.bai_filename)
            print(samfile.get_index_statistics())
            sam_stats = samfile.get_index_statistics()
        chrom_sizes = self.genome.chrom_sizes()
        segments = []
        index = 0
        for stat in sam_stats:
            contig = stat.contig
            contig_name = contig
            if contig_name not in chrom_sizes:
                contig_name = "chr{}".format(contig)
            if contig_name not in chrom_sizes:
                continue
            size = chrom_sizes[contig_name]['size']
            count = size // segment_len
            count = max(count, 0)
            segment_size = size // count
            for seg_index in range(count):
                begin_pos = seg_index * segment_size
                end_pos = (seg_index + 1) * segment_size - 1
                if seg_index + 1 == count:
                    end_pos = size
                index += 1
                segments.append((index, (contig, begin_pos, end_pos)))
        return segments

    PROGRESS_STEP = 50_000
    FLUSH_STEP = 500_000

    def _process_segment(self, segment, region):
        with pysam.AlignmentFile(self.bam_filename, 'rb') as samfile:
            assert samfile.check_index(), \
                (self.bam_filename, self.bai_filename)
            chrom, begin, end = region

            mapped = 0
            cell_records = defaultdict(lambda: [])
            for rec in samfile.fetch(chrom, begin, end):
                if not rec.has_tag('CB'):
                    continue
                mapped += 1
                barcode = rec.get_tag('CB')
                cell_records[barcode].append(rec)
                if len(cell_records[barcode]) > self.FLUSH_STEP \
                        and barcode in self.barcodes:
                    self._write_to_fastq(
                        barcode, segment, cell_records[barcode])
                    cell_records[barcode] = []

            for barcode, records in cell_records.items():
                if barcode in self.barcodes:
                    self._write_to_fastq(barcode, segment, records)

    def _write_to_fastq(self, barcode, segment, records):
        filename = "C{:0>5}_{:0>5}.fastq.gz".format(
            self.barcodes[barcode], segment)
        filename = os.path.join(self.fastq_dirname, filename)
        with gzip.open(filename, "at") as outfile:
            for rec in records:
                assert barcode == rec.get_tag('CB')
                outfile.write('@')
                outfile.write(rec.query_name)
                outfile.write('\n')
                outfile.write(rec.query_sequence)
                outfile.write('\n')
                outfile.write('+\n')
                outfile.write(''.join(
                    [chr(ch + 33) for ch in rec.query_qualities]))
                outfile.write('\n')

    def _split_once(self, dry_run, segment):
        segment_index, region = segment
        contig, begin, end = region
        print(colored(
            "processing 10x bam file {} segment {}, {}:{}-{}".format(
                self.bam_filename,
                segment_index,
                contig,
                begin,
                end
            ),
            "green"))
        if dry_run:
            return
        self._process_segment(segment_index, region)

    def _merge_cell(self, dry_run, cell_id):
        pattern = "C{:0>5}_*.fastq.gz".format(cell_id)
        pattern = os.path.join(
            self.fastq_dirname,
            pattern
        )
        print(colored(
            "merging fastq cell files {} ".format(
                pattern
            ),
            "green"))
        if dry_run:
            return []

        filenames = glob.glob(pattern)
        if not filenames:
            print(colored(
                "can't find segments for cell {}: {}".format(
                    cell_id, pattern
                ), 'red'
            ))
            return []
        print(filenames)
        outfile = "C{:0>5}.fastq.gz".format(cell_id)
        outfile = os.path.join(self.fastq_dirname, outfile)

        command = "cat {} > {}".format(
            " ".join(filenames),
            outfile
        )
        print(colored(command, "yellow"))
        os.system(command)
        command = "rm -f {}".format(" ".join(filenames))
        print(colored(command, "yellow"))
        os.system(command)
        return filenames

    def _bowtie_stage(self, cellname):
        reportfile = os.path.join(
            self.mapping_dirname,
            "{}.bowtie_report.log".format(cellname)
        )

        bowtie_opts = self.config.mapping_10x.\
            mapping_10x_bowtie_opts.split(' ')
        return [
            'bowtie',
            '-S', '-t', '-v', '0', '-m', '1',
            '--best', '--strata', '--chunkmbs', '256',
            *bowtie_opts,
            self.config.genome_index_filename(),
            '-',
            '2>',
            reportfile,
        ]

    def _samtools_view_store_stage(self, cellname):
        outfile = os.path.join(
            self.mapping_dirname,
            "{}.rmdup.bam".format(cellname)
        )
        return [
            'samtools',
            'view',
            '-b',
            '-o',
            outfile,
            '-',
        ]

    def _samtools_index_bam(self, cellname):
        outfile = os.path.join(
            self.mapping_dirname,
            "{}.rmdup.bam".format(cellname)
        )
        return [
            'samtools'
            'index',
            outfile,
        ]

    def _mapping_once(self, dry_run, cell_id):
        cellname = "C{:0>5}".format(cell_id)
        fastq_filename = os.path.join(
            self.fastq_dirname,
            "{}.fastq.gz".format(cellname))
        pipeline = [
            *MappingPipeline.unarchive_stage(fastq_filename),
            '|',
            # *MappingPipeline.head_stage(cellname, lines=40000),
            # '|',
            *self._bowtie_stage(cellname),
            '|',
            *MappingPipeline.samtools_view_stage(cellname),
            '|',
            *MappingPipeline.samtools_view_remove_unmappted_stage(cellname),
            '|',
            *MappingPipeline.samtools_sort_stage(cellname),
            '|',
            *MappingPipeline.samtools_rmdup_stage(cellname),
            '|',
            *self._samtools_view_store_stage(cellname),
        ]
        command = ' '.join(pipeline)
        print(colored(command, "green"))
        if not dry_run:
            res = subprocess.check_call(
                command,
                # stdout=subprocess.DEVNULL,
                # stderr=subprocess.DEVNULL,
                shell=True)
            print(res)

    def run(self, dask_client):
        try:
            self.config.check_nonempty_workdir(self.fastq_dirname)
        except NonEmptyWorkDirectory:
            return
        try:
            self.config.check_nonempty_workdir(
                self.mapping_dirname
            )
        except NonEmptyWorkDirectory:
            return

        command = "rm -rf {}/*".format(self.fastq_dirname)
        print(colored(command, 'yellow'))
        os.system(command)

        command = "rm -rf {}/*".format(self.mapping_dirname)
        print(colored(command, 'yellow'))
        os.system(command)

        segments = self._build_segment_regions(segment_len=50_000_000)
        delayed_tasks = dask_client.map(
            lambda segment: self._split_once(self.config.dry_run, segment),
            segments
        )
        distributed.wait(delayed_tasks)

        cells = self.barcodes.values()
        delayed_tasks = dask_client.map(
            lambda cell_id: self._merge_cell(self.config.dry_run, cell_id),
            cells
        )
        distributed.wait(delayed_tasks)

        cells = self.barcodes.values()
        delayed_tasks = dask_client.map(
            lambda cell_id: self._mapping_once(self.config.dry_run, cell_id),
            cells
        )
        distributed.wait(delayed_tasks)
