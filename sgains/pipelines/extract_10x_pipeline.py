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


class Base10xPipeline(object):

    def __init__(self, config):
        self.config = config
        self.summary_filename = self.config.data_10x.data_10x_cell_summary
        self.reads_dirname = self.config.reads.reads_dir
        self.bam_filename = self.config.data_10x.data_10x_bam
        self.bai_filename = self.config.data_10x.data_10x_bai

        assert os.path.exists(self.summary_filename), self.summary_filename
        assert os.path.exists(self.bam_filename), self.bam_filename
        assert os.path.exists(self.bai_filename), self.bai_filename

        self.summary_df = pd.read_csv(self.summary_filename, sep=',')
        print(self.summary_df.head())

        self.barcodes = {
            k: v for (k, v) in
            self.summary_df[['barcode', 'cell_id']].to_records(index=False)
        }
        self.genome = Genome(self.config)
        assert self.genome is not None


class Extract10xPipeline(Base10xPipeline):

    def __init__(self, config):
        super(Extract10xPipeline, self).__init__(config)

    def _build_segment_regions(self, segment_len=50_000_000):
        
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
            count = max(count, 1)
            segment_size = size // count
            for seg_index in range(count):
                begin_pos = seg_index * segment_size
                end_pos = (seg_index + 1) * segment_size - 1
                if seg_index + 1 == count:
                    end_pos = size
                index += 1
                segments.append((index, (contig, begin_pos, end_pos)))
        return segments

    PROGRESS_STEP = 10_000
    FLUSH_STEP = 2_000

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
                if barcode not in self.barcodes:
                    continue

                cell_records[barcode].append(rec)
                if len(cell_records[barcode]) > self.FLUSH_STEP:
                    if barcode in self.barcodes:
                        self._write_to_fastq(
                            barcode, segment, cell_records[barcode])
                    cell_records[barcode] = []

            for barcode, records in cell_records.items():
                if barcode in self.barcodes:
                    self._write_to_fastq(barcode, segment, records)
                cell_records[barcode] = []

    @staticmethod
    def _cell_name(cell_id):
        cell_name = f"C{cell_id:0>8}"
        return cell_name

    def _cell_segment_filename(self, cell_id, segment):
        cell_name = self._cell_name(cell_id)
        segment_name = f"{segment:0>8}"
        filename = os.path.join(
            self.reads_dirname,
            cell_name,
            f"{cell_name}_{segment_name}{self.config.reads.reads_suffix}"
        )
        return filename

    def _cell_fragments(self, cell_id):
        cell_name = self._cell_name(cell_id)

        pattern = f"{cell_name}_*{self.config.reads.reads_suffix}"
        pattern = os.path.join(
            self.reads_dirname,
            cell_name,
            pattern
        )
        print(colored(
            "merging fastq cell files {} ".format(
                pattern
            ),
            "green"))
        filenames = glob.glob(pattern)
        return filenames

    def _write_to_fastq(self, barcode, segment, records):

        cell = self.barcodes[barcode]
        filename = self._cell_segment_filename(cell, segment)
        os.makedirs(os.path.dirname(filename), exist_ok=True)

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
        filenames = self._cell_fragments(cell_id)

        if not filenames:
            print(colored(
                f"can't find segments for cell {cell_id}", 'red'))
            return []
        cell_name = self._cell_name(cell_id)
        outfile = os.path.join(
            self.reads_dirname,
            f"{cell_name}{self.config.reads.reads_suffix}")

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

    def run(self, dask_client):
        try:
            self.config.check_nonempty_workdir(self.reads_dirname)
        except NonEmptyWorkDirectory:
            return

        command = "rm -rf {}/*".format(self.reads_dirname)
        print(colored(command, 'yellow'))
        os.system(command)
        os.makedirs(self.reads_dirname, exist_ok=True)

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

