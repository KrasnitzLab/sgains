from collections import defaultdict, namedtuple
from sgains.genome import Genome
import os
import time
import pandas as pd
import numpy as np
import pysam
from dask.distributed import as_completed, Queue, worker_client, wait
import gzip
from io import BytesIO


class Varbin10xPipeline(object):

    def __init__(self, config):
        self.config = config
        self.hg = Genome(config)

        self.summary_filename = config.build_data_10x_summary()
        self.bam_filename = config.build_data_10x_bam()
        self.bai_filename = config.build_data_10x_bai()

        assert os.path.exists(self.summary_filename), self.summary_filename
        assert os.path.exists(self.bam_filename), self.bam_filename
        assert os.path.exists(self.bai_filename), self.bai_filename

        self.summary_df = pd.read_csv(self.summary_filename, sep=',')
        self.barcodes = {
            k: v for (k, v) in
            self.summary_df[['barcode', 'cell_id']].to_records(index=False)
        }
        self.bins_filename = self.config.bins_boundaries_filename()
        assert os.path.exists(self.bins_filename), self.bins_filename

        self.bins_df = pd.read_csv(self.bins_filename, sep='\t')
        self.chrom2contig, self.contig2chrom = self._chrom2contig_mapping()
        self.chrom_sizes = self.hg.chrom_sizes()
    
    def _chrom2contig_mapping(self):
        with pysam.AlignmentFile(self.bam_filename, 'rb') as samfile:
            assert samfile.check_index(), \
                (self.bam_filename, self.bai_filename)
            sam_stats = samfile.get_index_statistics()
        chroms = set(self.bins_df['bin.chrom'].values)
        chrom2contig = {}
        contig2chrom = {}
        for stat in sam_stats:
            contig = stat.contig
            chrom_name = contig
            if chrom_name not in chroms:
                chrom_name = "chr{}".format(contig)
            if chrom_name not in chroms:
                continue
            chrom2contig[chrom_name] = contig
            contig2chrom[contig] = chrom_name
        return chrom2contig, contig2chrom

    Region = namedtuple('Region', ['chrom', 'contig', 'start', 'end'])

    def split_bins(self, bins_step, bins_region=None):
        total_bins = len(self.bins_df)
        regions = []
        bin_start = 0
        bin_end = total_bins
        if bins_region is not None:
            bin_start, bin_end = bins_region
            bin_end = min(bin_end, total_bins)

        index = bin_start
        while index < bin_end:
            start = index
            end = index + bins_step - 1
            if end >= total_bins:
                end = total_bins - 1

            chrom_start = self.bins_df.iloc[start, 0]
            chrom_end = self.bins_df.iloc[end, 0]
            while chrom_end != chrom_start:
                end -= 1
                assert end >= start
                chrom_end = self.bins_df.iloc[end, 0]

            pos_start = self.bins_df.iloc[start, 1]
            pos_end = self.bins_df.iloc[end, 3]

            regions.append((chrom_start, pos_start, pos_end))
            index = end + 1

        return [
            self.Region(chrom, self.chrom2contig[chrom], start, end)
            for chrom, start, end in regions
        ]

    @staticmethod
    def compress_reads(reads):
        if not reads:
            return None
        df = pd.DataFrame(reads, columns=['cell_id', 'chrom', 'pos'])
        buff = BytesIO()
        df.to_feather(buff)

        return str(buff.getbuffer(), encoding='latin1')
    
    @staticmethod
    def decompress_reads(data):
        if data is None:
            return None
        buff = BytesIO(bytes(data, encoding='latin1'))
        df = pd.read_feather(buff)
        return df

    Read = namedtuple('Read', ['cell_id', 'chrom', 'pos'])

    def process_region_reads(self, region):
        print(f"started region {region}")
        with pysam.AlignmentFile(self.bam_filename, 'rb') as samfile:
            assert samfile.check_index(), \
                (self.bam_filename, self.bai_filename)
            cells_reads = []
            mapped = 0
            for rec in samfile.fetch(region.contig, region.start, region.end):
                if not rec.has_tag('CB'):
                    continue
                mapped += 1
                barcode = rec.get_tag('CB')
                if barcode not in self.barcodes:
                    continue
                cell_id = self.barcodes[barcode]
                contig = rec.reference_name
                if contig not in self.contig2chrom:
                    continue
                read = self.Read(
                    cell_id, self.contig2chrom[contig], rec.reference_start)
                cells_reads.append(read)
            print(f"done region {region}; reads processed {len(cells_reads)}")
            return self.compress_reads(cells_reads)

    def find_bin_index(self, abspos, bins):
        index = np.searchsorted(
            abspos, bins, side='right')

        index = index - 1
        return index

    def find_bin_index_binsearch(self, bins, abspos):
        index_up = len(bins)
        index_down = 0
        index_mid = int((index_up - index_down) / 2.0)

        while True:
            if abspos >= int(bins[index_mid]):
                index_down = index_mid + 0
                index_mid = int((index_up - index_down) / 2.0) + index_mid
            else:
                index_up = index_mid + 0
                index_mid = int((index_up - index_down) / 2.0) + index_down

            if index_up - index_down < 2:
                break

        return index_down

    def varbin_cell_reads(self, cell_reads):
        assert self.bins_df is not None
        cell_reads = sorted(cell_reads, key=lambda r: f"{r.chrom}:{r.pos:15}")
        count = 0
        dups = 0
        total_reads = 0

        prev_pos = 0
        bin_counts = defaultdict(int)

        bins = self.bins_df['bin.start.abspos'].values

        for read in cell_reads:
            total_reads += 1

            abspos = self.chrom_sizes[read.chrom].abspos + read.pos
            if prev_pos == abspos:
                dups += 1
                continue
            count += 1
            index = self.find_bin_index_binsearch(bins, abspos)

            bin_counts[index] += 1
            prev_pos = abspos

        number_of_reads_per_bin = float(count) / len(self.bins_df)
        result = []
        for index, row in self.bins_df.iterrows():
            bin_count = bin_counts[index]
            ratio = float(bin_count) / number_of_reads_per_bin
            result.append(
                [
                    row['bin.chrom'],
                    row['bin.start'],
                    row['bin.start.abspos'],
                    bin_count,
                    ratio
                ]
            )
        df = pd.DataFrame.from_records(
            result,
            columns=[
                'chrom',
                'chrompos',
                'abspos',
                'bincount',
                'ratio',
            ])
        df.sort_values(by=['abspos'], inplace=True)
        return df

    def varbin_cell_reads_save(self, cell_id, data, outdir):
        reads_df = self.decompress_reads(data)
        cell_reads = []
        for read in reads_df.to_dict(orient='records'):
            cell_reads.append(
                self.Read(
                    cell_id=read['cell_id'],
                    chrom= read['chrom'],
                    pos=read['pos']
                )
            )
        reads_df = None

        df = self.varbin_cell_reads(cell_reads)
        outfile = self.config.varbin_filename(str(cell_id))

        # outfile = os.path.join(outdir, filename)

        df.to_csv(outfile, sep='\t', index=False)
        return cell_id

    def _process_reads_producer(
            self, task_queue, regions):

        with worker_client() as client:
            for index, region in enumerate(regions):
                print(
                    f"submitting task {index} out of {len(regions)} "
                    f"for region {region}")

                task = client.submit(
                    self.process_region_reads, region
                )
                task_queue.put(task)

    def _varbins_cleaner(self, task_queue, total):
        with worker_client() as client:
            processed_tasks = 0
            while task_queue.qsize() > 0 or processed_tasks < total:
                print(f"processed tasks {processed_tasks} out of {total}")
                task = task_queue.get()
                processed_tasks += 1
                cell_id = task.result()
                print(f"processing {cell_id} done")
                client.cancel(task)

    def run(self, dask_client, bins_step=20, bins_region=None, outdir='.'):
        parallel = int(self.config.parallel)
        maxsize = int(2 * parallel)
        task_queue = Queue(maxsize=maxsize)

        regions = self.split_bins(bins_step=bins_step, bins_region=bins_region)

        producer = dask_client.submit(
            self._process_reads_producer, task_queue, 
            regions
        )

        print("producer started...")
        while task_queue.qsize() == 0:
            print("task queue empty. sleeping...")
            time.sleep(10)

        reads_data = []
        processed_tasks = 0
        while task_queue.qsize() > 0 or processed_tasks < len(regions):

            task = task_queue.get()
            result = dask_client.gather(task)
            start = time.time()
            print(f"task finished...")
            result = task.result()
            # result = dask_client.gather(task)
            elapsed = time.time() - start
            print(f"results collected in {elapsed}...")

            reads_data.append(result)
            processed_tasks += 1
            print(f'processed {processed_tasks} merge reads')
            dask_client.cancel(task)

        print("about to wait for producer")
        wait(producer)
        print("[done]")

        cells_reads = defaultdict(list)
        while len(reads_data) > 0:
            print(f"reads data processing {len(reads_data)}")
            data = reads_data.pop()
            df = self.decompress_reads(data)
            if df is None:
                continue
            for rec in df.to_dict(orient='records'):
                read = self.Read(rec['cell_id'], rec['chrom'], rec['pos'])
                cells_reads[read.cell_id].append(read)
            df = None
            data = None
        
        cleaner = dask_client.submit(
            self._varbins_cleaner, task_queue, 
            len(cells_reads)
        )

        for cell_id, cell_reads in cells_reads.items():
            print(f"about to start varbins for cell {cell_id}")
            data = self.compress_reads(cell_reads)
            task = dask_client.submit(
                self.varbin_cell_reads_save, cell_id, data, outdir
            )
            print(f"starting processing varbins for cell {cell_id}")
            task_queue.put(task)

        print("about to wait for cleaner")
        wait(cleaner)
        print("[done]")
