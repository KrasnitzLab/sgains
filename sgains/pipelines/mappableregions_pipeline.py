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


class InputGeneratorThread(Thread):

    def __init__(self, control_queue, aligner_input, input_function_generator):
        super(InputGeneratorThread, self).__init__()

        self.control_queue = control_queue
        self.aligner_input = aligner_input
        self.input_function_generator = input_function_generator

    def run(self):
        for rec in self.input_function_generator:
            out = Genome.to_fasta_string(rec)

            self.aligner_input.write(out)
            self.aligner_input.flush()

        self.aligner_input.close()


class AlignerOutputProcessingThread(Thread):

    def __init__(self, control_queue, aligner_output, output_process_function):
        super(AlignerOutputProcessingThread, self).__init__()

        self.control_queue = control_queue
        self.aligner_output = aligner_output
        self.output_process_function = output_process_function

    def run(self):
        prev = None
        state = MappableState.OUT

        while True:
            line = self.aligner_output.readline()
            if not line:
                break

            line = line.decode()
            if line[0] == '@':
                # comment
                continue
            mapping = Mapping.parse_sam(line)

            if state == MappableState.OUT:
                if mapping.acceptable():
                    prev = MappableRegion(mapping)
                    state = MappableState.IN
            else:
                if mapping.acceptable():
                    if mapping.chrom == prev.chrom:
                        if mapping.start < prev.end:
                            print(
                                "WARN: prev=", prev,
                                "; mapping=", mapping,
                                "; line=", line)
                        prev.extend(mapping)
                    else:

                        if prev.start >= prev.end:
                            print("ERROR: prev=", prev, "; line=", line)
                        # assert prev.start < prev.end
                        self.output_process_function(prev)
                        prev = MappableRegion(mapping)
                else:
                    if prev.start >= prev.end:
                        print("ERROR: prev=", prev, "; line=", line)
                    # assert prev.start < prev.end
                    self.output_process_function(prev)
                    state = MappableState.OUT

        if state == MappableState.IN:
            if prev.start >= prev.end:
                print("ERROR: prev=", prev, "; line=", line)
            # assert prev.start < prev.end
            self.output_process_function(prev)

        self.aligner_output.close()
        self.control_queue.put("out_done")


class MappableRegionsPipeline(object):

    def __init__(self, config, aligner=None):
        self.config = config
        self.genome = Genome(self.config)
        if aligner is not None:
            self.aligner = aligner
        else:
            assert self.genome.aligner is not None
            self.aligner = self.genome.aligner
        assert self.aligner is not None

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
                seq_record = self.genome.load_chrom(chrom)
                for i in range(len(seq_record) - read_length + 1):
                    seq = seq_record.seq[i: i + read_length]
                    out_record = SeqRecord(
                        seq,
                        id="{}.{}".format(chrom, i + 1),
                        description="generated_read"
                    )
                    # if 'N' in seq:
                    #     print('skipping: ', out_record)
                    #     continue
                    yield out_record
        finally:
            pass

    def generate_mappable_regions(
            self, chroms, read_length,
            outfile=None, aligner_options=[]):

        if outfile is None:
            outfile = sys.stdout

        reads_generator = self.generate_reads(chroms, read_length)

        def aligner_output_process_function(line):
            outfile.write(str(line))
            outfile.write("\n")

        aligner_command = self.aligner.build_mappable_regions_command(
            options=aligner_options
        )
        print('aligner command', ' '.join(aligner_command))

        with Popen(aligner_command, stdout=PIPE, stdin=PIPE) as proc:

            control_queue = queue.Queue()
            input_thread = InputGeneratorThread(
                control_queue,
                proc.stdin,
                reads_generator)
            output_thread = AlignerOutputProcessingThread(
                control_queue,
                proc.stdout,
                aligner_output_process_function)

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

    def mappable_regions_chrom_filename(self, chrom):
        mname = "{}_{}".format(
            chrom, self.config.mappable_regions.mappable_file)
        filename = os.path.join(
            self.config.mappable_regions.mappable_dir,
            mname
        )
        return filename

    def mappable_regions_filename(self):
        mname = self.config.mappable_regions.mappable_file
        filename = os.path.join(
            self.config.mappable_regions.mappable_dir,
            mname
        )
        return filename

    def run_once(self, chrom):
        outfilename = self.mappable_regions_chrom_filename(chrom)
        with open(outfilename, "w") as outfile:
            self.generate_mappable_regions(
                [chrom], read_length=50,
                outfile=outfile)
        return outfilename

    def concatenate_all_chroms(self):
        dst = self.mappable_regions_filename()
        if os.path.exists(dst) and not self.config.force:
            print(colored(
                "destination mappable regions file already exists"
                "use --force to overwrite", "red"))
            raise ValueError("destination file exists... use --force")

        if not self.config.dry_run:
            with open(dst, 'wb') as output:
                for chrom in self.genome.version.CHROMS:
                    src = self.mappable_regions_chrom_filename(chrom)
                    print(colored(
                        "appending {} to {}".format(src, dst),
                        "green"))
                    with open(src, 'rb') as src:
                        if not self.config.dry_run:
                            shutil.copyfileobj(src, output, 1024 * 1024 * 10)

    def run(self, dask_client):
        outfilename = self.mappable_regions_filename()
        print(colored(
            "going to generate mappable regions with length {} "
            "from genome {} into {}".format(
                self.config.mappable_regions.mappable_read_length,
                self.config.genome.genome_dir,
                outfilename
            ),
            "green"))

        if os.path.exists(outfilename) and not self.config.force:
            print(colored(
                "output file {} already exists; "
                "use --force to overwrite".format(outfilename),
                "red"))
            raise ValueError("output file already exists")

        genome_index_filenames = self.aligner.genome_index_filenames
        if not os.path.exists(genome_index_filenames[0]):
            print(colored(
                "genome index file {} not found".format(
                    genome_index_filenames),
                "red"))
            raise ValueError("genome index file not found")

        if self.config.dry_run:
            return

        os.makedirs(self.config.mappable_regions.mappable_dir, exist_ok=True)

        assert dask_client

        delayed_tasks = dask_client.map(
                self.run_once, self.genome.version.CHROMS)

        distributed.wait(delayed_tasks)

        for fut in delayed_tasks:
            print("fut done:", fut.done())
            print("fut exception:", fut.exception())
            print("fut traceback:", fut.traceback())
            print("fut result:", fut.result())

            # if fut.traceback() is not None:
            #     traceback.print_tb(fut.traceback())
            # if fut.exception() is None:
            #     print(fut.result())

        self.concatenate_all_chroms()
