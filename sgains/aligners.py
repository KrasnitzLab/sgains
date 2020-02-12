import os

from typing import List, Optional

from sgains.config import Config


class Aligner:

    def __init__(self, config, genome_version):
        self.config = config
        self.genome_version = genome_version

    @classmethod
    def from_config(cls, config, genome_version) -> Optional['Aligner']:
        assert config is not None
        if config.genome.aligner == 'bowtie':
            return Bowtie(config, genome_version)
        elif config.genome.aligner == 'hisat2':
            return Hisat2(config, genome_version)

        assert False, f'Unsupported aligner {config.genome.aligner}'

    @property
    def name(self):
        raise NotImplementedError()

    @property
    def default_genome_sequence(self) -> str:
        return self.genome_version.sequence_filename

    @property
    def default_genome_index_prefix(self) -> str:
        return self.genome_version.index_prefix

    @property
    def genome_index_filenames(self) -> List[str]:
        raise NotImplementedError()

    def build_index_command(
            self, sequence_filename: str,
            index_prefix: str) -> List[str]:
        raise NotImplementedError()

    def build_align_reads_command(self, options: List[str] = []) -> List[str]:
        raise NotImplementedError()

    def build_mapping_command(
            self, fastq_filename: str, options: List[str] = []) -> List[str]:
        raise NotImplementedError()

    def build_mappable_regions_command(
            self, options: List[str] = []) -> List[str]:
        raise NotImplementedError()

    def build_postmapping_filter(self) -> List[str]:
        raise NotImplementedError()

    def build_fastq_stream_command(
            self, fastq_filename: str, num_lines: int = 0) -> List[str]:
        _base, ext = os.path.splitext(fastq_filename)
        if ext == '.gz':
            result = ['gunzip', '-c', fastq_filename]
        elif ext == '.fastq':
            result = ['cat', fastq_filename]
        elif ext == '.bz2':
            result = ['bzcat', '-c', fastq_filename]
        else:
            assert False, f'unexptected archive extention: ' \
                f'{ext}, {fastq_filename}'
        if num_lines > 0:
            result.extend([
                '|', 'head', '-n', str(num_lines)
            ])
        return result

    def build_report_filename(self, fastq_filename: str) -> str:
        cellname = Config.cellname(fastq_filename)
        reportfile = self.config.abspath(os.path.join(
            self.config.mapping_dirname(),
            f"{cellname}.aligner_report.log"
        ))
        return reportfile

    def build_output_samfile(self, fastq_filename: str) -> str:
        cellname = Config.cellname((fastq_filename))
        outfile = os.path.join(
            self.config.mapping_dirname(),
            "{}.rmdup.bam".format(cellname)
        )
        outfile = self.config.abspath(outfile)
        return outfile

    def build_samtools_store_command(self, fastq_filename: str) -> List[str]:
        outfile = self.build_output_samfile(fastq_filename)
        return [
            'samtools',
            'view',
            '-b',
            '-o',
            outfile,
        ]

    def build_samtools_index_command(self, fastq_filename: str) -> List[str]:
        outfile = self.build_output_samfile(fastq_filename)
        return [
            'samtools',
            'index',
            outfile
        ]

    def build_mapping_pipeline(
            self, fastq_filename: str,
            num_lines: int = 0,
            options: List[str] = []) -> List[str]:

        pipeline = [
            self.build_fastq_stream_command(
                fastq_filename, num_lines=num_lines),
            self.build_mapping_command(fastq_filename, options=options),
            # post mapping filter
            self.build_postmapping_filter(),
            # sort
            ['samtools', 'sort'],
            # rmdup
            ['samtools', 'rmdup', '-s', '-', '-'],
            self.build_samtools_store_command(fastq_filename),
        ]
        pipeline = [
            ' '.join(command) for command in pipeline
        ]
        return pipeline


class Bowtie(Aligner):

    def __init__(self, config, genome_version):
        super(Bowtie, self).__init__(config, genome_version)

    @property
    def name(self):
        return 'bowtie'

    def build_index_command(
            self, sequence_filename: Optional[str] = None,
            index_prefix: Optional[str] = None,
            options: List[str] = []) -> List[str]:

        if sequence_filename is None:
            sequence_filename = self.genome_version.sequence_filename
        if index_prefix is None:
            index_prefix = self.genome_version.index_prefix

        result = [
            "bowtie-build", *options,
            "-f", sequence_filename,
            index_prefix,
        ]
        return result

    def build_align_reads_command(
            self, options: List[str] = []) -> List[str]:
        return []

    @property
    def genome_index_filenames(
            self,
            genome_index_prefix: Optional[str] = None) -> List[str]:
        if genome_index_prefix is None:
            genome_index_prefix = self.default_genome_index_prefix

        index_files = [
            f"{genome_index_prefix}.{index}.ebwt"
            for index in range(1, 5)
        ]
        rev_index_files = [
            f"{genome_index_prefix}.rev.{index}.ebwt"
            for index in range(1, 3)
        ]
        return [index_files, *rev_index_files]

    def build_mapping_command(
            self, fastq_filename: str, options: List[str] = []) -> List[str]:
        reportfile = self.build_report_filename(fastq_filename)
        return [
            'bowtie',
            '-S', '-t', '-v', '0', '-m', '1',
            '--best', '--strata', '--chunkmbs', '256',
            *options,
            self.default_genome_index_prefix,
            '-',
            '2>',
            reportfile,
        ]

    def build_mappable_regions_command(
            self, options: List[str] = []) -> List[str]:
        return [
            'bowtie', '-S', '-t', '-v', '0', '-m', '1',
            *options,
            '-f', self.default_genome_index_prefix,
            '-',
        ]

    def build_postmapping_filter(self) -> List[str]:
        return [
            'samtools',
            'view',
            '-bu',
            '-F',
            '4',
        ]


class Hisat2(Aligner):

    def __init__(self, config, genome_version):
        super(Hisat2, self).__init__(config, genome_version)

    @property
    def name(self):
        return 'hisat2'

    def build_index_command(
            self, sequence_filename: str,
            index_prefix: str) -> List[str]:

        if sequence_filename is None:
            sequence_filename = self.genome_version.sequence_filename
        if index_prefix is None:
            index_prefix = self.genome_version.index_prefix

        result = [
            "hisat2-build",
            self.genome_version.sequence_filename,
            self.genome_version.index_prefix,
        ]
        return result

    @property
    def genome_index_filenames(
            self,
            genome_index_prefix: Optional[str] = None) -> List[str]:
        if genome_index_prefix is None:
            genome_index_prefix = self.default_genome_index_prefix

        result = [
            f"{genome_index_prefix}.{index}.ht2"
            for index in range(1, 9)
        ]
        return result

    def build_mapping_command(
            self, fastq_filename: str, options: List[str] = []) -> List[str]:
        reportfile = self.build_report_filename(fastq_filename)
        return [
            'hisat2',
            '-X', '1500',
            '--no-spliced-alignment',
            *options,
            '-x', self.default_genome_index_prefix,
            '-U', '-',
            '2>', reportfile,
        ]

    def build_mappable_regions_command(
            self, options: List[str] = []) -> List[str]:
        return [
            'hisat2',
            '-X', '1500',
            '--no-spliced-alignment',
            *options,
            '-x', self.default_genome_index_prefix,
            '-f', '-',
        ]

    def build_postmapping_filter(self) -> List[str]:
        return [
            'samtools',
            'view',
            '-bu',
            '-q', '30',
            '-F', '0xff00',
        ]