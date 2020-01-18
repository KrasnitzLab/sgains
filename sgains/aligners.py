from typing import List, Optional


class Aligner:

    def __init__(self, genome_version):
        self.genome_version = genome_version

    @classmethod
    def from_config(cls, config, genome_version) -> Optional['Aligner']:
        assert config is not None
        if config.genome.aligner == 'bowtie':
            return Bowtie(genome_version)
        elif config.genome.aligner == 'hist2':
            return Hisat2(genome_version)
        else:
            return None

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

    def build_mapping_command(self, options: List[str] = []) -> List[str]:
        raise NotImplementedError()


class Bowtie(Aligner):

    def __init__(self, genome_version):
        super(Bowtie, self).__init__(genome_version)

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
        raise NotImplementedError()

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

    def build_mapping_command(self, options: List[str] = []) -> List[str]:
        return [
            'bowtie',
            '-S', '-t', '-v', '0', '-m', '1',
            '--best', '--strata', '--chunkmbs', '256',
            *options,
            self.default_genome_index_prefix,            
        ]


class Hisat2(Aligner):

    def __init__(self, genome_version):
        super(Bowtie, self).__init__(genome_version)

    @property
    def name(self):
        return 'hisat2'

    def build_index_command(
            self, sequence_filename: str,
            index_prefix: str) -> List[str]:

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

    def build_mapping_command(self, options: List[str] = []) -> List[str]:
        return [
            'hisat2',
            '-X', '1500',
            '--no-spliced-alignment',
            *options,
            self.default_genome_index_prefix,            
        ]
