'''
Created on Aug 2, 2017

@author: lubo
'''
from config import Config


class OptionsBase(object):

    def __init__(self, config):
        self.config = config
        self.subconfig = None
        self.parser = None
        self.subparser = None

    @staticmethod
    def common_options(parser):
        parser.add_argument(
            "-v", "--verbose",
            dest="verbose",
            action="count",
            help="set verbosity level [default: %(default)s]",
            default=0
        )
        parser.add_argument(
            "-c", "--config",
            dest="config",
            help="configuration file",
            metavar="path"
        )

        parser.add_argument(
            "-n", "--dry-run",
            dest="dry_run",
            action="store_true",
            help="perform a trial run with no changes made",
            default=False
        )

        parser.add_argument(
            "--force", "-F",
            dest="force",
            action="store_true",
            help="allows overwriting nonempty results directory",
            default=False
        )

    def common_updates(self, args):
        if args.config is not None:
            config = Config.load(args.config)
            print("config loaded: {}".format(config))
            self.config = config
            print("config updated: {}".format(self.config))

        if args.dry_run is not None:
            self.config.dry_run = args.dry_run
        if args.force is not None:
            self.config.force = args.force


class DataDirMixin(object):
    __slots__ = ()

    def data_dir_options(self, config=None, glob=False):
        assert self.subparser is not None
        if config is None:
            config = self.subconfig

        group = self.subparser.add_argument_group(
            "input data options")
        group.add_argument(
            "--data-dir", "-i",
            dest="data_dir",
            help="input data directory where the input data is located",
            default=config.data_dir
        )
        if glob:
            group.add_argument(
                "--glob", "-g",
                dest="data_glob",
                help="glob pattern for finding input data",
                default=config.data_glob)

        return group

    def data_dir_update(self, args, config=None, glob=False):
        assert self.subparser is not None
        if config is None:
            config = self.subconfig

        if args.data_dir is not None:
            config.data_dir = args.data_dir
        if glob:
            if args.data_glob is not None:
                config.data_glob = args.data_glob


class WorkDirMixin(object):
    __slots__ = ()

    def work_dir_options(self, config=None):
        assert self.subparser is not None
        if config is None:
            config = self.subconfig

        group = self.subparser.add_argument_group(
            "output data options")
        group.add_argument(
            "--work-dir", "-o",
            dest="work_dir",
            help="output directory where results from processing are stored",
            default=config.work_dir
        )
        return group

    def work_dir_update(self, args, config=None):
        if config is None:
            config = self.subconfig

        if args.work_dir is not None:
            config.work_dir = args.work_dir


class GenomeIndexMixin(object):
    __slots__ = ()

    def genome_index_options(self, input_dir=False):
        assert self.subconfig is not None
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "genome index options")
        group.add_argument(
            "--genome-index", "-G",
            dest="genome_index",
            help="genome index name",
            default=self.config.genome.index)

        group.add_argument(
            "--genome-dir",
            dest="genome_dir",
            help="genome index directory",
            default=self.config.genome.work_dir)

        group.add_argument(
            "--genome-version",
            dest="genome_version",
            help="version of reference genome in use (supports only hg19)",
            default=self.config.genome.version)
        if input_dir:
            group.add_argument(
                "--genome-pristine",
                dest="genome_pristine",
                help="directory where clean copy of reference genome "
                "is located",
                default=self.config.genome.data_dir)
        return group

    def genome_index_update(self, args, input_dir=False):
        if args.genome_index is not None:
            self.config.genome.index = args.genome_index
        if args.genome_dir is not None:
            self.config.genome.work_dir = args.genome_dir
        if args.genome_version is not None:
            self.config.genome.version = args.genome_version

        if input_dir:
            if args.genome_pristine is not None:
                self.config.genome.data_dir = args.genome_pristine


class BinsBoundariesMixin(object):
    __slots__ = ()

    def bins_boundaries_options(self, bins_count=False):
        group = self.subparser.add_argument_group(
            "bins boundaries")
        group.add_argument(
            "--bins-boundaries", "-B",
            dest="bins_boundaries",
            help="bins boundaries filename",
            default=self.config.bins.bins_boundaries
        )
        group.add_argument(
            "--bins-dir",
            dest="bins_dir",
            help="bins working directory",
            default=self.config.bins.work_dir
        )
        if bins_count:
            group.add_argument(
                "--bins-count", "-C",
                dest="bins_count",
                type=int,
                help="number of bins",
                default=self.config.bins.bins_count
            )

    def bins_boundaries_updates(self, args, bins_count=False):
        if args.bins_boundaries:
            self.config.bins.bins_boundaries = args.bins_boundaries
        if args.bins_dir:
            self.config.bins.work_dir = args.bins_dir
        if bins_count:
            if args.bins_count:
                self.config.bins.bins_count = args.bins_count


class MappableRegionsMixin(object):
    __slots__ = ()

    def mappable_regions_options(self, read_length=False):
        group = self.subparser.add_argument_group(
            "mappable regions options")
        group.add_argument(
            "--mappable-dir", "-m",
            dest="mappable_dir",
            help="directory where mappable regions file is stroed",
            default=self.config.mappable_regions.work_dir
        )
        group.add_argument(
            "--mappable-regions", "-M",
            dest="mappable_regions",
            help="filename where mappable regions are stored",
            default=self.config.mappable_regions.mappable_regions
        )
        if read_length:
            group.add_argument(
                "--read-length", "-l",
                dest="length",
                type=int,
                help="read length to use for generation of mappable regions",
                default=self.config.mappable_regions.length
            )
            group.add_argument(
                "--bowtie-opts",
                dest="bowtie_opts",
                help="additional bowtie options",
                default=self.config.mappable_regions.bowtie_opts
            )
        return group

    def mappable_regions_update(self, args, read_length=False):
        if args.mappable_dir is not None:
            self.config.mappable_regions.work_dir = args.mappable_dir
        if args.mappable_regions:
            self.config.mappable_regions.mappable_regions = \
                args.mappable_regions
        if read_length:
            if args.bowtie_opts:
                self.config.mappable_regions.bowtie_opts = args.bowtie_opts
            if args.length:
                self.config.mappable_regions.length = args.length
