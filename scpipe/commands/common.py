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

    def common_updates(self, args):
        if args.config is not None:
            config = Config.load(args.config)
            self.config.update(config)

        if args.dry_run is not None:
            self.config.dry_run = args.dry_run
        if args.force is not None:
            self.config.force = args.force


class DataDirMixin(object):
    __slots__ = ()

    def data_dir_options(self, glob=False):
        assert self.subconfig is not None
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "input data options")
        group.add_argument(
            "--data-dir", "-i",
            dest="data_dir",
            help="input data directory where the input data is located",
            default=self.subconfig.data_dir
        )
        if glob:
            group.add_argument(
                "--glob", "-g",
                dest="data_glob",
                help="glob pattern for finding input data",
                default=self.subconfig.data_glob)

        return group

    def data_dir_update(self, args, glob=False):
        assert self.subconfig is not None
        assert self.subparser is not None

        if args.data_dir is not None:
            self.subconfig.data_dir = args.data_dir
        if glob:
            if args.data_glob is not None:
                self.subconfig.data_glob = args.data_glob


class WorkDirMixin(object):
    __slots__ = ()

    def work_dir_options(self):
        assert self.subconfig is not None
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "output data options")
        group.add_argument(
            "--work-dir", "-o",
            dest="work_dir",
            help="output directory where results from processing are stored",
            default=self.subconfig.work_dir
        )
        return group

    def work_dir_update(self, args):
        if args.work_dir is not None:
            self.subconfig.work_dir = args.work_dir


class GenomeIndexMixin(object):
    __slots__ = ()

    def genome_index_options(self, genome_dir=False):
        assert self.subconfig is not None
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "genome index options")
        group.add_argument(
            "--genome-index", "-G",
            dest="genome_index",
            help="genome index name",
            default=self.config.genome.index)

        if genome_dir:
            group.add_argument(
                "--genome-dir",
                dest="genome_dir",
                help="genome index directory",
                default=self.config.genome.work_dir)
        return group

    def genome_index_update(self, args, genome_dir=False):
        if args.genome_index is not None:
            self.config.genome.index = args.genome_index
        if genome_dir:
            if args.genome_dir is not None:
                self.config.genome.work_dir = args.genome_dir


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
            "--mappable-dir", "-o",
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
