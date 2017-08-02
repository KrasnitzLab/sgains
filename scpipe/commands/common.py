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

    def bins_boundaries_options(self):
        group = self.subparser.add_argument_group(
            "bins boundaries")
        group.add_argument(
            "--bins-boundaries", "-B",
            dest="bins_boundaries",
            help="bins boundaries filename",
            default=self.subconfig.bins_boundaries)

        group.add_argument(
            "--bins-dir",
            dest="bins_dir",
            help="bins working directory",
            default=self.subconfig.work_dir)

    def bins_boundaries_updates(self, args):
        if args.bins_boundaries:
            self.subconfig.bins_boundaries = args.bins_boundaries
        if args.bins_dir:
            self.subconfig.work_dir = args.bins_dir
