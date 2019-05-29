'''
Created on Aug 2, 2017

@author: lubo
'''
from distributed import Client, LocalCluster
from contextlib import closing

from sgains.config import Config


class OptionsBase(object):

    def __init__(self):
        self.config = None
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

        parser.add_argument(
            "--parallel", "-p",
            dest="parallel",
            help="number of task to run in parallel",
            type=int,
            default=1
        )

        parser.add_argument(
            "--sge",
            dest="sge",
            action="store_true",
            help="parallelilizes commands using SGE cluster manager",
            default=False
        )


    def create_local_cluster(self):
        workers = self.config.parallel
        threads_per_worker = 1
        print("workers=", workers, " threads_per_worker=", threads_per_worker)
        cluster = LocalCluster(
            n_workers=workers, threads_per_worker=threads_per_worker)
        return cluster

    def create_sge_cluster(self):
        from dask_jobqueue import SGECluster
        
        workers = self.config.parallel
        threads_per_worker = 1
        print("workers=", workers, " threads_per_worker=", threads_per_worker)
        queue = "all.q@wigclust1.cshl.edu,all.q@wigclust3.cshl.edu,all.q@wigclust4.cshl.edu,"\
            "all.q@wigclust5.cshl.edu,all.q@wigclust6.cshl.edu,all.q@wigclust7.cshl.edu,"\
            "all.q@wigclust8.cshl.edu,all.q@wigclust10.cshl.edu,all.q@wigclust11.cshl.edu,"\
            "all.q@wigclust12.cshl.edu,all.q@wigclust13.cshl.edu,all.q@wigclust14.cshl.edu,"\
            "all.q@wigclust15.cshl.edu,all.q@wigclust16.cshl.edu,all.q@wigclust17.cshl.edu,"\
            "all.q@wigclust18.cshl.edu,all.q@wigclust19.cshl.edu"

        cluster = SGECluster(
            queue=queue,
            walltime="02:00:00",
            processes=1,   # we request 10 processes per worker
            memory='4GB',  # for memory requests, this must be specified
            resource_spec='m_mem_free=8G',  # for memory requests, this also needs to be specified
            cores=2)        
        cluster.scale_up(n=workers)
        return cluster

    def create_dask_cluster(self):
        if self.config.sge:
            return self.create_sge_cluster()
        else:
            return self.create_local_cluster()

    def create_dask_client(self, dask_cluster):
        client = Client(dask_cluster)
        return client

    def run_pipeline(self, pipeline):
        dask_cluster = self.create_dask_cluster()
        with closing(dask_cluster) as cluster:
            dask_client = self.create_dask_client(cluster)
            with closing(dask_client) as client:

                pipeline.run(dask_client=client)

    def common_updates(self, args):
        if args.config is not None:
            config = Config.load(args.config)
        else:
            try:
                config = Config.load('./sgains.yml')
            except AssertionError:
                config = Config.default()

        self.config = config

        if args.dry_run is not None:
            self.config.dry_run = args.dry_run
        if args.force is not None:
            self.config.force = args.force
        if args.sge is not None:
            self.config.sge = args.sge

        if args.parallel is not None:
            self.config.parallel = args.parallel


class DataDirMixin(object):
    __slots__ = ()

    def data_dir_options(self, config, glob=False):
        assert self.subparser is not None

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

    def data_dir_update(self, args, config, glob=False):
        assert self.subparser is not None

        if args.data_dir is not None:
            config.data_dir = args.data_dir
        if glob:
            if args.data_glob is not None:
                config.data_glob = args.data_glob


class WorkDirMixin(object):
    __slots__ = ()

    def work_dir_options(self, config):
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "output data options")
        group.add_argument(
            "--work-dir", "-o",
            dest="work_dir",
            help="output directory where results from processing are stored",
            default=config.work_dir
        )
        return group

    def work_dir_update(self, args, config):
        if args.work_dir is not None:
            config.work_dir = args.work_dir


class GenomeIndexMixin(object):
    __slots__ = ()

    def genome_index_options(self, config, input_dir=False):
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "genome index options")
        group.add_argument(
            "--genome-index", "-G",
            dest="genome_index",
            help="genome index name",
            default=config.genome.index)

        group.add_argument(
            "--genome-dir",
            dest="genome_dir",
            help="genome index directory",
            default=config.genome.work_dir)

        group.add_argument(
            "--genome-version",
            dest="genome_version",
            help="version of reference genome to use",
            default=config.genome.version)
        if input_dir:
            group.add_argument(
                "--genome-pristine",
                dest="genome_pristine",
                help="directory where clean copy of reference genome "
                "is located",
                default=config.genome.data_dir)
        return group

    def genome_index_updates(self, args, input_dir=False):
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

    def bins_boundaries_options(self, config, bins_count=False):
        group = self.subparser.add_argument_group(
            "bins boundaries")
        group.add_argument(
            "--bins-boundaries", "-B",
            dest="bins_boundaries",
            help="bins boundaries filename",
            default=config.bins.bins_boundaries
        )
        group.add_argument(
            "--bins-dir",
            dest="bins_dir",
            help="bins working directory",
            default=config.bins.bins_dir
        )
        if bins_count:
            group.add_argument(
                "--bins-count", "-C",
                dest="bins_count",
                type=int,
                help="number of bins",
                default=config.bins.bins_count
            )

    def bins_boundaries_updates(self, args, bins_count=False):
        if args.bins_boundaries:
            self.config.bins.bins_boundaries = args.bins_boundaries
        if args.bins_dir:
            self.config.bins.bins_dir = args.bins_dir
        if bins_count:
            if args.bins_count:
                self.config.bins.bins_count = args.bins_count


class MappableRegionsMixin(object):
    __slots__ = ()

    def mappable_regions_options(self, config, read_length=False):
        group = self.subparser.add_argument_group(
            "mappable regions options")
        group.add_argument(
            "--mappable-dir", "-m",
            dest="mappable_dir",
            help="directory where mappable regions file is stroed",
            default=config.mappable_regions.work_dir
        )
        group.add_argument(
            "--mappable-regions", "-M",
            dest="mappable_regions",
            help="filename where mappable regions are stored",
            default=config.mappable_regions.mappable_regions
        )
        if read_length:
            group.add_argument(
                "--read-length", "-l",
                dest="length",
                type=int,
                help="read length to use for generation of mappable regions",
                default=config.mappable_regions.length
            )
            group.add_argument(
                "--bowtie-opts",
                dest="bowtie_opts",
                help="additional bowtie options",
                default=config.mappable_regions.bowtie_opts
            )
        return group

    def mappable_regions_updates(self, args, read_length=False):
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


class MappingMixin(object):
    __slots__ = ()

    def reads_dir_options(self, config):
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "sequencing reads options")
        group.add_argument(
            "--reads-dir", "-R",
            dest="reads_dir",
            help="data directory where sequencing reads are located",
            default=config.mapping.reads_dir
        )
        group.add_argument(
            "--reads-suffix",
            dest="reads_suffix",
            help="reads files suffix pattern",
            default=config.mapping.reads_suffix)

        return group

    def reads_dir_updates(self, args):
        assert self.subparser is not None

        if args.reads_dir is not None:
            self.config.mapping.reads_dir = args.reads_dir
        if args.reads_suffix is not None:
            self.config.mapping.reads_suffix = args.reads_suffix

    def mapping_dir_options(self, config):
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "mapping files options")
        group.add_argument(
            "--mapping-dir", "-M",
            dest="mapping_dir",
            help="data directory where mapping files are located",
            default=config.mapping.mapping_dir
        )
        group.add_argument(
            "--mapping-suffix",
            dest="mapping_suffix",
            help="mapping files suffix pattern",
            default=config.mapping.mapping_suffix)
        return group

    def mapping_dir_updates(self, args):
        assert self.subparser is not None

        if args.mapping_dir is not None:
            self.config.mapping.mapping_dir = args.mapping_dir
        if args.mapping_suffix is not None:
            self.config.mapping.mapping_suffix = args.mapping_suffix

    def mapping_bowtie_opts(self, config):
        self.subparser.add_argument(
            "--mapping-bowtie-opts",
            dest="mapping_bowtie_opts",
            help="bowtie mapping options",
            default=config.mapping.mapping_bowtie_opts
        )

    def mapping_bowtie_updates(self, args):
        if args.mapping_bowtie_opts:
            self.config.mapping.mapping_bowtie_opts = \
                args.mapping_bowtie_opts

    def mapping_options(self, config):
        self.reads_dir_options(config=config)
        self.mapping_dir_options(config=config)
        self.mapping_bowtie_opts(config=config)

    def mapping_updates(self, args):
        self.reads_dir_updates(args)
        self.mapping_dir_updates(args)
        self.mapping_bowtie_updates(args)


class SCclustMixin(object):
    __slots__ = ()

    def scclust_options(self, config):
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "SCclust options")
        group.add_argument(
            "--scclust-dir", "-S",
            dest="scclust_dir",
            help="SCGV directory",
            default=config.scclust.scclust_dir
        )
        group.add_argument(
            "--case-name",
            dest="case_name",
            help="case name",
            default=config.scclust.case_name)

        return group

    def scclust_updates(self, args):
        assert self.subparser is not None

        if args.scclust_dir is not None:
            self.config.scclust.scclust_dir = args.scclust_dir
        if args.case_name is not None:
            self.config.scclust.case_name = args.case_name
