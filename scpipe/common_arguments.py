'''
Created on Jun 12, 2017

@author: lubo
'''
import os
from config import Config
from hg19 import HumanGenome19
from box import Box


class Parser:

    def __init__(self, parser):
        self.parser = parser

    @staticmethod
    def from_argument_parser(parser):
        parser.add_argument(
            "-v", "--verbose",
            dest="verbose",
            action="count",
            help="set verbosity level [default: %(default)s]"
        )

        parser.add_argument(
            "-c", "--config",
            dest="config",
            help="configuration file",
            metavar="path"
        )

        parser.add_argument(
            "-p", "--threads",
            dest="threads",
            help="number of threads to use for bowtie",
            type=int
        )

        parser.add_argument(
            "-o", "--output",
            dest="output",
            help="output file to store results into. "
            "To use standart output you need to specify '-o -'",
            metavar="path"
        )

        parser.add_argument(
            "-C", "--chroms",
            dest="chroms",
            help="list of chromosomes to work with. Default: all chromosomes",
        )

        parser.add_argument(
            "-g", "--genome",
            dest="genome",
            help="genome version. Default is 'hg19'",
            metavar="path"
        )

        parser.add_argument(
            "--genome-dir",
            dest="genome_dir",
            help="directory where genome and genome index are located",
            metavar="path"
        )

        parser.add_argument(
            "--genome-index",
            dest="genome_index",
            help="genome index name",
        )

        parser.add_argument(
            "-l", "--length",
            dest="length",
            help="read lengths to use for region mappings. Default: 100",
            type=int,
        )

        parser.add_argument(
            "--reads-dir",
            dest="reads_dir",
            help="directory where mappable regions are stored",
            metavar="path"
        )

        parser.add_argument(
            "-b", "--bins",
            dest="bins",
            type=int,
            help="number of bins to use: Default: 100"
        )

        parser.add_argument(
            "--bins-dir",
            dest="bins_dir",
            help="directory where bin boundaries are stored",
            metavar="path"
        )

        result = Parser(parser)
        return result

    def parse_arguments(self, argv):
        self.args = self.parser.parse_args(argv)
        args = self.args

        if args.config:
            assert os.path.exists(args.config)
            config = Config.load(args.config)
        else:
            config = Config.default()

        if args.output:
            if os.path.isabs(args.output):
                config.output = args.output
            elif '-' == args.output:
                config.output = '-'
            else:
                output = os.path.join(
                    config.dirname,
                    args.output
                )
                config.output = os.path.abspath(output)
        else:
            config.output = None

        if args.genome:
            config.genome.version = args.genome
        if args.genome_dir:
            config.genome.work_dir = args.genome_dir
        if args.genome_index:
            config.genome.index = args.genome_index

        if args.chroms is None:
            config.chroms = HumanGenome19.CHROMS
        else:
            chroms = [c.strip() for c in args.chroms.split(',')]
            config.chroms = chroms

        if args.length:
            config.reads.length = args.length
        if args.reads_dir:
            config.reads.work_dir = args.reads_dir

        if args.bins:
            config.bins.bins_count = args.bins
        if args.bins_dir:
            config.bins.work_dir = args.bins_dir

        if args.threads:
            config.threads = args.threads
        else:
            config.threads = 1

        result = Config(Box(config.to_dict(), frozen_box=True))
        reads_cache = result.abspath(result.reads.work_dir)
        if not os.path.exists(reads_cache):
            os.makedirs(reads_cache)

        bins_cache = result.abspath(result.bins.work_dir)
        if not os.path.exists(bins_cache):
            os.makedirs(bins_cache)

        return result
