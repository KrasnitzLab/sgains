'''
Created on Jun 12, 2017

@author: lubo
'''
import os
from config import Config
from hg19 import HumanGenome19
from box import Box


class Parser:
    DEFAULT_CONFIG = {
        "genome": {
            "version": "hg19",
            "pristine": "hg19_safe",
            "cache_dir": "data/safe",
            "index": "genomeindex",
        },
        "reads": {
            "length": 100,
            "cache_dir": "data/R100",
            "mappable_regions": "mappable_regions.tsv",
            "mappable_positions_count": "mappable_positions_count.yml",
            "chrom_sizes": "chrom_sizes.yml"
        },
        "bins": {
            "bins_count": 10000,
            "cache_dir": "data/R100_B10k",
            "bin_boundaries": "bin_boundaries.tst",
        }
    }

    def __init__(self, parser):
        self.parser = parser

    @staticmethod
    def from_argument_parser(parser):
        parser.add_argument(
            "-v", "--verbose",
            dest="verbose",
            action="count",
            help="set verbosity level [default: %(default)s]")

        parser.add_argument(
            "-c", "--config",
            dest="config",
            help="configuration file",
            metavar="path")

        parser.add_argument(
            "-C", "--chroms",
            dest="chroms",
            help="list of chromosomes to work with. Default: all chromosomes",
        )

        parser.add_argument(
            "-g", "--genome",
            dest="genome",
            help="genome version. Default is 'hg19'",
            metavar="path")

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

        result = Parser(parser)
        return result

    def parse_arguments(self, argv):
        args = self.parser.parse_args(argv)
        if args.config:
            assert os.path.exists(args.config)
            config = Config.load(args.config)
        else:
            config = Box(self.DEFAULT_CONFIG, default_box=True)
            config.filename = None
            config.dirname = os.getcwd()

        if args.genome:
            config.genome.version = args.genome
        if args.genome_dir:
            config.genome.cache_dir = args.genome_dir
        if args.genome_index:
            config.genome.index = args.genome_index

        if args.chroms is None:
            config.chroms = HumanGenome19.CHROMS
        else:
            chroms = [c.strip() for c in args.chroms.split(',')]
            config.chroms = chroms

        result = Box(config.to_dict(), frozen_box=True)
        return result
