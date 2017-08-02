'''
Created on Aug 2, 2017

@author: lubo
'''
from commands.common import WorkDirMixin, OptionsBase, GenomeIndexMixin
import argparse
from mappableregions_pipeline import MappableRegionsPipeline


class MappableRegionsCommand(
        WorkDirMixin,
        GenomeIndexMixin,
        OptionsBase):

    def __init__(self, config, parser, subparsers):
        super(MappableRegionsCommand, self).__init__(config)
        self.subconfig = config.genome
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="mappable-regions",
            help="finds all mappable regions in specified genome",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.parser.set_defaults(func=self.run)

    def mappable_regions_options(self):
        group = self.subparser.add_argument_group(
            "mappable regions options")
        group.add_argument(
            "--work-dir", "-o",
            dest="work_dir",
            help="output directory where results from processing are stored",
            default=self.config.mappable_regions.work_dir
        )
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

    def mappable_regions_update(self, args):
        if args.bowtie_opts:
            self.config.mappable_regions.bowtie_opts = args.bowtie_opts
        if args.length:
            self.config.mappable_regions.length = args.length
        if args.work_dir is not None:
            self.config.mappable_regions.work_dir = args.work_dir

    def add_options(self):
        self.mappable_regions_options()
        self.genome_index_options(genome_dir=True)

    def process_args(self, args):
        self.common_updates(args)
        self.genome_index_update(args, genome_dir=True)
        self.mappable_regions_update(args)

    def run(self, args):
        print("mappable-regions subcommand called with args: {}".format(args))
        self.process_args(args)
        pipeline = MappableRegionsPipeline(self.config)
        pipeline.run()
