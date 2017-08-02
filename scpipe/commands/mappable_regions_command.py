'''
Created on Aug 2, 2017

@author: lubo
'''
from commands.common import WorkDirMixin, OptionsBase, GenomeIndexMixin,\
    MappableRegionsMixin
import argparse
from mappableregions_pipeline import MappableRegionsPipeline


class MappableRegionsCommand(
        WorkDirMixin,
        GenomeIndexMixin,
        MappableRegionsMixin,
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

    def add_options(self):
        self.mappable_regions_options(read_length=True)
        self.genome_index_options(genome_dir=True)

    def process_args(self, args):
        self.common_updates(args)
        self.genome_index_update(args, genome_dir=True)
        self.mappable_regions_update(args, read_length=True)

    def run(self, args):
        print("mappable-regions subcommand called with args: {}".format(args))
        self.process_args(args)
        pipeline = MappableRegionsPipeline(self.config)
        pipeline.run()
