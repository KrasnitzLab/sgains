'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse

from sgains.commands.common import OptionsBase, GenomeIndexMixin,\
    MappableRegionsMixin
from sgains.pipelines.mappableregions_pipeline import MappableRegionsPipeline


class MappableRegionsCommand(
        GenomeIndexMixin,
        MappableRegionsMixin,
        OptionsBase):

    def __init__(self, parser, subparsers):
        super(MappableRegionsCommand, self).__init__()
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="mappable-regions",
            help="finds all mappable regions in specified genome",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser.set_defaults(func=self.run)

    def add_options(self, config):
        self.mappable_regions_options(config=config, read_length=True)
        self.genome_index_options(config=config, input_dir=False)

    def process_args(self, args):
        self.common_updates(args)
        self.genome_index_updates(args, input_dir=False)
        self.mappable_regions_updates(args, read_length=True)

    def run(self, args):
        print("mappable-regions subcommand called with args: {}".format(args))
        self.process_args(args)
        pipeline = MappableRegionsPipeline(self.config)
        pipeline.run()
