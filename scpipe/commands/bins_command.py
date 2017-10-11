'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse

from pipelines.bins_pipeline import BinsPipeline
from commands.common import OptionsBase,\
    MappableRegionsMixin, BinsBoundariesMixin, GenomeIndexMixin


class BinsCommand(
        GenomeIndexMixin,
        MappableRegionsMixin,
        BinsBoundariesMixin,
        OptionsBase):

    def __init__(self, parser, subparsers):
        super(BinsCommand, self).__init__()
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="bins",
            help="calculates all bins boundaries for specified bins count "
            "and read length",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser.set_defaults(func=self.run)

    def add_options(self, config):
        self.mappable_regions_options(config=config, read_length=False)
        self.bins_boundaries_options(config=config, bins_count=True)
        self.genome_index_options(config=config, input_dir=False)

    def process_args(self, args):
        self.common_updates(args)
        self.mappable_regions_updates(args, read_length=False)
        self.bins_boundaries_updates(args, bins_count=True)
        self.genome_index_updates(args, input_dir=False)

    def run(self, args):
        print("bins subcommand called with args: {}".format(args))
        self.process_args(args)
        pipeline = BinsPipeline(self.config)
        pipeline.run()
