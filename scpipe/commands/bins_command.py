'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse

from bins_pipeline import BinsPipeline
from commands.common import OptionsBase,\
    MappableRegionsMixin, BinsBoundariesMixin


class BinsCommand(
        MappableRegionsMixin,
        BinsBoundariesMixin,
        OptionsBase):

    def __init__(self, config, parser, subparsers):
        super(BinsCommand, self).__init__(config)
        self.subconfig = config.genome
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="bins",
            help="calculates all bins boundaries for specified bins count "
            "and read length",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.parser.set_defaults(func=self.run)

    def add_options(self):
        self.mappable_regions_options(read_length=False)
        self.bins_boundaries_options(bins_count=True)

    def process_args(self, args):
        self.common_updates(args)
        self.mappable_regions_update(args, read_length=False)
        self.bins_boundaries_updates(args, bins_count=True)

    def run(self, args):
        print("bins subcommand called with args: {}".format(args))
        self.process_args(args)
        pipeline = BinsPipeline(self.config)
        pipeline.run()
