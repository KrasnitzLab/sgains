'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse

from pipelines.bins_pipeline import BinsPipeline
from commands.common import OptionsBase,\
    MappableRegionsMixin, BinsBoundariesMixin


class BinsCommand(
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

    def process_args(self, args):
        self.common_updates(args)
        self.mappable_regions_updates(args, read_length=False)
        self.bins_boundaries_updates(args, bins_count=True)

    def run(self, args):
        print("bins subcommand called with args: {}".format(args))
        self.process_args(args)
        pipeline = BinsPipeline(self.config)
        pipeline.run()
