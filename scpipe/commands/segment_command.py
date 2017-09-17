'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse

from commands.common import OptionsBase, BinsBoundariesMixin, SegmentMixin
from pipelines.r_pipeline import Rpipeline
from commands.varbin_command import VarbinMixin


class SegmentCommand(
        VarbinMixin,
        SegmentMixin,
        BinsBoundariesMixin,
        OptionsBase):

    def __init__(self, parser, subparsers):
        super(SegmentCommand, self).__init__()
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="segment",
            help="segments bin counts and prepares the SCGV input data",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser.set_defaults(func=self.run)

    def add_options(self, config):
        self.varbin_dir_options(config)
        self.segment_options(config)
        self.bins_boundaries_options(config, bins_count=False)

    def process_args(self, args):
        self.common_updates(args)
        self.varbin_dir_updates(args)
        self.segment_updates(args)
        self.bins_boundaries_updates(args, bins_count=False)

    def run(self, args):
        print("segment subcommand called with args: {}".format(args))
        self.process_args(args)
        pipeline = Rpipeline(self.config)
        pipeline.run()
