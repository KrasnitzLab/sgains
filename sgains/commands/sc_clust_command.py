'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse
from termcolor import colored

from sgains.commands.common import OptionsBase, BinsBoundariesMixin, SCclustMixin
from sgains.pipelines.r_pipeline import Rpipeline
from sgains.commands.varbin_command import VarbinMixin


class SCclustCommand(
        VarbinMixin,
        SCclustMixin,
        BinsBoundariesMixin,
        OptionsBase):

    def __init__(self, parser, subparsers):
        super(SCclustCommand, self).__init__()
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="scclust",
            help="segmentation and clustering based bin counts and "
            "preparation of the SCGV input data",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser.set_defaults(func=self.run)

    def add_options(self, config):
        self.varbin_dir_options(config)
        self.scclust_options(config)
        self.bins_boundaries_options(config, bins_count=False)

    def process_args(self, args):
        self.common_updates(args)
        self.varbin_dir_updates(args)
        self.scclust_updates(args)
        self.bins_boundaries_updates(args, bins_count=False)

    def run(self, args):
        print(colored(
            "scclust subcommand called with args: {}".format(args),
            "yellow"))
        self.process_args(args)
        pipeline = Rpipeline(self.config)
        pipeline.run()
