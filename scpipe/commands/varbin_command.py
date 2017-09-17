'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse

from commands.common import OptionsBase, BinsBoundariesMixin
from pipelines.varbin_pipeline import VarbinPipeline
from commands.mapping_command import MappingMixin


class VarbinMixin(object):

    def varbin_dir_options(self, config):
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "varbin options")
        group.add_argument(
            "--varbin-dir", "-V",
            dest="varbin_dir",
            help="varbin directory",
            default=config.varbin.varbin_dir
        )
        group.add_argument(
            "--varbin-suffix",
            dest="varbin_suffix",
            help="varbin files suffix pattern",
            default=config.varbin.varbin_suffix)

        return group

    def varbin_dir_updates(self, args):
        assert self.subparser is not None

        if args.varbin_dir is not None:
            self.config.varbin.varbin_dir = args.varbin_dir
        if args.varbin_suffix is not None:
            self.config.varbin.varbin_suffix = args.varbin_suffix


class VarbinCommand(
        VarbinMixin,
        BinsBoundariesMixin,
        MappingMixin,
        OptionsBase):

    def __init__(self, parser, subparsers):
        super(VarbinCommand, self).__init__()
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="varbin",
            help="applies varbin algorithm to count read mappings in each bin",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser.set_defaults(func=self.run)

    def add_options(self, config):
        self.mapping_dir_options(config)
        self.varbin_dir_options(config)
        self.bins_boundaries_options(config, bins_count=False)

    def process_args(self, args):
        self.common_updates(args)
        self.mapping_dir_updates(args)
        self.varbin_dir_updates(args)
        self.bins_boundaries_updates(args, bins_count=False)

    def run(self, args):
        print("varbin subcommand called with args: {}".format(args))
        self.process_args(args)
        pipeline = VarbinPipeline(self.config)
        pipeline.run()
