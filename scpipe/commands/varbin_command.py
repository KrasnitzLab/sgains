'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse

from commands.common import OptionsBase, DataDirMixin,\
    WorkDirMixin, BinsBoundariesMixin
from pipelines.varbin_pipeline import VarbinPipeline


class VarbinMixin(DataDirMixin, WorkDirMixin, BinsBoundariesMixin):

    def mapping_options(self, config):
        self.data_dir_options(config=config.varbin, glob=True)
        group = self.work_dir_options(config=config.varbin)
        group.add_argument(
            "--suffix", "-s",
            help="suffix for output files",
            dest="suffix",
            default=self.config.varbin.suffix
        )
        self.bins_boundaries_options(bins_count=False)

    def mapping_updates(self, args):
        self.common_updates(args)
        self.work_dir_update(args, config=self.config.varbin)
        self.data_dir_update(args, config=self.config.varbin, glob=True)
        self.bins_boundaries_updates(args, bins_count=False)
        if args.suffix is not None:
            self.config.varbin.suffix = args.suffix


class VarbinCommand(
        VarbinMixin,
        OptionsBase):

    def __init__(self, config, parser, subparsers):
        super(VarbinCommand, self).__init__(config)
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="varbin",
            help="applies varbin algorithm to count read mappings in each bin",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser.set_defaults(func=self.run)

    def add_options(self, config):
        self.mapping_options(config)

    def process_args(self, args):
        self.common_updates(args)
        self.mapping_updates(args)

    def run(self, args):
        print("varbin subcommand called with args: {}".format(args))
        self.process_args(args)
        pipeline = VarbinPipeline(self.config)
        pipeline.run()
