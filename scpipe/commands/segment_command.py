'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse

from commands.common import OptionsBase, DataDirMixin,\
    WorkDirMixin, BinsBoundariesMixin
from pipelines.r_pipeline import Rpipeline


class SegmentMixin(DataDirMixin, WorkDirMixin, BinsBoundariesMixin):

    def segment_options(self, config):
        self.data_dir_options(config=config.segment, glob=True)
        group = self.work_dir_options(config=config.segment)
        group.add_argument(
            "--study-name", "-s",
            help="study name",
            dest="study_name",
            default=self.config.segment.study_name
        )
        self.bins_boundaries_options(bins_count=False)

    def segment_updates(self, args):
        self.common_updates(args)
        self.work_dir_update(args, config=self.config.segment)
        self.data_dir_update(args, config=self.config.segment, glob=True)
        self.bins_boundaries_updates(args, bins_count=False)
        if args.study_name is not None:
            self.config.segment.study_name = args.study_name


class SegmentCommand(
        SegmentMixin,
        OptionsBase):

    def __init__(self, config, parser, subparsers):
        super(SegmentCommand, self).__init__(config)
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="segment",
            help="segments bin counts and prepares the SCGV input data",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser.set_defaults(func=self.run)

    def add_options(self, config):
        self.segment_options(config)

    def process_args(self, args):
        self.common_updates(args)
        self.segment_updates(args)

    def run(self, args):
        print("segment subcommand called with args: {}".format(args))
        self.process_args(args)
        pipeline = Rpipeline(self.config)
        pipeline.run()
