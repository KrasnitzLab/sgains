'''
Created on Aug 3, 2017

@author: lubo
'''
from commands.common import OptionsBase, \
    GenomeIndexMixin, BinsBoundariesMixin, MappingMixin, SegmentMixin

import argparse
import os
from config import Config
from pipelines.mapping_pipeline import MappingPipeline
from pipelines.varbin_pipeline import VarbinPipeline
from pipelines.r_pipeline import Rpipeline
from termcolor import colored


class ProcessCommand(
        GenomeIndexMixin,
        BinsBoundariesMixin,
        MappingMixin,
        SegmentMixin,
        OptionsBase):

    def __init__(self, parser, subparsers):
        super(ProcessCommand, self).__init__()
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="process",
            help="combines mapping, varbin and segment subcommands into "
            "single command",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser.set_defaults(func=self.run)

    def add_options(self, config):
        self.reads_dir_options(config=config)
        self.mapping_bowtie_opts(config=config)
        self.segment_options(config=config)

        self.genome_index_options(config=config, input_dir=False)
        self.bins_boundaries_options(config=config, bins_count=False)

    def process_args(self, args):
        self.common_updates(args)
        self.reads_dir_updates(args)
        self.mapping_bowtie_updates(args)
        self.segment_updates(args)
        self.genome_index_updates(args, input_dir=False)
        self.bins_boundaries_updates(args, bins_count=False)

    def run(self, args):
        print(colored(
            "process subcommand called with args: {}".format(args),
            "yellow"))
        self.process_args(args)

        mapping_workdir = os.path.join(
            self.config.segment_dirname(),
            'mapping')
        varbin_workdir = os.path.join(
            self.config.segment_dirname(),
            'varbin')
        segment_workdir = os.path.join(
            self.config.segment_dirname(),
            'segment')

        mapping_config = Config.copy(self.config)
        mapping_config.mapping.mapping_dir = mapping_workdir

        varbin_config = Config.copy(self.config)
        varbin_config.mapping.mapping_dir = mapping_workdir
        varbin_config.varbin.varbin_dir = varbin_workdir

        segment_config = Config.copy(self.config)
        segment_config.varbin.varbin_dir = varbin_workdir
        segment_config.segment.segment_dir = segment_workdir

        pipeline = MappingPipeline(mapping_config)
        pipeline.run()

        pipeline = VarbinPipeline(varbin_config)
        pipeline.run()

        pipeline = Rpipeline(segment_config)
        pipeline.run()
