'''
Created on Aug 3, 2017

@author: lubo
'''
from commands.common import OptionsBase, \
    GenomeIndexMixin, BinsBoundariesMixin, MappingMixin

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
        OptionsBase):

    def output_options(self, config):
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "process output options")
        group.add_argument(
            "--output-dir", "-o",
            dest="output_dir",
            help="output directory",
            default=config.scclust.scclust_dir
        )
        group.add_argument(
            "--case-name",
            dest="case_name",
            help="case name",
            default=config.scclust.case_name)

        return group

    def output_updates(self, args):
        assert self.subparser is not None

        if args.output_dir is not None:
            self.config.scclust.scclust_dir = args.output_dir
        if args.case_name is not None:
            self.config.scclust.case_name = args.case_name

    def __init__(self, parser, subparsers):
        super(ProcessCommand, self).__init__()
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="process",
            help="combines mapping, varbin and scclust subcommands into "
            "single command",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser.set_defaults(func=self.run)

    def add_options(self, config):
        self.reads_dir_options(config=config)
        self.mapping_bowtie_opts(config=config)
        self.output_options(config=config)

        self.genome_index_options(config=config, input_dir=False)
        self.bins_boundaries_options(config=config, bins_count=False)

    def process_args(self, args):
        self.common_updates(args)
        self.reads_dir_updates(args)
        self.mapping_bowtie_updates(args)

        self.output_updates(args)
        self.genome_index_updates(args, input_dir=False)
        self.bins_boundaries_updates(args, bins_count=False)

    def run(self, args):
        print(colored(
            "process subcommand called with args: {}".format(args),
            "yellow"))
        self.process_args(args)

        mapping_workdir = os.path.join(
            self.config.scclust_dirname(),
            'mapping')
        varbin_workdir = os.path.join(
            self.config.scclust_dirname(),
            'varbin')
        scclust_workdir = os.path.join(
            self.config.scclust_dirname(),
            'scclust')

        mapping_config = Config.copy(self.config)
        mapping_config.mapping.mapping_dir = mapping_workdir

        varbin_config = Config.copy(self.config)
        varbin_config.mapping.mapping_dir = mapping_workdir
        varbin_config.varbin.varbin_dir = varbin_workdir

        segment_config = Config.copy(self.config)
        segment_config.varbin.varbin_dir = varbin_workdir
        segment_config.scclust.scclust_dir = scclust_workdir

        if not os.path.exists(mapping_workdir):
            os.makedirs(mapping_workdir)
        if not os.path.exists(varbin_workdir):
            os.makedirs(varbin_workdir)
        if not os.path.exists(scclust_workdir):
            os.makedirs(scclust_workdir)

        pipeline = MappingPipeline(mapping_config)
        # pipeline.run()

        pipeline = VarbinPipeline(varbin_config)
        # pipeline.run()

        pipeline = Rpipeline(segment_config)
        # pipeline.run()
