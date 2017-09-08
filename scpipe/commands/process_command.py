'''
Created on Aug 3, 2017

@author: lubo
'''
from commands.common import OptionsBase, DataDirMixin, WorkDirMixin,\
    GenomeIndexMixin, BinsBoundariesMixin
import argparse
import os
from config import Config
from pipelines.mapping_pipeline import MappingPipeline
from pipelines.varbin_pipeline import VarbinPipeline
from pipelines.r_pipeline import Rpipeline


class ProcessCommand(
        GenomeIndexMixin,
        BinsBoundariesMixin,
        DataDirMixin,
        WorkDirMixin,
        OptionsBase):

    def __init__(self, config, parser, subparsers):
        super(ProcessCommand, self).__init__(config)
        self.subconfig = config.varbin
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="process",
            help="combines mapping, varbin and segment subcommands into "
            "single command",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser.set_defaults(func=self.run)

    def add_options(self):
        self.data_dir_options(config=self.config.mapping, glob=True)
        group = self.work_dir_options(config=self.config.segment)
        group.add_argument(
            "--study-name", "-s",
            help="study name",
            dest="study_name",
            default=self.config.segment.study_name
        )
        group = self.genome_index_options(input_dir=False)
        group.add_argument(
            "--bowtie-opts",
            dest="bowtie_opts",
            help="additional bowtie options",
            default=self.config.mapping.bowtie_opts
        )
        self.bins_boundaries_options(bins_count=False)

    def process_args(self, args):
        self.common_updates(args)
        self.data_dir_update(args, config=self.config.mapping, glob=True)
        self.work_dir_update(args, config=self.config.segment)
        self.genome_index_update(args, input_dir=False)
        self.bins_boundaries_updates(args, bins_count=False)
        if args.study_name is not None:
            self.config.segment.study_name = args.study_name
        if args.bowtie_opts:
            self.config.mapping.bowtie_opts = args.bowtie_opts

    def run(self, args):
        print("process subcommand called with args: {}".format(args))
        self.process_args(args)

        mapping_workdir = os.path.join(
            self.config.segment_work_dirname(),
            'mapping')
        varbin_workdir = os.path.join(
            self.config.segment_work_dirname(),
            'varbin')
        segment_workdir = os.path.join(
            self.config.segment_work_dirname(),
            'segment')

        mapping_config = Config.copy(self.config)
        mapping_config.mapping.work_dir = mapping_workdir

        varbin_config = Config.copy(self.config)
        varbin_config.varbin.data_dir = mapping_workdir
        varbin_config.varbin.work_dir = varbin_workdir

        segment_config = Config.copy(self.config)
        segment_config.segment.data_dir = varbin_workdir
        segment_config.segment.work_dir = segment_workdir

        pipeline = MappingPipeline(mapping_config)
        pipeline.run()

        pipeline = VarbinPipeline(varbin_config)
        pipeline.run()

        pipeline = Rpipeline(segment_config)
        pipeline.run()
