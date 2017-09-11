'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse

from commands.common import GenomeIndexMixin, OptionsBase, DataDirMixin,\
    WorkDirMixin
from pipelines.mapping_pipeline import MappingPipeline


class MappingMixin(DataDirMixin, WorkDirMixin, GenomeIndexMixin):

    def mapping_options(self, config):
        self.genome_index_options(config=config, input_dir=False)
        group = self.data_dir_options(config=config.mapping, glob=True)
        group.add_argument(
            "--bowtie-opts",
            dest="bowtie_opts",
            help="additional bowtie options",
            default=config.mapping.bowtie_opts
        )
        self.work_dir_options(config=config.mapping)

    def mapping_updates(self, args):
        self.common_updates(args)
        self.work_dir_update(args, config=self.config.mapping)
        self.data_dir_update(args, config=self.config.mapping, glob=True)
        self.genome_index_update(args, input_dir=False)
        if args.bowtie_opts:
            self.config.mapping.bowtie_opts = args.bowtie_opts


class MappingCommand(
        MappingMixin,
        OptionsBase):

    def __init__(self, parser, subparsers):
        super(MappingCommand, self).__init__()
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="mapping",
            help="performs mapping of cell reads to reference genome",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser.set_defaults(func=self.run)

    def add_options(self, config):
        self.mapping_options(config)

    def process_args(self, args):
        self.common_updates(args)
        self.mapping_updates(args)

    def run(self, args):
        print("mapping subcommand called with args: {}".format(args))
        self.process_args(args)
        pipeline = MappingPipeline(self.config)
        pipeline.run()
