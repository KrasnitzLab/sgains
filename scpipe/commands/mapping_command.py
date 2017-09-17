'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse

from commands.common import GenomeIndexMixin, OptionsBase, MappingMixin
from pipelines.mapping_pipeline import MappingPipeline
from termcolor import colored


class MappingCommand(
        MappingMixin,
        GenomeIndexMixin,
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
        self.genome_index_options(config, input_dir=False)

    def process_args(self, args):
        self.common_updates(args)
        self.mapping_updates(args)
        self.genome_index_updates(args, input_dir=False)

    def run(self, args):
        print(colored(
            "mapping subcommand called with args: {}".format(args),
            "yellow"))
        self.process_args(args)
        pipeline = MappingPipeline(self.config)
        pipeline.run()
