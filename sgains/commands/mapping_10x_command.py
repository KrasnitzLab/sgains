'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse

from termcolor import colored

from sgains.commands.common import GenomeIndexMixin, OptionsBase, \
     Mapping10xMixin
from sgains.pipelines.mapping_10x_pipeline import Mapping10xPipeline


class Mapping10xCommand(
        Mapping10xMixin,
        GenomeIndexMixin,
        OptionsBase):

    def __init__(self, parser, subparsers):
        super(Mapping10xCommand, self).__init__()
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="mapping-10x",
            help="performs mapping of reads from 10xGenomics dataset"
            " to reference genome",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser.set_defaults(func=self.run)

    def add_options(self, config):
        self.mapping_10x_options(config)
        self.genome_index_options(config, input_dir=False)

    def process_args(self, args):
        self.common_updates(args)
        self.mapping_10x_updates(args)
        self.genome_index_updates(args, input_dir=False)

    def run(self, args):
        print(colored(
            "mapping_10x subcommand called with args: {}".format(args),
            "yellow"))
        self.process_args(args)
        pipeline = Mapping10xPipeline(self.config)
        self.run_pipeline(pipeline)
