'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse
from commands.common import OptionsBase,\
    GenomeIndexMixin
from genomeindex_pipeline import GenomeIndexPipeline


class GenomeIndexCommand(
        GenomeIndexMixin,
        OptionsBase):

    def __init__(self, config, parser, subparsers):
        super(GenomeIndexCommand, self).__init__(config)
        self.subconfig = config.genome
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="genomeindex",
            help="build appropriate genome index",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.parser.set_defaults(func=self.run)

    def add_options(self):
        self.genome_index_options(input_dir=True)

    def process_args(self, args):
        self.common_updates(args)
        self.genome_index_update(args, input_dir=True)

        return self.config

    def run(self, args):
        print("genomeindex subcommand runned with args: {}".format(args))
        self.process_args(args)
        pipeline = GenomeIndexPipeline(self.config)
        pipeline.run()
