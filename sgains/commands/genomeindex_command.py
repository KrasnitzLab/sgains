'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse
from sgains.commands.common import OptionsBase,\
    GenomeIndexMixin
from sgains.pipelines.genomeindex_pipeline import GenomeIndexPipeline


class GenomeIndexCommand(
        GenomeIndexMixin,
        OptionsBase):

    def __init__(self, parser, subparsers):
        super(GenomeIndexCommand, self).__init__()
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="genomeindex",
            help="builds appropriate bowtie index for the reference genome",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser.set_defaults(func=self.run)

    def add_options(self, config):
        self.genome_index_options(config=config, input_dir=True)

    def process_args(self, args):
        self.common_updates(args)
        self.genome_index_updates(args, input_dir=True)

        return self.config

    def run(self, args):
        print("genomeindex subcommand runned with args: {}".format(args))
        self.process_args(args)
        pipeline = GenomeIndexPipeline(self.config)
        pipeline.run()
