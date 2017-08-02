'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse
from commands.common import OptionsBase, WorkDirMixin, DataDirMixin,\
    GenomeIndexMixin
from genomeindex_pipeline import GenomeIndexPipeline


class GenomeIndexCommand(
        WorkDirMixin,
        DataDirMixin,
        GenomeIndexMixin,
        OptionsBase):

    def __init__(self, config, parser, subparsers):
        super(GenomeIndexCommand, self).__init__(config)
        self.config = config
        self.subconfig = config.genome
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="genomeindex",
            help="build appropriate genome index",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.parser.set_defaults(func=self.run)

    def add_options(self):
        self.data_dir_options(glob=False)
        self.work_dir_options()
        self.genome_index_options(genome_dir=False)

    def process_args(self, args):
        self.common_updates(args)
        self.data_dir_update(args, glob=False)
        self.work_dir_update(args)
        self.genome_index_update(args, genome_dir=False)

        return self.config

    def run(self, args):
        print("genomeindex subcommand runned with args: {}".format(args))
        self.process_args(args)
        pipeline = GenomeIndexPipeline(self.config)
        pipeline.run()
