'''
Created on Aug 2, 2017

@author: lubo
'''
from commands.common import GenomeIndexMixin, OptionsBase, \
    MappableRegionsMixin,\
    BinsBoundariesMixin
import argparse


class PrepareCommand(
        GenomeIndexMixin,
        MappableRegionsMixin,
        BinsBoundariesMixin,
        OptionsBase):

    def __init__(self, config, parser, subparsers):
        super(PrepareCommand, self).__init__(config)
        self.subconfig = config.genome
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="prepare",
            help="combines all preparation steps "
            "(genomeindex, mappable_regions, bins) into single command",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.parser.set_defaults(func=self.run)

    def add_options(self):
        self.genome_index_options(input_dir=True)
        self.mappable_regions_options(read_length=True)
        self.bins_boundaries_options(bins_count=True)

    def process_args(self, args):
        self.common_updates(args)
        self.genome_index_update(args, input_dir=True)
        self.mappable_regions_update(args, read_length=True)
        self.bins_boundaries_updates(args, bins_count=True)

        return self.config

    def run(self, args):
        print("prepare subcommand runned with args: {}".format(args))
        self.process_args(args)
