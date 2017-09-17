'''
Created on Aug 2, 2017

@author: lubo
'''
from commands.common import GenomeIndexMixin, OptionsBase, \
    MappableRegionsMixin,\
    BinsBoundariesMixin
import argparse
from pipelines.genomeindex_pipeline import GenomeIndexPipeline
from pipelines.mappableregions_pipeline import MappableRegionsPipeline
from pipelines.bins_pipeline import BinsPipeline


class PrepareCommand(
        GenomeIndexMixin,
        MappableRegionsMixin,
        BinsBoundariesMixin,
        OptionsBase):

    def __init__(self, parser, subparsers):
        super(PrepareCommand, self).__init__()
        self.parser = parser
        self.subparser = subparsers.add_parser(
            name="prepare",
            help="combines all preparation steps "
            "(genomeindex, mappable_regions, bins) into single command",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        self.subparser.set_defaults(func=self.run)

    def add_options(self, config):
        self.genome_index_options(config=config, input_dir=True)
        self.mappable_regions_options(config=config, read_length=True)
        self.bins_boundaries_options(config=config, bins_count=True)

    def process_args(self, args):
        self.common_updates(args)
        self.genome_index_updates(args, input_dir=True)
        self.mappable_regions_updates(args, read_length=True)
        self.bins_boundaries_updates(args, bins_count=True)

        return self.config

    def run(self, args):
        print("prepare subcommand runned with args: {}".format(args))
        self.process_args(args)

        pipeline = GenomeIndexPipeline(self.config)
        pipeline.run()

        pipeline = MappableRegionsPipeline(self.config)
        pipeline.run()

        pipeline = BinsPipeline(self.config)
        pipeline.run()
