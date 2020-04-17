'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse

from sgains.commands.common import CommandBase
from sgains.pipelines.mappableregions_pipeline import MappableRegionsPipeline


class MappableRegionsCommand(CommandBase):

    def __init__(self, config):
        super(MappableRegionsCommand, self).__init__(config)

    def run(self, args):
        print("mappable-regions subcommand called with args: {}".format(args))
        pipeline = MappableRegionsPipeline(self.config)
        self.run_pipeline(pipeline)
