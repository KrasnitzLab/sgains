'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse
from sgains.commands.common import CommandBase
from sgains.pipelines.genomeindex_pipeline import GenomeIndexPipeline


class GenomeIndexCommand(CommandBase):

    def __init__(self, config):
        super(GenomeIndexCommand, self).__init__(config)

    def run(self, args):
        print("genomeindex subcommand runned with args: {}".format(args))
        pipeline = GenomeIndexPipeline(self.config)
        pipeline.run()
