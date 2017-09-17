'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse

from commands.common import GenomeIndexMixin, OptionsBase
from pipelines.mapping_pipeline import MappingPipeline


class MappingMixin(GenomeIndexMixin):

    def reads_dir_options(self, config):
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "sequencing reads options")
        group.add_argument(
            "--reads-dir", "-R",
            dest="reads_dir",
            help="data directory where sequencing reads are located",
            default=config.mapping.reads_dir
        )
        group.add_argument(
            "--reads-suffix",
            dest="reads_suffix",
            help="reads files suffix pattern",
            default=config.mapping.reads_suffix)

        return group

    def reads_dir_updates(self, args):
        assert self.subparser is not None

        if args.reads_dir is not None:
            self.config.mapping.reads_dir = args.reads_dir
        if args.reads_suffix is not None:
            self.config.mapping.reads_suffix = args.reads_suffix

    def mapping_dir_options(self, config):
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "mapping files options")
        group.add_argument(
            "--mapping-dir", "-M",
            dest="mapping_dir",
            help="data directory where mapping files are located",
            default=config.mapping.mapping_dir
        )
        group.add_argument(
            "--mapping-suffix",
            dest="mapping_suffix",
            help="mapping files suffix pattern",
            default=config.mapping.reads_suffix)

        return group

    def mapping_dir_updates(self, args):
        assert self.subparser is not None

        if args.mapping_dir is not None:
            self.config.mapping.mapping_dir = args.mapping_dir
        if args.mapping_suffix is not None:
            self.config.mapping.mapping_suffix = args.mapping_suffix

    def mapping_options(self, config):
        self.reads_dir_options(config=config)
        group = self.mapping_dir_options(config=config)
        group.add_argument(
            "--mapping-bowtie-opts",
            dest="mapping_bowtie_opts",
            help="bowtie mapping options",
            default=config.mapping.mapping_bowtie_opts
        )

    def mapping_updates(self, args):
        self.reads_dir_updates(args)
        self.mapping_dir_updates(args)
        if args.mapping_bowtie_opts:
            self.config.mapping.mapping_bowtie_opts = \
                args.mapping_bowtie_opts


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
        self.genome_index_update(args, input_dir=False)

    def run(self, args):
        print("mapping subcommand called with args: {}".format(args))
        self.process_args(args)
        pipeline = MappingPipeline(self.config)
        pipeline.run()
