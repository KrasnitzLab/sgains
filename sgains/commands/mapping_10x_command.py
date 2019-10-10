'''
Created on Aug 2, 2017

@author: lubo
'''
import argparse

from termcolor import colored

from sgains.commands.common import GenomeIndexMixin, OptionsBase, MappingMixin
from sgains.pipelines.mapping_10x_pipeline import Mapping10xPipeline



class Mapping10xMixin(object):
    __slots__ = ()

    def reads_dir_options(self, config):
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "sequencing reads options")
        group.add_argument(
            "--data-10x-dir",
            dest="data_10x_dir",
            help="data directory where 10xGenomics dataset is located",
            default=config.build_data_10x_dir()
        )
        group.add_argument(
            "--data-10x-prefix",
            dest="data_10x_prefix",
            help="10xGenomics dataset common prefix",
            default=config.mapping_10x.data_10x_prefix)

        group.add_argument(
            "--data-10x-summary",
            dest="data_10x_summary",
            help="10xGenomics dataset per cell summary filename",
            default=config.build_data_10x_summary())

        group.add_argument(
            "--data-10x-bam",
            dest="data_10x_bam",
            help="10xGenomics dataset possition sorted BAM file",
            default=config.build_data_10x_bam())

        group.add_argument(
            "--data-10x-bai",
            dest="data_10x_bai",
            help="10xGenomics dataset possition sorted BAM file index",
            default=config.build_data_10x_bai())

        return group

    def reads_dir_updates(self, args):
        assert self.subparser is not None

        if args.data_10x_dir is not None:
            self.config.mapping_10x.data_10x_dir = args.data_10x_dir
        if args.data_10x_prefix is not None:
            self.config.mapping_10x.data_10x_prefix = args.data_10x_prefix

        if args.data_10x_summary is not None:
            self.config.mapping_10x.data_10x_summary = args.data_10x_summary
        if args.data_10x_bam is not None:
            self.config.mapping_10x.data_10x_bam = args.data_10x_bam
        if args.data_10x_bai is not None:
            self.config.mapping_10x.data_10x_bai = args.data_10x_bai


    def mapping_options(self, config):
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "mapping files options")
        group.add_argument(
            "--mapping-10x-dir", "-M",
            dest="mapping_10x_dir",
            help="data directory where mapping files are located",
            default=config.mapping_10x.mapping_10x_dir
        )
        group.add_argument(
            "--mapping-10x-suffix",
            dest="mapping_10x_suffix",
            help="mapping files suffix pattern",
            default=config.mapping_10x.mapping_10x_suffix)

        group.add_argument(
            "--mapping-10x-bowtie-opts",
            dest="mapping_10x_bowtie_opts",
            help="additional bowtie options for 10xGenomics dataset",
            default=config.mapping_10x.mapping_10x_bowtie_opts)

        return group

    def mapping_updates(self, args):
        assert self.subparser is not None

        if args.mapping_10x_dir is not None:
            self.config.mapping_10x.mapping_10x_dir = args.mapping_10x_dir
        if args.mapping_10x_suffix is not None:
            self.config.mapping_10x.mapping_10x_suffix = args.mapping_10x_suffix
        if args.mapping_10x_bowtie_opts is not None:
            self.config.mapping_10x.mapping_10x_bowtie_opts = \
                args.mapping_10x_bowtie_opts

    def mapping_10x_options(self, config):
        self.reads_dir_options(config=config)
        self.mapping_options(config=config)

    def mapping_10x_updates(self, args):
        self.reads_dir_updates(args)
        self.mapping_updates(args)


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
