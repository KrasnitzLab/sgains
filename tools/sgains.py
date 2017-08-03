#!/usr/bin/env python
# encoding: utf-8

'''
sgains - sparse genomic analysis of individual nuclei by sequencing pipeline
Created on Jul 20, 2017

@author: lubo
'''
from argparse import ArgumentParser,\
    RawDescriptionHelpFormatter
import functools
import os
import sys
import traceback

from cli_commands import parser_varbin_options, \
    parser_varbin_updates, parser_common_options,\
    parser_segment_options, parser_segment_updates, parser_process_options,\
    parser_process_updates
from commands.bins_command import BinsCommand
from commands.genomeindex_command import GenomeIndexCommand
from commands.mappable_regions_command import MappableRegionsCommand
from config import Config
import config
from mapping_pipeline import MappingPipeline
from r_pipeline import Rpipeline
from varbin_pipeline import VarbinPipeline
from commands.prepare_command import PrepareCommand
from commands.mapping_command import MappingCommand
from commands.varbin_command import VarbinCommand
from commands.segment_command import SegmentCommand


class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''

    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


def do_process(defaults_config, args):
    if args.config is not None:
        config = Config.load(args.config)
        defaults_config.update(config)
    defaults_config = parser_process_updates(args, defaults_config)

    mapping_workdir = os.path.join(
        defaults_config.segment_work_dirname(),
        'mapping')
    varbin_workdir = os.path.join(
        defaults_config.segment_work_dirname(),
        'varbin')
    segment_workdir = os.path.join(
        defaults_config.segment_work_dirname(),
        'segment')

    mapping_config = Config.copy(defaults_config)
    mapping_config.mapping.work_dir = mapping_workdir

    varbin_config = Config.copy(defaults_config)
    varbin_config.varbin.data_dir = mapping_workdir
    varbin_config.varbin.work_dir = varbin_workdir

    segment_config = Config.copy(defaults_config)
    segment_config.segment.data_dir = varbin_workdir
    segment_config.segment.work_dir = segment_workdir

    pipeline = MappingPipeline(mapping_config)
    pipeline.run()

    pipeline = VarbinPipeline(varbin_config)
    pipeline.run()

    pipeline = Rpipeline(segment_config)
    pipeline.run()


def main(argv=None):
    if argv is None:
        argv = sys.argv
    else:
        argv.extend(sys.argv)

    # Setup argument parser
    program_name = os.path.basename(sys.argv[0])
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_description = '''%s

USAGE
''' % (program_shortdesc, )

    try:
        defaults_config = Config.load("sgains.yml")

        argparser = ArgumentParser(
            description=program_description,
            formatter_class=RawDescriptionHelpFormatter)

        parser_common_options(argparser)

        # parser = Parser.from_argument_parser(argparser)
        subparsers = argparser.add_subparsers(
            title="subcommands"
        )

        genomeindex_command = GenomeIndexCommand(
            defaults_config, argparser, subparsers)
        genomeindex_command.add_options()

        mappable_regions_command = MappableRegionsCommand(
            defaults_config, argparser, subparsers)
        mappable_regions_command.add_options()

        bins_command = BinsCommand(
            defaults_config, argparser, subparsers)
        bins_command.add_options()

        prepare_command = PrepareCommand(
            defaults_config, argparser, subparsers)
        prepare_command.add_options()

        mapping_command = MappingCommand(
            defaults_config, argparser, subparsers)
        mapping_command.add_options()

        varbin_command = VarbinCommand(
            defaults_config, argparser, subparsers)
        varbin_command.add_options()

        segment_command = SegmentCommand(
            defaults_config, argparser, subparsers)
        segment_command.add_options()

        process_parser = parser_process_options(subparsers, defaults_config)
        process_parser.set_defaults(
            func=functools.partial(do_process, defaults_config))

        args = argparser.parse_args(argv[1:])
        args.func(args)

        # bin_boundaries(hg, config, outfile)

        return 0
    except KeyboardInterrupt:
        return 0
    except Exception as e:
        traceback.print_exc()

        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        sys.stderr.write('\n')
        return 2


if __name__ == "__main__":
    sys.exit(main())
