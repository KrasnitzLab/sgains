#!/usr/bin/env python
# encoding: utf-8

'''
sgains - sparse genomic analysis of individual nuclei by sequencing pipeline
Created on Jul 20, 2017

@author: lubo
'''
from argparse import ArgumentParser,\
    RawDescriptionHelpFormatter
import os
import sys
import traceback

from commands.bins_command import BinsCommand
from commands.genomeindex_command import GenomeIndexCommand
from commands.mappable_regions_command import MappableRegionsCommand
from commands.mapping_command import MappingCommand
from commands.prepare_command import PrepareCommand
from commands.process_command import ProcessCommand
from commands.segment_command import SegmentCommand
from commands.varbin_command import VarbinCommand
from config import Config
from commands.common import OptionsBase


class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''

    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


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

        OptionsBase.common_options(argparser)

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

        process_command = ProcessCommand(
            defaults_config, argparser, subparsers)
        process_command.add_options()

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
