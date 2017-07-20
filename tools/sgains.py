#!/usr/bin/env python
# encoding: utf-8

'''
sgains - sparse genomic analysis of individual nuclei by sequencing pipeline
Created on Jul 20, 2017

@author: lubo
'''
import sys
import os
from argparse import ArgumentParser,\
    RawDescriptionHelpFormatter
import config
import traceback
from config import Config
from commands import parser_mapping_options, parser_mapping_updates
import functools
from mapping_pipeline import MappingPipeline


class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''

    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


def bin_boundaries(hg, config, outfile):
    regions_df = hg.load_mappable_regions()
    bins_df = hg.calc_bin_boundaries(
        config.chroms,
        regions_df
    )
    df = hg.calc_bins_gc_content(config.chroms, bins_df)

    df.to_csv(outfile, sep='\t', index=False)


def do_mapping(defaults_config, args):
    if args.config is not None:
        config = Config.load(args.config)
        defaults_config.update(config)

    defaults_config = parser_mapping_updates(args, defaults_config)
    pipeline = MappingPipeline(defaults_config)
    pipeline.run()


def main(argv=None):
    if argv is None:
        argv = sys.argv
    else:
        argv.extend(sys.argv)

    print(argv)

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

        argparser.add_argument(
            "-v", "--verbose",
            dest="verbose",
            action="count",
            help="set verbosity level [default: %(default)s]",
            default=0
        )
        argparser.add_argument(
            "-c", "--config",
            dest="config",
            help="configuration file",
            metavar="path"
        )

        argparser.add_argument(
            "-n", "--dry-run",
            dest="dry_run",
            action="store_true",
            help="perform a trial run with no changes made",
            default=False
        )

        argparser.add_argument(
            "--force",
            dest="force",
            action="store_true",
            help="allows overwriting nonempty results directory",
            default=False
        )

        # parser = Parser.from_argument_parser(argparser)
        subparsers = argparser.add_subparsers(
            title="subcommands"
        )

        mapping_parser = parser_mapping_options(subparsers, defaults_config)
        mapping_parser.set_defaults(
            func=functools.partial(do_mapping, defaults_config))

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
