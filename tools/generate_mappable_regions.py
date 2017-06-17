#!/usr/bin/env python
# encoding: utf-8

'''
generate_mappable_regions -- generates mappable regions for whole genome

Created on Jun 12, 2017

@author: lubo
'''

import sys
import os
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import config
from hg19 import HumanGenome19
import common_arguments
import traceback


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
        parser = ArgumentParser(
            description=program_description,
            formatter_class=RawDescriptionHelpFormatter)
        common_arguments.genome_arguments(parser)

        parser.add_argument(
            "-C", "--chrom",
            dest="chrom",
            help="chromosome for which to generate mappable regions")

        parser.add_argument(
            "-o", "--outfile",
            dest="outfile",
            help="output file to write generated reads. "
            "stdout if not specified",
            metavar="path")

        parser.add_argument(
            "-l", "--length",
            dest="length",
            type=int,
            help="read length to generate")

        # process arguments
        args = parser.parse_args(argv[1:])
        config = common_arguments.process_genome_agrments(args)

        hg = None
        if config.genome.version == 'hg19':
            hg = HumanGenome19(config)

        if hg is None:
            raise CLIError("wrong genome version")

        chrom = args.chrom
        if chrom is not None:
            chroms = [chrom]
        else:
            chroms = hg.CHROMS

        length = args.length
        assert length is not None

        outfile = None
        if args.outfile is None:
            os.makedirs(
                config.abspath(config.bins.cache_dir)
            )
            filename = os.path.join(
                config.bins.cache_dir,
                config.bins.mappable_regions
            )
            outfile = open(config.abspath(filename), 'w')
        elif args.outfile == '-':
            outfile = sys.stdout
        else:
            outfile = os.path.abspath(args.outfile)
            outfile = open(outfile, "w")

        mappings_generator = hg.mappings_generator(chroms, length)
        generator = hg.mappable_regions_generator(mappings_generator)
        for mappable_region in generator:
            outfile.write(str(mappable_region))
            outfile.write('\n')

        outfile.close()
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
