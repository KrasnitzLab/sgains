#!/usr/bin/env python
# encoding: utf-8
'''
mappable_regions -- masks pseudoautosomal region of Y chromsome

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
            "-m", "--mapped",
            dest="mapped",
            help="file with mapped reads. stdin if not specified",
            metavar="path")

        parser.add_argument(
            "-o", "--outfile",
            dest="outfile",
            help="output file to write mappable regions. "
            "stdout if not specified",
            metavar="path")

        # process arguments
        args = parser.parse_args(argv[1:])
        config = common_arguments.process_genome_agrments(args)
        mapped = args.mapped
        outfile = args.outfile

        if mapped is None:
            mapped = sys.stdin
        if outfile is None:
            outfile = sys.stdout
        else:
            outfile = open(outfile, "w")

        generator = None
        if config.genome.version == 'hg19':
            generator = HumanGenome19(config)

        if generator is None:
            raise CLIError("wrong genome version")

        for mapping in generator.mappable_generator(mapped):
            outfile.write(str(mapping))
            outfile.write('\n')

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