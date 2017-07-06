#!/usr/bin/env python
# encoding: utf-8

'''
generate_reads -- generates mappable regions for whole genome

Created on Jun 12, 2017

@author: lubo
'''

import sys
import os
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import config
from hg19 import HumanGenome19
import traceback
from common_arguments import Parser


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
        argparser = ArgumentParser(
            description=program_description,
            formatter_class=RawDescriptionHelpFormatter)
        parser = Parser.from_argument_parser(argparser)

        config = parser.parse_arguments(argv[1:])

        hg = None
        if config.genome.version == 'hg19':
            hg = HumanGenome19(config)

        if hg is None:
            raise CLIError("wrong genome version")

        if config.output is None:
            filename = config.bin_boundaries_filename()
            outfile = open(config.abspath(filename), 'w')
        elif config.output == '-':
            outfile = sys.stdout
        else:
            outfile = config.abspath(config.output)
            outfile = open(outfile, "w")

        for rec in hg.generate_reads(config.chroms, config.reads.length):
            outfile.write(">{}\n".format(rec.id))
            outfile.write(str(rec.seq).upper())
            outfile.write("\n")

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
