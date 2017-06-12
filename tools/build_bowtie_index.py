#!/usr/bin/env python
# encoding: utf-8
'''
build_bowtie_index -- builds bowtie index for genome specified

Created on Jun 12, 2017

@author: lubo
'''

import sys
import os

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import common_arguments
import glob

__all__ = []
__version__ = 0.1
__date__ = '2017-06-12'
__updated__ = '2017-06-12'

DEBUG = 1
TESTRUN = 0
PROFILE = 0


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

        # process arguments
        args = parser.parse_args(argv[1:])
        config = common_arguments.process_genome_agrments(args)

        chrom_files = glob.glob("{}/*.fa".format(config.genome.src))
        command = "cp {} {}".format(
            " ".join(chrom_files),
            config.genome.dst
        )
        os.system(command)

        chrom_y_file = "{}/chrY.fa".format(config.genome.dst)
        command = "rm -f {}".format(chrom_y_file)
        os.system(command)

        chrom_files = glob.glob(
            "{}/chr*.fa".format(config.genome.dst)
        )
        command = "cat {} > {}/genome.fa".format(
            " ".join(chrom_files),
            config.genome.dst
        )
        os.system(command)

        os.chdir(config.genome.dst)
        command = "bowtie2-build -f genome.fa genomeindex"
        os.system(command)

        return 0
    except KeyboardInterrupt:
        return 0
    except Exception as e:

        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        sys.stderr.write('\n')
        return 2


if __name__ == "__main__":
    sys.exit(main())
