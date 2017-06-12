#!/usr/bin/env python
# encoding: utf-8
'''
mask_chr_y -- masks pseudoautosomal region of Y chromsome

Created on Jun 12, 2017

@author: lubo
'''
import sys
import os
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from config import load_config
import config
from box import Box
from hg19 import HumanGenome19


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
        parser.add_argument(
            "-v", "--verbose",
            dest="verbose",
            action="count",
            help="set verbosity level [default: %(default)s]")
        parser.add_argument(
            "-g", "--genome",
            dest="genome",
            help="genome version",
            metavar="path")
        parser.add_argument(
            "-s", "--src",
            dest="src",
            help="source directory",
            metavar="path")
        parser.add_argument(
            "-d", "--dst",
            dest="dst",
            help="destination output directory",
            metavar="path")
        parser.add_argument(
            "-c", "--config",
            dest="config",
            help="configuration file",
            metavar="path")

        # process arguments
        args = parser.parse_args(argv[1:])

        config = args.config
        genome = args.genome
        dst = args.dst
        src = args.src

        if config:
            assert os.path.exists(config)
            config = load_config(config)
        else:
            config = Box(default_box=True)

        if genome:
            config.genome.version = genome
        if dst:
            config.genome.dst = dst
        if src:
            config.genome.src = src

        print(config.to_dict())

        generator = None
        if config.genome.version == 'hg19':
            generator = HumanGenome19(config)

        if generator is None:
            raise CLIError("wrong genome version")

        generator.mask_pseudoautosomal_chrY()

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
