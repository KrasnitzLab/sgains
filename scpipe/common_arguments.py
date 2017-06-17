'''
Created on Jun 12, 2017

@author: lubo
'''
import os
from config import Config
from box import Box


def genome_arguments(parser):
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

    return parser


def process_genome_agrments(args):
    config = args.config
    genome = args.genome
    dst = args.dst
    src = args.src

    if config:
        assert os.path.exists(config)
        config = Config.load(config)
    else:
        config = Box(default_box=True)

    if genome:
        config.genome.version = genome
    if dst:
        config.genome.dst = dst
    if src:
        config.genome.src = src

    return config
