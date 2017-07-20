'''
Created on Jul 20, 2017

@author: lubo
'''
import argparse


def parser_input_data_options(parser, config):
    group = parser.add_argument_group(
        "input data options")
    group.add_argument(
        "--data-dir", "-i",
        dest="data_dir",
        help="input data directory where the input data is located",
        default=config.data_dir
    )
    group.add_argument(
        "--glob", "-g",
        dest="data_glob",
        help="glob pattern for finding input data",
        default=config.data_glob)


def parser_output_data_options(parser, config):
    group = parser.add_argument_group(
        "output data options")
    group.add_argument(
        "--work-dir", "-o",
        dest="work_dir",
        help="output directory where results from processing are stored",
        default=config.work_dir
    )


def parser_genome_options(parser, defaults_config):
    group = parser.add_argument_group(
        "genome index options")
    group.add_argument(
        "--genome-index", "-G",
        dest="genome_index",
        help="genome index name",
        default=defaults_config.index)

    group.add_argument(
        "--genome-dir",
        dest="genome_dir",
        help="genome index directory",
        default=defaults_config.work_dir)


def parser_mapping_options(subparsers, defaults_config):
    mapping_parser = subparsers.add_parser(
        name="mapping",
        help="performs actual mapping of cell reads",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    mapping_parser.add_argument(
        "--bowtie-opts",
        dest="bowtie_opts",
        help="additional bowtie options",
        default=defaults_config.mapping.bowtie_opts)

    parser_input_data_options(mapping_parser, defaults_config.mapping)
    parser_output_data_options(mapping_parser, defaults_config.mapping)
    parser_genome_options(mapping_parser, defaults_config.genome)
    return mapping_parser


def parser_mapping_updates(args, defaults_config):
    if args.dry_run is not None:
        defaults_config.dry_run = args.dry_run
    if args.force is not None:
        defaults_config.force = args.force

    if args.data_dir is not None:
        defaults_config.mapping.data_dir = args.data_dir
    if args.work_dir is not None:
        defaults_config.mapping.work_dir = args.work_dir
    if args.data_glob is not None:
        defaults_config.mapping.data_glob = args.data_glob
    if args.genome_index is not None:
        defaults_config.genome.index = args.genome_index
    if args.genome_dir is not None:
        defaults_config.genome.work_dir = args.genome_dir

    return defaults_config
