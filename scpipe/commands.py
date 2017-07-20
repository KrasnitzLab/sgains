'''
Created on Jul 20, 2017

@author: lubo
'''
import argparse


def parser_common_options(parser):
    parser.add_argument(
        "-v", "--verbose",
        dest="verbose",
        action="count",
        help="set verbosity level [default: %(default)s]",
        default=0
    )
    parser.add_argument(
        "-c", "--config",
        dest="config",
        help="configuration file",
        metavar="path"
    )

    parser.add_argument(
        "-n", "--dry-run",
        dest="dry_run",
        action="store_true",
        help="perform a trial run with no changes made",
        default=False
    )

    parser.add_argument(
        "--force",
        dest="force",
        action="store_true",
        help="allows overwriting nonempty results directory",
        default=False
    )


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


def parser_bins_options(parser, defaults_config):
    group = parser.add_argument_group(
        "bins boundaries")
    group.add_argument(
        "--bins-boundaries", "-B",
        dest="bins_boundaries",
        help="bins boundaries filename",
        default=defaults_config.bins_boundaries)

    group.add_argument(
        "--bins-dir",
        dest="bins_dir",
        help="bins working directory",
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


def _args_common_updates(args, defaults_config):
    if args.dry_run is not None:
        defaults_config.dry_run = args.dry_run
    if args.force is not None:
        defaults_config.force = args.force


def _args_input_data_updates(args, defaults_config):
    if args.data_dir is not None:
        defaults_config.data_dir = args.data_dir
    if args.data_glob is not None:
        defaults_config.data_glob = args.data_glob


def _args_output_data_updates(args, defaults_config):
    if args.work_dir is not None:
        defaults_config.work_dir = args.work_dir


def _args_genome_updates(args, defaults_config):
    if args.genome_index is not None:
        defaults_config.genome.index = args.genome_index
    if args.genome_dir is not None:
        defaults_config.genome.work_dir = args.genome_dir


def parser_mapping_updates(args, defaults_config):
    _args_common_updates(args, defaults_config)
    _args_input_data_updates(args, defaults_config.mapping)
    _args_output_data_updates(args, defaults_config.mapping)
    _args_genome_updates(args, defaults_config)

    return defaults_config


def parser_varbin_options(subparsers, defaults_config):
    varbin_parser = subparsers.add_parser(
        name="varbin",
        help="applies varbin algorithm to count read mappings in each bin",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser_input_data_options(varbin_parser, defaults_config.varbin)
    parser_output_data_options(varbin_parser, defaults_config.varbin)
    output_group = varbin_parser.add_argument_group("output suffix")
    output_group.add_argument(
        "--suffix", "-s",
        help="suffix for output files",
        dest="suffix",
        default=defaults_config.varbin.suffix
    )

    parser_bins_options(varbin_parser, defaults_config.bins)

    return varbin_parser


def parser_varbin_updates(args, defaults_config):
    _args_common_updates(args, defaults_config)
    _args_input_data_updates(args, defaults_config.varbin)
    _args_output_data_updates(args, defaults_config.varbin)

    if args.suffix is not None:
        defaults_config.varbin.suffix = args.suffix

    return defaults_config


def parser_segment_options(subparsers, defaults_config):
    segment_parser = subparsers.add_parser(
        name="segment",
        help="prepares the SCGV input data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser_input_data_options(segment_parser, defaults_config.segment)
    parser_output_data_options(segment_parser, defaults_config.segment)
    output_group = segment_parser.add_argument_group("study name")
    output_group.add_argument(
        "--study-name", "-s",
        help="study name",
        dest="study_name",
        default=defaults_config.segment.study_name
    )

    parser_bins_options(segment_parser, defaults_config.bins)

    return segment_parser


def parser_segment_updates(args, defaults_config):
    _args_common_updates(args, defaults_config)
    _args_input_data_updates(args, defaults_config.segment)
    _args_output_data_updates(args, defaults_config.segment)

    if args.study_name is not None:
        defaults_config.segment.study_name = args.study_name

    return defaults_config
