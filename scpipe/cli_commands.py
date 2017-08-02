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
        "--force", "-F",
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
    return group


def _parser_bins_options(parser, defaults_config):
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

    parser_input_data_options(mapping_parser, defaults_config.mapping)
    parser_output_data_options(mapping_parser, defaults_config.mapping)
    group = parser_genome_options(mapping_parser, defaults_config.genome)
    group.add_argument(
        "--bowtie-opts",
        dest="bowtie_opts",
        help="additional bowtie options",
        default=defaults_config.mapping.bowtie_opts)

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

    if args.bowtie_opts:
        defaults_config.mapping.bowtie_ops = args.bowtie_opts

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

    _parser_bins_options(varbin_parser, defaults_config.bins)

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
        help="segments bin counts and prepares the SCGV input data",
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

    _parser_bins_options(segment_parser, defaults_config.bins)

    return segment_parser


def parser_segment_updates(args, defaults_config):
    _args_common_updates(args, defaults_config)
    _args_input_data_updates(args, defaults_config.segment)
    _args_output_data_updates(args, defaults_config.segment)

    if args.study_name is not None:
        defaults_config.segment.study_name = args.study_name

    return defaults_config


def parser_process_options(subparsers, defaults_config):
    process_parser = subparsers.add_parser(
        name="process",
        help="combines mapping, varbin and segment subcommands into "
        "single command",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    output_group = process_parser.add_argument_group("study name")
    output_group.add_argument(
        "--study-name", "-s",
        help="study name",
        dest="study_name",
        default=defaults_config.segment.study_name
    )

    parser_input_data_options(process_parser, defaults_config.mapping)

    parser_output_data_options(process_parser, defaults_config.segment)

    group = parser_genome_options(process_parser, defaults_config.genome)
    group.add_argument(
        "--bowtie-opts",
        dest="bowtie_opts",
        help="additional bowtie options",
        default=defaults_config.mapping.bowtie_opts)

    _parser_bins_options(process_parser, defaults_config.bins)

    return process_parser


def parser_process_updates(args, defaults_config):
    _args_common_updates(args, defaults_config)
    _args_input_data_updates(args, defaults_config.mapping)
    _args_output_data_updates(args, defaults_config.segment)
    _args_genome_updates(args, defaults_config)

    if args.study_name is not None:
        defaults_config.segment.study_name = args.study_name
    if args.bowtie_opts:
        defaults_config.mapping.bowtie_opts = args.bowtie_opts

    return defaults_config


def parser_genomeindex_updates(args, defaults_config):
    _args_common_updates(args, defaults_config)
    if args.data_dir is not None:
        defaults_config.genome.data_dir = args.data_dir
    _args_output_data_updates(args, defaults_config.genome)
    if args.genome_index is not None:
        defaults_config.genome.index = args.genome_index

    return defaults_config


def parser_genomeindex_options(subparsers, defaults_config):
    genomeindex_parser = subparsers.add_parser(
        name="genomeindex",
        help="build appropriate genome index",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    group = genomeindex_parser.add_argument_group(
        "input data options")
    group.add_argument(
        "--data-dir", "-i",
        dest="data_dir",
        help="input data directory where the input data is located",
        default=defaults_config.genome.data_dir
    )

    parser_output_data_options(genomeindex_parser, defaults_config.genome)

    group = genomeindex_parser.add_argument_group(
        "genome index options")
    group.add_argument(
        "--genome-index", "-G",
        dest="genome_index",
        help="genome index name",
        default=defaults_config.genome.index)

    return genomeindex_parser


def parser_mappable_regions_options(subparsers, defaults_config):
    mappable_regions_parser = subparsers.add_parser(
        name="mappable-regions",
        help="finds all mappable regions in specified genome",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    group = parser_genome_options(
        mappable_regions_parser, defaults_config.genome)
    group.add_argument(
        "--bowtie-opts",
        dest="bowtie_opts",
        help="additional bowtie options",
        default=defaults_config.mappable_regions.bowtie_opts)

    group = mappable_regions_parser.add_argument_group(
        "mappable regions options")
    group.add_argument(
        "--work-dir", "-o",
        dest="work_dir",
        help="output directory where results from processing are stored",
        default=defaults_config.mappable_regions.work_dir
    )
    group.add_argument(
        "--read-length", "-l",
        dest="length",
        type=int,
        help="read length to use for generation of mappable regions",
        default=defaults_config.mappable_regions.length
    )

    return mappable_regions_parser


def parser_mappable_regions_updates(args, defaults_config):
    _args_common_updates(args, defaults_config)
    _args_output_data_updates(args, defaults_config.mappable_regions)
    _args_genome_updates(args, defaults_config)

    if args.bowtie_opts:
        defaults_config.mappable_regions.bowtie_opts = args.bowtie_opts

    return defaults_config


def parser_bins_updates(args, defaults_config):
    _args_common_updates(args, defaults_config)
    _args_output_data_updates(args, defaults_config.bins)
    if args.bins_count:
        defaults_config.bins.bins_count = args.bins_count
    if args.bins_boundaries:
        defaults_config.bins.bins_boundaries = args.bins_boundaries
    if args.data_dir:
        defaults_config.mappable_regions.work_dir = args.data_dir
    if args.mappable_regions:
        defaults_config.mappable_regions.mappable_regions = \
            args.mappable_regions

    return defaults_config


def parser_bins_options(subparsers, defaults_config):
    bins_parser = subparsers.add_parser(
        name="bins",
        help="finds all bins boundaries for specified bins count "
        "and read length",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    group = bins_parser.add_argument_group(
        "bins options")
    group.add_argument(
        "--work-dir", "-o",
        dest="work_dir",
        help="output directory where results from processing are stored",
        default=defaults_config.bins.work_dir
    )
    group.add_argument(
        "--bins-count", "-C",
        dest="bins_count",
        type=int,
        help="number of bins",
        default=defaults_config.bins.bins_count
    )
    group.add_argument(
        "--bins-boundaries", "-B",
        dest="bins_boundaries",
        help="filename where bins boundaries are stored",
        default=defaults_config.bins.bins_boundaries)

    group = bins_parser.add_argument_group(
        "mappable regions options")
    group.add_argument(
        "--data-dir", "-i",
        dest="data_dir",
        help="directory where to find mappable regions file",
        default=defaults_config.mappable_regions.work_dir
    )
    group.add_argument(
        "--mappable-regions", "-M",
        dest="mappable_regions",
        help="filename where mappable regions are stored",
        default=defaults_config.mappable_regions.mappable_regions)

    return bins_parser