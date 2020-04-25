import os
import sys
from copy import deepcopy

import traceback
import functools
from collections import defaultdict

import yaml

from argparse import ArgumentParser,\
    RawDescriptionHelpFormatter, ArgumentDefaultsHelpFormatter

from sgains.configuration.parser import SgainsValidator, Config
from sgains.configuration.schema import sgains_schema

from sgains.commands.common import Command

from sgains.pipelines.mappableregions_pipeline import MappableRegionsPipeline
from sgains.pipelines.genomeindex_pipeline import GenomeIndexPipeline
from sgains.pipelines.bins_pipeline import BinsPipeline
from sgains.pipelines.mapping_pipeline import MappingPipeline
from sgains.pipelines.extract_10x_pipeline import Extract10xPipeline

from sgains.pipelines.varbin_pipeline import VarbinPipeline
from sgains.pipelines.r_pipeline import Rpipeline


SGAINS_COMMANDS = {
    "genome": {
        "config_groups": ["aligner", "genome"],
        "help": "builds appropriate bowtie index for the reference genome",
    },
    "mappable_regions": {
        "config_groups": ["aligner", "genome", "mappable_regions", "sge"],
        "help": "finds all mappable regions in specified genome",
    },
    "bins": {
        "config_groups": ["genome", "mappable_regions", "bins", "sge"],
        "help": "calculates all bins boundaries for specified bins count "
        "and read length",
    },
    "mapping": {
        "config_groups": ["aligner", "genome", "reads", "mapping", "sge"],
        "help": "performs mapping of cells reads to the reference genome",
    },
    "extract_10x": {
        "config_groups": [
            "data_10x", "reads", "sge"],
        "help": "extracts cells reads from 10x Genomics datasets",
    },
    "varbin": {
        "config_groups": ["bins", "mapping", "varbin", "sge"],
        "help": "applies varbin algorithm to count read mappings in each bin",
    },
    "scclust": {
        "config_groups": ["bins", "varbin", "scclust"],
        "help": "segmentation and clustering based bin counts and "
        "preparation of the SCGV input data"
    },
}


def build_common_options(parser):
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

    parser.add_argument(
        "--parallel", "-p",
        dest="parallel",
        help="number of task to run in parallel",
        type=int,
        default=1
    )

    parser.add_argument(
        "--sge",
        dest="sge",
        action="store_true",
        help="parallelilizes commands using SGE cluster manager",
        default=False
    )


def _get_config_value(config, group_name, name):

    if config is None:
        return None
    group = config.config.get(group_name)
    if group is None:
        return None
    result = getattr(group, name)
    print("getting config for:", group_name, name, "->", result)
    return result


def build_cli_options(argparser, command=None, config=None):
    work_dirname = os.getcwd()
    if config is not None:
        work_dirname = config.work_dir

    validator = SgainsValidator(
        deepcopy(sgains_schema), work_dirname=work_dirname)

    if command is None:
        config_groups = list(validator.schema.keys())
    else:
        assert command in SGAINS_COMMANDS

        command = SGAINS_COMMANDS[command]
        config_groups = command["config_groups"]

    for group_name in config_groups:

        group = validator.schema.get(group_name)
        group_parser = argparser.add_argument_group(f"{group_name} group:")
        assert group["type"] == "dict", (group_name, group)

        group_schema = group["schema"]

        for arg_name, arg_spec in group_schema.items():
            name = f"--{arg_name.replace('_', '-')}"
            arg_type = str
            arg_type = arg_spec.get("type", "string")
            if arg_type == "string":
                arg_type = str
            elif arg_type == "integer":
                arg_type = int
            elif arg_type == "float":
                arg_type = float
            elif arg_type == "list":
                arg_type = list
            else:
                raise ValueError(f"wrong argument type {arg_type}")

            help_data = None
            meta_data = arg_spec.get("meta")
            if meta_data is not None:
                help_data = meta_data.get("help")

            arg_default = _get_config_value(config, group_name, arg_name)
            if arg_default is None:
                arg_default = arg_spec.get("default")

            group_parser.add_argument(
                name,
                help=help_data,
                dest=arg_name,
                type=arg_type,
                default=arg_default)

    return argparser


def parse_cli_options(args):
    config_dict = defaultdict(dict)
    work_dirname = os.getcwd()

    if args.config is not None:
        assert os.path.exists(args.config), args.config

        with open(args.config, "r") as infile:
            config_dict = yaml.safe_load(infile)
        work_dirname = os.path.dirname(args.config)

    validator = SgainsValidator(
        deepcopy(sgains_schema), work_dirname=work_dirname)

    result = defaultdict(dict)

    config_groups = list(validator.schema.keys())

    for group_name in config_groups:

        group = validator.schema.get(group_name)
        group_schema = group.get("schema")
        if group_schema is None:
            continue

        group_result = {}
        for arg_name in group_schema.keys():
            arg_value = getattr(args, arg_name, None)
            if arg_value is not None:
                group_result[arg_name] = arg_value
            else:
                config_value = config_dict.get(group_name, None)
                if config_value is not None:
                    config_value = config_value.get(arg_name, None)
                if config_value is not None:
                    group_result[arg_name] = config_value
        if group_result:
            result[group_name] = group_result

    config = Config.from_dict(result, work_dirname)
    config.verbose = args.verbose
    config.config_file = args.config
    config.dry_run = args.dry_run
    config.force = args.force
    config.parallel = args.parallel
    config.sge = args.sge

    return config


def main(argv=sys.argv[1:]):
    program_name = os.path.basename(sys.argv[0])
    program_shortdesc = \
        'sgains - sparse genomic analysis of individual nuclei by ' \
        'sequencing pipeline'
    program_description = '''%s

USAGE
''' % (program_shortdesc, )

    try:

        config = Config.parse_argv(argv)

        argparser = ArgumentParser(
            description=program_description,
            formatter_class=ArgumentDefaultsHelpFormatter)
        
        build_common_options(argparser)
        subparsers = argparser.add_subparsers(
            title="sGAINS subcommands"
        )

        for command in SGAINS_COMMANDS:
            command_name = command.replace("_", "-")
            command_help = SGAINS_COMMANDS[command].get("help", "")
            subparser = subparsers.add_parser(
                name=command_name,
                help=command_help,
                formatter_class=ArgumentDefaultsHelpFormatter
            )

            build_cli_options(subparser, command, config)
            subparser.set_defaults(func=functools.partial(execute, command))

        args = argparser.parse_args(argv)
        args.func(args)

    except KeyboardInterrupt:
        traceback.print_exc()
        return 0
    except Exception as e:
        traceback.print_exc()

        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        sys.stderr.write('\n')
        return 2

def create_pipeline(command, config):
    if command == "genome":
        return GenomeIndexPipeline(config)
    elif command == "mappable_regions":
        return MappableRegionsPipeline(config)
    elif command == "bins":
        return BinsPipeline(config)
    elif command == "mapping":
        return MappingPipeline(config)
    elif command == "varbin":
        return VarbinPipeline(config)
    elif command == "scclust":
        return Rpipeline(config)
    elif command == "extract_10x":
        return Extract10xPipeline(config)

    raise ValueError(f"Unexpected command: {command}")

def execute(command, args):
    print("EXECUTE!!!!")
    print(command, args)
    config = parse_cli_options(args)
    print(config.config)

    pipeline = create_pipeline(command, config)
    assert pipeline is not None, command

    executor = Command(config)
    executor.run_pipeline(pipeline)


if __name__ == "__main__":
    sys.exit(main())
