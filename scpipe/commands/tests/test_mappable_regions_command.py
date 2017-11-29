'''
Created on Aug 2, 2017

@author: lubo
'''


def test_mappable_regions_long(
        argparser, tests_config, mappable_regions_command):
    mappable_regions_command.add_options(tests_config)

    argv = [
        "--dry-run", "--force",
        "mappable-regions",
        "--mappable-dir", "data/proba",
        "--genome-index", "genomeindex",
        "--genome-dir", "data/hg19/",
        "--read-length", "200",
        "--bowtie-opts", "-1 -2 -3",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    config = mappable_regions_command.config
    assert config is not None

    assert config.force
    assert config.dry_run

    assert config.genome.work_dir == "data/hg19/"
    assert config.genome.index == "genomeindex"

    assert config.mappable_regions.length == 200
    assert config.mappable_regions.work_dir == "data/proba"
    assert config.mappable_regions.bowtie_opts == "-1 -2 -3"


def test_mappable_regions_short(
        argparser, tests_config, mappable_regions_command):
    mappable_regions_command.add_options(tests_config)

    argv = [
        "-n", "-F",
        "mappable-regions",
        "-m", "data/proba",
        "-G", "genomeindex",
        "--genome-dir", "data/hg19/",
        "-l", "200",
        "--bowtie-opts", "-1 -2 -3",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    config = mappable_regions_command.config
    assert config is not None

    assert config.force
    assert config.dry_run

    assert config.genome.work_dir == "data/hg19/"
    assert config.genome.index == "genomeindex"

    assert config.mappable_regions.length == 200
    assert config.mappable_regions.work_dir == "data/proba"
    assert config.mappable_regions.bowtie_opts == "-1 -2 -3"
