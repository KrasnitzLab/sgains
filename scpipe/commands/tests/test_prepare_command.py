'''
Created on Aug 2, 2017

@author: lubo
'''


def test_prepare_long(
        argparser, tests_config, prepare_command):
    prepare_command.add_options()

    argv = [
        "--dry-run", "--force",
        "prepare",
        "--mappable-dir", "data/proba",
        "--genome-index", "probaindex",
        "--genome-dir", "data/hg19_safe/",
        "--read-length", "200",
        "--bowtie-opts", "-1 -2 -3",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    assert tests_config.force
    assert tests_config.dry_run

    assert tests_config.genome.work_dir == "data/hg19_safe/"
    assert tests_config.genome.index == "probaindex"

    assert tests_config.mappable_regions.length == 200
    assert tests_config.mappable_regions.work_dir == "data/proba"
    assert tests_config.mappable_regions.bowtie_opts == "-1 -2 -3"
