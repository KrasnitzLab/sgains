'''
Created on Aug 2, 2017

@author: lubo
'''


def test_genomeindex_long(argparser, tests_config, genomeindex_command):
    genomeindex_command.add_options()

    argv = [
        "--dry-run", "--force",
        "genomeindex",
        "--genome-dir", "data/proba",
        "--genome-index", "probaindex",
        "--genome-pristine", "data/hg19_safe/",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    assert tests_config.force
    assert tests_config.dry_run

    assert tests_config.genome.work_dir == "data/proba"
    assert tests_config.genome.index == "probaindex"
    assert tests_config.genome.data_dir == "data/hg19_safe/"


def test_genomeindex_short(argparser, tests_config, genomeindex_command):
    genomeindex_command.add_options()

    argv = [
        "-n", "-F",
        "genomeindex",
        "--genome-dir", "data/proba",
        "-G", "probaindex",
        "--genome-pristine", "data/hg19_safe/",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    assert tests_config.force
    assert tests_config.dry_run

    assert tests_config.genome.work_dir == "data/proba"
    assert tests_config.genome.index == "probaindex"
    assert tests_config.genome.data_dir == "data/hg19_safe/"
