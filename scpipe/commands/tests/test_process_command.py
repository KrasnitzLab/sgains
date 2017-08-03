'''
Created on Aug 2, 2017

@author: lubo
'''


def test_process_long(
        argparser, tests_config, process_command):
    process_command.add_options()

    argv = [
        "--dry-run", "--force",
        "process",
        "--bins-boundaries", "test_bins_boundaries.txt",
        "--bins-dir", "data/test_study/bins",

        "--data-dir", "data/test_study/raw",
        "--glob", "*.fastq.gz",
        "--bowtie-opts", "-1 -2 -3",

        "--genome-index", "probaindex",
        "--genome-dir", "data/hg19_safe",

        "--work-dir", "data/proba",
        "--study-name", "test_study",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    assert tests_config.force
    assert tests_config.dry_run

    assert tests_config.bins.work_dir == "data/test_study/bins"
    assert tests_config.bins.bins_boundaries == "test_bins_boundaries.txt"

    assert tests_config.genome.work_dir == "data/hg19_safe"
    assert tests_config.genome.index == "probaindex"

    assert tests_config.mapping.data_dir == "data/test_study/raw"
    assert tests_config.mapping.data_glob == "*.fastq.gz"
    assert tests_config.mapping.bowtie_opts == "-1 -2 -3"

    assert tests_config.segment.study_name == "test_study"
    assert tests_config.segment.work_dir == "data/proba"
