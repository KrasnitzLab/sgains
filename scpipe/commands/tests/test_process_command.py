'''
Created on Aug 2, 2017

@author: lubo
'''


def test_process_long(
        argparser, tests_config, process_command):
    process_command.add_options(tests_config)

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

    config = process_command.config
    assert config is not None

    assert config.force
    assert config.dry_run

    assert config.bins.work_dir == "data/test_study/bins"
    assert config.bins.bins_boundaries == "test_bins_boundaries.txt"

    assert config.genome.work_dir == "data/hg19_safe"
    assert config.genome.index == "probaindex"

    assert config.mapping.data_dir == "data/test_study/raw"
    assert config.mapping.data_glob == "*.fastq.gz"
    assert config.mapping.bowtie_opts == "-1 -2 -3"

    assert config.segment.study_name == "test_study"
    assert config.segment.work_dir == "data/proba"
