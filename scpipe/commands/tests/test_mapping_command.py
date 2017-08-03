'''
Created on Aug 2, 2017

@author: lubo
'''


def test_bins_long(
        argparser, tests_config, mapping_command):
    mapping_command.add_options()

    argv = [
        "--dry-run", "--force",
        "mapping",
        "--genome-index", "probaindex",
        "--genome-dir", "data/hg19_safe",
        "--bowtie-opts", "-1 -2 -3",
        "--data-dir", "data/test_study/raw",
        "--glob", "*.fastq.gz",
        "--work-dir", "data/proba/bams",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    assert tests_config.force
    assert tests_config.dry_run

    assert tests_config.genome.work_dir == "data/hg19_safe"
    assert tests_config.genome.index == "probaindex"

    assert tests_config.mapping.work_dir == "data/proba/bams"
    assert tests_config.mapping.data_dir == "data/test_study/raw"
    assert tests_config.mapping.data_glob == "*.fastq.gz"
    assert tests_config.mapping.bowtie_opts == "-1 -2 -3"
