'''
Created on Aug 2, 2017

@author: lubo
'''


def test_mapping_long(
        argparser, tests_config, mapping_command):
    mapping_command.add_options(tests_config)

    argv = [
        "--dry-run", "--force",
        "--parallel", "10",
        "mapping",
        "--genome-index", "probaindex",
        "--genome-dir", "data/hg19_safe",
        "--mapping-bowtie-opts", "-1 -2 -3",
        "--reads-dir", "data/test_study/raw",
        "--reads-suffix", ".fastq.gz",
        "--mapping-dir", "data/proba/bams",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    config = mapping_command.config
    assert config is not None

    assert config.force
    assert config.dry_run
    assert config.parallel == 10

    assert config.genome.work_dir == "data/hg19_safe"
    assert config.genome.index == "probaindex"

    assert config.mapping.mapping_dir == "data/proba/bams"
    assert config.mapping.reads_dir == "data/test_study/raw"
    assert config.mapping.reads_suffix == ".fastq.gz"
    assert config.mapping.mapping_bowtie_opts == "-1 -2 -3"
