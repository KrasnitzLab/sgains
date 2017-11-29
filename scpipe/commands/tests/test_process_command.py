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

        "--reads-dir",
        "data/Navin2011/T10_small/reads/",
        "--reads-suffix", ".fastq.gz",
        "--mapping-bowtie-opts", "-1 -2 -3",

        "--genome-index", "genomeindex",
        "--genome-dir", "data/hg19",

        "--study-name", "test_study",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    config = process_command.config
    assert config is not None

    assert config.force
    assert config.dry_run

    assert config.bins.bins_dir == "data/test_study/bins"
    assert config.bins.bins_boundaries == "test_bins_boundaries.txt"

    assert config.genome.work_dir == "data/hg19"
    assert config.genome.index == "genomeindex"

    assert config.mapping.reads_dir == "data/Navin2011/T10_small/reads/"
    assert config.mapping.reads_suffix == ".fastq.gz"
    assert config.mapping.mapping_bowtie_opts == "-1 -2 -3"

    assert config.segment.study_name == "test_study"
