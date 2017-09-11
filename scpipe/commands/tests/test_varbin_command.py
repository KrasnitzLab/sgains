'''
Created on Aug 2, 2017

@author: lubo
'''


def test_varbin_long(
        argparser, tests_config, varbin_command):
    varbin_command.add_options(tests_config)

    argv = [
        "--dry-run", "--force",
        "varbin",
        "--bins-boundaries", "test_bins_boundaries.txt",
        "--bins-dir", "data/test_study/bins",
        "--data-dir", "data/test_study/bam",
        "--glob", "*.rmdup.bam",
        "--work-dir", "data/proba/varbin",
        "--suffix", "varbin.txt",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    config = varbin_command.config
    assert config is not None

    assert config.force
    assert config.dry_run

    assert config.bins.work_dir == "data/test_study/bins"
    assert config.bins.bins_boundaries == "test_bins_boundaries.txt"
    assert config.varbin.data_dir == "data/test_study/bam"
    assert config.varbin.data_glob == "*.rmdup.bam"
