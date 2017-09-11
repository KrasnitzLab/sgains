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

    assert tests_config.force
    assert tests_config.dry_run

    assert tests_config.bins.work_dir == "data/test_study/bins"
    assert tests_config.bins.bins_boundaries == "test_bins_boundaries.txt"
    assert tests_config.varbin.data_dir == "data/test_study/bam"
    assert tests_config.varbin.data_glob == "*.rmdup.bam"
