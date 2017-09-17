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
        "--mapping-dir", "data/test_study/bam",
        "--mapping-suffix", ".rmdup.bam",
        "--varbin-dir", "data/proba/varbin",
        "--varbin-suffix", "varbin.txt",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    config = varbin_command.config
    assert config is not None

    assert config.force
    assert config.dry_run

    assert config.bins.bins_dir == "data/test_study/bins"
    assert config.bins.bins_boundaries == "test_bins_boundaries.txt"
    assert config.mapping.mapping_dir == "data/test_study/bam"
    assert config.mapping.mapping_suffix == ".rmdup.bam"
