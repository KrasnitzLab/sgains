'''
Created on Aug 2, 2017

@author: lubo
'''


def test_sc_clust_long(
        argparser, tests_config, sc_clust_command, mocker):

    mocker.patch("os.path.exists", return_value=True)
    mocker.patch("os.listdir", return_value=[])
    mocker.patch("sgains.config.Config.mapping_reads_filenames")

    sc_clust_command.add_options(tests_config)

    argv = [
        "--dry-run", "--force",
        "--config", "sgains/tests/data/scpipe_tests.yml",
        "scclust",
        "--bins-boundaries", "test_bins_boundaries.txt",
        "--bins-dir", "data/test_study/bins",
        "--varbin-dir", "data/test_study/results",
        "--varbin-suffix", ".varbin.txt",
        "--scclust-dir", "data/proba",
        "--case-name", "test_study",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    config = sc_clust_command.config
    assert config is not None

    assert config.force
    assert config.dry_run

    assert config.bins.bins_dir == "data/test_study/bins"
    assert config.bins.bins_boundaries == "test_bins_boundaries.txt"

    assert config.varbin.varbin_dir == "data/test_study/results"
    assert config.varbin.varbin_suffix == ".varbin.txt"
    assert config.scclust.case_name == "test_study"
    assert config.scclust.scclust_dir == "data/proba"
