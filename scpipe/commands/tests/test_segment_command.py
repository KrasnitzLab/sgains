'''
Created on Aug 2, 2017

@author: lubo
'''


def test_segment_long(
        argparser, tests_config, segment_command):
    segment_command.add_options(tests_config)

    argv = [
        "--dry-run", "--force",
        "segment",
        "--bins-boundaries", "test_bins_boundaries.txt",
        "--bins-dir", "data/test_study/bins",
        "--data-dir", "data/test_study/results",
        "--glob", "*.varbin.txt",
        "--work-dir", "data/proba",
        "--study-name", "test_study",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    config = segment_command.config
    assert config is not None

    assert config.force
    assert config.dry_run

    assert config.bins.work_dir == "data/test_study/bins"
    assert config.bins.bins_boundaries == "test_bins_boundaries.txt"

    assert config.segment.data_dir == "data/test_study/results"
    assert config.segment.data_glob == "*.varbin.txt"
    assert config.segment.study_name == "test_study"
    assert config.segment.work_dir == "data/proba"
