'''
Created on Aug 2, 2017

@author: lubo
'''


def test_bins_long(
        argparser, tests_config, bins_command):
    bins_command.add_options(tests_config)

    argv = [
        "--dry-run", "--force",
        "bins",
        "--bins-dir", "data/proba",
        "--bins-boundaries", "test_bins_boundaries.txt",
        "--bins-count", "33",
        "--mappable-dir", "data/R100",
        "--mappable-regions", "mappable_regions.tsv",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    config = bins_command.config
    assert config is not None

    assert config.force
    assert config.dry_run

    assert config.mappable_regions.work_dir == "data/R100"
    assert config.mappable_regions.mappable_regions == \
        "mappable_regions.tsv"

    assert config.bins.bins_dir == "data/proba"
    assert config.bins.bins_boundaries == "test_bins_boundaries.txt"
    assert config.bins.bins_count == 33
