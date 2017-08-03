'''
Created on Aug 2, 2017

@author: lubo
'''


def test_bins_long(
        argparser, tests_config, bins_command):
    bins_command.add_options()

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

    assert tests_config.force
    assert tests_config.dry_run

    assert tests_config.mappable_regions.work_dir == "data/R100"
    assert tests_config.mappable_regions.mappable_regions == \
        "mappable_regions.tsv"

    assert tests_config.bins.work_dir == "data/proba"
    assert tests_config.bins.bins_boundaries == "test_bins_boundaries.txt"
    assert tests_config.bins.bins_count == 33
