'''
Created on Aug 2, 2017

@author: lubo
'''


def test_genomeindex_long(
        argparser, tests_config, genomeindex_command, mocker):

    mocker.patch("os.path.exists", return_value=True)
    mocker.patch("os.listdir", return_value=[])

    genomeindex_command.add_options(tests_config)
    argv = [
        "--dry-run", "--force",
        "genomeindex",
        "--genome-dir", "data/proba",
        "--genome-index", "probaindex",
        "--genome-pristine", "data/hg19_safe/",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    config = genomeindex_command.config
    assert config is not None

    assert config.force
    assert config.dry_run

    assert config.genome.work_dir == "data/proba"
    assert config.genome.index == "probaindex"
    assert config.genome.data_dir == "data/hg19_safe/"


def test_genomeindex_short(
        argparser, tests_config, genomeindex_command, mocker):

    mocker.patch("os.path.exists", return_value=True)
    mocker.patch("os.listdir", return_value=[])

    genomeindex_command.add_options(tests_config)

    argv = [
        "-n", "-F",
        "genomeindex",
        "--genome-dir", "data/proba",
        "-G", "probaindex",
        "--genome-pristine", "data/hg19_safe/",
    ]

    args = argparser.parse_args(argv)
    args.func(args)

    config = genomeindex_command.config
    assert config is not None

    assert config.force
    assert config.dry_run

    assert config.genome.work_dir == "data/proba"
    assert config.genome.index == "probaindex"
    assert config.genome.data_dir == "data/hg19_safe/"
