import os
import pytest
import pickle

from sgains.configuration.parser import Config


def relative_to_this_fixtures_folder(path):
    return os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "fixtures", path)


@pytest.fixture(scope="session")
def fixture_filename(request):
    def builder(relpath):
        return relative_to_this_fixtures_folder(relpath)
    return builder


@pytest.mark.parametrize("filename", [
    "test_sgains_config.yml",
    "test_sgains_config_10x.yml",
    "test_sgains_config_cshl.yml",
])
def test_config_parser_simple(fixture_filename, filename, mocker):
    filename = fixture_filename(filename)
    assert os.path.exists(filename)

    mocker.patch("os.path.exists", return_value=True)

    result = Config.parse(filename)
    assert result is not None

    print(result)
    print(type(result))

    assert result.aligner.aligner_name == "hisat2"

@pytest.mark.parametrize("filename", [
    "test_sgains_config.yml",
    "test_sgains_config_10x.yml",
    "test_sgains_config_cshl.yml",
])
def test_config_parser_pickle(fixture_filename, filename, mocker):
    filename = fixture_filename(filename)
    assert os.path.exists(filename)

    mocker.patch("os.path.exists", return_value=True)

    config = Config.parse(filename)
    assert config is not None

    print(config)

    buffer = pickle.dumps(config)
    assert buffer is not None

    result = pickle.loads(buffer)

    print(result)
    assert result is not None
