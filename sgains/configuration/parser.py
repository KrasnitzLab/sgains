import os
from collections import namedtuple

import yaml

from cerberus import Validator

from sgains.configuration.schema import sgains_schema


def _dict_to_namedtuple(input_dict, dict_name="root"):
    CONFIG_TUPLE = namedtuple(dict_name, input_dict.keys())

    for key, value in input_dict.items():
        if isinstance(value, dict):
            input_dict[key] = _dict_to_namedtuple(value, key)
        elif isinstance(value, list):
            input_dict[key] = [
                _dict_to_namedtuple(item)
                if isinstance(item, dict)
                else item
                for item in value
            ]

    return CONFIG_TUPLE(*input_dict.values())  # type: ignore


class SgainsValidator(Validator):
    def _normalize_coerce_abspath(self, value: str) -> str:
        directory = self._config["work_dirname"]
        if not os.path.isabs(value):
            value = os.path.join(directory, value)
        return os.path.normpath(value)


class Config:

    def __init__(self, config_tuple):
        self.config = config_tuple

    @staticmethod
    def parse_argv(argv):
        if '-c' in argv:
            index = argv.index('-c')
        elif '--config' in argv:
            index = argv.index('--config')
        else:
            return None

        index += 1
        if index < 0 or index >= len(argv):
            raise ValueError('config filename not found')

        filename = argv[index]
        config = Config.parse(filename)
        return config

    @staticmethod
    def parse(filename):
        assert os.path.exists(filename)

        with open(filename, "r") as infile:
            result = yaml.safe_load(infile)
        
        conf_dirname = os.path.dirname(filename)

        validator = SgainsValidator(
            sgains_schema, work_dirname=conf_dirname)
        assert validator.validate(result), validator.errors

        return Config(_dict_to_namedtuple(validator.document, "sgains"))
