import os
from collections import namedtuple

from box import Box
import yaml
from termcolor import colored

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


def _dict_to_box(input_dict):
    return Box(input_dict, frozen_box=True)


class SgainsValidator(Validator):
    def _normalize_coerce_abspath(self, value: str) -> str:
        directory = self._config["work_dirname"]
        if not os.path.isabs(value):
            value = os.path.join(directory, value)
        return os.path.normpath(value)


class Config:

    def __init__(self, config):
        self.config = config
        self.verbose = 0
        self.config_file = None
        self.dry_run = False
        self.force = False
        self.parallel = 1

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
            config_dict = yaml.safe_load(infile)
        
        conf_dirname = os.path.dirname(filename)
        return Config.from_dict(config_dict, conf_dirname)

    @staticmethod
    def from_dict(config_dict, work_dirname):
        config_dict["work_dir"] = os.path.abspath(work_dirname)

        validator = SgainsValidator(
            sgains_schema, work_dirname=work_dirname)
        assert validator.validate(config_dict), validator.errors

        return Config(_dict_to_box(validator.document))

        # return Config(
        #     _dict_to_namedtuple(validator.document, "sgains"),
        #     validator)

    @property
    def schema(self):
        return sgains_schema

    def check_nonempty_workdir(self, dirname):
        if not os.path.exists(dirname):
            return
        if len(os.listdir(dirname)) and \
                not self.force and not self.dry_run:
            print(colored(
                "ERROR: non-empty output directory and no --force option",
                "red"))
            raise ValueError(f"Non empyt directory {dirname}")

    def mappable_regions_filename(self, chrom=None):
        mname = self.config.mappable_regions.mappable_file
        if chrom:
            mname = "{}_{}".format(
                chrom, self.config.mappable_regions.mappable_file)
        filename = os.path.join(
            self.mappable_regions.mappable_dir,
            mname
        )
        return filename

    def bins_boundaries_filename(self, chrom=None):
        bname = self.config.bins.bins_file
        if chrom:
            bname = "{}_{}".format(
                chrom, self.config.bins.bins_file)
        filename = os.path.join(
            self.config.bins.bins_dir,
            bname
        )
        return filename

    def __getattr__(self, attr_name):
        # FIXME Temporary hack to enable default values
        # only for public attributes
        if attr_name[0:2] == "__":
            raise AttributeError()

        if attr_name not in self.schema.keys():
            raise ValueError(f"Unexpected attribute {attr_name}")
        return getattr(self.config, attr_name)

    @staticmethod
    def cellname(filename):
        return os.path.basename(filename).split(os.extsep, 1)[0]
