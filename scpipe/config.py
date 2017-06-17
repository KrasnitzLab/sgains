'''
Created on Jun 10, 2017

@author: lubo
'''
from box import Box
import os


class Config(Box):

    def __init__(self, data, **kwargs):
        super(Config, self).__init__(
            data,
            **kwargs)

    @staticmethod
    def load(filename):
        assert os.path.exists(filename)

        with open(filename, 'r') as infile:
            config = Box.from_yaml(infile)
            config.filename = os.path.abspath(filename)
            config.dirname = os.path.dirname(config.filename)

            return Config(
                config.to_dict(),
                default_box=True,
                camel_case_killer=True
            )

    def abspath(self, filename):
        return os.path.join(
            self.dirname,
            filename
        )
