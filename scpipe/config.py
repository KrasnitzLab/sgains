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

    def mappable_regions_filename(self):
        filename = os.path.join(
            self.bins.cache_dir,
            self.bins.mappable_regions
        )
        return self.abspath(filename)

    def mappable_positions_count_filename(self):
        filename = os.path.join(
            self.bins.cache_dir,
            self.bins.mappable_positions_count
        )
        return self.abspath(filename)

    def chrom_sizes_filename(self):
        filename = os.path.join(
            self.bins.cache_dir,
            self.bins.chrom_sizes
        )
        filename = self.abspath(filename)
        return filename
