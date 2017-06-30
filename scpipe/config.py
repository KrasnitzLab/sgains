'''
Created on Jun 10, 2017

@author: lubo
'''
from box import Box
import os


class Config(Box):
    DEFAULT_CONFIG = {
        "genome": {
            "version": "hg19",
            "pristine": "data/hg19_safe",
            "cache_dir": "data/hg19",
            "index": "genomeindex",
        },
        "reads": {
            "length": 100,
            "cache_dir": "data/R100",
            "mappable_regions": "mappable_regions.tsv",
            "mappable_positions_count": "mappable_positions_count.yml",
            "chrom_sizes": "chrom_sizes.yml"
        },
        "bins": {
            "bins_count": 10000,
            "cache_dir": "data/R100_B10k",
            "bin_boundaries": "bin_boundaries.tst",
        }
    }

    def __init__(self, data, **kwargs):
        super(Config, self).__init__(
            data,
            **kwargs)

    @staticmethod
    def default():
        config = Box(Config.DEFAULT_CONFIG, default_box=True)
        config.filename = None
        config.dirname = os.getcwd()
        return Config(
            config.to_dict(),
            default_box=True,
        )

    @staticmethod
    def load(filename):
        default = Config.default()

        filename = os.path.abspath(filename)
        assert os.path.exists(filename), filename

        with open(filename, 'r') as infile:
            config = Box.from_yaml(infile)
            config.filename = os.path.abspath(filename)
            config.dirname = os.path.dirname(config.filename)

            default.update(config.to_dict())
            return Config(
                default.to_dict(),
                default_box=True,
            )

    def abspath(self, filename):
        return os.path.join(
            self.dirname,
            filename
        )

    def genome_index_filename(self):
        filename = os.path.join(
            self.genome.cache_dir,
            self.genome.index
        )
        return self.abspath(filename)

    def mappable_regions_filename(self):
        filename = os.path.join(
            self.reads.cache_dir,
            self.reads.mappable_regions
        )
        return self.abspath(filename)

    def mappable_positions_count_filename(self):
        filename = os.path.join(
            self.reads.cache_dir,
            self.reads.mappable_positions_count
        )
        return self.abspath(filename)

    def chrom_sizes_filename(self):
        filename = os.path.join(
            self.reads.cache_dir,
            self.reads.chrom_sizes
        )
        filename = self.abspath(filename)
        return filename

    def bin_boundaries_filename(self):
        filename = os.path.join(
            self.bins.cache_dir,
            self.bins.bin_boundaries
        )
        return self.abspath(filename)

    def chrom_filename(self, chrom, pristine=False):
        if pristine:
            cache_dir = self.genome.pristine
        else:
            cache_dir = self.genome.cache_dir
        filename = os.path.join(
            cache_dir,
            "{}.fa".format(chrom)
        )
        return self.abspath(filename)
