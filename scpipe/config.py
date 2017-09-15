'''
Created on Jun 10, 2017

@author: lubo
'''
from box import Box, _get_box_config
import os
import glob
from termcolor import colored
import copy


class NonEmptyWorkDirectory(Exception):

    def __init__(self, dirname):
        super(NonEmptyWorkDirectory, self).__init__(dirname)


class Config(Box):
    DEFAULT_CONFIG = {
        "force": False,
        "dry_run": False,
        "parallel": 1,

        "genome": {
            "version": "hg19",
            "data_dir": "data/hg19_safe",
            "work_dir": "data/hg19",
            "index": "genomeindex",
        },
        "mappable_regions": {
            "length": 100,
            "work_dir": "data/R100",
            "mappable_regions": "hg19_R100_mappable_regions.txt",
            "chrom_sizes": "chrom_sizes.yml",
            "bowtie_opts": "",
        },
        "bins": {
            "bins_count": 10000,
            "work_dir": "data/R100_B10k",
            "bins_boundaries": "bins_boundaries.tst",
        },
        "mapping": {
            "data_dir": "",
            "work_dir": "",
            "bowtie_opts": "",
        },
        "varbin": {
            "data_dir": "",
            "data_glob": "*.rmdup.bam",
            "work_dir": "",
            "suffix": ".varbin.txt",
        },
        "segment": {
            "data_dir": ".",
            "data_glob": "*.varbin.txt",
            "work_dir": ".",
            "study_name": "test",
        }
    }

    def __new__(cls, *args, **kwargs):
        """
        Due to the way pickling works in python 3, we need to make sure
        the box config is created as early as possible.
        """
        obj = super(Box, cls).__new__(cls, *args, **kwargs)
        obj._box_config = _get_box_config(cls, kwargs)
        return obj

    def __init__(self, *args, **kwargs):
        super(Config, self).__init__(
            *args,
            **kwargs)

    @staticmethod
    def copy(config):
        new_box = Box(copy.deepcopy(config.to_dict()))
        return Config(new_box)

    @staticmethod
    def cellname(filename):
        return os.path.basename(filename).split(os.extsep, 1)[0]

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
    def parse_args(argv):
        if '-c' in argv:
            index = argv.index('-c')
        elif '--config' in argv:
            index = argv.index('--config')
        else:
            config = Config.load('./sgains.yml')
            return config
        index += 1
        if index < 0 or index >= len(argv):
            raise ValueError('config filename not found')

        filename = argv[index]
        config = Config.load(filename)
        return config

    @staticmethod
    def load(filename):
        default = Config.default()

        filename = os.path.abspath(filename)
        assert os.path.exists(filename), filename

        with open(filename, 'r') as infile:
            config = Box.from_yaml(infile)
            config.filename = os.path.abspath(filename)
            # config.dirname = os.path.dirname(config.filename)
            config.dirname = os.curdir

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
            self.genome.work_dir,
            self.genome.index
        )
        return self.abspath(filename)

    def mappable_regions_filename(self):
        filename = os.path.join(
            self.mappable_regions.work_dir,
            self.mappable_regions.mappable_regions
        )
        return self.abspath(filename)

    def mappable_positions_count_filename(self):
        filename = os.path.join(
            self.bins.work_dir,
            "B{}_mappable_positions_count.yaml".format(self.bins.bins_count)
        )
        return self.abspath(filename)

    def chrom_sizes_filename(self):
        filename = os.path.join(
            self.mappable_regions.work_dir,
            self.mappable_regions.chrom_sizes
        )
        filename = self.abspath(filename)
        return filename

    def bins_boundaries_filename(self):
        filename = os.path.join(
            self.bins.work_dir,
            self.bins.bins_boundaries
        )
        return self.abspath(filename)

    def chrom_filename(self, chrom, pristine=False):
        if pristine:
            cache_dir = self.genome.data_dir
        else:
            cache_dir = self.genome.work_dir
        filename = os.path.join(
            cache_dir,
            "{}.fa".format(chrom)
        )
        return self.abspath(filename)

    def varbin_data_filenames(self):
        assert os.path.exists(self.varbin.data_dir)
        dirname = self.abspath(self.varbin.data_dir)
        pattern = os.path.join(
            dirname,
            self.varbin.data_glob
        )
        filenames = glob.glob(pattern)
        return filenames

    def varbin_work_dirname(self):
        dirname = self.abspath(self.varbin.work_dir)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        return dirname

    def segment_data_dirname(self):
        dirname = self.abspath(self.segment.data_dir)
        return self.abspath(dirname)

    def segment_work_dirname(self):
        dirname = self.abspath(self.segment.work_dir)
        dirname = os.path.join(dirname, self.segment.study_name)
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        return dirname

    def varbin_work_filename(self, cellname):
        outfile = os.path.join(
            self.varbin_work_dirname(),
            "{}.{}".format(cellname, self.varbin.suffix)
        )
        outfile = self.abspath(outfile)
        return outfile

    def varbin_work_filenames(self):
        dirname = self.varbin_work_dirname()
        pattern = os.path.join(
            dirname,
            "*.{}".format(self.varbin.suffix)
        )
        return glob.glob(pattern)

    def mapping_data_dirname(self):
        assert os.path.exists(self.mapping.data_dir)
        dirname = self.abspath(self.mapping.data_dir)
        return dirname

    def mapping_fastq_filenames(self):
        pattern = os.path.join(
            self.mapping_data_dirname(),
            "*{}".format(self.mapping.data_glob)
        )
        filenames = glob.glob(pattern)
        return filenames

    def check_nonempty_workdir(self, dirname):
        if os.path.exists(dirname) and len(os.listdir(dirname)) and \
                not self.force:
            print(colored(
                "ERROR: non-empty output directory and no --force option",
                "red"))
            raise NonEmptyWorkDirectory(dirname)

    def mapping_work_dirname(self):
        dirname = self.abspath(self.mapping.work_dir)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        return dirname

    def __getstate__(self):
        state = self.to_dict()
        state['_box_config'] = copy.copy(self._box_config)

    def __setstate__(self, state):
        self._box_config = state['_box_config']
        self.__dict__.update(state)
