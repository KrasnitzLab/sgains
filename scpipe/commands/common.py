'''
Created on Aug 2, 2017

@author: lubo
'''


class OptionsBase(object):

    def __init__(self, config):
        self.config = config
        self.subconfig = None
        self.parser = None

    def common_options(self):
        pass


class DataDirMixin(object):
    __slots__ = ()

    def data_dir_options(self, glob=False):
        assert self.subconfig is not None
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "input data options")
        group.add_argument(
            "--data-dir", "-i",
            dest="data_dir",
            help="input data directory where the input data is located",
            default=self.subconfig.data_dir
        )
        if glob:
            group.add_argument(
                "--glob", "-g",
                dest="data_glob",
                help="glob pattern for finding input data",
                default=self.subconfig.data_glob)

        return group

    def data_dir_update(self, args, glob=False):
        assert self.subconfig is not None
        assert self.subparser is not None

        if args.data_dir is not None:
            self.subconfig.data_dir = args.data_dir
        if glob:
            if args.data_glob is not None:
                self.subconfig.data_glob = args.data_glob


class WorkDirMixin(object):
    __slots__ = ()

    def work_dir_options(self):
        assert self.subconfig is not None
        assert self.subparser is not None

        group = self.subparser.add_argument_group(
            "output data options")
        group.add_argument(
            "--work-dir", "-o",
            dest="work_dir",
            help="output directory where results from processing are stored",
            default=self.subconfig
        )
        return group

    def work_dir_update(self, args):
        if args.work_dir is not None:
            self.subconfig.work_dir = args.work_dir
