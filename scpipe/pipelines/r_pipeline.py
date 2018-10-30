'''
Created on Jul 9, 2017

@author: lubo
'''
import subprocess
from config import Config, NonEmptyWorkDirectory
import os
from termcolor import colored
import shutil


class Rpipeline(object):

    def __init__(self, config):
        self.config = config

    def run(self):
        print(self.config)

        bin_boundaries_filename = os.path.abspath(
            self.config.bins_boundaries_filename())
        varbin_dir = self.config.varbin_dirname()
        varbin_suffix = self.config.varbin.varbin_suffix

        scclust_dirname = self.config.scclust_dirname()
        case_name = self.config.scclust.case_name

        cytoband = os.path.abspath(self.config.scclust.cytoband)
        assert os.path.exists(cytoband), cytoband

        nsim = self.config.scclust.nsim
        sharemin = self.config.scclust.sharemin

        print(colored(
            "processing study {}; files from {} with suffix {} to {}".format(
                case_name,
                varbin_dir,
                varbin_suffix,
                scclust_dirname),
            "green"))
        print(colored(
            "using: bin boundaries {}".format(
                bin_boundaries_filename),
            "green"))
        print(colored(
            "using: cytoband {}; nsim {}; sharemin: {}".format(
                cytoband, nsim, sharemin),
            "green"))

        if os.path.exists(scclust_dirname) and \
                len(os.listdir(scclust_dirname)) > 0:
            if self.config.force:
                if not self.config.dry_run:
                    shutil.rmtree(scclust_dirname)
                    os.makedirs(scclust_dirname)
            else:
                print(colored(
                    "results directory {} is not empty; "
                    "use --force to overwrite".format(scclust_dirname),
                    "red"))
                raise NonEmptyWorkDirectory(scclust_dirname)

        basedir = os.path.dirname(__file__)
        rscript = os.path.abspath(
            os.path.join(
                basedir,
                '../../'
                'scripts/pipeline.R'
            ))
        print(colored("executing Rscript with: {}".format(rscript), "yellow"))

        if not self.config.dry_run:
            with open(os.devnull, 'w') as shutup:

                assert os.path.exists(rscript)
                subprocess.check_call(
                    [
                        'Rscript', rscript,
                        scclust_dirname,
                        case_name,
                        varbin_dir, varbin_suffix,
                        bin_boundaries_filename,
                        cytoband,
                        str(nsim),
                        str(sharemin)
                    ],
                    shell=False,
                    # stdout=sys.stdout, stderr=sys.stdout,
                    stdout=shutup, stderr=shutup,
                )


if __name__ == "__main__":
    config = Config.load("sgains.yml")
    print(config)

    pipeline = Rpipeline(config)

    pipeline.run()
