'''
Created on Jul 9, 2017

@author: lubo
'''
import os
import sys

import shutil
import subprocess

from termcolor import colored


class Rpipeline(object):

    def __init__(self, config):
        self.config = config

    def run(self, dask_client):
        print(self.config)

        bin_boundaries_filename = self.config.bins_boundaries_filename()
        varbin_dir = self.config.varbin.varbin_dir
        varbin_suffix = self.config.varbin.varbin_suffix

        scclust_dirname = self.config.scclust.scclust_dir
        case_name = self.config.scclust.scclust_case

        cytoband = self.config.scclust.scclust_cytoband_file
        assert os.path.exists(cytoband), cytoband

        nsim = self.config.scclust.scclust_nsim
        sharemin = self.config.scclust.scclust_sharemin
        fdrthres = self.config.scclust.scclust_fdrthres
        nshare = self.config.scclust.scclust_nshare
        climbtoshare = self.config.scclust.scclust_climbtoshare

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
            "using: cytoband {}; nsim {}; sharemin: {}; "
            "fdrthres: {}; nshare: {}; climbtoshare: {}".format(
                cytoband, nsim, sharemin, fdrthres, nshare, climbtoshare),
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
                raise ValueError(f"non empty directory {scclust_dirname}")

        basedir = os.path.dirname(__file__)
        rscript = os.path.abspath(
            os.path.join(
                basedir,
                '../'
                'scripts/pipeline.R'
            ))
        print(colored("executing Rscript with: {}".format(rscript), "yellow"))

        if not self.config.dry_run:
            with open(os.devnull, 'w') as shutup:  # noqa

                assert os.path.exists(rscript), rscript
                subprocess.check_call(
                    [
                        'Rscript', rscript,
                        scclust_dirname,
                        case_name,
                        varbin_dir, varbin_suffix,
                        bin_boundaries_filename,
                        cytoband,
                        str(nsim),
                        str(sharemin),
                        str(fdrthres),
                        str(nshare),
                        str(climbtoshare),
                    ],
                    shell=False,
                    stdout=sys.stdout, stderr=sys.stdout,
                    # stdout=shutup, stderr=shutup,
                )
