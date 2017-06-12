'''
Created on Jun 10, 2017

@author: lubo
'''
from box import Box
import os


def load_config(filename):
    assert os.path.exists(filename)

    with open(filename, 'r') as infile:
        config = Box.from_yaml(infile)
        config.filename = os.path.abspath(filename)
        return config
