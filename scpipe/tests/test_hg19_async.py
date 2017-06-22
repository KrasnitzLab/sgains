'''
Created on Jun 22, 2017

@author: lubo
'''


def test_mappings_generator1(hg):
    reads_generator = hg.generate_reads(['chrM'], 100)
    hg.mappings_generator1(reads_generator)
