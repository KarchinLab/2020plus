# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '..'))

import src.features.python.features as feat
import src.utils.python.util as _utils

def test_features():
    # setup input
    example_og = os.path.join(file_dir, 'data/og_example.txt')
    example_tsg = os.path.join(file_dir, 'data/tsg_example.txt')
    example_summary = os.path.join(file_dir, 'data/summary_example.txt')
    covar_file = os.path.join(file_dir, '../data/mutsigcv_gene_features.txt')
    biogrid_file = os.path.join(file_dir, '../data/biogrid_stats.txt')

    # test feature generation with biogrid permutation
    opts = {
        'og_test': example_og,
        'tsg_test': example_tsg,
        'summary': example_summary,
        'covariates': covar_file,
        'biogrid': biogrid_file,
        'permute_biogrid': True,
        'random_seed': 71,
        'output': 'tests/data/feature_test.txt'
    }
    feat.main(opts)

    # test feature generation
    opts['permute_biogrid'] = False
    feat.main(opts)
