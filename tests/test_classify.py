# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '..'))

import src.classify.python.classifier as clf
import src.utils.python.util as _utils

def test_trained_classifier():
    trained_clf = os.path.join(file_dir, 'data/test_train.Rdata')
    example_features = os.path.join(file_dir, 'data/example_sim_features.txt')
    null_dist = os.path.join(file_dir, 'data/example_null_dist.txt')
    opts = {
        'trained_classifier': trained_clf,
        'features': example_features,
        'null_distribution': null_dist,
        'simulated': True,
        'min_count': 0,
        'driver_rate': .7,
        'other_ratio': 1.0,
        'ntrees': 500,
        'random_seed': 71
    }
    clf.main(opts)


def test_regular_classifier():
    example_features = os.path.join(file_dir, 'data/features_pancan_subset.txt')
    null_dist = os.path.join(file_dir, 'data/example_null_dist.txt')
    opts = {
        'trained_classifier': None,
        'simulated': False,
        'features': example_features,
        'min_count': 0,
        'driver_rate': .7,
        'null_distribution': None,
        'other_ratio': 1.0,
        'ntrees': 200,
        'random_seed': 71
    }
    # setup directory
    _utils.make_result_dir('tests/data/result')

    # classify without empirical null
    clf.main(opts)

    # now with an empirical null
    opts['null_distribution'] = null_dist
    clf.main(opts)
