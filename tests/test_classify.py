# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '..'))

import src.classify.python.classifier as clf

def test_trained_classifier():
    trained_clf = os.path.join(file_dir, 'data/pancan.Rdata')
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
