# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '..'))

import src.train.python.train as train
import src.utils.python.util as _utils


def test_training():
    """Test the training of a random forest model."""
    example_features = os.path.join(file_dir, 'data/features_pancan_subset.txt')
    test_train = os.path.join(file_dir, 'data/test_train.Rdata')
    opts = {
        'features': example_features,
        'driver_rate': .7,
        'other_ratio': 1.0,
        'ntrees': 200,
        'min_count': 0,
        'random_seed': 71,
        'output': test_train
    }
    train.main(opts)
