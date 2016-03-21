# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '..'))

def test_run():
    import run
    import src.classify.python.classifier
    import src.features.python.features
    import src.savedb.python.gene_tsv
    import src.savedb.python.gene_features
    import src.savedb.python.gene_maf
    import src.savedb.python.merge_mutations
    import src.train.python.train
    import rpy2
    import pandas
