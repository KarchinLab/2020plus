from sklearn.naive_bayes import MultinomialNB
import sklearn.metrics as metrics
from generic_classifier import GenericClassifier
import numpy as np
import logging

class MultinomialNaiveBayes(GenericClassifier):

    def __init__(self, df, min_ct=10):
        self.logger = logging.getLogger(__name__)
        super(MultinomialNaiveBayes, self).__init__()  # call base constructor
        self.set_min_count(min_ct)

        # process data
        df = self._filter_rows(df)  # filter out low count rows
        row_sums = df.sum(axis=1).astype(float)
        df = df.div(row_sums, axis=0)  # normalize each row
        # df = df.mul(100)
        df.to_csv('tmp.nbclf.txt', sep='\t')
        self.x, self.y = self._randomize(df)

        # setup classifier
        self.clf = MultinomialNB(alpha=1,  # laplacian smooth, i.e. pseudocounts
                                 fit_prior=True)  # use data for prior class probs

