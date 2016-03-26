from sklearn.naive_bayes import MultinomialNB
from src.classify.python.generic_classifier import GenericClassifier
import src.features.python.features as features
import logging

class MultinomialNaiveBayes(GenericClassifier):

    def __init__(self, df, weight=True, min_ct=0, total_iter=5):
        self.logger = logging.getLogger(__name__)
        super(MultinomialNaiveBayes, self).__init__(total_iterations=total_iter)  # call base constructor
        #self.set_min_count(min_ct)
        self.is_weighted_sample = weight

        # process data
        #df = self._filter_rows(df)  # filter out low count rows
        # row_sums = df.sum(axis=1).astype(float)
        # df = df.div(row_sums, axis=0)  # normalize each row
        # df = df.mul(100)
        # df.to_csv('tmp.nbclf.txt', sep='\t')
        df = df.fillna(df.mean())
        total = df['total']
        df = df[['recurrent missense', 'recurrent indel', 'frame shift',
                 'nonsense', 'missense', 'synonymous', 'inframe indel', 'no protein',
                 'lost stop', 'splicing mutation']]
        df = df.mul(total, axis=0).astype(int)  # get back counts instead of pct
        self.x, self.y = features.randomize(df)

        # setup classifier
        self.clf = MultinomialNB(alpha=1,  # laplacian smooth, i.e. pseudocounts
                                 fit_prior=True)  # use data for prior class probs

