from sklearn.dummy import DummyClassifier
from generic_classifier import GenericClassifier
import features.python.features as features
import logging


class DummyClf(GenericClassifier):

    def __init__(self, df,
                 strategy='most_frequent',
                 weight=False,
                 min_ct=0):
        self.logger = logging.getLogger(__name__)
        super(DummyClf, self).__init__()  # call base constructor
        #self.set_min_count(min_ct)
        self.is_weighted_sample = False

        # process data
        #df = self._filter_rows(df)  # filter out low count rows
        df = df.fillna(df.mean())
        self.x, self.y = features.randomize(df)

        # setup classifier
        self.clf = DummyClassifier(strategy=strategy)
