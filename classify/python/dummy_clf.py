from sklearn.dummy import DummyClassifier
from generic_classifier import GenericClassifier
import logging


class DummyClf(GenericClassifier):

    def __init__(self, df,
                 strategy='most_frequent',
                 min_ct=10):
        self.logger = logging.getLogger(__name__)
        super(DummyClf, self).__init__()  # call base constructor
        self.set_min_count(min_ct)

        # process data
        df = self._filter_rows(df)  # filter out low count rows
        self.x, self.y = self._randomize(df)

        # setup classifier
        self.clf = DummyClassifier(strategy=strategy)
