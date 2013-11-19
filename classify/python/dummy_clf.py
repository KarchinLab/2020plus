from sklearn.dummy import DummyClassifier
from generic_classifier import GenericClassifier

class DummyClf(GenericClassifier):

    def __init__(self, df,
                 strategy='most_frequent'):
        super(DummyClf, self).__init__()  # call base constructor

        # process data
        df = self._filter_rows(df)  # filter out low count rows
        self.x = self._random_sort(df)  # randomly sort data
        self.y = self.x.index.to_series().apply(self._label_gene)  # get gene labels

        # setup classifier
        self.clf = DummyClassifier(strategy=strategy)
