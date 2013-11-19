from sklearn.ensemble import RandomForestClassifier
from generic_classifier import GenericClassifier


class RandomForest(GenericClassifier):

    def __init__(self, df,
                 min_ct=10):
        super(RandomForest, self).__init__()  # call base constructor
        self.set_min_count(min_ct)

        # process data
        df = self._filter_rows(df)  # filter out low count rows
        self.x = self._random_sort(df)  # randomly sort data
        # row_sums = self.x.sum(axis=1)
        # self.x = self.x.div(row_sums, axis=0)  # normalize each row
        # self.x['total'] = row_sums  # add a column for total counts
        self.y = self.x.index.to_series().apply(self._label_gene)  # get gene labels

        # setup classifier
        self.clf = RandomForestClassifier(n_estimators=100)
