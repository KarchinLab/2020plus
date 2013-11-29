from __future__ import division
from sklearn.ensemble import RandomForestClassifier
from generic_classifier import GenericClassifier
import pandas as pd


class RandomForest(GenericClassifier):

    def __init__(self, df,
                 min_ct=10):
        super(RandomForest, self).__init__()  # call base constructor
        self.set_min_count(min_ct)

        self.feature_importance = []
        self.num_pred = 0  # no predictions have been made yet

        # process data
        df = self._filter_rows(df)  # filter out low count rows
        recurrent_mutation = df['recurrent missense'] + df['recurrent indel']
        deleterious_mutation = df['lost stop'] + df['nonsense'] + df['frame shift'] + df['no protein']
        row_sums = df.sum(axis=1).astype(float)
        df = df.div(row_sums, axis=0)  # normalize each row
        df['recurrent count'] = recurrent_mutation
        df['deleterious count'] = deleterious_mutation
        df.to_csv('tmp.rclf.txt', sep='\t')
        #self.x = self._random_sort(df)  # randomly sort data
        # self.x['total'] = row_sums  # add a column for total counts
        #self.y = self.x.index.to_series().apply(self._label_gene)  # get gene labels
        self.x, self.y = self._randomize(df)

        # setup classifier
        self.clf = RandomForestClassifier(n_estimators=100)

    def _update_metrics(self):
        self.feature_importance.append(self.clf.feature_importances_)
        self.num_pred += 1

    def _on_finish(self):
        self.feature_importance = pd.DataFrame(self.feature_importance,
                                               columns=self.x.columns)
        self.mean_importance = self.feature_importance.mean()
        self.std_importance = self.feature_importance.std()
