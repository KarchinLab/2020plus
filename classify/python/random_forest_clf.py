from __future__ import division
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from generic_classifier import GenericClassifier
import pandas as pd
import logging


class RandomForest(GenericClassifier):

    def __init__(self, df,
                 total_iter=5,
                 weight=True,
                 min_ct=0):
        self.logger = logging.getLogger(__name__)
        super(RandomForest, self).__init__(total_iter)  # call base constructor
        # self.set_min_count(min_ct)
        self.is_weighted_sample = weight

        # process data
        #df = self._filter_rows(df)  # filter out low count rows
        #recurrent_mutation = df['recurrent missense'] + df['recurrent indel']
        #deleterious_mutation = df['lost stop'] + df['nonsense'] + df['frame shift'] + df['no protein']
        #row_sums = df.sum(axis=1).astype(float)
        #df = df.div(row_sums, axis=0)  # normalize each row
        #df['recurrent count'] = recurrent_mutation
        #df['deleterious count'] = deleterious_mutation
        #df.to_csv('tmp.rclf.txt', sep='\t')
        df = df.fillna(df.mean())
        self.x, self.y = self._randomize(df)

        # setup classifier
        self.clf = RandomForestClassifier(n_estimators=250)
        #self.clf = ExtraTreesClassifier(n_estimators=1000,
                                        #max_features=2,
                                        #n_jobs=4)

    def _update_onco_metrics(self, y_true, y_pred, prob):
        super(RandomForest, self)._update_onco_metrics(y_true, y_pred, prob)

        # evaluate feature importance for random forest
        self.feature_importance.append(self.clf.feature_importances_)

    def _update_tsg_metrics(self, y_true, y_pred, prob):
        super(RandomForest, self)._update_tsg_metrics(y_true, y_pred, prob)

        # evaluate feature importance for random forest
        self.feature_importance.append(self.clf.feature_importances_)

    def _on_finish(self):
        super(RandomForest, self)._on_finish()
        self.feature_importance = pd.DataFrame(self.feature_importance,
                                               columns=self.x.columns)
        self.mean_importance = self.feature_importance.mean()
        self.std_importance = self.feature_importance.std()
