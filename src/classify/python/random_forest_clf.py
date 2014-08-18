from __future__ import division
from sklearn.ensemble import RandomForestClassifier
from generic_classifier import GenericClassifier
import features.python.features as features
import pandas as pd
import logging


class RandomForest(GenericClassifier):

    def __init__(self, df,
                 ntrees=250,
                 total_iter=5,
                 weight=True,
                 min_ct=0):
        self.logger = logging.getLogger(__name__)
        super(RandomForest, self).__init__(total_iter)  # call base constructor
        self.is_weighted_sample = weight

        # process data
        if 'total' in df.columns:
            df = df.drop('total', axis=1)
        df = df.fillna(df.mean())
        #myfeatures = {'nonsense': df['nonsense'],
                      #'no protein': df['no protein'],
                      #'lost stop': df['lost stop'],
                      #'frame shift': df['frame shift'],
                       #'deleterious percent': df['nonsense'] + df['no protein'] + df['lost stop'] + df['frame shift'],
                      #'recurrent missense': df['recurrent missense'],
                      #'recurrent indel': df['recurrent indel'],
                       #'recurrent percent': df['recurrent missense'] + df['recurrent indel'],
                      #'missense': df.missense,
                      #'indel': df.indel,
                      #'recurrent count': df['recurrent count'],
                      #'deleterious count': df['deleterious count'],
                      #'synonymous': df.synonymous,
                      #'missense position entropy': df['missense position entropy'],
                      #'mutation position entropy': df['mutation position entropy']}
        #tmpdf = pd.DataFrame(myfeatures)

        # optional columns from mutsigcv paper
        #if 'gene_length' in df.columns:
            #tmpdf['gene_length'] = df['gene_length']
        #if 'noncoding_mutation_rate' in df.columns:
            #tmpdf['noncoding_mutation_rate'] = df['noncoding_mutation_rate']
        #if 'expression' in df.columns:
            #tmpdf['expression'] = df['expression']
        #if 'replication_time' in df.columns:
            #tmpdf['replication_time'] = df['replication_time']

        self.x, self.y = features.randomize(df)

        # setup classifier
        self.clf = RandomForestClassifier(n_estimators=ntrees)

    def _update_metrics(self, y_true, y_pred, onco_prob, tsg_prob):
        super(RandomForest, self)._update_metrics(y_true, y_pred, onco_prob, tsg_prob)

        # evaluate feature importance for random forest
        self.feature_importance.append(self.clf.feature_importances_)

    def _update_onco_metrics(self, y_true, y_pred, prob):
        super(RandomForest, self)._update_onco_metrics(y_true, y_pred, prob)

        # evaluate feature importance for random forest
        # self.feature_importance.append(self.clf.feature_importances_)

    def _update_tsg_metrics(self, y_true, y_pred, prob):
        super(RandomForest, self)._update_tsg_metrics(y_true, y_pred, prob)

        # evaluate feature importance for random forest
        # self.feature_importance.append(self.clf.feature_importances_)

    def _on_finish(self):
        super(RandomForest, self)._on_finish()
        self.feature_importance = pd.DataFrame(self.feature_importance,
                                               columns=self.x.columns)
        self.mean_importance = self.feature_importance.mean()
        self.std_importance = self.feature_importance.std()
