import utils.python.util as _utils
import numpy as np
from scipy import interp
from sklearn import cross_validation
import sklearn.metrics as metrics


class GenericClassifier(object):

    def __init__(self):
        self.min_count = 10  # min mutations for a gene

        # integer codes for classes
        self.other_num = 0
        self.onco_num = 1
        self.tsg_num = 1

        # genes categorized as oncogenes/tsg by vogelstein's
        # science paper
        self.vogelsteins_oncogenes = _utils.oncogene_set
        self.vogelsteins_tsg = _utils.tsg_set

    def _filter_rows(self, df):
        """Filter out rows with counts less than the minimum."""
        row_sums = df.T.sum()
        filtered_df = df[row_sums > self.min_count]
        return filtered_df

    def _label_gene(self, gene):
        """Label a gene according to Vogelstein's list of oncogenes
        and tsg."""
        if gene in self.vogelsteins_oncogenes:
            return self.onco_num
        elif gene in self.vogelsteins_tsg:
            return self.tsg_num
        else:
            return self.other_num

    def _random_sort(self, df):
        """Randomly sort a data frame.

        This method is an attempt to prevent order bias in the input data.
        Non-random ordering could effect cross validation results.

        Args:
            df (pd.DataFrame): pandas dataframe

        Returns:
            pd.DataFrame: randomly sorted data frame
        """
        random_indices = np.random.choice(df.index.values,  # sample from 'genes'
                                          len(df),  # number of samples
                                          replace=False)  # sample without replacement
        random_df = df.reindex(random_indices)  # change order of df
        return random_df

    def _randomize(self, df):
        x = self._random_sort(df)  # randomly sort data
        y = x.index.to_series().apply(self._label_gene)  # get gene labels
        return x, y

    def kfold_validation(self, k=5):
        # generate indices for kfold cross validation
        k_fold = cross_validation.KFold(n=len(self.x),  # len of df
                                        n_folds=k,  # k fold
                                        indices=True)  # return indices
        self.num_pred = 0  # number of predictions

        mean_tpr = 0.0
        mean_fpr = np.linspace(0, 1, 100)

        for i in range(k):
            self.x, self.y = self._randomize(self.x)  # randomize for another round

            # evaluate k-fold cross validation
            for train_ix, test_ix in k_fold:
                self.clf.fit(self.x.ix[train_ix], self.y.ix[train_ix])
                # print self.clf.score(self.x.ix[test_ix], self.y.ix[test_ix])

                # print 'Feature Importance', self.clf.feature_importances_
                y_pred = self.clf.predict(self.x.ix[test_ix])
                proba_ = self.clf.predict_proba(self.x.ix[test_ix])
                fpr, tpr, thresholds = metrics.roc_curve(self.y.ix[test_ix], proba_[:,1])
                mean_tpr += interp(mean_fpr, fpr, tpr)
                mean_tpr[0] = 0.0
                self._update_metrics()  # update info about classification
                #print self.clf.feature_importances_
                # print 'ROC AUC', metrics.auc(fpr, tpr)

        mean_tpr /= self.num_pred  # divide by number of folds squared
        mean_tpr[-1] = 1.0  # it always ends at 1
        mean_auc = metrics.auc(mean_fpr, mean_tpr)

        self._on_finish()  # update info for kfold cross-validation

        return mean_tpr, mean_fpr, mean_auc

    def set_min_count(self, ct):
        if ct >= 0:
            self.min_count = ct
