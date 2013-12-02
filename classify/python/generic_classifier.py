import utils.python.util as _utils
import numpy as np
from numpy import interp
from scipy.stats import itemfreq
from sklearn import cross_validation
import sklearn.metrics as metrics
import matplotlib.pyplot as plt


class GenericClassifier(object):

    def __init__(self):
        self.min_count = 10  # min mutations for a gene

        # integer codes for classes
        self.other_num = 0
        self.onco_num = 1
        self.tsg_num = 0
        self.num_classes = 2  # total number of classes to predict

        self._init_metrics()  # initialize metrics

        # genes categorized as oncogenes/tsg by vogelstein's
        # science paper
        self.vogelsteins_oncogenes = _utils.oncogene_set
        self.vogelsteins_tsg = _utils.tsg_set

    def _init_metrics(self):
        self.feature_importance = []
        self.confusion_matrix = np.zeros((self.num_classes, self.num_classes))
        self.f1_score = 0.0
        self.precision = 0.0
        self.recall = 0.0
        self.num_pred = 0  # no predictions have been made yet
        self.mean_tpr = 0.0
        self.mean_fpr = np.linspace(0, 1, 100)
        self.mean_precision = 0.0
        #self.mean_recall = 0.0
        self.mean_recall = np.linspace(0, 1, 100)

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

    def _update_metrics(self, y_true, y_pred, prob):
        self.num_pred += 1

        # estabilish confusion matrix
        tmp_confusion_matrix = metrics.confusion_matrix(y_true, y_pred)
        self.confusion_matrix += tmp_confusion_matrix

        # compute metrics for classification
        prec, recall, fscore, support = metrics.precision_recall_fscore_support(y_true, y_pred)
        self.precision += prec
        self.recall += recall
        self.f1_score += fscore
        self.logger.debug('Iter %d: Precission=%s, Recall=%s, f1_score=%s' % (
                          self.num_pred, str(prec), str(recall), str(fscore)))

        # compute ROC curve metrics
        fpr, tpr, thresholds = metrics.roc_curve(y_true, prob)
        self.mean_tpr += interp(self.mean_fpr, fpr, tpr)
        self.mean_tpr[0] = 0.0

        # compute Precision-Recall curve metrics
        p, r, thresh = metrics.precision_recall_curve(y_true, prob)
        p, r = p[::-1], r[::-1]  # reverse order of results
        plt.plot(r, p)
        plt.savefig('tmp.pr.png')
        self.mean_precision += interp(self.mean_recall, r, p)
        #print interp(self.mean_recall, r, p)
        #print p
        #print r
        # self.mean_precision += p
        # self.mean_recall += r

    def _on_finish(self):
        self.confusion_matrix /= self.num_pred
        self.f1_score /= self.num_pred
        self.precision /= self.num_pred
        self.recall /= self.num_pred

        # ROC curve metrics
        self.mean_tpr /= self.num_pred  # divide by number of folds squared
        self.mean_tpr[-1] = 1.0  # it always ends at 1
        self.mean_roc_auc = metrics.auc(self.mean_fpr, self.mean_tpr)

        # Precision-Recall curve metrics
        self.mean_precision /= self.num_pred
        # self.mean_recall /= self.num_pred
        self.mean_pr_auc = metrics.auc(self.mean_recall, self.mean_precision)

    def kfold_validation(self, k=5):
        # generate indices for kfold cross validation
        k_fold = cross_validation.KFold(n=len(self.x),  # len of df
                                        n_folds=k,  # k fold
                                        indices=True)  # return indices
        self.num_pred = 0  # number of predictions

        for i in range(k):
            self.x, self.y = self._randomize(self.x)  # randomize for another round

            # evaluate k-fold cross validation
            for train_ix, test_ix in k_fold:
                self.clf.fit(self.x.ix[train_ix], self.y.ix[train_ix])
                y_pred = self.clf.predict(self.x.ix[test_ix])
                proba_ = self.clf.predict_proba(self.x.ix[test_ix])
                self._update_metrics(self.y.ix[test_ix], y_pred, proba_[:,1])  # update info about classification

        self._on_finish()  # update info for kfold cross-validation

    def get_roc_metrics(self):
        return self.mean_tpr, self.mean_fpr, self.mean_roc_auc

    def get_pr_metrics(self):
        return self.mean_precision, self.mean_recall, self.mean_pr_auc

    def set_min_count(self, ct):
        if ct >= 0:
            self.min_count = ct
