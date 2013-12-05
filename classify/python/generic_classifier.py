import utils.python.util as _utils
import numpy as np
from numpy import interp
from sklearn import cross_validation
import sklearn.metrics as metrics
import pandas as pd


class GenericClassifier(object):

    def __init__(self):
        self.min_count = 0  # min mutations for a gene

        # set integer codes for classes
        self.set_classes(oncogene=True, tsg=True)

        self._init_metrics()  # initialize metrics

        # genes categorized as oncogenes/tsg by vogelstein's
        # science paper
        self.vogelsteins_oncogenes = _utils.oncogene_set
        self.vogelsteins_tsg = _utils.tsg_set

    def _init_metrics(self):
        """Initialize classification diagnostric metrics."""
        self.feature_importance = []
        self.confusion_matrix = np.zeros((self.num_classes, self.num_classes))
        self.num_pred = 0  # no predictions have been made yet

        # oncogene metrics
        self.onco_f1_score = 0.0
        self.onco_precision = 0.0
        self.onco_recall = 0.0
        self.onco_mean_tpr = 0.0
        self.onco_mean_fpr = np.linspace(0, 1, 100)
        self.onco_mean_precision = 0.0
        self.onco_mean_recall = np.linspace(0, 1, 100)

        # tsg metrics
        self.tsg_f1_score = 0.0
        self.tsg_precision = 0.0
        self.tsg_recall = 0.0
        self.tsg_mean_tpr = 0.0
        self.tsg_mean_fpr = np.linspace(0, 1, 100)
        self.tsg_mean_precision = 0.0
        self.tsg_mean_recall = np.linspace(0, 1, 100)

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

    def _update_onco_metrics(self, y_true, y_pred, prob):

        # estabilish confusion matrix
        #tmp_confusion_matrix = metrics.confusion_matrix(y_true, y_pred)
        #self.confusion_matrix += tmp_confusion_matrix

        # compute metrics for classification
        prec, recall, fscore, support = metrics.precision_recall_fscore_support(y_true, y_pred)
        self.onco_precision += prec
        self.onco_recall += recall
        self.onco_f1_score += fscore
        self.logger.debug('Onco Iter %d: Precission=%s, Recall=%s, f1_score=%s' % (
                          self.num_pred, str(prec), str(recall), str(fscore)))

        # compute ROC curve metrics
        fpr, tpr, thresholds = metrics.roc_curve(y_true, prob)
        self.onco_mean_tpr += interp(self.onco_mean_fpr, fpr, tpr)
        self.onco_mean_tpr[0] = 0.0

        # compute Precision-Recall curve metrics
        p, r, thresh = metrics.precision_recall_curve(y_true, prob)
        p, r, thresh = p[::-1], r[::-1], thresh[::-1]  # reverse order of results
        self.onco_mean_precision += interp(self.onco_mean_recall, r, p)

    def _update_tsg_metrics(self, y_true, y_pred, prob):

        # estabilish confusion matrix
        #tmp_confusion_matrix = metrics.confusion_matrix(y_true, y_pred)
        #self.confusion_matrix += tmp_confusion_matrix

        # compute metrics for classification
        prec, recall, fscore, support = metrics.precision_recall_fscore_support(y_true, y_pred)
        self.tsg_precision += prec
        self.tsg_recall += recall
        self.tsg_f1_score += fscore
        self.logger.debug('Tsg Iter %d: Precission=%s, Recall=%s, f1_score=%s' % (
                          self.num_pred, str(prec), str(recall), str(fscore)))

        # compute ROC curve metrics
        fpr, tpr, thresholds = metrics.roc_curve(y_true, prob)
        self.tsg_mean_tpr += interp(self.tsg_mean_fpr, fpr, tpr)
        self.tsg_mean_tpr[0] = 0.0

        # compute Precision-Recall curve metrics
        p, r, thresh = metrics.precision_recall_curve(y_true, prob)
        p, r, thresh = p[::-1], r[::-1], thresh[::-1]  # reverse order of results
        self.tsg_mean_precision += interp(self.tsg_mean_recall, r, p)

    def _on_finish(self):
        self.confusion_matrix /= self.num_pred
        self.onco_f1_score /= self.num_pred
        self.onco_precision /= self.num_pred
        self.onco_recall /= self.num_pred
        self.tsg_f1_score /= self.num_pred
        self.tsg_precision /= self.num_pred
        self.tsg_recall /= self.num_pred

        # ROC curve metrics
        self.onco_mean_tpr /= self.num_pred  # divide by number of folds squared
        self.onco_mean_tpr[-1] = 1.0  # it always ends at 1
        self.onco_mean_roc_auc = metrics.auc(self.onco_mean_fpr, self.onco_mean_tpr)
        self.tsg_mean_tpr /= self.num_pred  # divide by number of folds squared
        self.tsg_mean_tpr[-1] = 1.0  # it always ends at 1
        self.tsg_mean_roc_auc = metrics.auc(self.tsg_mean_fpr, self.tsg_mean_tpr)

        # Precision-Recall curve metrics
        self.onco_mean_precision /= self.num_pred
        self.onco_mean_pr_auc = metrics.auc(self.onco_mean_recall, self.onco_mean_precision)
        self.tsg_mean_precision /= self.num_pred
        self.tsg_mean_pr_auc = metrics.auc(self.tsg_mean_recall, self.tsg_mean_precision)

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
                # do prediction
                self.clf.fit(self.x.iloc[train_ix], self.y.iloc[train_ix])
                y_pred = self.clf.predict(self.x.iloc[test_ix])
                proba_ = self.clf.predict_proba(self.x.iloc[test_ix])

                # update information
                self.num_pred += 1
                true_onco = (self.y.iloc[test_ix]==self.onco_num).astype(int)  # true oncogenes
                pred_onco = (y_pred==self.onco_num).astype(int)  # predicted oncogenes
                self._update_onco_metrics(true_onco,
                                          pred_onco,
                                          proba_[:, self.onco_num])
                true_tsg = (self.y.iloc[test_ix]==self.tsg_num).astype(int)  # true oncogenes
                pred_tsg = (y_pred==self.tsg_num).astype(int)  # predicted oncogenes
                self._update_tsg_metrics(true_tsg,
                                         pred_tsg,
                                         proba_[:, self.tsg_num])

        self._on_finish()  # update info for kfold cross-validation

    def kfold_prediction(self, k=5):
        # generate indices for kfold cross validation
        k_fold = cross_validation.KFold(n=len(self.x),  # len of df
                                        n_folds=k,  # k fold
                                        indices=True)  # return indices
        self.num_pred = 0  # number of predictions
        self.x, self.y = self._randomize(self.x)  # randomize data

        prediction = pd.Series(index=self.y.index)  # predicted class
        onco_prob = pd.Series(index=self.y.index).fillna(0)
        tsg_prob = pd.Series(index=self.y.index).fillna(0)

        for i in range(k):
            self.x, self.y = self._randomize(self.x)  # randomize for another round
            # obtain predictions from single round of kfold validation
            for train_ix, test_ix in k_fold:
                # retreive indices from pandas dataframe using row number
                tmp_train_ix = self.x.iloc[train_ix].index
                tmp_test_ix = self.x.iloc[test_ix].index

                # predict test data in kfold validation
                self.clf.fit(self.x.ix[tmp_train_ix], self.y.ix[tmp_train_ix])
                tmp_prob = self.clf.predict_proba(self.x.ix[tmp_test_ix])
                onco_prob.ix[tmp_test_ix] += tmp_prob[:, self.onco_num]
                tsg_prob.ix[tmp_test_ix] += tmp_prob[:, self.tsg_num]
            self.num_pred += 1
        onco_prob /= self.num_pred
        tsg_prob /= self.num_pred

        # return prediction.astype(int), prob
        return onco_prob, tsg_prob

    def get_onco_roc_metrics(self):
        return self.onco_mean_tpr, self.onco_mean_fpr, self.onco_mean_roc_auc

    def get_tsg_roc_metrics(self):
        return self.tsg_mean_tpr, self.tsg_mean_fpr, self.tsg_mean_roc_auc

    def get_onco_pr_metrics(self):
        return self.onco_mean_precision, self.onco_mean_recall, self.onco_mean_pr_auc

    def get_tsg_pr_metrics(self):
        return self.tsg_mean_precision, self.tsg_mean_recall, self.tsg_mean_pr_auc

    def set_classes(self, oncogene, tsg):
        """Sets the integers used to represent classes in classification."""
        if not oncogene and not tsg:
            raise ValueError('Classification needs at least two classes')
        elif oncogene and tsg:
            self.other_num = 0
            self.onco_num = 1
            self.tsg_num = 2
            self.num_classes = 3
        else:
            self.other_num = 0
            self.num_classes = 2
            self.onco_num = 1 if oncogene else 0
            self.tsg_num = 1 if tsg else 0

    def set_min_count(self, ct):
        if ct >= 0:
            self.min_count = ct
