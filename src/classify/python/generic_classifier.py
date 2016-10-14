"""The generic_classifier module contains code to do validation
of the classifier predictions."""

import src.utils.python.util as _utils
import src.utils.python.math as mymath
import src.features.python.feature_utils as futils
import numpy as np
from numpy import interp
from sklearn import cross_validation
import sklearn.metrics as metrics
import pandas as pd


class GenericClassifier(object):

    def __init__(self,
                 total_iterations=5,
                 classify_oncogene=True,
                 classify_tsg=True,
                 rseed=None):
        self.min_count = 0  # min mutations for a gene

        # set the the state of the prng
        self.prng = np.random.RandomState(rseed)
        self.rseed = rseed

        # set integer codes for classes
        self.set_classes(oncogene=classify_oncogene, tsg=classify_tsg)

        self.set_total_iter(total_iterations)  # initialize total number of simulations
        self._init_metrics()  # initialize metrics

        # genes categorized as oncogenes/tsg by vogelstein's
        # science paper
        self.vogelsteins_oncogenes = _utils.oncogene_set
        self.vogelsteins_tsg = _utils.tsg_set


    def set_total_iter(self, myiterations):
        self.total_iter = myiterations

    def _filter_rows(self, df):
        """Filter out rows with counts less than the minimum."""
        row_sums = df.T.sum()
        filtered_df = df[row_sums > self.min_count]
        return filtered_df

    def train(self):
        """Train classifier on entire data set provided."""
        self.x, self.y = futils.randomize(self.x, self.prng)
        self.clf.fit(self.x, self.y)

    def train_cv(self, k=10):
        """Train classifier on entire data set provided, but done in cross-validation."""
        # generate indices for kfold cross validation
        self.num_pred = 0  # number of predictions

        for i in range(self.total_iter):
            # randomize for another round
            self.x, self.y = futils.randomize(self.x, self.prng)

            # set up stratified kfold iterator
            k_fold = cross_validation.StratifiedKFold(self.y,
                                                      n_folds=k)

            # obtain predictions from single round of kfold validation
            for train_ix, test_ix in k_fold:
                # retreive indices from pandas dataframe using row number
                tmp_train_ix = self.x.iloc[train_ix].index

                if self.is_weighted_sample:
                    # figure out sample weights
                    num_train = len(train_ix)
                    sample_weight = np.zeros(num_train)
                    onco_ix = np.nonzero(self.y.ix[tmp_train_ix]==self.onco_num)[0]
                    tsg_ix = np.nonzero(self.y.ix[tmp_train_ix]==self.tsg_num)[0]
                    other_ix = np.nonzero(self.y.ix[tmp_train_ix]==self.other_num)[0]
                    sample_weight[onco_ix] = 1. / len(onco_ix)
                    sample_weight[tsg_ix] = 1. / len(tsg_ix)
                    sample_weight[other_ix] = 1. / len(other_ix)

                    # do training with sample weighting
                    self.clf.fit(self.x.ix[tmp_train_ix].copy(),
                                 self.y.ix[tmp_train_ix].copy(),
                                 sample_weight=sample_weight)
                else:
                    # do training without weighting
                    self.clf.fit(self.x.ix[tmp_train_ix].copy(),
                                 self.y.ix[tmp_train_ix].copy())
                self.clf.append_fold_result()  # add the training result from each fold
            self.clf.append_cv_result()  # add the training result for a single CV to the R variable

            self.num_pred += 1


    def predict(self):
        """Predict after using the train method."""
        self.num_pred = 1  # only one prediction

        # initialize predicted results variables
        num_genes = len(self.y)

        # do predictions
        y_pred = self.clf.predict(self.x)
        proba_ = self.clf.predict_proba(self.x)

        # update information
        onco_prob = proba_[:, self.onco_num] # predicted oncogenes
        tsg_prob = proba_[:, self.tsg_num] # predicted oncogenes
        other_prob = 1 - (onco_prob + tsg_prob)

        return onco_prob, tsg_prob, other_prob

    def kfold_validation(self, k=10):
        """Records the performance in terms of ROC and PR AUC for cross-validation.

        Params
        ------
        k : int (10)
            Number of cross-validation folds
        """
        self.num_pred = 0  # number of predictions

        for i in range(self.total_iter):
            # randomize for another round
            self.x, self.y = futils.randomize(self.x, self.prng)

            # initialize predicted results variables
            num_genes = len(self.y)
            onco_pred = np.zeros(num_genes)
            onco_prob = np.zeros(num_genes)
            tsg_pred = np.zeros(num_genes)
            tsg_prob = np.zeros(num_genes)
            overall_pred = np.zeros(num_genes)

            # set up stratified kfold iterator
            k_fold = cross_validation.StratifiedKFold(self.y,
                                                      n_folds=k)

            # evaluate k-fold cross validation
            for train_ix, test_ix in k_fold:
                if self.is_weighted_sample:
                    # weight classes by using sample weights
                    num_train = len(train_ix)
                    sample_weight = np.zeros(num_train)
                    onco_ix = np.nonzero(self.y[train_ix]==self.onco_num)[0]
                    tsg_ix = np.nonzero(self.y[train_ix]==self.tsg_num)[0]
                    other_ix = np.nonzero(self.y[train_ix]==self.other_num)[0]
                    sample_weight[onco_ix] = 1. / len(onco_ix)
                    sample_weight[tsg_ix] = 1. / len(tsg_ix)
                    sample_weight[other_ix] = 1. / len(other_ix)

                    # do training
                    self.clf.fit(self.x.iloc[train_ix].copy(),
                                 self.y.iloc[train_ix].copy(),
                                 sample_weight=sample_weight)
                else:
                    # do training without sample weights
                    self.clf.fit(self.x.iloc[train_ix].copy(),
                                 self.y.iloc[train_ix].copy())

                # do prediction
                y_pred = self.clf.predict(self.x.iloc[test_ix])
                proba_ = self.clf.predict_proba(self.x.iloc[test_ix])

                # update information
                overall_pred[test_ix] = y_pred  # prediction including all classes
                onco_pred[test_ix] = (y_pred==self.onco_num).astype(int)  # predicted oncogenes
                onco_prob[test_ix] = proba_[:, self.onco_num] # predicted oncogenes
                tsg_pred[test_ix] = (y_pred==self.tsg_num).astype(int)  # predicted oncogenes
                tsg_prob[test_ix] = proba_[:, self.tsg_num] # predicted oncogenes

            # update information
            true_onco = (self.y==self.onco_num).astype(int)
            self._update_onco_metrics(true_onco,
                                      onco_pred,
                                      onco_prob)
            true_tsg = (self.y==self.tsg_num).astype(int)  # true oncogenes
            self._update_tsg_metrics(true_tsg,
                                     tsg_pred,
                                     tsg_prob)
            self._update_metrics(self.y,
                                 overall_pred,
                                 onco_prob,
                                 tsg_prob)
            self.num_pred += 1

        self._on_finish()  # update info for kfold cross-validation

    def kfold_prediction(self, k=10):
        # generate indices for kfold cross validation
        self.num_pred = 0  # number of predictions

        prediction = pd.Series(index=self.y.index)  # predicted class
        onco_prob = pd.Series(index=self.y.index).fillna(0)
        tsg_prob = pd.Series(index=self.y.index).fillna(0)

        for i in range(self.total_iter):
            # randomize for another round
            self.x, self.y = futils.randomize(self.x, self.prng)

            # set up stratified kfold iterator
            k_fold = cross_validation.StratifiedKFold(self.y,
                                                      n_folds=k)

            # obtain predictions from single round of kfold validation
            for train_ix, test_ix in k_fold:
                # retreive indices from pandas dataframe using row number
                tmp_train_ix = self.x.iloc[train_ix].index
                tmp_test_ix = self.x.iloc[test_ix].index

                if self.is_weighted_sample:
                    # figure out sample weights
                    num_train = len(train_ix)
                    sample_weight = np.zeros(num_train)
                    onco_ix = np.nonzero(self.y.ix[tmp_train_ix]==self.onco_num)[0]
                    tsg_ix = np.nonzero(self.y.ix[tmp_train_ix]==self.tsg_num)[0]
                    other_ix = np.nonzero(self.y.ix[tmp_train_ix]==self.other_num)[0]
                    sample_weight[onco_ix] = 1. / len(onco_ix)
                    sample_weight[tsg_ix] = 1. / len(tsg_ix)
                    sample_weight[other_ix] = 1. / len(other_ix)

                    # do training with sample weighting
                    self.clf.fit(self.x.ix[tmp_train_ix].copy(),
                                 self.y.ix[tmp_train_ix].copy(),
                                 sample_weight=sample_weight)
                else:
                    # do training without weighting
                    self.clf.fit(self.x.ix[tmp_train_ix].copy(),
                                 self.y.ix[tmp_train_ix].copy())

                # predict test data in kfold validation
                tmp_prob = self.clf.predict_proba(self.x.ix[tmp_test_ix])
                onco_prob.ix[tmp_test_ix] += tmp_prob[:, self.onco_num]
                tsg_prob.ix[tmp_test_ix] += tmp_prob[:, self.tsg_num]

            self.num_pred += 1

        # convert number of trees to fraction of trees
        onco_prob /= self.num_pred
        tsg_prob /= self.num_pred
        other_prob = 1 - (onco_prob + tsg_prob)

        # return prediction.astype(int), prob
        return onco_prob, tsg_prob, other_prob

    def _init_metrics(self):
        """Initialize classification diagnostric metrics."""
        self.feature_importance = []
        self.confusion_matrix = np.zeros((self.num_classes, self.num_classes))
        self.num_pred = 0  # no predictions have been made yet

        # oncogene metrics
        num_points = 100  # number of points to calculate metrics
        self.onco_f1_score = np.zeros(self.total_iter)
        self.onco_precision = np.zeros(self.total_iter)
        self.onco_recall = np.zeros(self.total_iter)
        self.onco_gene_count = np.zeros(self.total_iter)
        self.onco_tpr_array = np.zeros((self.total_iter, num_points))
        self.onco_fpr_array = np.linspace(0, 1, num_points)
        self.onco_precision_array = np.zeros((self.total_iter, num_points))
        self.onco_recall_array = np.linspace(0, 1, num_points)
        self.onco_threshold_array = np.zeros((self.total_iter, num_points))

        # tsg metrics
        self.tsg_f1_score = np.zeros(self.total_iter)
        self.tsg_precision = np.zeros(self.total_iter)
        self.tsg_recall = np.zeros(self.total_iter)
        self.tsg_gene_count = np.zeros(self.total_iter)
        self.tsg_tpr_array = np.zeros((self.total_iter, num_points))
        self.tsg_fpr_array = np.linspace(0, 1, num_points)
        self.tsg_precision_array = np.zeros((self.total_iter, num_points))
        self.tsg_recall_array = np.linspace(0, 1, num_points)

        # driver metrics
        self.driver_precision = np.zeros(self.total_iter)
        self.driver_recall = np.zeros(self.total_iter)
        self.driver_tpr_array = np.zeros((self.total_iter, num_points))
        self.driver_fpr_array = np.linspace(0, 1, num_points)
        self.driver_precision_array = np.zeros((self.total_iter, num_points))
        self.driver_recall_array = np.linspace(0, 1, num_points)
        self.driver_threshold_array = np.zeros((self.total_iter, num_points))
        self.cancer_gene_count = np.zeros(self.total_iter)

        # overall metrics
        self.f1_score = np.zeros(self.total_iter)
        self.precision = np.zeros(self.total_iter)
        self.recall = np.zeros(self.total_iter)

    def _update_metrics(self, y_true, y_pred,
                        onco_prob, tsg_prob):
        # record which genes were predicted what
        self.driver_gene_pred = pd.Series(y_pred, self.y.index)
        self.driver_gene_score = pd.Series(onco_prob+tsg_prob, self.y.index)

        # evaluate performance
        prec, recall, fscore, support = metrics.precision_recall_fscore_support(y_true, y_pred,
                                                                                average='macro')
        cancer_gene_pred = ((onco_prob + tsg_prob)>.5).astype(int)
        self.cancer_gene_count[self.num_pred] = np.sum(cancer_gene_pred)
        self.precision[self.num_pred] = prec
        self.recall[self.num_pred] = recall
        self.f1_score[self.num_pred] = fscore

        # compute Precision-Recall curve metrics
        driver_prob = onco_prob + tsg_prob
        driver_true = (y_true > 0).astype(int)
        p, r, thresh = metrics.precision_recall_curve(driver_true, driver_prob)
        p, r, thresh = p[::-1], r[::-1], thresh[::-1]  # reverse order of results
        thresh = np.insert(thresh, 0, 1.0)
        self.driver_precision_array[self.num_pred, :] = interp(self.driver_recall_array, r, p)
        self.driver_threshold_array[self.num_pred, :] = interp(self.driver_recall_array, r, thresh)

        # calculate prediction summary statistics
        prec, recall, fscore, support = metrics.precision_recall_fscore_support(driver_true, cancer_gene_pred)
        self.driver_precision[self.num_pred] = prec[1]
        self.driver_recall[self.num_pred] = recall[1]

        # save driver metrics
        fpr, tpr, thresholds = metrics.roc_curve(driver_true, driver_prob)
        self.driver_tpr_array[self.num_pred, :] = interp(self.driver_fpr_array, fpr, tpr)

    def _update_onco_metrics(self, y_true, y_pred, prob):
        self.onco_gene_pred = pd.Series(y_pred, self.y.index)
        self.onco_gene_score = pd.Series(prob, self.y.index)

        # compute metrics for classification
        self.onco_gene_count[self.num_pred] = sum(y_pred)
        prec, recall, fscore, support = metrics.precision_recall_fscore_support(y_true, y_pred)
        self.onco_precision[self.num_pred] = prec[self.onco_num]
        self.onco_recall[self.num_pred] = recall[self.onco_num]
        self.onco_f1_score[self.num_pred] = fscore[self.onco_num]
        self.logger.debug('Onco Iter %d: Precission=%s, Recall=%s, f1_score=%s' % (
                          self.num_pred + 1, str(prec), str(recall), str(fscore)))

        # compute ROC curve metrics
        fpr, tpr, thresholds = metrics.roc_curve(y_true, prob)
        self.onco_tpr_array[self.num_pred, :] = interp(self.onco_fpr_array, fpr, tpr)
        #self.onco_mean_tpr[0] = 0.0

        # compute Precision-Recall curve metrics
        p, r, thresh = metrics.precision_recall_curve(y_true, prob)
        p, r, thresh = p[::-1], r[::-1], thresh[::-1]  # reverse order of results
        thresh = np.insert(thresh, 0, 1.0)
        self.onco_precision_array[self.num_pred, :] = interp(self.onco_recall_array, r, p)
        self.onco_threshold_array[self.num_pred, :] = interp(self.onco_recall_array, r, thresh)

    def _update_tsg_metrics(self, y_true, y_pred, prob):
        self.tsg_gene_pred = pd.Series(y_pred, self.y.index)
        self.tsg_gene_score = pd.Series(prob, self.y.index)

        # compute metrics for classification
        self.tsg_gene_count[self.num_pred] = sum(y_pred)
        prec, recall, fscore, support = metrics.precision_recall_fscore_support(y_true, y_pred)
        tsg_col = 1  # column for metrics relate to tsg
        self.tsg_precision[self.num_pred] = prec[tsg_col]
        self.tsg_recall[self.num_pred] = recall[tsg_col]
        self.tsg_f1_score[self.num_pred] = fscore[tsg_col]
        self.logger.debug('Tsg Iter %d: Precission=%s, Recall=%s, f1_score=%s' % (
                          self.num_pred + 1, str(prec), str(recall), str(fscore)))

        # compute ROC curve metrics
        fpr, tpr, thresholds = metrics.roc_curve(y_true, prob)
        self.tsg_tpr_array[self.num_pred, :] = interp(self.tsg_fpr_array, fpr, tpr)
        #self.tsg_tpr_array[0] = 0.0

        # compute Precision-Recall curve metrics
        p, r, thresh = metrics.precision_recall_curve(y_true, prob)
        p, r, thresh = p[::-1], r[::-1], thresh[::-1]  # reverse order of results
        self.tsg_precision_array[self.num_pred, :] = interp(self.tsg_recall_array, r, p)

    def _on_finish(self):
        """Record the mean ROC and PR AUC after all of the cross-validation is finished
        in kfold_validation method.
        """
        # ROC curve metrics
        self.onco_mean_roc_auc = float(metrics.auc(self.onco_fpr_array,
                                                   np.mean(self.onco_tpr_array, axis=0)))
        self.tsg_mean_roc_auc = float(metrics.auc(self.tsg_fpr_array,
                                                  np.mean(self.tsg_tpr_array, axis=0)))
        self.driver_mean_roc_auc = float(metrics.auc(self.driver_fpr_array,
                                                     np.mean(self.driver_tpr_array, axis=0)))

        # Precision-Recall curve metrics
        self.onco_mean_pr_auc = float(metrics.auc(self.onco_recall_array,
                                                  np.mean(self.onco_precision_array, axis=0)))
        self.tsg_mean_pr_auc = float(metrics.auc(self.tsg_recall_array,
                                                 np.mean(self.tsg_precision_array, axis=0)))
        self.driver_mean_pr_auc = float(metrics.auc(self.driver_recall_array,
                                                    np.mean(self.driver_precision_array, axis=0)))

        # log info on classifier predictions
        self.logger.debug('TSG: Precision=%s, Recall=%s, Fscore=%s' % (
                          np.mean(self.tsg_precision), np.mean(self.tsg_recall), np.mean(self.tsg_f1_score)))
        self.logger.debug('Oncogene: Precision=%s, Recall=%s, Fscore=%s' % (
                          np.mean(self.onco_precision), np.mean(self.onco_recall), np.mean(self.onco_f1_score)))
        self.logger.debug('Driver: Precision=%s, Recall=%s' % (
                          np.mean(self.driver_precision), np.mean(self.driver_recall)))


    def get_onco_roc_metrics(self):
        """Simple get method for oncogene ROC metrics."""
        return self.onco_tpr_array, self.onco_fpr_array, self.onco_mean_roc_auc

    def get_tsg_roc_metrics(self):
        """Simple get method for tumor supressor ROC metrics."""
        return self.tsg_tpr_array, self.tsg_fpr_array, self.tsg_mean_roc_auc

    def get_onco_pr_metrics(self):
        """Simple get method for oncogene Precision-Recall metrics."""
        return self.onco_precision_array, self.onco_recall_array, self.onco_mean_pr_auc

    def get_tsg_pr_metrics(self):
        """Simple get method for tumor supressor Precision-Recall metrics"""
        return self.tsg_precision_array, self.tsg_recall_array, self.tsg_mean_pr_auc

    def get_driver_pr_metrics(self):
        """Simple get method for driver gene Precision-Recall metrics"""
        return self.driver_precision_array, self.driver_recall_array, self.driver_mean_pr_auc

    def get_driver_roc_metrics(self):
        """Simple get method for driver gene ROC metrics."""
        return self.driver_tpr_array, self.driver_fpr_array, self.driver_mean_roc_auc

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

        # significantly mutated genes
        self.smg_num = 1

    def set_min_count(self, ct):
        """Set a minimum count threshold."""
        if ct >= 0:
            self.min_count = ct
