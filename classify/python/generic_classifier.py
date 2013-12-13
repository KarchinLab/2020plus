import utils.python.util as _utils
import numpy as np
from numpy import interp
from sklearn import cross_validation
import sklearn.metrics as metrics
import pandas as pd


class GenericClassifier(object):

    def __init__(self,
                 total_iterations=5):
        self.min_count = 0  # min mutations for a gene

        # set integer codes for classes
        self.set_classes(oncogene=True, tsg=True)

        self.set_total_iter(total_iterations)  # initialize total number of simulations
        self._init_metrics()  # initialize metrics

        # genes categorized as oncogenes/tsg by vogelstein's
        # science paper
        self.vogelsteins_oncogenes = _utils.oncogene_set
        self.vogelsteins_tsg = _utils.tsg_set

        # gene's categorized by the significantly mutate gene
        # list from the pan-cancer nature paper
        self.smg_list = _utils.smg_list

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
        self.onco_tpr_array = np.zeros((self.total_iter, num_points))
        self.onco_fpr_array = np.linspace(0, 1, num_points)
        self.onco_precision_array = np.zeros((self.total_iter, num_points))
        self.onco_recall_array = np.linspace(0, 1, num_points)
        self.onco_threshold_array = np.zeros((self.total_iter, num_points))

        # tsg metrics
        self.tsg_f1_score = np.zeros(self.total_iter)
        self.tsg_precision = np.zeros(self.total_iter)
        self.tsg_recall = np.zeros(self.total_iter)
        self.tsg_tpr_array = np.zeros((self.total_iter, num_points))
        self.tsg_fpr_array = np.linspace(0, 1, num_points)
        self.tsg_precision_array = np.zeros((self.total_iter, num_points))
        self.tsg_recall_array = np.linspace(0, 1, num_points)

        # overall metrics
        self.f1_score = np.zeros(self.total_iter)
        self.precision = np.zeros(self.total_iter)
        self.recall = np.zeros(self.total_iter)

    def set_total_iter(self, myiterations):
        self.total_iter = myiterations

    def _filter_rows(self, df):
        """Filter out rows with counts less than the minimum."""
        row_sums = df.T.sum()
        filtered_df = df[row_sums > self.min_count]
        return filtered_df

    def _label_gene(self, gene,
                    kind='vogelstein'):
        """Label a gene according to Vogelstein's list of oncogenes
        and tsg."""
        if kind == 'vogelstein':
            if gene in self.vogelsteins_oncogenes:
                return self.onco_num
            elif gene in self.vogelsteins_tsg:
                return self.tsg_num
            else:
                return self.other_num
        elif kind == 'smg':
            if gene in self.smg_list:
                return self.smg_num
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

    def _update_metrics(self, y_true, y_pred):
        prec, recall, fscore, support = metrics.precision_recall_fscore_support(y_true, y_pred,
                                                                                average='macro')
        self.precision[self.num_pred] = prec
        self.recall[self.num_pred] = recall
        self.f1_score[self.num_pred] = fscore

    def _update_onco_metrics(self, y_true, y_pred, prob):

        # estabilish confusion matrix
        #tmp_confusion_matrix = metrics.confusion_matrix(y_true, y_pred)
        #self.confusion_matrix += tmp_confusion_matrix

        # compute metrics for classification
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

        # estabilish confusion matrix
        #tmp_confusion_matrix = metrics.confusion_matrix(y_true, y_pred)
        #self.confusion_matrix += tmp_confusion_matrix

        # compute metrics for classification
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
        #self.confusion_matrix /= self.num_pred
        #self.onco_f1_score /= self.num_pred
        #self.onco_precision /= self.num_pred
        #self.onco_recall /= self.num_pred
        #self.tsg_f1_score /= self.num_pred
        #self.tsg_precision /= self.num_pred
        #self.tsg_recall /= self.num_pred

        # ROC curve metrics
        #self.onco_mean_tpr /= self.num_pred  # divide by number of folds squared
        #self.onco_tpr_array[:, -1] = 1.0  # it always ends at 1
        #self.onco_tpr_array[:, 0] = 0.0
        self.onco_mean_roc_auc = metrics.auc(self.onco_fpr_array,
                                             np.mean(self.onco_tpr_array, axis=0))
        #self.tsg_mean_tpr /= self.num_pred  # divide by number of folds squared
        #self.tsg_tpr_array[:, -1] = 1.0  # it always ends at 1
        #self.tsg_tpr_array[:, 0] = 0.0  # it always begins as 0
        self.tsg_mean_roc_auc = metrics.auc(self.tsg_fpr_array,
                                            np.mean(self.tsg_tpr_array, axis=0))

        # Precision-Recall curve metrics
        #self.onco_mean_precision /= self.num_pred
        self.onco_mean_pr_auc = metrics.auc(self.onco_recall_array,
                                            np.mean(self.onco_precision_array, axis=0))
        #self.onco_mean_threshold /= self.num_pred
        #self.tsg_mean_precision /= self.num_pred
        self.tsg_mean_pr_auc = metrics.auc(self.tsg_recall_array,
                                           np.mean(self.tsg_precision_array, axis=0))

    def kfold_validation(self, k=5):
        # generate indices for kfold cross validation
        #k_fold = cross_validation.KFold(n=len(self.x),  # len of df
                                        #n_folds=k,  # k fold
                                        #indices=True)  # return indices
        self.num_pred = 0  # number of predictions

        for i in range(self.total_iter):
            self.x, self.y = self._randomize(self.x)  # randomize for another round

            # initialize predicted results variables
            num_genes = len(self.y)
            onco_pred = np.zeros(num_genes)
            onco_prob = np.zeros(num_genes)
            tsg_pred = np.zeros(num_genes)
            tsg_prob = np.zeros(num_genes)
            overall_pred = np.zeros(num_genes)

            # set up stratified kfold iterator
            k_fold = cross_validation.StratifiedKFold(self.y,
                                                      n_folds=k,
                                                      indices=True)

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
                    self.clf.fit(self.x.iloc[train_ix],
                                 self.y.iloc[train_ix],
                                 sample_weight=sample_weight)
                else:
                    # do training without sample weights
                    self.clf.fit(self.x.iloc[train_ix], self.y.iloc[train_ix])

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
                                 overall_pred)
            self.num_pred += 1

        self._on_finish()  # update info for kfold cross-validation

    def kfold_prediction_old(self, k=5):
        # generate indices for kfold cross validation
        k_fold = cross_validation.KFold(n=len(self.x),  # len of df
                                        n_folds=k,  # k fold
                                        indices=True)  # return indices
        self.num_pred = 0  # number of predictions
        self.x, self.y = self._randomize(self.x)  # randomize data

        prediction = pd.Series(index=self.y.index)  # predicted class
        onco_prob = pd.Series(index=self.y.index).fillna(0)
        tsg_prob = pd.Series(index=self.y.index).fillna(0)

        for i in range(self.total_iter):
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
        other_prob = 1 - (onco_prob + tsg_prob)

        # return prediction.astype(int), prob
        return onco_prob, tsg_prob, other_prob

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

        for i in range(self.total_iter):
            self.x, self.y = self._randomize(self.x)  # randomize for another round
            # obtain predictions from single round of kfold validation
            for train_ix, test_ix in k_fold:
                # retreive indices from pandas dataframe using row number
                tmp_train_ix = self.x.iloc[train_ix].index
                tmp_test_ix = self.x.iloc[test_ix].index

                if self.is_weighted_sample:
                    num_train = len(train_ix)
                    sample_weight = np.zeros(num_train)
                    onco_ix = np.nonzero(self.y.ix[tmp_train_ix]==self.onco_num)[0]
                    tsg_ix = np.nonzero(self.y.ix[tmp_train_ix]==self.tsg_num)[0]
                    other_ix = np.nonzero(self.y.ix[tmp_train_ix]==self.other_num)[0]
                    sample_weight[onco_ix] = 1. / len(onco_ix)
                    sample_weight[tsg_ix] = 1. / len(tsg_ix)
                    sample_weight[other_ix] = 1. / len(other_ix)

                    # do training with sample weighting
                    self.clf.fit(self.x.ix[tmp_train_ix],
                                 self.y.ix[tmp_train_ix],
                                 sample_weight=sample_weight)
                else:
                    # do training without weighting
                    self.clf.fit(self.x.ix[tmp_train_ix], self.y.ix[tmp_train_ix])

                # predict test data in kfold validation
                tmp_prob = self.clf.predict_proba(self.x.ix[tmp_test_ix])
                onco_prob.ix[tmp_test_ix] += tmp_prob[:, self.onco_num]
                tsg_prob.ix[tmp_test_ix] += tmp_prob[:, self.tsg_num]

            self.num_pred += 1
        onco_prob /= self.num_pred
        tsg_prob /= self.num_pred
        other_prob = 1 - (onco_prob + tsg_prob)

        # return prediction.astype(int), prob
        return onco_prob, tsg_prob, other_prob

    def get_onco_roc_metrics(self):
        return self.onco_tpr_array, self.onco_fpr_array, self.onco_mean_roc_auc

    def get_tsg_roc_metrics(self):
        return self.tsg_tpr_array, self.tsg_fpr_array, self.tsg_mean_roc_auc

    def get_onco_pr_metrics(self):
        return self.onco_precision_array, self.onco_recall_array, self.onco_mean_pr_auc

    def get_tsg_pr_metrics(self):
        return self.tsg_precision_array, self.tsg_recall_array, self.tsg_mean_pr_auc

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
        if ct >= 0:
            self.min_count = ct
