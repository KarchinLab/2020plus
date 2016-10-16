"""Random Forest classifier using R's randomForest library.
RPy2 is used to interface with R."""
import readline  # hopefully fix libreadline error
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import pandas.rpy.common as com
import pandas as pd
from src.classify.python.generic_classifier import GenericClassifier
import src.features.python.feature_utils as futils
import os
import logging

class MyClassifier(object):
    """
    The MyClassifier class is intended to only be used by RRandomForest.
    Essentially MyClassifier is a wrapper around R's randomForest libary,
    to make it compatible with the methods used by scikit learn (predict,
    fit, predic_proba). The advantage of R's version of random forest is
    that it allows the sampling rate to be specified.
    """

    def __init__(self,
                 ntrees=200,
                 other_sample_ratio=1.,
                 driver_sample=.7):
        self.ntrees = ntrees
        # self.other_sample_rate = other_sample
        self.other_sample_ratio = other_sample_ratio
        self.driver_sample_rate = driver_sample

        # Code for R's random forest using rpy2
        ro.r("suppressPackageStartupMessages(library(randomForest))")  # load randomForest library

        # set up lists holding trained models and what genes
        # were trained on
        ro.r("trained.models <- list()")
        ro.r("tmp.trained.models <- list()")

        # R function for fitting a random forest
        ro.r('''rf_fit <- function(df, ntree, sampSize){
                df$true_class <- as.factor(df$true_class)
                rf_clf <<- randomForest(true_class~.,
                                        data=df,
                                        replace=TRUE,
                                        ntree=ntree,
                                        classwt=1/sampSize,
                                        sampsize=sampSize,
                                        importance=TRUE)
                return(rf_clf)
             }''')
        self.rf_fit = ro.r['rf_fit']

        # R function for getting feature importance
        ro.r('''rf_imp <- function(rf){
                myimp <- importance(rf)
                return(myimp[,5])
             }''')
        self.rf_imp = ro.r['rf_imp']

        # R function for predicting class probability
        ro.r('''rf_pred_prob <- function(rf, xtest){
                prob <- predict(rf, xtest, type="prob")
                return(prob)
             }''')
        self.rf_pred_prob = ro.r['rf_pred_prob']

        # R function for predicting class
        ro.r('''rf_pred <- function(rf, xtest){
                prob <- predict(rf, xtest)
                return(prob)
             }''')
        self.rf_pred = ro.r['rf_pred']

    def set_sample_size(self, sampsize):
        # sampsize[0] *= self.other_sample_rate
        if self.is_onco_pred and self.is_tsg_pred:
            sampsize[1] *= self.driver_sample_rate
            sampsize[2] *= self.driver_sample_rate
            tmp_total_driver = sampsize[1] + sampsize[2]
        else:
            sampsize[1] *= self.driver_sample_rate
            tmp_total_driver = sampsize[1]
        sampsize[0] = int(self.other_sample_ratio * tmp_total_driver)
        self.sample_size = ro.IntVector(sampsize)

    def fit(self, xtrain, ytrain):
        """The fit method trains R's random forest classifier.

        NOTE: the method name ("fit") and method signature were choosen
        to be consistent with scikit learn's fit method.

        Parameters
        ----------
        xtrain : pd.DataFrame
            features for training set
        ytrain : pd.DataFrame
            true class labels (as integers) for training set
        """
        label_counts = ytrain.value_counts()
        if self.is_onco_pred and self.is_tsg_pred:
            sampsize = [label_counts[self.other_num],
                        label_counts[self.onco_num],
                        label_counts[self.tsg_num]]
        elif self.is_onco_pred:
            sampsize = [label_counts[self.other_num],
                        label_counts[self.onco_num]]
        elif self.is_tsg_pred:
            sampsize = [label_counts[self.other_num],
                        label_counts[self.tsg_num]]

        self.set_sample_size(sampsize)
        ytrain.index = xtrain.index  # ensure indexes match
        xtrain['true_class'] = ytrain
        r_xtrain = com.convert_to_r_dataframe(xtrain)
        #ro.globalenv['trainData'] = r_xtrain
        #r_xtrain = pandas2ri.py2ri(xtrain)
        self.rf = self.rf_fit(r_xtrain, self.ntrees, self.sample_size)
        r_imp = self.rf_imp(self.rf)  # importance dataframe in R
        self.feature_importances_ = com.convert_robj(r_imp)
        #self.feature_importances_ = pandas2ri.ri2py(r_imp)

    def save(self, path):
        """Save random forest model as a Rdata file."""
        ro.r('''save(rf_clf, rf_pred_prob, rf_pred,
                     rf_imp, rf_fit, file="{0}")'''.format(path))

    def save_cv(self, path):
        """Save random forest model as a Rdata file."""
        ro.r('''save(rf_pred_prob, rf_pred, rf_imp, rf_fit, cvFoldDf,
                     trained.models, file="{0}")'''.format(path))

    def load(self, path):
        set_wd_str = 'setwd("{0}")'.format(os.getcwd())
        ro.r(set_wd_str)
        ro.r('load("{0}")'.format(path))
        self.rf = ro.r["rf_clf"]

    def load_cv(self, path):
        set_wd_str = 'setwd("{0}")'.format(os.getcwd())
        ro.r(set_wd_str)
        ro.r('load("{0}")'.format(path))
        self.rf_cv = ro.r["trained.models"]
        self.cv_folds = com.convert_robj(ro.r["cvFoldDf"])

    def set_model(self, num_iter, num_fold):
        """Set which random forest model is currently active.

        Whith cross-validation there may be many models saved.
        There will be one for each iteration and fold for cross-validation.

        Note: due to rpy2 indexing, indexes should start from 1

        Parameters
        ----------
        num_iter : int
            Iteration index amongst repeated cross-validations
        num_fold : int
            The particular fold for each cross-validation
        """
        self.rf = self.rf_cv.rx2(num_iter).rx2(num_fold)

    def append_cv_result(self):
        """Append result for cross-validation."""
        ro.r("trained.models <- append(trained.models, list(tmp.trained.models)); tmp.trained.models <- list()")

    def append_fold_result(self):
        """Append result for each cross-validation fold."""
        ro.r("tmp.trained.models <- append(tmp.trained.models, list(rf_clf))")

    def set_cv_fold(self, df):
        """Send which genes are valid test sets for each CV fold."""
        r_df = com.convert_to_r_dataframe(df)
        ro.globalenv['cvFoldDf'] = r_df

    def set_classes(self, oncogene, tsg):
        """Sets the integers used to represent classes in classification."""
        if not oncogene and not tsg:
            raise ValueError('Classification needs at least two classes')
        self.is_onco_pred = oncogene
        self.is_tsg_pred = tsg
        if oncogene and tsg:
            self.other_num = 0
            self.onco_num = 1
            self.tsg_num = 2
            self.num_classes = 3
        else:
            self.other_num = 0
            self.num_classes = 2
            self.onco_num = 1 if oncogene else 0
            self.tsg_num = 1 if tsg else 0

    def set_seed(self, seed):
        if seed is not None:
            ro.r('set.seed({0})'.format(seed))

    def predict(self, xtest):
        """Predicts class via majority vote.

        Parameters
        ----------
        xtest : pd.DataFrame
            features for test set
        """
        r_xtest = com.convert_to_r_dataframe(xtest)
        #r_xtest = pandas2ri.py2ri(xtest)
        pred = self.rf_pred(self.rf, r_xtest)
        py_pred = com.convert_robj(pred)
        #py_pred = pandas2ri.ri2py(pred)
        genes, pred_class = zip(*py_pred.items())
        tmp_df = pd.DataFrame({'pred_class': pred_class},
                              index=genes)
        tmp_df = tmp_df.reindex(xtest.index)
        tmp_df -= 1  # for some reason the class numbers start at 1
        return tmp_df

    def predict_proba(self, xtest):
        """Predicts the probability for each class.

        Parameters
        ----------
        xtest : pd.DataFrame
            features for test set
        """
        r_xtest = com.convert_to_r_dataframe(xtest)
        #r_xtest = pandas2ri.ri2py(xtest)
        pred_prob = self.rf_pred_prob(self.rf, r_xtest)
        py_pred_prob = com.convert_robj(pred_prob)
        #py_pred_prob = pandas2ri.ri2py(pred_prob)
        return py_pred_prob.values


class RRandomForest(GenericClassifier):

    def __init__(self, df,
                 total_iter=5,
                 weight=False,
                 ntrees=200,
                 other_sample_ratio=3.,
                 driver_sample=.7,
                 seed=None):
        self.logger = logging.getLogger(__name__)
        onco_flag, tsg_flag = True, True  # classes to actually classify
        super(RRandomForest, self).__init__(total_iter,
                                            classify_oncogene=onco_flag,
                                            classify_tsg=tsg_flag,
                                            rseed=seed)  # call base constructor
        self.is_weighted_sample = weight

        if 'total' in df.columns:
            # hack to get rid of total mutation count column
            df = df.drop('total', axis=1)
        df = df.fillna(df.mean())

        # randomization is mostly done in prediciton methods
        self.x, self.y = futils.randomize(df, self.prng)

        # use the MyClassifier wrapper class around R
        self.clf = MyClassifier(ntrees=ntrees,
                                driver_sample=driver_sample,
                                other_sample_ratio=other_sample_ratio)
        self.clf.set_classes(onco_flag, tsg_flag)
        self.clf.set_seed(seed)

    def _update_metrics(self, y_true, y_pred, onco_prob, tsg_prob):
        super(RRandomForest, self)._update_metrics(y_true, y_pred, onco_prob, tsg_prob)

        # evaluate feature importance for random forest
        self.feature_importance.append(self.clf.feature_importances_.tolist())

    def _on_finish(self):
        super(RRandomForest, self)._on_finish()
        self.feature_importance = pd.DataFrame(self.feature_importance,
                                               columns=self.x.columns)
        self.mean_importance = self.feature_importance.mean()
        self.std_importance = self.feature_importance.std()
