import rpy2.robjects as ro
import pandas.rpy.common as com
from generic_classifier import GenericClassifier
import logging


class MyClassifier(object):

    def __init__(self,
                 ntree=200,
                 other_sample=.3,
                 driver_sample=.7):
        self.ntree = ntree
        self.other_sample_rate = other_sample
        self.driver_sample_rate = driver_sample
        ro.r("library(randomForest)")
        ro.r('''rf_fit <- function(df, ntree, sampSize){
                df$true_class <- as.factor(df$true_class)
                rf <- randomForest(true_class~.,
                                   data=df,
                                   replace=True,
                                   ntree=ntree,
                                   classwt=1/sampleSize,
                                   sampsize=sampleSize)
                return(rf)
             }''')
        self.rf_fit = ro.r['rf_fit']
        ro.r('''rf_pred_prob <- function(rf, xtest){
                prob <- predict(rf, xtest, type="prob")
                return(prob)
             }''')
        self.rf_pred_prob = ro.r['rf_pred_prob']
        ro.r('''rf_pred <- function(rf, xtest){
                prob <- predict(rf, xtest)
                return(prob)
             }''')
        self.rf_pred = ro.r['rf_pred']

    def set_sample_size(self, sampsize):
        sampsize[0] *= self.other_sample_rate
        sampsize[1] *= self.driver_sample_rate
        sampsize[2] *= self.driver_sample_rate
        self.sample_size = ro.FloatVector(sampsize)

    def fit(self, xtrain, ytrain):
        label_counts = ytrain.value_counts()
        sampsize = [label_counts['other'],
                    label_counts['oncogene'],
                    label_counts['tsg']]
        self.set_sample_size(sampsize)
        xtrain['true_class'] = ytrain
        r_xtrain = com.convert_to_r_dataframe(xtrain)
        self.rf = self.rf_fit(r_xtrain, self.ntree, self.sample_size)

    def predict(self, xtest):
        r_xtest = com.convert_to_r_dataframe(xtest)
        pred = self.rf_pred(self.rf, r_xtest)
        py_pred = com.convert_robj(pred)
        return py_pred

    def predict_proba(self, xtest):
        r_xtest = com.convert_to_r_dataframe(xtest)
        pred_prob = self.rf_pred_prob(self.rf, r_xtest)
        py_pred_prob = com.convert_robj(pred_prob)
        return py_pred_prob


class RRandomForest(GenericClassifier):

    def __init__(self, df,
                 total_iter=5,
                 weight=True,
                 min_ct=5):
        self.logger = logging.getLogger(__name__)
        super(RRandomForest, self).__init__(total_iter)  # call base constructor
        # self.set_min_count(min_ct)
        self.is_weighted_sample = weight

        df = df.fillna(df.mean())
        self.x, self.y = self._randomize(df)

        self.clf = MyClassifier()
