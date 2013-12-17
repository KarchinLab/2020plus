import rpy2.robjects as ro
import pandas.rpy.common as com
import pandas as pd
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
                                   replace=TRUE,
                                   ntree=ntree,
                                   classwt=1/sampSize,
                                   sampsize=sampSize)
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
        self.sample_size = ro.IntVector(sampsize)

    def fit(self, xtrain, ytrain):
        label_counts = ytrain.value_counts()
        sampsize = [label_counts[self.other_num],
                    label_counts[self.onco_num],
                    label_counts[self.tsg_num]]
        self.set_sample_size(sampsize)
        xtrain['true_class'] = ytrain
        r_xtrain = com.convert_to_r_dataframe(xtrain)
        self.rf = self.rf_fit(r_xtrain, self.ntree, self.sample_size)

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

    def predict(self, xtest):
        r_xtest = com.convert_to_r_dataframe(xtest)
        pred = self.rf_pred(self.rf, r_xtest)
        py_pred = com.convert_robj(pred)
        genes, pred_class = zip(*py_pred.items())
        tmp_df = pd.DataFrame({'pred_class': pred_class},
                              index=genes)
        tmp_df = tmp_df.reindex(xtest.index)
        tmp_df -= 1  # for some reason the class numbers start at 1
        return tmp_df

    def predict_proba(self, xtest):
        r_xtest = com.convert_to_r_dataframe(xtest)
        pred_prob = self.rf_pred_prob(self.rf, r_xtest)
        py_pred_prob = com.convert_robj(pred_prob)
        return py_pred_prob.values


class RRandomForest(GenericClassifier):

    def __init__(self, df,
                 total_iter=5,
                 weight=False,
                 min_ct=5):
        self.logger = logging.getLogger(__name__)
        onco_flag, tsg_flag = True, True  # classes to actually classify
        super(RRandomForest, self).__init__(total_iter,
                                            classify_oncogene=onco_flag,
                                            classify_tsg=tsg_flag)  # call base constructor
        # self.set_min_count(min_ct)
        self.is_weighted_sample = weight

        df = df.fillna(df.mean())
        self.x, self.y = self._randomize(df)

        self.clf = MyClassifier()
        self.clf.set_classes(onco_flag, tsg_flag)
