from sklearn.naive_bayes import MultinomialNB
from generic_classifier import GenericClassifier


class MultinomialNaiveBayes(GenericClassifier):

    def __init__(self, df, min_ct=10):
        super(MultinomialNaiveBayes, self).__init__()  # call base constructor
        self.set_min_count(min_ct)

        # process data
        df = self._filter_rows(df)  # filter out low count rows
        row_sums = df.sum(axis=1).astype(float)
        df = df.div(row_sums, axis=0)  # normalize each row
        #self.x = self._random_sort(df)  # randomly sort data
        #self.y = self.x.index.to_series().apply(self._label_gene)  # get gene labels
        self.x, self.y = self._randomize(df)

        # setup classifier
        self.clf = MultinomialNB(alpha=1,  # laplacian smooth, i.e. pseudocounts
                                 fit_prior=True)  # use data for prior class probs

    def _update_metrics(self):
        self.num_pred += 1

    def _on_finish(self):
        pass
