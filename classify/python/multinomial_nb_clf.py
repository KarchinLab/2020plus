from sklearn.naive_bayes import MultinomialNB
from generic_classifier import GenericClassifier


class MultinomialNaiveBayes(GenericClassifier):

    def __init__(self, df, min_ct=10):
        super(MultinomialNaiveBayes, self).__init__()  # call base constructor
        self.set_min_count(min_ct)

        # process data
        df = self._filter_rows(df)  # filter out low count rows
        self.x = self._random_sort(df)  # randomly sort data
        self.y = self.x.index.to_series().apply(self._label_gene)  # get gene labels

        # setup classifier
        self.clf = MultinomialNB(alpha=1,  # laplacian smooth, i.e. pseudocounts
                                 fit_prior=True)  # use data for prior class probs
