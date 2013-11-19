from __future__ import division

class VogelsteinClassifier(object):
    """Oncogene and TSG classifier based on the 20/20 rule.

    Essentially the 20/20 rule states that oncogenes have at least
    20% recurrent missense mutations. While tumor suppressor genes
    have atleast 20% deleterius mutations. This is a simple rule-
    based classifier. To reduce errors for genes with low mutation
    counts, vogelstein et al. manually curated genes with between
    7 and 18 mutations. This class can not reproduce their manual
    curation but can give an estimate on the accuracy of a naive
    implementation of the 20/20 rule. The threshold of 20% is also
    changeable.

    Information on the 20/20 rule from vogelstein's science paper:
        http://www.sciencemag.org/content/339/6127/1546.full
    """

    def __init__(self, onco_threshold=.2, tsg_threshold=.2,
                 min_count=10):
        # check valid percentage
        if not 0 < onco_threshold < 1:
            raise ValueError("Oncogene threshold is invalid")
        if not 0 < tsg_threshold < 1:
            raise ValueError("TSG threshold is invalid")

        # assign percentage thresholds
        self.onco_threshold = onco_threshold
        self.tsg_threshold = tsg_threshold

        # set min count to classify gene
        self.min_count = min_count

        # labels to classify genes as
        self.onco_label = "oncogene"
        self.tsg_label = "tsg"
        self.other_label = "other"

    def predict_list(self, input_list, kind='count'):
        """Predict a list of inputs as either oncogene/tsg/other."""
        gene_class_list = []
        if kind == 'count':
            for recur_ct, del_ct, total_ct in input_list:
                tmp_gene_class = self.predict_by_cts(recur_ct,
                                                     del_ct,
                                                     total_ct)
                gene_class_list.append(tmp_gene_class)
        else:
            for recur_pct, del_pct, total_cts in input_list:
                tmp_gene_class = self.predict_by_pct(recur_pct,
                                                     del_pct,
                                                     total_cts)
                gene_class_list.append(tmp_gene_class)
        return gene_class_list

    def predict_by_cts(self, recurrent, deleterious, total):
        """Predicts oncogene/tsg/other by gene mutation counts."""
        if total < self.min_count:
            # too few mutations case
            return self.other_label

        # sufficient number of counts
        recur_perc = recurrent / total
        del_perc = deleterious / total
        gene_class = self.predict_by_pct(recur_perc,
                                         del_perc,
                                         total)
        return gene_class

    def predict_by_pct(self, recur_pct, del_pct, total=10):
        """The actual 20/20 rule logic to classify genes."""
        if total < self.min_count:
            # too few mutations case
            return self.other_label
        elif recur_pct >= self.onco_threshold:
            # high number of recurrent missense case
            if recur_pct >= del_pct:
                return self.onco_label
            else:
                return self.tsg_label
        elif del_pct >= self.tsg_threshold:
            # high number of deleterious mutations case
            return self.tsg_label
        else:
            # doesn't classify as oncogene or tsg
            return self.other_label

    def set_onco_threshold(self, threshold):
        """Setter for percentage threshold for recurrent missense mutations
        to call it an oncogene."""
        if 0 < threshold < 1:
            self.onco_threshold = threshold

    def set_tsg_threshold(self, threshold):
        """Setter for percentage threshold for deleterious mutations to
        call it a tsg."""
        if 0 < threshold < 1:
            self.tsg_threshold = threshold

    def set_min_count(self, count):
        """Setter for minimum count that can be classified for either a
        oncogene or tsg."""
        if count > 0:
            self.min_count = count
