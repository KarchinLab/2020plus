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

    def __init__(self,
                 onco_threshold=.2,
                 tsg_threshold=.2,
                 kind='vogelstein',
                 min_count=0,
                 tsg_min=7,
                 onco_min=10,
                 db_size=404863):  # db size is as reported from Cancer Genome Landscapes paper
        # check valid percentage
        if not 0 < onco_threshold < 1:
            raise ValueError("Oncogene threshold is invalid")
        if not 0 < tsg_threshold < 1:
            raise ValueError("TSG threshold is invalid")

        self.kind = kind  # either 'vogelstein' or 'min'

        # set parameters as reported in Cancer genome landscapes paper
        self.db_size = db_size
        self.db_tsg_min = tsg_min
        self.db_onco_min = onco_min

        # assign percentage thresholds
        self.onco_threshold = onco_threshold
        self.tsg_threshold = tsg_threshold

        # set min count to classify gene
        self.min_count = min_count
        self.tsg_min = tsg_min
        self.onco_min = onco_min

        # labels to classify genes as
        self.onco_label = "oncogene"
        self.tsg_label = "tsg"
        self.other_label = "other"

    def predict_list(self, input_list,
                     kind='count', scale_type='linear'):
        """Predict a list of inputs as either oncogene/tsg/other."""
        # scale count thresholds
        all_cts = sum([x[-1] for x in input_list])
        if scale_type:
            self.tsg_min = self.db_tsg_min * float(all_cts)/self.db_size
            self.onco_min = self.db_onco_min * float(all_cts)/self.db_size
        else:
            self.tsg_min = self.db_tsg_min
            self.onco_min = self.db_onco_min

        # perform prediction
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
        recur_perc = recurrent / float(total)
        del_perc = deleterious / float(total)
        gene_class = self.predict_by_pct(recur_perc,
                                         del_perc,
                                         total)
        return gene_class

    def predict_by_pct(self, recur_pct, del_pct, total):
        """The actual 20/20 rule logic to classify genes."""
        # calc counts
        recur_ct = recur_pct * total
        del_ct = del_pct * total

        # 20/20 rule logic
        if self.kind == 'vogelstein':
            if recur_pct >= self.onco_threshold and recur_ct >= self.onco_min:
                if del_pct <= .05:
                    return self.onco_label
                elif del_ct >= self.tsg_min:
                    return self.tsg_label
                else:
                    return self.other_label
            elif del_pct >= self.tsg_threshold and del_ct >= self.tsg_min:
                return self.tsg_label
            else:
                return self.other_label
        elif self.kind == 'min':
            if total < self.min_count:
                # too few mutations case
                return self.other_label
            # if recur_pct >= self.onco_threshold and (total*recur_pct) >= self.min_count:
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
