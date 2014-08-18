from ..src.classify.python.vogelstein_classifier import VogelsteinClassifier


def test_oncogene_counts():
    """Test if VogelsteinClassifier classifies oncogenes by cts correctly."""
    vclf = VogelsteinClassifier()  # get 20/20 rule classifier

    # example data
    recur1, del1, total1 = 4, 3, 7
    recur2, del2, total2 = 3, 4, 7
    test_list = [[5, 1, 8], [12, 0, 13]]

    # case 1 -- should recognize as oncogene
    gene_class = vclf.predict_by_cts(recur1, del1, total1)
    assert gene_class == vclf.onco_label

    # case 2 -- low counts
    vclf.set_min_count(8)
    gene_class = vclf.predict_by_cts(recur1, del1, total1)
    assert gene_class == vclf.other_label

    # case 3 -- higher deleterious counts despite >20%
    vclf.set_min_count(7)
    gene_class = vclf.predict_by_cts(recur2, del2, total2)
    assert gene_class == vclf.tsg_label

    # case 4 -- input a list of counts
    gene_class_list = vclf.predict_list(test_list)
    assert gene_class_list == [vclf.onco_label] * 2


def test_oncogene_pct():
    """Test if VogelsteinClassifier classifies oncogenes by pct correctly."""
    vclf = VogelsteinClassifier()

    # example data
    recur1, del1, total1 = .3, .1, 10
    recur2, del2, total2 = .3, .4, 10

    # case 1 -- should recognize as oncogene
    gene_class = vclf.predict_by_pct(recur1, del1, total1)
    assert gene_class == vclf.onco_label

    # case 2 -- low counts
    vclf.set_min_count(11)
    gene_class = vclf.predict_by_pct(recur1, del1, total1)
    assert gene_class == vclf.other_label

    # case 3 -- higher deleterious counts despite >20%
    vclf.set_min_count(7)
    gene_class = vclf.predict_by_pct(recur2, del2, total2)
    assert gene_class == vclf.tsg_label


def test_tsg_counts():
    """Test if VogelsteinClassifier classifies tsg by cts correctly."""
    vclf = VogelsteinClassifier()  # get 20/20 rule classifier

    # example data
    recur1, del1, total1 = 3, 4, 7
    recur2, del2, total2 = 4, 3, 7
    test_list = [[1, 5, 8], [0, 12, 13]]

    # case 1 -- should recognize as tsg
    gene_class = vclf.predict_by_cts(recur1, del1, total1)
    assert gene_class == vclf.tsg_label

    # case 2 -- low counts
    vclf.set_min_count(8)
    gene_class = vclf.predict_by_cts(recur1, del1, total1)
    assert gene_class == vclf.other_label

    # case 3 -- higher recurrent counts despite >20%
    vclf.set_min_count(7)
    gene_class = vclf.predict_by_cts(recur2, del2, total2)
    assert gene_class == vclf.onco_label

    # case 4 -- input a list of counts
    gene_class_list = vclf.predict_list(test_list)
    assert gene_class_list == [vclf.tsg_label] * 2


def test_tsg_pct():
    """Test if VogelsteinClassifier classifies tsg by pct correctly."""
    vclf = VogelsteinClassifier()

    # example data
    recur1, del1, total1 = .1, .3, 10
    recur2, del2, total2 = .4, .3, 10

    # case 1 -- should recognize as tsg
    gene_class = vclf.predict_by_pct(recur1, del1, total1)
    assert gene_class == vclf.tsg_label

    # case 2 -- low counts
    vclf.set_min_count(11)
    gene_class = vclf.predict_by_pct(recur1, del1, total1)
    assert gene_class == vclf.other_label

    # case 3 -- higher recurrent counts despite >20%
    vclf.set_min_count(7)
    gene_class = vclf.predict_by_pct(recur2, del2, total2)
    assert gene_class == vclf.onco_label

