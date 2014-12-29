import numpy as np


def rank_genes(s1, s2,
               mask1=None, mask2=None,
               thresh=None):
    #s1.sort(ascending=False)
    #s2.sort(ascending=False)
    # use score threshold if provided
    if thresh:
        all_ixs = list(set(s1[s1>thresh].index) | set(s2[s2>thresh].index))
    elif mask1 is not None and mask2 is not None:
        all_ixs = list(set(s1[mask1].index) | set(s2[mask2].index))
    else:
        all_ixs = list(set(s1.index) | set(s2.index))

    top_s1 = s1[all_ixs]
    top_s2 = s2[all_ixs]
    top_rank1 = top_s1.fillna(0).rank()
    top_rank2 = top_s2.fillna(0).rank()
    top_rank2 = top_rank2[top_rank1.index]  # match ordering of index
    return top_rank1, top_rank2


def overlap(s1, s2,
            mask1=None, mask2=None,
            thresh=None,
            depth=None):
    # get values of interest, either use a threshold, mask array, or specified
    # depth
    if thresh:
        s1, s2 = s1[s1>thresh], s2[s2>thresh]
    elif mask1 is not None and mask2 is not None:
        s1, s2 = s1[mask1], s2[mask2]
    elif depth is not None:
        s1, s2 = s1[:depth], s2[:depth]

    # genes are supposed to be the index of the series
    s1_genes = set(s1.index)
    s2_genes = set(s2.index)

    # calculate jaccard index
    num_intersect = len(s1_genes & s2_genes)
    num_total = len(s1_genes)
    if num_total:
        # provided series are not empty
        ov_sim = num_intersect / float(num_total)
    else:
        # empty series case
        ov_sim = 0
    return ov_sim


def weighted_overlap(s1, s2,
                     max_depth,
                     step_size,
                     weight_factor):
    # calculate jaccard index at specified intervals
    num_depths = (max_depth) // step_size
    num_depths_total = (len(s2)) // step_size
    ov = np.zeros(num_depths)
    ov_all = np.zeros(num_depths_total)
    for i, depth in enumerate(range(step_size, num_depths_total+1, step_size)):
        if depth <= max_depth:
            #ov_tmp = overlap(s1.iloc[:depth].copy(), s2, depth=max_depth)
            ov_tmp = overlap(s1.iloc[:depth].copy(), s2.iloc[:depth+100].copy())
            ov[i] = ov_tmp
            ov_all[i] = ov_tmp
        else:
            #ov_all[i] = overlap(s1.iloc[:max_depth], s2, depth=depth)
            ov_all[i] = overlap(s1.iloc[:max_depth], s2.iloc[:depth+100].copy())

    # calculate the weighting for jaccard index
    p = weight_factor ** (1./(num_depths-1))
    w = p*np.ones(num_depths_total)
    w[0] = 1
    w = np.cumprod(w)  # calculate geometric weights
    w = w / w.sum()  # normalize weights to 1

    weighted_mean_ov = np.dot(w, ov_all)
    mean_ov = np.mean(ov)

    return ov, mean_ov, weighted_mean_ov


def jaccard_index(s1, s2,
                  mask1=None, mask2=None,
                  thresh=None,
                  depth=None):
    # get values of interest, either use a threshold, mask array, or specified
    # depth
    if thresh:
        s1, s2 = s1[s1>thresh], s2[s2>thresh]
    elif mask1 is not None and mask2 is not None:
        s1, s2 = s1[mask1], s2[mask2]
    elif depth is not None:
        s1, s2 = s1[:depth], s2[:depth]

    # genes are supposed to be the index of the series
    s1_genes = set(s1.index)
    s2_genes = set(s2.index)

    # calculate jaccard index
    num_intersect = len(s1_genes & s2_genes)
    num_union = len(s1_genes | s2_genes)
    if num_union:
        # provided series are not empty
        jaccard_sim = num_intersect / float(num_union)
    else:
        # empty series case
        jaccard_sim = 0
    return jaccard_sim


def weighted_jaccard_index(s1, s2,
                           max_depth,
                           step_size,
                           weight_factor):
    # calculate jaccard index at specified intervals
    num_depths = (max_depth) // step_size
    num_depths_total = (len(s2)) // step_size
    ji = np.zeros(num_depths)
    ji_all = np.zeros(num_depths_total)
    for i, depth in enumerate(range(step_size, num_depths_total+1, step_size)):
        if depth <= max_depth:
            ji_tmp = jaccard_index(s1.iloc[:depth].copy(), s2, depth=max_depth)
            ji[i] = ji_tmp
            ji_all[i] = ji_tmp
        else:
            ji_all[i] = jaccard_index(s1.iloc[:max_depth], s2, depth=depth)

    # calculate the weighting for jaccard index
    p = weight_factor ** (1./(num_depths-1))
    w = p*np.ones(num_depths_total)
    w[0] = 1
    w = np.cumprod(w)  # calculate geometric weights
    w = w / w.sum()  # normalize weights to 1

    weighted_mean_ji = np.dot(w, ji_all)
    mean_ji = np.mean(ji)

    return ji, mean_ji, weighted_mean_ji
