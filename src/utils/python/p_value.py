import numpy as np
import pandas as pd
import bisect

# genes to be removed from MLFC calc
mlfc_remove_genes = set([
    'PLCG1', 'CRLF2', 'SMARCD1', 'SH2B3', 'STK11', 'MEN1', 'IKBKB', 'AKT1', 'B2M', 'MLH1',
    'USP28', 'TSHR', 'FGFR4', 'GPS2', 'CDC73', 'PIK3CA', 'MAP3K1', 'CACNA1D', 'FGFR3', 'TSC1',
    'ZC3H13', 'CBFB', 'ARHGAP26', 'CDH1', 'JAK2', 'ABCD2', 'NFKBIE', 'CUL1', 'GNAS', 'BRCA1',
    'ERBB3', 'BRCA2', 'DDX3X', 'PIK3R1', 'IDH2', 'IDH1', 'JAK1', 'KLF4', 'STAT5B', 'IWS1',
    'RPL10', 'SMAD2', 'CSF3R', 'KRAS', 'SETD2', 'FGFR2', 'GATA3', 'UBR5', 'ALK', 'KAT8',
    'SOCS1', 'MAX', 'PTPN11', 'ASXL1', 'POLE', 'TCF7L2', 'KIT', 'FOXA1', 'DDX5', 'FAT1',
    'RUNX1', 'CYLD', 'CD79B', 'ARID1B', 'APC', 'DAXX', 'ARID1A', 'CTCF', 'KDM5C', 'IL7R',
    'JAK3', 'CD1D', 'ACVR1B', 'CDKN1B', 'HLA-A', 'RPL5', 'COPS4', 'GGCT', 'TNFRSF14', 'FBXO11',
    'SETBP1', 'TRAF7', 'DNM2', 'CCAR1', 'ARID4B', 'XPO1', 'SPEN', 'KANSL1', 'CHD4', 'WDR47',
    'U2AF1', 'TTR', 'GK2', 'SPOP', 'PPP6C', 'NF2', 'NF1', 'PCBP1', 'AMPH', 'MYD88',
    'GRM3', 'MYH2', 'ATP1A1', 'STAG2', 'ARID2', 'RNF43', 'TUBA3C', 'TRRAP', 'PAX5', 'CUX1',
    'HRAS', 'RAD21', 'HMCN1', 'ING1', 'DNER', 'MGA', 'TP53', 'GNAQ', 'ESR1', 'FAM47B',
    'MPL', 'CBL', 'STK31', 'KRT15', 'PRDM1', 'NFE2L2', 'NSD1', 'MYCN', 'AGTR1', 'CREBBP',
    'ZRSR2', 'PDGFRA', 'GIGYF2', 'SMAD4', 'GATA1', 'ATM', 'EPHA2', 'POT1', 'SMAD3', 'SMO',
    'TBX3', 'CBLB', 'ATR', 'BIRC3', 'ABL1', 'AMER1', 'FAM46C', 'RAF1', 'FLT3', 'NCOR1',
    'CD79A', 'NCOA2', 'TP63', 'STAT3', 'H3F3B', 'MYH9', 'ECT2L', 'CDKN2A(p14)', 'RET', 'MAPK1',
    'KDM6A', 'BRWD3', 'PENK', 'CDKN1A', 'KDR', 'FBXW7', 'TGFBR2', 'FUBP1', 'PDYN', 'RQCD1',
    'GATA2', 'EGFR', 'TET2', 'ZFP36L1', 'ZFP36L2', 'MTOR', 'PRKAR1A', 'BCOR', 'ATRX', 'EP300',
    'TNFAIP3', 'DICER1', 'TBL1XR1', 'KCNJ5', 'MAP2K2', 'COL2A1', 'ALB', 'MAP2K1', 'KEAP1', 'EZH2',
    'CDK4', 'RAC1', 'PBRM1', 'CMTR2', 'BRE', 'RHEB', 'LPAR4', 'AMOT', 'CIC', 'PPP2R1A',
    'ACVR1', 'WT1', 'ZNF318', 'MSH2', 'SF3B1', 'MSH6', 'CTNNB1', 'VHL', 'USP9X', 'SOX9',
    'NOTCH2', 'MAP2K4', 'ELF3', 'SMARCA4', 'H3F3A', 'CEBPA', 'AXIN2', 'AXIN1', 'TWIST1', 'FAS',
    'NRAS', 'RB1', 'CDKN2A', 'KLF6', 'MED12', 'HNF1A', 'ETNK1', 'ATG5', 'CNOT3', 'NRG3',
    'TERT', 'AJUBA', 'NT5C2', 'BRAF', 'KMT2C', 'KMT2B', 'KMT2A', 'DNMT3A', 'SMARCB1', 'KMT2D',
    'PTEN', 'RBM10', 'CARD11', 'GNA11', 'RHOA', 'PTPRB', 'MAP3K13', 'HIST1H3B', 'PHF6', 'TP53BP1',
    'TSC2', 'SUFU', 'DACH1', 'TRIP12', 'PHOX2B', 'NPM1', 'RASA1', 'MYOD1', 'PTCH1', 'ERBB2',
    'CALR', 'SRSF2', 'DKK2', 'NOTCH1', 'CASP8', 'CDK12', 'GRIN2A', 'PTPRC', 'ARHGAP35', 'CHD8',
    'FOXL2', 'BAP1', 'BCL6', 'MET'
])


def compute_p_value(scores, null_p_values):
    """Get the p-value for each score by examining the list null distribution
    where scores are obtained by a certain probability.

    NOTE: uses score2pval function

    Parameters
    ----------
    scores : pd.Series
        series of observed scores
    null_p_values: pd.Series
        Empirical null distribution, index are scores and values are p values

    Returns
    -------
    pvals : pd.Series
        Series of p values for scores
    """
    num_scores = len(scores)
    pvals = pd.Series(np.zeros(num_scores))
    null_p_val_scores = list(reversed(null_p_values.index.tolist()))
    #null_p_values = null_p_values.ix[null_p_val_scores].copy()
    null_p_values.sort_values(inplace=True, ascending=False)
    pvals = scores.apply(lambda x: score2pval(x, null_p_val_scores, null_p_values))
    return pvals


def score2pval(score, null_scores, null_pvals):
    """Looks up the P value from the empirical null distribution based on the provided
    score.

    NOTE: null_scores and null_pvals should be sorted in ascending order.

    Parameters
    ----------
    score : float
        score to look up P value for
    null_scores : list
        list of scores that have a non-NA value
    null_pvals : pd.Series
        a series object with the P value for the scores found in null_scores

    Returns
    -------
    pval : float
        P value for requested score
    """
    # find position in simulated null distribution
    pos = bisect.bisect_right(null_scores, score)

    # if the score is beyond any simulated values, then report
    # a p-value of zero
    if pos == null_pvals.size and score > null_scores[-1]:
        return 0
    # condition needed to prevent an error
    # simply get last value, if it equals the last value
    elif pos == null_pvals.size:
        return null_pvals.iloc[pos-1]
    # normal case, just report the corresponding p-val from simulations
    else:
        return null_pvals.iloc[pos]


def cummin(x):
    """A python implementation of the cummin function in R"""
    for i in range(1, len(x)):
        if x[i-1] < x[i]:
            x[i] = x[i-1]
    return x


def bh_fdr(pval):
    """A python implementation of the Benjamani-Hochberg FDR method.

    This code should always give precisely the same answer as using
    p.adjust(pval, method="BH") in R.

    Parameters
    ----------
    pval : list or array
        list/array of p-values

    Returns
    -------
    pval_adj : np.array
        adjusted p-values according the benjamani-hochberg method
    """
    pval_array = np.array(pval)
    sorted_order = np.argsort(pval_array)
    original_order = np.argsort(sorted_order)
    pval_array = pval_array[sorted_order]

    # calculate the needed alpha
    n = float(len(pval))
    pval_adj = np.zeros(n)
    i = np.arange(1, n+1, dtype=float)[::-1]  # largest to smallest
    pval_adj = np.minimum(1, cummin(n/i * pval_array[::-1]))[::-1]
    return pval_adj[original_order]


def mean_log_fold_change(data, genes):
    """Mean log fold change function

    Parameters
    ----------
    data : pd.Series
        a series of p-values

    Returns
    -------
    mlfc : float
        mean log fold change.
    """
    tmp = data.copy()
    tmp = tmp[~genes.isin(mlfc_remove_genes)]
    tmp.sort_values(ascending=True, inplace=True)
    tmp[tmp==0] = tmp[tmp>0].min()  # avoid infinity in log by avoiding zero pvals
    dist_quant = np.arange(1, len(tmp)+1)/float(len(tmp))
    mlfc = np.mean(np.abs(np.log2(tmp/dist_quant)))
    return mlfc
