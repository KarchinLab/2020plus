import pandas as pd
import numpy as np
from numpy import interp
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
import argparse


def plot_pr_curve(lasso_prec, lasso_auc,
                  tuson_prec, tuson_auc,
                  recall,
                  save_path,
                  title):
    plot_df = pd.DataFrame({'Davoli \\textit{et al.} LASSO (AUC = %.3f)' % lasso_auc: lasso_prec,
                            'Davoli \\textit{et al.} TUSON (AUC = %.3f)' % tuson_auc: tuson_prec},
                           index=recall)
    plot_df.plot(kind='line')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(title)
    plt.savefig(save_path)


def calc_pr_metrics(truth_df, score_df):
    recall_array = np.linspace(0, 1, 100)
    p, r, thresh = metrics.precision_recall_curve(truth_df, score_df)
    p, r, thresh = p[::-1], r[::-1], thresh[::-1]  # reverse order of results
    thresh = np.insert(thresh, 0, 1.0)
    precision_array = interp(recall_array, r, p)
    threshold_array = interp(recall_array, r, thresh)
    pr_auc = metrics.auc(recall_array, precision_array)
    return precision_array, recall_array, pr_auc


def read_gene_list(file_path):
    gene_list = []
    with open(file_path) as handle:
        for line in handle:
            gene_list.append(line.strip())
    return gene_list


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-op', '--operformance',
                        type=str, required=True,
                        help='File containing Davoli p-values and'
                        ' Lasso probabilities for oncogenes')
    parser.add_argument('-og', '--oncogenes',
                        type=str, required=True,
                        help='List of "known" oncogenes')
    parser.add_argument('-oo', '--oncogene-output',
                        type=str, required=True,
                        help='Oncogene performance plot')
    parser.add_argument('-tp', '--tperformance',
                        type=str, required=True,
                        help='File containing Davoli p-values and'
                        ' Lasso probabilities for tsg')
    parser.add_argument('-tg', '--tsg',
                        type=str, required=True,
                        help='List of "known" tumor suppressor genes')
    parser.add_argument('-to', '--tsg-output',
                        type=str, required=True,
                        help='TSG performance plot')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in oncogene information
    df = pd.read_csv(opts['operformance'], sep='\t')
    known_og = read_gene_list(opts['oncogenes'])
    df['known og'] = df['Gene'].apply(lambda x: 1 if x in known_og else 0)

    # get lasso oncogene pr curve
    lasso_precision_array, recall_array, lasso_pr_auc = calc_pr_metrics(df['known og'],
                                                                             df['OG_Probability_LASSO'])
    # get tuson oncogene pr curve
    tuson_precision_array, recall_array, tuson_pr_auc = calc_pr_metrics(df['known og'], -df['TUSON_q_value_OG'])

    # plot results
    plot_pr_curve(lasso_precision_array, lasso_pr_auc,
                  tuson_precision_array, tuson_pr_auc,
                  recall_array,
                  save_path=opts['oncogene_output'],
                  title=r'Davoli \textit{et al.} Oncogene Precision-Recall Curve')

    # read in tumor suppressor
    df = pd.read_csv(opts['tperformance'], sep='\t')
    known_tsg = read_gene_list(opts['tsg'])
    df['known tsg'] = df['Gene'].apply(lambda x: 1 if x in known_tsg else 0)

    # get lasso oncogene pr curve
    lasso_precision_array, recall_array, lasso_pr_auc = calc_pr_metrics(df['known tsg'],
                                                                        df['TSG_Probability_LASSO'])
    # get tuson oncogene pr curve
    tuson_precision_array, recall_array, tuson_pr_auc = calc_pr_metrics(df['known tsg'], -df['TUSON_q_value_TSG'])

    plot_pr_curve(lasso_precision_array, lasso_pr_auc,
                  tuson_precision_array, tuson_pr_auc,
                  recall_array,
                  save_path=opts['tsg_output'],
                  title=r'Davoli \textit{et al.} TSG Precision-Recall Curve')



if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
