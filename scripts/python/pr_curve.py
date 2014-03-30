import pandas as pd
import numpy as np
from numpy import interp
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
import re
import argparse


def plot_pr_curve(prec_df,
                  recall,
                  save_path,
                  title):
    prec_df = prec_df.set_index(recall)
    prec_df.plot(kind='line')
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


def calc_all_pr_metrics(truth, perf_df, ptypes):
    all_pr_metrics = []
    for i, col in enumerate(perf_df.columns):
        myseries = perf_df[col].dropna()
        all_pr_metrics.append(calc_pr_metrics(truth.ix[myseries.index],
                                              ptypes[i] * myseries))
    return all_pr_metrics


def construct_performance_df(perf_files, header_names, names_list):
    perf_df = None
    for i, f in enumerate(perf_files):
        tmp_df = pd.read_csv(f, sep='\t', index_col=0)
        if not perf_df:
            # initialize df if not yet
            perf_df = pd.DataFrame(index=tmp_df.index)
        perf_df[names_list[i]] = tmp_df[header_names[i]]
    return perf_df


def parse_arguments():
    parser = argparse.ArgumentParser()
    help_str = ('Comma separated list of files containing '
                'performance measures')
    parser.add_argument('-p', '--performance-files',
                        type=str, required=True,
                        help=help_str)
    help_str = ('Use a column name to grab correct column '
                'from performance file (comma separated)')
    parser.add_argument('-hn', '--header-names',
                        type=str, required=True,
                        help=help_str)
    help_str = 'file path only containing true classification'
    parser.add_argument('-t', '--truth',
                        type=str, required=True,
                        help=help_str)
    help_str = ('Name of method used for each performance file'
                ' (comma separated)')
    parser.add_argument('-n', '--names',
                        type=str, required=True,
                        help=help_str)
    parser.add_argument('-s', '--save-path',
                        type=str, required=True,
                        help='path to save pr curve')
    help_str = ('comma separated list of "-1" and "+1", which indicate'
                'performance metric is a p/q-value vs probability, respectively')
    parser.add_argument('-pt', '--performance-type',
                        type=str, required=True,
                        help=help_str)
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # get info from CLI
    perf_files = re.split('\s*,\s*', opts['performance_files'])
    header_names = re.split('\s*,\s*', opts['header_names'])
    names_list = re.split('\s*,\s*', opts['names'])
    with open(opts['truth']) as handle:
        truth = [line.strip() for line in handle.readlines()]
    perf_types = [int(x) for x in re.split('\s*,\s*', opts['performance_type'])]

    # construct all the performance columns to make pr curves
    perf_df = construct_performance_df(perf_files, header_names, names_list)

    # 1/0 series for being in "truth" list (i.e. training list)
    truth_df = perf_df.index.to_series().apply(lambda x: 1 if x in truth else 0)
    truth_df = truth_df.reindex(perf_df.index)

    all_pr_metrics = calc_all_pr_metrics(truth_df, perf_df, perf_types)
    recall = all_pr_metrics[0][1]
    prec_df = pd.DataFrame(index=recall)
    for i, m in enumerate(all_pr_metrics):
        orig_name = perf_df.columns[i]
        new_name = orig_name + ' (AUC={0:.3f})'.format(all_pr_metrics[i][2])
        prec_df[new_name] = all_pr_metrics[i][0]

    plot_pr_curve(prec_df, recall, opts['save_path'], 'Precision-Recall Curve')


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
