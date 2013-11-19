from __future__ import division
from vogelstein_classifier import VogelsteinClassifier
import utils.python.util as _utils
import plot_data
import pandas as pd
import numpy as np
import glob
import re


def calc_onco_info(df, onco_pct, tsg_pct):
    # calculate the number of genes classified as oncogene
    vclf = VogelsteinClassifier(onco_pct, tsg_pct)
    df['total'] = df.T.sum()
    input_list = ((row['recurrent missense'],
                   row['frame shift'] + row['nonsense'],
                   row['total'])
                  for i, row in df.iterrows())
    df['2020 class'] = vclf.predict_list(input_list)
    class_cts = df['2020 class'].value_counts()

    # calculate the pct of known oncogenes found
    df['curated class'] = [_utils.classify_gene(gene)
                           for gene in df.index.tolist()]
    tmpdf = df.copy()  # prevent wierd behavior
    known_onco = tmpdf[tmpdf['curated class']=='oncogene']
    num_onco_found = len(known_onco[known_onco['2020 class']=='oncogene'])
    total_onco = len(known_onco)  # total number of oncogenes with counts
    pct_onco_found = num_onco_found / total_onco

    return class_cts['oncogene'], pct_onco_found


def num_onco_by_recurrent_mutations(onco_pct, tsg_pct):
    """Count number of oncogenes while varying the definition of recurrency"""
    # calculate counts for oncogenes/tsg with varying the required the number
    # of mutations to define a recurrent position
    file_match_pattern = './data_analysis/results/genes/gene_design_matrix.r*.txt'
    gene_design_matrix_paths = glob.glob(file_match_pattern)
    onco_ct_list, onco_pct_list = [], []  # list of cts/pct for oncogenes
    for file_path in gene_design_matrix_paths:
        tmp_df = pd.read_csv(file_path, sep='\t', index_col=0)
        tmp_ct, tmp_pct = calc_onco_info(tmp_df, onco_pct=onco_pct, tsg_pct=tsg_pct)
        onco_ct_list.append(tmp_ct)
        onco_pct_list.append(tmp_pct)

    # extract the '-r' parameter from the file name
    recur_param_pattern = '\d+'
    recur_param_list = [int(re.search(recur_param_pattern, mypath).group())
                        for mypath in gene_design_matrix_paths]

    # return dataframe with counts for each use of a recurrent mutation counts
    mycts = pd.Series(onco_ct_list, index=recur_param_list)
    mypct = pd.Series(onco_pct_list, index=recur_param_list)
    return mycts, mypct


def num_onco_by_pct_threshold():
    cts, pct = num_onco_by_recurrent_mutations(.2, .2)
    df_ct = pd.DataFrame(index=cts.index)
    df_pct = pd.DataFrame(index=pct.index)
    for threshold in np.arange(.15, .5, .05):
        tmp_ct, tmp_pct = num_onco_by_recurrent_mutations(threshold, threshold)
        df_ct[str(threshold)] = tmp_ct
        df_pct[str(threshold)] = tmp_pct
    return df_ct, df_pct


def main():
    cfg_opts = _utils.get_output_config('classifier')
    input_opts = _utils.get_output_config('stats')

    # get oncogene info
    count_df, pct_df = num_onco_by_pct_threshold()
    count_df = count_df.sort_index()  # sort df by recurrent mutation cts
    pct_df = pct_df.sort_index()  # sort df by recurrent mutation cts

    # save results
    count_df.to_csv(_utils.clf_result_dir + cfg_opts['oncogene_parameters_ct'], sep='\t')
    pct_df.to_csv(_utils.clf_result_dir + cfg_opts['oncogene_parameters_pct'], sep='\t')

    # plot results
    # plot number of predicted oncogenes while varying parameters
    tmp_save_path = _utils.clf_plot_dir + cfg_opts['number_oncogenes_plot']
    tmp_title = 'Vogelstein\'s Classifier Predicted Oncogenes'
    tmp_ylabel = 'Number of Oncogenes'
    tmp_xlabel = 'Number of Mutations Required for Recurrency'
    plot_data.onco_mutations_parameter(count_df,
                                       tmp_save_path,
                                       title=tmp_title,
                                       ylabel=tmp_ylabel,
                                       xlabel=tmp_xlabel)
    # plot percentage of vogelstein's oncogenes recovered
    tmp_title = 'Percentage of Vogelstein\'s Oncogenes Recovered'
    tmp_ylabel = 'Fraction of Oncogenes'
    tmp_xlabel = 'Number of Mutations Required for Recurrency'
    tmp_save_path = _utils.clf_plot_dir + cfg_opts['pct_oncogenes_plot']
    plot_data.onco_mutations_parameter(pct_df,
                                       tmp_save_path,
                                       title=tmp_title,
                                       ylabel=tmp_ylabel,
                                       xlabel=tmp_xlabel)


if __name__ == "__main__":
    main()
