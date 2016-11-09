import src.utils.python.plot as myplt
import src.utils.python.util as _utils
import src.utils.python.p_value as pval
import matplotlib
matplotlib.use('agg', warn=False)
import matplotlib.pyplot as plt
import pandas as pd
import logging
import numpy as np

logger = logging.getLogger(__name__)

# genes to remove for MLFC calculation


def feature_importance_barplot(mean_df,
                               std_df,
                               save_path):
    """Creates a feature import barplot based on the mean decrease in the
    gini index."""
    logger.info('Plotting feature importances from random forest . . .')

    # rename columns so latex doesn't complain about '_'
    rename_dict = {col: col.replace('_', ' ') for col in mean_df.index}
    mean_df.rename(rename_dict, inplace=True)
    std_df.rename(rename_dict, inplace=True)

    # rename to fit Table S1
    rename_dict = {
        'silent': 'silent fraction',
        'nonsense': 'nonsense fraction',
        'splice site': 'splice site fraction',
        'missense': 'missense fraction',
        'recurrent missense': 'recurrent missense fraction',
        'frameshift indel': 'frameshift indel fraction',
        'inframe indel': 'inframe indel fraction',
        'lost start and stop': 'lost start and stop fraction',
        'normalized missense position entropy': 'normalized missense position entropy',
        'missense to silent': 'missense to silent',
        'non-silent to silent': 'non-silent to silent',
        'normalized mutation entropy': 'normalized mutation entropy',
        'Mean Missense MGAEntropy': 'mean missense MGAEntropy',
        'Mean VEST Score': 'mean VEST score',
        'inactivating p-value': 'inactivating SNV fraction p-value',
        'entropy p-value': 'missense entropy p-value',
        'vest p-value': 'missense VEST p-value',
        'combined p-value': 'missense combined p-value',
        'gene degree': 'gene degree',
        'gene betweenness': 'gene betweenness centrality',
        'gene length': 'gene length',
        'expression CCLE': 'gene expression CCLE',
        'replication time': 'replication time',
        'HiC compartment': 'HiC compartment',
    }
    mean_df.rename(rename_dict, inplace=True)
    std_df.rename(rename_dict, inplace=True)

    mean_df.sort_values(inplace=True, ascending=True)  # sort with most important features first
    std_df = std_df.ix[mean_df.index]  # match ordering in bar plot
    myplt.barplot(mean_df,
                  save_path,
                  kind='barh',
                  xerr=std_df,
                  ecolor='#00008B',  # dark blue
                  title='Feature Importance in Random Forest',
                  xlabel='Feature Importance')
    logger.info('Finished plotting feature importance bar plot.')


def precision_recall_curve(df,
                           save_path,
                           style,
                           sem=None,
                           title='Precision-Recall Curve',
                           xlabel='Recall',
                           ylabel='Precision'):
    logger.info('Plotting precision-recall curve (%s) ...' % save_path)
    if not sem:
        # normal line plot
        myplt.line(df, save_path,
                   style=style,
                   title=title,
                   xlabel=xlabel,
                   ylabel=ylabel)
    else:
        # include standard error of the mean (sem) in the plot
        myplt.line_fill_between(df, sem, save_path,
                                style=style,
                                title=title,
                                xlabel=xlabel,
                                ylabel=ylabel)
    logger.info('Finished plotting PR curve.')


def receiver_operator_curve(df,
                            save_path,
                            style,
                            sem=None,
                            title='ROC Curve',
                            xlabel='False-Positive Rate',
                            ylabel='True-Positive Rate'):
    logger.info('Plotting receiver operator curve (%s) ...' % save_path)
    if not sem:
        # normal line plot
        myplt.line(df, save_path,
                   style=style,
                   title='ROC Curve',
                   xlabel='False-Positive Rate',
                   ylabel='True-Positive Rate')
    else:
        # include standard error of the mean (sem) in plot
        myplt.line_fill_between(df, sem, save_path,
                                title=title,
                                xlabel=xlabel,
                                ylabel=ylabel)
    logger.info('Finished plotting ROC curve.')


def prob_kde(df, col_name, save_path,
             title,
             xlabel='Probability'):
    df[df['training list class']==_utils.tsg_label][col_name].plot(kind='kde', label='TSG')
    df[df['training list class']==_utils.onco_label][col_name].plot(kind='kde', label='Oncogenes')
    df[df['training list class']==_utils.other_label][col_name].plot(kind='kde', label='Other genes')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.legend()
    plt.savefig(save_path)
    plt.close()


def prob_scatter(df, plot_path, title):
    """
    Scatter plot of oncogene versus tsg probabilities.

    **Parameters**

    df : pd.DataFrame
        results from random forest prediction
    plot_path : str
        path to save scatter plot of tsg/onco probabilities
    title : str
        title of scatter plot
    """
    # scatter plot of oncogene/tsg probabilities
    myplt.scatter(df['oncogene score'],
                  df['tsg score'],
                  plot_path,
                  xlabel='Oncogene Score',
                  ylabel='TSG Score',
                  title=title,
                  colors='#348ABD')


def sample_boxplot(pred_onco,
                   pred_tsg,
                   pred_driver,
                   save_path_type,
                   save_path_driver,
                   xlabel='',
                   ylabel='',
                   title=''):
    """Create a box plot for distribution of percentage of tumor samples
    containing a non-silent mutation in different categories of genes (ie
    oncogenes, tsgs, and drivers).

    Parameters
    ----------
    pred_onco : list
        list of genes predicted as oncogenes
    pred_tsg : list
        list of genes predicted as tsgs
    pred_driver : list
        list of genes predicted as drivers
    save_path_type : str
        path to save figure for comparing oncogenes, tsgs, and other
    save_path_driver : str
        path to save figure for comparing drivers vs other
    xlabel : str
        x-axis label
    ylabel : str
        y-axis label
    title : str
        title of figures
    """
    cfg = _utils.get_output_config('sample')
    df = pd.read_csv(_utils.result_dir + cfg['max_gene_pct_sample_out'],
                     sep='\t', index_col=0)
    df['Predicted Type'] = [("oncogene" if g in pred_onco else "other") for g in df.index]
    df.ix[pred_tsg, 'Predicted Type'] = 'TSG'
    df['Predicted Driver'] = [("driver" if g in pred_driver else "other") for g in df.index]

    # set figure labels
    if not xlabel:
        xlabel = 'Predicted Type'
    if not ylabel:
        ylabel = 'Maximum Pct of Samples for a Tumor Type'
    if not title:
        title = 'Percentage of Samples with Non-Silent Mutation'

    # plot with oncogenes, tsgs, and other
    myplt.boxplot(df,
                  by='Predicted Type',
                  column=['all mutation sample pct', 'non-silent sample pct'],
                  save_path=save_path_type,
                  xlabel=xlabel,
                  ylabel=ylabel,
                  title=title)

    # plot with drivers vs other
    myplt.boxplot(df,
                  by='Predicted Driver',
                  column=['all mutation sample pct', 'non-silent sample pct'],
                  save_path=save_path_driver,
                  xlabel=xlabel,
                  ylabel=ylabel,
                  title=title)


def qqplot(data,
           ax=None, log=False, title=None,
           use_xlabel=True, use_ylabel=True,
           **kwargs):
    """Function for qq-plot with uniform distribution.

    Parameters
    ----------
    data : pd.Series
        p-values for a method
    ax : matplotlib axis
        provided matplotlib axis object to plot on
    log : bool
        indicator to use log scale
    """
    # sort p-values
    tmp = data.copy()
    tmp.sort_values(inplace=True)

    # expected p-values
    dist_quant = np.arange(1, len(tmp)+1)/float(len(tmp)+1)
    if log:
        log_quant = -np.log10(dist_quant)
        if ax is None:
            plt.plot(log_quant, -np.log10(tmp),'o', markersize=3, **kwargs)
            plt.plot([0, log_quant[0]], [0, log_quant[0]], ls="-", color='red')
        else:
            ax.plot(log_quant, -np.log10(tmp),'o', markersize=3, **kwargs)
            ax.plot([0, log_quant[0]], [0, log_quant[0]], ls="-", color='red')
        # set axis labels
        if use_xlabel:
            if ax is None: plt.xlabel('Theoretical ($-log_{10}(p)$)')
            else: ax.set_xlabel('Theoretical ($-log_{10}(p)$)')
        if use_ylabel:
            if ax is None: plt.ylabel('Observed ($-log_{10}(p)$)')
            else: ax.set_ylabel('Observed ($-log_{10}(p)$)')
    else:
        if ax is None:
            plt.plot(dist_quant, tmp,'o', markersize=3, **kwargs)
            plt.plot([0, 1], [0, 1], ls="-", color='red')
        else:
            ax.plot(dist_quant, tmp,'o', markersize=3, **kwargs)
            ax.plot([0, 1], [0, 1], ls="-", color='red')
            ax.set_ylabel('p-value')
        if use_xlabel:
            if ax is None: plt.xlabel('Theoretical p-value')
            else: ax.set_xlabel('Theoretical p-value')
        if use_ylabel:
            if ax is None: plt.ylabel('Observed p-value')
            else: ax.set_ylabel('Observed p-value')
    if title:
        ax.set_title(title)

    return ax


def create_qqplots(pval_df, pval_col, save_path):
    """Create qq plots for oncogene, tsg, and driver p-value."""
    keep_cols = ['gene', pval_col]
    plot_df = pval_df[keep_cols].copy()
    plot_df = plot_df[~plot_df['gene'].isin(pval.mlfc_remove_genes)]
    fig, ax = plt.subplots(1, 1)
    qqplot(pval_df[pval_col], ax)
    mlfc = pval.mean_log_fold_change(plot_df[pval_col], plot_df['gene'])
    ax.text(.025, .95, 'MLFC = {0:.2f}'.format(mlfc))
    ax.set_xlim((0, 1))
    ax.set_ylim((0, 1))
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()
