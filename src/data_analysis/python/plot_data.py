"""
Performs necessary tasks prior to plotting data via the
utils/python/plot.py module. In many cases the plot functions
read the data from file so that data observed in the results
directory is always consistent with the actual plots.
"""

from __future__ import division  # prevents integer division
import pandas as pd
import src.utils.python
import numpy as np
import src.utils.python.plot as myplt
import src.utils.python.util as _utils
from matplotlib.mlab import PCA
import matplotlib.pyplot as plt
import logging

logger = logging.getLogger(__name__)  # logger obj for this module

def aa_missense_heatmap(file_path, save_path):
    """Plot a heatmap for missense mutations.

    Rows are normalize in order to sum to 1. Each cell in the heatmap represents
    the mutation transition probability. The y-axis represents the initial amino
    acid and the x-axis represents the mutated amino acid.

    **Parameters**

    file_path : str
        file to data containing missense mutation counts
    save_path : str
        file name of heatmap figure
    """
    logger.info('reading in %s ...' % file_path)
    df = pd.read_csv(file_path, sep='\t')  # read in data
    logger.info('finished reading.')

    # pivot data to create a mutation count matrix
    ptable = pd.pivot_table(df,
                            values='count',
                            rows='initial',
                            cols='mutated',
                            aggfunc=sum)

    # normalize rows to sum to 1 (transition prob. matrix)
    ptable_norm = (ptable.T / ptable.T.sum()).T
    ptable_norm.fillna(0)  # fill missing with 0 probability

    # reorder rows/columns to be in chemically meaningful order for AA
    order = ['A', 'C', 'G', 'I', 'L', 'M', 'F', 'P', 'W', 'V',  # nonpolar
             'N', 'Q', 'S', 'T', 'Y',  # neutral polar
             'R', 'H', 'K',  # basic polar
             'D', 'E']  # acidic polar
    ptable_norm = ptable_norm.ix[order]
    ptable_norm = ptable_norm[order]

    # plot and save heatmap figure
    logger.info('Plotting missense mutation heatmap (%s) ...' % save_path)
    myplt.heatmap(ptable_norm,
                  file_path=save_path,
                  xlabel='Mutated AA',
                  ylabel='Initial AA')
    logger.info('Finished plotting heatmap.')


def nuc_substitution_heatmap(file_path, save_path, title=''):
    """Plot a heatmap for DNA substiution mutations.

    Rows are normalize in order to sum to 1 (legal probability). Each cell in the
    heatmap represents the mutation transition probability. The y-axis represents the
    initial nucleotide and the x-axis represents the mutated nucleotide.

    **Parameters**

    file_path : str
        file to data containing substiution mutation counts
    save_path : str
        file name of heatmap figure
    title : str, (default='')
        plot title
    """
    logger.info('reading in %s ...' % file_path)
    df = pd.read_csv(file_path, sep='\t')  # read in data
    logger.info('finished reading.')

    # pivot data to create a mutation count matrix
    ptable = pd.pivot_table(df,
                            values='count',
                            rows='initial',
                            cols='mutated',
                            aggfunc=sum)

    # normalize rows to sum to 1 (transition prob. matrix)
    ptable_norm = (ptable.T / ptable.T.sum()).T
    ptable_norm.fillna(0)  # fill missing with 0 probability

    # plot and save heatmap figure
    logger.info('Plotting substitution mutation heatmap (%s) ...' % save_path)
    myplt.heatmap(ptable_norm,
                  file_path=save_path,
                  title=title,
                  xlabel='Mutated Base',
                  ylabel='Initial Base')
    logger.info('Finished plotting heatmap.')


def aa_property_heatmap(file_path, save_path):
    """Plot a heatmap for mutation changes in chemical properties.
    """
    df = _utils.read_aa_properties(file_path)

    # normalize rows to sum to 1 (transition prob. matrix)
    df_norm = (df.T / df.T.sum()).T
    df_norm.fillna(0)  # fill missing with 0 probability

    # reorder rows/columns to go from non-polar to polar
    order = ['nonpolar', 'polar', 'basic polar', 'acidic polar']
    df_norm = df_norm.ix[order]
    df_norm = df_norm[order]

    # plot and save heatmap figure
    logger.info('Plotting change in chemical property heatmap (%s) ...' % save_path)
    myplt.heatmap(df_norm,
                  file_path=save_path,
                  xlabel='Mutated AA Properties',
                  ylabel='Initial AA Properties')
    logger.info('Finished plotting heatmap of AA chemical properties.')


def aa_property_barplot(file_path, save_path):
    df = _utils.read_aa_properties(file_path)
    logger.info('Plotting change in chemical property barplot (%s) ...' % save_path)
    myplt.barplot(df,
                  file_path=save_path,
                  title='Amino Acid Missense Mutations by Property',
                  ylabel='Counts')
    logger.info('Finished plotting heatmap of AA chemical barplot.')


def nuc_substitution_barplot(file_path, save_path,
                             title='DNA Substitution Mutations'):
    df = pd.read_csv(file_path, sep='\t')

    # pivot data to create a mutation count matrix
    ptable = pd.pivot_table(df,
                            values='count',
                            rows='initial',
                            cols='mutated',
                            aggfunc=sum)

    logger.info('Plotting substitution barplot (%s) ...' % save_path)
    myplt.barplot(ptable,
                  file_path=save_path,
                  title=title,
                  ylabel='Counts',
                  stacked=True)
    logger.info('Finished plotting bar plot.')


def mutation_types_barplot(mutation_cts,
                           save_path,
                           title='Mutations by Type'):
    """Create a barplot graphing counts of amino acid/DNA mutation types.

    Currently synonymous, missense, nonsense, frame shift, and indels
    are plotted for amino acids in the bar graph.

    **Parameters**

    mutation_cts : pd.Series
        unique counts for mutation types
    save_path : str
        path to save barplot
    title : str
        title for plot
    """
    logger.info('Plotting mutation type counts barplot (%s) . . .' % save_path)
    myplt.barplot(mutation_cts,
                  save_path,
                  title=title,
                  ylabel='Counts')
    logger.info('Finished plotting barplot of mutation types.')


def gene_mutation_histogram(gene_cts,
                            save_path,
                            title='Gene Mutation Histogram'):
    logger.info('Plotting gene mutation histogram (%s) . . .' % save_path)
    myplt.histogram(gene_cts,
                    save_path,
                    bins=range(0, 500, 10),  # not many genes >300
                    log=True,  # log scale y-axis
                    title=title,
                    ylabel='Counts (log)')
    logger.info('Finished plotting gene mutation histogram.')


def cumulative_gene_mutation(gene_cts,
                             save_path,
                             title='Cumulative Gene Mutations'):
    logger.info('Plotting cumulative gene mutations (%s) . . .' % save_path)
    df = pd.DataFrame(gene_cts)
    df['pseudo_count'] = 1  # each gene only counts once
    my_counts = df.groupby('count')['pseudo_count'].sum()  # numr of genes for each mutation count
    cumulative_cts = my_counts.cumsum()
    myplt.line(cumulative_cts,
               save_path,
               logx=True,
               title='Cumulative Gene Mutations',
               ylabel='Number of Genes',
               xlabel='Number of Gene Mutations (log)',
               vlines=[7, 18])  # vogelstein curates at between 7-18 counts
    logger.info('Finished plotting cumulative gene mutations.')


def pca_plot(file_path,
             save_path,
             norm_class=False,
             low_count_filter=1,
             title='Gene Mutation PCA'):
    """Create a PCA plot using features from genes.

    Each point represents one gene. The features for each gene
    are reduced to two components for plotting.

    Parameters
    ----------
    file_path : str
        path to file for gene feature matrix
    save_path : str
        path to save plot
    norm_class : int, (Default=3)
        ratio of driver to passenger genes for subsampling genes to even unbalanced classes.
    low_count_filter : int, (Default=1)
        Genes should have at least low_count_filter number of mutations
    title : str
        title for plot
    """
    logger.info('Plotting PCA of gene mutations (%s) . . .' % save_path)

    # normalize counts
    df = pd.read_csv(file_path,  # path
                     sep='\t',  # tab delim
                     index_col=0)  # index df by gene name

    if norm_class:
        # non-driver genes (i.e. not oncogenes or tsg) heavily
        # out numbered driver genes. Thus PCA will tend to explain
        # variance for essentially only the non-driver genes. tolist
        # counter this effect if "norm_class" is True then non-driver
        # genes will be sub-sampled
        driver_list = list(_utils.oncogene_list + _utils.tsg_list)
        driver_df = df.ix[driver_list]
        non_driver_list = list(set(df.index) - set(driver_list))
        other_df = df.ix[non_driver_list]
        len_driver = len(driver_df)
        sub_sample = other_df.loc[np.random.choice(other_df.index, norm_class*len_driver, replace=False)]
        df = pd.concat([driver_df, sub_sample])  # new df with more equal classes

    # plot oncogenes and tumor suppressor genes as different colors
    oncogenes = set(_utils.read_oncogenes())  # get oncogenes
    tsgs = set(_utils.read_tsgs())  # get tumor suppressor genes
    colors = []
    for g in df.index.tolist():
        if g in oncogenes:
            colors.append('red')
        elif g in tsgs:
            colors.append('purple')
        else:
            colors.append('blue')

    # normalize data by row for PCA
    row_sums = df.sum(axis=1)
    df = df.div(row_sums.astype(float), axis=0)
    keep_rows = row_sums >= low_count_filter
    df = df.ix[keep_rows]  # filter low mutation ct genes

    # get marker size for scatter plot
    MAX_SIZE = 300  # some genes take up to much space
    scatter_size = [size if size < MAX_SIZE else MAX_SIZE for size in row_sums]

    # drop columns that are potentially all zeros
    # since this will make SVD not converge.
    # all zeros could come from not using COSMIC
    # data
    drop_cols = []
    for c in df.columns:
        if df[c].sum() == 0:
            drop_cols.append(c)
    df = df.drop(drop_cols, axis=1)

    # perform PCA
    results = PCA(df)
    first_eigen_value, second_eigen_value = results.fracs[:2]
    xy_data = [[item[0], item[1]] for item in results.Y]  # first two components
    x, y = zip(*xy_data)
    myplt.scatter(x, y,
                  save_path,
                  colors=colors,
                  size=scatter_size,
                  title=title,
                  xlabel='1st component (%f)' % first_eigen_value,
                  ylabel='2nd component (%f)' % second_eigen_value)
    logger.info('Finished PCA plot.')


def all_mut_type_barplot(df,
                         save_path,
                         title='Protein Mutation Types by Gene Label'):
    logger.info('Plotting protein mutation types by gene type (%s) . . .' % save_path)
    myplt.barplot(df,
                  save_path,
                  title='Protein Mutation Type by Gene Type',
                  ylabel='Counts',
                  stacked=True)
    logger.info('Finished plotting protein mutation types by gene type.')


def non_silent_ratio_kde(df, save_path,
                         xlim=None,
                         title='',
                         xlabel='',
                         ylabel=''):
    labels = df['label'].unique()
    for label in labels:
        df[df['label']==label]['non-silent/silent'].plot(kind='kde', label=label)

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc='best')
    plt.savefig(save_path)
    plt.close()


def recurrent_missense_pos_line(df,
                                save_path,
                                xlim=(2, 25)):
    logger.info('Plotting number of recurrent missense positions (%s) ...' % save_path)
    # filter data
    df = df.fillna(0)  # fill absent counts with zeros
    df = df[(df.index>=xlim[0]) & (df.index<xlim[1])]  # only include certain range

    # normalized based on number of genes
    num_onco = len(_utils.oncogene_list)
    num_tsg = len(_utils.tsg_list)
    num_other = 19140
    df['oncogene'] = df['oncogene'] / num_onco
    df['tsg'] = df['tsg'] / num_tsg
    df['other'] = df['other'] / num_other

    myplt.line(df,
               save_path,
               title='Number of Recurrent Missense Positions per Gene',
               xlabel='Number of Mutations per Position',
               ylabel='Number of Recurrent Missense Positions per Gene')
    logger.info('Finished plotting number of recurrent missense positions')


def entropy_kde(df,
                column,
                save_path,
                title='',
                xlabel='',
                ylabel='Density'):
    df[df['true class']==_utils.tsg_label][column].dropna().plot(kind='kde', label='TSG')
    df[df['true class']==_utils.onco_label][column].dropna().plot(kind='kde', label='Oncogenes')
    df[df['true class']==_utils.other_label][column].dropna().plot(kind='kde', label='Other genes')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc='best')
    plt.savefig(save_path)
    plt.close()


def sample_barplot(df, save_path, title='', xlabel='', ylabel=''):
    logger.info('plotting number of mutation give sample size . . .')
    df['category'] = ''
    tmp_ix = df[df['MutationCounts']>100].index
    df['category'].ix[tmp_ix] = '100+'
    tmp_ix = df[(df['MutationCounts']>75) & (df['MutationCounts']<=100)].index
    df['category'].ix[tmp_ix] = '75-100'
    tmp_ix = df[(df['MutationCounts']>50) & (df['MutationCounts']<=75)].index
    df['category'].ix[tmp_ix] = '50-75'
    tmp_ix = df[(df['MutationCounts']>25) & (df['MutationCounts']<=50)].index
    df['category'].ix[tmp_ix] = '25-50'
    tmp_ix = df[(df['MutationCounts']>=1) & (df['MutationCounts']<=25)].index
    df['category'].ix[tmp_ix] = '1-25'
    category_cts = df.groupby('category').aggregate(sum)
    myplt.barplot(category_cts,
                  save_path,
                  title=title,
                  ylabel='Counts')
    logger.info('Finished plotting number of mutation given sample size.')


def sample_kde(df, save_path,
               xlabel='',
               ylabel='',
               title=''):
    logger.info('Plotting KDE of percent of samples with a non-silent mutation . . .')

    df = df.rename(columns={'non-silent sample pct': 'sample_pct'}).copy()  # plotting expects col named 'sample_pct'

    # categorize genes into onco/tsg
    df['true class'] = df.index.to_series().apply(_utils.classify_gene)
    df['true class'] = df['true class'].apply(lambda x: _utils.class_to_label[x])

    # plot kde
    df[df['true class']==_utils.tsg_label]['sample_pct'].dropna().plot(kind='kde', label='TSG')
    df[df['true class']==_utils.onco_label]['sample_pct'].dropna().plot(kind='kde', label='Oncogenes')
    df[df['true class']==_utils.other_label]['sample_pct'].dropna().plot(kind='kde', label='Other genes')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc='best')
    plt.savefig(save_path)
    plt.close()
    logger.info('Finished plotting KDE of samples.')


def sample_boxplot(df, save_path,
                   xlabel='',
                   ylabel='',
                   title=''):
    logger.info('Plotting box plot of percent of samples with specific '
                'mutations ({0}) . . .'.format(save_path))
    # categorize genes into onco/tsg
    df['true class'] = df.index.to_series().apply(_utils.classify_gene)

    # specify text for plot
    if not xlabel:
        xlabel = 'Category of Gene in Training'
    if not ylabel:
        ylabel = 'Maximum Pct of Samples in a Tumor Type'
    if not title:
        title = 'Pct of Samples containing a Non-silent Mutation'

    # plot boxplot
    myplt.boxplot(df,
                  by='true class',
                  column=[],  # use default column for boxplot
                  save_path=save_path,
                  xlabel=xlabel,
                  ylabel=ylabel,
                  title=title)
    logger.info('Finished box plot.')


def entropy_sample_correlation(x, y,
                               save_path,
                               xlabel='',
                               ylabel='',
                               title=''):
    logger.info('Plotting correlation between sample prevalence and '
                ' oncogene position entropy (%s) ...' % save_path)

    # define plot text if not provided
    if not xlabel:
        xlabel = 'Maximum percent of samples a gene is mutated in a given Tumor Type'
    if not ylabel:
        ylabel = 'Percent of Uniform Mutation Entropy'
    if not title:
        title = 'Oncogene comparison between prevalence and positional distribution'

    # make scatter plot
    myplt.correlation_plot(x, y,
                           save_path,
                           title=title,
                           xlabel=xlabel,
                           ylabel=ylabel)
    logger.info('Finished plotting correlation between sample prevalence and '
                ' oncogene position entropy.')


def non_silent_tumor_type_barplot(df, save_path,
                                  xlabel='',
                                  ylabel='',
                                  title=''):
    # get rid of underscores
    df.rename(columns=lambda x: x.replace('_', ' '), inplace=True)
    df.rename(index=lambda x: x.replace('_', ' '), inplace=True)

    # specify labels
    if not xlabel:
        xlabel='Tumor Types'
    if not ylabel:
        ylabel='Number of non-silent mutations'
    if not title:
        title='Number of non-silent mutations for each tumor type'
    # make bar plot
    myplt.barplot(df,
                  save_path,
                  xlabel=xlabel,
                  ylabel=ylabel,
                  title=title)


def js_distance_kde(df, save_path,
                    xlabel='',
                    ylabel='',
                    title=''):
    # initialize labels
    if not xlabel:
        xlabel = 'Jensen-Shannon Distance'
    if not ylabel:
        ylabel = 'Density'
    if not title:
        title = 'JS distance from Database Tumor Type distribution'

    # plot kde
    df[df['true class']=='tsg']['JS distance'].dropna().plot(kind='kde', label='TSG')
    df[df['true class']=='oncogene']['JS distance'].dropna().plot(kind='kde', label='Oncogenes')
    df[df['true class']=='other']['JS distance'].dropna().plot(kind='kde', label='Other genes')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc='best')
    plt.savefig(save_path)
    plt.close()
