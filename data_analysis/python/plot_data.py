from __future__ import division  # prevents integer division
import pandas as pd
import utils.python
import utils.python.plot as myplt
import utils.python.util as _utils
from matplotlib.mlab import PCA
import logging


def aa_missense_heatmap(file_path='data_analysis/results/aa_change.missense.txt',
                        save_path='data_analysis/plots/aa_missense.heatmap.png'):
    """Plot a heatmap for missense mutations.

    Rows are normalize in order to sum to 1. Each cell in the heatmap represents
    the mutation transition probability. The y-axis represents the initial amino
    acid and the x-axis represents the mutated amino acid.

    Kwargs:
        file_path (str): file to data containing missense mutation counts
        save_path (str): file name of heatmap figure
    """
    logger = logging.getLogger(__name__)
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


def aa_property_heatmap(file_path='data_analysis/results/aa_change.properties.txt',
                        save_path='data_analysis/plots/aa_property.heatmap.png'):
    """Plot a heatmap for mutation changes in chemical properties.

    """
    logger = logging.getLogger(__name__)
    # logger.info('reading in %s ...' % file_path)
    # df = pd.read_csv(file_path, sep='\t')  # read in data
    # df = df.set_index('initial_prop')  # set rows as initial property
    # logger.info('finished reading.')
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


def aa_property_barplot(file_path='data_analysis/results/aa_change.properties.txt',
                        save_path='data_analysis/plots/aa_property.barplot.png'):
    logger = logging.getLogger(__name__)
    df = _utils.read_aa_properties(file_path)
    logger.info('Plotting change in chemical property barplot (%s) ...' % save_path)
    myplt.barplot(df,
                  file_path=save_path,
                  title='Amino Acid Missense Mutations by Property',
                  ylabel='Counts')
    logger.info('Finished plotting heatmap of AA chemical barplot.')


def aa_mutation_types_barplot(mutation_cts,
                              save_path='data_analysis/plots/aa_mut_types.barplot.png',
                              title='Protein Mutations by Type'):
    """Create a barplot graphing counts of amino acid mutation types.

    Currently synonymous, missense, nonsense, frame shift, and indels
    are plotted in the bar graph.

    Args:
        mutation_cts (pd.Series): unique counts for mutation types

    Kwargs:
        save_path (str): path to save barplot
    """
    logger = logging.getLogger(__name__)
    logger.info('Plotting AA mutation type counts barplot (%s) . . .' % save_path)
    myplt.barplot(mutation_cts,
                  save_path,
                  title=title,
                  ylabel='Counts')
    logger.info('Finished plotting barplot of AA mutation types.')


def gene_mutation_histogram(gene_cts,
                            save_path='data_analysis/plots/gene_mutations.histogram.png',
                            title='Gene Mutation Histogram'):
    logger = logging.getLogger(__name__)
    logger.info('Plotting gene mutation histogram (%s) . . .' % save_path)
    myplt.histogram(gene_cts,
                    save_path,
                    title=title,
                    ylabel='Counts (log)')
    logger.info('Finished plotting gene mutation histogram.')


def cumulative_gene_mutation(gene_cts,
                             save_path='data_analysis/plots/gene_mutation.cumulative.png',
                             title='Cumulative Gene Mutations'):
    logger = logging.getLogger(__name__)
    logger.info('Plotting cumulative gene mutations (%s) . . .' % save_path)
    df = pd.DataFrame(gene_cts)
    df['pseudo_count'] = 1  # each gene only counts once
    my_counts = df.groupby('count')['pseudo_count'].sum()  # numr of genes for each mutation count
    cumulative_cts = my_counts.cumsum()
    myplt.line(cumulative_cts,
               'data_analysis/plots/gene_mutations.cumulative.png',
               logx=True,
               title='Cumulative Gene Mutations',
               ylabel='Number of Genes',
               xlabel='Number of Gene Mutations (log)',
               vlines=[7, 18])  # vogelstein curates at between 7-18 counts
    logger.info('Finished plotting cumulative gene mutations.')


def pca_plot(file_path='data_analysis/results/gene_design_matrix.txt',
             save_path='data_analysis/plots/gene_pca.scatter.png',
             title='Gene Mutation PCA'):
    logger = logging.getLogger(__name__)
    logger.info('Plotting PCA of gene mutations (%s) . . .' % save_path)

    # normalize counts
    df = pd.read_csv(file_path, sep='\t')
    oncogenes = set(_utils.read_oncogenes())
    tsgs = set(_utils.read_tsgs())
    colors = []
    for g in df['gene']:
        if g in oncogenes:
            colors.append('red')
        elif g in tsgs:
            colors.append('purple')
        else:
            colors.append('blue')
    df = df[df.columns.tolist()[1:]]  # remove gene column
    tmp_total_cts = df.T.sum()
    df = (df.T / tmp_total_cts).T
    df['total_cts'] = tmp_total_cts.T

    results = PCA(df)
    xy_data = [[item[0], item[1]] for item in results.Y]  # first two components
    x, y = zip(*xy_data)
    myplt.scatter(x, y,
                  save_path,
                  colors=colors,
                  title='Mutation PCA',
                  xlabel='1st component',
                  ylabel='2nd component')


def all_mut_type_barplot(df,
                         save_path='data_analysis/plots/all_mut_type.barplot.png',
                         title='Protein Mutation Types by Gene Label'):
    logger = logging.getLogger(__name__)
    logger.info('Plotting protein mutation types by gene type (%s) . . .' % save_path)

    myplt.barplot(df,
                  save_path,
                  title='Protein Mutation Type by Gene Type',
                  ylabel='Counts',
                  stacked=True)

    logger.info('Finished plotting protein mutation types by gene type.')
