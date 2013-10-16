from __future__ import division  # prevents integer division
import pandas as pd
import utils.python
import utils.python.plot as myplt
import utils.python.utils as utils
import logging


def plot_aa_missense_heatmap(file_path='data_analysis/results/aa_change.missense.txt',
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


def plot_aa_property_heatmap(file_path='data_analysis/results/aa_change.properties.txt',
                             save_path='data_analysis/plots/aa_property.heatmap.png'):
    """Plot a heatmap for mutation changes in chemical properties.

    """
    logger = logging.getLogger(__name__)
    # logger.info('reading in %s ...' % file_path)
    # df = pd.read_csv(file_path, sep='\t')  # read in data
    # df = df.set_index('initial_prop')  # set rows as initial property
    # logger.info('finished reading.')
    df = utils.read_aa_properties(file_path)

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


def plot_aa_property_barplot(file_path='data_analysis/results/aa_change.properties.txt',
                             save_path='data_analysis/plots/aa_property.barplot.png'):
    logger = logging.getLogger(__name__)
    df = utils.read_aa_properties(file_path)
    logger.info('Plotting change in chemical property barplot (%s) ...' % save_path)
    myplt.barplot(df,
                  file_path=save_path,
                  title='Amino Acid Missense Mutations by Property',
                  ylabel='Counts')
    logger.info('Finished plotting heatmap of AA chemical barplot.')


def plot_aa_mutation_types_barplot(mutation_cts,
                                   save_path='data_analysis/plots/aa_mut_types.barplot.png'):
    """Create a barplot graphing amino acid mutation type counts.

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
                  title='Amino Acid Mutations by Type',
                  ylabel='Counts')
    logger.info('Finished plotting barplot of AA mutation types.')
