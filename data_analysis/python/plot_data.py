import pandas as pd
import utils.python.plot as myplt
import logging


def plot_aa_missense_change(file_path='data_analysis/results/aa_change.missense.txt',
                            save_path='data_analysis/plots/aa_heatmap.missense.png'):
    """Plot a heatmap for missense mutations.

    Rows are normalize in order to sum to 1. Each cell in the heatmap represents
    the mutation transition probability.

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
