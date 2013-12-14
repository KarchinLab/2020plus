from utils.python.amino_acid import AminoAcid
import utils.python
import utils.python.util as _utils
import plot_data
import pandas as pd
import pandas.io.sql as psql
import logging
import csv

logger = logging.getLogger(__name__)

def count_missense(conn, table):
    """Count amino acid changes.

    Parameters
    ----------
    conn : mysqldb connection object
        either could be sqlite or mysql

    Returns
    -------
    aa_change_counter : dict
        contains counts eg. {('aa1', 'aa2'): 4}
    """
    logger.info('Starting to count amino acid changes . . .')

    # perform query
    if table == 'cosmic_aa':
        sql = 'SELECT aachange, occurrences FROM cosmic_aa'
    elif table == 'nucleotide':
        sql = 'SELECT AminoAcid FROM nucleotide'
    df = psql.frame_query(sql, con=conn)

    # count amino acid missense mutations
    aa_change_counter = {}
    for i, row in df.iterrows():
        # get amino acid information
        if table == 'cosmic_aa' :
            aa = AminoAcid(hgvs=row['aachange'],
                           occurrence=row['occurrences'])
        elif table == 'nucleotide':
            aa = AminoAcid(hgvs=row['AminoAcid'])

        # count amino acid
        if aa.is_missense and aa.is_valid and not aa.is_missing_info:
            aa_change_counter.setdefault((aa.initial, aa.mutated), 0)
            aa_change_counter[(aa.initial, aa.mutated)] += aa.occurrence
    logger.info('Finished counting amino acid changes.')
    return aa_change_counter


def save_missense(aacounter,
                  missense_save,
                  property_save):
    """Saves missense mutation counts to file.
    """
    logger.info('Saving protein missense mutation information . . .')

    # save missense mutation counts into a file
    header = [['initial', 'mutated', 'count']]
    aa_list = sorted([[key[0], key[1], val]
                      for key, val in aacounter.iteritems() if "*" not in key])
    csv.writer(open(missense_save, 'wb'),
               delimiter='\t').writerows(header + aa_list)

    # re-slice the mutation data
    df = pd.read_csv(missense_save, sep='\t')
    # add properties of initial/mutated amino acids
    df['initial_prop'] = df['initial'].apply(lambda x: utils.python.letter_to_prop[x])
    df['mutated_prop'] = df['mutated'].apply(lambda x: utils.python.letter_to_prop[x])
    ptable = pd.pivot_table(df,
                            values='count',
                            rows='initial_prop',
                            cols='mutated_prop',
                            aggfunc=sum)
    ptable.to_csv(property_save, sep='\t')

    logger.info('Finished saving protein missense information.')


def main(conn, table):
    cfg_opts = _utils.get_output_config('missense')
    result_dir = _utils.result_dir
    plot_dir = _utils.plot_dir

    # handle missense mutation data
    aa_counter = count_missense(conn, table=table)
    save_missense(aa_counter,
                  missense_save=result_dir + cfg_opts['missense'],
                  property_save=result_dir + cfg_opts['property'])
    plot_data.aa_missense_heatmap(result_dir + cfg_opts['missense'],  # read in path
                                  plot_dir + cfg_opts['missense_heatmap'])  # plot path
    plot_data.aa_property_heatmap(result_dir + cfg_opts['property'],
                                  plot_dir + cfg_opts['property_heatmap'])
    plot_data.aa_property_barplot(result_dir + cfg_opts['property'],
                                  plot_dir + cfg_opts['property_barplot'])
