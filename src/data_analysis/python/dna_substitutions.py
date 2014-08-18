"""This module counts mutational changes at the DNA level.
"""
import utils.python.util as _utils
from utils.python.nucleotide import Nucleotide
import plot_data
import pandas.io.sql as psql
import csv
import itertools as it
import logging

logger = logging.getLogger(__name__)  # module logger

def count_substitutions(hgvs_iterable):
    """Count nucleotide substitutions from HGVS strings.

    Parameters
    ----------
    hgvs_iterable : iterable
        iterable container of HGVS DNA strings

    Returns
    -------
    nuc_counter : dict
        keys are substitutions and values are counts
    """
    # count DNA substitutions
    nuc_counter = {}  # {(initial, mutated): count,...}
    valid_letters = ['A', 'C', 'G', 'T']

    # initiliaze all substitutions to zero
    for possible_sub in it.permutations(valid_letters, 2):
        nuc_counter.setdefault(possible_sub, 0)

    # count actual occurences of substitutions
    for nuc_change in hgvs_iterable:
        nuc = Nucleotide(nuc_change)
        if nuc.is_valid and not nuc.is_missing_info and nuc.is_substitution:
            # only take valid substitution events
            if len(nuc.initial) == 1 and len(nuc.mutated) == 1:
                # restrict attention to single nucleotide substitutions
                if nuc.initial in valid_letters and nuc.mutated in valid_letters:
                    # unfortunately, there are a couple counts from non DNA
                    # which need to be filtered
                    tmp_key = (nuc.initial, nuc.mutated)  # reference => mutated
                    nuc_counter.setdefault(tmp_key, 0)
                    nuc_counter[tmp_key] += 1
    return nuc_counter


def count_nuc_substitutions(conn):
    """Count specific single nucleotide substitutions."""
    logger.info('Starting to count DNA substitutions . . .')

    # query `nucleotide` table
    sql = "SELECT DNA_Change as Nucleotide FROM mutation"
    df = psql.frame_query(sql, con=conn)

    # count DNA substitutions
    nuc_counter = count_substitutions(df['Nucleotide'])
    logger.info('Finished counting DNA substitutions.')
    return nuc_counter


def save_nuc_substitutions(nuc_counter, save_as):
    """Saves DNA substitutions to a file."""
    header = [['initial', 'mutated', 'count']]
    nuc_list = sorted([[key[0], key[1], val]
                       for key, val in nuc_counter.iteritems()])
    csv.writer(open(save_as, 'wb'),
               delimiter='\t').writerows(header + nuc_list)


def count_oncogene_substitutions(conn):
    """Count DNA substitutions in oncogenes."""
    logger.info('Counting oncogene substitution mutations . . .')

    # prepare sql statement
    oncogenes = _utils.oncogene_list
    sql = "SELECT DNA_Change as Nucleotide FROM mutation WHERE Gene in " + str(oncogenes)
    logger.debug('Oncogene SQL statement: ' + sql)

    df = psql.frame_query(sql, con=conn)  # execute query
    nuc_count = count_substitutions(df['Nucleotide'])
    logger.info('Finished counting oncogene substitutions.')
    return nuc_count


def count_tsg_substitutions(conn):
    """Count DNA substitutions in tumor suppressor genes."""
    logger.info('Counting tsg substitution mutations . . .')

    # prepare sql statement
    tsgs = _utils.tsg_list
    sql = "SELECT DNA_Change as Nucleotide FROM mutation WHERE Gene in " + str(tsgs)
    logger.debug('Oncogene SQL statement: ' + sql)

    df = psql.frame_query(sql, con=conn)  # execute query
    nuc_count = count_substitutions(df['Nucleotide'])
    logger.info('Finished counting tsg substitutions.')
    return nuc_count


def main(conn):
    out_dir = _utils.result_dir  # text file output directory
    cfg_opts = _utils.get_output_config('dna_substitutions')  # get config

    # handle all DNA substitutions
    nuc_ctr = count_nuc_substitutions(conn)  # get substitution counts
    count_save_path = out_dir + cfg_opts['sub_out']
    save_nuc_substitutions(nuc_ctr, count_save_path)

    # count oncogene substitutions
    onco_nuc_ctr = count_oncogene_substitutions(conn)
    onco_save_path = out_dir + cfg_opts['onco_out']
    save_nuc_substitutions(onco_nuc_ctr, onco_save_path)

    # count oncogene substitutions
    tsg_nuc_ctr = count_tsg_substitutions(conn)
    tsg_save_path = out_dir + cfg_opts['tsg_out']
    save_nuc_substitutions(tsg_nuc_ctr, tsg_save_path)

    # plot results
    all_heatmap = _utils.plot_dir + cfg_opts['sub_heatmap']
    plot_data.nuc_substitution_heatmap(count_save_path, all_heatmap)
    all_barplot = _utils.plot_dir + cfg_opts['sub_barplot']
    plot_data.nuc_substitution_barplot(count_save_path, all_barplot,
                                       title='All DNA Substitution Mutations')
    onco_barplot = _utils.plot_dir + cfg_opts['onco_barplot']
    plot_data.nuc_substitution_barplot(onco_save_path, onco_barplot,
                                       title='Oncogene DNA Substitution Mutations')
    tsg_barplot = _utils.plot_dir + cfg_opts['tsg_barplot']
    plot_data.nuc_substitution_barplot(tsg_save_path, tsg_barplot,
                                       title='TSG DNA Substitution Mutations')
