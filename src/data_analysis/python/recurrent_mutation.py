"""
This module handles counting recurrent mutations.
"""

from src.utils.python.amino_acid import AminoAcid
import src.utils.python.util as _utils
import pandas as pd
import pandas.io.sql as psql
import plot_data
import numpy as np
from collections import Counter
import logging

logger = logging.getLogger(__name__)

def _count_recurrent_missense(hgvs_iterable,
                              bin_width=1):
    """Counts the total missense mutations and stratifies missense
    counts according to position.

    Parameters
    ----------
    hgvs_iterable : iterable
        container object with HGVS protein strings
    bin_width : int
        size (length) to still consider mutation in same bin.
        Default is 1, meaning mutations have to occur at the
        exact same position.

    Returns
    -------
    gene_pos_ctr : dict
        counts missense mutations by position. {mutation_position: count, ...}
    total_missense_ctr : int
        total missense mutation count

    NOTE: This function requires the input HGVS container to contain
    mutations only for a single gene.
    """
    gene_pos_ctr = {}
    total_missense_ctr = 0
    for hgvs in hgvs_iterable:
        aa = AminoAcid(hgvs)
        if aa.mutation_type == 'missense':
            gene_pos_ctr.setdefault(aa.pos, 0)
            gene_pos_ctr[aa.pos] += 1  # add 1 to dict of pos
            total_missense_ctr += 1  # add 1 to total missense

    # now aggregate mutations using bin width
    unique_bins = set([int(key/bin_width)
                       for key in gene_pos_ctr.keys()])
    gene_pos_bin_ctr = {mybin: 0 for mybin in unique_bins}
    for p in gene_pos_ctr.keys():
        tmp_bin = int(p/bin_width)
        gene_pos_bin_ctr[tmp_bin] += gene_pos_ctr[p]

    # return gene_pos_ctr, total_missense_ctr
    return gene_pos_bin_ctr, total_missense_ctr


def count_missense_types(hgvs_iterable,
                         recurrent_min=2,
                         recurrent_max=float('inf'),
                         bin_width=1):
    """Count the number of recurrent missense and regular missense
    mutations given a valid recurrency range.

    Parameters
    ----------
    hgvs_iterable : iterable
        contains HGVS protein strings for a single gene.
    recurrent_min : int, (default=2)
        minimum number of missense mutations to define a recurrent position.
    recurrent_max : int, (default=infinity)
        maximum number of missense mutations to define a recurrent position.
        This is just a purely debugging parameter.

    Returns
    -------
    recurrent_cts : int
        number of recurrent missense mutations in a gene
    missense_cts : int
        number of regular missense mutations in a gene
    """
    pos_count, total_missense = _count_recurrent_missense(hgvs_iterable, bin_width)
    recurrent_cts = sum([cts for cts in pos_count.values()  # sum mutations passing theshold
                         if cts >= recurrent_min and cts <= recurrent_max])
    missense_cts = total_missense - recurrent_cts  # subtract reccurent from total
    return recurrent_cts, missense_cts


def count_recurrent_by_number(hgvs_iterable):
    """Count the number of recurrent positions with a given number recurrent
    counts."""
    pos_count, total_missense = _count_recurrent_missense(hgvs_iterable)
    mycounter = Counter(pos_count.values())
    return mycounter


def unique_missense_positions(conn):
    """This function appears to not be used. Possible bug in calling
    count_recurrent_by_number function."""
    logger.info('Calculating unique missense positions. . .')

    # query database
    sql = "SELECT Gene, AminoAcid FROM mutation"  # get everything from table
    df = psql.frame_query(sql, con=conn)
    gene_to_indexes = df.groupby('Gene').groups

    # calculate missense position entropy by gene
    gene_items = gene_to_indexes.items()
    gene_list, _ = zip(*gene_items)
    result_df = pd.DataFrame(np.zeros(len(gene_list)), columns=['unique missense position'], index=gene_list)
    for gene, indexes in gene_items:
        tmp_df = df.ix[indexes]
        myct = count_recurrent_by_number(tmp_df['AminoAcid'])
        pos_ct = np.array(myct.values())  # convert to numpy array
        unique_missense = len(pos_ct)
        result_df.ix[gene, 'unique missense position'] = unique_missense

    logger.info('Finsihed counting unique missense positions.')
    return result_df


def count_recurrent(conn):
    """Counts the number of recurrent mutations by how many recurrent
    mutations were at a particular position.

    All counts are stratified by gene type (oncogene, TSG, other).

    Parameters
    ----------
    conn : db connection
        connection to 20/20+ database
    """
    logger.info('Counting the number of recurrent positions . . .')

    # query database
    sql = "SELECT Gene, Protein_Change as AminoAcid FROM mutation"  # get everything from table
    df = psql.frame_query(sql, con=conn)
    gene_to_indexes = df.groupby('Gene').groups

    # counters
    onco_ct = Counter()  # counter for oncogenes
    tsg_ct = Counter()  # counter for tumor suppressor genes
    all_ct = Counter()  # counter for all genes
    other_ct = Counter()  # counter for gene other than tsg/onco

    # count recurrent missense
    for gene, indexes in gene_to_indexes.iteritems():
        tmp_df = df.ix[indexes]
        myct = count_recurrent_by_number(tmp_df['AminoAcid'])
        all_ct += myct
        if gene in _utils.oncogene_set:
            onco_ct += myct
        elif gene in _utils.tsg_set:
            tsg_ct += myct
        else:
            other_ct += myct

    logger.info('Finished counting recurent.')
    return all_ct, onco_ct, tsg_ct, other_ct


def main(conn):
    cfg_opts = _utils.get_output_config('recurrent_mutation')

    # get count info
    all_counts, onco_counts, tsg_counts, other_counts = count_recurrent(conn)

    # extract count data
    all_indx, all_data = zip(*list(all_counts.iteritems()))
    onco_indx, onco_data = zip(*list(onco_counts.iteritems()))
    tsg_indx, tsg_data = zip(*list(tsg_counts.iteritems()))
    other_indx, other_data = zip(*list(other_counts.iteritems()))

    # create temporary data frames
    all_df = pd.DataFrame({'all': all_data},
                          index=all_indx)
    onco_df = pd.DataFrame({'oncogene': onco_data},
                           index=onco_indx)
    tsg_df = pd.DataFrame({'tsg': tsg_data},
                          index=tsg_indx)
    other_df = pd.DataFrame({'other': other_data},
                            index=other_indx)
    all_df.index.name = 'Number of Mutations'
    other_df.index.name = 'Number of Mutations'
    onco_df.index.name = 'Number of Mutations'
    tsg_df.index.name = 'Number of Mutations'

    # write output to file
    all_df.to_csv(_utils.result_dir + cfg_opts['all_recurrent'], sep='\t')
    onco_df.to_csv(_utils.result_dir + cfg_opts['onco_recurrent'], sep='\t')
    tsg_df.to_csv(_utils.result_dir + cfg_opts['tsg_recurrent'], sep='\t')
    other_df.to_csv(_utils.result_dir + cfg_opts['other_recurrent'], sep='\t')

    # add data to single dataframe for plotting purposes
    other_df['oncogene'] = onco_df['oncogene']
    other_df['tsg'] = tsg_df['tsg']

    # plot recurrent missense positions per gene
    tmp_path = _utils.plot_dir + cfg_opts['recur_missense_line']
    plot_data.recurrent_missense_pos_line(other_df,
                                          tmp_path)
