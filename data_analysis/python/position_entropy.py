import utils.python.util as _utils
import recurrent_mutation
import utils.python.math as mymath
from utils.python.amino_acid import AminoAcid
import plot_data
import pandas as pd
import numpy as np
import pandas.io.sql as psql
import logging

logger = logging.getLogger(name=__file__)

def _count_mutation_position(hgvs_iterable):
    """Counts the total mutations and stratifies mutation
    counts according to position.

    **Parameters**

    hgvs_iterable : iterable
        container object with HGVS protein strings

    **Returns**

    gene_pos_ctr : dict
        maps positions to mutation counts
    total_mutation_ctr : int
        total number of mutations

    NOTE: This function requires the input HGVS container to contain
    mutations only for a single gene.
    """
    gene_pos_ctr = {}
    total_mutation_ctr = 0
    for hgvs in hgvs_iterable:
        aa = AminoAcid(hgvs)
        if aa.is_valid and aa.pos:
            gene_pos_ctr.setdefault(aa.pos, 0)
            gene_pos_ctr[aa.pos] += 1  # add 1 to dict of pos
            total_mutation_ctr += 1  # add 1 to total missense
    return gene_pos_ctr, total_mutation_ctr


def query_amino_acid(conn):
    sql = "SELECT Gene, Protein_Change as AminoAcid FROM mutation"  # get everything from table
    df = psql.frame_query(sql, con=conn)
    return df


def mutation_position_entropy(df):
    logger.info('Calculating mutation position entropy . . .')

    gene_to_indexes = df.groupby('Gene').groups

    # calculate missense position entropy by gene
    gene_items = gene_to_indexes.items()
    gene_list, _ = zip(*gene_items)
    result_df = pd.DataFrame(np.zeros((len(gene_list), 2)),
                             columns=['mutation position entropy',
                                      'pct of uniform mutation entropy'],
                             index=gene_list)
    for gene, indexes in gene_items:
        tmp_df = df.ix[indexes]
        myct, total_ct = _count_mutation_position(tmp_df['AminoAcid'])
        if total_ct > 0:
            pos_ct = np.array(myct.values())  # convert to numpy array
            p = pos_ct / float(total_ct)  # normalize to a probability
            mutation_entropy = mymath.shannon_entropy(p)  # calc shannon entropy
            if total_ct > 1:
                max_ent = mymath.max_shannon_entropy(total_ct) if total_ct > 1 else 1
                percent_mutation_entropy = mutation_entropy / max_ent  # percent of max entropy
            else:
                percent_mutation_entropy = 0
        else:
            mutation_entropy = 0
            percent_mutation_entropy = 0
        result_df.ix[gene, 'mutation position entropy'] = mutation_entropy  # store result
        result_df.ix[gene, 'pct of uniform mutation entropy'] = percent_mutation_entropy  # store result
    logger.info('Finsihed calculating mutation position entropy.')
    return result_df


def missense_position_entropy(df):
    logger.info('Calculating missense position entropy . . .')

    gene_to_indexes = df.groupby('Gene').groups

    # calculate missense position entropy by gene
    gene_items = gene_to_indexes.items()
    gene_list, _ = zip(*gene_items)
    result_df = pd.DataFrame(np.zeros((len(gene_list), 2)),
                             columns=['missense position entropy',
                                      'pct of uniform missense entropy'],
                             index=gene_list)
    gene_len_df = _utils.get_gene_length()['gene_length']
    gene_len_mean = gene_len_df.mean()
    known_genes = set(gene_len_df.index.tolist())
    mybin_width = 1
    for gene, indexes in gene_items:
        tmp_df = df.ix[indexes]
        # myct = recurrent_mutation.count_recurrent_by_number(tmp_df['AminoAcid'])
        myct, total_missense = recurrent_mutation._count_recurrent_missense(tmp_df['AminoAcid'],
                                                                            bin_width=mybin_width)
        if total_missense > 0:

            # process positional info
            pos_ct = np.array(myct.values())  # convert to numpy array
            total_missense = np.sum(pos_ct)  # total number of missense
            p = pos_ct / float(total_missense)  # normalize to a probability
            missense_entropy = mymath.shannon_entropy(p)  # calc shannon entropy

            # get gene length in case number of mutations is greater than length
            tmp_gene_len = gene_len_df.ix[gene] if gene in known_genes else gene_len_mean
            # in sum cases gene length may be smaller than what a mutation
            # occurs at. In this case, extend gene_length to allow for mutation.
            len_from_bins = (max(myct.keys()) + 1) * mybin_width
            tmp_gene_len = len_from_bins if len_from_bins > tmp_gene_len else tmp_gene_len

            # calculate pct of maximum entropy
            if (total_missense > 1 and len(pos_ct) > 1) and tmp_gene_len > mybin_width:
                tmp_num_bins = int(tmp_gene_len/mybin_width)
                max_bins = tmp_num_bins if not tmp_gene_len % mybin_width else tmp_num_bins + 1

                # find the effective max number of counts at unique positions/bins
                max_count = total_missense if total_missense < max_bins else max_bins

                # calc max entropy
                # max_ent = mymath.max_shannon_entropy(total_missense)
                max_ent = mymath.max_shannon_entropy(max_count)
                percent_missense_entropy = missense_entropy / max_ent
            else:
                percent_missense_entropy = 0
        else:
            missense_entropy = 0
            percent_missense_entropy = 0
        result_df.ix[gene, 'missense position entropy'] = missense_entropy  # store result
        result_df.ix[gene, 'pct of uniform missense entropy'] = percent_missense_entropy  # store result

    logger.info('Finsihed calculating missense position entropy.')
    return result_df


def main(conn):
    cfg_opts = _utils.get_output_config('position_entropy')

    # query for amino acid mutations
    aa_df = query_amino_acid(conn)

    # get information about missense position entropy
    entropy_df = missense_position_entropy(aa_df)
    entropy_df['true class'] = 0
    onco_mask = [True if gene in _utils.oncogene_set else False for gene in entropy_df.index]
    tsg_mask = [True if gene in _utils.tsg_set else False for gene in entropy_df.index]
    entropy_df.ix[onco_mask, 'true class'] = 1
    entropy_df.ix[tsg_mask, 'true class'] = 2
    entropy_df.to_csv(_utils.result_dir + cfg_opts['missense_pos_entropy'], sep='\t')

    # plot distribution of missense position entropy
    plot_data.entropy_kde(entropy_df,
                          'missense position entropy',  # use the missense entropy column
                          _utils.plot_dir + cfg_opts['missense_pos_entropy_dist'],
                          title='Distribution of Missense Position Entropy',
                          xlabel='Missense Position Entropy (bits)')
    plot_data.entropy_kde(entropy_df,
                          'pct of uniform missense entropy',  # use the missense entropy column
                          _utils.plot_dir + cfg_opts['pct_missense_pos_entropy_dist'],
                          title='Distribution of Percentage of Maximum Missense Position Entropy',
                          xlabel='Percentage of Maximum Missense Position Entropy (bits)')


    # get information about mutation position entropy
    entropy_df = mutation_position_entropy(aa_df)
    entropy_df['true class'] = 0
    onco_mask = [True if gene in _utils.oncogene_set else False for gene in entropy_df.index]
    tsg_mask = [True if gene in _utils.tsg_set else False for gene in entropy_df.index]
    entropy_df.ix[onco_mask, 'true class'] = 1
    entropy_df.ix[tsg_mask, 'true class'] = 2
    entropy_df.to_csv(_utils.result_dir + cfg_opts['mutation_pos_entropy'], sep='\t')

    # plot distribution of mutation position entropy
    plot_data.entropy_kde(entropy_df,
                          'mutation position entropy',  # use the mutation entropy column
                          _utils.plot_dir + cfg_opts['mutation_pos_entropy_dist'],
                          title='Distribution of Mutation Position Entropy',
                          xlabel='Mutation Position Entropy (bits)')
    plot_data.entropy_kde(entropy_df,
                          'pct of uniform mutation entropy',  # use the missense entropy column
                          _utils.plot_dir + cfg_opts['pct_mutation_pos_entropy_dist'],
                          title='Distribution of Percentage of Maximum Mutation Position Entropy',
                          xlabel='Percentage of Maximum Mutation Position Entropy (bits)')
