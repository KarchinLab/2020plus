"""The features module generates a feature matrix each row is a gene
and features are specified in different columns. The generate_feature_matrix
function specifies how the feature matrix is constructed.
"""
from src.utils.python.amino_acid import AminoAcid
import recurrent_mutation as recur
import plot_data
import pandas.io.sql as psql
import src.utils.python.util as _utils
from collections import OrderedDict
import shutil
import csv
import logging

logger = logging.getLogger(__name__)

def generate_feature_matrix(df, recurrency_threshold,
                            recurrency_cap=float('inf')):
    """Generate a feature matrix potentially useful for classifying genes.

    Parameters
    ----------
    df : pd.DataFrame
        data frame consisting of cosmic mutations from db
    recurrency_threshold : int
        minimum number of mutations to define recurrency
    recurrency_cap : int
        maximum number of mutations to define recurrency.
        This is a debug parameter, so generally leave it
        at default setting of infinity.

    Returns
    -------
    features : list
        feature matrix
    """
    logger.info('Creating design matrix . . .')

    # query database
    mtypes = _utils.get_mutation_types(df['AminoAcid'],
                                       df['Nucleotide'],
                                       known_type=df['Variant_Classification'])
    mtypes.index = df.index  # make sure indexes match the data frame
    df['mut_types'] = mtypes  # add mutation types to SQL output
    # gene_to_indexes = df.groupby('Gene').groups
    gene_to_indexes = df.groupby('Gene').groups

    # aggregate info
    design_matrix = []
    not_used_types = ['not valid',
                      'missing',
                      'unknown effect',
                      'no protein']  # only include known mutations
    for gene, indexes in gene_to_indexes.iteritems():
        tmp_df = df.ix[indexes]
        gene_pos_counter = {}
        identical_indel = {}
        mut_type_ctr = OrderedDict([['Missense_Mutation', 0],
                                    ['Frame_Shift_Indel', 0],
                                    ['Silent', 0],
                                    #['not valid', 0],
                                    ['In_Frame_Indel', 0],
                                    #['no protein', 0],
                                    ['Nonstop_Mutation', 0],
                                    ['Translation_Start_Site', 0],
                                    ['Splice_Site', 0],
                                    #['missing', 0],
                                    ['Nonsense_Mutation', 0]])
                                    #['unknown effect', 0]])
        # count identical indels
        #for i, hgvs in enumerate(tmp_df['AminoAcid']):
            #aa = AminoAcid(hgvs)
            #if aa.mutation_type not in not_used_types:
                # do not use 'missing', 'unkown effect' or 'not valid'
                # if aa.mutation_type == 'missense':
                    # keep track of missense pos for recurrency
                #    gene_pos_counter.setdefault(aa.pos, 0)
                #    gene_pos_counter[aa.pos] += 1
                #if aa.mutation_type == 'In_Frame_Indel' and tmp_df['mut_types'].iloc[i] != 'Splice_Site':
                    # keep track of missense pos for recurrency
                    #identical_indel.setdefault(aa.hgvs_original, 0)
                    #identical_indel[aa.hgvs_original] += 1

        # count mutation types
        for mt in tmp_df['mut_types']:
            if mt not in not_used_types:
                mut_type_ctr[mt] += 1

        # combine lost start and lost stop
        lost_stop_ct = mut_type_ctr.pop('Nonstop_Mutation', 0)
        lost_start_ct = mut_type_ctr.pop('Translation_Start_Site', 0)
        mut_type_ctr['Nonstop_Mutation+Translation_Start_Site'] = lost_start_ct + lost_stop_ct

        recur_ct, missense_ct = recur.count_missense_types(tmp_df['AminoAcid'],
                                                           recurrency_threshold,
                                                           recurrency_cap)

        # needs to have at least one count
        if sum(mut_type_ctr.values()):
            #mut_type_ctr['missense'] = missense_ct
            #recurrent_cts = sum([cts for cts in gene_pos_counter.values()
            #                     if cts >= recurrency_threshold])
            identical_cts = sum([cts for cts in identical_indel.values()
                                 if cts >= recurrency_threshold])
            #mut_type_ctr['missense'] -= recurrent_cts  # subtract off the recurrent missense
            # mut_type_ctr['In_Frame_Indel'] -= identical_cts  # subtract off the recurrent missense
            design_matrix.append([gene, recur_ct] + list(mut_type_ctr.values()))
    header = [['gene', 'recurrent missense'] + list(mut_type_ctr)]
    logger.info('Finished creating feature matrix.')
    return header + design_matrix


def main(recurrent, recurrent_max, conn):
    cfg_opts = _utils.get_output_config('feature_matrix')  # get config

    # query db for all mutations
    sql = ("SELECT Gene, DNA_Change as Nucleotide, "
           "       Protein_Change as AminoAcid, Variant_Classification "
           "FROM mutation")
    df = psql.frame_query(sql, con=conn)  # get all

    # generate features
    feature_matrix = generate_feature_matrix(df, recurrent,
                                             recurrency_cap=recurrent_max)

    # save features
    feature_path = _utils.result_dir + cfg_opts['gene_feature_matrix']
    with open(feature_path, 'wb') as handle:
        csv.writer(handle, delimiter='\t').writerows(feature_matrix)
    copy_path = feature_path.strip('txt') + 'r%d.txt' % recurrent
    shutil.copy(feature_path, copy_path)  # record a second file with reccurent param in name

    # PCA plots
    # unnormalized PCA plot
    plot_data.pca_plot(_utils.result_dir + cfg_opts['gene_feature_matrix'],
                       _utils.plot_dir + cfg_opts['pca_plot'],
                       title='Protein Mutation Type Composition PCA')
    # normalized PCA by removing class imbalance
    plot_data.pca_plot(_utils.result_dir + cfg_opts['gene_feature_matrix'],
                       _utils.plot_dir + cfg_opts['pca_plot_rand_3'],
                       norm_class=3,
                       low_count_filter=10,
                       title='Protein Mutation Type Composition PCA Subsampled by '
                       'Gene Type (3:1)')
    plot_data.pca_plot(_utils.result_dir + cfg_opts['gene_feature_matrix'],
                       _utils.plot_dir + cfg_opts['pca_plot_rand_1'],
                       norm_class=1,
                       low_count_filter=10,
                       title='Protein Mutation Type Composition PCA Subsampled by '
                       'Gene Type (1:1)')
