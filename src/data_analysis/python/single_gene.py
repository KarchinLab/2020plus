import src.utils.python.util as _utils
import recurrent_mutation as recur
import pandas as pd
import pandas.io.sql as psql
import logging
import src.utils.python.plot as myplt
import os

logger = logging.getLogger(__name__)

def count_primary_tissues(gene, conn):
    """Count the number of mutations in a single gene for each primary tissue.

    **Parameters**

    gene : str
        gene name
    conn : sqlite/MySQL connection
        connection to database with `nucleotide` table
    """
    logger.info('Counting %s mutations for each primary tissue ...' % gene)
    sql = ("SELECT PrimaryTissue, COUNT(*) as Counts"
           " FROM cosmic_mutation cm"
           " WHERE cm.Gene='%s'"
           " GROUP BY PrimaryTissue"
           " ORDER BY Counts Desc;" % gene)
    df = psql.frame_query(sql, con=conn)
    df.PrimaryTissue = df.PrimaryTissue.str.replace('_', ' ')
    df = df.set_index('PrimaryTissue')
    logger.info('Finished counting %s mutations.' % gene)
    return df


def count_types_primary_tissue(gene, recurrent_min, conn):
    """Count mutation types in a single gene for each primary tissue.

    **Parameters**

    gene : str
        gene name
    conn : sqlite/MySQL connection
        connection to database with `nucleotide` table
    """
    sql = ("SELECT PrimaryTissue, Nucleotide, AminoAcid "
           "FROM cosmic_mutation cm "
           "WHERE cm.Gene='%s' "
           "ORDER BY PrimaryTissue" % gene)
    df = psql.frame_query(sql, con=conn)

    sql = ("SELECT PrimaryTissue, COUNT(*) as Counts "
           "FROM cosmic_mutation cm "
           "WHERE cm.Gene='%s' "
           "GROUP BY PrimaryTissue "
           "ORDER BY Counts Desc" % gene)
    sorted_df = psql.frame_query(sql, con=conn)
    sorted_df['PrimaryTissue'] = sorted_df['PrimaryTissue'].apply(lambda x: x.replace('_', ' '))

    groups = df.groupby('PrimaryTissue').groups
    aa, nuc = {}, {}
    for k, v in groups.iteritems():
        aa_types = _utils.count_mutation_types(df.ix[v]['AminoAcid'],
                                               df.ix[v]['Nucleotide'])
        recur_ct, missense_ct = recur.count_missense_types(df.ix[v]['AminoAcid'],
                                                           recurrent_min=recurrent_min)
        aa_types = aa_types.set_value('missense', missense_ct)
        aa_types = aa_types.set_value('recurrent missense', recur_ct)
        nuc_types = _utils.count_mutation_types(df.ix[v]['Nucleotide'],
                                                kind='nucleotide')
        k = k.replace('_', ' ')
        aa[k] = aa_types
        nuc[k] = nuc_types
    aa_df = pd.DataFrame(aa).fillna(0).astype(int).T.ix[sorted_df['PrimaryTissue']]
    nuc_df = pd.DataFrame(nuc).fillna(0).astype(int).T.ix[sorted_df['PrimaryTissue']]
    return aa_df, nuc_df


def main(gene_name, recurrent_count, conn):
    cfg_opts = _utils.get_output_config('single_gene')

    # set up directory
    gene_result_dir = _utils.result_dir + "single_gene/" + gene_name + "/"
    gene_plot_dir = _utils.plot_dir + "single_gene/" + gene_name + "/"
    if not os.path.isdir(gene_result_dir):
        os.mkdir(gene_result_dir)
    if not os.path.isdir(gene_plot_dir):
        os.mkdir(gene_plot_dir)

    tissue_cts = count_primary_tissues(gene_name, conn)
    tissue_cts.to_csv(gene_result_dir + cfg_opts['primary_tissue'],
                      sep='\t')
    myplt.barplot(tissue_cts,
                  gene_plot_dir + cfg_opts['primary_tissue_barplot'],
                  ylabel='Counts',
                  title='%s Mutation Counts in Primary Tissues' % gene_name.replace('_',' '))

    # stratify by mutation types
    aa_df, nuc_df = count_types_primary_tissue(gene_name, recurrent_count, conn)
    aa_df.to_csv(gene_result_dir + cfg_opts['aa_type_primary_tissue'],
                 sep='\t')
    nuc_df.to_csv(gene_result_dir + cfg_opts['nuc_type_primary_tissue'],
                  sep='\t')
    myplt.barplot(aa_df,
                  gene_plot_dir + cfg_opts['aa_type_primary_tissue_barplot'],
                  stacked=True,
                  ylabel='Counts',
                  title='%s Mutation Counts in Primary Tissues' % gene_name.replace('_', ' '))
    myplt.barplot(nuc_df,
                  gene_plot_dir + cfg_opts['nuc_type_primary_tissue_barplot'],
                  stacked=True,
                  ylabel='Counts',
                  title='%s Mutation Counts in Primary Tissues' % gene_name.replace('_', ' '))
