import utils.python.util as _utils
from utils.python.cosmic_db import get_cosmic_db
import pandas as pd
import pandas.io.sql as psql
import logging
import utils.python.plot as myplt

logger = logging.getLogger(__name__)

def count_primary_tissues(conn):
    """Count the number of mutations in TP53 for each primary tissue."""
    logger.info('Counting TP53 mutations for each primary tissue . . .')
    sql = ("SELECT PrimaryTissue, COUNT(*) as Counts"
           " FROM `nucleotide` nuc"
           " WHERE nuc.Gene='TP53'"
           " GROUP BY PrimaryTissue"
           " ORDER BY Counts Desc;")
    df = psql.frame_query(sql, con=conn)
    df.PrimaryTissue = df.PrimaryTissue.str.replace('_', ' ')
    df = df.set_index('PrimaryTissue')
    logger.info('Finished counting TP53 mutations.')
    return df


def count_types_primary_tissue(conn):
    sql = ("SELECT PrimaryTissue, Nucleotide, AminoAcid "
           "FROM `nucleotide` nuc "
           "WHERE nuc.Gene='TP53' "
           "ORDER BY PrimaryTissue")
    df = psql.frame_query(sql, con=conn)

    sql = ("SELECT PrimaryTissue, COUNT(*) as Counts "
           "FROM `nucleotide` nuc "
           "WHERE nuc.Gene='TP53' "
           "GROUP BY PrimaryTissue "
           "ORDER BY Counts Desc")
    sorted_df = psql.frame_query(sql, con=conn)
    sorted_df['PrimaryTissue'] = sorted_df['PrimaryTissue'].apply(lambda x: x.replace('_', ' '))

    groups = df.groupby('PrimaryTissue').groups
    aa, nuc = {}, {}
    for k, v in groups.iteritems():
        aa_types = _utils.count_mutation_types(df.ix[v]['AminoAcid'])
        nuc_types = _utils.count_mutation_types(df.ix[v]['Nucleotide'],
                                                kind='nucleotide')
        k = k.replace('_', ' ')
        aa[k] = aa_types
        nuc[k] = nuc_types
    aa_df = pd.DataFrame(aa).fillna(0).astype(int).T.ix[sorted_df['PrimaryTissue']]
    nuc_df = pd.DataFrame(nuc).fillna(0).astype(int).T.ix[sorted_df['PrimaryTissue']]
    return aa_df, nuc_df


def main():
    cfg_opts = _utils.get_output_config('tp53')
    result_dir = _utils.result_dir
    plot_dir = _utils.plot_dir

    conn = get_cosmic_db()  # open connection to COSMIC_nuc
    tissue_cts = count_primary_tissues(conn)
    tissue_cts.to_csv(result_dir + cfg_opts['primary_tissue'],
                      sep='\t')
    myplt.barplot(tissue_cts,
                  plot_dir + cfg_opts['primary_tissue_barplot'],
                  ylabel='Counts',
                  title='TP53 Mutation Counts in Primary Tissues')

    # stratify by mutation types
    aa_df, nuc_df = count_types_primary_tissue(conn)
    aa_df.to_csv(result_dir + cfg_opts['aa_type_primary_tissue'],
                 sep='\t')
    nuc_df.to_csv(result_dir + cfg_opts['nuc_type_primary_tissue'],
                  sep='\t')
    myplt.barplot(aa_df,
                  plot_dir + cfg_opts['aa_type_primary_tissue_barplot'],
                  stacked=True,
                  ylabel='Counts',
                  title='TP53 Mutation Counts in Primary Tissues')
    myplt.barplot(nuc_df,
                  plot_dir + cfg_opts['nuc_type_primary_tissue_barplot'],
                  stacked=True,
                  ylabel='Counts',
                  title='TP53 Mutation Counts in Primary Tissues')

    conn.close()
