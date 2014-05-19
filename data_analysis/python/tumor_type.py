"""Calculates stats related to tumor types."""
import utils.python.util as _utils
from utils.python.amino_acid import AminoAcid
import pandas.io.sql as psql
import pandas as pd
import plot_data
import utils.python.math as mymath
import numpy as np
import logging
import IPython

logger = logging.getLogger(__name__)

def filter_silent(df):
    """Filter out mutations that are silent.

    Parameters
    ----------
    df : pd.DataFrame
        dataframe from get_all_mutations function

    Results
    -------
    df : pd.DataFrame
        contains no rows with silent mutations
    """
    non_silent = map(lambda x: AminoAcid(x).is_non_silent,
                     df['AminoAcid'])
    df['aa_non_silent'] = non_silent
    is_splice_site = df['Variant_Classification'] == 'Splice_Site'
    df['non_silent'] = (df['aa_non_silent'] | is_splice_site).astype(int)
    df = df[df['non_silent']==1]  # only keep non-silent mutations
    return df


def ct_ns_tumor_types(df):
    """Count number of non-silent mutations per tumor type.

    Parameters
    ----------
    df : pd.DataFrame
        dataframe from get_all_mutations function

    Results
    -------
    non_silent_type_df : pd.DataFrame
        Dataframe aggregating number of non-silent mutations
        per tumor type
    """
    logger.info('Calculating non-silent mutations for tumor types . . .')
    # identify non-silent mutations
    df = filter_silent(df)

    # sum non-silent mutations for each tumor type
    non_silent_type_df = df[['Tumor_Type', 'non_silent']].groupby('Tumor_Type').sum()
    logger.info('Finished calculating non-silent mutations for tumor types.')
    return non_silent_type_df


def get_all_mutations(conn):
    """Get all mutations from the mutation table.

    Parameters
    ----------
    conn : sqlite or mysql db connection
        database contains mutation table of all mutations

    Returns
    -------
    df : pd.DataFrame
        Query results of selected cols from mutation table
    """
    sql = ('SELECT Gene, Tumor_Type, Tumor_Sample, Protein_Change as AminoAcid, '
           '       DNA_Change as Nucleotide, Variant_Classification '
           'FROM mutation')
    df = psql.frame_query(sql, con=conn)
    return df


def gene_ns_tumor_types(df):
    """Stratifies pct of non-silent mutations for each tumor type for genes.

    Parameters
    ----------
    df : pd.DataFrame
        df with silent mutations already filtered.

    Returns
    -------
    gene_mut_table : pd.DataFrame
        pivoted table with rows as genes and columns as tumor types.
        Each element is the percent of non-silent mutations for a given
        gene that are in a specific tumor type.
    """
    gene_mut_table = pd.pivot_table(df,
                                    values='non_silent',
                                    cols='Tumor_Type',
                                    rows='Gene',
                                    aggfunc=np.sum)
    gene_mut_table = gene_mut_table.fillna(0)  # NAs represent 0 counts
    gene_mut_table = gene_mut_table + .01  # add a pseudo count
    gene_mut_table = gene_mut_table.div(gene_mut_table.sum(axis=1).astype(float),
                                        axis=0)  # row normalize to 1 (i.e. a probability)
    return gene_mut_table


def js_dist_ttype(df, ttype_df):
    """Get the Jensen-Shannon distance of each genes tumor type
    distribution from the aggregate database distribution.

    Parameters
    ----------
    df : pd.DataFrame
        df from mutation table except with columns name `AminoAcid`
        and `Nucleotide`
    ttype_df : pd.Series
        indices are tumor types, values are # of non-silent mutations

    Returns
    -------
    js_df : pd.DataFrame
        df with 'JS distance' containing Jensen-Shannon distance and
        'true class' based on training labels
    """
    logger.info('Calculating JS distance of tumor type distribution . . .')
    pct_ttype_df = ttype_df / float(ttype_df.sum())
    non_silent_df = filter_silent(df)
    gene_mut_table = gene_ns_tumor_types(non_silent_df)
    # pct_ttype_df['non silent'] = 1/float(len(pct_ttype_df))
    # db_ttypes = pct_ttype_df['non silent'].reindex(gene_mut_table.columns)  # make sure indices align
    db_ttypes = pct_ttype_df.reindex(gene_mut_table.columns)
    db_ttypes = db_ttypes['non silent'].fillna(0.0)  # get series object
    js_df = pd.DataFrame({'JS distance': gene_mut_table.apply(mymath.js_distance,
                                                              args=(db_ttypes,),
                                                              axis=1)})
    js_df['true class'] = js_df.index.to_series().apply(_utils.classify_gene)
    logger.info('Finished calculating JS distance')
    return js_df


def main(conn):
    out_dir = _utils.result_dir  # output directory for text files
    plot_dir = _utils.plot_dir
    opts = _utils.get_output_config('tumor_type')

    # get mutations
    df = get_all_mutations(conn)  # get data from mutations table

    # get number of non-silent mutations for each tumor type
    ttype_df = ct_ns_tumor_types(df)
    ttype_df.to_csv(out_dir + opts['non_silent_count'], sep='\t')

    # plot barplot of non-silent mutations for each tumor type
    ttype_df.sort('non_silent', inplace=True)
    plot_data.non_silent_tumor_type_barplot(ttype_df,
                                            plot_dir + opts['non_silent_barplot'])

    # save non-silent mutations for each gene stratified by tumor type
    ns_df = filter_silent(df)
    ns_table = gene_ns_tumor_types(ns_df)
    ns_table.to_csv(out_dir + opts['gene_ns_ttype'], sep='\t')

    # plot js distance for tumor types
    js_df = js_dist_ttype(df, ttype_df)
    plot_data.js_distance_kde(js_df, plot_dir + opts['js_distance_kde'])
