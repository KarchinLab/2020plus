import utils.python.util as _utils
from utils.python.amino_acid import AminoAcid
import pandas.io.sql as psql
import pandas as pd
import plot_data
import utils.python.math as mymath
import numpy as np
import logging

logger = logging.getLogger(__name__)

def get_non_silent(df):
    non_silent = map(lambda x: AminoAcid(x).is_non_silent,
                     df['AminoAcid'])
    df['aa_non_silent'] = non_silent
    is_splice_site = df['Variant_Classification'] == 'Splice_Site'
    df['non_silent'] = (df['aa_non_silent'] | is_splice_site).astype(int)
    df = df[df['non_silent']==1]  # only keep non-silent mutations
    return df


def count_non_silent_tumor_types(df):
    # identify non-silent mutations
    df = get_non_silent(df)

    # sum non-silent mutations for each tumor type
    non_silent_type_df = df[['Tumor_Type', 'non_silent']].groupby('Tumor_Type').sum()
    return non_silent_type_df


def main(conn):
    out_dir = _utils.result_dir  # output directory for text files
    plot_dir = _utils.plot_dir
    opts = _utils.get_output_config('tumor_type')

    # get mutations
    sql = ('SELECT Gene, Tumor_Type, Tumor_Sample, Protein_Change as AminoAcid, '
           '       DNA_Change as Nucleotide, Variant_Classification '
           'FROM mutation')
    df = psql.frame_query(sql, con=conn)

    # get number of non-silent mutations for each tumor type
    ttype_df = count_non_silent_tumor_types(df)
    ttype_df.to_csv(out_dir + opts['non_silent_count'], sep='\t')

    # plot barplot of non-silent mutations for each tumor type
    ttype_df.rename(columns=lambda x: x.replace('_', ' '), inplace=True)
    ttype_df.rename(index=lambda x: x.replace('_', ' '), inplace=True)
    ttype_df.sort('non silent', inplace=True)
    plot_data.non_silent_tumor_type_barplot(ttype_df,
                                            plot_dir + opts['non_silent_barplot'])

    # plot js distance
    pct_ttype_df = ttype_df / float(ttype_df.sum())
    non_silent_df = get_non_silent(df)
    gene_mut_table = pd.pivot_table(non_silent_df,
                                    values='non_silent',
                                    cols='Tumor_Type',
                                    rows='Gene',
                                    aggfunc=np.sum)
    gene_mut_table = gene_mut_table.fillna(0)  # NAs represent 0 counts
    gene_mut_table = gene_mut_table + 1/5.  # add a pseudo count
    gene_mut_table = gene_mut_table.div(gene_mut_table.sum(axis=1).astype(float),
                                        axis=0)  # row normalize to 1 (i.e. a probability)
    pct_ttype_df['non silent'] = 1/float(len(pct_ttype_df))
    db_ttypes = pct_ttype_df['non silent'].reindex(gene_mut_table.columns)  # make sure indices align
    js_df = pd.DataFrame({'JS distance': gene_mut_table.apply(mymath.js_distance, args=(db_ttypes,), axis=1)})
    js_df['true class'] = js_df.index.to_series().apply(_utils.classify_gene)
    plot_data.js_distance_kde(js_df, plot_dir + opts['js_distance_kde'])
