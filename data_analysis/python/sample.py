"""
The sample module analyzes characteristics of samples in
the COSMIC database by using the cosmic_mutation.
"""

import utils.python.util as _utils
import pandas.io.sql as psql
import plot_data


def count_mutated_genes(conn):
    """Count the number of genes that are mutated in each sample.

    **Parameters**

    conn : MySQLdb connection
        connection to COSMIC_nuc

    **Returns**

    df : pd.DataFrame
        two column data frame of sample names and gene cts
    """
    sql = ('SELECT x.SampleName, SUM(x.gene_indicator) as GeneCounts'
          ' FROM ('
          '     SELECT SampleName, Gene, 1 as gene_indicator'
          '     FROM cosmic_mutation'
          '     GROUP BY SampleName, Gene'
          ' ) x GROUP BY SampleName'
          ' ORDER BY GeneCounts Desc;')
    df = psql.frame_query(sql, con=conn)
    df.GeneCounts = df.GeneCounts.astype(int)  # pandas is not auto detecting
    return df


def count_mutations(conn):
    """Count the number of mutations in each sample.

    **Parameters**

    conn : MySQLdb connection
        connection to COSMIC_nuc

    **Returns**

    df : pd.DataFrame
        two column data frame of sample names and mutation cts
    """
    sql = ('SELECT x.SampleName, SUM(x.mut_indicator) as MutationCounts'
          ' FROM ('
          '     SELECT SampleName, 1 as mut_indicator'
          '     FROM cosmic_mutation'
          ' ) x GROUP BY SampleName'
          ' ORDER BY MutationCounts Desc;')
    df = psql.frame_query(sql, con=conn)
    df.MutationCounts = df.MutationCounts.astype(int)  # pandas is not auto detecting
    return df


def main(conn):
    out_dir = _utils.result_dir  # output directory for text files
    cfg_opts = _utils.get_output_config('sample')

    # get info about sample names
    sample_gene_cts = count_mutated_genes(conn)
    sample_mutation_cts = count_mutations(conn)
    sample_gene_cts.to_csv(out_dir + cfg_opts['gene_out'],
                           sep='\t', index=False)
    sample_mutation_cts.to_csv(out_dir + cfg_opts['mutation_out'],
                               sep='\t', index=False)

    # plot results
    plot_data.sample_barplot(sample_mutation_cts,
                             out_dir + cfg_opts['sample_barplot'],
                             title='Composition of database from size of samples',
                             xlabel='Size of sample',
                             ylabel='Number of Mutations')
