import utils.python
from utils.python.cosmic_db import get_cosmic_db
import dna_substitutions
import utils.python.util as _utils
import plot_data
import mutation_types
import single_gene
import recurrent_mutation
import position_entropy
import feature_matrix
import tables.sample_name as sample_name
import tables.cosmic_aa as cosmic_aa
import tables.cosmic_genomic as cosmic_genomic
import missense
import pandas.io.sql as psql
import sqlite3
import logging

logger = logging.getLogger(__name__)

def count_gene_mutations(conn):
    """Count the total number of mutations for each gene.

    **Parameters**

    conn : mysql/sqlite connection
        connection to db with nucleotide table

    **Returns**

    df : pd.DataFrame
        two column data frame with gene and counts
    """
    logger.info('Counting number of mutations for each gene . . .')

    sql = ('SELECT Gene, COUNT(*) as count'
          ' FROM cosmic_mutation GROUP BY Gene'
          ' ORDER BY count DESC;')
    logger.debug('Gene mutation count SQL statement: ' + sql)

    df = psql.frame_query(sql, con=conn)  # execute query
    logger.info('Finished getting gene mutation counts.')
    return df


def main(recurrent, recurrent_max, db, classify_only):
    cfg_opts = _utils.get_output_config('stats')  # get config

    if db == 'cosmic_nuc':
        conn = get_cosmic_db()  # connect to COSMIC_nuc

        # skip re-runing entire data_analysis if user specified
        if not classify_only:
            missense.main(conn, 'cosmic_aa')  # plot missense info

            # check info about COSMIC_nuc tables
            cosmic_aa.main()  # check cosmic_aa table
            cosmic_genomic.main()  # check cosmic_genomic table
    else:
        # connect to sqlite db at data/genes.db
        genes_db_path = _utils.get_db_config('champ')['db']
        conn = sqlite3.connect(genes_db_path)

        if not classify_only:
            missense.main(conn, 'cosmic_mutation')  # specify different table name

    feature_matrix.main(recurrent, recurrent_max, conn)  # generate mutation count feature matrix
    position_entropy.main(conn)

    # user can specify a flag to prevent complete updates of the
    # data_analysis results. if classify_only is specified only
    # the pertinent information for classification is updated.
    if not classify_only:

        # look at individual genes
        with open(_utils.config_dir + 'single_gene.txt') as handle:
            for row in handle:
                gene = row.strip()
                try:
                    single_gene.main(gene, recurrent, conn)
                except:
                    raise
                    # be careful this catches all exceptions which might
                    # hide other errors
                    logger.debug('(Problem) Gene not found: %s' % gene)

        recurrent_mutation.main(conn)

        sample_name.main(conn)  # info related to mutations for each sample

        # handle DNA substitutions
        dna_substitutions.main(conn)

        # handle stats related to mutation types
        mutation_types.main(conn)

        # gene mutation counts
        gene_ct_df = count_gene_mutations(conn)
        gene_ct_df.set_index('Gene').to_csv(_utils.result_dir + cfg_opts['gene_mutation_counts'],
                                            sep='\t')
        plot_data.gene_mutation_histogram(gene_ct_df['count'],
                                        _utils.plot_dir + cfg_opts['gene_mutation_histogram'])
        plot_data.cumulative_gene_mutation(gene_ct_df['count'],
                                        _utils.plot_dir + cfg_opts['cumulative_gene_mutation'])

    conn.close()  # close data/genes.db connection
