import utils.python
from utils.python.cosmic_db import get_cosmic_db
from utils.python.amino_acid import AminoAcid
import utils.python.util as _utils
import plot_data
import pandas as pd
import pandas.io.sql as psql
import csv
from collections import OrderedDict
import logging


def count_mutations(cursor):
    """Count the number of entries"""
    cursor.execute("""SELECT COUNT(COSMICSampleID)
                   FROM `nucleotide`""")
    return cursor.fetchone()[0]  # COUNT query returns a tuple


def count_aa_mutation_types(cursor):
    """Count the amino acid mutation types (missense, indel, etc.).
    """
    df = psql.frame_query("""SELECT * FROM `nucleotide`""", con=cursor)
    unique_cts = _utils.count_mutation_types(df['AminoAcid'])
    return unique_cts


def generate_design_matrix(conn):
    """Generate a design matrix potentially useful for classifying genes.

    Args:
        conn (MySQLdb connection): database connection to mysql

    Returns:
        pd.DataFrame: Design matrix
    """
    logger = logging.getLogger(__name__)
    logger.info('Creating design matrix . . .')

    df = psql.frame_query("SELECT * FROM `nucleotide`", con=conn)  # get all
    mutation_types = _utils.get_mutation_types(df['AminoAcid'])
    df['mut_types'] = mutation_types  # add mutation types to SQL output
    gene_to_indexes = df.groupby('Gene').groups

    # aggregate info
    design_matrix = []
    for gene, indexes in gene_to_indexes.iteritems():
        tmp_df = df.ix[indexes]
        gene_pos_counter = {}
        mut_type_ctr = OrderedDict([['missense', 0],
                                    ['frame shift', 0],
                                    ['synonymous', 0],
                                    ['not valid', 0],
                                    ['indel', 0],
                                    ['missing', 0],
                                    ['nonsense', 0]])
        for hgvs in tmp_df['AminoAcid']:
            aa = AminoAcid(hgvs)
            if aa.is_missense:
                gene_pos_counter.setdefault(aa.pos, 0)
                gene_pos_counter[aa.pos] += 1
            mut_type_ctr[aa.mutation_type] += 1

        recurrent_cts = sum([cts for cts in gene_pos_counter.values() if cts > 1])
        mut_type_ctr['missense'] -= recurrent_cts  # subtract off the recurrent missense
        design_matrix.append([gene, recurrent_cts] + list(mut_type_ctr.values()))
    header = [['gene', 'recurrent missense'] + list(mut_type_ctr)]
    return header + design_matrix


def count_onco_aa_mut_types(conn):
    logger = logging.getLogger(__name__)
    logger.info('Counting oncogene mutation types . . .')

    # prepare sql statement
    oncogenes = _utils.read_oncogenes()
    sql = "SELECT * FROM `nucleotide` WHERE Gene in " + str(oncogenes)
    logger.debug('Oncogene SQL statement: ' + sql)

    df = psql.frame_query(sql, con=conn)  # execute query

    # count mutation types
    counts = _utils.count_mutation_types(df['AminoAcid'])
    logger.info('Finished counting oncogene mutation types.')
    return counts


def count_tsg_aa_mut_types(conn):
    logger = logging.getLogger(__name__)
    logger.info('Counting tumor suppressor gene mutation types . . .')

    # prepare sql statement
    tsgs = _utils.read_tsgs()
    sql = "SELECT * FROM `nucleotide` WHERE Gene in " + str(tsgs)
    logger.debug('Oncogene SQL statement: ' + sql)

    df = psql.frame_query(sql, con=conn)  # execute query

    # count mutation types
    counts = _utils.count_mutation_types(df['AminoAcid'])
    logger.info('Finished counting tumor suppressor gene mutation types.')
    return counts


def count_gene_mutations(conn):
    logger = logging.getLogger(__name__)
    logger.info('Counting number of mutations for each gene . . .')

    sql = """SELECT Gene, COUNT(*) as count
          FROM `nucleotide` GROUP BY Gene
          ORDER BY count DESC;"""
    logger.debug('Gene mutation count SQL statement: ' + sql)

    df = psql.frame_query(sql, con=conn)  # execute query
    logger.info('Finished getting gene mutation counts.')
    return df


def count_aa_missense_changes(cursor):
    """Count amino acid changes.

    Args:
        cursor: mysqldb cursor object

    Returns:
        dict. containing counts eg. {('aa1', 'aa2'): 4}
    """
    logger = logging.getLogger(name=__name__)
    logger.info('Starting to count amino acid changes . . .')
    cursor.execute("""SELECT aachange, occurrences
                   FROM `cosmic_aa`""")
    aa_change_counter = {}
    for aachange, occurrences in cursor.fetchall():
        aa = AminoAcid(hgvs=aachange,
                       occurrence=occurrences)
        if aa.is_valid and not aa.is_missing_info:
            aa_change_counter.setdefault((aa.initial, aa.mutated), 0)
            aa_change_counter[(aa.initial, aa.mutated)] += aa.occurrence
    logger.info('Finished counting amino acid changes.')
    return aa_change_counter


def save_aa_missense_counts(aacounter):
    """Saves missense mutation counts to file.

    """
    # save missense mutation counts into a file
    file_path = 'data_analysis/results/aa_change.missense.txt'  # save file
    header = [['initial', 'mutated', 'count']]
    aa_list = sorted([[key[0], key[1], val]
                      for key, val in aacounter.iteritems() if "*" not in key])
    csv.writer(open(file_path, 'wb'),
               delimiter='\t').writerows(header + aa_list)

    # re-slice the mutation data
    new_file_path = 'data_analysis/results/aa_change.properties.txt'
    df = pd.read_csv(file_path, sep='\t')
    # add properties of initial/mutated amino acids
    df['initial_prop'] = df['initial'].apply(lambda x: utils.python.letter_to_prop[x])
    df['mutated_prop'] = df['mutated'].apply(lambda x: utils.python.letter_to_prop[x])
    ptable = pd.pivot_table(df,
                            values='count',
                            rows='initial_prop',
                            cols='mutated_prop',
                            aggfunc=sum)
    ptable.to_csv(new_file_path, sep='\t')


def mut_type_cts_by_gene_type(file_path='data_analysis/results/gene_design_matrix.txt'):
    """Returns protein mutation type counts by gene type (oncogenes, tsg, other).

    Kwargs:
        file_path (str): path to mutation type cts by gene file

    Returns:
        pd.DataFrame: mutation type counts by gene type
    """
    df = pd.read_csv(file_path, sep='\t')
    df['gene_type'] = df['gene'].apply(_utils.classify_gene)
    mut_ct_df = df.iloc[:,1:]  # remove the "gene" column
    mut_ct_df = mut_ct_df.groupby('gene_type').sum()  # get counts for each gene type
    return mut_ct_df


def main():
    # count mutation types
    conn = get_cosmic_db()
    mut_cts = count_aa_mutation_types(conn)  # all mutation cts
    mut_cts.to_csv('data_analysis/results/aa_mut_type_cts.txt', sep='\t')
    plot_data.aa_mutation_types_barplot(mut_cts)
    onco_mut_cts = count_onco_aa_mut_types(conn)  # oncogene mutation cts
    onco_mut_cts.to_csv('data_analysis/results/aa_onco_mut_type_cts.txt', sep='\t')
    plot_data.aa_mutation_types_barplot(onco_mut_cts,
                                        save_path='data_analysis/plots'
                                        '/aa_onco_mut_types.barplot.png',
                                        title='Oncogene Protein Mutations'
                                        ' By Type')
    tsg_mut_cts = count_tsg_aa_mut_types(conn)
    tsg_mut_cts.to_csv('data_analysis/results/aa_tsg_mut_type_cts.txt', sep='\t')
    plot_data.aa_mutation_types_barplot(tsg_mut_cts,
                                        save_path='data_analysis/plots'
                                        '/aa_tsg_mut_types.barplot.png',
                                        title='Tumor Suppressor Protein '
                                        'Mutations By Type')

    # gene mutation counts
    gene_ct_df = count_gene_mutations(conn)
    gene_ct_df.set_index('Gene').to_csv('data_analysis/results/gene_mutation_counts.txt', sep='\t')
    plot_data.gene_mutation_histogram(gene_ct_df['count'])
    plot_data.cumulative_gene_mutation(gene_ct_df['count'])

    # get design matrix
    design_matrix = generate_design_matrix(conn)
    with open('data_analysis/results/gene_design_matrix.txt', 'wb') as handle:
        csv.writer(handle, delimiter='\t').writerows(design_matrix)
    plot_data.pca_plot()

    # plot protein mutation type counts by gene type
    tmp_mut_df = mut_type_cts_by_gene_type()
    tmp_mut_df.to_csv('data_analysis/results/gene_mutation_counts_by_gene_type.txt', sep='\t')
    plot_data.all_mut_type_barplot(tmp_mut_df)
    conn.close()

    with get_cosmic_db() as cursor:
        # handle missense mutation data
        aa_counter = count_aa_missense_changes(cursor)
        save_aa_missense_counts(aa_counter)
        plot_data.aa_missense_heatmap()
        plot_data.aa_property_heatmap()
        plot_data.aa_property_barplot()


if __name__=="__main__":
    main()
