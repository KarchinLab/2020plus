import utils.python
from utils.python.cosmic_db import get_cosmic_db
from utils.python.amino_acid import AminoAcid
from utils.python.nucleotide import Nucleotide
import utils.python.util as _utils
import plot_data
import mutation_types as mut_type
import pandas as pd
import pandas.io.sql as psql
import csv
from collections import OrderedDict
import logging


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
            if aa.mutation_type == 'missense':
                gene_pos_counter.setdefault(aa.pos, 0)
                gene_pos_counter[aa.pos] += 1
            mut_type_ctr[aa.mutation_type] += 1

        recurrent_cts = sum([cts for cts in gene_pos_counter.values() if cts > 1])
        mut_type_ctr['missense'] -= recurrent_cts  # subtract off the recurrent missense
        design_matrix.append([gene, recurrent_cts] + list(mut_type_ctr.values()))
    header = [['gene', 'recurrent missense'] + list(mut_type_ctr)]
    return header + design_matrix


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


def count_nuc_changes(conn):
    sql = "SELECT Nucleotide FROM `nucleotide`"
    df = psql.frame_query(sql, con=conn)
    for nuc_hgvs in df['Nucleotide']:
        nuc = Nucleotide(nuc_hgvs)


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


def main():
    # count mutation types
    conn = get_cosmic_db()
    # handle DNA nucleotides
    mut_nuc_cts = mut_type.count_nucleotides(conn)
    mut_nuc_cts.to_csv(_utils.result_dir + 'nuc_mut_type_cts.txt', sep='\t')
    tmp_plot_path = _utils.plot_dir + 'nuc_mut_types.barplot.png'  # plot path
    plot_data.mutation_types_barplot(mut_nuc_cts,
                                     save_path= tmp_plot_path,
                                     title='DNA Mutations by Type')

    # handle amino acids
    mut_cts = mut_type.count_amino_acids(conn)  # all mutation cts
    mut_cts.to_csv(_utils.result_dir + 'aa_mut_type_cts.txt', sep='\t')
    plot_data.mutation_types_barplot(mut_cts,
                                     title='Protein Mutations by Type')

    # handle oncogene mutation types
    onco_aa_cts, onco_nuc_cts = mut_type.count_oncogenes(conn)  # oncogene mutation cts
    onco_aa_cts.to_csv(_utils.result_dir + 'aa_onco_mut_type_cts.txt', sep='\t')
    onco_nuc_cts.to_csv(_utils.result_dir + 'nuc_onco_mut_type_cts.txt', sep='\t')
    plot_data.mutation_types_barplot(onco_aa_cts,
                                     save_path=_utils.plot_dir + \
                                     '/aa_onco_mut_types.barplot.png',
                                     title='Oncogene Protein Mutations'
                                     ' By Type')
    plot_data.mutation_types_barplot(onco_nuc_cts,
                                     save_path=_utils.plot_dir + \
                                     '/nuc_onco_mut_types.barplot.png',
                                     title='Oncogene DNA Mutations'
                                     ' By Type')

    # handle tumor suppressor mutation types
    tsg_aa_cts, tsg_nuc_cts = mut_type.count_tsg(conn)
    tsg_aa_cts.to_csv(_utils.result_dir + 'aa_tsg_mut_type_cts.txt', sep='\t')
    tsg_nuc_cts.to_csv(_utils.result_dir + 'nuc_tsg_mut_type_cts.txt', sep='\t')
    plot_data.mutation_types_barplot(tsg_aa_cts,
                                     save_path=_utils.plot_dir + \
                                     '/aa_tsg_mut_types.barplot.png',
                                     title='Tumor Suppressor Protein '
                                     'Mutations By Type')
    plot_data.mutation_types_barplot(tsg_nuc_cts,
                                     save_path=_utils.plot_dir + \
                                     '/nuc_tsg_mut_types.barplot.png',
                                     title='Tumor Suppressor DNA '
                                     'Mutations By Type')

    # gene mutation counts
    gene_ct_df = count_gene_mutations(conn)
    gene_ct_df.set_index('Gene').to_csv('data_analysis/results/gene_mutation_counts.txt', sep='\t')
    plot_data.gene_mutation_histogram(gene_ct_df['count'])
    plot_data.cumulative_gene_mutation(gene_ct_df['count'])

    # get design matrix
    design_matrix = generate_design_matrix(conn)
    with open(_utils.result_dir + 'gene_design_matrix.txt', 'wb') as handle:
        csv.writer(handle, delimiter='\t').writerows(design_matrix)
    plot_data.pca_plot()

    # plot protein mutation type counts by gene type
    tmp_mut_df = mut_type.count_gene_types()
    tmp_mut_df.to_csv(_utils.result_dir + 'gene_mutation_counts_by_gene_type.txt',
                      sep='\t')
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
