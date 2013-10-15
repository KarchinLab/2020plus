from utils.python.cosmic_db import get_cosmic_db
from utils.python.amino_acid import AminoAcid
import plot_data
import csv
import logging


def count_mutations(cursor):
    """Count the number of entries"""
    cursor.execute("""SELECT COUNT(COSMICSampleID)
                   FROM `nucleotide`""")
    return cursor.fetchone()[0]  # COUNT query returns a tuple


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


def main():
    conn = get_cosmic_db()
    cursor = conn.cursor()
    aa_counter = count_aa_missense_changes(cursor)
    header = [['initial', 'mutated', 'count']]
    aa_list = sorted([[key[0], key[1], val] for key, val in aa_counter.iteritems() if "*" not in key])
    csv.writer(open('data_analysis/results/aa_change.missense.txt', 'wb'),
               delimiter='\t').writerows(header + aa_list)
    plot_data.plot_aa_missense_change()
    conn.close()

if __name__=="__main__":
    main()
