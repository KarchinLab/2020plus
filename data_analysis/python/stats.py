from utils.python.cosmic_db import get_cosmic_db
from utils.python.amino_acid import AminoAcid
import matplotlib.pyplot as plt
import csv

def count_mutations(cursor):
    """Count the number of entries"""
    cursor.execute("""SELECT COUNT(COSMICSampleID)
                   FROM `nucleotide`""")
    return cursor.fetchone()[0]  # COUNT query returns a tuple

def count_aa_changes(cursor):
    """Count amino acid change counts.

    Args:
        cursor: mysqldb cursor object

    Returns:
        dict. containing counts eg. {('aa1', 'aa2'): 4}
    """
    cursor.execute("""SELECT aachange, occurrences
                   FROM `cosmic_aa`""")
    aa_change_counter = {}
    for aachange, occurrences in cursor.fetchall():
        aa = AminoAcid(amino_acid=aachange,
                       occurrence=occurrences)
        if aa.is_valid:
            aa_change_counter.setdefault((aa.initial, aa.mutated), 0)
            aa_change_counter[(aa.initial, aa.mutated)] += aa.occurrence
    return aa_change_counter

def main():
    conn = get_cosmic_db()
    cursor = conn.cursor()
    # mycount = count_mutations(cursor)  # count mutations
    aa_counter = count_aa_changes(cursor)
    header = [['initial', 'mutated', 'count']]
    aa_list = sorted([[key[0], key[1], val] for key, val in aa_counter.iteritems()])
    # aa_list = sorted([['=>'.join(key), val] for key, val in aa_counter.iteritems()])
    csv.writer(open('data_analysis/results/aa_change.txt','wb'),
               delimiter='\t').writerows( header + aa_list)
    # change, count = zip(*aa_list)
    # x = range(len(change))
    # plt.bar(x, count)
    # plt.xticks(x, change)
    # plt.savefig('test.png')
    conn.close()

if __name__=="__main__":
    main()
