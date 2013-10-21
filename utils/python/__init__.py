import csv
import os

# parent directory of this file
file_dir = os.path.dirname(os.path.abspath(__file__))

# columns of aa_properties.txt
NAME_COL = 0
LETTER_COL = 1
THREE_LETTER_COL = 2
PROPERTY_COL = 3

# read in properties of amino acids
aa_prop = list(csv.reader(open(os.path.join(file_dir, 'aa_properties.txt'), 'r'),
                          delimiter='\t'))[1:]

# convenience dicts for mapping amino acid letter to other atributes
one_to_three_letter = dict([[row[LETTER_COL], row[THREE_LETTER_COL]]
                            for row in aa_prop])
letter_to_prop = dict([[row[LETTER_COL], row[PROPERTY_COL]]
                       for row in aa_prop])
letter_to_name = dict([[row[LETTER_COL], row[NAME_COL]]
                       for row in aa_prop])
