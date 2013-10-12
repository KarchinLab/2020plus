import re

class AminoAcid(object):
    """ The AminoAcid class represents aa changes in the Cosmic Database.
    """

    def __init__(self, amino_acid='', occurrence=0):
        if re.search('^[A-Za-z]\d+[A-Za-z]$', amino_acid):
            self.set_amino_acid(amino_acid)
            if occurrence:
                self.occurrence = occurrence
            self.is_valid = True
        else:
            self.is_valid = False

    def set_occurrence(self, occur):
        self.occurrence = occur

    def set_amino_acid(self, aa):
        """Set amino acid change and position."""
        self.initial, self.pos, self.mutated = self.__read_amino_acid_change(aa)

    def __read_amino_acid_change(self, aa_change):
        """Convert A207T syntax for amino acid change to tuple ("A", 207, "T").

        Args:
            aa_change (str): amino acid string ("A207T")

        Returns:
            tuple. The order::

                0 -- initial aa
                1 -- aa position
                2 -- mutated aa
        """
        try:
            initial_aa = aa_change[0].upper()  # starting aa
            pos = int(aa_change[1:-1])
            mutated_aa = aa_change[-1].upper()  # mutated aa
        except:
            print aa_change
            raise
        return initial_aa, pos, mutated_aa
