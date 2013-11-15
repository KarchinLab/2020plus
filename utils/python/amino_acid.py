import re
import logging


class AminoAcid(object):
    """ The AminoAcid class represents aa changes in the Cosmic Database.

    The AminoAcid class follows the syntax of HGVS
    (http://www.hgvs.org/mutnomen/recs-prot.html). Although the parsing
    generally follows the HGVS syntax, it does have slight variations where
    the COSMIC database uses idiosyncratic conventions.

    The AminoAcid class is only intended to be used in the following way::

        >>> aa = AminoAcid('p.A267C')
        >>> aa.pos
        267
        >>> aa.is_missense
        True

    Namely, the constructor parses the necessary HGVS string and then extracts
    attributes that can be used.
    """

    def __init__(self, hgvs='', occurrence=1):
        self.logger = logging.getLogger(__name__)
        if not (type(hgvs) is str or type(hgvs) is type(u'')):
            # catches cases where wierd non-string input is used
            self.is_valid = False
            self.set_mutation_type()
            print type(hgvs)
        else:
            # expected case of string
            self.hgvs_original = hgvs
            hgvs = hgvs.upper().replace('>', '')  # convert everything to upper case
            self.hgvs = hgvs
            self.occurrence = occurrence
            self.set_amino_acid(hgvs)
            self.set_mutation_type()

    def set_mutation_type(self, mut_type=''):
        """Sets the mutation type attribute to a single label based on
        attribute flags.

        Kwargs:
            mut_type (str): value to set self.mut_type
        """
        if mut_type:
            # user specifies a mutation type
            self.mutation_type = mut_type
        else:
            # mutation type is taken from object attributes
            if not self.is_valid:
                # does not correctly fall into a category
                self.mutation_type = 'not valid'
            elif self.is_missing_info:
                # mutation has a ?
                self.mutation_type = 'missing'
            else:
                # valid mutation type to be counted
                if self.is_synonymous:
                    # synonymous must go before missense since mutations
                    # can be categorized as synonymous and missense. Although
                    # in reality such cases are actually synonymous and not
                    # missense mutations.
                    self.mutation_type = 'synonymous'
                elif self.is_missense:
                    self.mutation_type = 'missense'
                elif self.is_indel:
                    self.mutation_type = 'indel'
                elif self.is_nonsense_mutation:
                    self.mutation_type = 'nonsense'
                elif self.is_frame_shift:
                    self.mutation_type = 'frame shift'

    def set_occurrence(self, occur):
        self.occurrence = occur

    def set_amino_acid(self, aa):
        """Set amino acid change and position."""
        aa = aa.upper()  # make sure it is upper case
        aa = aa[2:] if aa.startswith('P.') else aa  # strip "p."
        self.__set_mutation_status()  # set flags detailing the type of mutation
        self.__parse_hgvs_syntax(aa)  # read in specific mutations

    def __set_mutation_status(self):
        # strip "p." from HGVS protein syntax
        hgvs_tmp = self.hgvs[2:] if self.hgvs.startswith("P.") else self.hgvs

        # set evidence status
        self.__set_unkown_effect(hgvs_tmp)  # unknown effect
        self.__set_no_protein(hgvs_tmp)  # no protein
        self.__set_mutation_type(hgvs_tmp)  # indel, missense, etc.

    def __set_mutation_type(self, hgvs_string):
        """Interpret the mutation type (missense, etc.) and set appropriate flags.

        Args:
            hgvs_string (str): hgvs syntax with "p." removed
        """
        self.__set_missense_status(hgvs_string)  # missense mutations
        self.__set_indel_status()  # indel mutations
        self.__set_frame_shift_status()  # check for fs
        self.__set_premature_stop_codon_status(hgvs_string)  # check for stops

    def __set_missense_status(self, hgvs_string):
        """Sets the self.is_missense flag."""
        # set missense status
        if re.search('^[A-Z?*]\d+[A-Z?]$', hgvs_string):
            self.is_missense = True
        else:
            self.is_missense = False

    def __set_frame_shift_status(self):
        """Check for frame shift and set the self.is_frame_shift flag."""
        if 'fs' in self.hgvs_original:
            self.is_frame_shift = True
        elif re.search('[*A-Z]\d+[A-Z]+\*', self.hgvs_original):
            # it looks like some mutations dont follow the convention
            # of using 'fs' to indicate frame shift
            self.is_frame_shift = True
        else:
            self.is_frame_shift = False

    def __set_premature_stop_codon_status(self, hgvs_string):
        """Set whether there is a premature stop codon."""
        if re.search('.+\*(\d+)?$', hgvs_string):
            self.is_premature_stop_codon = True

            # check if it is also a nonsense mutation
            if hgvs_string.endswith('*'):
                self.is_nonsense_mutation = True
            else:
                self.is_nonsense_mutation = False
        else:
            self.is_premature_stop_codon = False
            self.is_nonsense_mutation = False

    def __set_indel_status(self):
        """Sets flags related to the mutation being an indel."""
        # set indel status
        if "ins" in self.hgvs_original:
            # mutation is insertion
            self.is_insertion = True
            self.is_deletion = False
            self.is_indel = True
        elif "del" in self.hgvs_original:
            # mutation is deletion
            self.is_deletion = True
            self.is_insertion = False
            self.is_indel = True
        else:
            # not an indel
            self.is_deletion = False
            self.is_insertion = False
            self.is_indel = False

    def __set_unkown_effect(self, hgvs_string):
        """Sets a flag for unkown effect according to HGVS syntax. The
        COSMIC database also uses unconventional questionmarks to denote
        missing information.

        Args:
            hgvs_string (str): hgvs syntax with "p." removed
        """
        # Standard use by HGVS of indicating unknown effect.
        unknown_effect_list = ['?', '(=)', '=']  # unknown effect symbols
        if hgvs_string in unknown_effect_list:
            self.unknown_effect = True
        elif "(" in hgvs_string:
            # parethesis in HGVS indicate expected outcomes
            self.unknown_effect = True
        else:
            self.unknown_effect = False

        # detect if there are missing information. commonly COSMIC will
        # have insertions with p.?_?ins? or deleteions with ?del indicating
        # missing information.
        if "?" in hgvs_string:
            self.is_missing_info = True
        else:
            self.is_missing_info = False

    def __set_no_protein(self, hgvs_string):
        """Set a flag for no protein expected.

        Args:
            hgvs_string (str): hgvs syntax with "p." removed
        """
        no_protein_list = ['0', '0?']  # no protein symbols
        if hgvs_string in no_protein_list:
            self.no_protein = True
        else:
            self.no_protein = False

    def __parse_hgvs_syntax(self, aa_hgvs):
        """Convert HGVS syntax for amino acid change into attributes.

        Specific details of the mutation are stored in attributes like
        self.intial (prior to mutation), sel.pos (mutation position),
        self.mutated (mutation), and self.stop_pos (position of stop codon,
        if any).

        Args:
            aa_hgvs (str): amino acid string following HGVS syntax
        """
        self.is_valid = True  # assume initially the syntax is legitimate
        self.is_synonymous = False  # assume not synonymous until proven
        if self.is_missense:
            self.initial = aa_hgvs[0]
            self.mutated = aa_hgvs[-1]
            self.pos = int(aa_hgvs[1:-1])
            self.stop_pos = None  # not a nonsense mutation
            if self.initial == self.mutated:
                self.is_synonymous = True
        elif self.is_indel:
            if self.is_insertion:
                if not self.is_missing_info:
                    self.initial = re.findall('([A-Z])\d+', aa_hgvs)[:2]  # first two
                    self.pos = map(int, re.findall('[A-Z](\d+)', aa_hgvs)[:2])  # first two
                    self.mutated = re.findall('(?<=INS)[A-Z0-9?*]+', aa_hgvs)[0]
                    self.mutated = self.mutated.strip('?')  # remove the missing info '?'
                else:
                    self.initial = ''
                    self.pos = []
                    self.mutated = ''
            elif self.is_deletion:
                if not self.is_missing_info:
                    self.initial = re.findall('([A-Z])\d+', aa_hgvs)
                    self.pos = map(int, re.findall('[A-Z](\d+)', aa_hgvs))
                    self.mutated = re.findall('(?<=DEL)[A-Z]*', aa_hgvs)[0]
                else:
                    self.initial = ''
                    self.pos = []
                    self.mutated = ''
        elif self.is_frame_shift:
            self.initial = aa_hgvs[0]
            self.mutated = ''
            try:
                self.pos = int(re.findall('[A-Z*](\d+)', aa_hgvs)[0])
                if self.is_premature_stop_codon:
                    self.stop_pos = int(re.findall('\*>?(\d+)$', aa_hgvs)[0])
                else:
                    self.stop_pos = None
            except IndexError:
                # unconventional usage of indicating frameshifts will cause
                # index errors. For example, in some cases 'fs' is not used.
                # In other cases, either amino acids were not included or
                # just designated as a '?'
                self.logger.debug('(Parsing-Problem) frame shift hgvs string: "%s"' % aa_hgvs)
                self.pos = None
                self.stop_pos = None
                self.is_missing_info = True
        elif self.is_nonsense_mutation:
            self.initial = aa_hgvs[0]
            self.mutated = '*'  # there is actually a stop codon
            self.stop_pos = 0  # indicates same position is stop codon
            try:
                self.pos = int(aa_hgvs[1:-1])
            except ValueError:
                # wierd error of p.E217>D*
                self.is_valid = False
                self.logger.debug('(Parsing-Problem) Invalid HGVS Amino Acid '
                                'syntax: ' + aa_hgvs)
            if self.initial == self.mutated:
                # classify nonsense-to-nonsense mutations as synonymous
                self.is_synonymous = True
        else:
            self.is_valid = False  # did not match any of the possible cases
            self.logger.debug('(Parsing-Problem) Invalid HGVS Amino Acid '
                              'syntax: ' + aa_hgvs)
