import re
import logging


class Nucleotide(object):
    """The Nucleotide class represents DNA changes in the COSMIC database.

    The Nucleotide class follows the syntax of HGVS
    (http://www.hgvs.org/mutnomen/recs-DNA.html).
    """

    def __init__(self, hgvs='', occurrence=1,
                 len5ss=2, len3ss=-2):
        self.logger = logging.getLogger(__name__)
        self.occurrence = occurrence
        self.len5ss = len5ss  # 5' splice site len
        self.len3ss = len3ss  # 3' splice site len

        self.hgvs_original = hgvs  # unmodified hgvs copy
        self.hgvs = hgvs if not hgvs.startswith('c.') else hgvs[2:]  # modified
        self.set_nucleotide()
        self.set_mutation_type()

    def set_nucleotide(self, hgvs_nuc=''):
        hgvs_tmp = hgvs_nuc if hgvs_nuc else self.hgvs
        self.__set_unknown_effect(hgvs_tmp)  # completely unknown
        self.__set_missing_info(hgvs_tmp)  # has missing information
        self.__set_nucleotide_mutation(hgvs_tmp)  # set mutation type
        self.__parse_hgvs_syntax(hgvs_tmp)

    def set_mutation_type(self, mut_type=''):
        """Sets a string designating the mutation type.

        Args:
            mut_type (str): name of mutation type
        """
        if mut_type:
            # specified mutation type
            self.mutation_type = mut_type
        else:
            # interpret mutation type from attributes
            if not self.is_valid:
                # does not correctly fall into a category
                self.mutation_type = 'not valid'
            elif self.unknown_effect:
                self.mutation_type = 'unknown effect'
            elif self.is_missing_info:
                self.mutation_type = 'missing'
            elif self.is_substitution:
                self.mutation_type = 'substitution'
            elif self.is_deletion:
                self.mutation_type = 'deletion'
            elif self.is_insertion:
                self.mutation_type = 'insertion'

        # check if mutation at splice site
        self.__set_splice_mutation()

    def __set_splice_mutation(self):
        """Set the is_splicing_mutation flag"""
        #len5ss = 6  # positive number since 5SS
        #len3ss = -20  # use negative syntax like HGVS
        if type(self.intron_pos) == int:
            # SNV case, only one position
            if self.len3ss <= self.intron_pos <= self.len5ss:
                self.is_splicing_mutation = True
            else:
                self.is_splicing_mutation = False
        elif type(self.intron_pos) == list:
            # deletion case, now have a window to check overlap
            if self.intron_pos[0]:
                first_in_splice = self.len3ss <= self.intron_pos[0] <= self.len5ss
                tmp_pos1 = self.intron_pos[0]
            else:
                first_in_splice = False
                tmp_pos1 = 0
            if self.intron_pos[1]:
                second_in_splice = self.len3ss <= self.intron_pos[1] <= self.len5ss
                tmp_pos2 = self.intron_pos[1]
            else:
                second_in_splice = False
                tmp_pos2 = 0

            # set splice site mutation flag
            if first_in_splice or second_in_splice:
                self.is_splicing_mutation = True
            elif (tmp_pos1 == 0 and tmp_pos2 > self.len5ss) or (tmp_pos1 < self.len3ss and tmp_pos2 == 0):
                self.is_splicing_mutation = True
            else:
                self.is_splicing_mutation = False
        else:
            self.is_splicing_mutation = False

    def __set_unknown_effect(self, hgvs_str):
        """Sets a flag for unkown effect (c.? or ?).

        Note: Unavailable information according to HGVS is usually
        marked with a c.?, ?, or parethesis.

        Args:
            hgvs_str (str): DNA HGVS string
        """
        unknown_effect_list = ['c.?', '?']
        if hgvs_str.lower() in unknown_effect_list:
            self.unknown_effect = True
        elif hgvs_str.startswith("("):
            self.unknown_effect = True
        else:
            self.unknown_effect = False

    def __set_missing_info(self, hgvs_str):
        """Sets a flag for missing data (? in HGVS syntax).

        Args:
            hgvs_str (str): DNA HGVS string
        """
        if '?' in hgvs_str:
            self.is_missing_info = True
        else:
            self.is_missing_info = False

    def __set_nucleotide_mutation(self, hgvs_str):
        """Interpret the HGVS syntax and set appropriate mutation type
        attributes (substitution, insertion, etc.).

        Args:
            hgvs_str (str): string representing HGVS DNA mutation (no "c.")
        """
        self.__set_substitution_status(hgvs_str)
        self.__set_indel_status(hgvs_str)

    def __set_substitution_status(self, hgvs_str):
        self.is_substitution = '>' in hgvs_str

    def __set_indel_status(self, hgvs_str):
        """Sets attribute flags for whether mutation is a insertion, deletion,
        and indel.

        """
        # set deletion status
        self.is_deletion = 'del' in hgvs_str

        # set insertion status
        self.is_insertion = 'ins' in hgvs_str

        # set indel status
        if self.is_insertion or self.is_deletion:
            self.is_indel = True
        else:
            self.is_indel = False

    def __parse_hgvs_syntax(self, hgvs_str):
        """Parse the HGVS DNA mutation string to set attributes.

        Look at tests/test_nucleotide.py for examples on how
        specific HGVS strings should be parsed.

        Args:
            hgvs_str (str): DNA HGVS string
        """
        self.is_valid = True  # assume initially the syntax is valid
        if self.is_substitution:
            sub_pattern = '(?:(\d+)([+-]\d+)?_)?(\d+)([+-]\d+)?([A-Z]+)>([A-Z]+)$'
            matches = re.findall(sub_pattern, hgvs_str)
            if matches:
                init_pos, init_intron, reg_pos, reg_intron, initial, mutated = matches[0]
                if not init_pos:
                    self.pos = int(reg_pos)
                    self.intron_pos = int(reg_intron) if reg_intron != '' else None
                    self.initial = initial
                    self.mutated = mutated
                else:
                    init_pos = init_pos.strip('_')  # remove separating underscore
                    self.pos = [int(init_pos), int(reg_pos)]
                    intron_tmp1 = int(init_intron) if init_intron != '' else None
                    intron_tmp2 = int(reg_intron) if reg_intron != '' else None
                    self.intron_pos = [intron_tmp1, intron_tmp2]
                    self.initial = initial
                    self.mutated = mutated
            else:
                self.is_valid = False
                self.intron_pos = None
                self.logger.debug('(Parsing-Problem) Invalid DNA Substitution: ' + hgvs_str)
                return
        elif self.is_deletion:
            del_pattern = '(?:([0-9?]+)([-+]\d+)?(?:_))?([0-9?]+)([-+]\d+)?del([A-Z?0-9]+)$'
            matches = re.findall(del_pattern, hgvs_str)
            if matches:
                init_pos, init_intron, reg_pos, reg_intron, del_nuc = matches[0]
                if not init_pos:
                    # only one nucleotide deleted
                    self.pos = int(reg_pos) if reg_pos != '?' else reg_pos
                    self.intron_pos = int(reg_intron) if reg_intron != '' else None
                    self.mutated = ''
                    self.initial = del_nuc
                else:
                    # more than one nucleotide deleted
                    init_pos = init_pos.strip('_')  # remove '_' because of regex
                    pos1 = int(init_pos) if init_pos != '?' else init_pos
                    pos2 = int(reg_pos) if reg_pos != '?' else reg_pos
                    self.pos = [pos1, pos2]
                    intron_tmp1 = int(init_intron) if init_intron != '' else None
                    intron_tmp2 = int(reg_intron) if reg_intron != '' else None
                    self.intron_pos = [intron_tmp1, intron_tmp2]
                    self.mutated = ''
                    self.initial = del_nuc
            else:
                self.intron_pos = False
        elif self.is_insertion:
            ins_pattern = '(?:([0-9?]+)([-+]\d+)?(?:_))?([0-9?]+)([-+]\d+)?ins([A-Z?0-9]+)$'
            matches = re.findall(ins_pattern, hgvs_str)
            if matches:
                init_pos, init_intron, reg_pos, reg_intron, ins_nuc = matches[0]
                if not init_pos:
                    # only one nucleotide inserted
                    self.pos = int(reg_pos) if reg_pos != '?' else reg_pos
                    self.intron_pos = int(reg_intron) if reg_intron != '' else None
                    self.initial = ''
                    self.mutated = ins_nuc
                else:
                    # more than one nucleotide inserted
                    init_pos = init_pos.strip('_')  # remove '_' because of regex
                    pos1 = int(init_pos) if init_pos != '?' else init_pos
                    pos2 = int(reg_pos) if reg_pos != '?' else reg_pos
                    self.pos = [pos1, pos2]
                    intron_tmp1 = int(init_intron) if init_intron != '' else None
                    intron_tmp2 = int(reg_intron) if reg_intron != '' else None
                    self.intron_pos = [intron_tmp1, intron_tmp2]
                    self.initial = ''
                    self.mutated = ins_nuc
            else:
                self.intron_pos = None
        elif self.unknown_effect:
            # unknown effect for mutation. usually denoted as c.?
            self.intron_pos = None
            return
        else:
            # mutation did not fall into any of the categories. thus it likely
            # has invalid syntax
            self.is_valid = False
            self.intron_pos = None
            self.logger.debug('(Parsing-Problem) Invalid HGVS DNA syntax: ' + hgvs_str)
            return
