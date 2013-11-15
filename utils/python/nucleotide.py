import re
import logging


class Nucleotide(object):
    """The Nucleotide class represents DNA changes in the COSMIC database.

    The Nucleotide class follows the syntax of HGVS
    (http://www.hgvs.org/mutnomen/recs-DNA.html).
    """

    def __init__(self, hgvs='', occurrence=1):
        self.logger = logging.getLogger(__name__)
        self.occurrence = occurrence

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

    def __set_unknown_effect(self, hgvs_str):
        unknown_effect_list = ['c.?', '?']
        if hgvs_str.lower() in unknown_effect_list:
            self.unknown_effect = True
        elif hgvs_str.startswith("("):
            self.unknown_effect = True
        else:
            self.unknown_effect = False

    def __set_missing_info(self, hgvs_str):
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

        """
        self.is_valid = True  # assume initially the syntax is valid
        if self.is_substitution:
            sub_pattern = '(?:(\d+)([+-]\d+)?_)?(\d+)([+-]\d+)?([A-Z]+)>([A-Z]+)$'
            # old_sub_pattern = '(\d+([+-]\d+)?_)?(\d+)([+-]\d+)?([A-Z]+)>([A-Z]+)$'
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
                self.logger.debug('(Parsing-Problem) Invalid DNA Substitution: ' + hgvs_str)
        elif self.is_deletion:
            del_pattern = '(?:([0-9?]+)([-+]\d+)?(?:_))?([0-9?]+)([-+]\d+)?del([A-Z?0-9]+)$'
            # old_del_pattern = '([0-9?]+([-+]\d+)?(?:_))?([0-9?]+)([-+]\d+)?del([A-Z?0-9]+)$'
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
        elif self.is_insertion:
            ins_pattern = '(?:([0-9?]+)([-+]\d+)?(?:_))?([0-9?]+)([-+]\d+)?ins([A-Z?0-9]+)$'
            # old_ins_pattern = '([0-9?]+([-+]\d+)?(?:_))?([0-9?]+)([-+]\d+)?ins([A-Z?0-9]+)$'
            matches = re.findall(ins_pattern, hgvs_str)
            if matches:
                init_pos, init_intron, reg_pos, reg_intron, ins_nuc = matches[0]
                if not init_pos:
                    # only one nucleotide inserted
                    self.pos = int(reg_pos) if reg_pos!= '?' else reg_pos
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
        elif self.unknown_effect:
            # unknown effect for mutation. usually denoted as c.?
            pass
        else:
            # mutation did not fall into any of the categories. thus it likely
            # has invalid syntax
            self.is_valid = False
            self.logger.debug('(Parsing-Problem) Invalid HGVS DNA syntax: ' + hgvs_str)

