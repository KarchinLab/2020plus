from ..utils.python.nucleotide import Nucleotide


def test_substitution():
    """Tests if the Nucleotide class can parse substitutions."""
    sub1 = 'c.1849G>T'
    sub2 = 'c.3082+1G>T'
    sub3 = 'c.1153_1159TGTAAAA>CGAACTTGTAAACTGTAA'

    # case 1 -- normal single substitution
    nuc = Nucleotide(sub1)
    assert nuc.is_substitution, 'Nucleotide did not detect substitution'
    assert nuc.pos == 1849
    assert nuc.initial == 'G'
    assert nuc.mutated == 'T'
    assert nuc.intron_pos is None

    # case 2 -- mutation in intron
    nuc = Nucleotide(sub2)
    assert nuc.is_substitution, 'Nucleotide did not detect substitution'
    assert nuc.pos == 3082
    assert nuc.intron_pos == 1, '+1 != %d' + nuc.intron_pos
    assert nuc.initial == 'G'
    assert nuc.mutated == 'T'

    # case 3 -- substitution of more than one base
    nuc = Nucleotide(sub3)
    assert nuc.is_substitution, 'Nucleotide did not detect substitution'
    assert nuc.pos == [1153, 1159]
    assert nuc.initial == 'TGTAAAA'
    assert nuc.mutated == 'CGAACTTGTAAACTGTAA'


def test_deletion():
    """Tests if the Nucleotide class can parse HGVS deletions."""
    del1 = 'c.3047_3048delGT'
    del2 = 'c.4248delC'
    del3 = 'c.?_?del?'
    del4 = 'c.2233_2247del15'

    # case 1 -- more than one nucleotide deletion
    nuc = Nucleotide(del1)
    assert nuc.is_deletion, 'Did not detect deletion'
    assert nuc.pos == [3047, 3048]
    assert nuc.initial == 'GT'
    assert nuc.mutated == ''

    # case 2 -- single nucleotide deletion
    nuc = Nucleotide(del2)
    assert nuc.is_deletion
    assert nuc.pos == 4248
    assert nuc.initial == 'C'
    assert nuc.mutated == ''
    assert nuc.intron_pos is None

    # case 3 -- missing info
    nuc = Nucleotide(del3)
    assert nuc.is_deletion
    assert nuc.is_missing_info
    assert nuc.pos == ['?', '?']
    assert nuc.initial == '?'

    # case 4 -- range deletion with length
    nuc = Nucleotide(del4)
    assert nuc.is_deletion
    assert nuc.initial == '15'


def test_insertion():
    """Tests if the Nucleotide class can parse HGVS insertions."""
    ins1 = 'c.3047_3048insGT'
    ins2 = 'c.4248insC'
    ins3 = 'c.?_?ins?'
    ins4 = 'c.2233_2247ins15'

    # case 1 -- more than one nucleotide insertion
    nuc = Nucleotide(ins1)
    assert nuc.is_insertion, 'Did not detect insertion'
    assert nuc.pos == [3047, 3048]
    assert nuc.mutated == 'GT'
    assert nuc.initial == ''

    # case 2 -- single nucleotide deletion
    nuc = Nucleotide(ins2)
    assert nuc.is_insertion
    assert nuc.pos == 4248
    assert nuc.mutated == 'C'
    assert nuc.initial == ''
    assert nuc.intron_pos is None

    # case 3 -- missing info
    nuc = Nucleotide(ins3)
    assert nuc.is_insertion
    assert nuc.is_missing_info
    assert nuc.pos == ['?', '?']
    assert nuc.mutated == '?'

    # case 4 -- range deletion with length
    nuc = Nucleotide(ins4)
    assert nuc.is_insertion
    assert nuc.mutated == '15'
