from ..utils.python.amino_acid import AminoAcid


def test_missense_mutation():
    """Tests if AminoAcid properly parses missense HGVS syntax."""
    # test strings
    missense1 = 'p.A391E'
    missense2 = 'a391e'
    missense3 = 'p.*757?'

    # case 1 -- normal missense mutation
    aa = AminoAcid(missense1)
    assert aa.is_missense, 'Did not correctly detect a missense mutation'
    assert aa.initial == 'A', 'Did not correctly parse initial AA'
    assert aa.pos == 391, 'Did not correctly parse missense location'
    assert aa.mutated == 'E', 'Did not correctly parse mutated AA'

    # case 2 -- missense without the p.
    aa.set_amino_acid(missense2)
    assert aa.is_missense, 'Likely issue with parsing non-standard HGVS'
    assert aa.initial == 'A', 'Likely issue with parsing non-standard HGVS'
    assert aa.pos == 391, 'Likely issue with parsing non-standard HGVS'
    assert aa.mutated == 'E', 'Likely issue with parsing non-standard HGVS'

    # case 3 -- nonsense to unknown missense
    aa = AminoAcid(missense3)
    assert aa.is_missense
    assert aa.initial == '*'
    assert aa.mutated == '?'
    assert aa.is_missing_info


def test_insertion_mutation():
    """Tests if AminoAcid properly parses insertions using HGVS syntax."""
    # test strings
    ins1 = 'p.L601_K602ins27'
    ins2 = 'p.k313_v314insS'
    ins3 = 'p.?_?ins?'

    # case 1 -- num of inserted bases
    aa = AminoAcid(ins1)
    assert aa.is_insertion, 'Did not detect an insertion.'
    assert aa.initial == ['L', 'K'], 'Initial AA: ["L", "K"] != ' + str(aa.initial)
    assert aa.pos == [601, 602], 'Insertion pos: [601, 602] != ' + str(aa.pos)
    assert aa.mutated == '27', 'Did not correctly parse inserstion mutation'

    # case 2 -- inserted base info
    aa = AminoAcid(ins2)
    assert aa.is_insertion, 'Did not detect an insertion.'
    assert aa.initial == ['K', 'V'], 'Uncorrectly handled lower case.'
    assert aa.pos == [313, 314], 'Did not correctly parse insertion position'
    assert aa.mutated == 'S', 'Did not correctly parse inserstion mutation'

    # case 3 -- mising information
    aa = AminoAcid(ins3)
    assert aa.is_missing_info, 'Problem with detecting missing information.'
    assert aa.is_insertion, 'Problem with missing values.'
    assert aa.initial == '', 'Problem with missing values.'
    assert aa.pos == [], 'Problem with missing values.'
    assert aa.mutated == '', 'Problem with missing values.'


def test_deletion_mutation():
    """Tests if AminoAcid properly parses deletions using HGVS syntax."""
    # test strings
    del1 = 'p.S16_Q20del'
    del2 = 'p.?del'
    del3 = 'p.S45del'
    del4 = 'p.T155_A161delTRVRAMA'

    # case 1 -- normal
    aa = AminoAcid(del1)
    assert aa.is_deletion, 'Problem with detecting deletions.'
    assert aa.initial == ['S', 'Q'], 'Problem with parsing initial AA.'
    assert aa.pos == [16, 20], 'Problem with parsing deletion positon.'
    assert aa.mutated == '', 'Problem with parsing self.mutated'

    # case 2 -- missing info
    aa = AminoAcid(del2)
    assert aa.is_missing_info, 'Problem detecting missing information'

    # case 3 -- single base deletion
    aa = AminoAcid(del3)
    assert aa.initial == ['S'], 'Problem parsing single bp deletion'
    assert aa.pos == [45], 'Problem parsing single bp deletion'

    # case 4 -- deleted base information
    aa = AminoAcid(del4)
    assert aa.mutated == 'TRVRAMA', 'Did not record deleted base information'


def test_frame_shift_mutation():
    """Tests if AminoAcid properly parses frame shifts using HGVS syntax."""
    # test strings
    fs1 = 'p.M1fs'
    fs2 = 'p.W288fs*12'

    # case 1 -- frame shift without premature stop codon
    aa = AminoAcid(fs1)
    assert aa.is_frame_shift, 'Problem with detecting frame shift events.'
    assert aa.initial == 'M'
    assert aa.mutated == ''
    assert aa.pos == 1
    assert aa.stop_pos == None

    # case 2 -- frame shift with premature stop codon
    aa = AminoAcid(fs2)
    assert aa.initial == 'W'
    assert aa.pos == 288
    assert aa.stop_pos == 12


def test_nonsense_mutation():
    """Tests if AminoAcid properly parses nonsense mutations."""
    # test strings
    non1 = 'p.R943*'  # normal nonsense mutation
    non2 = 'p.*829*'  # synonymous nonsense mutation

    # case 1 -- normal nonsense mutation
    aa = AminoAcid(non1)
    assert aa.is_nonsense_mutation, 'Should detect nonsense mutation'
    assert aa.initial == 'R'
    assert aa.mutated == '*'
    assert aa.pos == 943, '943 != ' + str(aa.pos)

    # case 2 -- synonymous nonsense mutation
    aa = AminoAcid(non2)
    assert aa.is_nonsense_mutation
    assert aa.initial == '*'
    assert aa.mutated == '*'
    assert aa.is_synonymous
    assert aa.pos == 829



