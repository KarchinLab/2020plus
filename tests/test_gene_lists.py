from ..utils.python import util as _utils
from ..utils.python.cosmic_db import get_cosmic_db
import pandas.io.sql as psql


def test_vogelstein_oncogenes():
    """Tests how many oncogenes in vogelstein's science paper is
    found in the COSMIC_nuc database"""
    onco_list = _utils.oncogene_list  # oncogenes according to paper
    num_oncogenes = len(onco_list)  # correspondig number of oncogenes

    # query COSMIC_nuc
    conn = get_cosmic_db()
    sql = ("SELECT Count(DISTINCT(Gene)) as NumOncoFound FROM `nucleotide` "
           "WHERE Gene in " + str(tuple(onco_list)))
    df = psql.frame_query(sql, con=conn)
    num_found_oncogenes = df['NumOncoFound'][0]
    conn.close()

    print 'ONCOGENES: Number found (%d) is less than total (%s)' % (num_found_oncogenes,
                                                                    num_oncogenes)
    assert num_oncogenes == num_found_oncogenes


def test_vogelstein_tsg():
    """Tests how many tsg in vogelstein's science paper are
    found in the COSMIC_nuc database"""
    tsg_list = _utils.tsg_list  # oncogenes according to paper
    num_tsg = len(tsg_list)  # correspondig number of oncogenes

    # query COSMIC_nuc
    conn = get_cosmic_db()
    sql = ("SELECT Count(DISTINCT(Gene)) as NumTsgFound FROM `nucleotide` "
           "WHERE Gene in " + str(tuple(tsg_list)))
    df = psql.frame_query(sql, con=conn)
    num_found_tsg = df['NumTsgFound'][0]
    conn.close()

    print 'TSG: Number found (%d) is less than total (%s)' % (num_found_tsg,
                                                                    num_tsg)
    assert num_tsg == num_found_tsg
