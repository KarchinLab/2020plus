from ..utils.python.cosmic_db import get_cosmic_db

def test_get_cosmic_db():
    """Test if database connection to COSMIC_nuc DB is working."""
    conn = get_cosmic_db()
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM `nucleotide` LIMIT 10")
    print list(cursor.fetchall())
    conn.close
