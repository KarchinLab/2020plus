import MySQLdb
import ConfigParser
import os

def get_cosmic_db():
    """Return a connection to the cosmic database."""
    cfg_opts = get_config()  # read config file for mysql db
    conn = MySQLdb.connect(host=cfg_opts['host'],
                           user=cfg_opts['user'],
                           passwd=cfg_opts['passwd'],
                           db=cfg_opts['db'])
    return conn

def get_config():
    """Returns content of the database config file as a dict."""
    cfg = ConfigParser.ConfigParser()
    cfg.read(os.path.join(os.path.dirname(__file__),
                          'db.cfg'))
    cfg_options = dict(cfg.items('database'))
    return cfg_options

