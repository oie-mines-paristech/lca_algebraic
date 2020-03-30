#
# This file defines several utility functions above brightway2 to be used by notebooks
#

import lca_algebraic.helpers
from .base_utils import _eprint
from .base_utils import *
from .helpers import *
from .lca import *
from .stats import *
from .params import *

def resetDb(db_name):
    """ Create or cleanup a user DB"""
    if db_name in bw.databases:
        _eprint("Db %s was here. Reseting it" % db_name)
        del bw.databases[db_name]
    db = bw.Database(db_name)
    db.write(dict())

def initDb(project_name) :
    '''Init brightway and detect version of existing installation of ecoinvent'''
    bw.projects.set_current(project_name)
    bw.bw2setup()


def importDb(dbname, path):
    '''Import eco invent DB'''

    if dbname in bw.databases:
        _eprint("Database '%s' has already been imported " % dbname)
    else:
        ei34 = bw.SingleOutputEcospold2Importer(path, dbname)
        ei34.apply_strategies()
        ei34.statistics()
        ei34.write_database()




