#
# This file defines several utility functions above brightway2 to be used by notebooks
#

import lca_algebraic.helpers
from .base_utils import *
from .helpers import *
from .lca import *
from .stats import *
from .params import *
from .io import *

def deleteDb(db_name) :
    del bw.databases[db_name]

def resetDb(db_name, foreground=True):
    """ Create or cleanup a user DB"""
    if db_name in bw.databases:
        error("Db %s was here. Reseting it" % db_name)
        del bw.databases[db_name]
    db = bw.Database(db_name)
    db.write(dict())
    if foreground :
        setForeground(db_name)
    else:
        setBackground(db_name)


def initProject(project_name) :
    '''Setup the project if not already done.'''
    bw.projects.set_current(project_name)
    bw.bw2setup()

def initDb(project_name) :
    '''Deprecated : use initProject(...) '''
    error("Deprecated : use initProject")
    initProject(project_name)

def importDb(dbname, path, parallel=False):
    '''Import eco invent DB'''

    if dbname in bw.databases:
        error("Database '%s' has already been imported " % dbname)
    else:
        ei34 = bw.SingleOutputEcospold2Importer(path, dbname, use_mp=parallel)
        ei34.apply_strategies()
        ei34.statistics()
        ei34.write_database()




