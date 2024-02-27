import brightway2 as bw
from typing_extensions import deprecated

from .helpers import setBackground, setForeground
from .log import logger


def deleteDb(db_name):
    del bw.databases[db_name]


def resetDb(db_name, foreground=True):
    """Create or cleanup a user DB"""
    if db_name in bw.databases:
        logger.warning("Db %s was here. Reseting it" % db_name)
        del bw.databases[db_name]

    db = bw.Database(db_name)
    db.write(dict())
    if foreground:
        setForeground(db_name)
    else:
        setBackground(db_name)


@deprecated("DEPRECATED : Use bw2io.import_ecoinvent_release() instead")
def initProject(project_name):
    """Setup the project if not already done."""
    bw.projects.set_current(project_name)
    bw.bw2setup()


@deprecated("DEPRECATED : Use the new bw2io.import_ecoinvent_release instead")
def importDb(dbname, path, parallel=False):
    """Import eco invent DB

    DEPRECATED : Use the new bw2io.import_ecoinvent_release instead
    """
    if dbname in bw.databases:
        logger.warning("Database '%s' has already been imported " % dbname)
    else:
        ei34 = bw.SingleOutputEcospold2Importer(path, dbname, use_mp=parallel)
        ei34.apply_strategies()
        ei34.statistics()
        ei34.write_database()
