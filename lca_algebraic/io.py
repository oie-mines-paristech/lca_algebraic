from os.path import basename, dirname

import bw2data
from bw2data.parameters import DatabaseParameter, ProjectParameter
from bw2io import BW2Package
import brightway2 as bw

from lca_algebraic import loadParams, error
from lca_algebraic.params import _listParams


__all__ = ['export_db', 'import_db']

def param_data(param) :
    """Return param data except id"""
    res =  {key:val for key, val in param.__data__.items() if key != "id"}
    data = res.pop("data")
    res.update(data)

    return res

def export_db(db_name, filename) :
    """Export Db and linked parameters"""
    db = bw.Database(db_name)
    db_params = DatabaseParameter.select().where(DatabaseParameter.database == db_name)

    # Export Db params
    db.metadata["database_parameters"] = [param_data(param) for param in db_params]

    # List of all project params used in this dataset
    used_project_params = list(param.name for param in _listParams(db_name) if param.dbname is None)

    if len(used_project_params) > 0 :
        error('Warning : this DB uses project parameters that are exported as well and might override project params at import time : ', used_project_params)

        proj_params = list(ProjectParameter.get(ProjectParameter.name==name) for name in used_project_params)

        db.metadata["project_parameters"] = [param_data(param) for param in proj_params]

    BW2Package._write_file(filename, [BW2Package._prepare_obj(db, False)])



def import_db(filename) :
    """Export Db and linked parameters"""

    db = BW2Package.import_file(filename)[0]
    if "database_parameters" in db.metadata :
        params = db.metadata["database_parameters"]
        bw.parameters.new_database_parameters(params, db.name)

    if "project_parameters" in db.metadata:
        params = db.metadata["project_parameters"]
        bw.parameters.new_project_parameters(params)

    # Reload the parameters
    loadParams()

    return db

