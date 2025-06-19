from bw2data import Database
from bw2data.parameters import DatabaseParameter, ParameterManager, ProjectParameter
from bw2io import BW2Package

from lca_algebraic.log import warn
from lca_algebraic.params import _listParams, loadParams

__all__ = ["export_db", "import_db"]


def _param_data(param):
    """Return param data except id"""
    res = {key: val for key, val in param.__data__.items() if key != "id"}
    data = res.pop("data")
    res.update(data)

    return res


def export_db(db_name, filename):
    """This function exports a database to the **BW2Package** format, including the definition of parameters"""
    db = Database(db_name)
    db_params = DatabaseParameter.select().where(DatabaseParameter.database == db_name)

    # Export Db params
    db.metadata["database_parameters"] = [_param_data(param) for param in db_params]

    # List of all project params used in this dataset
    used_project_params = list(
        param.name for param in _listParams(db_name) if param.dbname is None
    )

    if len(used_project_params) > 0:
        warn(
            "Warning : this DB uses project parameters that are exported as well "
            "and might override project params at import time : ",
            used_project_params,
        )

        proj_params = list(
            ProjectParameter.get(ProjectParameter.name == name)
            for name in used_project_params
        )

        db.metadata["project_parameters"] = [
            _param_data(param) for param in proj_params
        ]

    BW2Package._write_file(filename, [BW2Package._prepare_obj(db, False)])


def import_db(filename):
    """Import Db from BW2Package with linked parameters (as produced by **export_db**)"""

    db = BW2Package.import_file(filename)[0]
    param_mgr = ParameterManager()
    if "database_parameters" in db.metadata:
        params = db.metadata["database_parameters"]
        param_mgr.new_database_parameters(params, db.name)

    if "project_parameters" in db.metadata:
        params = db.metadata["project_parameters"]
        param_mgr.new_project_parameters(params)

    # Reload the parameters
    loadParams()

    return db
