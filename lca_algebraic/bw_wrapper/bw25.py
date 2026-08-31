from typing import Callable, List, Tuple

import bw2calc
import bw2data
import numpy as np
import pandas as pd
from bw2data import Database
from bw2data.backends import (
    Activity,
    ActivityDataset,
    Exchange,
    ExchangeDataset,
    SQLiteBackend,
    sqlite3_lci_db,
)
from bw2data.backends.utils import dict_as_exchangedataset
from bw2data.configuration import labels
from bw2data.errors import UnknownObject
from bw2data.parameters import ParameterManager
from bw2data.proxies import ActivityProxyBase


class BW25Backend:
    name = "bw25"

    Activity = Activity
    ActivityDataset = ActivityDataset
    Exchange = Exchange
    ExchangeDataset = ExchangeDataset
    ActivityProxyBase = ActivityProxyBase
    LCIBackend = SQLiteBackend
    sqlite3_lci_db = sqlite3_lci_db
    dict_as_exchangedataset = staticmethod(dict_as_exchangedataset)

    @property
    def databases(self):
        return bw2data.databases

    @property
    def projects(self):
        return bw2data.projects

    @property
    def methods(self):
        return bw2data.methods

    def Database(self, name: str):
        return Database(name)

    def Method(self, key: Tuple):
        return bw2data.Method(key)

    def get_activity(self, key):
        return bw2data.get_activity(key)

    def iter_database(self, db_name: str):
        return iter(Database(db_name))

    def multi_lca(
        self,
        activities: List[dict],
        methods: List[Tuple],
        act_name_fn: Callable,
        method_name_fn: Callable,
    ) -> pd.DataFrame:
        meth_cfg = {"impact_categories": methods}
        fu = {act_name_fn(act): {act.id: 1} for act_amount in activities for act in act_amount}

        bw2data.databases.clean()
        data_objs = bw2data.get_multilca_data_objs(fu, meth_cfg)
        lca = bw2calc.MultiLCA(demands=fu, method_config=meth_cfg, data_objs=data_objs)
        lca.lci()
        lca.lcia()

        rows = [method_name_fn(method) for method in methods]
        cols = [act_name_fn(act) for act_amount in activities for act in act_amount]
        results = lca.scores

        return pd.DataFrame(
            np.array([[results[m, a] for a in cols] for m in methods]),
            index=rows,
            columns=cols,
        )

    def is_output_exchange(self, exc) -> bool:
        return exc.get("type") == labels.production_edge_default

    def exchange_type_for_sub_act(self, sub_act) -> str:
        if sub_act.get("type") in labels.lci_node_types:
            return labels.consumption_edge_default
        return labels.biosphere_edge_default

    def process_node_default(self) -> str:
        return labels.process_node_default

    def production_edge_default(self) -> str:
        return labels.production_edge_default

    def activity_not_found(self):
        return UnknownObject

    def new_database_parameters(self, params, db_name: str) -> None:
        ParameterManager().new_database_parameters(params, db_name)

    def new_project_parameters(self, params) -> None:
        ParameterManager().new_project_parameters(params)

    def copy_activity_skip_keys(self) -> set:
        return {"database", "code", "id"}

    def bw2setup(self) -> None:
        raise NotImplementedError("bw2setup is only supported with brightway2")

    def SingleOutputEcospold2Importer(self, path, dbname, use_mp=False):
        raise NotImplementedError("SingleOutputEcospold2Importer is only supported with brightway2")
