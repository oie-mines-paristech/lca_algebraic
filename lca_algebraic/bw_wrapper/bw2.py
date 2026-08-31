from typing import Callable, List, Tuple

import brightway2 as bw
import pandas as pd
from bw2data.backends import LCIBackend
from bw2data.backends.peewee import (
    Activity,
    ActivityDataset,
    Exchange,
    ExchangeDataset,
    sqlite3_lci_db,
)
from bw2data.backends.peewee.utils import dict_as_exchangedataset
from bw2data.proxies import ActivityProxyBase
from peewee import DoesNotExist

_TECHNO_TYPES = [None, "process", "product", "processwithreferenceproduct", "multifunctional"]


class BW2Backend:
    name = "bw2"

    Activity = Activity
    ActivityDataset = ActivityDataset
    Exchange = Exchange
    ExchangeDataset = ExchangeDataset
    ActivityProxyBase = ActivityProxyBase
    LCIBackend = LCIBackend
    sqlite3_lci_db = sqlite3_lci_db
    dict_as_exchangedataset = staticmethod(dict_as_exchangedataset)

    @property
    def databases(self):
        return bw.databases

    @property
    def projects(self):
        return bw.projects

    @property
    def methods(self):
        return bw.methods

    def Database(self, name: str):
        return bw.Database(name)

    def Method(self, key: Tuple):
        return bw.Method(key)

    def get_activity(self, key):
        return bw.get_activity(key)

    def iter_database(self, db_name: str):
        return iter(bw.Database(db_name))

    def multi_lca(
        self,
        activities: List[dict],
        methods: List[Tuple],
        act_name_fn: Callable,
        method_name_fn: Callable,
    ) -> pd.DataFrame:
        bw.calculation_setups["process"] = {"inv": activities, "ia": methods}
        lca = bw.MultiLCA("process")
        cols = [act_name_fn(act) for act_amount in activities for act in act_amount]
        return pd.DataFrame(lca.results.T, index=[method_name_fn(method) for method in methods], columns=cols)

    def is_output_exchange(self, exc) -> bool:
        return exc.get("type") == "production"

    def exchange_type_for_sub_act(self, sub_act) -> str:
        if sub_act.get("type") in _TECHNO_TYPES:
            return "technosphere"
        return "biosphere"

    def process_node_default(self) -> str:
        return "process"

    def production_edge_default(self) -> str:
        return "production"

    def activity_not_found(self):
        return DoesNotExist

    def new_database_parameters(self, params, db_name: str) -> None:
        bw.parameters.new_database_parameters(params, db_name)

    def new_project_parameters(self, params) -> None:
        bw.parameters.new_project_parameters(params)

    def copy_activity_skip_keys(self) -> set:
        return {"database", "code"}

    def bw2setup(self):
        bw.bw2setup()

    def SingleOutputEcospold2Importer(self, path, dbname, use_mp=False):
        return bw.SingleOutputEcospold2Importer(path, dbname, use_mp=use_mp)
