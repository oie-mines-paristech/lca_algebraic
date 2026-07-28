from typing import Any, Callable, Iterable, List, Protocol, Tuple, Type

import pandas as pd


class BrightwayBackend(Protocol):
    """Brightway 2 / 2.5 compatibility surface used by lca_algebraic."""

    Activity: Any
    Exchange: Any
    ExchangeDataset: Any
    ActivityProxyBase: Any
    LCIBackend: Any
    sqlite3_lci_db: Any
    dict_as_exchangedataset: Callable[..., Any]

    name: str

    @property
    def databases(self) -> Any:
        ...

    @property
    def projects(self) -> Any:
        ...

    @property
    def methods(self) -> Iterable:
        ...

    def Database(self, name: str) -> Any:
        ...

    def Method(self, key: Tuple) -> Any:
        ...

    def get_activity(self, key) -> Any:
        ...

    def iter_database(self, db_name: str) -> Iterable:
        ...

    def multi_lca(
        self,
        activities: List[dict],
        methods: List[Tuple],
        act_name_fn: Callable,
        method_name_fn: Callable,
    ) -> pd.DataFrame:
        ...

    def is_output_exchange(self, exc) -> bool:
        ...

    def exchange_type_for_sub_act(self, sub_act) -> str:
        ...

    def process_node_default(self) -> str:
        ...

    def production_edge_default(self) -> str:
        ...

    def activity_not_found(self) -> Type[Exception]:
        ...

    def new_database_parameters(self, params, db_name: str) -> None:
        ...

    def new_project_parameters(self, params) -> None:
        ...

    def copy_activity_skip_keys(self) -> set:
        ...

    def bw2setup(self) -> None:
        ...

    def SingleOutputEcospold2Importer(self, path, dbname, use_mp=False) -> Any:
        ...
