from lca_algebraic import copyActivity, getActByCode, newActivity, resetDb
from test.conftest import BG_DB, USER_DB


def _links(act):
    """Inputs of `act`, excluding its own production : {input_key: amount}"""
    return {tuple(ex["input"]): ex["amount"] for ex in act.exchanges() if tuple(ex["input"]) != tuple(act.key)}


def _productions(act):
    return [ex for ex in act.exchanges() if tuple(ex["input"]) == tuple(act.key)]


def test_overwrite_preserves_incoming_links(data):
    target = copyActivity(USER_DB, data.bg_act1, code="target")
    newActivity(USER_DB, "c1", "kg", {target: 1})
    newActivity(USER_DB, "c2", "kg", {target: 2})

    # Overwrite target : consumers must keep their links
    copyActivity(USER_DB, data.bg_act2, code="target")

    assert _links(getActByCode(USER_DB, "c1")) == {(USER_DB, "target"): 1}
    assert _links(getActByCode(USER_DB, "c2")) == {(USER_DB, "target"): 2}


def test_overwrite_preserves_links_from_other_db(data):
    resetDb("fg2")
    target = copyActivity(USER_DB, data.bg_act1, code="target")
    other = newActivity("fg2", "other", "kg", {target: 3})

    copyActivity(USER_DB, data.bg_act2, code="target")

    assert _links(getActByCode("fg2", "other")) == {(USER_DB, "target"): 3}
    assert _links(other) == {(USER_DB, "target"): 3}


def test_overwrite_replaces_content(data):
    act = copyActivity(USER_DB, data.bg_act1, code="target")
    act["custom_flag"] = "stale"
    act.save()

    copyActivity(USER_DB, data.bg_act2, code="target")

    act = getActByCode(USER_DB, "target")
    assert _links(act) == {(BG_DB, "bio2"): 1}
    assert "custom_flag" not in act
    assert tuple(act["inherited_from"]) == data.bg_act2.key

    prods = _productions(act)
    assert len(prods) == 1
    assert prods[0]["amount"] == 1


def test_double_copy_is_idempotent(data):
    copyActivity(USER_DB, data.bg_act1, code="target")
    copyActivity(USER_DB, data.bg_act1, code="target")

    act = getActByCode(USER_DB, "target")
    assert _links(act) == {(BG_DB, "bio1"): 1}
    assert len(_productions(act)) == 1


def test_new_activity_overwrite_preserves_links(data):
    act = newActivity(USER_DB, "act", "kg", {data.bg_act1: 1})
    newActivity(USER_DB, "cons", "kg", {act: 1})

    newActivity(USER_DB, "act", "kg", {data.bio1: 2})

    assert _links(getActByCode(USER_DB, "cons")) == {(USER_DB, "act"): 1}

    act = getActByCode(USER_DB, "act")
    assert _links(act) == {(BG_DB, "bio1"): 2}
    assert len(_productions(act)) == 1
