import sys

from lca_algebraic.methods import findMethods, method_name

# The package re-exports bw 'methods', which hides the submodule attribute
methods_module = sys.modules["lca_algebraic.methods"]

# Classic bw2 / Activity Browser layout : (family, category, indicator)
CLASSIC = ("EF v3.1", "climate change", "global warming potential (GWP100)")

# bw2io.import_ecoinvent_release layout : (ecoinvent version, family, category, indicator)
ECOINVENT_RELEASE = ("ecoinvent-3.11", "EF v3.1", "climate change", "global warming potential (GWP100)")

# Some methods have only two elements, such as ReCiPe 2016
SHORT = ("ReCiPe 2016 v1.03, midpoint (H)", "climate change")


def test_find_methods_main_cat(monkeypatch):
    monkeypatch.setattr(methods_module, "methods", [CLASSIC, ECOINVENT_RELEASE, SHORT])

    assert findMethods("climate", mainCat="EF v3.1") == [CLASSIC, ECOINVENT_RELEASE]
    assert findMethods("climate", mainCat="ReCiPe 2016 v1.03, midpoint (H)") == [SHORT]
    assert findMethods("climate") == [CLASSIC, ECOINVENT_RELEASE, SHORT]
    assert findMethods(mainCat="EF v3.1") == [CLASSIC, ECOINVENT_RELEASE]


def test_method_name():
    assert method_name(CLASSIC) == "climate change - global warming potential (GWP100)"
    assert method_name(ECOINVENT_RELEASE) == "climate change - global warming potential (GWP100)"
    assert method_name(SHORT) == "ReCiPe 2016 v1.03, midpoint (H) - climate change"
