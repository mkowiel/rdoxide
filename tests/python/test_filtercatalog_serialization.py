"""FilterCatalog / FilterCatalogEntry serialization (rdoxide's own binary format)."""

import pickle

from rdoxide import Chem
from rdoxide.Chem import FilterCatalog as FC


def test_can_serialize():
    assert FC.FilterCatalogCanSerialize()


def test_builtin_catalog_round_trip():
    params = FC.FilterCatalogParams()
    params.AddCatalog(FC.FilterCatalogParams.FilterCatalogs.PAINS_A)
    cat = FC.FilterCatalog(params)
    for other in (FC.FilterCatalog(cat.Serialize()), pickle.loads(pickle.dumps(cat))):
        assert other.GetNumEntries() == cat.GetNumEntries()
        assert other.GetEntryWithIdx(3).GetDescription() == cat.GetEntryWithIdx(3).GetDescription()
        mol = Chem.MolFromSmiles("O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2")
        assert [e.GetDescription() for e in other.GetMatches(mol)] == [e.GetDescription() for e in cat.GetMatches(mol)]


def test_entry_with_combined_matchers():
    a = FC.SmartsMatcher("benzene", "c1ccccc1", 1)
    b = FC.SmartsMatcher("chlorine", "[Cl]", 1)
    excl = FC.ExclusionList()
    excl.AddPattern(b)
    matcher = FC.FilterMatchOps.And(a, FC.FilterMatchOps.Not(b))
    entry = FC.FilterCatalogEntry("clean benzene", matcher)
    entry.SetProp("Scope", "test")
    for e2 in (FC.FilterCatalogEntry(entry.Serialize()), pickle.loads(pickle.dumps(entry))):
        assert e2.GetDescription() == "clean benzene"
        assert e2.GetProp("Scope") == "test"
        assert e2.HasFilterMatch(Chem.MolFromSmiles("c1ccccc1"))
        assert not e2.HasFilterMatch(Chem.MolFromSmiles("c1ccccc1Cl"))


def test_functional_group_hierarchy_round_trip():
    hier = FC.GetFunctionalGroupHierarchy()
    again = FC.FilterCatalog(hier.Serialize())
    mol = Chem.MolFromSmiles("CBr")
    assert again.GetNumEntries() == hier.GetNumEntries()
    names = lambda cat: [m.filterMatch.GetName() for e in cat.GetMatches(mol) for m in e.GetFilterMatches(mol)]  # noqa: E731
    assert names(again) == names(hier) != []
