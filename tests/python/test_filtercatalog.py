import pytest

from rdoxide import Chem
from rdoxide.Chem import FilterCatalog
from rdoxide.Chem.FilterCatalog import FilterCatalogParams, FilterMatchOps


def test_builtin_catalog_sizes():
    for cat, n in [(FilterCatalogParams.FilterCatalogs.PAINS_A, 16),
                   (FilterCatalogParams.FilterCatalogs.BRENK, 105),
                   (FilterCatalogParams.FilterCatalogs.ZINC, 50)]:
        assert FilterCatalog.FilterCatalog(cat).GetNumEntries() == n


def test_entry_properties_and_matches():
    params = FilterCatalogParams()
    assert params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    catalog = FilterCatalog.FilterCatalog(params)
    mol = Chem.MolFromSmiles('O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2')
    entry = catalog.GetFirstMatch(mol)
    assert entry.GetDescription() == 'hzone_phenol_A(479)'
    assert entry.GetProp('Scope') == 'PAINS filters (family A)'
    assert catalog.HasMatch(mol)
    assert len(catalog.GetMatches(mol)) >= 1
    with pytest.raises(KeyError):
        entry.GetProp('nope')


def test_python_matcher_and_operators():
    class Always(FilterCatalog.FilterMatcher):
        def IsValid(self):
            return True

        def HasMatch(self, mol):
            return True

    always = Always("always")
    assert always.GetName() == "always"
    not_ = FilterMatchOps.Not(always)
    assert not not_.HasMatch(Chem.MolFromSmiles("C"))
    assert str(not_) == "(Not always)"


def test_run_filter_catalog_threads():
    catalog = FilterCatalog.FilterCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    smiles = ["CC(=O)Cl", "CCCC", "bad smiles"]
    one = FilterCatalog.RunFilterCatalog(catalog, smiles, numThreads=1)
    many = FilterCatalog.RunFilterCatalog(catalog, smiles, numThreads=3)
    assert [[e.GetDescription() for e in r] for r in one] == \
        [[e.GetDescription() for e in r] for r in many]
    assert one[2][0].GetDescription() == "no valid RDKit molecule"


def test_catalog_not_picklable():
    import pickle
    with pytest.raises(RuntimeError):
        pickle.dumps(FilterCatalog.FilterCatalog())
