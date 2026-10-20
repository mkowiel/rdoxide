from rdoxide import Chem


def test_complex_query_bonds_do_not_add_to_valence():
    # the implicit SMARTS bond ("single or aromatic") is a complex bond type query
    m = Chem.MolFromSmarts("[C:1](=[O:2])O")
    m.UpdatePropertyCache(False)
    assert [a.GetValence(Chem.ValenceType.EXPLICIT) for a in m.GetAtoms()] == [2, 2, 0]
    Chem.SetConjugation(m)
    Chem.SetHybridization(m)
    assert [a.GetHybridization() for a in m.GetAtoms()] == [Chem.HybridizationType.SP2] * 3
    # plain bond type queries still count
    m = Chem.MolFromSmarts("C-O")
    m.UpdatePropertyCache(False)
    assert [a.GetValence(Chem.ValenceType.EXPLICIT) for a in m.GetAtoms()] == [1, 1]
