"""SubstructMatchParameters callbacks; expectations hold for RDKit as well."""

import pytest
from rdoxide import Chem


def params(**kw):
    ps = Chem.SubstructMatchParameters()
    for k, v in kw.items():
        setattr(ps, k, v)
    return ps


def test_atom_check_sees_query_then_molecule_atoms():
    m = Chem.MolFromSmiles("CCOCC")
    q = Chem.MolFromSmiles("CO")
    seen = []

    def check(qa, ma):
        seen.append(
            (
                qa.GetOwningMol().GetNumAtoms() == q.GetNumAtoms(),
                ma.GetOwningMol().GetNumAtoms() == m.GetNumAtoms(),
            )
        )
        return True

    ps = params()
    ps.setExtraAtomCheckFunc(check)
    assert m.GetSubstructMatches(q, ps) == ((1, 2), (3, 2))
    assert seen and all(s == (True, True) for s in seen)


def test_atom_check_override_and_bond_check():
    m = Chem.MolFromSmiles("CCOCC")
    q = Chem.MolFromSmiles("CO")
    ps = params(extraAtomCheckOverridesDefaultCheck=True)
    ps.setExtraAtomCheckFunc(lambda qa, ma: qa.GetIdx() == 0 and ma.GetIdx() == 0 or qa.GetIdx() == 1)
    assert m.GetSubstructMatches(q, ps) == ((0, 1),)

    m = Chem.MolFromSmiles("CC=CC")
    q = Chem.MolFromSmiles("CC")
    ps = params()
    ps.setExtraBondCheckFunc(lambda qb, mb: mb.GetBondTypeAsDouble() == 1.0)
    assert m.GetSubstructMatches(q, ps) == ((0, 1), (2, 3))


def test_final_check():
    m = Chem.MolFromSmiles("CCOCC")
    p = Chem.MolFromSmiles("CCO")
    ps = params()
    ps.setExtraFinalCheck(lambda mol, vect: False)
    assert m.GetSubstructMatches(p, ps) == ()
    assert not m.HasSubstructMatch(p, ps)
    ps.setExtraFinalCheck(lambda mol, vect: sum(vect) > 5)
    assert m.GetSubstructMatches(p, ps) == ((4, 3, 2),)


def test_callback_exception_propagates():
    m = Chem.MolFromSmiles("CCO")
    ps = params()

    def boom(qa, ma):
        raise ValueError("boom")

    ps.setExtraAtomCheckFunc(boom)
    with pytest.raises(ValueError, match="boom"):
        m.GetSubstructMatches(Chem.MolFromSmiles("CO"), ps)


def test_coords_matcher():
    m = Chem.MolFromMolBlock(Chem.MolToMolBlock(Chem.MolFromSmiles("CC")))
    assert m is not None
    m = Chem.MolFromSmiles("CCC |(0,0,0;1,0,0;2,0,0)|")
    q = Chem.MolFromSmiles("CC |(1,0,0.05;2,0,0)|")
    ps = params()
    ps.setExtraAtomCheckFunc(Chem.AtomCoordsMatcher(tol=0.1))
    assert m.GetSubstructMatches(q, ps) == ((1, 2),)
    ps.setExtraAtomCheckFunc(Chem.AtomCoordsMatcher(tol=0.01))
    assert m.GetSubstructMatches(q, ps) == ()
    # no conformers: nothing matches
    assert Chem.MolFromSmiles("CCC").GetSubstructMatches(q, ps) == ()
