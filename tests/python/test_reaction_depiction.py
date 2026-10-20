import pytest

from rdoxide.Chem import rdChemReactions as rc

SMARTS = "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]"


def _coords(rxn):
    out = []
    for m in list(rxn.GetReactants()) + list(rxn.GetProducts()):
        c = m.GetConformer()
        out.append([v for i in range(m.GetNumAtoms()) for v in (c.GetAtomPosition(i).x, c.GetAtomPosition(i).y)])
    return out


def test_templates_are_laid_out_side_by_side():
    rxn = rc.ReactionFromSmarts(SMARTS)
    rc.Compute2DCoordsForReaction(rxn)
    got = _coords(rxn)
    # reference values from RDKit 2026.09.1
    assert got[0] == pytest.approx([1.299038, -0.5, 2.598076, 0.25, 0.0, 0.25], abs=1e-5)
    assert got[1] == pytest.approx([3.598076, 0.0], abs=1e-5)
    assert got[2] == pytest.approx([5.897114, -0.5, 4.598076, 0.25, 7.196152, 0.25], abs=1e-5)


def test_spacing_and_bond_length_options():
    rxn = rc.ReactionFromSmarts(SMARTS)
    rc.Compute2DCoordsForReaction(rxn, spacing=3.0, canonOrient=False, bondLength=2.0)
    got = _coords(rxn)
    assert got[0] == pytest.approx([1.732051, 0.0, 3.464102, 1.0, 0.0, 1.0], abs=1e-5)
    assert got[2] == pytest.approx([11.196152, 0.0, 9.464102, 1.0, 12.928203, 1.0], abs=1e-5)
    # the default bond length is restored afterwards
    rxn2 = rc.ReactionFromSmarts(SMARTS)
    rc.Compute2DCoordsForReaction(rxn2)
    assert _coords(rxn2)[0][2] == pytest.approx(2.598076, abs=1e-5)
