"""Wildman-Crippen logP/MR and AddHs (expectations from RDKit 2026.09.1)."""

import pytest

from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors
from rdkit.Geometry import Point3D


def test_calc_crippen_descriptors():
    logp, mr = rdMolDescriptors.CalcCrippenDescriptors(Chem.MolFromSmiles("c1ccccc1O"))
    assert logp == pytest.approx(1.3922)
    assert mr == pytest.approx(28.1068)
    # with includeHs=True the implicit Hs are counted, with includeHs=False they are not
    assert rdMolDescriptors.CalcCrippenDescriptors(Chem.MolFromSmiles("CCO")) == pytest.approx((-0.0014, 12.7598))
    assert rdMolDescriptors.CalcCrippenDescriptors(Chem.MolFromSmiles("CCO"), False) == pytest.approx((-0.3487, 6.0798))
    assert rdMolDescriptors.CalcCrippenDescriptors(Chem.MolFromSmiles("C[N+](C)(C)C.[Cl-]")) == pytest.approx((-2.6736, 23.9564))


def test_morgan_wrappers_crippen():
    assert Crippen.MolLogP(Chem.MolFromSmiles("OC(=O)c1ccccc1")) == pytest.approx(1.3848)
    assert Crippen.MolMR(Chem.MolFromSmiles("OC(=O)c1ccccc1")) == pytest.approx(33.4013)


def test_crippen_atom_contribs():
    m = Chem.MolFromSmiles("CCO")
    types, labels = [0] * 3, [""] * 3
    contribs = rdMolDescriptors._CalcCrippenContribs(m, False, types, labels)
    assert [tuple(c) for c in contribs] == [
        pytest.approx((0.1441, 2.503)),
        pytest.approx((-0.2035, 2.753)),
        pytest.approx((-0.2893, 0.8238)),
    ]
    assert types == [1, 6, 69]
    assert labels == ["C1", "C3", "O2"]
    # without the optional lists the contributions are the same
    assert rdMolDescriptors._CalcCrippenContribs(m)[1] == contribs[1]


def test_crippen_atom_contribs_cache():
    m = Chem.MolFromSmiles("CCO")
    rdMolDescriptors.CalcCrippenDescriptors(m, False)
    types, labels = [0] * 3, [""] * 3
    contribs = rdMolDescriptors._CalcCrippenContribs(m, False, types, labels)
    # served from the cache stored on the molecule: atom types are not recomputed
    assert [c[0] for c in contribs] == [pytest.approx(0.1441), pytest.approx(-0.2035), pytest.approx(-0.2893)]
    assert types == [0, 0, 0]
    assert labels == ["", "", ""]
    types, labels = [0] * 3, [""] * 3
    rdMolDescriptors._CalcCrippenContribs(m, True, types, labels)
    assert types == [1, 6, 69]


def test_crippen_atom_contribs_bad_lists():
    m = Chem.MolFromSmiles("CCO")
    with pytest.raises(ValueError):
        rdMolDescriptors._CalcCrippenContribs(m, False, [0, 0])
    with pytest.raises(ValueError):
        rdMolDescriptors._CalcCrippenContribs(m, False, [0] * 3, [""] * 2)


def test_addhs_defaults_and_options():
    m = Chem.MolFromSmiles("CCO")
    assert Chem.AddHs(m).GetNumAtoms() == 9
    explicit = Chem.AddHs(m, explicitOnly=True)
    assert explicit.GetNumAtoms() == 3
    assert Chem.MolToSmiles(explicit) == "CCO"
    assert [a.GetSymbol() for a in Chem.AddHs(m, onlyOnAtoms=[0]).GetAtoms()] == ["C", "C", "O", "H", "H", "H"]
    # added implicit Hs are flagged, the input molecule is not changed
    flags = [a.GetPropsAsDict(True, True).get("isImplicit") for a in Chem.AddHs(m).GetAtoms()]
    assert flags == [None, None, None] + [1] * 6
    assert m.GetNumAtoms() == 3


def test_addhs_parameters():
    params = Chem.AddHsParameters()
    params.explicitOnly = True
    assert Chem.AddHs(Chem.MolFromSmiles("CCO"), params).GetNumAtoms() == 3
    assert Chem.AddHs(Chem.MolFromSmiles("CCO"), params, onlyOnAtoms=[0]).GetNumAtoms() == 3


def test_addhs_coordinates():
    m = Chem.MolFromSmiles("CCO")
    conf = Chem.Conformer(3)
    for i, xyz in enumerate([(0.0, 0.0, 0.0), (1.5, 0.2, 0.1), (2.2, 1.4, -0.3)]):
        conf.SetAtomPosition(i, Point3D(*xyz))
    m.AddConformer(conf, assignId=True)
    mh = Chem.AddHs(m, addCoords=True)
    got = [(p.x, p.y, p.z) for p in (mh.GetConformer().GetAtomPosition(i) for i in range(mh.GetNumAtoms()))]
    expected = [
        (0.0, 0.0, 0.0),
        (1.5, 0.2, 0.1),
        (2.2, 1.4, -0.3),
        (-0.430439, -0.057392, 1.010657),
        (-0.447471, 0.84643, -0.541595),
        (-0.210067, -0.934102, -0.541595),
        (1.746079, 0.055869, 1.16239),
        (1.84817, -0.502318, -0.671704),
        (2.489911, 1.89699, 0.505638),
    ]
    assert len(got) == len(expected)
    for g, e in zip(got, expected):
        assert g == pytest.approx(e, abs=1e-6)
