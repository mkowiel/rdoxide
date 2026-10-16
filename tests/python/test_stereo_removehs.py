"""RemoveHs / RemoveAllHs / RemoveHsParameters, legacy AssignStereochemistry
and the default MolFromSmiles stereo perception. Expected values checked
against RDKit 2026.09.1 (this file also passes with the official wheel);
cases adapted from upstream Code/GraphMol/Wrap/rough_test.py (BSD-3-Clause)."""

import pytest

from rdoxide import Chem


def cip_codes(m):
    return [a.GetProp("_CIPCode") if a.HasProp("_CIPCode") else None for a in m.GetAtoms()]


@pytest.mark.parametrize(
    "smi,expected",
    [
        ("C[C@H](O)F", [None, "R", None, None]),
        ("C[C@@H](O)F", [None, "S", None, None]),
        ("F[C@H](Cl)Br", [None, "R", None, None]),
        ("C[C@H](N)C(=O)O", [None, "S", None, None, None, None]),
    ],
)
def test_cip_codes_from_smiles(smi, expected):
    assert cip_codes(Chem.MolFromSmiles(smi)) == expected


def test_non_stereocenter_tag_is_cleaned():
    m = Chem.MolFromSmiles("C[C@H](C)O")
    assert m.GetAtomWithIdx(1).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
    assert not m.GetAtomWithIdx(1).HasProp("_CIPCode")
    # the H that was only there for the stereo spec is dropped (Issue 194)
    assert m.GetAtomWithIdx(1).GetNumExplicitHs() == 0


def test_ring_stereo_kept_without_cip_code():
    m = Chem.MolFromSmiles("C[C@H]1CC[C@@H](C)CC1")
    assert [a.GetIdx() for a in m.GetAtoms() if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED] == [1, 4]
    assert list(m.GetAtomWithIdx(1).GetPropsAsDict(True, True)["_ringStereoAtoms"]) == [5]


@pytest.mark.parametrize(
    "smi,stereo,atoms",
    [
        ("F/C=C/F", Chem.BondStereo.STEREOE, [0, 3]),
        ("F/C=C\\F", Chem.BondStereo.STEREOZ, [0, 3]),
        ("C/C=C(/F)Cl", Chem.BondStereo.STEREOZ, [0, 4]),
        ("C1/C=C/CCCCCC1", Chem.BondStereo.STEREOE, [0, 3]),
        ("C1/C=C/CCCC1", Chem.BondStereo.STEREONONE, []),  # 7-ring
        ("C/C=C(/C)C", Chem.BondStereo.STEREONONE, []),
    ],
)
def test_double_bond_stereo(smi, stereo, atoms):
    m = Chem.MolFromSmiles(smi)
    b = m.GetBondWithIdx(1)
    assert b.GetStereo() == stereo
    assert list(b.GetStereoAtoms()) == atoms


def test_conflicting_directions():
    m = Chem.MolFromSmiles("F/C(\\Cl)=C/F")
    assert m.GetBondWithIdx(2).GetStereo() == Chem.BondStereo.STEREONONE
    assert all(b.GetBondDir() == Chem.BondDir.NONE for b in m.GetBonds())


def test_assign_stereochemistry_flags():
    m = Chem.MolFromSmiles("CC(O)F")
    assert not m.GetAtomWithIdx(1).HasProp("_CIPCode")
    assert m.GetAtomWithIdx(1).HasProp("_ChiralityPossible")
    m.GetAtomWithIdx(1).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    Chem.AssignStereochemistry(m)  # already done: no-op without force
    assert not m.GetAtomWithIdx(1).HasProp("_CIPCode")
    Chem.AssignStereochemistry(m, force=True)
    assert m.GetAtomWithIdx(1).GetProp("_CIPCode") == "S"
    Chem.AssignStereochemistry(m, cleanIt=True, force=True, flagPossibleStereoCenters=False)
    assert not m.GetAtomWithIdx(1).HasProp("_ChiralityPossible")


def test_compute_atom_cip_ranks():
    m = Chem.MolFromSmiles("OC(F)Cl")
    assert Chem.ComputeAtomCIPRanks(m) == (1, 0, 2, 3)
    assert m.GetAtomWithIdx(3).GetUnsignedProp("_CIPRank") == 3


def test_molfromsmiles_removes_hs_and_keeps_chirality():
    m = Chem.MolFromSmiles("[H][C@](F)(Cl)Br")
    assert m.GetNumAtoms() == 4
    assert m.GetAtomWithIdx(0).GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW
    assert m.GetAtomWithIdx(0).GetProp("_CIPCode") == "S"


def test_remove_hs_defaults():
    m = Chem.MolFromSmiles("[H]C([H])([H])O[2H]", sanitize=False)
    res = Chem.RemoveHs(m)
    assert m.GetNumAtoms() == 6  # input untouched
    assert [a.GetAtomicNum() for a in res.GetAtoms()] == [6, 8, 1]
    assert res.GetAtomWithIdx(0).GetNumExplicitHs() == 0


def test_remove_hs_kept_cases():
    for smi in ("[H][H]", "[H]", "*[H]", "C[H-]C", "[H]/N=C/C"):
        m = Chem.MolFromSmiles(smi, sanitize=False)
        assert Chem.RemoveHs(m, sanitize=False).GetNumAtoms() == m.GetNumAtoms(), smi


def test_remove_hs_overloads():
    m = Chem.MolFromSmiles("[H]C([H])([H])[H]", sanitize=False)
    assert Chem.RemoveHs(m).GetNumAtoms() == 1
    assert Chem.RemoveHs(m, updateExplicitCount=True).GetAtomWithIdx(0).GetNumExplicitHs() == 4
    assert Chem.RemoveHs(m, True).GetNumAtoms() == 5  # implicitOnly: nothing marked implicit
    ps = Chem.RemoveHsParameters()
    assert ps.removeMapped and ps.removeWithWedgedBond and not ps.removeIsotopes
    assert Chem.RemoveHs(m, ps).GetNumAtoms() == 1
    assert Chem.RemoveHs(m, ps, False).GetNumAtoms() == 1
    assert Chem.RemoveHs(m, params=ps, sanitize=False).GetNumAtoms() == 1


def test_remove_hs_parameters_isotopes():
    m = Chem.MolFromSmiles("[2H]C([2H])O", sanitize=False)
    ps = Chem.RemoveHsParameters()
    assert Chem.RemoveHs(m, ps).GetNumAtoms() == 4
    ps.removeIsotopes = True
    assert Chem.RemoveHs(m, ps).GetNumAtoms() == 2
    ps.removeIsotopes = False
    ps.removeAndTrackIsotopes = True
    res = Chem.RemoveHs(m, ps)
    assert res.GetNumAtoms() == 2
    assert list(res.GetAtomWithIdx(0).GetPropsAsDict(True, True)["_isotopicHs"]) == [2, 2]


def test_remove_all_hs():
    m = Chem.MolFromSmiles("[H][H].[2H]C([H])O*[H]", sanitize=False)
    assert [a.GetAtomicNum() for a in Chem.RemoveAllHs(m).GetAtoms()] == [6, 8, 0]


def test_remove_hs_sanitize_error_propagates():
    m = Chem.MolFromSmiles("[H]C(C)(C)(C)(C)C", sanitize=False)
    with pytest.raises(Chem.AtomValenceException):
        Chem.RemoveHs(m)
    assert Chem.RemoveHs(m, sanitize=False).GetNumAtoms() == 6
