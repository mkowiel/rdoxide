"""AdjustQueryProperties; verified against RDKit 2026.09.1."""

import pytest
from rdoxide import Chem
from rdoxide.Chem import rdmolops


def adjusted(smiles, **attrs):
    params = rdmolops.AdjustQueryParameters()
    for name, value in attrs.items():
        setattr(params, name, value)
    return rdmolops.AdjustQueryProperties(Chem.MolFromSmiles(smiles), params)


def describe(atom):
    return atom.DescribeQuery().replace("\n", "/")


def test_defaults_add_degree_queries_to_ring_atoms_and_dummy_queries():
    m = rdmolops.AdjustQueryProperties(Chem.MolFromSmiles("*c1ccccc1"))
    assert describe(m.GetAtomWithIdx(0)) == "AtomNull/"
    # ring atoms get degree queries, chain atoms do not
    assert "AtomExplicitDegree 3" in describe(m.GetAtomWithIdx(1))
    assert "AtomExplicitDegree 2" in describe(m.GetAtomWithIdx(2))
    assert not Chem.MolFromSmiles("CC").GetAtomWithIdx(0).HasQuery()


def test_input_is_not_modified():
    mol = Chem.MolFromSmiles("*c1ccccc1")
    rdmolops.AdjustQueryProperties(mol)
    assert not any(a.HasQuery() for a in mol.GetAtoms())


def test_no_adjustments_changes_nothing():
    params = rdmolops.AdjustQueryParameters.NoAdjustments()
    m = rdmolops.AdjustQueryProperties(Chem.MolFromSmiles("*c1ccccc1"), params)
    assert not any(a.HasQuery() for a in m.GetAtoms())


def test_ring_count_ring_chain_and_heavy_degree():
    m = adjusted(
        "C1CC1C",
        adjustDegree=False,
        adjustRingCount=True,
        adjustRingCountFlags=rdmolops.ADJUST_IGNORENONE,
    )
    assert [describe(a) for a in m.GetAtoms()][0].startswith("AtomAnd")
    assert "AtomInNRings 1" in describe(m.GetAtomWithIdx(0))
    assert "AtomInNRings 0" in describe(m.GetAtomWithIdx(3))
    m = adjusted("C1CC1C", adjustDegree=False, adjustRingChain=True)
    assert "AtomInRing" in describe(m.GetAtomWithIdx(0))
    m = adjusted("C1CC1C", adjustDegree=False, adjustHeavyDegree=True)
    assert "AtomHeavyAtomDegree 3" in describe(m.GetAtomWithIdx(2))


def test_generic_atoms_and_bonds():
    m = adjusted("CC", adjustDegree=False, makeAtomsGeneric=True, makeBondsGeneric=True)
    assert all(describe(a) == "AtomNull/" for a in m.GetAtoms())
    assert m.GetBondWithIdx(0).HasQuery()


def test_flags_select_atoms():
    ignore_rings = adjusted(
        "C1CC1C",
        adjustDegree=True,
        adjustDegreeFlags=rdmolops.ADJUST_IGNORERINGS,
    )
    assert not ignore_rings.GetAtomWithIdx(0).HasQuery()
    assert ignore_rings.GetAtomWithIdx(3).HasQuery()


def test_single_bonds_next_to_aromatic_atoms():
    m = adjusted(
        "c1ccccc1-c1ccccc1",
        adjustDegree=False,
        adjustSingleBondsBetweenAromaticAtoms=True,
    )
    assert m.GetBondBetweenAtoms(5, 6).HasQuery()
    m = adjusted(
        "c1ccccc1C",
        adjustDegree=False,
        adjustSingleBondsToDegreeOneNeighbors=True,
    )
    assert m.GetBondBetweenAtoms(5, 6).HasQuery()


def test_conjugated_five_rings():
    m = adjusted("C1=CC=CC1", adjustDegree=False, adjustConjugatedFiveRings=True)
    assert all(b.HasQuery() for b in m.GetBonds())


def test_stereo_care_for_bonds():
    m = adjusted("F/C=C/F", adjustDegree=False, useStereoCareForBonds=True)
    assert m.GetBondWithIdx(1).GetStereo() == Chem.BondStereo.STEREONONE
    m = adjusted("F/C=C/F", adjustDegree=False)
    assert m.GetBondWithIdx(1).GetStereo() == Chem.BondStereo.STEREOE


def test_flags_enum_and_parameters_object():
    assert rdmolops.ADJUST_IGNORERINGS == 4
    assert int(rdmolops.ADJUST_IGNOREALL) == 0xFFFFFFF
    params = rdmolops.AdjustQueryParameters()
    assert params.adjustDegree and params.makeDummiesQueries and params.aromatizeIfPossible
    assert params.adjustDegreeFlags == (
        rdmolops.ADJUST_IGNOREDUMMIES | rdmolops.ADJUST_IGNORECHAINS
    )
    assert not rdmolops.AdjustQueryParameters.NoAdjustments().adjustDegree


def test_bad_arguments():
    mol = Chem.MolFromSmiles("CC")
    with pytest.raises(TypeError):
        rdmolops.AdjustQueryProperties(mol, 3)
    with pytest.raises(Exception) as err:
        rdmolops.AdjustQueryProperties(3)
    assert type(err.value).__name__ == "ArgumentError"
