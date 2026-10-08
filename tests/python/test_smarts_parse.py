"""MolFromSmarts / AtomFromSmarts / BondFromSmarts / SmartsParserParams and
query descriptions. Expected values checked against RDKit 2026.09.1 (this
file also passes with the official wheel); cases adapted from upstream
Code/GraphMol/Wrap/rough_test.py (BSD-3-Clause)."""

import pytest

from rdkit import Chem


def test_basic_smarts():
    m = Chem.MolFromSmarts("[C;H2]-,=[N,O]")
    assert m.GetNumAtoms() == 2 and m.GetNumBonds() == 1
    a0, a1 = m.GetAtoms()
    assert a0.HasQuery() and a1.HasQuery()
    assert a0.DescribeQuery() == "AtomAnd\n  AtomType 6 = val\n  AtomHCount 2 = val\n"
    assert a1.DescribeQuery() == "AtomOr\n  AtomType 7 = val\n  AtomType 8 = val\n"
    assert m.GetBondWithIdx(0).DescribeQuery() == "BondOr\n  BondOrder 1 = val\n  BondOrder 2 = val\n"


def test_unspecified_bond_is_single_or_aromatic():
    m = Chem.MolFromSmarts("cC")
    b = m.GetBondWithIdx(0)
    assert b.DescribeQuery() == "SingleOrAromaticBond 1 = val\n"
    assert b.GetBondType() == Chem.BondType.SINGLE


def test_recursive_and_negation():
    a = Chem.MolFromSmarts("[!$(C=O)]").GetAtomWithIdx(0)
    assert a.DescribeQuery() == "RecursiveStructure val not in )\n"
    a = Chem.MolFromSmarts("[!#6;!#1]").GetAtomWithIdx(0)
    assert a.DescribeQuery() == "AtomAnd\n  AtomAtomicNum 6 != val\n  AtomAtomicNum 1 != val\n"


def test_ranges_and_ring_queries():
    describe = lambda s: Chem.MolFromSmarts(s).GetAtomWithIdx(0).DescribeQuery()
    assert describe("[D{2-3}]") == "range_AtomExplicitDegree 2 val 3\n"
    assert describe("[D{-2}]") == "greater_AtomExplicitDegree 2 >= \n"
    assert describe("[k5]") == "AtomRingSize 5 = val\n"
    assert describe("[R]") == "AtomInNRings -1 = val\n"
    assert describe("[r]") == "AtomInRing 1 = val\n"


def test_query_type_labels():
    assert Chem.MolFromSmarts("*").GetAtomWithIdx(0).GetQueryType() == ""
    assert Chem.MolFromSmarts("[C]").GetAtomWithIdx(0).GetQueryType() == ""


def test_hydrogens_and_charges():
    m = Chem.MolFromSmarts("[2H+:3]")
    a = m.GetAtomWithIdx(0)
    assert (a.GetAtomicNum(), a.GetIsotope(), a.GetFormalCharge(), a.GetAtomMapNum()) == (1, 2, 1, 3)
    a = Chem.MolFromSmarts("[CH2]").GetAtomWithIdx(0)
    assert a.GetNumExplicitHs() == 2 and a.GetNoImplicit()


def test_atom_and_bond_from_smarts():
    a = Chem.AtomFromSmarts("[C,N]")
    assert a.DescribeQuery() == "AtomOr\n  AtomType 6 = val\n  AtomType 7 = val\n"
    b = Chem.BondFromSmarts("!@")
    assert b.DescribeQuery() == "BondInRing 1 != val\n"
    assert Chem.AtomFromSmarts("") is None


def test_failures_return_none():
    for s in ["[C", "C)", "[C;H{1-2}]", "C1CC", "[$(C]"]:
        assert Chem.MolFromSmarts(s) is None, s


def test_params_and_replacements():
    ps = Chem.SmartsParserParams()
    assert ps.parseName and ps.allowCXSMILES and not ps.mergeHs
    assert Chem.MolFromSmarts("CC", ps).GetNumAtoms() == 2
    m = Chem.MolFromSmarts("C{X}", replacements={"{X}": "[O,N]"})
    assert m.GetNumAtoms() == 2
