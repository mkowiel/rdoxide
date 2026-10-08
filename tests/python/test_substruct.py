"""HasSubstructMatch / GetSubstructMatch / GetSubstructMatches and
SubstructMatchParameters. Expected values checked against RDKit 2026.09.1
(this file also passes with the official wheel); cases adapted from upstream
Code/GraphMol/Wrap/rough_test.py and Substruct/catch_tests.cpp (BSD-3-Clause)."""

import pytest

from rdkit import Chem


def test_basic_matching():
    m = Chem.MolFromSmiles("OC(=O)c1ccccc1C(=O)O")
    q = Chem.MolFromSmarts("C(=O)[OH]")
    assert m.HasSubstructMatch(q)
    assert m.GetSubstructMatch(q) == (1, 2, 0)
    assert m.GetSubstructMatches(q) == ((1, 2, 0), (9, 10, 11))
    assert not m.HasSubstructMatch(Chem.MolFromSmarts("[#7]"))
    assert m.GetSubstructMatch(Chem.MolFromSmarts("[#7]")) == ()
    assert m.GetSubstructMatches(Chem.MolFromSmarts("[#7]")) == ()


def test_uniquify_and_max_matches():
    m = Chem.MolFromSmiles("c1ccccc1")
    q = Chem.MolFromSmiles("c1ccccc1")
    assert len(m.GetSubstructMatches(q)) == 1
    allm = m.GetSubstructMatches(q, uniquify=False)
    assert len(allm) == 12
    assert allm[:2] == ((0, 1, 2, 3, 4, 5), (0, 5, 4, 3, 2, 1))
    assert len(m.GetSubstructMatches(q, uniquify=False, maxMatches=3)) == 3


def test_recursive_queries():
    m = Chem.MolFromSmiles("OC(=O)c1ccccc1N")
    assert m.GetSubstructMatches(Chem.MolFromSmarts("[$(c-C(=O)O)]")) == ((3,),)
    assert m.GetSubstructMatches(Chem.MolFromSmarts("[c;!$(c-[#7]);!$(c-C)]")) == ((4,), (5,), (6,), (7,))
    ps = Chem.SubstructMatchParameters()
    ps.recursionPossible = False
    assert m.GetSubstructMatches(Chem.MolFromSmarts("[$(c-C(=O)O)]"), ps) == ()


def test_chirality():
    m = Chem.MolFromSmiles("C[C@H](N)C(=O)O")
    q_same = Chem.MolFromSmiles("C[C@H](N)C=O")
    q_other = Chem.MolFromSmiles("C[C@@H](N)C=O")
    assert m.HasSubstructMatch(q_other)
    assert m.HasSubstructMatch(q_same, useChirality=True)
    assert not m.HasSubstructMatch(q_other, useChirality=True)
    ps = Chem.SubstructMatchParameters()
    ps.useChirality = True
    assert not m.GetSubstructMatches(q_other, ps)
    m2 = Chem.MolFromSmiles("CC(N)C(=O)O")
    assert not m2.HasSubstructMatch(q_same, useChirality=True)
    ps.specifiedStereoQueryMatchesUnspecified = True
    assert m2.HasSubstructMatch(q_same, ps)


def test_double_bond_stereo():
    m = Chem.MolFromSmiles("F/C=C/F")
    assert m.HasSubstructMatch(Chem.MolFromSmiles("F/C=C/F"), useChirality=True)
    assert not m.HasSubstructMatch(Chem.MolFromSmiles("F/C=C\\F"), useChirality=True)
    assert m.HasSubstructMatch(Chem.MolFromSmiles("F/C=C\\F"))


def test_params_defaults():
    ps = Chem.SubstructMatchParameters()
    assert ps.uniquify and ps.recursionPossible and not ps.useChirality
    assert ps.maxMatches == 1000 and ps.maxRecursiveMatches == 1000 and ps.numThreads == 1
    assert not ps.aromaticMatchesConjugated and not ps.useQueryQueryMatches


def test_aromatic_matches_conjugated():
    m = Chem.MolFromSmiles("C1=CC=CC=C1")
    Chem.Kekulize(m, clearAromaticFlags=True)
    q = Chem.MolFromSmiles("c1ccccc1")
    assert not m.HasSubstructMatch(q)
    ps = Chem.SubstructMatchParameters()
    ps.aromaticMatchesConjugated = True
    assert m.HasSubstructMatch(q, ps)


def test_query_query_matches():
    target = Chem.MolFromSmarts("[C,N]")
    q = Chem.MolFromSmarts("[C]")
    assert not target.HasSubstructMatch(q)
    assert target.HasSubstructMatch(q, useQueryQueryMatches=True)


def test_atom_properties():
    m = Chem.MolFromSmiles("CCO")
    q = Chem.MolFromSmiles("CO")
    m.GetAtomWithIdx(1).SetProp("tag", "a")
    q.GetAtomWithIdx(0).SetProp("tag", "a")
    ps = Chem.SubstructMatchParameters()
    ps.atomProperties = ["tag"]
    assert m.GetSubstructMatches(q, ps) == ((1, 2),)
