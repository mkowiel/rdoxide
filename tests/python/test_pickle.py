"""Mol pickling: Mol.ToBinary / Mol(pkl) / pickle module support and the
default pickle properties. Expected values checked against RDKit 2026.09.1
(this file also passes with the official wheel); cases adapted from upstream
Code/GraphMol/Wrap/rough_test.py (BSD-3-Clause)."""

import copy
import pickle

import pytest

from rdoxide import Chem


def test_to_binary_header():
    b = Chem.MolFromSmiles("CCO").ToBinary()
    assert b[:4] == bytes.fromhex("efbeadde")
    assert b[8:20] == bytes([16, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0])


def test_roundtrip():
    m = Chem.MolFromSmiles("C[C@H](O)c1ccccc1")
    m2 = Chem.Mol(m.ToBinary())
    assert Chem.MolToSmiles(m2) == "C[C@H](O)c1ccccc1"
    assert m2.GetRingInfo().NumRings() == 1


def test_python_pickle_keeps_type_and_dict():
    m = Chem.MolFromSmiles("CCO")
    m.foo = 3
    m2 = pickle.loads(pickle.dumps(m))
    assert type(m2) is Chem.Mol and m2.foo == 3
    rw = pickle.loads(pickle.dumps(Chem.RWMol(m)))
    assert type(rw) is Chem.RWMol
    assert m.__reduce__()[2] == ({"foo": 3},)


def test_property_flags():
    m = Chem.MolFromSmiles("CCO")
    m.SetProp("_Name", "ethanol")
    m.SetIntProp("count", 3)
    m.GetAtomWithIdx(0).SetProp("atag", "x")
    assert not Chem.Mol(m.ToBinary()).HasProp("count")
    P = Chem.PropertyPickleOptions
    m2 = Chem.Mol(m.ToBinary(P.AllProps))
    assert m2.GetIntProp("count") == 3 and m2.GetProp("_Name") == "ethanol"
    assert m2.GetAtomWithIdx(0).GetProp("atag") == "x"
    m3 = Chem.Mol(m.ToBinary(P.AllProps), P.NoProps)
    assert not m3.HasProp("count")
    assert Chem.GetDefaultPickleProperties() == 0
    try:
        Chem.SetDefaultPickleProperties(P.MolProps)
        assert pickle.loads(pickle.dumps(m)).GetIntProp("count") == 3
    finally:
        Chem.SetDefaultPickleProperties(P.NoProps)


def test_query_pickles():
    q = Chem.MolFromSmarts("[$(c-C(=O)O)]~[c;!R0]")
    q2 = Chem.Mol(q.ToBinary())
    assert Chem.MolToSmarts(q2) == Chem.MolToSmarts(q)
    m = Chem.MolFromSmiles("OC(=O)c1ccccc1")
    assert m.GetSubstructMatches(q2) == m.GetSubstructMatches(q)


def test_bad_pickle():
    with pytest.raises(RuntimeError):
        Chem.Mol(b"xx")


def test_copy_still_works():
    m = Chem.MolFromSmiles("CCO")
    assert Chem.MolToSmiles(copy.deepcopy(m)) == "CCO"
