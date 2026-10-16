"""Ported from upstream Code/GraphMol/Wrap/rough_test.py::test1Table and
Code/GraphMol/catch_periodictable.cpp (BSD-3-Clause, RDKit contributors)."""

import pytest

from rdoxide import Chem


def feq(a, b, tol=1e-4):
    return abs(a - b) <= tol


def test_table():
    tbl = Chem.GetPeriodicTable()
    assert tbl
    assert feq(tbl.GetAtomicWeight(6), 12.011)
    assert feq(tbl.GetAtomicWeight("C"), 12.011)
    assert tbl.GetAtomicNumber("C") == 6
    assert feq(tbl.GetRvdw(6), 1.7)
    assert feq(tbl.GetRvdw("C"), 1.7)
    assert feq(tbl.GetRcovalent(6), 0.76)
    assert feq(tbl.GetRcovalent("C"), 0.76)
    assert tbl.GetDefaultValence(6) == 4
    assert tbl.GetDefaultValence("C") == 4
    assert tuple(tbl.GetValenceList(6)) == (4,)
    assert tuple(tbl.GetValenceList("C")) == (4,)
    assert tuple(tbl.GetValenceList(16)) == (2, 4, 6)
    assert tuple(tbl.GetValenceList("S")) == (2, 4, 6)
    assert tbl.GetNOuterElecs(6) == 4
    assert tbl.GetNOuterElecs("C") == 4
    assert tbl.GetMostCommonIsotope(6) == 12
    assert tbl.GetMostCommonIsotope("C") == 12
    assert tbl.GetMostCommonIsotopeMass(6) == 12.0
    assert tbl.GetMostCommonIsotopeMass("C") == 12.0
    assert tbl.GetAbundanceForIsotope(6, 12) == 98.93
    assert tbl.GetAbundanceForIsotope("C", 12) == 98.93
    assert feq(tbl.GetRb0(6), 0.77)
    assert feq(tbl.GetRb0("C"), 0.77)
    assert tbl.GetElementSymbol(6) == "C"
    assert tbl.GetElementName(6) == "Carbon"
    assert tbl.GetRow(6) == 2
    assert tbl.GetRow("C") == 2


def test_singleton_and_module():
    assert Chem.GetPeriodicTable() is Chem.GetPeriodicTable()
    assert Chem.PeriodicTable.__module__ == "rdoxide.Chem.rdchem"
    with pytest.raises(TypeError):
        Chem.PeriodicTable()


def test_errors():
    tbl = Chem.GetPeriodicTable()
    with pytest.raises(RuntimeError):
        tbl.GetAtomicNumber("Xx")
    with pytest.raises(RuntimeError):
        tbl.GetAtomicWeight(500)


def test_legacy_symbols():
    tbl = Chem.GetPeriodicTable()
    assert tbl.GetAtomicNumber("Uut") == 113
    assert tbl.GetElementSymbol(113) == "Nh"
    assert tbl.GetMaxAtomicNumber() == 118
