"""MergeQueryHs / HasQueryHs / MolFromSmarts(mergeHs=True); verified against RDKit 2026.09.1."""

import pytest
from rdkit import Chem


def merged(smarts, *args):
    return Chem.MolToSmarts(Chem.MergeQueryHs(Chem.MolFromSmarts(smarts), *args))


def test_explicit_hs_become_negated_h_counts():
    assert merged("C[H]") == "[C&!H0]"
    assert merged("C([H])[H]") == "[C&!H0&!H1]"
    assert merged("[CH3][H]") == "[C&H3&!H0]"
    assert merged("[#6][#1]") == "[#6&!H0]"


def test_hydrogens_without_a_heavy_neighbor_stay():
    assert Chem.MergeQueryHs(Chem.MolFromSmarts("[H][H]")).GetNumAtoms() == 2
    assert Chem.MergeQueryHs(Chem.MolFromSmarts("[H]")).GetNumAtoms() == 1


def test_mapped_and_isotopic_hs():
    q = Chem.MolFromSmarts("[H:1]C")
    assert Chem.MergeQueryHs(q).GetNumAtoms() == 1
    assert Chem.MergeQueryHs(q, mergeUnmappedOnly=True).GetNumAtoms() == 2
    iso = Chem.MolFromSmarts("[2H]C")
    assert Chem.MergeQueryHs(iso).GetNumAtoms() == 2
    assert Chem.MergeQueryHs(iso, mergeIsotopes=True).GetNumAtoms() == 1


def test_recursive_queries_are_merged():
    q = Chem.MergeQueryHs(Chem.MolFromSmarts("[$(C[H])]C"))
    assert q.GetNumAtoms() == 2
    assert "H" in Chem.MolToSmarts(q)


def test_input_is_not_modified():
    q = Chem.MolFromSmarts("C[H]")
    Chem.MergeQueryHs(q)
    assert q.GetNumAtoms() == 2


def test_has_query_hs():
    assert Chem.HasQueryHs(Chem.MolFromSmarts("C[H]")) == (True, False)
    assert Chem.HasQueryHs(Chem.MolFromSmarts("CC")) == (False, False)
    assert Chem.HasQueryHs(Chem.MolFromSmarts("[#1,#6]C")) == (True, True)
    assert Chem.HasQueryHs(Chem.MolFromSmarts("[$(C[H])]C")) == (True, False)


def test_smarts_parser_merges_on_request():
    assert Chem.MolFromSmarts("C[H]", mergeHs=True).GetNumAtoms() == 1
    params = Chem.SmartsParserParams()
    params.mergeHs = True
    assert Chem.MolToSmarts(Chem.MolFromSmarts("C([H])[H]", params)) == "[C&!H0&!H1]"
    assert Chem.MolFromSmarts("C[H]").GetNumAtoms() == 2
