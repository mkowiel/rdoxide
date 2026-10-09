"""rdkit.Chem.rdMMPA; expected values come from reference RDKit."""

import pytest
from rdkit import Chem
from rdkit.Chem import rdMMPA


def test_single_and_double_cuts():
    m = Chem.MolFromSmiles("c1ccccc1OC")
    assert rdMMPA.FragmentMol(m, resultsAsMols=False) == (
        ("", "CO[*:1].c1ccc([*:1])cc1"),
        ("O([*:1])[*:2]", "C[*:1].c1ccc([*:2])cc1"),
        ("", "C[*:1].c1ccc(O[*:1])cc1"),
    )


def test_min_max_cuts_overload():
    m = Chem.MolFromSmiles("CCOCCN")
    res = rdMMPA.FragmentMol(m, minCuts=1, maxCuts=2, maxCutBonds=10, resultsAsMols=False)
    assert len(res) == 15
    assert res[0] == ("", "C[*:1].NCCOC[*:1]")


def test_bonds_to_cut():
    m = Chem.MolFromSmiles("CCOCCN")
    assert rdMMPA.FragmentMol(m, [1, 2], minCuts=2, maxCuts=2, resultsAsMols=False) == (
        ("O([*:1])[*:2]", "CC[*:1].NCC[*:2]"),
    )


def test_results_as_mols():
    core, chains = rdMMPA.FragmentMol(Chem.MolFromSmiles("c1ccccc1OC"))[1]
    assert Chem.MolToSmiles(core) == "O([*:1])[*:2]"
    assert Chem.MolToSmiles(chains) == "C[*:1].c1ccc([*:2])cc1"
    assert rdMMPA.FragmentMol(Chem.MolFromSmiles("O")) == ()


def test_too_many_cut_bonds_gives_nothing():
    m = Chem.MolFromSmiles("CCOCCN")
    assert rdMMPA.FragmentMol(m, maxCutBonds=1) == ()


def test_bad_cut_counts():
    m = Chem.MolFromSmiles("CCO")
    with pytest.raises(ValueError):
        rdMMPA.FragmentMol(m, minCuts=0, maxCuts=2, maxCutBonds=5)
    with pytest.raises(ValueError):
        rdMMPA.FragmentMol(m, minCuts=3, maxCuts=2, maxCutBonds=5)
