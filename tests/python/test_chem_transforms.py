"""ReplaceSubstructs, ReplaceSidechains, ReplaceCore and MurckoDecompose (expectations from RDKit 2026.09.1)."""

import json
import os

import pytest

from rdkit import Chem
from rdkit.Chem import rdmolops

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "mols3d.json")) as f:
    BLOCKS = json.load(f)

M = Chem.MolFromSmiles
Q = Chem.MolFromSmarts
S = Chem.MolToSmiles


def smiles(mols):
    return [S(m) for m in mols]


def test_replace_substructs():
    assert smiles(rdmolops.ReplaceSubstructs(M("COCCOC"), Q("O[CH3]"), M("NC"))) == ["CNCCOC", "CNCCOC"]
    assert smiles(rdmolops.ReplaceSubstructs(M("CCOC"), Q("[CH3]O"), M("NC"))) == ["CC.CN"]
    assert smiles(rdmolops.ReplaceSubstructs(M("COCCOC"), Q("O[CH3]"), M("CN"), True, 1)) == ["CNCCNC"]
    # no match: a copy of the molecule
    assert smiles(rdmolops.ReplaceSubstructs(M("CCCC"), Q("N"), M("O"))) == ["CCCC"]


def test_replace_substructs_bad_connection_point():
    with pytest.raises(RuntimeError, match="bad replacementConnectionPoint"):
        rdmolops.ReplaceSubstructs(M("CCOC"), Q("O"), M("N"), False, 5)


def test_replace_sidechains():
    assert S(rdmolops.ReplaceSidechains(M("CCC1CCC1"), M("C1CCC1"))) == "[1*]C1CCC1"
    assert S(rdmolops.ReplaceSidechains(M("C1CC2C1CCC2"), M("C1CCC1"))) == "[1*]C1CCC1[2*]"
    assert rdmolops.ReplaceSidechains(M("CCC1CC1"), M("C1CCC1")) is None


def test_replace_core():
    assert S(rdmolops.ReplaceCore(M("CCC1CCC1"), M("C1CCC1"))) == "[1*]CC"
    assert S(rdmolops.ReplaceCore(M("C1CC2C1CCC2"), M("C1CCC1"))) == "[1*]CCC[2*]"
    assert S(rdmolops.ReplaceCore(M("CCN1CCC1"), M("C1CCN1"), labelByIndex=True)) == "[3*]CC"
    assert rdmolops.ReplaceCore(M("CCC1CC1"), M("C1CCC1")) is None
    assert rdmolops.ReplaceCore(M("CC1CC1"), M("C1CC1"), requireDummyMatch=True) is None


def test_replace_core_dummies_and_matches():
    mol, core = M("C1CNCC1"), Q("[*]N[*]")
    assert S(rdmolops.ReplaceCore(mol, core)) == "[1*]CC[2*]"
    assert S(rdmolops.ReplaceCore(mol, core, replaceDummies=False)) == "[1*]CCCC[2*]"
    # molecule atom indices in core order, and (core, molecule) pairs
    assert S(rdmolops.ReplaceCore(mol, core, mol.GetSubstructMatch(core))) == "[1*]CC[2*]"
    assert S(rdmolops.ReplaceCore(mol, core, ((0, 1), (1, 2), (2, 3)))) == "[1*]CC[2*]"


def test_replace_core_bad_matches():
    mol, core = M("C1CNCC1"), Q("[*]N[*]")
    with pytest.raises(ValueError, match="out of bounds of the core molecule"):
        rdmolops.ReplaceCore(mol, core, ((0, 1), (1, 2), (5, 3)))
    with pytest.raises(ValueError, match="out of bounds of the target molecule"):
        rdmolops.ReplaceCore(mol, core, ((0, 1), (1, 2), (2, 30)))
    with pytest.raises(ValueError, match="supplied core requires 3 entries"):
        rdmolops.ReplaceCore(mol, core, (1, 2))
    with pytest.raises(ValueError, match="Input not a vector"):
        rdmolops.ReplaceCore(mol, core, ((0, 1, 2), (1, 2, 3), (2, 3, 4)))
    with pytest.raises(TypeError, match="did not match C\\+\\+ signature"):
        rdmolops.ReplaceCore(mol, "x")


def test_replace_core_keeps_conformer_with_dummy_positions():
    mol = Chem.MolFromMolBlock(BLOCKS["ethanol"], removeHs=False)
    res = rdmolops.ReplaceCore(mol, Q("CO"), labelByIndex=True)
    assert res.GetNumConformers() == 1
    assert res.GetNumAtoms() == 11  # 7 sidechain atoms (C and 6 H) + 4 dummies
    dummies = [a for a in res.GetAtoms() if a.GetAtomicNum() == 0]
    assert len(dummies) == 4


@pytest.mark.parametrize(
    "smi, scaffold",
    [
        ("c1ccccc1CCc1ccncc1", "c1ccc(CCc2ccncc2)cc1"),
        ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "c1ccccc1"),
        ("CCO", ""),
        ("C[C@H](F)C1CCC(=O)CC1", "O=C1CCCCC1"),
        ("O=C1CCCCC1C(=O)C1CC1", "O=C1CCCCC1C(=O)C1CC1"),
        ("c1ccc2[nH]ccc2c1CCN", "c1ccc2[nH]ccc2c1"),
    ],
)
def test_murcko_decompose(smi, scaffold):
    assert S(rdmolops.MurckoDecompose(M(smi))) == scaffold


