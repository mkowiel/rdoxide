"""rdMolDescriptors.CalcOxidationNumbers (values taken from RDKit 2026.03.6)."""

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as md

CASES = [
    ("CO", [-2, -2]),
    ("CC(=O)O", [-3, 3, -2, -2]),
    ("c1ccccc1N(=O)=O", [-1, -1, -1, -1, -1, 1, 3, -2, -2]),
    ("[Na+].[Cl-]", [1, -1]),
    ("C[S](C)=O", [-2, 0, -2, -2]),
    ("c1ccncc1", [-1, -1, 1, -3, 0, -1]),
    ("OO", [-1, -1]),
    ("[O-][N+](=O)C", [-2, 3, -2, -2]),
]


@pytest.mark.parametrize("smiles,expected", CASES)
def test_oxidation_numbers(smiles, expected):
    m = Chem.MolFromSmiles(smiles)
    md.CalcOxidationNumbers(m)
    assert [a.GetIntProp("OxidationNumber") for a in m.GetAtoms()] == expected
    # the input stays aromatic: the calculation uses a kekulized copy
    assert Chem.MolToSmiles(m) == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
