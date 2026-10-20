import json
from pathlib import Path

import pytest
from rdoxide import Chem
from rdoxide.Chem import rdMolDescriptors

# values computed with RDKit 2026.09.1 (`GetAtomFeatures`)
REFERENCE = json.loads((Path(__file__).parent / "data_atom_features.json").read_text())


@pytest.mark.parametrize("smiles", list(REFERENCE))
def test_matches_rdkit(smiles):
    m = Chem.MolFromSmiles(smiles)
    for add_chiral, expected in zip((False, True), REFERENCE[smiles]):
        for i, row in enumerate(expected):
            got = rdMolDescriptors.GetAtomFeatures(m, i, add_chiral)
            assert len(got) == (52 if add_chiral else 49)
            assert got == pytest.approx(row)


def test_bad_atom_index():
    m = Chem.MolFromSmiles("CC")
    with pytest.raises(RuntimeError):
        rdMolDescriptors.GetAtomFeatures(m, 5)
    assert rdMolDescriptors._GetAtomFeatures_version == "1.0.0"
