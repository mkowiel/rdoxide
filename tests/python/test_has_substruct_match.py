"""HasSubstructMatch (existence search) agrees with GetSubstructMatch."""

import pytest
from rdoxide import Chem

SMILES = [
  "CC(=O)Nc1ccccc1", "OC(=O)c1ccccc1C(=O)O", "c1ccc2ccccc2c1", "c1ccc2c(c1)oc1ccccc12",
  "CCCCCCS", "CCN(CC)CCOC(=O)c1ccc(N)cc1", "C1CC1C(F)(F)F", "[Na+].[O-]C(=O)C",
  "C[C@H](N)C(=O)O", "ClC(Cl)(Cl)Cl", "O=C1CCCCC1", "c1ccncc1", "CC", "C",
]
SMARTS = [
  "c1ccccc1", "C(=O)[OH]", "[NX3;H2,H1;!$(NC=O)]", "[$(C=O)]N", "c1ccc2c(c1)oc1ccccc12",
  "[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#16]", "*~*~*~*~*~*~*", "C.C", "[Cl]C([Cl])[Cl]", "[C@H]",
  "[Na+]", "[!#6]~[#6]~[!#6]", "[R2]", "[#7]~[#6]~[#6]~[#6]~[#6]~[#6]~[#15]",
]


@pytest.mark.parametrize("smarts", SMARTS)
@pytest.mark.parametrize("useChirality", [False, True])
def test_has_match_agrees_with_get_match(smarts, useChirality):
  q = Chem.MolFromSmarts(smarts)
  for smi in SMILES:
    m = Chem.MolFromSmiles(smi)
    assert m.HasSubstructMatch(q, useChirality=useChirality) == bool(
      m.GetSubstructMatch(q, useChirality=useChirality)), (smi, smarts)
