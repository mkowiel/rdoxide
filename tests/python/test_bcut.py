"""BCUT2D descriptors (expectations from RDKit 2026.09.1, eigenvalues to 1e-9)."""

import pytest

from rdoxide import Chem
from rdoxide.Chem import Descriptors, rdMolDescriptors as D

EXPECTED = {
  "CC(=O)Nc1ccc(O)cc1": [
    16.254688714, 10.28019086, 1.928036059, -1.994554074, 2.125247579, -2.044011313, 5.880954785,
    -0.113735783
  ],
  "CCO": [
    16.248991329, 10.908168035, 1.340293519, -1.522358693, 1.296242638, -1.568363547, 3.835802336,
    0.318405613
  ],
  "C": [12.011, 12.011, -0.077557891, -0.077557891, 0.1441, 0.1441, 2.503, 2.503],
}


@pytest.mark.parametrize("smi", EXPECTED)
def test_bcut2d(smi):
  assert D.BCUT2D(Chem.MolFromSmiles(smi)) == pytest.approx(EXPECTED[smi], abs=1e-9)


def test_hydrogens_are_ignored():
  with_hs = Chem.AddHs(Chem.MolFromSmiles("CCO"))
  assert D.BCUT2D(with_hs) == pytest.approx(EXPECTED["CCO"], abs=1e-9)


def test_empty_molecule():
  assert D.BCUT2D(Chem.MolFromSmiles("")) == [0.0] * 8


def test_explicit_atom_properties():
  mol = Chem.MolFromSmiles("CCO")
  expected = (3.7323843332905975, 0.26828233325261325)
  assert D.BCUT2D(mol, [1.0, 2.0, 3.0]) == pytest.approx(expected, abs=1e-12)
  assert D.BCUT2D(mol, (1.0, 2.0, 3.0)) == pytest.approx(expected, abs=1e-12)
  for atom in mol.GetAtoms():
    atom.SetDoubleProp("x", atom.GetIdx() + 1.0)
  assert D.BCUT2D(mol, "x") == pytest.approx(expected, abs=1e-12)
  with pytest.raises(KeyError):
    D.BCUT2D(mol, "missing")


def test_version_and_descriptor_names():
  assert D._BCUT2D_version == "1.0.0"
  names = [n for n, _ in Descriptors._descList if n.startswith("BCUT2D")]
  assert names == [f"BCUT2D_{s}" for s in
                   ("MWHI", "MWLOW", "CHGHI", "CHGLO", "LOGPHI", "LOGPLOW", "MRHI", "MRLOW")]
  assert Descriptors.BCUT2D_MWHI(Chem.MolFromSmiles("CCO")) == pytest.approx(16.248991329, abs=1e-9)
