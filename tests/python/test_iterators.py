"""Atom and bond iterators fail when their molecule is edited meanwhile."""

import pytest
from rdoxide import Chem


def test_iterators_fail_when_the_molecule_is_edited():
  mol = Chem.RWMol(Chem.MolFromSmiles("*C*"))
  atoms = iter(mol.GetAtoms())
  mol.RemoveAtom(next(atoms).GetIdx())
  with pytest.raises(RuntimeError):
    next(atoms)

  mol = Chem.RWMol(Chem.MolFromSmiles("*C*"))
  bonds = iter(mol.GetBonds())
  bond = next(bonds)
  mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
  with pytest.raises(RuntimeError):
    next(bonds)
