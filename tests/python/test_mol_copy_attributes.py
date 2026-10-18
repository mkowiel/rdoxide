"""Python attributes of molecules follow copy.copy (shared values) and copy.deepcopy (copied)."""

import copy

from rdoxide import Chem


def test_mol_attributes_survive_copies():
  m = Chem.MolFromSmiles("CC")
  m.foo = [1]
  shallow, deep = copy.copy(m), copy.deepcopy(m)
  shallow.foo.append(4)
  assert m.foo == [1, 4]
  assert deep.foo == [1]
  assert not hasattr(Chem.RWMol(m), "foo")


def test_rwmol_attributes_survive_copies():
  m = Chem.RWMol(Chem.MolFromSmiles("CC"))
  m.foo = [1]
  shallow, deep = copy.copy(m), copy.deepcopy(m)
  shallow.foo.append(4)
  assert m.foo == [1, 4]
  assert deep.foo == [1]
