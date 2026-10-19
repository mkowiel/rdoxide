"""Atom and bond property lists written on the molecule."""

from rdoxide import Chem


def test_atom_and_bond_lists():
  m = Chem.MolFromSmiles("CCO")
  m.GetAtomWithIdx(0).SetIntProp("n", 3)
  m.GetAtomWithIdx(2).SetIntProp("n", 5)
  Chem.CreateAtomIntPropertyList(m, "n")
  assert m.GetProp("atom.iprop.n") == "3 n/a 5"
  Chem.CreateAtomIntPropertyList(m, "n", missingValueMarker="?")
  assert m.GetProp("atom.iprop.n") == "[?] 3 ? 5"

  m.GetBondWithIdx(1).SetProp("tag", "x")
  Chem.CreateBondStringPropertyList(m, "tag")
  assert m.GetProp("bond.prop.tag") == "n/a x"


def test_long_lists_wrap():
  m = Chem.MolFromSmiles("C" * 12)
  for atom in m.GetAtoms():
    atom.SetDoubleProp("q", 0.5)
  Chem.CreateAtomDoublePropertyList(m, "q", lineSize=20)
  lines = m.GetProp("atom.dprop.q").split("\n")
  assert len(lines) > 1 and all(len(line) < 20 for line in lines)
  assert " ".join(lines).split() == ["0.5"] * 12
