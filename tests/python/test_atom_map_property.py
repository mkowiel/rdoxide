"""A map number set as a string property counts as the atom map number."""

from rdoxide import Chem


def test_string_map_number_property():
  m = Chem.MolFromSmiles("CC*")
  m.GetAtomWithIdx(2).SetProp("molAtomMapNumber", "30")
  assert m.GetAtomWithIdx(2).GetAtomMapNum() == 30
  assert Chem.MolToSmiles(m) == "CC[*:30]"
