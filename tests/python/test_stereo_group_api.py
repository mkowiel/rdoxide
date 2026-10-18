"""StereoGroup objects follow their molecule; CreateStereoGroup validates its arguments."""

import pytest
from rdoxide import Chem


def test_write_id_is_shared_with_the_molecule():
  m = Chem.MolFromSmiles("C[C@H](O)Cl |o5:1|")
  group = m.GetStereoGroups()[0]
  assert (group.GetReadId(), group.GetWriteId()) == (5, 0)
  group.SetWriteId(7)
  assert Chem.MolToCXSmiles(m) == "C[C@H](O)Cl |o7:1|"
  copy = Chem.RWMol(m)
  Chem.ForwardStereoGroupIds(m)
  assert group.GetWriteId() == 5
  assert Chem.MolToCXSmiles(m) == "C[C@H](O)Cl |o5:1|"
  assert copy.GetStereoGroups()[0].GetWriteId() == 7


def test_create_stereo_group_checks_its_arguments():
  m = Chem.RWMol(Chem.MolFromSmiles("F[C@@H](Br)[C@H](F)Cl"))
  kind = Chem.StereoGroupType.STEREO_OR
  with pytest.raises(ValueError):
    Chem.CreateStereoGroup(kind, m, [], [])
  with pytest.raises(ValueError):
    Chem.CreateStereoGroup(kind, m, [100])
  with pytest.raises(TypeError):
    Chem.CreateStereoGroup(kind, None, [1])
  with pytest.raises(TypeError):
    Chem.CreateStereoGroup(kind, m, [1, "text"])
