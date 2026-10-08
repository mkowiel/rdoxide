"""Property registry of rdMolDescriptors (behaviour checked against RDKit 2026.09.1)."""

import pytest

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors as D

PARACETAMOL = "CC(=O)Nc1ccc(O)cc1"


@pytest.fixture
def mol():
  return Chem.MolFromSmiles(PARACETAMOL)


def test_available_properties_keep_upstream_order():
  names = list(D.Properties.GetAvailableProperties())
  assert len(names) >= 43
  assert names[:4] == ["exactmw", "amw", "lipinskiHBA", "lipinskiHBD"]
  assert names[-3:] == ["kappa2", "kappa3", "Phi"]


def test_compute_subset(mol):
  props = D.Properties(["exactmw", "NumHBD", "lipinskiHBA"])
  assert list(props.GetPropertyNames()) == ["exactmw", "NumHBD", "lipinskiHBA"]
  values = props.ComputeProperties(mol)
  assert type(values).__name__ == "_vectd"
  assert list(values) == [pytest.approx(151.063328528), 2.0, 3.0]


def test_annotate(mol):
  props = D.Properties(["exactmw", "NumHBD"])
  props.AnnotateProperties(mol)
  assert mol.GetDoubleProp("NumHBD") == 2.0
  assert "exactmw" in mol.GetPropsAsDict()
  other = Chem.MolFromSmiles("CCO")
  props.ComputeProperties(other)
  assert "exactmw" not in other.GetPropsAsDict()
  props.ComputeProperties(other, True)
  assert other.GetDoubleProp("exactmw") == pytest.approx(46.041864812)


def test_unknown_property_names():
  with pytest.raises(KeyError, match="nope"):
    D.Properties(["nope"])
  with pytest.raises(KeyError, match="nope"):
    D.Properties.GetProperty("nope")
  with pytest.raises(KeyError, match="nope"):
    D.MakePropertyRangeQuery("nope", 0, 1)


def test_builtin_functor(mol):
  functor = D.Properties.GetProperty("exactmw")
  assert (functor.GetName(), functor.GetVersion()) == ("exactmw", "1.1.0")
  assert functor(mol) == pytest.approx(151.063328528)


def test_range_query_is_open_interval(mol):
  assert D.MakePropertyRangeQuery("NumHBD", 0, 3).Match(mol)
  assert not D.MakePropertyRangeQuery("NumHBD", 0, 2).Match(mol)
  assert not D.MakePropertyRangeQuery("NumHBD", 2, 3).Match(mol)


def test_python_functor_registration(mol):
  class DoubleAtoms(Descriptors.PropertyFunctor):

    def __init__(self):
      Descriptors.PropertyFunctor.__init__(self, "double_atoms_test", "1.0")

    def __call__(self, m):
      return m.GetNumAtoms() * 2.0

  functor = DoubleAtoms()
  assert (functor.GetName(), functor.GetVersion()) == ("double_atoms_test", "1.0")
  first = D.Properties.RegisterProperty(functor)
  # registering the same name again replaces the entry in place
  assert D.Properties.RegisterProperty(DoubleAtoms()) == first
  assert "double_atoms_test" in list(D.Properties.GetAvailableProperties())
  assert list(D.Properties(["double_atoms_test"]).ComputeProperties(mol)) == [22.0]


def test_python_functor_without_data_function_cannot_be_a_query(mol):

  class Plain(Descriptors.PropertyFunctor):

    def __init__(self):
      Descriptors.PropertyFunctor.__init__(self, "plain_test", "1.0")

  D.Properties.RegisterProperty(Plain())
  query = D.MakePropertyRangeQuery("plain_test", 0, 1)
  with pytest.raises(RuntimeError):
    query.Match(mol)


def test_descriptors_module_computes_most_descriptors(mol):
  assert len(Descriptors._descList) >= 200
  values = {}
  for name, fn in Descriptors._descList:
    if name == "SPS":
      continue  # needs the new stereo perception, which is not ported yet
    values[name] = fn(mol)
  assert values["MolWt"] == pytest.approx(151.165)
  assert values["NumHDonors"] == 2
  assert values["fr_phenol"] == 1
  assert values["TPSA"] == pytest.approx(49.33)
