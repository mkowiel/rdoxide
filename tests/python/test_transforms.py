"""DeleteSubstructs and FindPotentialStereoBonds (expectations from RDKit 2026.09.1)."""

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdDepictor


def delete(smi, query, **kwargs):
  res = Chem.DeleteSubstructs(Chem.MolFromSmiles(smi), Chem.MolFromSmarts(query), **kwargs)
  return Chem.MolToSmiles(res)


@pytest.mark.parametrize("smi, query, kwargs, expected", [
  ("CCOC", "OC", {}, "C"),
  ("CCOC", "OC", dict(onlyFrags=True), "CCOC"),
  ("CCOCCl.Cl", "Cl", dict(onlyFrags=True), "CCOCCl"),
  ("CCOCCl.Cl", "Cl", {}, "CCOC"),
  ("CC", "N", {}, "CC"),
])
def test_delete_substructs(smi, query, kwargs, expected):
  assert delete(smi, query, **kwargs) == expected


def test_original_is_not_modified():
  mol = Chem.MolFromSmiles("CCOC")
  Chem.DeleteSubstructs(mol, Chem.MolFromSmarts("OC"))
  assert mol.GetNumAtoms() == 4


def test_conformers_follow_the_kept_atoms():
  mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
  rdDepictor.Compute2DCoords(mol)
  res = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts("[OX2H]"))
  assert res.GetNumAtoms() == mol.GetNumAtoms() - 1
  kept = mol.GetConformer().GetPositions()
  assert res.GetConformer().GetPositions().tolist() == kept[:2].tolist() + kept[3:].tolist()


def stereo_bonds(smi, **kwargs):
  mol = Chem.MolFromSmiles(smi)
  Chem.FindPotentialStereoBonds(mol, **kwargs)
  return [(str(b.GetStereo()), list(b.GetStereoAtoms()))
          for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE]


def test_potential_stereo_bonds_get_stereo_atoms():
  assert stereo_bonds("CC=CC") == [("STEREOANY", [0, 3])]
  assert stereo_bonds("CC=CC=CC")[0][0] == "STEREOANY"


@pytest.mark.parametrize("smi", ["C=C", "CC=C(C)C", "C1CC=CCC1"])
def test_bonds_without_distinguishable_ends_are_left_alone(smi):
  assert all(stereo == "STEREONONE" for stereo, _ in stereo_bonds(smi))


def test_mol_is_marked_and_second_call_is_a_no_op():
  mol = Chem.MolFromSmiles("CC=CC")
  Chem.FindPotentialStereoBonds(mol)
  assert mol.HasProp("_BondsPotentialStereo")
  Chem.FindPotentialStereoBonds(mol, cleanIt=True)
  assert mol.GetBondWithIdx(1).GetStereo() == Chem.BondStereo.STEREOANY
