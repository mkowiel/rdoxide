from pathlib import Path

import pytest
from rdoxide import Chem

DATA = (
    Path(__file__).resolve().parents[2]
    / "upstream/rdkit/Code/GraphMol/FileParsers/test_data/macromols"
)
needs_data = pytest.mark.skipif(not DATA.is_dir(), reason="the upstream checkout is not there")


def test_params_defaults_and_setters():
    p = Chem.MolFromSCSRParams()
    assert p.includeLeavingGroups is True
    assert p.scsrTemplateNames == Chem.SCSRTemplateNames.AsEntered
    assert p.scsrBaseHbondOptions == Chem.SCSRBaseHbondOptions.UseSapAll
    p.scsrBaseHbondOptions = Chem.SCSRBaseHbondOptions.Auto
    p.scsrTemplateNames = Chem.SCSRTemplateNames.UseSecondName
    p.includeLeavingGroups = False
    assert p.scsrBaseHbondOptions == Chem.SCSRBaseHbondOptions.Auto
    assert p.scsrTemplateNames == Chem.SCSRTemplateNames.UseSecondName
    assert p.includeLeavingGroups is False
    with pytest.raises(ValueError):
        p.scsrBaseHbondOptions = 17


@needs_data
def test_leaving_groups_are_optional():
    path = str(DATA / "Triplet.mol")
    p = Chem.MolFromSCSRParams()
    full = Chem.MolFromSCSRFile(path, True, True, p)
    p.includeLeavingGroups = False
    bare = Chem.MolFromSCSRFile(path, True, True, p)
    assert (full.GetNumAtoms(), bare.GetNumAtoms()) == (30, 27)


@needs_data
def test_bad_input():
    assert Chem.MolFromSCSRBlock("not a mol block") is None
    with pytest.raises(OSError):
        Chem.MolFromSCSRFile(str(DATA / "does_not_exist.mol"))
    assert Chem.MolFromSCSRFile(str(DATA / "badMissingTemplate.mol")) is None


@needs_data
def test_expanded_residues_are_marked_with_groups():
    mol = Chem.MolFromSCSRFile(str(DATA / "Triplet.mol"))
    assert mol.GetNumAtoms() == 30
    assert len(Chem.GetMolSubstanceGroups(mol)) == 6
