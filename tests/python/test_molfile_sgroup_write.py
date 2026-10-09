"""SubstanceGroups in written mol blocks; expectations validated against RDKit 2026.09.1."""

import pytest

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Geometry import Point3D


def depicted(smiles="CCOCC"):
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    return mol


def sru_mol():
    mol = depicted()
    group = Chem.CreateMolSubstanceGroup(mol, "SRU")
    for a in (1, 2):
        group.AddAtomWithIdx(a)
    for b in (0, 2):
        group.AddBondWithIdx(b)
    group.SetProp("LABEL", "n")
    group.SetProp("CONNECT", "HT")
    group.SetProp("SUBTYPE", "ALT")
    group.SetUnsignedProp("ID", 1)
    group.SetProp("BRKTYP", "PAREN")
    group.SetProp("ESTATE", "E")
    group.AddBracket([Point3D(1, 1, 0), Point3D(1, 2, 0)])
    group.AddBracket([Point3D(3, 1, 0), Point3D(3, 2, 0)])
    return mol


def lines(block, prefix):
    return [line for line in block.split("\n") if line.startswith(prefix)]


def test_v2000_group_lines():
    block = Chem.MolToMolBlock(sru_mol())
    assert "  5  4  0  1  0  0  0  0  0  0999 V2000" in block
    expected = [
        "M  STY  1   1 SRU", "M  SLB  1   1   1", "M  SST  1   1 ALT", "M  SCN  1   1 HT ",
        "M  SDS EXP  1   1", "M  SBT  1   1   1", "M  SAL   1  2   2   3",
        "M  SBL   1  2   1   3",
        "M  SDI   1  4    1.0000    1.0000    1.0000    2.0000",
        "M  SDI   1  4    3.0000    1.0000    3.0000    2.0000", "M  SMT   1 n",
    ]
    assert lines(block, "M  S") == expected


def test_v3000_group_lines():
    block = Chem.MolToMolBlock(sru_mol(), forceV3000=True)
    assert "M  V30 COUNTS 5 4 1 0 0" in block
    assert ("M  V30 1 SRU 1 ATOMS=(2 2 3) XBONDS=(2 1 3) SUBTYPE=ALT CONNECT=HT "
            "LABEL=n -") in block


def test_groups_survive_a_round_trip():
    for v3000 in (False, True):
        block = Chem.MolToMolBlock(sru_mol(), forceV3000=v3000)
        mol = Chem.MolFromMolBlock(block)
        (group,) = Chem.GetMolSubstanceGroups(mol)
        assert (list(group.GetAtoms()), list(group.GetBonds())) == ([1, 2], [0, 2])
        assert group.GetProp("LABEL") == "n" and len(group.GetBrackets()) == 2
        assert group.GetProp("BRKTYP") == "PAREN" and group.GetProp("ESTATE") == "E"


def test_long_data_fields_and_quoting():
    mol = depicted("CC(C)OCCNCC")
    group = Chem.CreateMolDataSubstanceGroup(mol, "fld", "x" * 150)
    group.AddAtomWithIdx(1)
    group.SetProp("FIELDINFO", "info(1)")
    block = Chem.MolToMolBlock(mol, forceV3000=True)
    assert 'FIELDNAME=fld FIELDINFO="info(1)" -' in block
    assert 'M  V30 ' + "x" * 72 + "-" in block.replace("-\nM  V30 ", "-\nM  V30 ")
    v2000 = Chem.MolToMolBlock(mol)
    assert len(lines(v2000, "M  SCD")) == 2 and len(lines(v2000, "M  SED")) == 1
    back = Chem.MolFromMolBlock(v2000)
    (group,) = Chem.GetMolSubstanceGroups(back)
    assert list(group.GetStringVectProp("DATAFIELDS")) == ["x" * 150]


def test_parent_and_attachment_lines():
    mol = depicted("CCOCC")
    first = Chem.CreateMolSubstanceGroup(mol, "MUL")
    for a in (1, 2, 3):
        first.AddAtomWithIdx(a)
    first.AddBondWithIdx(0)
    first.AddBondWithIdx(3)
    first.SetProp("MULT", "3")
    second = Chem.CreateMolSubstanceGroup(mol, "GEN")
    second.AddAtomWithIdx(0)
    second.SetUnsignedProp("PARENT", 1)
    second.SetUnsignedProp("COMPNO", 2)
    block = Chem.MolToMolBlock(mol)
    assert "M  STY  2   1 MUL   2 GEN" in block
    assert lines(block, "M  SPL") == ["M  SPL  1   2   1"]
    assert lines(block, "M  SNC") == ["M  SNC  1   2   2"]
    assert "M  SMT   1 3" in block


def test_problems_are_reported():
    mol = depicted()
    Chem.CreateMolSubstanceGroup(mol, "SRU").SetProp("BRKTYP", "X")
    with pytest.raises(RuntimeError, match="Invalid BRKTYP value 'X' for SGroup 1"):
        Chem.MolToMolBlock(mol)
    mol = depicted()
    Chem.CreateMolDataSubstanceGroup(mol, "f", "y" * 250)
    with pytest.raises(RuntimeError, match="longer than limit of 200"):
        Chem.MolToMolBlock(mol)


def test_query_atoms_and_zero_order_bonds_become_groups():
    block = Chem.MolToMolBlock(Chem.MolFromSmarts("[N;H1]-[O]-*"))
    assert "M  SED   1 [N&H1]" in block
    mol = Chem.RWMol(Chem.MolFromSmiles("[Fe].[NH3+]"))
    mol.AddBond(0, 1, Chem.BondType.ZERO)
    block = Chem.MolToMolBlock(mol.GetMol())
    for expected in ("M  ZBO  1   1   0", "M  HYD  2   1   0   2   3", "M  ZCH  1   2   1",
                     "M  SED   2 0;3", "M  SED   3 0;1"):
        assert expected in block
    back = Chem.MolFromMolBlock(block, sanitize=False)
    assert back.GetBondWithIdx(0).GetBondType() == Chem.BondType.ZERO
    assert back.GetAtomWithIdx(1).GetFormalCharge() == 1
