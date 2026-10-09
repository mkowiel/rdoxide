"""PDB reading/writing (rdmolfiles). Written to pass on upstream RDKit as well."""

import io
import os

import pytest
from rdkit import Chem

GLY = """\
HEADER    TEST PROTEIN                            01-JAN-00   1ABC
COMPND    MOL_ID: 1;
COMPND   2 MOLECULE: GLYCINE TEST;
ATOM      1  N   GLY A   1      -0.966   0.493   1.500  1.00  0.00           N
ATOM      2  CA  GLY A   1       0.257   0.418   0.692  1.00  0.00           C
ATOM      3  C   GLY A   1      -0.094   0.017  -0.716  1.00  0.00           C
ATOM      4  O   GLY A   1      -1.056  -0.682  -0.923  1.00  0.00           O
END
"""


def info(a):
    i = a.GetPDBResidueInfo()
    return (i.GetName(), i.GetResidueName(), i.GetResidueNumber(), i.GetChainId(), i.GetSerialNumber())


def test_basic_block():
    m = Chem.MolFromPDBBlock(GLY)
    assert m.GetNumAtoms() == 4
    assert m.GetProp("_Name") == "GLYCINE TEST"
    assert Chem.MolToSmiles(m) == "NCC=O"
    assert info(m.GetAtomWithIdx(1)) == (" CA ", "GLY", 1, "A", 2)
    assert m.GetNumConformers() == 1
    assert m.GetConformer().Is3D()
    # C=O is a standard residue double bond
    assert m.GetBondBetweenAtoms(2, 3).GetBondType() == Chem.BondType.DOUBLE


def test_bytes_and_file(tmp_path):
    assert Chem.MolFromPDBBlock(GLY.encode()).GetNumAtoms() == 4
    p = tmp_path / "x.pdb"
    p.write_text(GLY)
    assert Chem.MolFromPDBFile(str(p)).GetNumAtoms() == 4
    with pytest.raises(OSError):
        Chem.MolFromPDBFile(str(tmp_path / "missing.pdb"))


def test_empty_and_garbage():
    assert Chem.MolFromPDBBlock("") is None
    assert Chem.MolFromPDBBlock("hello\nworld\n") is None


def test_proximity_bonding_flag():
    m = Chem.MolFromPDBBlock(GLY, proximityBonding=False)
    assert m.GetNumBonds() == 0


def test_alt_locs_and_flavor():
    block = GLY.replace(
        "ATOM      3  C   GLY A   1",
        "ATOM      3  C  BGLY A   1",
    )
    assert Chem.MolFromPDBBlock(block, sanitize=False).GetNumAtoms() == 3
    assert Chem.MolFromPDBBlock(block, sanitize=False, flavor=1).GetNumAtoms() == 4


def test_charges_and_isotopes():
    block = (
        "HETATM    1  D   UNL     1       0.000   0.000   0.000  1.00  0.00           D  \n"
        "HETATM    2  N   UNL     1       1.500   0.000   0.000  1.00  0.00           N1+\n"
        "HETATM    3  O   UNL     1       3.000   0.000   0.000  1.00  0.00           O1-\n"
    )
    m = Chem.MolFromPDBBlock(block, sanitize=False, proximityBonding=False)
    assert m.GetAtomWithIdx(0).GetIsotope() == 2
    assert m.GetAtomWithIdx(1).GetFormalCharge() == 1
    assert m.GetAtomWithIdx(2).GetFormalCharge() == -1
    assert m.GetAtomWithIdx(1).GetPDBResidueInfo().GetIsHeteroAtom()


def test_conect_orders():
    block = (
        "HETATM    1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  UNL     1       1.200   0.000   0.000  1.00  0.00           O\n"
        "CONECT    1    2    2\n"
        "CONECT    2    1    1\n"
    )
    m = Chem.MolFromPDBBlock(block, proximityBonding=False)
    assert m.GetBondBetweenAtoms(0, 1).GetBondType() == Chem.BondType.DOUBLE


def test_multi_model():
    block = GLY.replace("END\n", "") + "ENDMDL\n" + "".join(
        l.replace("0.", "1.", 1) for l in GLY.splitlines(True) if l.startswith("ATOM")
    ) + "ENDMDL\nEND\n"
    m = Chem.MolFromPDBBlock(block)
    assert m.GetNumConformers() == 2
    assert m.GetNumAtoms() == 4
    # flavor 2 stops at the first ENDMDL
    assert Chem.MolFromPDBBlock(block, flavor=2).GetNumConformers() == 1


def test_crlf():
    m = Chem.MolFromPDBBlock(GLY.replace("\n", "\r\n"))
    assert m.GetNumAtoms() == 4


def test_write_block_roundtrip():
    m = Chem.MolFromPDBBlock(GLY)
    block = Chem.MolToPDBBlock(m)
    lines = block.splitlines()
    assert lines[0] == "COMPND    GLYCINE TEST"
    assert lines[1].startswith("ATOM      1  N   GLY A   1      -0.966   0.493   1.500  1.00  0.00           N")
    assert lines[-1] == "END"
    m2 = Chem.MolFromPDBBlock(block)
    assert Chem.MolToSmiles(m2) == Chem.MolToSmiles(m)


def test_write_flavors():
    m = Chem.MolFromPDBBlock(GLY)
    block = Chem.MolToPDBBlock(m, flavor=1 | 16 | 32)
    assert "TER       5      GLY A   1" in block
    assert block.splitlines()[-2].startswith("MASTER")
    assert "CONECT" not in Chem.MolToPDBBlock(m, flavor=2 | 8)
    both = Chem.MolToPDBBlock(m, flavor=4)
    assert both.count("CONECT") > Chem.MolToPDBBlock(m).count("CONECT")


def test_write_without_residue_info():
    m = Chem.MolFromSmiles("CC(=O)[O-]")
    block = Chem.MolToPDBBlock(m)
    assert "HETATM    1  C1  UNL     1" in block
    assert "O1-" in block


def test_writer_class(tmp_path):
    m = Chem.MolFromPDBBlock(GLY)
    buf = io.StringIO()
    w = Chem.PDBWriter(buf, flavor=1)
    w.write(m)
    w.write(m)
    assert w.NumMols() == 2
    w.close()
    text = buf.getvalue()
    assert text.count("MODEL") == 2 and text.count("ENDMDL") == 2
    with pytest.raises(RuntimeError):
        w.write(m)
    path = str(tmp_path / "o.pdb")
    with Chem.PDBWriter(path) as w2:
        w2.write(m)
    assert os.path.getsize(path) > 0
    Chem.MolToPDBFile(m, str(tmp_path / "p.pdb"))
    assert Chem.MolFromPDBFile(str(tmp_path / "p.pdb")).GetNumAtoms() == 4


def test_upstream_quirks():
    atoms = [l for l in GLY.splitlines() if l.startswith("ATOM")]
    # a '-' followed by a digit in the charge columns gives a *positive* charge
    a = atoms[0][:78] + "-3"
    m = Chem.MolFromPDBBlock(a + "\n", sanitize=False)
    assert m.GetAtomWithIdx(0).GetFormalCharge() == 3
    # the line-end search starts at the second character, so a blank line
    # swallows the record that follows it
    block = "\n".join([atoms[0], "", atoms[1], atoms[2]]) + "\n"
    assert Chem.MolFromPDBBlock(block, sanitize=False).GetNumAtoms() == 2
