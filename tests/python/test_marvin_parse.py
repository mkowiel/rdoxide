"""Marvin (MRV) reader; cases from upstream MarvinParse/testMrvToMol.cpp."""
import os

import pytest
from rdkit import Chem

DATA = os.path.join(
    os.path.dirname(__file__), "..", "..", "upstream", "rdkit", "Code", "GraphMol",
    "MarvinParse", "test_data")

needs_data = pytest.mark.skipif(not os.path.isdir(DATA), reason="upstream checkout not available")

SIMPLE = """<cml><MDocument><MChemicalStruct><molecule molID="m1">
<atomArray><atom id="a1" elementType="C" x2="0" y2="0"/><atom id="a2" elementType="O" x2="1.5" y2="0"/></atomArray>
<bondArray><bond id="b1" atomRefs2="a1 a2" order="2"/></bondArray>
</molecule></MChemicalStruct></MDocument></cml>"""


def test_block():
    m = Chem.MolFromMrvBlock(SIMPLE)
    assert Chem.MolToSmiles(m) == "C=O"
    assert m.GetNumConformers() == 1


def test_bad_block_gives_none():
    assert Chem.MolFromMrvBlock("<cml><MDocument><MChemicalStruct/></MDocument></cml>") is not None
    assert Chem.MolFromMrvBlock("<cml><MDocument></cml>") is None


def test_missing_file_is_ioerror():
    with pytest.raises(IOError):
        Chem.MolFromMrvFile("/does/not/exist.mrv")


@needs_data
@pytest.mark.parametrize("name,natoms,nbonds", [
    ("aspirin.mrv", 13, 13),
    ("MultipleSgroup.mrv", 123, 122),
    ("EmbeddedSgroupMUL_MUL.mrv", 141, 140),
    ("MarvinOldSuperGroupTest.mrv", 89, 93),
    ("DataSgroup.mrv", 7, 6),
])
def test_files(name, natoms, nbonds):
    m = Chem.MolFromMrvFile(os.path.join(DATA, name), removeHs=False)
    assert m.GetNumAtoms() == natoms
    assert m.GetNumBonds() == nbonds


@needs_data
def test_substance_groups():
    m = Chem.MolFromMrvFile(os.path.join(DATA, "MultipleSgroup.mrv"), removeHs=False)
    assert any(sg.GetProp("TYPE") == "MUL" for sg in Chem.GetMolSubstanceGroups(m))


@needs_data
def test_bad_files_give_none():
    assert Chem.MolFromMrvFile(os.path.join(DATA, "MarvinBadX2.mrv")) is None


def test_write_round_trip():
    m = Chem.MolFromMrvBlock(SIMPLE)
    block = Chem.MolToMrvBlock(m)
    assert block.startswith('<?xml version="1.0" encoding="utf-8"?>')
    assert 'role=' not in block
    m2 = Chem.MolFromMrvBlock(block)
    assert Chem.MolToSmiles(m2) == "C=O"
    assert "\n\t<MDocument>" in Chem.MolToMrvBlock(m, prettyPrint=True)


def test_write_file(tmp_path):
    m = Chem.MolFromSmiles("CC(=O)O")
    path = str(tmp_path / "a.mrv")
    Chem.MolToMrvFile(m, path)
    assert Chem.MolToSmiles(Chem.MolFromMrvFile(path)) == "CC(=O)O"


def test_unlabelled_query_atoms_become_stars():
    # (as upstream: only labelled queries such as Q or X are rejected)
    m = Chem.MolFromSmarts("[C,N]C")
    assert 'elementType="*"' in Chem.MolToMrvBlock(m)
