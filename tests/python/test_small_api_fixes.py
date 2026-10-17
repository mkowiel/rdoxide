"""Small API details: deprecated GetNumAtoms argument, error text, ring families,
_HasSubstructMatchStr and text-mode writer streams."""

import gzip
import io
import tempfile

import pytest
from rdoxide import Chem


def test_get_num_atoms_deprecated_argument_selects_explicit_count():
    mol = Chem.MolFromSmiles("n1ccccc1")
    assert mol.GetNumAtoms() == 6
    assert mol.GetNumAtoms(1) == 6
    # the int is handed on as onlyExplicit, as upstream does
    assert mol.GetNumAtoms(0) == 11


def test_range_error_names_the_upstream_source():
    mol = Chem.MolFromSmiles("C")
    with pytest.raises(RuntimeError) as e:
        mol.GetAtomWithIdx(3)
    text = str(e.value)
    assert "Code/GraphMol/ROMol.cpp" in text
    assert "Failed Expression: 3 < 1" in text


def test_ring_families_are_initialized_by_sanitization():
    assert Chem.MolFromSmiles("C1CC2C1CC2").GetRingInfo().AreRingFamiliesInitialized()
    assert not Chem.MolFromSmiles("C1CC2C1CC2", sanitize=False).GetRingInfo().AreRingFamiliesInitialized()


def test_has_substruct_match_str():
    query = Chem.MolFromSmarts("[n,p]1ccccc1")
    assert Chem._HasSubstructMatchStr(Chem.MolFromSmiles("n1ccccc1").ToBinary(), query)
    assert not Chem._HasSubstructMatchStr(Chem.MolFromSmiles("S1=CC=CC=C1").ToBinary(), query)
    with pytest.raises(ValueError):
        Chem._HasSubstructMatchStr(b"not a pickle", query)


def test_writers_need_text_mode_streams():
    with tempfile.TemporaryFile() as tmp, gzip.open(tmp) as binary:
        with pytest.raises(ValueError):
            Chem.SDWriter(binary)
    text = io.StringIO()
    writer = Chem.SDWriter(text)
    writer.write(Chem.MolFromSmiles("C"))
    writer.close()
    assert "$$$$" in text.getvalue()
