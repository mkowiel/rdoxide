import io

import pytest
from rdoxide import Chem

# expectations from RDKit 2026.09.1 (`MaeWriter.GetText`)
ETHANE_NAMED = """f_m_ct {
  i_m_ct_stereo_status
  s_m_title
  s_rdkit_note
  :::
  1
  "a \\"b\\""
  "hello world"
  m_atom[2] {
    # First column is Index #
    r_m_x_coord
    r_m_y_coord
    r_m_z_coord
    i_m_atomic_number
    i_m_formal_charge
    s_rdkit_label
    :::
    1 0.000000 0.000000 0.000000 6 0 <>
    2 1.299038 0.750000 0.000000 6 0 L1
    :::
  }
  m_bond[1] {
    # First column is Index #
    i_m_from
    i_m_order
    i_m_to
    :::
    1 1 1 2
    :::
  }
}

"""


def ethane():
    m = Chem.MolFromSmiles("CC")
    m.SetProp("_Name", 'a "b"')
    m.SetProp("note", "hello world")
    m.GetAtomWithIdx(1).SetProp("label", "L1")
    return m


def test_get_text_layout():
    assert Chem.MaeWriter.GetText(ethane()) == ETHANE_NAMED


def test_file_header_and_blocks():
    out = io.StringIO()
    with Chem.MaeWriter(out) as w:
        w.write(Chem.MolFromSmiles("CC"))
        w.write(Chem.MolFromSmiles("O"))
        assert w.NumMols() == 2
    text = out.getvalue()
    assert text.startswith(" {\n  s_m_m2io_version\n  :::\n  2.0.0\n}\n\nf_m_ct {")
    assert text.count("f_m_ct {") == 2
    assert text.endswith("}\n\n")
    assert Chem.MaeWriter.GetText(Chem.MolFromSmiles("CC")) in text


def test_no_mols_writes_nothing():
    out = io.StringIO()
    Chem.MaeWriter(out).close()
    assert out.getvalue() == ""


def test_selected_props():
    m = ethane()
    m.SetProp("drop", "x")
    text = Chem.MaeWriter.GetText(m, props_list=["note"])
    assert "s_rdkit_note" in text
    assert "drop" not in text and "label" not in text


def test_stereo_dative_and_properties_types():
    m = Chem.MolFromSmiles("C[C@H](F)Cl")
    assert "s_st_Chirality_1" in Chem.MaeWriter.GetText(m)
    assert "2_ANS_1_3_4" in Chem.MaeWriter.GetText(m) or "2_ANR_1_3_4" in Chem.MaeWriter.GetText(m)
    m = Chem.MolFromSmiles("CC")
    m.SetBoolProp("flag", True)
    m.SetDoubleProp("x", 0.1234567)
    m.SetIntProp("n", 3)
    m.GetAtomWithIdx(0).SetBoolProp("abool", False)
    text = Chem.MaeWriter.GetText(m)
    assert "b_rdkit_flag" in text and "r_rdkit_x" in text and "0.123457" in text
    assert "b_rdkit_abool" in text


def test_errors():
    with pytest.raises(ValueError, match="without atoms"):
        Chem.MaeWriter.GetText(Chem.Mol())
    m = Chem.RWMol(Chem.MolFromSmiles("CC"))
    m.GetBondWithIdx(0).SetBondType(Chem.BondType.QUADRUPLE)
    with pytest.raises(ValueError, match="Bond 0 has a type that is not supported"):
        Chem.MaeWriter.GetText(m)
    with pytest.raises(ValueError, match="Bad Conformer Id"):
        m2 = Chem.AddHs(Chem.MolFromSmiles("C"))
        Chem.rdDepictor.Compute2DCoords(m2)
        Chem.MaeWriter.GetText(m2, confId=5)
    m3 = Chem.MolFromSmiles("CC")
    m3.GetAtomWithIdx(0).SetIntProp("s_user_bad", 1)
    with pytest.raises(RuntimeError, match="not an integer value"):
        Chem.MaeWriter.GetText(m3)


def test_closed_writer_and_text_stream_required(tmp_path):
    out = io.StringIO()
    w = Chem.MaeWriter(out)
    w.close()
    with pytest.raises(RuntimeError, match="no output stream"):
        w.write(Chem.MolFromSmiles("C"))
    with pytest.raises(ValueError):
        Chem.MaeWriter(io.BytesIO())
    with pytest.raises(IOError):
        Chem.MaeWriter(str(tmp_path / "missing" / "x.mae"))


def test_roundtrip_with_the_supplier(tmp_path):
    smiles = "C[C@H](F)c1ccccc1O"
    m = Chem.MolFromSmiles(smiles)
    m.SetProp("_Name", "test mol")
    path = tmp_path / "out.mae"
    with Chem.MaeWriter(str(path)) as w:
        w.write(m)
    back = next(Chem.MaeMolSupplier(str(path)))
    assert Chem.MolToSmiles(back) == Chem.MolToSmiles(m)
    assert back.GetProp("_Name") == "test mol"
