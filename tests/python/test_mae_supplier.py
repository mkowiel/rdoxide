import gzip
import io

import pytest
from rdoxide import Chem

MAE = """{
 s_m_m2io_version
 :::
 2.0.0
}

f_m_ct {
 s_m_title
 i_m_ct_format
 b_sd_flag
 r_sd_weight
 s_st_Chirality_1
 :::
 "my \\"mol\\""
 2
 1
 12.5
 2_S_1_3_4_5
 m_atom[5] {
  # First column is atom index #
  i_m_mmod_type
  r_m_x_coord
  r_m_y_coord
  r_m_z_coord
  i_m_atomic_number
  i_m_formal_charge
  s_m_label
  :::
  1 3 0.0 0.0 0.0 9 0 F
  2 3 1.0 0.0 0.0 6 0 <>
  3 3 1.5 1.0 0.0 17 0 "C l"
  4 3 1.5 -0.5 0.9 35 0 <>
  5 3 1.8 -0.5 -0.9 1 0 <>
  :::
 }
 m_bond[4] {
  # First column is bond index #
  i_m_from
  i_m_to
  i_m_order
  :::
  1 1 2 1
  2 2 3 1
  3 2 4 1
  4 2 5 1
  :::
 }
}
"""


def test_properties_and_structure():
    suppl = Chem.MaeMolSupplier()
    suppl.SetData(MAE, removeHs=False)
    assert len(suppl) == 1
    m = next(suppl)
    assert m.GetProp("_Name") == 'my "mol"'
    assert m.GetNumAtoms() == 5
    assert m.GetBoolProp("b_sd_flag") is True
    assert m.GetDoubleProp("r_sd_weight") == 12.5
    assert m.GetIntProp("i_m_ct_format") == 2
    assert m.GetAtomWithIdx(0).GetProp("s_m_label") == "F"
    assert m.GetAtomWithIdx(2).GetProp("s_m_label") == "C l"
    assert not m.GetAtomWithIdx(1).HasProp("s_m_label")
    assert m.GetConformer().Is3D()
    assert m.GetAtomWithIdx(1).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    with pytest.raises(StopIteration):
        next(suppl)
    assert suppl.atEnd()
    suppl.reset()
    assert not suppl.atEnd()
    assert Chem.MolToSmiles(suppl[0]) == Chem.MolToSmiles(m)


def test_stream_and_file(tmp_path):
    path = tmp_path / "a.mae"
    path.write_text(MAE)
    assert len(Chem.MaeMolSupplier(str(path))) == 1
    assert len(Chem.MaeMolSupplier(io.BytesIO(MAE.encode()))) == 1
    gz = tmp_path / "a.maegz"
    gz.write_bytes(gzip.compress(MAE.encode()))
    with Chem.MaeMolSupplier(gzip.open(gz)) as suppl:
        assert next(suppl).GetNumAtoms() == 4  # the H is removed
    with pytest.raises(IOError):
        Chem.MaeMolSupplier(str(tmp_path))
    with pytest.raises(ValueError):
        Chem.MaeMolSupplier(io.StringIO(MAE))


def test_indexing():
    suppl = Chem.MaeMolSupplier()
    suppl.SetData(MAE + MAE)
    assert len(suppl) == 2
    assert suppl[-1].GetNumAtoms() == 4
    with pytest.raises(IndexError):
        suppl[2]
    with pytest.raises(IndexError):
        suppl[-3]


@pytest.mark.parametrize(
    "text, message",
    [
        ("f_m_ct {\n s_m_title\n :::\n \"a\"\n", "Line 5, column 1: Missing '}' for block."),
        ("f_m_ct {\n s_m_title\n :::\n \"a", "Line 4, column 4: Unterminated quoted string at EOF."),
        ("fm_ct { }", "Line 1, column 2: Bad format for outer block name"),
        ("f_m_ct { s_m_title i_m_foo ::: x }", "Line 1, column 34: Unexpected character."),
        ("f_m_ct { s_m_ ::: x }", "Line 1, column 14: Bad format for property"),
        ("f_m_ct s_m_title ::: x }", "Line 1, column 8: Missing '{' for outer block."),
    ],
)
def test_read_errors(text, message):
    suppl = Chem.MaeMolSupplier()
    with pytest.raises(RuntimeError) as exc:
        suppl.SetData(text)
    assert message in str(exc.value)


def test_build_errors():
    suppl = Chem.MaeMolSupplier()
    suppl.SetData("f_m_ct { s_m_title ::: x }")
    with pytest.raises(RuntimeError, match="Indexed block not found: m_atom"):
        next(suppl)
    atoms = (
        "f_m_ct { s_m_title ::: x m_atom[1] { i_m_atomic_number r_m_x_coord "
        "r_m_y_coord r_m_z_coord ::: 1 %s ::: } }"
    )
    suppl.SetData(atoms % "6 0.0abc 0 0")
    with pytest.raises(RuntimeError, match="Bad floating point representation."):
        next(suppl)
    suppl.SetData(atoms % "99999 0.0 0 0")
    with pytest.raises(RuntimeError, match="Atomic number not found"):
        next(suppl)


def test_no_data_is_an_error():
    suppl = Chem.MaeMolSupplier()
    with pytest.raises(RuntimeError):
        next(suppl)
