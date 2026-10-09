"""ForwardSDMolSupplier. Expected values were checked against RDKit 2026.09.1."""
import gzip
import io
import os

import pytest

from rdkit import Chem

SDF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "molfiles", "first_200.props.sdf")


@pytest.fixture(scope="module")
def raw():
    with open(SDF, "rb") as f:
        return f.read()


def atom_counts(suppl):
    return [None if m is None else m.GetNumAtoms() for m in suppl]


def test_sources(raw, tmp_path):
    assert len(atom_counts(Chem.ForwardSDMolSupplier(io.BytesIO(raw)))) == 200
    assert len(atom_counts(Chem.ForwardSDMolSupplier(SDF))) == 200
    with open(SDF, "rb") as f:
        assert len(atom_counts(Chem.ForwardSDMolSupplier(f))) == 200
    gz = tmp_path / "a.sdf.gz"
    gz.write_bytes(gzip.compress(raw))
    assert len(atom_counts(Chem.ForwardSDMolSupplier(gzip.open(gz)))) == 200


def test_input_errors(tmp_path):
    with pytest.raises(ValueError, match="binary mode"):
        Chem.ForwardSDMolSupplier(io.StringIO("x"))
    with pytest.raises(OSError, match="Bad input file"):
        Chem.ForwardSDMolSupplier(str(tmp_path / "missing.sdf"))
    (tmp_path / "empty.sdf").write_bytes(b"")
    assert atom_counts(Chem.ForwardSDMolSupplier(str(tmp_path / "empty.sdf"))) == []
    assert atom_counts(Chem.ForwardSDMolSupplier(io.BytesIO(b""))) == []


def test_iteration_state(raw):
    s = Chem.ForwardSDMolSupplier(io.BytesIO(raw))
    assert not s.atEnd() and not s.GetEOFHitOnRead()
    m = next(s)
    assert m.GetNumAtoms() == 9 and m.GetProp("AMW") == "122.12344"
    assert s.GetProcessPropertyLists()
    s.SetProcessPropertyLists(False)
    assert not s.GetProcessPropertyLists()
    for _ in range(199):
        next(s)
    with pytest.raises(StopIteration, match="End of supplier hit"):
        next(s)
    assert s.atEnd() and s.GetEOFHitOnRead()
    with pytest.raises(StopIteration):
        next(s)
    assert not hasattr(s, "__len__") and not hasattr(s, "__getitem__")


def test_damaged_records_yield_none(raw):
    first, second = raw.split(b"$$$$\n")[:2]
    bad = first + b"$$$$\n" + b"garbage\nmore\n\n\n$$$$\n" + second + b"$$$$\n"
    assert atom_counts(Chem.ForwardSDMolSupplier(io.BytesIO(bad))) == [9, None, 20]
    assert atom_counts(Chem.ForwardSDMolSupplier(io.BytesIO(first + b"$$$$\n" + second + b"$$$$"))) == [9, 20]
    assert atom_counts(Chem.ForwardSDMolSupplier(io.BytesIO(first + b"$$$$\n\n\n\n\n\n\n"))) == [9, None]


def test_context_manager_closes(raw):
    with Chem.ForwardSDMolSupplier(io.BytesIO(raw)) as s:
        assert next(s).GetNumAtoms() == 9
    for call in (s.atEnd, lambda: next(s)):
        with pytest.raises(RuntimeError, match="no stream"):
            call()
