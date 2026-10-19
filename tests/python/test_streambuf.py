import io

from rdoxide import Chem, rdBase


def test_streambuf_wraps_a_file_object():
    data = b"\n  RDKit          2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n"
    suppl = Chem.ForwardSDMolSupplier(rdBase.streambuf(io.BytesIO(data)))
    assert sum(1 for m in suppl if m is not None) == 1
