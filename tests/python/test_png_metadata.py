import struct
import zlib

import pytest
from rdoxide import Chem


def _png():
    def chunk(kind, data):
        body = kind + data
        return struct.pack(">I", len(data)) + body + struct.pack(">I", zlib.crc32(body))

    ihdr = struct.pack(">IIBBBBB", 1, 1, 8, 2, 0, 0, 0)
    idat = zlib.compress(b"\x00\xff\x00\x00")
    return (
        b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) + chunk(b"IDAT", idat) + chunk(b"IEND", b"")
    )


def test_round_trip_all_formats():
    m = Chem.MolFromSmiles("C[C@H](F)Cl")
    ps = Chem.PNGMetadataParams()
    ps.includeMol = True
    png = Chem.MolMetadataToPNGString(m, _png(), ps)
    meta = Chem.MetadataFromPNGString(png)
    assert {k.split()[0] for k in meta} == {"rdkitPKL", "SMILES", "MOL"}
    for tag in ("rdkitPKL", "SMILES", "MOL"):
        (m2,) = Chem.MolsFromPNGString(png, tag=tag)
        assert Chem.MolToSmiles(m2) == Chem.MolToSmiles(m)
    assert Chem.MolToSmiles(Chem.MolFromPNGString(png)) == Chem.MolToSmiles(m)
    assert len(Chem.MetadataFromPNGString(png, asList=True)) == 3


def test_flags_and_plain_metadata():
    m = Chem.MolFromSmiles("CCO")
    png = Chem.MolMetadataToPNGString(m, _png(), includePkl=False, includeSmiles=True)
    assert [k.split()[0] for k in Chem.MetadataFromPNGString(png)] == ["SMILES"]
    png = Chem.AddMetadataToPNGString({"note": "hello"}, png)
    assert Chem.MetadataFromPNGString(png)["note"] == b"hello"
    # the image data are untouched
    assert png.startswith(_png()[:8])


def test_files(tmp_path):
    f = tmp_path / "a.png"
    f.write_bytes(_png())
    m = Chem.MolFromSmiles("c1ccccc1O")
    f.write_bytes(Chem.MolMetadataToPNGFile(m, str(f)))
    assert Chem.MolToSmiles(Chem.MolFromPNGFile(str(f))) == "Oc1ccccc1"
    assert len(Chem.MolsFromPNGFile(str(f))) == 1
    f.write_bytes(Chem.AddMetadataToPNGFile({"k": "v"}, str(f)))
    assert Chem.MetadataFromPNGFile(str(f))["k"] == b"v"
    with pytest.raises(IOError):
        Chem.MolFromPNGFile(str(tmp_path / "missing.png"))


def test_bad_png():
    assert Chem.MolFromPNGString(b"not a png") is None
    assert Chem.MolFromPNGString(_png()) is None  # no molecule metadata
    with pytest.raises(Exception):
        Chem.MetadataFromPNGString(b"not a png")


def test_params_semantics():
    ps = Chem.PNGMetadataParams()
    assert ps.includePkl and ps.includeSmiles and not ps.includeMol
    ps.smilesWriteParams.canonical = False
    assert ps.smilesWriteParams.canonical is False
    with pytest.raises(AttributeError):
        ps.notAnOption = 1


def test_zlib_streams():
    # a dynamic-Huffman zTXt chunk, as written by zlib itself
    text = b"CC(=O)O " * 200
    body = b"SMILES\x00\x00" + zlib.compress(text, 9)
    png = _png()
    end = png[-12:]
    blk = b"zTXt" + body
    new = png[:-12] + struct.pack(">I", len(body)) + blk + struct.pack(">I", zlib.crc32(blk)) + end
    assert Chem.MetadataFromPNGString(new)["SMILES"] == text
