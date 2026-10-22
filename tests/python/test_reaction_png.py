"""Reactions in PNG metadata (`ReactionMetadataToPNG*`, `ReactionFromPNG*`)."""

import struct
import zlib

import pytest

from rdoxide import rdBase
from rdoxide.Chem import rdChemReactions


def _blank_png():
    def chunk(kind, body):
        crc = zlib.crc32(kind + body) & 0xFFFFFFFF
        return struct.pack(">I", len(body)) + kind + body + struct.pack(">I", crc)

    ihdr = struct.pack(">IIBBBBB", 1, 1, 8, 2, 0, 0, 0)
    idat = zlib.compress(b"\x00\x00\x00\x00")
    return b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) + chunk(b"IDAT", idat) + chunk(b"IEND", b"")


RXN = "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]"


def test_flag():
    assert rdBase._iostreamsEnabled


@pytest.mark.parametrize(
    "opts",
    [
        {"includePkl": True, "includeSmiles": False},
        {"includePkl": False, "includeSmiles": True},
        {"includePkl": False, "includeSmiles": False, "includeSmarts": True},
        {"includePkl": False, "includeSmiles": False, "includeRxn": True},
    ],
)
def test_round_trip(opts):
    rxn = rdChemReactions.ReactionFromSmarts(RXN)
    png = rdChemReactions.ReactionMetadataToPNGString(rxn, _blank_png(), **opts)
    assert png.startswith(b"\x89PNG")
    back = rdChemReactions.ReactionFromPNGString(png)
    assert back.GetNumReactantTemplates() == 2
    assert back.GetNumProductTemplates() == 1


def test_file(tmp_path):
    rxn = rdChemReactions.ReactionFromSmarts(RXN)
    f = tmp_path / "r.png"
    f.write_bytes(_blank_png())
    png = rdChemReactions.ReactionMetadataToPNGFile(rxn, str(f))
    f.write_bytes(png)
    back = rdChemReactions.ReactionFromPNGFile(str(f))
    assert rdChemReactions.ReactionToSmarts(back) == rdChemReactions.ReactionToSmarts(rxn)


def test_no_metadata():
    with pytest.raises(Exception):
        rdChemReactions.ReactionFromPNGString(_blank_png())
    with pytest.raises(Exception):
        rdChemReactions.ReactionFromPNGFile("/nonexistent/file.png")
