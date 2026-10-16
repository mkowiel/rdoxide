"""DataStructs.InitFromDaylightString (expectations from the reference RDKit)."""

import pytest
from rdoxide import DataStructs


@pytest.mark.parametrize("cls", [DataStructs.ExplicitBitVect, DataStructs.SparseBitVect])
def test_decodes_bits(cls):
    bv = cls(32)
    DataStructs.InitFromDaylightString(bv, "+.ab3")
    assert list(bv.GetOnBits()) == [5, 12, 15, 16, 18, 21, 22, 23]


def test_pad_marker_trims_trailing_bytes():
    bv = DataStructs.ExplicitBitVect(32)
    DataStructs.InitFromDaylightString(bv, "+.ab2")
    assert list(bv.GetOnBits()) == [5, 12, 15]


def test_trailing_newline_and_reset():
    bv = DataStructs.ExplicitBitVect(32)
    bv.SetBit(31)
    DataStructs.InitFromDaylightString(bv, "+.ab3\n")
    assert 31 not in bv.GetOnBits()
    assert bv.GetNumOnBits() == 8


@pytest.mark.parametrize("bad", ["", "abcd", "abcd4"])
def test_bad_strings(bad):
    with pytest.raises(Exception):
        DataStructs.InitFromDaylightString(DataStructs.ExplicitBitVect(64), bad)
