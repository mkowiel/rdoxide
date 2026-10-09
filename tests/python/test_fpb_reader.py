"""FPBReader / MultiFPBReader on files copied from RDKit (expectations from testFPB.py)."""

from pathlib import Path

import pytest
from rdkit import DataStructs

DATA = Path(__file__).parent.parent / "data" / "fpb"


@pytest.fixture(params=[False, True], ids=["memory", "lazy"])
def reader(request):
    r = DataStructs.FPBReader(str(DATA / "zim.head100.fpb"), request.param)
    r.Init()
    return r


def test_basics(reader):
    assert (len(reader), reader.GetNumBits()) == (100, 2048)
    assert reader.GetId(3) == "ZINC04803506"
    fp, name = reader[0]
    assert name == "ZINC00902219" and fp.GetNumOnBits() == 17
    assert len(reader.GetBytes(0)) == 256


def test_searches(reader):
    bv = reader.GetBytes(0)
    assert reader.GetTanimoto(1, bv) == pytest.approx(0.3704, abs=1e-4)
    nbrs = reader.GetTanimotoNeighbors(bv, threshold=0.3)
    assert [n[1] for n in nbrs[:2]] == [0, 1] and len(nbrs) == 5
    assert reader.GetTverskyNeighbors(bv, 1, 1, threshold=0.3) == nbrs
    assert reader.GetContainingNeighbors(reader.GetBytes(1)) == (1, 2, 3, 4)


def test_errors():
    with pytest.raises(RuntimeError, match="Bad input file"):
        DataStructs.FPBReader(str(DATA / "missing.fpb"))
    r = DataStructs.FPBReader(str(DATA / "zim.head100.fpb"))
    with pytest.raises(RuntimeError):
        len(r)
    r.Init()
    with pytest.raises(ValueError):
        r.GetId(100)


def test_multi_reader_threads():
    multi = DataStructs.MultiFPBReader()
    for i in (1, 2):
        assert multi.AddReader(DataStructs.FPBReader(str(DATA / f"zinc_random200.{i}.patt.fpb"))) == i
    multi.Init()
    assert (len(multi), multi.GetNumBits()) == (2, 1024)
    bv = multi.GetReader(1).GetBytes(3)
    serial = multi.GetTanimotoNeighbors(bv, 0.5)
    assert (1.0, 3, 1) in serial
    for threads in (2, 4, 0):
        assert multi.GetTanimotoNeighbors(bv, 0.5, numThreads=threads) == serial
    assert (3, 1) in multi.GetContainingNeighbors(bv, numThreads=2)
