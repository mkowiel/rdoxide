"""SubstructLibrary bindings; cases follow upstream rdSubstructLibrary rough_test.py."""
import io
import pickle

import pytest

from rdoxide import Chem
from rdoxide.Chem import rdSubstructLibrary as ssl

HOLDERS = [ssl.MolHolder, ssl.CachedMolHolder, ssl.CachedSmilesMolHolder,
           ssl.CachedTrustedSmilesMolHolder]
SMILES = ["c1ccccc1", "CCO", "c1ccccc1O", "CCCC", "C1CCCCC1", "c1ccncc1", "OCCO", "CC(=O)O"]


def make_library(holder_cls, fp_cls=None, key_cls=None):
    holder = holder_cls()
    fps = fp_cls() if fp_cls else None
    keys = key_cls() if key_cls else None
    if fps is not None and keys is not None:
        lib = ssl.SubstructLibrary(holder, fps, keys)
    elif fps is not None:
        lib = ssl.SubstructLibrary(holder, fps)
    elif keys is not None:
        lib = ssl.SubstructLibrary(holder, keys)
    else:
        lib = ssl.SubstructLibrary(holder)
    for idx, smi in enumerate(SMILES):
        m = Chem.MolFromSmiles(smi)
        m.SetProp("_Name", str(idx))
        assert lib.AddMol(m) == idx
    return lib


@pytest.mark.parametrize("holder_cls", HOLDERS)
@pytest.mark.parametrize("fp_cls", [None, ssl.PatternHolder, ssl.TautomerPatternHolder])
@pytest.mark.parametrize("key_cls", [None, ssl.KeyFromPropHolder])
def test_search_and_pickle(holder_cls, fp_cls, key_cls):
    lib = make_library(holder_cls, fp_cls, key_cls)
    assert len(lib) == 8
    for copy in (lib, pickle.loads(pickle.dumps(lib)),
                 ssl.SubstructLibrary(lib.Serialize())):
        assert len(copy) == 8
        q = Chem.MolFromSmarts("c1ccccc1")
        assert list(copy.GetMatches(q)) == [0, 2]
        assert copy.CountMatches(q) == 2
        assert copy.HasMatch(q)
        assert not copy.HasMatch(Chem.MolFromSmiles("CCCCCCCCCC"))
        assert list(copy.GetMatches(q, 1, 8)) == [2]
        assert copy.CountMatches(q, 1, 8) == 1
        assert copy.HasMatch(q, 1, 8)
        assert list(copy.GetMatches(q, numThreads=2, maxResults=1)) == [0]
        if key_cls:
            assert copy.GetKeyHolder().GetKey(3) == "3"
            assert list(copy.GetKeyHolder().GetKeys([2, 0])) == ["2", "0"]
        else:
            assert copy.GetKeyHolder() is None
        assert (copy.GetFpHolder() is None) == (fp_cls is None)


def test_holder_types_round_trip():
    lib = make_library(ssl.CachedSmilesMolHolder, ssl.TautomerPatternHolder, ssl.KeyFromPropHolder)
    assert isinstance(lib.GetMolHolder(), ssl.CachedSmilesMolHolder)
    assert isinstance(lib.GetFpHolder(), ssl.TautomerPatternHolder)
    assert isinstance(lib.GetFpHolder(), ssl.PatternHolder)
    assert isinstance(lib.GetKeyHolder(), ssl.KeyFromPropHolder)
    assert len(lib.GetMolHolder()) == len(lib.GetFpHolder()) == 8
    # holders are shared with the library
    holder = lib.GetMolHolder()
    holder.AddSmiles("CCN")
    assert len(lib) == 9
    assert Chem.MolToSmiles(lib.GetMol(8)) == "CCN"


def test_parameters_overload():
    lib = ssl.SubstructLibrary()
    for smi in ("C[C@H](F)Cl", "C[C@@H](F)Cl", "CC(F)Cl"):
        lib.AddMol(Chem.MolFromSmiles(smi))
    q = Chem.MolFromSmiles("C[C@H](F)Cl")
    params = Chem.SubstructMatchParameters()
    assert sorted(lib.GetMatches(q, params)) == [0, 1, 2]
    params.useChirality = True
    assert sorted(lib.GetMatches(q, params)) == [0]
    assert lib.CountMatches(q, params) == 1
    assert lib.HasMatch(q, params)
    assert list(lib.GetMatches(q, 1, 3, params)) == []
    assert sorted(lib.GetMatches(q, useChirality=False)) == [0, 1, 2]
    assert list(lib.GetMatches(q)) == [0]
    with pytest.raises(TypeError):
        lib.GetMatches(q, bogus=1)


def test_search_order():
    for keyholder in (None, ssl.KeyFromPropHolder()):
        lib = ssl.SubstructLibrary(ssl.MolHolder(), keyholder)
        for idx, smi in enumerate(("CCCOC", "CCCCOCC", "CCOC", "COC", "CCCCCOC")):
            m = Chem.MolFromSmiles(smi)
            m.SetProp("_Name", str(idx))
            lib.AddMol(m)
        lib.SetSearchOrder((3, 2, 0, 1, 4))
        assert lib.GetSearchOrder() == (3, 2, 0, 1, 4)
        qm = Chem.MolFromSmiles("COC")
        assert list(lib.GetMatches(qm, maxResults=2)) == [3, 2]
        if keyholder:
            assert keyholder.GetPropName() == "_Name"
            assert list(lib.GetKeyHolder().GetKeys(lib.GetMatches(qm, maxResults=2))) == ["3", "2"]
        lib.SetSearchOrder(None)
        assert lib.GetSearchOrder() == ()
        lib.SetSearchOrder([3, 2, 0, 1, 4])
        lib.SetSearchOrder([])
        assert lib.GetSearchOrder() == ()
        with pytest.raises(IndexError):
            lib.SetSearchOrder([7])


def test_key_prop_holder():
    for propname in (None, "foo"):
        keys = ssl.KeyFromPropHolder() if propname is None else ssl.KeyFromPropHolder(propname)
        lib = ssl.SubstructLibrary(ssl.MolHolder(), keys)
        m = Chem.MolFromSmiles("CCC")
        m.SetProp(propname or "_Name", "Z11234")
        lib.AddMol(m)
        assert list(lib.GetKeyHolder().GetKeys(lib.GetMatches(m))) == ["Z11234"]
        lib.AddMol(Chem.MolFromSmiles("CCC"))
        assert lib.GetKeyHolder().GetKey(1) == "LIBIDX-1"
        with pytest.raises(IndexError):
            lib.GetKeyHolder().GetKey(2)


def test_add_patterns():
    for fp_cls in (ssl.PatternHolder, ssl.TautomerPatternHolder):
        lib = make_library(ssl.CachedMolHolder)
        assert lib.GetFpHolder() is None
        holder = fp_cls(1024)
        ssl.AddPatterns(lib, holder, 2)
        assert len(lib.GetFpHolder()) == 8
        assert len(holder) == 8
        assert list(lib.GetMatches(Chem.MolFromSmiles("c1ccccc1"))) == [0, 2]
        with pytest.raises(RuntimeError):
            ssl.AddPatterns(lib, 1)
    lib = make_library(ssl.CachedSmilesMolHolder)
    ssl.AddPatterns(lib, numThreads=-1)
    assert isinstance(lib.GetFpHolder(), ssl.PatternHolder)
    assert lib.GetFpHolder().GetFingerprint(0).GetNumBits() == 2048


def test_pattern_holder():
    for cls in (ssl.PatternHolder, ssl.TautomerPatternHolder):
        holder = cls()
        mols = ssl.CachedTrustedSmilesMolHolder()
        lib = ssl.SubstructLibrary(mols, holder)
        for smi in ("c1ccccc1", "CCO"):
            lib.AddMol(Chem.MolFromSmiles(smi))
        query = holder.MakeFingerprint(Chem.MolFromSmarts("c1ccccc1"))
        assert holder.PassesFilter(0, query)
        assert not holder.PassesFilter(1, query)
        fp = holder.GetFingerprint(0)
        assert holder.AddFingerprint(fp) == 2
        with pytest.raises(IndexError):
            holder.GetFingerprint(3)


def test_bad_smiles_never_match():
    for cls in (ssl.CachedSmilesMolHolder, ssl.CachedTrustedSmilesMolHolder):
        holder = cls()
        holder.AddSmiles("&CCS(=O)")
        holder.AddSmiles("CCO")
        lib = ssl.SubstructLibrary(holder)
        assert list(lib.GetMatches(Chem.MolFromSmiles("C"))) == [1]
        assert holder.GetMol(0) is None
        ssl.AddPatterns(lib, -1)
        assert len(lib.GetFpHolder()) == 2


def test_streams():
    lib = make_library(ssl.CachedSmilesMolHolder, ssl.PatternHolder)
    out = io.StringIO()
    lib.ToStream(out)
    lib2 = ssl.SubstructLibrary()
    lib2.InitFromStream(io.BytesIO(out.getvalue().encode()))
    assert len(lib2) == 8
    assert list(lib2.GetMatches(Chem.MolFromSmiles("CCO"))) == [1, 6, 7]
    assert ssl.SubstructLibraryCanSerialize()


def test_thread_counts_agree():
    lib = ssl.SubstructLibrary(ssl.CachedMolHolder(), ssl.PatternHolder())
    for smi in SMILES * 40:
        lib.AddMol(Chem.MolFromSmiles(smi))
    q = Chem.MolFromSmarts("[OX2H]")
    expected = list(lib.GetMatches(q, numThreads=1, maxResults=-1))
    for threads in (2, 3, 8):
        assert list(lib.GetMatches(q, numThreads=threads, maxResults=-1)) == expected
        assert list(lib.GetMatches(q, numThreads=threads, maxResults=7)) == expected[:7]


def test_errors():
    lib = make_library(ssl.MolHolder)
    with pytest.raises(IndexError):
        lib.GetMol(8)
    with pytest.raises(TypeError):
        lib.GetMatches("CC")
    with pytest.raises(TypeError):
        ssl.SubstructLibrary(ssl.MolHolder(), 5)
