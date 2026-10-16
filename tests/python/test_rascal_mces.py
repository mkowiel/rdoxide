"""rdRascalMCES: API surface not covered by the upstream tests."""
import pytest
from rdkit import Chem
from rdkit.Chem import rdRascalMCES


def test_result_properties():
    m1 = Chem.MolFromSmiles("c1ccccc1Cl")
    m2 = Chem.MolFromSmiles("c1ccccc1F")
    (res,) = rdRascalMCES.FindMCES(m1, m2)
    assert res.numFragments == 1
    assert res.largestFragmentSize == 6
    assert res.timedOut is False
    assert 0.0 < res.similarity <= 1.0
    assert res.tier1Sim >= 0.7 and res.tier2Sim >= 0.7


def test_options_are_checked():
    opts = rdRascalMCES.RascalOptions()
    with pytest.raises(AttributeError):
        opts.noSuchOption = 1
    opts.timeout = 5
    assert opts.timeout == 5


def test_cluster_option_b_is_upstreams_a():
    # upstream binds "b" to the member "a"
    co = rdRascalMCES.RascalClusterOptions()
    co.b = 0.5
    assert co.a == 0.5


def test_cluster_rejects_none():
    with pytest.raises(ValueError):
        rdRascalMCES.RascalCluster([Chem.MolFromSmiles("CC"), None])


def test_cluster_is_independent_of_thread_count():
    mols = [Chem.MolFromSmiles(s) for s in
            ["c1ccccc1CCN", "c1ccccc1CCO", "c1ccccc1CCCN", "CCCCCC", "CCCCCCC", "c1ccncc1CCN"]]
    results = []
    for n in (1, 2, 4):
        co = rdRascalMCES.RascalClusterOptions()
        co.numThreads = n
        co.similarityCutoff = 0.5
        results.append(rdRascalMCES.RascalCluster(mols, co))
    assert results[0] == results[1] == results[2]
