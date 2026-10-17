"""Thread control (RDOXIDE_NUM_THREADS, rdBase.SetMaxThreads) and the
guarantee that results do not depend on the number of threads."""

import os
import random
import subprocess
import sys

import pytest
from rdoxide import Chem, DataStructs, rdBase
from rdoxide.Chem import rdDistGeom


@pytest.fixture
def restore_threads():
    before = rdBase.GetMaxThreads()
    yield
    rdBase.SetMaxThreads(before)


def _max_threads_with_env(value):
    env = dict(os.environ)
    if value is None:
        env.pop("RDOXIDE_NUM_THREADS", None)
    else:
        env["RDOXIDE_NUM_THREADS"] = value
    out = subprocess.run(
        [sys.executable, "-c", "from rdoxide import rdBase; print(rdBase.GetMaxThreads())"],
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    return int(out.stdout)


def test_environment_variable_sets_the_ceiling():
    cores = os.cpu_count()
    assert _max_threads_with_env(None) == cores
    assert _max_threads_with_env("all") == cores
    assert _max_threads_with_env("0") == cores
    assert _max_threads_with_env("3") == 3
    assert _max_threads_with_env("-1") == max(cores - 1, 1)
    # unusable text falls back to all cores
    assert _max_threads_with_env("many") == cores


def test_set_max_threads(restore_threads):
    rdBase.SetMaxThreads(2)
    assert rdBase.GetMaxThreads() == 2
    rdBase.SetMaxThreads(0)
    assert rdBase.GetMaxThreads() == os.cpu_count()


def test_embed_parameters_read_one_but_run_on_all_cores():
    params = rdDistGeom.ETKDGv3()
    assert params.numThreads == 1  # upstream's visible default
    params.numThreads = 3
    assert params.numThreads == 3


def test_embedding_does_not_depend_on_the_thread_count(restore_threads):
    mol = Chem.AddHs(Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O"))
    positions = {}
    for threads in (1, 4):
        rdBase.SetMaxThreads(threads)
        params = rdDistGeom.ETKDGv3()
        params.randomSeed = 42
        m = Chem.Mol(mol)
        rdDistGeom.EmbedMultipleConfs(m, 12, params)
        positions[threads] = [c.GetPositions().tolist() for c in m.GetConformers()]
    assert positions[1] == positions[4]


def test_distance_matrix_does_not_depend_on_the_thread_count(restore_threads):
    mol = Chem.MolFromSmiles("C1CCCCC1" * 1 + "C" * 400)
    results = []
    for threads in (1, 6):
        rdBase.SetMaxThreads(threads)
        results.append(Chem.GetDistanceMatrix(mol, useBO=False, force=True).tolist())
    assert results[0] == results[1]
    n = mol.GetNumAtoms()
    assert results[0][0][n - 1] == pytest.approx(float(n - 6) + 1.0, abs=1)  # chain plus ring link


def test_weighted_distance_matrix_does_not_depend_on_the_thread_count(restore_threads):
    mol = Chem.MolFromSmiles("c1ccccc1" + "C" * 600 + "C=O")
    results = []
    for threads in (1, 6):
        rdBase.SetMaxThreads(threads)
        results.append(Chem.GetDistanceMatrix(mol, useBO=True, force=True).tolist())
    assert results[0] == results[1]


def test_bulk_similarity_does_not_depend_on_the_thread_count(restore_threads):
    rng = random.Random(7)
    bvs = []
    for _ in range(25000):
        bv = DataStructs.ExplicitBitVect(512)
        for b in rng.sample(range(512), 40):
            bv.SetBit(b)
        bvs.append(bv)
    query = bvs[0]
    results = []
    for threads in (1, 6):
        rdBase.SetMaxThreads(threads)
        results.append(
            (
                DataStructs.BulkTanimotoSimilarity(query, bvs),
                DataStructs.TanimotoSimilarityNeighbors(bvs[:3], bvs),
            )
        )
    assert results[0] == results[1]
    assert results[0][0][0] == 1.0
