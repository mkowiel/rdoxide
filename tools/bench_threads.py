#!/usr/bin/env python3
"""Times a few heavy operations with one thread and with all cores.

    python tools/bench_threads.py rdoxide   # sweeps rdBase.SetMaxThreads
    python tools/bench_threads.py rdkit     # numThreads where RDKit has it
Prints {operation: {"1": ms, "all": ms}} as JSON.
"""

import json
import os
import random
import sys
import time

lib = sys.argv[1] if len(sys.argv) > 1 else "rdoxide"
if lib == "rdkit":
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdDistGeom
    rdBase = None
else:
    from rdoxide import Chem, DataStructs, rdBase
    from rdoxide.Chem import rdDistGeom

chain = Chem.MolFromSmiles("C" * 1500)
weighted = Chem.MolFromSmiles("c1ccccc1" + "C" * 700 + "C=O")
steroid = Chem.AddHs(Chem.MolFromSmiles("C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O"))
rng = random.Random(1)
bvs = []
for _ in range(100000):
    bv = DataStructs.ExplicitBitVect(2048)
    for b in rng.sample(range(2048), 60):
        bv.SetBit(b)
    bvs.append(bv)


def embed(threads):
    p = rdDistGeom.ETKDGv3()
    p.randomSeed = 42
    p.numThreads = threads
    rdDistGeom.EmbedMultipleConfs(Chem.Mol(steroid), 80, p)


OPS = {
    "topological distance matrix, 1500-atom chain": lambda t: Chem.GetDistanceMatrix(chain, force=True),
    "bond-order distance matrix, 700 atoms": lambda t: Chem.GetDistanceMatrix(weighted, useBO=True, force=True),
    "bulk Tanimoto, 1 vs 100000": lambda t: DataStructs.BulkTanimotoSimilarity(bvs[0], bvs),
    "Tanimoto neighbours, 50 x 100000": lambda t: DataStructs.TanimotoSimilarityNeighbors(bvs[:50], bvs),
    "ETKDG embedding, 80 conformers": embed,
}


def best_ms(fn, repeats=3):
    best = float("inf")
    for _ in range(repeats):
        t = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t)
    return best * 1e3


out = {"cores": os.cpu_count()}
for name, op in OPS.items():
    row = {}
    for label, threads in (("1", 1), ("all", 0)):
        if rdBase is not None:
            rdBase.SetMaxThreads(threads)  # 0 = all cores
        # RDKit only parallelizes embedding (numThreads); the rest stays serial
        row[label] = best_ms(lambda: op(threads))
    out[name] = row
print(json.dumps(out))
