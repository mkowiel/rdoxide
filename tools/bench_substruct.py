#!/usr/bin/env python3
"""Microbenchmark of HasSubstructMatch / GetSubstructMatches.

Run it with the interpreter of the library to time:
    python tools/bench_substruct.py rdoxide
    python tools/bench_substruct.py rdkit      # original RDKit, other venv
It prints one JSON object {query: milliseconds per 1000 molecule checks}.
"""

import json
import sys
import time
from pathlib import Path

lib = sys.argv[1] if len(sys.argv) > 1 else "rdoxide"
if lib == "rdkit":
    from rdkit import Chem
else:
    from rdoxide import Chem

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "upstream" / "rdkit" / "Docs" / "Book" / "data"
mols = []
for name in ("solubility.train.sdf", "bzr.sdf", "actives_5ht3.sdf"):
    mols += [m for m in Chem.SDMolSupplier(str(DATA / name)) if m is not None]
mols = mols[:1500]

QUERIES = {
    "benzene c1ccccc1": "c1ccccc1",
    "carboxylic acid C(=O)[OH]": "C(=O)[OH]",
    "amide C(=O)N": "C(=O)N",
    "primary/secondary amine": "[NX3;H2,H1;!$(NC=O)]",
    "phenol [OH]c": "[OH]c",
    "naphthalene": "c1ccc2ccccc2c1",
    "6-atom chain C~C~C~C~C~C": "[#6]~[#6]~[#6]~[#6]~[#6]~[#6]",
    "recursive [$(C=O)]N": "[$(C=O)]N",
    "anilide CC(=O)Nc1ccccc1": "CC(=O)Nc1ccccc1",
    "no-match dibenzofuran": "c1ccc2c(c1)oc1ccccc12",
    "any atom *": "*",
    "halogen-substituted aryl c[F,Cl,Br,I]": "c[F,Cl,Br,I]",
}


def best_ms(fn, repeats=5):
    best = float("inf")
    for _ in range(repeats):
        t = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t)
    return best * 1e3


out = {"molecules": len(mols)}
for label, text in QUERIES.items():
    q = Chem.MolFromSmarts(text)
    if q is None:
        q = Chem.MolFromSmiles(text)
    has = lambda: [m.HasSubstructMatch(q) for m in mols]
    allm = lambda: [m.GetSubstructMatches(q, uniquify=True, maxMatches=1000) for m in mols]
    h = best_ms(has) / len(mols) * 1000
    a = best_ms(allm) / len(mols) * 1000
    out[label] = {"has": h, "matches": a, "hits": sum(has())}
print(json.dumps(out))
