#!/usr/bin/env python3
"""Builds tests/data/coordgen_frag.txt: molecules (RDKit 2D coordinates, some atoms fixed or constrained) and the fragments
CoordgenFragmenter of coordgenlibs 3.0.2 makes of them (driver tools/coordgen_drv/frag.cpp).

    .venv-ref/bin/python tools/make_coordgen_frag_corpus.py DRIVER > tests/data/coordgen_frag.txt
"""
import importlib.util
import os
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("m", os.path.join(HERE, "make_coordgen_model_corpus.py"))
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def flags(k, n):
    fixed, constrained = set(), set()
    if k % 4 == 1:
        fixed = {0}
    elif k % 4 == 2:
        constrained = {i for i in range(n) if i % 5 == 0 and i < max(1, n // 3)}
    elif k % 7 == 3:
        fixed = {n - 1}
        constrained = {0}
    return fixed, constrained


CHAINS = (
    ["c1ccccc1" + "C" * n for n in range(2, 22, 2)]
    + ["C" * n for n in (4, 7, 10, 14)]
    + ["c1ccccc1" + "C" * n for n in (4, 5)]
    + ["c1ccc2ccccc2c1" + "C" * n for n in (7, 8)]
    + ["c1ccc2cc3ccccc3cc2c1" + "C" * n for n in (9, 10)]
    + ["c1ccccc1" + "C" * n + "c1ccccc1" for n in (3, 6, 9, 13)]
    + ["C1CCCCC1" + "CO" * n for n in (2, 4, 6)]
    + ["c1ccc2ccccc2c1" + "C" * n for n in (5, 8, 11)]
    + ["c1ccc2cc3ccccc3cc2c1" + "C" * n for n in (6, 9, 13, 15)]
    + ["CC(C)CC(C)CC(C)(C)CCCC", "CCC(CCCCCC)(CCCCC)CCCC", "CCCCCCCCCCCC=CCCCCC", "CCCCCCCC#CCCCCCC"]
)


def main():
    driver = sys.argv[1]
    blocks = []
    n_plain = len(base.smiles())
    for k, s in enumerate(base.smiles() + CHAINS):
        full = Chem.MolFromSmiles(s)
        if full is None or full.GetNumAtoms() == 0 or full.GetNumAtoms() > 60:
            continue
        # the sketcher works on one connected molecule at a time
        for m in Chem.GetMolFrags(full, asMols=True):
            Chem.Kekulize(m, clearAromaticFlags=True)
            AllChem.Compute2DCoords(m)
            fixed, constrained = flags(k, m.GetNumAtoms()) if k < n_plain else (set(), set())
            lines = base.block(m).split("\n")
            for i in range(m.GetNumAtoms()):
                f = ("f" if i in fixed else "") + ("c" if i in constrained else "")
                if f:
                    lines[1 + i] += " " + f
            blocks.append("\n".join(lines))
    inp = "".join(blocks)
    out = subprocess.run([driver], input=inp.encode(), capture_output=True, timeout=600).stdout.decode()
    recs = out.split("$$\n")[:-1]
    assert len(recs) == len(blocks), (len(recs), len(blocks))
    for b, r in zip(blocks, recs):
        sys.stdout.write(b + r + "$$\n")


if __name__ == "__main__":
    main()
