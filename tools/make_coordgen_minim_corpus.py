#!/usr/bin/env python3
"""Builds tests/data/coordgen_minim.txt: molecules with randomly displaced RDKit 2D coordinates and what
CoordgenMinimizer::minimizeMolecule of coordgenlibs 3.0.2 makes of them (driver tools/coordgen_drv/minim.cpp).

    .venv-ref/bin/python tools/make_coordgen_minim_corpus.py DRIVER > tests/data/coordgen_minim.txt
"""
import importlib.util
import os
import random
import struct
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("m", os.path.join(HERE, "make_coordgen_model_corpus.py"))
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

EXTRA = ["C/C=C/C", "C/C=C\\C", "C1CCCCCCCCCC1", "C1CCCCCCCCCCCCCCC1", "C/C=C/CCCC/C=C\\CCCC1", "C1=CCCCCCC/C=C/CCCCC1"]


def main():
    driver = sys.argv[1]
    rnd = random.Random(11)
    blocks = []
    for s in base.smiles() + EXTRA:
        full = Chem.MolFromSmiles(s)
        if full is None or full.GetNumAtoms() == 0 or full.GetNumAtoms() > 40:
            continue
        for m in Chem.GetMolFrags(full, asMols=True):
            if m.GetNumAtoms() < 2:
                continue
            Chem.Kekulize(m, clearAromaticFlags=True)
            AllChem.Compute2DCoords(m)
            lines = base.block(m).split("\n")
            for i in range(m.GetNumAtoms()):
                f = lines[1 + i].split(" ")
                x = struct.unpack("<f", struct.pack("<I", int(f[3], 16)))[0] + rnd.uniform(-8, 8)
                y = struct.unpack("<f", struct.pack("<I", int(f[4], 16)))[0] + rnd.uniform(-8, 8)
                f[3], f[4] = base.fb(x), base.fb(y)
                lines[1 + i] = " ".join(f)
            blocks.append("\n".join(lines))
    inp = "".join(b if b.endswith("\n") else b + "\n" for b in blocks)
    out = subprocess.run([driver], input=inp.encode(), capture_output=True, timeout=1200).stdout.decode()
    recs = out.split("$$\n")[:-1]
    assert len(recs) == len(blocks), (len(recs), len(blocks))
    for b, r in zip(blocks, recs):
        sys.stdout.write((b if b.endswith("\n") else b + "\n") + r + "$$\n")


if __name__ == "__main__":
    main()
