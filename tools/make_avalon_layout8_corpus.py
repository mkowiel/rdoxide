#!/usr/bin/env python3
"""Builds tests/data/avalon_layout8.jsonl: molfiles with shortcut (R + atext) atoms and what Avalon's layout.c
IsShortcutLikeAtom, PerceiveShortcutClasses and LayoutShortcutChains produce (driver: .work/avalon_drv/lay8.c).

    .venv-ref/bin/python tools/make_avalon_layout8_corpus.py DRIVER > tests/data/avalon_layout8.jsonl
"""
import json
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

LABELS = ["Ph", "Me", "Et", "Boc", "OMe", "Ac", "Bn", "tBu", "CF3", "R1", "(CH2)3", "Ts"]
SMILES = [
    "*-*-*-*-*", "*-*(-*)-*", "C*-*C", "*C*C*C*", "c1ccccc1-*-*", "*-*-*-*-*-*-*-*-*-*", "C1CC(*)C(*)C1(*)*", "*-*-*-*-*-*", "CC(*)(*)C(*)*", "*1-*-*-*-*-1",
    "*-*-*-*(-*-*)-*-*", "C-*-*-C-*-*-C", "N-*-*", "O*-*", "*(*)(*)*", "F-*-*-Cl", "CC*", "*.*-*", "C1CCC1-*-*-*", "*-*-*-*-*-*-*-*", "C(*-*)(*-*)*-*",
    "*-*-C-*-*-C-*-*", "c1ccc(cc1)-*-*-*-c1ccccc1",
]


def block(smi, rot):
    m = Chem.MolFromSmiles(smi)
    AllChem.Compute2DCoords(m)
    lines = Chem.MolToMolBlock(m).split("\n")
    out, k = [], rot
    alias = []
    for i, a in enumerate(m.GetAtoms()):
        if a.GetAtomicNum() == 0:
            alias.append("A  %3d\n%s" % (i + 1, LABELS[k % len(LABELS)]))
            k += 1
    end = lines.index("M  END")
    return "\n".join(lines[:end] + alias + lines[end:])


def main():
    recs = []
    for r in range(3):
        for smi in SMILES:
            recs.append((smi + " #%d" % r, block(smi, r * 5)))
    out = subprocess.run([sys.argv[1]], input=("$$$$\n".join(b for _, b in recs) + "\n").encode(), capture_output=True, timeout=600).stdout.decode()
    parts = out.split("@@\n")[1:]
    assert len(parts) == len(recs), (len(parts), len(recs))
    for (s, b), o in zip(recs, parts):
        print(json.dumps({"input": b, "name": s, "output": o}))


if __name__ == "__main__":
    main()
