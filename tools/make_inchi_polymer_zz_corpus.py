#!/usr/bin/env python3
"""Builds tests/data/inchi_polymers_zz.jsonl: polymer units capped by two Zz atoms (the input that
needs CRU folding and frame shift) with the results of the official RDKit 2026.09.1.

    .venv-ref/bin/python tools/make_inchi_polymer_zz_corpus.py > tests/data/inchi_polymers_zz.jsonl
"""
import json
import random

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdinchi

OPTIONS = ["", "-Polymers", "-Polymers -NPZz", "-Polymers -FoldCRU", "-Polymers105"]
BASE = [
    "CCOC", "CC(C)C", "CC(c1ccccc1)C", "C(=O)CCCCC(=O)NCCCCCCN", "CC(F)(F)C", "CCCC", "OCCO",
    "C[Si](C)(C)OC", "CC(O)CC(O)C", "CCCCCCCC", "c1ccc(cc1)C", "NCC(=O)N", "CC(=O)OC", "C1CCCCC1",
    "OCCOCCOCCOCCO", "NCCNCCN", "CC(C)CC(C)C", "C(C)C(C)", "CC(Cl)CC(Cl)C", "OC(=O)CCC(=O)OCC",
    "CC(C)(C)CC", "N1CCNCC1",
]


def build(smi, rng, stype, con, lab):
    m = Chem.MolFromSmiles("*" + smi + "*")
    rdDepictor.Compute2DCoords(m)
    n = m.GetNumAtoms()
    atoms = list(range(1, n - 1))
    inside = set(atoms)
    xb = [b.GetIdx() for b in m.GetBonds() if (b.GetBeginAtomIdx() in inside) != (b.GetEndAtomIdx() in inside)]
    sg = Chem.CreateMolSubstanceGroup(m, stype)
    for a in atoms:
        sg.AddAtomWithIdx(a)
    for b in xb:
        sg.AddBondWithIdx(b)
    sg.SetProp("LABEL", lab)
    sg.SetProp("CONNECT", con)
    mb = Chem.MolToMolBlock(m)
    return mb.replace(" R# ", " Zz ").replace(" R  ", " Zz ").replace(" *  ", " Zz ")


def one(mb):
    res = {}
    for opt in OPTIONS:
        r = rdinchi.MolBlockToInchi(mb, opt)
        res[opt] = [r[0], r[1], r[2], r[3], r[4]]
    return res


def main():
    import subprocess
    import sys
    if len(sys.argv) > 1:  # child: one molblock on stdin
        print(json.dumps(one(sys.stdin.read())))
        return
    rng = random.Random(7)
    for smi in BASE:
        for stype in ("SRU", "CRU", "MON"):
            for con in ("HT", "HH", "EU"):
                mb = build(smi, rng, stype, con, rng.choice(["n", "3", "2-5"]))
                p = subprocess.run([sys.executable, __file__, "child"], input=mb, capture_output=True, text=True)
                if p.returncode != 0:  # the reference InChI library crashes on some input: skip
                    continue
                print(json.dumps({"molblock": mb, "results": json.loads(p.stdout)}))


if __name__ == "__main__":
    main()
