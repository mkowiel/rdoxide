#!/usr/bin/env python3
"""Builds tests/data/avalon_layout10.jsonl: SMILES (read with SMIToMOL) and molfiles with shortcut (R + atext)
atoms and what Avalon's LayoutMolecule makes of them (driver: .work/avalon_drv/lay10.c, -ffp-contract=off).

    .venv-ref/bin/python tools/make_avalon_layout10_corpus.py DRIVER > tests/data/avalon_layout10.jsonl
"""
import importlib.util
import json
import os
import subprocess
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    smis = []
    for name in ("make_avalon_layout1_corpus", "make_avalon_layout5_corpus", "make_avalon_layout6_corpus", "make_avalon_layout7_corpus"):
        smis += load(name).SMILES
    smis += ["CC(=O)Nc1ccc(O)cc1", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "O=C1NC(=O)C=C1", "C1CCC(=O)NCCCCCCC1", "C1CC(=O)NCCCCCCCCCCC1", "C1CCCCCCCCCCCCNC(=O)1", "OC(=O)CCCCCC(=O)NC1CCCCC1"]
    seen, uniq = set(), []
    for s in smis:
        if s not in seen:
            seen.add(s)
            uniq.append(s)
    recs = [("SMILES:" + s, s) for s in uniq]
    l8 = load("make_avalon_layout8_corpus")
    for r in range(2):
        for smi in l8.SMILES:
            recs.append((l8.block(smi, r * 5), smi + " #%d" % r))
    out = subprocess.run([sys.argv[1]], input=("$$$$\n".join(b for b, _ in recs) + "\n").encode(), capture_output=True, timeout=900).stdout.decode()
    parts = out.split("@@\n")[1:]
    assert len(parts) == len(recs), (len(parts), len(recs))
    for (b, n), o in zip(recs, parts):
        print(json.dumps({"input": b, "name": n, "output": o}))


if __name__ == "__main__":
    main()
