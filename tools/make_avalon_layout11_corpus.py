#!/usr/bin/env python3
"""Builds tests/data/avalon_layout11.jsonl: SMIToMOL with DO_LAYOUT (and optional coordinate strings after '|')
as run by the C code (driver: .work/avalon_drv/smil.c, -ffp-contract=off).

    .venv-ref/bin/python tools/make_avalon_layout11_corpus.py DRIVER > tests/data/avalon_layout11.jsonl
"""
import importlib.util
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DO_LAYOUT, TRUE_DB_STEREO, DROP_TRIVIAL_HYDROGENS = 1, 0x80, 8


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


COORDS = [
    "CCO|0,0,1.5,0", "CCCC|0,0,1.5,0,3,0,4.5,0", "c1ccccc1C|0,0,1.3,0.75,1.3,2.25,0,3,-1.3,2.25,-1.3,0.75", "C1CCCCC1O|0,0,1.3,0.75,1.3,2.25,0,3,-1.3,2.25,-1.3,0.75",
    "CC(C)(C)O|,,1,2,,,3,4", "CCN|0,0", "OCCO|1,1,2,2,x,y,5,5", "C[C@H](N)C(=O)O|0,0,1.5,0,3,0,4.5,0,6,0", "C/C=C/C|0,0,1.3,0.75,2.6,0,3.9,0.75",
]


def main():
    smis = []
    for name in ("make_avalon_layout1_corpus", "make_avalon_layout5_corpus", "make_avalon_layout6_corpus"):
        smis += load(name).SMILES
    seen, uniq = set(), []
    for s in smis:
        if s not in seen:
            seen.add(s)
            uniq.append(s)
    recs = []
    for i, s in enumerate(uniq):
        recs.append((DO_LAYOUT | (TRUE_DB_STEREO if i % 2 else 0), s))
    for c in COORDS:
        recs.append((DO_LAYOUT | TRUE_DB_STEREO, c))
        recs.append((DO_LAYOUT | TRUE_DB_STEREO | DROP_TRIVIAL_HYDROGENS, c))
    out = subprocess.run([sys.argv[1]], input="".join("%d\t%s\n" % r for r in recs).encode(), capture_output=True, timeout=900).stdout.decode()
    parts = out.split("@@\n")[1:]
    assert len(parts) == len(recs), (len(parts), len(recs))
    for (f, s), o in zip(recs, parts):
        print(json.dumps({"flags": f, "input": s, "output": o}))


if __name__ == "__main__":
    main()
