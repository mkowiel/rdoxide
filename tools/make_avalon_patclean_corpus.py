#!/usr/bin/env python3
"""Builds tests/data/avalon_patclean.jsonl and the template files in tests/data/avalon_patterns/: molecules sent
through Avalon's struchk with template cleaning (-pc), rotation (-pr) and forced stereochemistry (-ps) and what the
C code (patclean.c via RunStruchk) answers (driver: .work/avalon_drv/sc.c). `{P}` in the options stands for the
directory with the template files.

    .venv-ref/bin/python tools/make_avalon_patclean_corpus.py DRIVER PATTERNDIR > tests/data/avalon_patclean.jsonl
"""
import json
import math
import os
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

CLEAN = ["c1ccccc1", "c1ccc2ccccc2c1", "c1ccncc1", "C1CCCCC1", "C(=O)N", "c1ccc(cc1)C(=O)O", "C1CC1", "c1ccc2[nH]ccc2c1",
         "CC(C)(C)c1ccccc1", "N1CCOCC1", "C1=CC=CC=C1", "O=C1CCCCC1"]
ROTATE = ["c1ccccc1", "c1ccc2ccccc2c1", "C1CCCCC1", "c1ccc(cc1)c1ccccc1", "CCCC", "C=CC", "c1ccncc1C"]
STEREO = ["C[C@H](N)C(=O)O", "C[C@@H](O)CC", "C[C@H]1CCCC[C@@H]1C", "N[C@@H](Cc1ccccc1)C(=O)O", "C[C@H](F)Cl", "O[C@H]1CCCC[C@@H]1O"]

SMILES = [
    "c1ccccc1", "Cc1ccccc1", "c1ccccc1C(=O)O", "c1ccc2ccccc2c1", "Cc1cccc2ccccc12", "c1ccncc1", "Cc1ccncc1", "CCc1ccncc1CC",
    "C1CCCCC1", "CC1CCCCC1C", "C1CC1C", "CC(=O)N", "CCC(=O)NC", "CC(C)(C)c1ccccc1", "O=C1CCCCC1", "CC1CCCCC1=O", "N1CCOCC1",
    "CN1CCOCC1", "c1ccc2[nH]ccc2c1", "Cc1c[nH]c2ccccc12", "CC(N)C(=O)O", "CCC(N)C(=O)O", "CC(O)CC", "CC(F)Cl", "NC(Cc1ccccc1)C(=O)O",
    "CC1CCCCC1C", "OC1CCCCC1O", "CC(N)C(=O)NC(C)C(=O)O", "c1ccc(cc1)c1ccccc1", "Cc1ccc(cc1)c1ccccc1C", "CCCC", "CCCCC", "C=CC",
    "CC=CC", "C(=O)(N)c1ccccc1", "Clc1ccc(cc1)C(=O)O", "NC(=O)c1ccncc1", "C1CCC2CCCCC2C1", "c1ccc(cc1)C1CCCCC1",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "CC[C@H](C)C(=O)O", "C[C@@H](O)c1ccccc1", "C1CCCCC1C(=O)N", "C1CC1c1ccccc1",
    "OC(=O)C1CCCCC1N", "CC(C)(C)c1ccccc1CC(N)C(=O)O", "c1ccccc1.CC(N)C(=O)O", "c1ccccc1.c1ccncc1", "CC(N)C(=O)O.CC(O)CC",
]

OPTS = [
    "-pc {P}/clean.sdf", "-pr {P}/rotate.sdf", "-ps {P}/stereo.sdf", "-pc {P}/clean.sdf|-ps {P}/stereo.sdf",
    "-pr {P}/rotate.sdf|-pc {P}/clean.sdf", "-ps {P}/stereo.sdf|-pr {P}/rotate.sdf|-pc {P}/clean.sdf|-cs",
]


def template_file(path, smis, stereo=False):
    out = []
    for s in smis:
        m = Chem.MolFromSmiles(s)
        AllChem.Compute2DCoords(m)
        if stereo:
            Chem.WedgeMolBonds(m, m.GetConformer())
        m.SetProp("_Name", s)
        out.append(Chem.MolToMolBlock(m) + "$$$$\n")
    open(path, "w").write("".join(out))


def transformed(m, angle, flip):
    c = m.GetConformer()
    ca, sa = math.cos(angle), math.sin(angle)
    for i in range(m.GetNumAtoms()):
        p = c.GetAtomPosition(i)
        x = -p.x if flip else p.x
        c.SetAtomPosition(i, (ca * x - sa * p.y + 3.0, sa * x + ca * p.y - 2.0, 0))


def main():
    driver, pdir = sys.argv[1], sys.argv[2]
    os.makedirs(pdir, exist_ok=True)
    template_file(os.path.join(pdir, "clean.sdf"), CLEAN)
    template_file(os.path.join(pdir, "rotate.sdf"), ROTATE)
    template_file(os.path.join(pdir, "stereo.sdf"), STEREO, True)
    inputs = [("SMILES:" + s, s) for s in SMILES]
    for k, s in enumerate(SMILES):
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        AllChem.Compute2DCoords(m)
        transformed(m, 0.4 + 0.7 * k, k % 2 == 1)
        inputs.append((Chem.MolToMolBlock(m), "2d " + s))
    recs = []
    for k, (body, name) in enumerate(inputs):
        for j in (k % len(OPTS), (k + 3) % len(OPTS)):
            recs.append((OPTS[j], body, name))
    # one process per record: the C code crashes when it frees the templates of a previous option set
    for o, b, n in recs:
        inp = o.replace("{P}", pdir) + "\n" + b + ("\n" if not b.endswith("\n") else "")
        out = subprocess.run([driver], input=inp.encode(), capture_output=True, timeout=600).stdout.decode()
        assert out.startswith("@@\n"), (o, n)
        print(json.dumps({"options": o, "input": b, "name": n, "output": out[3:]}))


if __name__ == "__main__":
    main()
