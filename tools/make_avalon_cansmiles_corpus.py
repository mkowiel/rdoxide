#!/usr/bin/env python3
"""Builds tests/data/avalon_cansmiles.jsonl: CanSmiles and RemoveStandardComponents as run by the C code
(driver: .work/avalon_drv/cs.c, -ffp-contract=off, glibc-style rand).

    .venv-ref/bin/python tools/make_avalon_cansmiles_corpus.py DRIVER > tests/data/avalon_cansmiles.jsonl
"""
import importlib.util
import json
import os
import subprocess
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
DB, CENTER, SCRAMBLE, COMPSET, MAIN, OCTET, ISO, MAP = 1, 2, 4, 8, 0x10, 0x20, 0x40, 0x80
TO_MONOMER, DY = 0x1000, 0x2000
FLAGS = [0, DB, DB | CENTER, DB | CENTER | MAIN | COMPSET, DB | CENTER | SCRAMBLE, OCTET | DB, ISO | MAP | CENTER, DY | DB | CENTER, MAIN | CENTER]


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


EXTRA = [
    "CC(=O)[O-].[Na+]", "[Na+].[Cl-]", "CCO.CCO", "c1ccccc1.c1ccccc1.O", "[2H]C([2H])([2H])O", "[13CH4]", "C[C@H](N)C(=O)O.Cl", "OC(=O)C(F)(F)F.CCN", "[Cu+2].[O-]S(=O)(=O)[O-]", "C[N+](C)(C)C.[Br-]",
    "[NH4+].[Cl-]", "C/C=C/C.C/C=C\\C", "F/C=C/F.F/C=C\\F", "[CH3:1][CH2:2]O", "CC(=O)OC1=CC=CC=C1C(O)=O", "O=[N+]([O-])c1ccccc1", "C[S+](C)[O-]", "[O-][N+](=O)C", "C=[N+]=[N-]", "CC#N", "N#N", "[H][H]",
    "O", "[Fe]", "*", "C*C", "{R1}C", "c1ccc2[nH]ccc2c1", "C1=CC=CC=C1", "c1cc[nH]c1", "O=c1cc[nH]cc1", "Oc1ccccn1", "C[C@@H]1CC[C@H](C)CC1", "N[C@@H](C)C(=O)N[C@@H](C)C(=O)O",
]
SALTS = ["[Na+]", "Cl", "[Cl-]", "O", "[K+]", "CCO", "Br"]
SALT_CASES = ["[Na+].CC(=O)[O-]", "CC(=O)[O-].[Na+]", "CCN.Cl.Cl", "Cl.CCN", "CCN.Cl", "O.O.CCN", "[Na+].[Cl-].[K+]", "CCN", "Cl.Cl", "Cl", "CCO.CCN.CCO", "CCN.[Na+].CCO.Br"]


def main():
    smis = list(EXTRA)
    for name in ("make_avalon_layout1_corpus", "make_avalon_layout5_corpus", "make_avalon_layout6_corpus", "make_avalon_layout7_corpus"):
        smis += load(name).SMILES
    seen, uniq = set(), []
    for s in smis:
        if s not in seen:
            seen.add(s)
            uniq.append(s)
    recs = []
    for i, s in enumerate(uniq):
        for k in range(2):
            recs.append(("%d" % FLAGS[(i + k * 4) % len(FLAGS)], s))
        m = Chem.MolFromSmiles(s)
        if m is not None:
            recs.append(("%d" % FLAGS[(i * 3) % len(FLAGS)], Chem.MolToSmiles(m, doRandom=True)))
    for s in uniq[:20]:
        recs.append(("%d" % TO_MONOMER, s))
    for s in SALT_CASES:
        recs.append(("S", s + "\t" + ",".join(SALTS)))
    inp = "".join("%s\t%s\n" % r for r in recs)
    out = subprocess.run([sys.argv[1]], input=inp.encode(), capture_output=True, timeout=1800).stdout.decode().split("\n")
    assert len(out) - 1 == len(recs), (len(out), len(recs))
    for r, o in zip(recs, out):
        print(json.dumps({"flags": r[0], "input": r[1], "output": o}))


if __name__ == "__main__":
    main()
