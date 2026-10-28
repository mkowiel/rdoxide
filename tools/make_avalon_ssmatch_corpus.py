#!/usr/bin/env python3
"""Builds tests/data/avalon_ssmatch.jsonl: molecule / query pairs (molfiles) and what Avalon's ssmatch.c
(SSMatch, PruneDuplicateAtomSets, RSSMatch, FGMatch, Refine) finds (driver: .work/avalon_drv/ss.c).

    .venv-ref/bin/python tools/make_avalon_ssmatch_corpus.py DRIVER > tests/data/avalon_ssmatch.jsonl
"""
import importlib.util
import json
import os
import random
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


QUERIES = ["c1ccccc1", "C(=O)O", "[N+](=O)[O-]", "C=C", "CC", "N", "C#N", "OC(=O)C", "C1CC1", "c1ccncc1", "C(=O)N", "CC(C)C", "CO", "C1CCCCC1", "S(=O)(=O)", "ClC", "c1ccc2ccccc2c1", "CN(C)C", "O=CC=C", "NC(=O)C"]


def molblock(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    AllChem.Compute2DCoords(m)
    return Chem.MolToMolBlock(m)


def mutate(mb, rng, kind):
    """query variants: atoms to A/Q, bonds to the query bond types 5..8"""
    lines = mb.split("\n")
    na, nb = int(lines[3][0:3]), int(lines[3][3:6])
    for i in range(na):
        ln = 4 + i
        if kind in (1, 3) and rng.random() < 0.4:
            lines[ln] = lines[ln][:31] + ("A" if rng.random() < 0.5 else "Q") + "  " + lines[ln][34:]
    for i in range(nb):
        ln = 4 + na + i
        if kind in (2, 3) and rng.random() < 0.5:
            lines[ln] = lines[ln][:6] + "%3d" % rng.choice([5, 6, 7, 8]) + lines[ln][9:]
    return "\n".join(lines)


def with_mapping(mb):
    """atom-atom mapping numbers for the non-terminal atoms (terminal ones keep 0 and must not have extra neighbours)"""
    lines = mb.split("\n")
    na, nb = int(lines[3][0:3]), int(lines[3][3:6])
    deg = [0] * na
    for i in range(nb):
        ln = lines[4 + na + i]
        deg[int(ln[0:3]) - 1] += 1
        deg[int(ln[3:6]) - 1] += 1
    for i in range(na):
        if deg[i] > 1:
            ln = lines[4 + i]
            lines[4 + i] = ln[:60] + "%3d" % (i + 1) + ln[63:]
    return "\n".join(lines)


def main():
    smis = []
    for name in ("make_avalon_layout1_corpus", "make_avalon_layout5_corpus", "make_avalon_layout6_corpus"):
        smis += load(name).SMILES
    seen, targets = set(), []
    for s in smis:
        if s not in seen and molblock(s):
            seen.add(s)
            targets.append(s)
    targets = targets[::2] + ["CC(=O)Oc1ccccc1C(=O)O", "O=C(O)c1ccccc1", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "c1ccc2ccccc2c1", "CCOC(=O)C=C", "NC(=O)CC(=O)N"]
    rng = random.Random(7)
    recs = []
    for i, t in enumerate(targets):
        tb = molblock(t)
        for k in range(5):
            q = QUERIES[(i * 3 + k * 7) % len(QUERIES)]
            qb = molblock(q)
            kind = k % 4
            if kind:
                qb = mutate(qb, rng, kind)
            else:
                qb = with_mapping(qb)
            recs.append((t + " | " + q + " #%d" % kind, tb, qb))
    inp = "$$$$\n".join(r[1] + "$$$$\n" + r[2] for r in recs) + "\n"
    out = subprocess.run([sys.argv[1]], input=inp.encode(), capture_output=True, timeout=1800).stdout.decode()
    parts = out.split("@@\n")[1:]
    assert len(parts) == len(recs), (len(parts), len(recs))
    for (n, tb, qb), o in zip(recs, parts):
        print(json.dumps({"name": n, "mol": tb, "query": qb, "output": o}))


if __name__ == "__main__":
    main()
