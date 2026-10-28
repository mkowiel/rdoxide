#!/usr/bin/env python3
"""Builds tests/data/avalon_geometry.jsonl: command lines for the functions of Avalon's geometry.c
(TP TransformPoints, PSM PointSetMatchTransformation, NGN NextGraphNeighbour, CUSP HasCusp,
PIR PointsIntoRing, NSP NextSubstituentPoint, TGS TransformGraphSubstituent) and the output of the C
code (driver: .work/avalon_drv/geom.c, built with -ffp-contract=off).

    .venv-ref/bin/python tools/make_avalon_geometry_corpus.py DRIVER > tests/data/avalon_geometry.jsonl
"""
import json
import random
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
rng = random.Random(21)

SMILES = ["CCO", "c1ccccc1O", "C[C@H](N)C(=O)O", "C1CC1", "CC(C)(C)C", "c1ccc2ccccc2c1", "CC(=O)Oc1ccccc1C(=O)O", "C12C3C4C1C5C2C3C45",
          "OC[C@@H](O)[C@H](O)CO", "CN1CCC[C@H]1c1cccnc1", "C=C=C", "FC(F)(F)C(F)(F)F", "C1CCCCC1", "CC12CCC3C(CCc4cc(O)ccc34)C1CCC2O", "[Na+].[Cl-]"]


def num(x):
    return repr(float(x))


def pts(p):
    return " ".join("%s %s" % (num(a), num(b)) for a, b in p)


def graph(smi, scale=1.5):
    m = Chem.MolFromSmiles(smi)
    AllChem.Compute2DCoords(m)
    c = m.GetConformer()
    p = [(c.GetAtomPosition(i).x * scale / 1.5, c.GetAtomPosition(i).y * scale / 1.5) for i in range(m.GetNumAtoms())]
    e = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()]
    return p, e


def edges_s(e):
    return " ".join("%d %d" % x for x in e)


def cases():
    graphs = [graph(s) for s in SMILES]
    for _ in range(150):
        n = rng.randint(1, 12)
        p = [(rng.uniform(-5, 5), rng.uniform(-5, 5)) for _ in range(n)]
        a = [(rng.uniform(-3, 3), rng.uniform(-3, 3)) for _ in range(4)]
        yield "TP %d %s %s" % (n, pts(a), pts(p)), n
        yield "TP %d %s %s" % (n, pts([a[0], a[0], a[2], a[3]]), pts(p)), n
        nm = rng.randint(1, n)
        yield "PSM %d %d %d %s %s %s" % (n, nm, rng.randint(0, 1), pts(p), pts(p[:nm]), pts([(rng.uniform(-5, 5), rng.uniform(-5, 5)) for _ in range(nm)])), n
        yield "CUSP %s %d %s" % (pts([(0, 0)]), n, pts([(rng.uniform(-2, 2), rng.uniform(-2, 2)) for _ in range(n)])), 1
    for p, e in graphs:
        n = len(p)
        for node in range(n):
            for last in (-1, 0, 2):
                yield "NGN %d %d %d %s" % (node, last, len(e), edges_s(e)), 1
        for seed in range(n):
            for flags in (0, 1, 2, 3, 4, 5):
                yield "NSP %d %d %s %d %s %d %s" % (seed, flags, pts([(-99.0, -99.0)]), n, pts(p), len(e), edges_s(e)), 1
            for _ in range(2):
                q = (p[seed][0] + rng.uniform(-2, 2), p[seed][1] + rng.uniform(-2, 2))
                yield "PIR %d %s %d %s %d %s" % (seed, pts([q]), n, pts(p), len(e), edges_s(e)), 1
            nb = [(p[b][0], p[b][1]) for a, b in e if a == seed] + [(p[a][0], p[a][1]) for a, b in e if b == seed]
            if nb:
                yield "CUSP %s %d %s" % (pts([p[seed]]), len(nb), pts(nb)), 1
    for _ in range(60):
        (sp, se), (mp, me) = rng.choice(graphs), rng.choice(graphs)
        apo, seed = rng.randrange(len(sp)), rng.randrange(len(mp))
        sp = [(x + rng.uniform(-0.3, 0.3), y) for x, y in sp]
        yield "TGS %d %s %d %s %d %d %s %d %s %d" % (len(sp), pts(sp), len(se), edges_s(se), apo, len(mp), pts(mp), len(me), edges_s(me), seed), len(sp)
    # single atom cases
    yield "TGS 1 %s 0  0 3 %s 2 0 1 1 2 1" % (pts([(0.3, 0.4)]), pts([(0, 0), (1.5, 0), (2, 1)])), 1


def main():
    cs = list(cases())
    inp = "\n".join(c for c, _ in cs) + "\n"
    out = subprocess.run([sys.argv[1]], input=inp.encode(), capture_output=True, timeout=300).stdout.decode().split("\n")
    pos = 0
    for c, n in cs:
        print(json.dumps({"cmd": c, "output": "\n".join(out[pos : pos + n])}))
        pos += n


if __name__ == "__main__":
    main()
