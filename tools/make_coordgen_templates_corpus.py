#!/usr/bin/env python3
"""Builds tests/data/coordgen_templates.txt: the normalized templates of coordgenlibs 3.0.2 ("T" lines), then molecules
(the templates themselves with permuted atom order, and RDKit molecules) with the Morgan scores, the template match and an
alignment matrix the C++ code gives (driver tools/coordgen_drv/tmpl.cpp).

    .venv-ref/bin/python tools/make_coordgen_templates_corpus.py DRIVER > tests/data/coordgen_templates.txt
"""
import importlib.util
import math
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
TEMPLATES = os.path.join(HERE, "..", "crates", "rdoxide", "src", "external", "coordgen", "data", "templates.txt")

EXTRA = [
    "C1CC2CCC1C2", "C12C3C4C1C5C2C3C45", "C1C2CC3CC1CC(C2)C3", "c1ccc2c(c1)ccc1ccccc12", "C1CCC2C(C1)CCC1C2CCC2CCCC21",
    "c1ccc2cc3ccccc3cc2c1", "C1CC2CCC1CC2", "C1CC1C1CC1", "C1CCC2(CC1)CCCC2", "c1cc2cccc3c4cccc5cccc(c(c1)c23)c54", "C1CCC2C(C1)C1CCC2C1",
    "C1C2CC3C1C2C3", "C12CC3CC(CC(C3)C1)C2", "c1ccc2c(c1)c1ccccc1c1ccccc21",
]


def fb(x):
    return "%08x" % struct.unpack("<I", struct.pack("<f", x))[0]


def template_blocks():
    lines = open(TEMPLATES).read().split("\n")
    i, out = 0, []
    rnd = random.Random(7)
    while i < len(lines) and lines[i].startswith("T "):
        na, nb = map(int, lines[i].split()[1:])
        atoms = [l.split() for l in lines[i + 1 : i + 1 + na]]
        bonds = [l.split() for l in lines[i + 1 + na : i + 1 + na + nb]]
        i += 1 + na + nb
        perm = list(range(na))
        rnd.shuffle(perm)
        ang = rnd.uniform(0, 6.28)
        flip = rnd.random() < 0.5
        inv = {old: new for new, old in enumerate(perm)}
        rows = []
        for new, old in enumerate(perm):
            z, x, y = atoms[old]
            x, y = float(x) * 50, float(y) * 50
            if flip:
                x = -x
            x, y = x * math.cos(ang) - y * math.sin(ang), x * math.sin(ang) + y * math.cos(ang)
            rows.append("A %s 0 %s %s" % (z, fb(x), fb(y)))
        brows = ["B %d %d %s" % (inv[int(a)], inv[int(b)], o) for a, b, o in bonds]
        # macrocycle double bonds need the CIP ranking, which is not ported yet
        rw = Chem.RWMol()
        for z, _, _ in atoms:
            rw.AddAtom(Chem.Atom(int(z)))
        for a, b, o in bonds:
            rw.AddBond(int(a), int(b), {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}[int(o)])
        Chem.FastFindRings(rw)
        ri = rw.GetRingInfo()
        if any(b.GetBondType() == Chem.BondType.DOUBLE and ri.MinBondRingSize(b.GetIdx()) >= 9 for b in rw.GetBonds()):
            continue
        out.append("MOL %d %d\n" % (na, nb) + "\n".join(rows + brows) + "\nEND\n")
    return out


def main():
    driver = sys.argv[1]
    dump = subprocess.run([driver, "dump"], capture_output=True).stdout.decode()
    blocks = template_blocks()
    for s in base.smiles()[:40] + EXTRA:
        m = Chem.MolFromSmiles(s)
        if m is None or m.GetNumAtoms() == 0 or m.GetNumAtoms() > 60 or "." in s:
            continue
        Chem.Kekulize(m, clearAromaticFlags=True)
        # double bonds in macrocycles need the CIP ranking, which is not ported yet
        ri = m.GetRingInfo()
        if any(b.GetBondType() == Chem.BondType.DOUBLE and ri.MinBondRingSize(b.GetIdx()) >= 9 for b in m.GetBonds()):
            continue
        AllChem.Compute2DCoords(m)
        blocks.append(base.block(m))
    out = subprocess.run([driver], input="".join(blocks).encode(), capture_output=True, timeout=600).stdout.decode()
    recs = out.split("$$\n")[:-1]
    assert len(recs) == len(blocks), (len(recs), len(blocks))
    sys.stdout.write(dump)
    for b, r in zip(blocks, recs):
        sys.stdout.write(b + r + "$$\n")


if __name__ == "__main__":
    main()
