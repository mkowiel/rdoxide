#!/usr/bin/env python3
"""Builds tests/data/pubchem_shape_align.txt: pairs of 3D molecules (atoms, radii, colour types) and the volumes, steric
frame classes and optimized overlay pubchem-align3d computes for them (driver tools/pcshape_drv/align.cpp).

    .venv-ref/bin/python tools/make_pubchem_shape_corpus.py DRIVER > tests/data/pubchem_shape_align.txt
"""
import struct
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
RAD = {1: 1.10, 6: 1.70, 7: 1.55, 8: 1.52, 9: 1.47, 15: 1.80, 16: 1.80, 17: 1.75, 35: 1.85}
SMILES = [
    "CC(=O)Oc1ccccc1C(=O)O", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "CC(C)NCC(O)c1ccc(O)c(O)c1",
    "OC(=O)C(N)Cc1ccccc1", "c1ccc2ccccc2c1", "CCCCCCCC", "OCC1OC(O)C(O)C(O)C1O", "CN1CCC[C@H]1c1cccnc1", "C1CCCCC1",
    "c1ccccc1", "CC(=O)Nc1ccc(O)cc1", "CS(=O)(=O)c1ccc(N)cc1", "O=C1NC(=O)C=C1", "ClC(Cl)(Cl)Cl",
]


def hexd(x):
    return "%016x" % struct.unpack(">Q", struct.pack(">d", x))[0]


def hexf(x):
    return "%08x" % struct.unpack(">I", struct.pack(">f", x))[0]


def atoms(m, hs):
    if not hs:
        m = Chem.RemoveHs(m)
    c = m.GetConformer()
    out = []
    for a in m.GetAtoms():
        p = c.GetAtomPosition(a.GetIdx())
        z = a.GetAtomicNum()
        color = 0
        if z == 7:
            color = 1
        elif z == 8:
            color = 2 if any(n.GetAtomicNum() == 1 for n in a.GetNeighbors()) else 3
        elif z == 6 and a.GetIsAromatic():
            color = 4
        out.append("%s %s %s %s %d" % (hexf(p.x), hexf(p.y), hexf(p.z), hexd(RAD.get(z, 1.8)), color))
    return out


def main():
    driver = sys.argv[1]
    mols = []
    for s in SMILES:
        m = Chem.AddHs(Chem.MolFromSmiles(s))
        if AllChem.EmbedMolecule(m, randomSeed=11) == 0:
            AllChem.MMFFOptimizeMolecule(m)
            mols.append(m)
    blocks = []
    for i, a in enumerate(mols):
        for j in (i, (i + 1) % len(mols), (i + 4) % len(mols)):
            b = mols[j]
            for hs, color, opt, pre, post in ((False, 1, 0.5, 3, 16), (False, 0, 1.0, 10, 30), (True, 1, 1.0, 10, 30)):
                ra, fb = atoms(a, hs), atoms(b, hs)
                if not any(x.endswith(' 0') for x in ra) or not any(x.endswith(' 0') for x in fb):
                    continue  # no volume atoms: upstream returns an uninitialized pose
                blocks.append("CASE %d %d %d %s %d %d\n" % (len(ra), len(fb), color, hexd(opt), pre, post) + "\n".join(ra + fb) + "\n")
    inp = "".join(blocks)
    out = subprocess.run([driver], input=inp.encode(), capture_output=True, timeout=3000).stdout.decode()
    recs = out.split("$$\n")[:-1]
    assert len(recs) == len(blocks), (len(recs), len(blocks))
    for b, r in zip(blocks, recs):
        sys.stdout.write(b + r + "$$\n")


if __name__ == "__main__":
    main()
