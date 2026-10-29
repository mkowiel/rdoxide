#!/usr/bin/env python3
"""Builds tests/data/freesasa_sasa.txt: coordinates and radii and what freesasa 2.0.3 makes of them
(driver tools/freesasa_drv/sasa.c).

    .venv-ref/bin/python tools/make_freesasa_corpus.py DRIVER > tests/data/freesasa_sasa.txt
"""
import struct
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
RAD = {1: 1.10, 6: 1.70, 7: 1.55, 8: 1.52, 9: 1.47, 15: 1.80, 16: 1.80, 17: 1.75, 35: 1.85, 53: 1.98}
SMILES = [
    "CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "OCC1OC(O)C(O)C(O)C1O", "CC(C)NCC(O)c1ccc(O)c(O)c1", "C1CCCCC1", "c1ccc2ccccc2c1", "CS(=O)(=O)c1ccc(N)cc1",
    "OC(=O)C(N)Cc1ccccc1", "CCCCCCCCCC", "NC(=O)c1ccccc1", "FC(F)(F)c1ccccc1Cl", "C1CC2CCC1C2", "CN1CCC[C@H]1c1cccnc1",
    "O=C1NC(=O)C=C1", "C1COCCOCCOCCOCCO1", "CC(C)(C)c1ccc(O)cc1", "Brc1ccccc1I",
]


def hexd(x):
    return "%016x" % struct.unpack(">Q", struct.pack(">d", x))[0]


def case(alg, pts, probe):
    out = ["CASE %d %d %s" % (alg, len(pts), hexd(probe))]
    for x, y, z, r in pts:
        out.append(" ".join(hexd(v) for v in (x, y, z, r)))
    return "\n".join(out) + "\n"


def main():
    driver = sys.argv[1]
    blocks = []
    for s in SMILES:
        m = Chem.AddHs(Chem.MolFromSmiles(s))
        if AllChem.EmbedMolecule(m, randomSeed=7) != 0:
            continue
        c = m.GetConformer()
        pts = [(c.GetAtomPosition(a.GetIdx()).x, c.GetAtomPosition(a.GetIdx()).y, c.GetAtomPosition(a.GetIdx()).z,
                RAD.get(a.GetAtomicNum(), 1.8)) for a in m.GetAtoms()]
        for alg in (0, 1):
            for probe in (1.4, 0.8):
                blocks.append(case(alg, pts, probe))
    # a protein: the atoms of 1d3z (ubiquitin) from the PDB file with element radii
    pdb = []
    for line in open("upstream/rdkit/External/FreeSASA/test_data/1d3z.pdb"):
        if line.startswith("ATOM") and line[17:20] != "HOH":
            el = line[76:78].strip() or line[12:14].strip()[:1]
            pdb.append((float(line[30:38]), float(line[38:46]), float(line[46:54]),
                        RAD.get(Chem.GetPeriodicTable().GetAtomicNumber(el.capitalize()), 1.8)))
        if len(pdb) >= 700 and line.startswith("ENDMDL"):
            break
    pdb = pdb[:700]
    for alg in (0, 1):
        blocks.append(case(alg, pdb, 1.4))
    inp = "".join(blocks)
    out = subprocess.run([driver], input=inp.encode(), capture_output=True, timeout=1800).stdout.decode()
    recs = out.split("$$\n")[:-1]
    assert len(recs) == len(blocks), (len(recs), len(blocks))
    for b, r in zip(blocks, recs):
        sys.stdout.write(b + r + "$$\n")


if __name__ == "__main__":
    main()
