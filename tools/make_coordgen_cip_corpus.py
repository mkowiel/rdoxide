#!/usr/bin/env python3
"""Builds tests/data/coordgen_cip.txt: stereo molecules (RDKit 2D coordinates and wedges) and what the CIP / stereo code
of coordgenlibs 3.0.2 answers (driver tools/coordgen_drv/cip.cpp).

    .venv-ref/bin/python tools/make_coordgen_cip_corpus.py DRIVER > tests/data/coordgen_cip.txt
"""
import importlib.util
import os
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("m", os.path.join(HERE, "make_coordgen_model_corpus.py"))
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

SMILES = [
    "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", "N[C@@H](Cc1ccccc1)C(=O)O", "OC[C@H](O)[C@@H](O)[C@H](O)CO", "F[C@](Cl)(Br)I",
    "C[C@H](F)[C@@H](Cl)Br", "C[S@](=O)CC", "C[N@+](CC)(CCC)CCCC", "CC(C)[C@H](N)C(C)C", "C[C@H]1CC[C@@H](C)CC1",
    "C[C@H]1CC[C@H](C)CC1", "C[C@@H]1CCC[C@H](C)C1", "C[C@H]1C[C@@H](C)C1", "C1C[C@H]2CC[C@@H]1C2", "C[C@H](O)[C@H](C)O",
    "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O",
    "C/C=C/C", "C/C=C\\C", "F/C=C/Cl", "F/C=C\\Cl", "C/C(F)=C/C", "CC/C=C/C(C)C", "C/C=C/C=C/C", "OC(=O)/C=C/c1ccccc1",
    "C1CCCCCC/C=C/CCCCC1", "C1CCCCCC/C=C\\CCCCC1", "C1CC/C=C/CCCCCCCC1", "C1CCCCCCCCCCCC=CC1", "CC(C)(C)C", "N[C@H](C)C(=O)N[C@@H](C)C(=O)O",
    "C[C@H](Cl)[C@H](C)Br", "CC[C@H](C)[C@@H](C)CC", "[2H][C@](F)(Cl)Br", "C[C@@H](c1ccccc1)N", "N[C@@H]1CCCC[C@H]1N", "CC[C@@]1(C)CCCC[C@H]1C",
    "C[C@H](CC)[C@H](CC)C", "OC[C@@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "CP(=O)(O)[C@H](C)N", "C[P@](=O)(CC)c1ccccc1",
]


def main():
    driver = sys.argv[1]
    blocks = []
    for s in SMILES + base.smiles()[::4]:
        m = Chem.MolFromSmiles(s)
        if m is None or m.GetNumAtoms() == 0 or m.GetNumAtoms() > 60:
            continue
        for frag in Chem.GetMolFrags(m, asMols=True):
            Chem.Kekulize(frag, clearAromaticFlags=True)
            AllChem.Compute2DCoords(frag)
            Chem.WedgeMolBonds(frag, frag.GetConformer())
            lines = base.block(frag).split("\n")
            for b in frag.GetBonds():
                d = b.GetBondDir()
                w = 1 if d == Chem.BondDir.BEGINWEDGE else 2 if d == Chem.BondDir.BEGINDASH else 0
                lines[1 + frag.GetNumAtoms() + b.GetIdx()] += " %d" % w
            blocks.append("\n".join(lines))
    inp = "".join(b if b.endswith("\n") else b + "\n" for b in blocks)
    out = subprocess.run([driver], input=inp.encode(), capture_output=True, timeout=600).stdout.decode()
    recs = out.split("$$\n")[:-1]
    assert len(recs) == len(blocks), (len(recs), len(blocks))
    for b, r in zip(blocks, recs):
        sys.stdout.write((b if b.endswith("\n") else b + "\n") + r + "$$\n")


if __name__ == "__main__":
    main()
