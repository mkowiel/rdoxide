#!/usr/bin/env python3
"""Builds tests/data/coordgen_sketch.txt: molecules (also salts and mixtures, with double bond stereo) and the coordinates
that coordgenlibs 3.0.2 makes of them the way RDKit's CoordGen::addCoords calls it (driver tools/coordgen_drv/sketch.cpp).

    .venv-ref/bin/python tools/make_coordgen_sketch_corpus.py DRIVER > tests/data/coordgen_sketch.txt
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

EXTRA = [
    "CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "C/C=C/C", "C/C=C\\C", "F/C=C/Cl", "OC(=O)/C=C/c1ccccc1", "CC(C)(C)c1ccccc1",
    "[Na+].[Cl-]", "CC(=O)[O-].[Na+]", "c1ccccc1.CCO", "CCN(CC)CC", "C1CCCCC1C(=O)O", "CS(=O)(=O)c1ccc(N)cc1", "OP(=O)(O)OCC",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "c1ccc2ccccc2c1", "c1ccc2c(c1)ccc1ccccc12", "C1CC2CCC1C2", "C12C3C4C1C5C2C3C45",
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "CC(C)NCC(O)c1ccc(O)c(O)c1", "N#Cc1ccccc1", "O=C1NC(=O)C=C1", "C1=CC=CC=C1C=CC1=CC=CC=C1",
    "[NH4+].[Cl-]", "OCC(O)C(O)C(O)C(O)CO", "CC=CC=CC", "C/C=C/C=C/C", "c1ccc(cc1)-c1ccccc1", "CCCCCCCCCC", "C1CCC2(CC1)CCCC2",
    "NC(=O)c1ccccc1.O", "FC(F)(F)C(F)(F)F", "CN1CCC[C@H]1c1cccnc1", "C[C@H](N)C(=O)O", "O=S(=O)(N)c1ccccc1",
]


def ok(m):
    ri = m.GetRingInfo()
    if any(len(r) >= 9 for r in ri.AtomRings()):
        return False
    if any(ri.NumBondRings(b.GetIdx()) > 2 for b in m.GetBonds()):
        return False
    return True


def main():
    driver = sys.argv[1]
    blocks = []
    for s in EXTRA + base.smiles()[::2]:
        m = Chem.MolFromSmiles(s)
        if m is None or m.GetNumAtoms() == 0 or m.GetNumAtoms() > 45 or not ok(m):
            continue
        stereo = []
        for b in m.GetBonds():
            if (
                b.GetBondType() == Chem.BondType.DOUBLE
                and b.GetStereo() in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOTRANS)
                and len(list(b.GetStereoAtoms())) >= 2
            ):
                a1, a2 = list(b.GetStereoAtoms())[:2]
                cis = b.GetStereo() in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS)
                stereo.append("S %d %d %d %s" % (b.GetIdx(), a1, a2, "c" if cis else "t"))
        Chem.Kekulize(m, clearAromaticFlags=True)
        AllChem.Compute2DCoords(m)
        lines = base.block(m).rstrip("\n").split("\n")
        assert lines[-1] == "END"
        lines = lines[:-1] + stereo + ["END"]
        blocks.append("\n".join(lines) + "\n")
    inp = "".join(blocks)
    out = subprocess.run([driver], input=inp.encode(), capture_output=True, timeout=1800).stdout.decode()
    recs = out.split("$$\n")[:-1]
    assert len(recs) == len(blocks), (len(recs), len(blocks))
    for b, r in zip(blocks, recs):
        sys.stdout.write(b + r + "$$\n")


if __name__ == "__main__":
    main()
