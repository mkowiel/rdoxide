#!/usr/bin/env python3
"""Builds tests/data/coordgen_macro.txt: molecules with rings of nine or more atoms (macrocycles, also fused and substituted, with double bond stereo) and the coordinates
that coordgenlibs 3.0.2 makes of them the way RDKit's CoordGen::addCoords calls it (driver tools/coordgen_drv/sketch.cpp).

    .venv-ref/bin/python tools/make_coordgen_macro_corpus.py DRIVER > tests/data/coordgen_macro.txt
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
    "C1CCCCCCCCC1", "C1CCCCCCCCCC1", "C1CCCCCCCCCCCC1", "C1CCCCCCCCCCCCCC1", "C1CCCCCCCCCCCCCCCCCC1",
    "C1COCCOCCOCCOCCOCCO1", "C1COCCOCCOCCOCCO1", "C1CNCCNCCNCCNCCN1", "C1CSCCSCCSCCSC1", "C1COCCOCCOCCOC1",
    "C1CCCCC=CCCCC1", "C1CCCCC/C=C/CCCC1", "C1CCCCC/C=C\\CCCC1", "C1CC/C=C/CC/C=C/CC1", "C1CC/C=C\\CC/C=C/CCCC1",
    "O=C1CCCCCCCCCCCC1", "O=C1CCCCCCCCCCCCCCO1", "CC1CCCCCCCCCCC(=O)O1", "OC1CCCCCCCCCCC1", "NC1CCCCCCCCCCC1",
    "C1CCCCCCCCCC1C", "CC1CCCCCCCCCC(C)C1", "c1ccc2c(c1)CCCCCCCCC2", "C1CCCCCCCCC2CCCCC2C1", "C1CCCCCCCCCCC2CCC12",
    "c1ccc2c(c1)OCCOCCOCCO2", "c1ccc2c(c1)OCCOCCOCCOCCO2", "c1cc2OCCOCCOCCOc2cc1", "C1CCCCCCCCC2(CC2)C1",
    "O=C1NCCCCCCCCCN1", "O=C1CCCCCCCCCCN1", "O=C1CCCCCCCCCCCCCCCCN1", "O=C1NC(C)C(=O)NC(C)C(=O)NC(C)C(=O)N1",
    "O=C1CNC(=O)CNC(=O)CNC(=O)CNC(=O)CN1", "CC1NC(=O)C(C)NC(=O)C(C)NC(=O)C(C)NC(=O)C(C)NC1=O",
    "O=C1NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)N1", "CC(C)C1NC(=O)C(Cc2ccccc2)NC(=O)C(C)NC(=O)C(CC(C)C)NC(=O)CNC1=O",
    "C1CCC(CC1)C1CCCCCCCCC1", "C1=CC=CC=CC=CC=CC=CC=C1", "c1ccc2cc3ccc(cc3cc2c1)-c1ccccc1",
    "CC1CC(C)C(O)C(C)C(=O)C(C)C(O)C(C)C(=O)OC(CC)C(C)C1=O", "CCC1OC(=O)C(C)C(O)C(C)C(O)C(C)C(=O)C(C)CC(C)(O)C1C",
    "O=C1CC(O)CC(O)CC(O)CC(O)CC(O)CCC(O)CCC(O)CC(O)CC=CC=CC=CC=CC=CC=CC(O)C(C)OC(=O)C1",
    "C1CCCCCCCCCCCC1.CCO", "C1CCCCCCCCCCCC1.[Na+].[Cl-]", "C1COCCOCCOCCOCCOCCO1.[K+]",
    "C1CCCCCCCCCNCCCCCCCCC1", "C1CCCCCCCCC1CCCCCCCCC", "C1CCCCCCCCCC1C1CCCCCCCCC1",
    "C1CCCCCCCCCCCC1(C)C", "CC1(C)CCCCCCCCCC1", "FC1CCCCCCCCC1F", "ClC1CCCCCCCCCCC1Cl",
    "C1CC2CCC1CCCCCCCC2", "C1CCCCCCCCC2CCCC12", "C1CCCCCCCCC(C1)C(=O)O", "C1CCCCCCCCCC=1", "C1=CCCCCCCCC1",
    "C1CCCCCCCCC=CC1", "C1CCCCCCCCCCCCC1", "C1CCCCCCCCCCCCCCC1", "C1CCCCCCCCCCCCCCCC1", "C1CCCCCCCCCCCCCCCCC1",
    "C1CCCC/C=C/CCC1", "C1CC/C=C/C=C/CCC1", "C1/C=C/C=C/C=C/CCCC1", "C1C/C=C\\C/C=C\\C/C=C\\C1", "C1CC/C=C/CC/C=C\\CC1",
    "O=C1/C=C/CCCCCCCCC1", "O=C1/C=C\\CCCCCCCCC1", "C1CCCC/C=C/C(=O)OCC1", "C1CCCCC(=O)OCCCCCC(=O)OC1",
    "c1ccc2c(c1)CC/C=C/CCCC2", "c1ccc2c(c1)CCCCCCCCCC2", "c1cc2cc1CCCCCCCC2", "c1ccc2c(c1)C(=O)NCCCCCCCCN2",
    "C1CCCCCCC2CCC3CCCC3C2CCC1", "C1CCCCCCCC2CCCCC2CC1", "C1CCCCCCCC2CC2CC1", "C1CCCCCCCC2OC2CC1",
    "C1CCCCCCCC2C1CCC2", "C1CCCCCCCC2C1C=CC2", "C1CCCCCCCCC(C1)C1CCCCCCCCC1", "C1CCC2(CC1)CCCCCCCCC2",
    "CCCC1CCCCCCCCCC1CCC", "CC1CCCC(C)CCCCCCC1C", "C1CCCCCCCCC(CCCCC)CC1", "OC1CCCCCCCCC(O)CCC1O", "NC1CCCCCCCC(N)CCC1N",
    "C1COCCOCCOCCOCCOC1", "C1COCCOCCNCCOCCOC1", "C1CCCCCCCCCCCCCCCCCCCCC1", "C1CCCCCCCCCCCCCCCCCCCCCCCCC1",
    "CC(C)C1CCCCCCCCC(C(C)C)C1", "C[C@H]1CCCCCCCC[C@@H](C)C1", "C[C@H]1CCCCCCCCCC(=O)O1", "O[C@H]1CCCCCCCCCC[C@@H]1O",
    "CN1C(=O)CN(C)C(=O)CN(C)C(=O)CN(C)C(=O)CN(C)C(=O)C1", "CC1=CC=CC=CC=CC=CC=CC=C1",
    "O=C1OC(C)CC=CC=CC(O)CC(O)CC(O)CCC(O)C(C)C(O)C1", "CC1CC=CC=CC(O)CC(=O)OC(C)CC=CC1",
    "C1=C/CC/C=C/CC/C=C/CC1", "C1CC=CCC=CCC=CC1", "C1=CC=CC=CC=CC=C1", "C1C#CCCCCCC1", "C1CCCCCC#CCCCC1",
    "CC1CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)C1", "CC1CC(CC)CC(CC)CCC(C)CC(CC)CC(C)C1", "C1CC(C(C)(C)C)CCC(C(C)(C)C)CCC(C(C)(C)C)CC1",
    "OC1CC(O)CC(O)CC(O)CC(O)CC(O)CC(O)C1", "NC1CC(N)CC(N)CC(N)CC(N)CC(N)C1", "CC1OC(C)CC(C)OC(C)CC(C)OC(C)C1",
    "C1CC(CCO)CCC(CCO)CCC(CCO)CC1", "O=C1CC(C)CC(=O)CC(C)CC(=O)CC(C)C1", "CC(C)C1CC(C(C)C)CC(C(C)C)CC(C(C)C)CC(C(C)C)C1",
    "CC1CCC(C)CCC(C)CCC(C)CCC(C)CCC1C", "CC1(C)CC(C)(C)CC(C)(C)CC(C)(C)CC(C)(C)C1", "COC1CC(OC)CC(OC)CC(OC)CC(OC)C1",
    "CSC1CC(SC)CC(SC)CC(SC)CC(SC)C1", "FC1CC(F)CC(F)CC(F)CC(F)C1", "CC1CNC(C)CNC(C)CNC(C)CNC(C)CN1",
    "CC1NC(=O)C(CC(C)C)NC(=O)C(Cc2ccccc2)NC(=O)C(CO)NC(=O)C(CC(=O)O)NC1=O", "CC1NC(=O)C(CCCN)NC(=O)C(CC(N)=O)NC(=O)C(C)NC(=O)C(Cc2c[nH]c3ccccc23)NC1=O",
    "CC1OC(=O)CC(O)C(C)C(O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(=O)C(C)CC1(C)O",
    "CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(N(C)C)C2O)C(C)(O)CC(C)C(=O)C(C)C(O)C1(C)O",
]


def ok(m):
    ri = m.GetRingInfo()
    if not any(len(r) >= 9 for r in ri.AtomRings()):
        return False
    if any(ri.NumBondRings(b.GetIdx()) > 2 for b in m.GetBonds()):
        return False
    return True


def main():
    driver = sys.argv[1]
    blocks = []
    for s in EXTRA + [x for x in base.smiles() if "1" in x]:
        m = Chem.MolFromSmiles(s)
        if m is None or m.GetNumAtoms() == 0 or m.GetNumAtoms() > 60 or not ok(m):
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
