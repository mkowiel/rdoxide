#!/usr/bin/env python3
"""Builds tests/data/avalon_layout2.jsonl: molfiles and what the second part of Avalon's layout.c
(GetColoredGraph, ColorStrain, ImproveBondByFlip, ImproveBondByStretch, ImproveMoleculeByBondFlip,
ImproveSPAtoms, CountColor) produces after the ring layout of part 1
(driver: .work/avalon_drv/lay2.c, built with -ffp-contract=off).

    .venv-ref/bin/python tools/make_avalon_layout2_corpus.py DRIVER > tests/data/avalon_layout2.jsonl
"""
import json
import os
import subprocess
import sys
import tempfile

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

SMILES = [
    "CC(C)(C)CC(C)(C)C", "CCCCCC(CC)CC", "CC#CC", "C#CC#C", "CC(C#N)C", "C=C=C", "CC=C=CC", "CC(=O)N(C)CC(C)C", "OCC(O)C(O)C(O)CO", "CN(C)CCOC(c1ccccc1)c1ccccc1", "CC(C)Cc1ccc(C(C)C(=O)O)cc1",

    "CCO", "c1ccccc1", "C1CC1", "C1CCCCC1", "c1ccc2ccccc2c1", "C12C3C4C1C5C2C3C45", "C1CC2CCC1C2", "C1CCC2(CC1)CCCC2", "CC12CCC3C(CCc4cc(O)ccc34)C1CCC2O",
    "c1ccc(cc1)-c1ccccc1", "C1CC1C1CC1", "c1ccc2c(c1)[nH]c1ccccc12", "C1CCC2CCCCC2C1", "O=C1CCCCC1", "C1=CC=CC=C1C=CC1=CC=CC=C1", "C/C=C/C=C/C",
    "F/C=C/C=C\\Cl", "CC(=O)Oc1ccccc1C(=O)O", "C1CCCCCCCC1", "c1cc2cccc3ccc4cccc1c4c23", "N[C@@H](Cc1ccccc1)C(=O)O", "C1CC2CC1CC2", "[Fe]", "C[Fe]C",
    "C1CC1.C1CCC1", "CC.CC", "c1ccccc1[Pt](Cl)Cl", "[Cu]1(N)(N)OCCO1", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "C1=CC2=CC=CC=C2C=C1",
    "C1CCC1C1CCC1", "[Na+].[Cl-]", "C12CC1C2", "C1CC11CC1", "c1ccc2[nH]ccc2c1", "CN1CCC[C@H]1c1cccnc1",
]


def blocks():
    for smi in SMILES:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        AllChem.Compute2DCoords(m)
        Chem.WedgeMolBonds(m, m.GetConformer())
        yield Chem.MolToMolBlock(m)
        mh = Chem.AddHs(m)
        AllChem.Compute2DCoords(mh)
        yield Chem.MolToMolBlock(mh)


def main():
    drv = sys.argv[1]
    for mb in blocks():
        with tempfile.NamedTemporaryFile("wb", suffix=".mol", dir=os.getcwd(), delete=False) as f:
            f.write(mb.encode())
        try:
            out = subprocess.run([drv, f.name], capture_output=True, timeout=60).stdout.decode()
        finally:
            os.unlink(f.name)
        print(json.dumps({"input": mb, "output": out}))


if __name__ == "__main__":
    main()
