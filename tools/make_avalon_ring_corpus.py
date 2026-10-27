#!/usr/bin/env python3
"""Builds tests/data/avalon_rings.jsonl: molfiles and what PerceiveRingBonds, PerceiveAromaticBonds,
PerceiveDYAromaticity and SetRingSizeFlags of the Avalon C code leave in the molecule (driver built from .work/avalon_drv/ring.c; its rand()
is the glibc generator).

    .venv-ref/bin/python tools/make_avalon_ring_corpus.py DRIVER > tests/data/avalon_rings.jsonl
"""
import json
import os
import subprocess
import sys
import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem

SMILES = [
    "c1ccccc1", "C1CC1", "C1CCC1", "c1ccc2ccccc2c1", "C12C3C4C1C5C2C3C45", "C1CC2CCC1C2",
    "C1CC2(C1)CCC2", "C1CCC2(CC1)CCCCC2", "CC12CCC3C(CCc4cc(O)ccc34)C1CCC2O", "c1ccc2cc3ccccc3cc2c1",
    "c1cc2ccc3cccc4ccc(c1)c2c34", "C1C2CC3CC1CC(C2)C3", "C1CC2CC1C1CCC21", "O=C1CCCCC1", "C1CCC1C1CCC1",
    "C1CC1C1CC1", "c1ccc(cc1)-c1ccccc1", "C1C2CC3C1C3C2", "C12CC1C2", "c1ccc2c(c1)[nH]c1ccccc12",
    "C1CC2C3CCC(C3)C2C1", "CC(C)C", "CCO", "C1=CC2=CC=CC=CC2=C1", "C1CC2CC3CC1CC(C2)C3",
    "C12C3C4C5C1C6C2C3C4C56", "O1CCOCC1", "C1CCC2C(C1)CCC1CCCCC12",
    "c1cc2cc3cc4cc5ccc6cc7ccc8c9c%10c(c%11c(c(c7c6c5c4c3c2c1)c89)cccc%11)cc%10", "C1CC2CCC3CCC1C2C3",
    "C1CC2CC2C1", "C1C2CC1C2", "C1CC11CC1", "C12CC3CC(C1)CC(C3)C2", "C1CCC2CC3CCCCC3CC2C1",
    "C1=CC=C2C(=C1)C3=CC=CC=C3C4=CC=CC=C24",
    "c1cc[nH]c1", "c1ccoc1", "c1ccsc1", "c1ccncc1", "c1ccc2[nH]ccc2c1", "c1nc2[nH]cnc2c(=O)[nH]1",
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "c1cnc[nH]1", "c1ccc2ncccc2c1", "O=c1cc[nH]cc1", "O=C1C=CC(=O)C=C1",
    "[cH-]1cccc1", "c1ccc2cccc2cc1", "C1=CC=CC=CC=C1", "[CH+]1C=CC=CC=C1", "c1cc[o+]cc1", "c1ccc[n+]([O-])c1",
    "c1cnoc1", "c1cscn1", "O=c1ccoc2ccccc12", "c1ccc2occc2c1", "n1ccnc2ccccc12", "C1=CC=CC=C1C1=CC=CC=C1",
    "C1=CNC=C1", "C1=CCC=C1", "C1=CC=CN1", "O=C1NC=CC=C1", "c1ccc2c(c1)ccc1ccccc12", "c1c[nH]c2ncccc12",
    "[nH]1cccc1-c1ccccc1", "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61", "S1C=CC=C1", "N1C=CC=C1C=O", "C=C1C=CC=CC1=O",
]
MODES = ["ring", "arom", "dy", "rsize"]


def main():
    drv = sys.argv[1]
    for smi in SMILES:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        AllChem.Compute2DCoords(m)
        mk = Chem.Mol(m)
        Chem.Kekulize(mk, clearAromaticFlags=True)
        for mb in (Chem.MolToMolBlock(m), Chem.MolToMolBlock(mk, kekulize=False), Chem.MolToV3KMolBlock(m)):
            for mode in MODES:
                with tempfile.NamedTemporaryFile("wb", suffix=".mol", dir=os.getcwd(), delete=False) as f:
                    f.write(mb.encode())
                try:
                    out = subprocess.run([drv, f.name, mode], capture_output=True, timeout=20).stdout.decode()
                finally:
                    os.unlink(f.name)
                print(json.dumps({"input": mb, "mode": mode, "output": out}))


if __name__ == "__main__":
    main()
