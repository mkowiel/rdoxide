#!/usr/bin/env python3
"""Builds tests/data/avalon_rings.jsonl: molfiles and the bond topography / ring size flags that
PerceiveRingBonds of the Avalon C code sets (driver built from .work/avalon_drv/ring.c; its rand()
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
]


def main():
    drv = sys.argv[1]
    for smi in SMILES:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        AllChem.Compute2DCoords(m)
        for mb in (Chem.MolToMolBlock(m), Chem.MolToV3KMolBlock(m)):
            with tempfile.NamedTemporaryFile("wb", suffix=".mol", dir=os.getcwd(), delete=False) as f:
                f.write(mb.encode())
            try:
                out = subprocess.run([drv, f.name], capture_output=True, timeout=20).stdout.decode()
            finally:
                os.unlink(f.name)
            print(json.dumps({"input": mb, "output": out}))


if __name__ == "__main__":
    main()
