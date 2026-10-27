#!/usr/bin/env python3
"""Builds tests/data/avalon_smiles.jsonl: molfiles, flags and the SMILES the Avalon C code's MOLToSMI
writes for them (driver: .work/avalon_drv/smi.c).

    .venv-ref/bin/python tools/make_avalon_smiles_corpus.py DRIVER > tests/data/avalon_smiles.jsonl
"""
import json
import os
import subprocess
import sys
import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem

SMILES = [
    "CCO", "c1ccccc1O", "C[N+](C)(C)C", "[O-]C(=O)C", "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", "C/C=C/C", "C/C=C\\C",
    "[13CH3]C", "[2H]C([2H])Cl", "[CH3]", "[Na+].[Cl-]", "C1CC1", "O=S(=O)(N)c1ccccc1", "[Fe+3]",
    "CC(C)(C)c1ccc(O)cc1", "C#CC", "[NH4+]", "C[Si](C)(C)C", "FC(F)(F)C(=O)O", "[14C]", "[18F]C",
    "OC[C@@H](O)[C@H](O)CO", "[CH2]C", "c1ccc2ccccc2c1", "[Se]", "[3H]O", "[N+](=O)[O-]",
    "C12C3C4C1C5C2C3C45", "CC12CCC3C(CCc4cc(O)ccc34)C1CCC2O", "c1cc[nH]c1", "c1ccoc1", "c1ccsc1", "c1ccncc1",
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "O=C1C=CC(=O)C=C1", "F/C=C/C=C/Cl", "F/C=C\\C=C/Cl", "C(=O)(O)[C@@H](N)Cc1ccccc1",
    "N[C@@H](C)C(=O)N[C@H](C(=O)O)Cc1c[nH]c2ccccc12", "C[C@H]1CC[C@@H](C)CC1", "CC(=O)Oc1ccccc1C(=O)O",
    "C=C=C", "CC=C=CC", "OC(=O)C(O)C(O)C(=O)O", "C1CCC1C1CCC1", "c1ccc(cc1)-c1ccccc1", "[Cl-].[Cl-].[Mg+2]",
    "CC[S+](C)C", "C[P](C)(C)=O", "CS(=O)C", "O=C1CCCN1", "CC(C)C[C@H](N)C(O)=O", "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O",
    "[C@H](F)(Cl)Br", "[C@@H](F)(Cl)Br", "C[C@](F)(Cl)Br", "N#Cc1ccccc1", "C1=CC=CC=C1C=O", "CC(C)=O", "[H][H]",
    "O=C(c1ccccc1)c1ccccc1", "CCN(CC)CC", "OCC(O)CO", "NC(=N)N", "C1CC2CCC1C2", "CN1CCC[C@H]1c1cccnc1",
    "[O-][n+]1ccccc1", "c1ccc2[nH]ccc2c1", "C[C@@H]1CCCC[C@H]1C", "C(/C=C/c1ccccc1)=O", "[nH]1cccc1",
]
FLAGS = [0, 4, 0x2000, 0x2004, 0x40, 0x44]


def blocks():
    for smi in SMILES:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        AllChem.Compute2DCoords(m)
        yield Chem.MolToMolBlock(m)
        yield Chem.MolToMolBlock(m, kekulize=False)
        mh = Chem.AddHs(m)
        AllChem.Compute2DCoords(mh)
        yield Chem.MolToMolBlock(mh)
        m3 = Chem.AddHs(m)
        if AllChem.EmbedMolecule(m3, randomSeed=11) == 0:
            yield Chem.MolToMolBlock(m3)


def main():
    drv = sys.argv[1]
    for mb in blocks():
        for fl in FLAGS:
            with tempfile.NamedTemporaryFile("wb", suffix=".mol", dir=os.getcwd(), delete=False) as f:
                f.write(mb.encode())
            try:
                out = subprocess.run([drv, f.name, str(fl)], capture_output=True, timeout=20).stdout.decode()
            finally:
                os.unlink(f.name)
            print(json.dumps({"input": mb, "flags": fl, "output": out}))


if __name__ == "__main__":
    main()
