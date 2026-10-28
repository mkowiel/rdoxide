#!/usr/bin/env python3
"""Builds tests/data/avalon_canon.jsonl: molfiles and the CanonicalNumbering result (numbers,
ordered neighbourhoods, molecule after the stereo clean-up) of the Avalon C code
(driver: .work/avalon_drv/canon.c).

    .venv-ref/bin/python tools/make_avalon_canon_corpus.py DRIVER > tests/data/avalon_canon.jsonl
"""
import json
import os
import random
import subprocess
import sys
import tempfile

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

SMILES = [
    "CCO", "c1ccccc1", "C1CCCCC1", "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", "OC[C@@H](O)[C@H](O)CO", "C/C=C/C", "C/C=C\\C",
    "F/C=C/C=C/Cl", "F/C=C\\C=C/Cl", "C[C@H]1CC[C@@H](C)CC1", "C[C@H]1CC[C@H](C)CC1", "[C@H](F)(Cl)Br", "C12C3C4C1C5C2C3C45",
    "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O", "CC(C)C[C@H](N)C(O)=O", "C=C=C", "CC=C=CC", "C/C=C=C/C",
    "OC(=O)C(O)C(O)C(=O)O", "O[C@H](C(=O)O)[C@@H](O)C(=O)O", "O[C@H](C(=O)O)[C@H](O)C(=O)O", "C1CC1C1CC1", "c1ccc2ccccc2c1",
    "N[C@@H](Cc1ccccc1)C(=O)O", "CC(C)(C)C", "FC(F)(F)C(F)(F)F", "[13CH3]C", "[2H]C([2H])Cl", "[3H]C", "C[N+](C)(C)C", "[O-]C(=O)C",
    "CN1CCC[C@H]1c1cccnc1", "C[C@@H]1CCCC[C@H]1C", "C[C@H]1CCCC[C@@H]1C", "ClC=CC=CCl", "C/C=C/C=C/C", "CC1=CC=CC=C1",
    "C[C@H](F)[C@@H](F)C", "C[C@H](Cl)C[C@@H](C)Cl", "C[C@H]1C[C@@H](C)C[C@H](C)C1", "OC1CC(O)CC(O)C1", "[Na+].[Cl-]",
    "C1CC2CCC1C2", "C1C2CC3CC1CC(C2)C3", "[H]C([H])([H])[H]", "O=C1CCCC1", "C1=CC=CC=C1C=CC1=CC=CC=C1",
]


def blocks():
    rng = random.Random(3)
    for smi in SMILES:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        AllChem.Compute2DCoords(m)
        Chem.WedgeMolBonds(m, m.GetConformer())
        yield Chem.MolToMolBlock(m)
        yield Chem.MolToMolBlock(m, kekulize=False)
        mh = Chem.AddHs(m)
        AllChem.Compute2DCoords(mh)
        Chem.WedgeMolBonds(mh, mh.GetConformer())
        yield Chem.MolToMolBlock(mh)
        order = list(range(m.GetNumAtoms()))
        rng.shuffle(order)
        yield Chem.MolToMolBlock(Chem.RenumberAtoms(m, order))


def main():
    drv = sys.argv[1]
    for mb in blocks():
        with tempfile.NamedTemporaryFile("wb", suffix=".mol", dir=os.getcwd(), delete=False) as f:
            f.write(mb.encode())
        try:
            out = subprocess.run([drv, f.name], capture_output=True, timeout=20).stdout.decode()
        finally:
            os.unlink(f.name)
        print(json.dumps({"input": mb, "output": out}))


if __name__ == "__main__":
    main()
