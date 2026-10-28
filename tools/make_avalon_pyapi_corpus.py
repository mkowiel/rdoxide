#!/usr/bin/env python3
"""Builds tests/data/avalon_pyapi.jsonl: what the official rdkit.Avalon.pyAvalonTools gives for GetCanonSmiles,
Generate2DCoords and the data string flavours of the fingerprint functions.

    .venv-ref/bin/python tools/make_avalon_pyapi_corpus.py COORDS_DRIVER > tests/data/avalon_pyapi.jsonl

The layout draws on the C library's `rand()`, which differs between platforms; the official wheel's coordinates
therefore depend on the platform it was built for. The coordinates of `Generate2DCoords(data, ...)` are taken from
the C oracle `.work/avalon_drv/coords.c` (glibc-style `rand`) instead, and for `Generate2DCoords(mol)` only
success or failure (which does not depend on `rand`) and the conformer id are recorded.
"""
import json
import os
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Avalon import pyAvalonTools as p
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
sys.path.insert(0, os.path.dirname(__file__))
from make_avalon_smiles_corpus import SMILES  # noqa: E402

EXTRA = [
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "c1ccc2c(c1)[nH]c1ccccc12",
    "C[C@@H]1CC[C@H]2[C@@H](C1)OC(=O)C2", "O=C(O)c1ccccc1O", "F/C=C/C=C\\Cl", "C[C@H](N)C(=O)O.Cl", "[Na+].CC(=O)[O-]", "C1CCCCC1C1CCCCC1", "[2H]C([2H])O",
    "CCO.CCO", "N[C@@H](Cc1ccccc1)C(=O)O", "c1ccccc1-c1ccccc1", "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O", "[H]C(=O)O",
]
FLAGS = [-1, 0, 1, 3, 0x1B, 0xFB, 0x2003, 0x1003]


def main():
    smis = SMILES[::2] + EXTRA
    coords = []
    for i, smi in enumerate(smis):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        mb = Chem.MolToMolBlock(mol)
        m2 = Chem.Mol(mol)
        AllChem.Compute2DCoords(m2)
        mb2 = Chem.MolToMolBlock(m2)
        for fl in (FLAGS[i % len(FLAGS)], -1):
            print(json.dumps({"kind": "canon_mol", "input": smi, "flags": fl, "out": p.GetCanonSmiles(mol, fl)}))
            print(json.dumps({"kind": "canon_mol2d", "input": mb2, "flags": fl, "out": p.GetCanonSmiles(m2, fl)}))
            print(json.dumps({"kind": "canon_smiles", "input": smi, "flags": fl, "out": p.GetCanonSmiles(smi, True, fl)}))
            print(json.dumps({"kind": "canon_molblock", "input": mb, "flags": fl, "out": p.GetCanonSmiles(mb, False, fl)}))
        coords.append(("coords_smiles", smi, "SMILES:" + smi))
        coords.append(("coords_molblock", mb, mb))
        m3 = Chem.Mol(mol)
        try:
            cid = p.Generate2DCoords(m3)
            print(json.dumps({"kind": "coords_mol", "input": smi, "cid": cid, "natoms": m3.GetNumAtoms()}))
        except RuntimeError:
            print(json.dumps({"kind": "coords_mol", "input": smi, "error": True}))
        mh = Chem.AddHs(mol)
        try:
            cid = p.Generate2DCoords(mh)
            print(json.dumps({"kind": "coords_molh", "input": smi, "cid": cid, "natoms": mh.GetNumAtoms()}))
        except RuntimeError:
            print(json.dumps({"kind": "coords_molh", "input": smi, "error": True}))
        for nb, q, fl in [(512, False, 15761407), (1024, True, 0x7FFF), (256, False, 0x10)]:
            fp = p.GetAvalonFP(smi, True, nb, q, False, fl)
            print(json.dumps({"kind": "fp_smiles", "input": smi, "args": [nb, q, fl], "out": list(fp.GetOnBits())}))
            print(json.dumps({"kind": "words_smiles", "input": smi, "args": [nb, q, fl], "out": list(p.GetAvalonFPAsWords(smi, True, nb, q, False, fl))}))
            cf = p.GetAvalonCountFP(smi, True, nb, q, fl)
            print(json.dumps({"kind": "count_smiles", "input": smi, "args": [nb, q, fl], "out": {str(k): v for k, v in cf.GetNonzeroElements().items()}}))
            fp = p.GetAvalonFP(mb, False, nb, q, False, fl)
            print(json.dumps({"kind": "fp_molblock", "input": mb, "args": [nb, q, fl], "out": list(fp.GetOnBits())}))
    recs = "$$$$\n".join(c[2] for c in coords) + "\n"
    out = subprocess.run([sys.argv[1]], input=recs.encode(), capture_output=True, timeout=900).stdout.decode()
    parts = out.split("@@\n")[1:]
    assert len(parts) == len(coords), (len(parts), len(coords))
    for (kind, inp, _), o in zip(coords, parts):
        print(json.dumps({"kind": kind, "input": inp, "out": o}))


if __name__ == "__main__":
    main()
