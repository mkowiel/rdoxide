#!/usr/bin/env python3
"""Builds tests/data/avalon_pyfp.jsonl: what the official rdkit.Avalon.pyAvalonTools gives for
GetAvalonFP / GetAvalonCountFP / GetAvalonFPAsWords on a set of molecules.

    .venv-ref/bin/python tools/make_avalon_pyfp_corpus.py > tests/data/avalon_pyfp.jsonl
"""
import json
import os
import sys

from rdkit import Chem
from rdkit.Avalon import pyAvalonTools as p
from rdkit.Chem import AllChem

sys.path.insert(0, os.path.dirname(__file__))
from make_avalon_smiles_corpus import SMILES  # noqa: E402

EXTRA = [
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "CC1=C(C(=O)NC(=O)N1)C", "c1ccc2c(c1)[nH]c1ccccc12", "C[C@@H]1CC[C@H]2[C@@H](C1)OC(=O)C2", "O=C(O)c1ccccc1O",
    "CCOC(=O)C1=C(C)NC(C)=C(C1c1ccccc1[N+](=O)[O-])C(=O)OC", "[Pt](Cl)(Cl)(N)N", "C[Sn](C)(C)Cl", "FC(F)(F)c1ccccc1",
    "O=P(O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1O", "CC(C)(C)OC(=O)NC(C)C(=O)O", "S=C=Nc1ccccc1", "C#N", "[2H]O[2H]",
]
CONFIGS = [(512, False, 15761407), (1024, False, 15761407), (512, True, 15761407), (2048, False, 15761407),
           (512, False, 0x7FFF), (512, False, 0xF00000), (256, False, 0x10), (4096, True, 0x7FFF), (1024, False, 0x2000)]


def main():
    for smi in SMILES[::2] + EXTRA:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        variants = [("plain", mol)]
        m2 = Chem.Mol(mol)
        AllChem.Compute2DCoords(m2)
        variants.append(("2d", m2))
        m3 = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(m3, randomSeed=5) == 0:
            variants.append(("3d_h", m3))
        for kind, m in variants:
            for k, (nbits, query, flags) in enumerate(CONFIGS):
                if (k + len(smi)) % 3 and kind != "plain":
                    continue
                rec = {"smiles": smi, "kind": kind, "molblock": Chem.MolToMolBlock(m) if kind != "plain" else None,
                       "nBits": nbits, "isQuery": query, "bitFlags": flags}
                fp = p.GetAvalonFP(m, nBits=nbits, isQuery=query, bitFlags=flags)
                rec["onBits"] = list(fp.GetOnBits())
                rec["counts"] = {i: c for i, c in enumerate(p.GetAvalonCountFP(m, nBits=nbits, isQuery=query, bitFlags=flags)) if c}
                rec["words"] = list(p.GetAvalonFPAsWords(m, nBits=nbits, isQuery=query, bitFlags=flags))
                print(json.dumps(rec))


if __name__ == "__main__":
    main()
