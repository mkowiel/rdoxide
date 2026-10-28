#!/usr/bin/env python3
"""Builds tests/data/avalon_layout5.jsonl: SMILES (read with SMIToMOL, so atom parities and double
bond marks come from the Avalon code) laid out up to the chain/fragment stage, and what the stereo
part of Avalon's layout.c (LayoutBondStereo, LayoutAtomStereo, ClearDBStereoInSmallRings,
ClearFlipFlags, MakeLandscape) produces (driver: .work/avalon_drv/lay5.c, -ffp-contract=off).

    .venv-ref/bin/python tools/make_avalon_layout5_corpus.py DRIVER > tests/data/avalon_layout5.jsonl
"""
import json
import subprocess
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

SMILES = [
    "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", "OC[C@@H](O)[C@H](O)CO", "C/C=C/C", "C/C=C\\C", "F/C=C/C=C/Cl", "F/C=C\\C=C/Cl", "C/C=C/C=C\\C",
    "C[C@H]1CC[C@@H](C)CC1", "C[C@H]1CC[C@H](C)CC1", "[C@H](F)(Cl)Br", "C[C@](F)(Cl)Br", "CC(C)C[C@H](N)C(O)=O", "N[C@@H](Cc1ccccc1)C(=O)O",
    "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O", "C1CCCC/C=C/CC1", "C1CCCC/C=C\\CC1", "C1CCCCCCC/C=C/CC1", "C1CCCCCCC/C=C\\CC1",
    "C1CC/C=C/CC1", "C1CC/C=C\\CC1", "C/C=C/C1CC1", "C[C@H](O)[C@@H](O)C", "C[C@H](Cl)C[C@@H](C)Cl", "N[C@@H](C)C(=O)N[C@H](C(=O)O)Cc1c[nH]c2ccccc12",
    "C[C@@H]1CCCC[C@H]1C", "CC[C@H](C)[C@@H](C)CC", "OC(=O)[C@H](O)[C@@H](O)C(=O)O", "F[C@](Cl)(Br)I", "C[C@@H](C(=O)O)[C@H](N)C", "C1C[C@@H]2CC[C@H]1C2",
    "C[C@H]1C[C@@H](C)C[C@H](C)C1", "CN1CCC[C@H]1c1cccnc1", "O=C(/C=C/c1ccccc1)O", "c1ccccc1/C=C\\C(=O)O", "CC/C=C/C=C/C(=O)O", "C/C(F)=C(/Cl)Br",
    "FC(F)(F)[C@@H](O)c1ccccc1", "C[C@H](N)[C@@H](O)c1ccccc1", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "C[C@@]1(O)CC[C@@H](C(C)C)CC1", "CC(C)=C[C@@H]1C[C@H]1C(=O)O",
    "C[C@H](Br)/C=C/C", "[C@@H]1(F)CCCC1", "[C@H]1(F)CCCC1Cl", "C[S@](=O)c1ccccc1", "C[P@](=O)(c1ccccc1)CC", "C=C=C", "CC=C=CC", "C/C=C=C/C",
    "CCCC/C=C/CCCC", "C1CCCCCCCCCC/C=C/1", "C/C=C/C(/C)=C/C", "OC(=O)/C=C\\C(=O)O", "OC(=O)/C=C/C(=O)O", "C[C@H](/C=C/C)[C@@H](C)O",
]


def main():
    smis = []
    for smi in SMILES:
        smis.append(smi)
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            for _ in range(2):
                smis.append(Chem.MolToSmiles(m, doRandom=True))
    seen, uniq = set(), []
    for s in smis:
        if s not in seen:
            seen.add(s)
            uniq.append(s)
    out = subprocess.run([sys.argv[1]], input=("\n".join(uniq) + "\n").encode(), capture_output=True, timeout=600).stdout.decode()
    parts = out.split("@@\n")[1:]
    assert len(parts) == len(uniq), (len(parts), len(uniq))
    for s, o in zip(uniq, parts):
        print(json.dumps({"input": s, "output": o}))


if __name__ == "__main__":
    main()
