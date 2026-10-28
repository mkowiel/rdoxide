#!/usr/bin/env python3
"""Builds tests/data/avalon_layout7.jsonl: SMILES laid out up to the fragment stage, then run through the
fragment bond flip parts of Avalon's layout.c (TopBottomSeparation, AttachmentCenterDistance,
ImproveFragmentByBondFlip, ImproveFragmentByAllBondFlips, FloodInvertColor, ChangeAtomColors; driver .work/avalon_drv/lay7.c).

    .venv-ref/bin/python tools/make_avalon_layout7_corpus.py DRIVER > tests/data/avalon_layout7.jsonl
"""
import json
import subprocess
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

SMILES = [
    "C1CCCCC1CC(C)CC1CCCCC1", "c1ccccc1-c1ccccc1", "c1ccc(cc1)C(C)(C)c1ccccc1", "C1CCC1C1CCCC1C1CCC1", "c1ccccc1CCc1ccccc1CCc1ccccc1", "C1CC1C(C)C1CCCCC1",
    "c1ccc2ccccc2c1-c1ccc2ccccc2c1", "C1CCCC1-C1CCCC1-C1CCCC1-C1CCCC1", "c1ccccc1C(=O)Nc1ccccc1C(=O)Nc1ccccc1", "OC1CCC(CC1)C1CCC(O)CC1", "C1CC2CCC1C2C1CC2CCC1C2",
    "c1ccccc1Oc1ccccc1Oc1ccccc1", "CC1CCC(CC1)C(C)(C)C1CCC(C)CC1", "c1cccnc1-c1ccccn1", "c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccccc1", "C1CCC1C2CCC2.C1CCCC1",
    "CC(C)(C)C1CCC(CC1)C(C)(C)C", "CN(C)c1ccc(cc1)C(c1ccc(cc1)N(C)C)c1ccc(cc1)N(C)C", "C1CCCCC1C1CCCCC1C1CCCCC1C1CCCCC1", "c1ccccc1C1CCCC(C1)c1ccccc1", "C1=CC=CC=C1C=CC1=CC=CC=C1",
    "CCCC", "C", "CC1CC1", "C1CCC1CC1CCC1", "O=C(c1ccccc1)c1ccccc1", "c1ccccc1S(=O)(=O)c1ccccc1", "c1ccc(cc1)-c1ccccc1-c1ccccc1-c1ccccc1",
]


def main():
    smis = []
    for smi in SMILES:
        smis.append(smi)
        m = Chem.MolFromSmiles(smi)
        if m is not None:
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
