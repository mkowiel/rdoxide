#!/usr/bin/env python3
"""Builds tests/data/avalon_layout6.jsonl: SMILES laid out up to the fragment stage, then run through the
collision / stretch / bridgehead parts of Avalon's layout.c (CountCollisions, ImproveCollisions,
ScaleByFixedFragments, FixLinearBridgeHeads, FixInwardHydrogens; driver .work/avalon_drv/lay6.c).

    .venv-ref/bin/python tools/make_avalon_layout6_corpus.py DRIVER > tests/data/avalon_layout6.jsonl
"""
import json
import subprocess
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

SMILES = [
    "CC(C)(C)C(C)(C)C(C)(C)C", "CC(C)(C)c1ccc(cc1)C(C)(C)C", "C1CC2CCC1CC2", "C1CC2CCC1C2", "C1CC2CC1C2", "C12C3C4C1C5C2C3C45", "C1C2CC3CC1CC(C2)C3",
    "c1ccc2c(c1)ccc1ccccc12", "C1CCC2(CC1)CCCCC2", "CC(C)(C)C(C)(C)C(C)(C)C(C)(C)C", "C(C(C)(C)C)(C(C)(C)C)C(C)(C)C", "CCCCCCCCCCCCCCCC",
    "N#CC(C#N)(C#N)C(C#N)(C#N)C#N", "C[Si](C)(C)C[Si](C)(C)C", "OC(c1ccccc1)(c1ccccc1)c1ccccc1", "CC1(C)C2CCC1(C)C(=O)C2", "C1CC2CC1CC2C1CCCC1",
    "[H][C@]12CC[C@@]([H])(CC1)C2", "[H][C@@]1(C)CCCC[C@]1([H])C", "[H][C@]12CCC[C@@]1([H])CCCC2", "[H][C@]12CC[C@]([H])(C1)CC2", "C[C@@]12CC[C@](C)(CC1)CC2",
    "[H][C@]12CC[C@@]3([H])CC[C@]([H])(C1)C23", "[H]C12CCC(CC1)C2", "[H][C@@]12C[C@@]3([H])C[C@@]([H])(C1)C[C@@]([H])(C3)C2", "C1CCC2(C1)CCC1(C2)CCCC1",
    "C1CCCCCCCCCCCC1", "C1CCCCCCCCCCCCCCCCC1", "C1=CC=CC=CC=C1", "C1C2C3C4C1C5C2C3C45", "O=C1CCC2(CC1)OCCO2", "CC12CCC3C(CCC4CC(O)CCC34C)C1CCC2O",
    "c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccccc1", "C(c1ccccc1)(c1ccccc1)(c1ccccc1)c1ccccc1", "CCC(CC)(CC)C(CC)(CC)CC", "[Fe](C1C=CC=C1)C1C=CC=C1", "C1CC1C1CC1C1CC1",
    "[H][C@]1(C)CCCC1", "[H][C@@]1(C)CC2CCC1C2", "[H][C@]1(CC)C2CCC1CC2", "[H]C1(C)C2CC3CC(C2)CC1C3", "[H][C@]1(O)C[C@@]2([H])CC[C@]1([H])C2", "[H]C1(N)C2CCC1CC2", "[H]C1(C)C2C3CC(C3)C12",
    "[H][C@]1(C)C[C@@]2(C)CCC[C@@]2(C)C1", "C[C@]1([H])CC[C@@]2(C)CCCC[C@]12[H]", "[H]C1(C)CC2CCC1C2C", "[H]C1(Cl)C2CCC1C2", "[H]C12CCC(C1)C(C)C2", "[H]C1(C)C2CCC3CCC1C23",
    "C12(CCCC1)CCCC2", "C1CC2(C1)CC1(C2)CC2(C1)CC1(C2)C1", "OC1C(O)C(O)C(O)C(O)C1O", "C(=O)(O)C(O)C(O)C(O)C(=O)O",
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
