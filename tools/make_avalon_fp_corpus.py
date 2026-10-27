#!/usr/bin/env python3
"""Builds tests/data/avalon_fp.jsonl: molfiles, feature selections and the hashed feature counts
SetFingerprintCountsWithFocus of the Avalon C code produces (driver: .work/avalon_drv/fp.c).

    .venv-ref/bin/python tools/make_avalon_fp_corpus.py DRIVER > tests/data/avalon_fp.jsonl
"""
import json
import os
import random
import subprocess
import sys
import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, os.path.dirname(__file__))
from make_avalon_smiles_corpus import SMILES  # noqa: E402

# (which_bits, as_query, fpflags, nbits)
BIT_SETS = [0x10, 0x04, 0x20, 0x400, 0x100, 0x534]
EXTRA = ["C(C)(C)(C)C", "CC(C)(C)C", "OC(=O)C", "O=[N+]([O-])c1ccccc1", "CS(=O)(=O)N", "ClC(Cl)(Cl)Cl", "C[Si](C)(C)C",
         "[Fe](C)(C)C", "Brc1ccc(I)cc1", "C1CC1(C)C", "C12(CCC1)CCC2", "[Na+].[Cl-]", "OP(=O)(O)O", "c1ccc2c(c1)oc1ccccc12"]


def with_sub(mb, rng):
    lines = mb.split("\n")
    i = lines.index("M  END")
    n = int(lines[3][:3])
    picks = rng.sample(range(1, n + 1), min(3, n))
    ents = "".join(f" {a:3d} {rng.choice([-2, -2, 1, 2, 3]):3d}" for a in picks)
    lines.insert(i, f"M  SUB{len(picks):3d}{ents}")
    return "\n".join(lines)


def main():
    drv = sys.argv[1]
    rng = random.Random(7)
    smis = SMILES[::2] + EXTRA
    for k, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        AllChem.Compute2DCoords(m)
        blocks = [Chem.MolToMolBlock(m)]
        if k % 3 == 0:
            blocks.append(Chem.MolToMolBlock(Chem.AddHs(m)))
        blocks.append(with_sub(Chem.MolToMolBlock(m), rng))
        for mb in blocks:
            for which in BIT_SETS:
                for as_query in (0, 1):
                    for fpflags in (0, 1):
                        nbits = rng.choice([512, 2048, 4096])
                        with tempfile.NamedTemporaryFile("wb", suffix=".mol", dir=os.getcwd(), delete=False) as f:
                            f.write(mb.encode())
                        try:
                            out = subprocess.run(
                                [drv, f.name, str(which), str(as_query), str(fpflags), str(nbits)],
                                capture_output=True, timeout=20,
                            ).stdout.decode()
                        finally:
                            os.unlink(f.name)
                        print(json.dumps({"input": mb, "which": which, "query": as_query, "fpflags": fpflags,
                                          "nbits": nbits, "output": out}))


if __name__ == "__main__":
    main()
