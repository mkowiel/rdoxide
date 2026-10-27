#!/usr/bin/env python3
"""Builds tests/data/avalon_molfile.jsonl: molfiles and what the Avalon C code (MolStr2Mol followed
by MolToMolStr) writes back for them.

    .venv-ref/bin/python tools/make_avalon_molfile_corpus.py DRIVER > tests/data/avalon_molfile.jsonl

DRIVER is a program that reads a molfile (path argument) and prints MolToMolStr(MolStr2Mol(text))
or NULL; it is built from the Avalon sources (.work/avalon_drv/drv.c).
"""
import json
import subprocess
import sys
import tempfile
import os

from rdkit import Chem
from rdkit.Chem import AllChem

SMILES = [
    "CCO", "c1ccccc1O", "C[N+](C)(C)C", "[O-]C(=O)C", "C[C@H](N)C(=O)O", "C/C=C/C", "C/C=C\\C",
    "[13CH3]C", "[2H]C([2H])Cl", "[CH3]", "[Na+].[Cl-]", "C1CC1", "O=S(=O)(N)c1ccccc1", "[Fe+3]",
    "CC(C)(C)c1ccc(O)cc1", "C#CC", "[NH4+]", "C[Si](C)(C)C", "FC(F)(F)C(=O)O", "[14C]", "[18F]C",
    "OC[C@@H](O)[C@H](O)CO", "[CH2]C", "c1ccc2ccccc2c1", "[Se]", "[3H]O", "[N+](=O)[O-]",
]


def blocks():
    for smi in SMILES:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        AllChem.Compute2DCoords(m)
        yield Chem.MolToMolBlock(m)
        yield Chem.MolToV3KMolBlock(m)
        m3 = Chem.AddHs(m)
        if AllChem.EmbedMolecule(m3, randomSeed=7) == 0:
            yield Chem.MolToMolBlock(m3)
    # hand-made edge cases
    yield "empty\n\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n"
    yield "x\n  progname 2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 R#  0  0\nM  END\n"
    yield "x\n\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0\n    1.0000    0.0000    0.0000 A   0  0\n  1  2  1  0\nM  END\n"
    yield "x\n\n\n  1  0  0  0  0  0            999 V2000\n    0.5000    0.5000    0.5000 L   0  0\n  1 F    2   6   7\nM  END\n"
    yield "x\n\n\n  1  0  0  0  0  0            999 V2000\n    0.5000    0.5000    0.5000 C   0  0\nV    1 1.5\nM  SUB  1   1   2\nM  RAD  1   1   2\nM  CHG  1   1   5\nM  END\n"
    yield "bad\n\n\n  x  y\n"
    yield "trunc\n\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0\n"
    yield "crlf\r\n\r\n\r\n  1  0  0  0  0  0  0  0  0  0999 V2000\r\n    0.0000    0.0000    0.0000 Cl  0  0\r\nM  END\r\n"
    yield "unk\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0\nM  FOO  something\nM  END\n"
    yield "iso\n\n\n  2  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0\n    0.0000    0.0000    0.0000 Cl  0  0\nM  ISO  2   1  13   2  36\nM  END\n"


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
