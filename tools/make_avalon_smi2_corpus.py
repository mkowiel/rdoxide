#!/usr/bin/env python3
"""Builds tests/data/avalon_smiles_atoms.jsonl (SmilesToMDLAtom and QuickCheck of the Avalon C code
for SMILES strings) and tests/data/avalon_smiles_valence.jsonl (ComputeSmilesValenceAndHCount for
molfiles). Driver: .work/avalon_drv/smi2.c.

    .venv-ref/bin/python tools/make_avalon_smi2_corpus.py DRIVER OUTDIR
"""
import itertools
import json
import os
import random
import re
import subprocess
import sys
import tempfile

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

SMILES = [
    "CCO", "c1ccccc1O", "C[N+](C)(C)C", "[O-]C(=O)C", "C[C@H](N)C(=O)O", "C/C=C/C", "[13CH3]C", "[2H]C([2H])Cl", "[CH3]",
    "[Na+].[Cl-]", "C1CC1", "O=S(=O)(N)c1ccccc1", "[Fe+3]", "C#CC", "[NH4+]", "[3H]O", "[N+](=O)[O-]", "c1cc[nH]c1",
    "[O-][n+]1ccccc1", "C12C3C4C1C5C2C3C45", "[C@H](F)(Cl)Br", "C[P](C)(C)=O", "CS(=O)C", "O=s1cccc1", "[Se]", "[14C]", "[18F]C",
    "c1ccc2ccccc2c1", "C[C@@]12CC[C@H]3C(CC[C@@H]4CCCC[C@@]34C)[C@@H]1CC[C@@H]2O",
]


def run(drv, args, stdin=None):
    return subprocess.run([drv] + args, capture_output=True, input=stdin, timeout=60).stdout.decode()


def atom_strings():
    out = set()
    for smi in SMILES:
        out.add(smi)
        out.update(re.findall(r"\[[^\]]*\]", smi))
    syms = ["C", "c", "N", "n", "O", "Fe", "Na", "H", "*", "R", "Cl", "Se", "Xx", "se", "Uuo", "D", "T"]
    for iso, sym, st, h, ch, mp in itertools.product(["", "13", "2", "3", "99", "100"], syms, ["", "@", "@@", "@2", "@@@"], ["", "H", "H2", "H9"], ["", "+", "-", "++", "--", "+2", "-3", "+12", "-10", "+-"], ["", ":5", ":12"]):
        out.add("[%s%s%s%s%s%s]" % (iso, sym, st, h, ch, mp))
    rng = random.Random(5)
    out.update(["[", "[]", "[C", "[Cx]", "[C+]]", "[C@@H+:3]CC", "[123]", "[1]", "[C@12]", "[-]", "[+]", "C", "", "[C]C", "[*:3]", "[*+3]", "[*-1]", "[*-2]", "[*-3]", "[12*]", "[R-1]", "[R7]", "[C++-]", "[C+--]"])
    rest = sorted(out)
    rng.shuffle(rest)
    return rest[:2500] + ["CCO)", "C(C", "C1CC", "C1CC1", "1CC", "C%12CC%12", "C%1CC", "C[", "C]", "C{C}", "C{C", "C=1CC1", "C-1CC1", "C:1CC1", "C(1)CC1", "C1CC1|0,0", "c1ccccc1|1,2", "[C-1]1CC1", "C1[CH2]C1", "C|", "C%99CC%99", "C%9CC", "C#1CC1", "C/1CC1", "C.1CC1"]


def with_hcount(mb, rng):
    lines = mb.split("\n")
    na = int(lines[3][0:3])
    for i in range(na):
        if rng.random() < 0.5:
            l = lines[4 + i]
            lines[4 + i] = l[:42] + str(rng.choice([0, 1, 2, 3, 4, 5])).rjust(3) + l[45:]
    return "\n".join(lines)


def main():
    drv, outdir = sys.argv[1:3]
    strs = atom_strings()
    outs = run(drv, ["atoms"], ("\n".join(strs) + "\n").encode()).split("\n")
    with open(os.path.join(outdir, "avalon_smiles_atoms.jsonl"), "w") as fo:
        for s, o in zip(strs, outs):
            fo.write(json.dumps({"input": s, "output": o}) + "\n")
    rng = random.Random(9)
    with open(os.path.join(outdir, "avalon_smiles_valence.jsonl"), "w") as fo:
        for smi in SMILES + ["O=[N+]([O-])c1ccccc1", "C[S+](C)C", "[O-][N+]=O", "C[Si](C)(C)C", "FC(F)(F)F", "[Cl-]", "[Cu+2]", "C=C=C", "O=P(O)(O)O", "OS(=O)(=O)O"]:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            AllChem.Compute2DCoords(m)
            for mb in (Chem.MolToMolBlock(m), Chem.MolToMolBlock(m, kekulize=False), Chem.MolToMolBlock(Chem.AddHs(m))):
                for variant in (mb, with_hcount(mb, rng), with_hcount(mb, rng)):
                    with tempfile.NamedTemporaryFile("wb", suffix=".mol", dir=os.getcwd(), delete=False) as f:
                        f.write(variant.encode())
                    try:
                        out = run(drv, [f.name])
                    finally:
                        os.unlink(f.name)
                    fo.write(json.dumps({"input": variant, "output": out}) + "\n")


if __name__ == "__main__":
    main()
