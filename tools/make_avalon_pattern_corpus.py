#!/usr/bin/env python3
"""Builds tests/data/avalon_pattern.jsonl (molfile -> what ConvertMoleculeToString, CheckAtoms,
TransformAugmentedAtoms and AlternativeAATransformation of the Avalon C code produce with the
STRUCHK rule files) and tests/data/avalon_aaparse.jsonl (augmented atom strings ->
StringToAugmentedAtom). Driver: .work/avalon_drv/pat.c.

    .venv-ref/bin/python tools/make_avalon_pattern_corpus.py DRIVER RULEDIR OUTDIR
"""
import json
import os
import re
import subprocess
import sys
import tempfile

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

SMILES = [
    "CCO", "C[N+](=O)[O-]", "C[N+](C)(C)C", "[Na+].[Cl-]", "[Na]OC", "CC(=O)[O-].[Na+]", "[O-]C(=O)C", "C[S+](C)C",
    "CS(=O)(=O)C", "CS(C)=O", "C[P+](C)(C)C", "c1ccccc1", "c1ccncc1", "c1cc[nH]c1", "[Fe+2]", "[Fe+3]", "[Cu+]", "[Mg+2]",
    "CC(=O)OC", "N#CC", "C=[N+]=[N-]", "CN=[N+]=[N-]", "[NH4+]", "C[NH3+]", "CC(=O)N", "O=C(O)c1ccccc1", "C1CC1",
    "OP(=O)(O)O", "[O-]P(=O)([O-])[O-]", "CC[Si](C)(C)C", "B(O)(O)C", "[B-](F)(F)(F)F", "ClC(Cl)(Cl)Cl", "[Cl-]",
    "C[Hg]Cl", "[Pt](Cl)(Cl)(N)N", "O=[N+]([O-])c1ccccc1", "C[N+]([O-])=O", "CN(=O)=O", "CC(=O)C=C", "C#C", "[CH2]C",
    "[C-]#[O+]", "O=C=O", "S=C=S", "C[Se]C", "[I-]", "I(=O)(=O)C", "CI(Cl)Cl", "C1CCCCC1", "c1ccc2ccccc2c1", "C=C(C)C=C",
    "[Li]C", "[K+].[O-]C", "[Al](C)(C)C", "[Sn](C)(C)(C)C", "[Zn]", "[Au+3]", "N[C@@H](C)C(=O)O", "C[N+](C)(C)CC(=O)[O-]",
]


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


def run(drv, args, stdin=None):
    return subprocess.run([drv] + args, capture_output=True, input=stdin, timeout=60).stdout.decode()


def aa_strings(ruledir):
    out = set()
    for fn in ("checkfgs.chk", "checkfgs.trn"):
        for line in open(os.path.join(ruledir, fn)):
            out.update(re.findall(r'"([^"]*)"', line))
    out.update([
        "", "C", "@C", "!@C", "@!C", "C+1", "C-3", "C+-", "C+9", "C.", "C|", "C:", "C:1", "C(-N)", "C(=O)(-N)", "C(-N", "C(-N)x",
        "C(*N)", "C(-=N)", "C(?N)", "C(-~N)", "C(~N)", "C(=~N)", "C(#N+1)", "C(-N.)", "C(-N:2)", "C(-N+1.:3)", "C(-N,O,S)",
        "N,O#", "1C", "C(-1)", "C()", "C(-N)(-N)(-N)(-N)(-N)", "O.(=C)", "Fe+2|(-N)(=O:1)", "!@C+1.(-N)", "C(-N)(", "C)",
    ])
    return sorted(out)


def main():
    drv, ruledir, outdir = sys.argv[1:4]
    with open(os.path.join(outdir, "avalon_pattern.jsonl"), "w") as fo:
        for mb in blocks():
            with tempfile.NamedTemporaryFile("wb", suffix=".mol", dir=os.getcwd(), delete=False) as f:
                f.write(mb.encode())
            try:
                out = run(drv, [f.name, ruledir])
            finally:
                os.unlink(f.name)
            fo.write(json.dumps({"input": mb, "output": out}) + "\n")
    strs = aa_strings(ruledir)
    out = run(drv, ["parse"], ("\n".join(strs) + "\n").encode()).split("\n")
    with open(os.path.join(outdir, "avalon_aaparse.jsonl"), "w") as fo:
        for s, o in zip(strs, out):
            fo.write(json.dumps({"input": s, "output": o}) + "\n")


if __name__ == "__main__":
    main()
