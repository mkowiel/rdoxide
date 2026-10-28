#!/usr/bin/env python3
"""Builds tests/data/avalon_strip.jsonl: molfiles and what SplitIonicBonds, ForceOctet,
RemoveInorganicFragments, StripMolecule and ScrambleMolecule of the Avalon C code produce
(driver: .work/avalon_drv/canon2.c).

    .venv-ref/bin/python tools/make_avalon_strip_corpus.py DRIVER > tests/data/avalon_strip.jsonl
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
    "CCO", "C[N+](=O)[O-]", "C[N+](C)(C)C", "[Na]OC", "C[Mg]Cl", "CC(=O)O[Na]", "[Na+].[Cl-]", "[Na]Cl.CC", "O=C([O-])C.[K+]",
    "C[S+](C)C", "CS(=O)(=O)C", "C[P+](C)(C)C", "c1ccccc1", "[Fe+2].[Cl-].[Cl-]", "C=[N+]=[N-]", "CN=[N+]=[N-]", "[NH4+]",
    "C[NH3+]", "CC(=O)N", "[B-](F)(F)(F)F.C[N+](C)(C)C", "ClC(Cl)(Cl)Cl", "O=[N+]([O-])c1ccccc1", "CN(=O)=O", "N#CC", "[C-]#[O+]",
    "O=C=O.CCO", "[Li]C", "[Al](C)(C)C.O", "N[C@@H](C)C(=O)O", "C[C@H](F)[C@@H](F)C", "[1*]CC", "[1*]C[2*]", "[N,O]C", "[Cl,Br]CC",
    "C[N+]([O-])=O", "O=[Cu]", "[Pt](Cl)(Cl)(N)N.CC", "[OH-].[Na+]", "C[O-].[Na+]", "[O-][O-]", "OP(=O)([O-])[O-].[Ca+2]",
    "N[C@H](C(=O)O)Cc1ccccc1", "C[C@H](N)[C@@H](O)C", "F[C@](Cl)(Br)I", "C1CC1[C@H](C)N", "[Si](C)(C)(C)C", "C[Mg]O",
]


def with_parities(mb, rng):
    lines = mb.split("\n")
    na, nb = int(lines[3][0:3]), int(lines[3][3:6])
    deg = [0] * na
    for k in range(4 + na, 4 + na + nb):
        a, b = int(lines[k][0:3]) - 1, int(lines[k][3:6]) - 1
        deg[a] += 1
        deg[b] += 1
    for i in range(na):
        if deg[i] >= 3 and rng.random() < 0.7:
            l = lines[4 + i]
            lines[4 + i] = l[:39] + str(rng.choice([1, 2])).rjust(3) + l[42:]
    return "\n".join(lines)


def blocks():
    rng = random.Random(11)
    for smi in SMILES:
        m = Chem.MolFromSmiles(smi) if "," not in smi and "*" not in smi else Chem.MolFromSmiles(smi) or Chem.MolFromSmarts(smi)
        if m is None:
            continue
        AllChem.Compute2DCoords(m)
        Chem.WedgeMolBonds(m, m.GetConformer())
        mb = Chem.MolToMolBlock(m)
        yield mb
        yield with_parities(mb, rng)
        mh = Chem.AddHs(m) if m.GetNumAtoms() < 40 else m
        AllChem.Compute2DCoords(mh)
        yield with_parities(Chem.MolToMolBlock(mh), rng)


def special_blocks():
    """R-group labels (M  RGP) and atom lists (M  ALS), which RDKit does not write for plain molecules"""
    for smi in ("[*]CCO", "[*]C(C)C[*]", "[*]C(F)(Cl)CC", "[*]CC(=O)O[Na]", "[*]C([*])CCCO"):
        m = Chem.MolFromSmiles(smi)
        AllChem.Compute2DCoords(m)
        mb = Chem.MolToMolBlock(m)
        lines = mb.split("\n")
        na = int(lines[3][0:3])
        rg = []
        for i in range(na):
            l = lines[4 + i]
            if l[31:34].strip() == "*":
                lines[4 + i] = l[:31] + "R# " + l[34:]
                rg.append(i + 1)
        props = "M  RGP%3d" % len(rg) + "".join("%4d%4d" % (a, k + 1) for k, a in enumerate(rg))
        yield "\n".join(lines).replace("M  END", props + "\nM  END")
        # atom list on the first carbon
        l2 = list(lines)
        for i in range(na):
            if l2[4 + i][31:34].strip() == "C":
                l2[4 + i] = l2[4 + i][:31] + "L  " + l2[4 + i][34:]
                als = "M  ALS %3d%3d F N   O   " % (i + 1, 2)
                yield "\n".join(l2).replace("M  END", als + "\nM  END")
                break


def main():
    drv = sys.argv[1]
    for mb in list(blocks()) + list(special_blocks()):
        with tempfile.NamedTemporaryFile("wb", suffix=".mol", dir=os.getcwd(), delete=False) as f:
            f.write(mb.encode())
        try:
            out = subprocess.run([drv, f.name], capture_output=True, timeout=20).stdout.decode()
        finally:
            os.unlink(f.name)
        print(json.dumps({"input": mb, "output": out}))


if __name__ == "__main__":
    main()
