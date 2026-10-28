#!/usr/bin/env python3
"""Builds tests/data/avalon_struchk.jsonl: InitCheckMol option sets applied to SMILES and molfiles and what Avalon's
struchk.c (RunStruchk) answers: result bits, message list and the checked molecule (driver: .work/avalon_drv/sc.c).
`{D}` in the option sets stands for the directory with checkfgs.chk / checkfgs.trn.

    .venv-ref/bin/python tools/make_avalon_struchk_corpus.py DRIVER DATADIR > tests/data/avalon_struchk.jsonl
"""
import importlib.util
import json
import os
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

OPTS = [
    "-ta {D}/checkfgs.trn|-ca {D}/checkfgs.chk|-cc|-cl 3|-cs|-cn 999|-or",
    "-tm|-ta {D}/checkfgs.trn|-tm|-or|-ca {D}/checkfgs.chk|-cc|-cl 3|-cs|-cn 999|-l /dev/null",
    "-da|-dg|-ds|-cc|-cs|-cn 20",
    "-cn 5|-tm",
    "-cc|-cl 40",
    "-cs|-ca {D}/checkfgs.chk",
]


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def with_props(mb, lines):
    ls = mb.split("\n")
    end = ls.index("M  END")
    return "\n".join(ls[:end] + lines + ls[end:])


def inputs():
    smis = []
    for name in ("make_avalon_layout1_corpus", "make_avalon_layout5_corpus"):
        smis += load(name).SMILES
    smis = smis[::2] + [
        "CC(=O)Oc1ccccc1C(=O)O.Cl", "[Na+].CC(=O)[O-]", "CCO.CCCCC", "O=C(O)C(F)(F)F.CCN", "C[N+](C)(C)C", "[Fe]", "CC(C)(C)c1ccccc1", "N#N", "[H][H]", "C=[N+]=[N-]", "OS(=O)(=O)O",
    ]
    out = []
    for s in smis:
        out.append(("SMILES:" + s, s))
    for s in smis[::3]:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        AllChem.Compute2DCoords(m)
        mb = Chem.MolToMolBlock(m)
        out.append((mb, "2d " + s))
        # collapsed coordinates: atoms clash
        m2 = Chem.Mol(m)
        conf = m2.GetConformer()
        for i in range(m2.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (p.x * 0.3, p.y * 0.3, 0))
        out.append((Chem.MolToMolBlock(m2), "shrunk " + s))
        mh = Chem.AddHs(m)
        if AllChem.EmbedMolecule(mh, randomSeed=3) == 0:
            out.append((Chem.MolToMolBlock(mh), "3d " + s))
        # wedges and either bonds
        m3 = Chem.Mol(m)
        for b in m3.GetBonds():
            if b.GetBondType() == Chem.BondType.DOUBLE and b.GetBeginAtom().GetDegree() > 1 and b.GetEndAtom().GetDegree() > 1:
                b.SetStereo(Chem.BondStereo.STEREOANY)
        out.append((Chem.MolToMolBlock(m3), "either " + s))
        out.append((with_props(mb, ["A    1", "OH-", "A    2", "r1", "G   1    0", "COMPONENT"]), "alias " + s))
    out.append((with_props(Chem.MolToMolBlock(Chem.MolFromSmiles("CCO")), ["A    1", "Xyz", "A    3", "^^+2"]), "alias bad"))
    return out


def main():
    driver, datadir = sys.argv[1], sys.argv[2]
    ins = inputs()
    recs = []
    for k, (body, name) in enumerate(ins):
        for j in (k % len(OPTS), (k + 2) % len(OPTS)):
            recs.append((OPTS[j], body, name))
    inp = "$$$$\n".join(o.replace("{D}", datadir) + "\n" + b + ("\n" if not b.endswith("\n") else "") for o, b, _ in recs)
    out = subprocess.run([driver], input=inp.encode(), capture_output=True, timeout=1800).stdout.decode()
    parts = out.split("@@\n")[1:]
    assert len(parts) == len(recs), (len(parts), len(recs))
    for (o, b, n), p in zip(recs, parts):
        print(json.dumps({"options": o, "input": b, "name": n, "output": p}))


if __name__ == "__main__":
    main()
