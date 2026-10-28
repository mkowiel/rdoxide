#!/usr/bin/env python3
"""Builds tests/data/coordgen_model.txt: molecules (RDKit 2D coordinates) and what the model classes of coordgenlibs
3.0.2 answer (driver tools/coordgen_drv/model.cpp). Each record is the input block, the driver output and "$$".

    .venv-ref/bin/python tools/make_coordgen_model_corpus.py DRIVER > tests/data/coordgen_model.txt
"""
import importlib.util
import os
import struct
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fb(x):
    return "%08x" % struct.unpack("<I", struct.pack("<f", x))[0]


EXTRA = [
    "C", "CC", "C1CC1", "C1CCCCC1", "c1ccccc1", "c1ccc2ccccc2c1", "C1CC2CCC1C2", "C1CCC2(CC1)CCCC2", "c1ccc2c(c1)[nH]c1ccccc12",
    "C1CCCCCCCCCC1", "C1CCCCCCCCCCCCCCC1", "O=C1CCCCCCCCCCCN1", "CC(C)(C)c1ccccc1", "[NH4+].[Cl-]", "CS(=O)(=O)O", "OP(=O)(O)O",
    "C1=CC=CC=C1C=CC1=CC=CC=C1", "C/C=C/C", "C/C=C\\C", "C1CC1C1CC1", "c1ccc(cc1)-c1ccccc1", "N#CC", "C[N+](C)(C)C", "O=C=O",
    "C1C2CC3CC1CC(C2)C3", "c1cc2cccc3c4cccc5cccc(c(c1)c23)c54", "CC1=C(C)C(=O)C(C)=C(C)C1=O", "[Fe](C)(C)(C)C", "FC(F)(F)C(F)(F)F",
    "OCC(O)C(O)C(O)C(O)CO", "C1CCC2C(C1)CCC1C2CCC2CCCC21", "N1C=CC=C1", "c1ccsc1", "c1ccoc1", "S1SSSSS1", "[SiH3]C", "O=S(=O)(N)c1ccccc1",
]


def smiles():
    out = list(EXTRA)
    for name in ("make_avalon_layout1_corpus", "make_avalon_layout5_corpus"):
        out += load(name).SMILES[::3]
    return out


def block(m):
    conf = m.GetConformer()
    lines = ["MOL %d %d" % (m.GetNumAtoms(), m.GetNumBonds())]
    for a in m.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append("A %d %d %s %s" % (a.GetAtomicNum(), a.GetFormalCharge(), fb(p.x * 50), fb(p.y * 50)))
    for b in m.GetBonds():
        o = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3, Chem.BondType.AROMATIC: 1}.get(b.GetBondType(), 1)
        lines.append("B %d %d %d" % (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), o))
    lines.append("END")
    return "\n".join(lines) + "\n"


def main():
    driver = sys.argv[1]
    blocks = []
    for s in smiles():
        m = Chem.MolFromSmiles(s)
        if m is None or m.GetNumAtoms() == 0:
            continue
        if m.GetNumAtoms() > 60:
            continue
        Chem.Kekulize(m, clearAromaticFlags=True)
        AllChem.Compute2DCoords(m)
        blocks.append(block(m))
    inp = "".join(blocks)
    out = subprocess.run([driver], input=inp.encode(), capture_output=True, timeout=600).stdout.decode()
    recs = out.split("$$\n")[:-1]
    assert len(recs) == len(blocks), (len(recs), len(blocks))
    for b, r in zip(blocks, recs):
        sys.stdout.write(b + r + "$$\n")


if __name__ == "__main__":
    main()
