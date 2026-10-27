#!/usr/bin/env python3
"""Builds tests/data/inchi_polymers.jsonl: polymer molfiles (SRU / CRU / MON / COP / ...) with the
InChI, return code, message and AuxInfo the official RDKit gives for several option sets.

Run with an interpreter that has the official `rdkit` (RDKit 2026.09.1):
    .venv-ref/bin/python tools/make_inchi_polymer_corpus.py > tests/data/inchi_polymers.jsonl
"""
import json
import random

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdinchi

OPTIONS = ["", "-Polymers", "-Polymers -NPZz", "-Polymers -FoldCRU", "-Polymers -NoFoldCRU", "-Polymers105"]

# (SMILES, [(sgroup type, atoms, bonds, connect)], comment)
BASE = [
    "CCCOC", "CCC(C)CC", "OCCOCCOCCO", "CC(c1ccccc1)CC", "C(=O)CCCCC(=O)NCCCCCCN", "CC(C(=O)OC)CC",
    "CCCCCC", "O=C(O)CCN", "CC(F)(F)CC", "c1ccc(cc1)CCC", "CCOCCOCC", "NCC(=O)NCC(=O)NCC(=O)O",
    "C[Si](C)(C)OC", "CCC(Cl)CC", "CC(O)CC(O)CC",
]


def build(smi, rng):
    m = Chem.MolFromSmiles(smi)
    rdDepictor.Compute2DCoords(m)
    n = m.GetNumAtoms()
    bonds = [b.GetIdx() for b in m.GetBonds()]
    # a contiguous run of heavy atoms of the backbone as the repeating unit
    size = rng.randint(1, max(1, min(3, n - 2)))
    start = rng.randint(1, max(1, n - size - 1))
    atoms = list(range(start, min(n - 1, start + size)))
    if not atoms:
        return None
    inside = set(atoms)
    xb = [b.GetIdx() for b in m.GetBonds() if (b.GetBeginAtomIdx() in inside) != (b.GetEndAtomIdx() in inside)]
    if len(xb) != 2:
        return None
    stype = rng.choice(["SRU", "SRU", "SRU", "MON", "CRU", "COP", "MER"])
    sg = Chem.CreateMolSubstanceGroup(m, stype if stype in ("SRU", "MON", "COP", "MER", "CRU") else "SRU")
    for a in atoms:
        sg.AddAtomWithIdx(a)
    for b in xb:
        sg.AddBondWithIdx(b)
    sg.SetProp("LABEL", rng.choice(["n", "m", "2-5", "3"]))
    sg.SetProp("CONNECT", rng.choice(["HT", "HH", "EU"]))
    return Chem.MolToMolBlock(m)


def main():
    rng = random.Random(20261027)
    seen = set()
    for smi in BASE:
        for _ in range(14):
            mb = build(smi, rng)
            if mb is None or mb in seen:
                continue
            seen.add(mb)
            res = {}
            for opt in OPTIONS:
                r = rdinchi.MolBlockToInchi(mb, opt)
                res[opt] = [r[0], r[1], r[2], r[3], r[4]]
            print(json.dumps({"molblock": mb, "results": res}))


if __name__ == "__main__":
    main()
