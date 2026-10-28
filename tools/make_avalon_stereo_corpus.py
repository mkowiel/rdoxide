#!/usr/bin/env python3
"""Builds tests/data/avalon_stereo.jsonl: molfiles and what the stereo checks of the Avalon C code
(CheckStereo, NoParityDefined, DubiousStereochemistry, AtomClash, CisTransPerception,
FixDubious3DMolecule, RemoveDubiousStereochemistry; driver: .work/avalon_drv/stereo.c) report for them.

    .venv-ref/bin/python tools/make_avalon_stereo_corpus.py DRIVER > tests/data/avalon_stereo.jsonl
"""
import json
import os
import random
import subprocess
import sys
import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem

SMILES = [
    "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", "OC[C@@H](O)[C@H](O)CO", "C/C=C/C", "C/C=C\\C", "F/C=C/C=C/Cl",
    "C[C@H]1CC[C@@H](C)CC1", "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O", "[C@H](F)(Cl)Br",
    "C[C@](F)(Cl)Br", "CC(C)C[C@H](N)C(O)=O", "C=C=C", "CC=C=CC", "C/C=C=C/C", "CC1(C)C2CCC1(C)C(=O)C2",
    "C[S@](=O)c1ccccc1", "C[P@](=O)(c1ccccc1)CC", "C[N+](C)(CC)CCC", "C[Si](C)(CC)CCC", "C1CC1", "CCO",
    "CN1CCC[C@H]1c1cccnc1", "N[C@@H](Cc1ccccc1)C(=O)O", "C[C@H](O)[C@@H](O)C", "CC(C)(C)C", "FC(F)(F)C(F)(F)F",
    "C12C3C4C1C5C2C3C45", "OC1C(O)C(O)C(O)C(O)C1O", "C#CC", "c1ccccc1C(C)(C)C", "C[N@](CC)CCC", "C(F)(Cl)(Br)I",
    "O=C1C=CC(=O)C=C1", "CC(=CC)C", "C/C(F)=C(/Cl)Br", "[O-][N+](=O)C(C)(C)C",
]


def mutate(mb, rng):
    """random stereo symbols on bonds to hit the illegal representations"""
    lines = mb.split("\n")
    na, nb = int(lines[3][0:3]), int(lines[3][3:6])
    for k in range(4 + na, 4 + na + nb):
        if rng.random() < 0.35:
            sym = rng.choice([1, 6, 4, 1, 6, 0])
            l = lines[k]
            lines[k] = l[:9] + str(sym).rjust(3) + l[12:]
    return "\n".join(lines)


def shift_atom(mb, rng):
    """moves one atom onto another atom or onto a bond"""
    lines = mb.split("\n")
    na = int(lines[3][0:3])
    if na < 3:
        return mb
    i, j = rng.sample(range(na), 2)
    xi, yi = float(lines[4 + j][0:10]), float(lines[4 + j][10:20])
    d = rng.choice([0.0, 0.05, 0.1, 0.3])
    lines[4 + i] = "%10.4f%10.4f" % (xi + d, yi) + lines[4 + i][20:]
    return "\n".join(lines)


def tilt(m, rng):
    """flat 2D coordinates with a random z, wedges kept"""
    c = m.GetConformer()
    for a in range(m.GetNumAtoms()):
        p = c.GetAtomPosition(a)
        c.SetAtomPosition(a, (p.x, p.y, rng.uniform(-0.8, 0.8)))
    return m


def blocks():
    rng = random.Random(7)
    for smi in SMILES:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        AllChem.Compute2DCoords(m)
        Chem.WedgeMolBonds(m, m.GetConformer())
        mb = Chem.MolToMolBlock(m)
        yield mb
        yield mb.replace("  0  0  0  0  0  0  0  0999 V2000", "  0  0  1  0  0  0  0  0999 V2000")
        for _ in range(4):
            yield mutate(mb, rng)
        yield shift_atom(mb, rng)
        yield shift_atom(mutate(mb, rng), rng)
        t = Chem.Mol(m)
        Chem.WedgeMolBonds(t, t.GetConformer())
        yield Chem.MolToMolBlock(tilt(t, rng))
        mh = Chem.AddHs(m)
        if AllChem.EmbedMolecule(mh, randomSeed=5) == 0:
            yield Chem.MolToMolBlock(mh)
            yield Chem.MolToMolBlock(Chem.RemoveHs(mh))
        m2 = Chem.Mol(m)
        if AllChem.EmbedMolecule(m2, randomSeed=9) == 0:
            Chem.AssignChiralTypesFromBondDirs(m2)
            Chem.WedgeMolBonds(m2, m2.GetConformer())
            yield Chem.MolToMolBlock(m2)


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
