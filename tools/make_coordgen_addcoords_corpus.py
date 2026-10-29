#!/usr/bin/env python3
"""Builds tests/data/coordgen_addcoords.txt: molecules and the coordinates RDKit's rdCoordGen.AddCoords makes of them
(plain, with a coordMap, and with a template molecule).

    .venv-ref/bin/python tools/make_coordgen_addcoords_corpus.py > tests/data/coordgen_addcoords.txt
"""
import importlib.util
import os

from rdkit import Chem, Geometry, RDLogger
from rdkit.Chem import rdCoordGen

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, path))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


sk = load("make_coordgen_sketch_corpus.py", "sk")
mc = load("make_coordgen_macro_corpus.py", "mc")


def emit(m, tag, extra=()):
    c = m.GetConformer()
    print("MOL", tag)
    for line in extra:
        print(line)
    for i in range(m.GetNumAtoms()):
        p = c.GetAtomPosition(i)
        print("P %d %r %r" % (i, p.x, p.y))
    print("END")


def main():
    smis = sk.EXTRA + mc.EXTRA[::3]
    n = 0
    for s in smis:
        m = Chem.MolFromSmiles(s)
        if m is None or m.GetNumAtoms() > 50:
            continue
        rdCoordGen.AddCoords(m)
        emit(m, s)
        n += 1
    # coordinate maps: the first three atoms of a chain fixed on a line
    for s in ["CCCCC", "c1ccccc1CCO", "CC(=O)Nc1ccccc1", "C1CCCCC1C(=O)O", "OCC1OC(O)C(O)C(O)C1O"]:
        m = Chem.MolFromSmiles(s)
        cm = {0: Geometry.Point2D(0.0, 0.0), 1: Geometry.Point2D(1.5, 0.0), 2: Geometry.Point2D(2.25, 1.3)}
        p = rdCoordGen.CoordGenParams()
        p.SetCoordMap(cm)
        rdCoordGen.AddCoords(m, p)
        emit(m, s, ["MAP %d %r %r" % (k, v.x, v.y) for k, v in cm.items()])
    # templates: a substructure with 2D coordinates
    for s, t in [("CC(=O)Nc1ccccc1", "c1ccccc1N"), ("OCCc1ccc(O)cc1", "c1ccccc1"), ("CCC1CCCCC1C", "C1CCCCC1")]:
        m = Chem.MolFromSmiles(s)
        tm = Chem.MolFromSmiles(t)
        rdCoordGen.AddCoords(tm)
        # rotate and shift the template so that it differs from the free layout
        conf = tm.GetConformer()
        for i in range(tm.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, Geometry.Point3D(-p.y + 3.0, p.x - 1.0, 0.0))
        p = rdCoordGen.CoordGenParams()
        p.SetTemplateMol(tm)
        rdCoordGen.AddCoords(m, p)
        extra = ["TPL " + t]
        for i in range(tm.GetNumAtoms()):
            q = conf.GetAtomPosition(i)
            extra.append("TP %d %r %r" % (i, q.x, q.y))
        emit(m, s, extra)


if __name__ == "__main__":
    main()
