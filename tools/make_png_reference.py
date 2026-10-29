#!/usr/bin/env python3
"""Writes tests/data/png_ref/*.png: what the official wheel's MolDraw2DCairo draws for the cases of
tests/python/test_png_reference.py (CASES is shared).

    PYTHONPATH=.work/pil .venv-ref/bin/python tools/make_png_reference.py
"""
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "python"))
from png_cases import CASES, draw  # noqa: E402
from rdkit.Chem.Draw import rdMolDraw2D  # noqa: E402
from rdkit import Chem  # noqa: E402

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "data", "png_ref")
os.makedirs(out, exist_ok=True)
for name in CASES:
    with open(os.path.join(out, name + ".png"), "wb") as f:
        f.write(draw(Chem, rdMolDraw2D, name))
