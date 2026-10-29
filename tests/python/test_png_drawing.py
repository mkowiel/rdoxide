"""Tests of the PNG drawer (rdMolDraw2D.MolDraw2DCairo, rendered with tiny-skia) and the PIL based Draw functions."""

import struct

import pytest

from rdoxide import Chem
from rdoxide.Chem import Draw
from rdoxide.Chem.Draw import rdMolDraw2D

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def png_size(png):
    assert png[:8] == PNG_SIGNATURE
    return struct.unpack(">II", png[16:24])


def test_cairo_drawer_makes_png():
    m = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    d = rdMolDraw2D.MolDraw2DCairo(300, 250)
    assert isinstance(d, rdMolDraw2D.MolDraw2D)
    d.DrawMolecule(m, legend="aspirin")
    d.FinishDrawing()
    png = d.GetDrawingText()
    assert isinstance(png, bytes)
    assert png_size(png) == (300, 250)
    assert len(png) > 2000  # not a blank canvas


def test_cairo_drawer_writes_a_file(tmp_path):
    m = Chem.MolFromSmiles("c1ccncc1")
    d = rdMolDraw2D.MolDraw2DCairo(200, 200)
    d.DrawMolecule(m)
    d.FinishDrawing()
    path = tmp_path / "pyridine.png"
    d.WriteDrawingText(str(path))
    assert path.read_bytes() == d.GetDrawingText()


def test_png_grid_and_highlights():
    ms = [Chem.MolFromSmiles(s) for s in ("CCO", "c1ccccc1", "CC(=O)O")]
    d = rdMolDraw2D.MolDraw2DCairo(600, 200, 200, 200)
    d.DrawMolecules(ms, highlightAtoms=[[0], [1, 2], []])
    d.FinishDrawing()
    assert png_size(d.GetDrawingText()) == (600, 200)


def test_png_is_deterministic():
    m = Chem.MolFromSmiles("CCN")

    def make():
        d = rdMolDraw2D.MolDraw2DCairo(150, 150)
        d.DrawMolecule(m)
        d.FinishDrawing()
        return d.GetDrawingText()

    assert make() == make()


def test_molecule_image_with_pil():
    pil = pytest.importorskip("PIL.Image")
    m = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    img = Draw.MolToImage(m, size=(200, 200))
    assert isinstance(img, pil.Image)
    assert img.size == (200, 200)
    grid = Draw.MolsToGridImage([m, Chem.MolFromSmiles("c1ccncc1")], molsPerRow=2, subImgSize=(200, 200))
    assert grid.size == (400, 200)
    # a drawing is mostly white with some coloured ink
    colours = img.convert("RGB").getcolors(maxcolors=100000)
    assert max(colours)[1] == (255, 255, 255)
    assert any(c[1][0] > 200 and c[1][1] < 80 for c in colours)  # the red oxygens


def test_png_carries_the_molecule():
    m = Chem.MolFromSmiles("CC(=O)Nc1ccccc1")
    d = rdMolDraw2D.MolDraw2DCairo(250, 200)
    d.DrawMolecule(m)
    d.FinishDrawing()
    nm = Chem.MolFromPNGString(d.GetDrawingText())
    assert Chem.MolToSmiles(nm) == Chem.MolToSmiles(m)
