"""Tests of rdoxide.Chem.Draw.rdMolDraw2D (SVG drawer, options, helpers)."""

import pytest

from rdoxide import Chem
from rdoxide.Chem.Draw import rdMolDraw2D
from rdoxide.Geometry import Point2D


def svg_of(mol, **kwargs):
    d = rdMolDraw2D.MolDraw2DSVG(300, 250)
    d.DrawMolecule(mol, **kwargs)
    d.FinishDrawing()
    return d.GetDrawingText()


def test_basic_drawing():
    m = Chem.MolFromSmiles("CCO")
    svg = svg_of(m, legend="ethanol")
    assert svg.startswith("<?xml version='1.0' encoding='iso-8859-1'?>")
    assert "width='300px' height='250px'" in svg
    assert svg.rstrip().endswith("</svg>")
    assert "class='bond-0 atom-0 atom-1'" in svg
    d = rdMolDraw2D.MolDraw2DSVG(300, 250)
    assert isinstance(d, rdMolDraw2D.MolDraw2D)
    assert (d.Width(), d.Height()) == (300, 250)


def test_highlights_both_overloads():
    m = Chem.MolFromSmiles("CCO")
    plain = svg_of(m)
    assert "<ellipse" not in plain
    atoms = svg_of(m, highlightAtoms=[0, 1])
    assert atoms.count("<ellipse") == 2
    d = rdMolDraw2D.MolDraw2DSVG(300, 250)
    d.DrawMolecule(m, [0, 1], [0], {0: (1, 0, 0)}, {0: (0, 1, 0)})
    d.FinishDrawing()
    svg = d.GetDrawingText()
    assert "#FF0000" in svg and "#00FF00" in svg
    with pytest.raises(ValueError):
        svg_of(m, highlightAtoms=[7])


def test_flexi_canvas_and_mol_size():
    m = Chem.MolFromSmiles("c1ccccc1CCO")
    d = rdMolDraw2D.MolDraw2DSVG(-1, -1)
    w, h = d.GetMolSize(m)
    d.DrawMolecule(m)
    assert (d.Width(), d.Height()) == (w, h)


def test_draw_molecules_grid():
    ms = [Chem.MolFromSmiles(s) for s in ("CCO", "c1ccccc1", "CCN")]
    d = rdMolDraw2D.MolDraw2DSVG(600, 300, 300, 300)
    d.DrawMolecules(ms, legends=["a", "b", "c"], highlightAtoms=[[0], [], [1]])
    d.FinishDrawing()
    svg = d.GetDrawingText()
    assert svg.count("<rect") == 1
    with pytest.raises(ValueError):
        d.DrawMolecules(ms, legends=["a"], highlightAtoms=[[0]])


def test_options():
    d = rdMolDraw2D.MolDraw2DSVG(200, 200)
    o = d.drawOptions()
    assert o.padding == 0.05 and o.addAtomIndices is False
    o.addAtomIndices = True
    o.padding = 0.1
    assert d.drawOptions().addAtomIndices is True
    assert d.drawOptions().padding == 0.1
    o.backgroundColour = (0.5, 0.5, 0.5)
    assert d.drawOptions().getBackgroundColour() == (0.5, 0.5, 0.5, 1.0)
    with pytest.raises(ValueError):
        o.setHighlightColour((2, 0, 0))
    with pytest.raises(AttributeError):
        o.notAnOption = 1
    o.atomLabels[1] = "foo"
    assert dict(d.drawOptions().atomLabels.items()) == {1: "foo"}
    assert len(o.atomLabels) == 1 and 1 in o.atomLabels
    del o.atomLabels[1]
    assert len(o.atomLabels) == 0
    o.legendPosition = rdMolDraw2D.LegendPosition.Top
    assert o.legendPosition == rdMolDraw2D.LegendPosition.Top
    assert str(o.legendPosition) == "Top"
    # standalone options can be copied over the drawer's
    other = rdMolDraw2D.MolDrawOptions()
    other.rotate = 90
    d.SetDrawOptions(other)
    assert d.drawOptions().rotate == 90
    o.useBWAtomPalette()
    assert d.drawOptions().getAtomPalette() == {-1: (0.0, 0.0, 0.0, 1.0)}
    assert o.fixedScale == -1.0


def test_modes():
    d = rdMolDraw2D.MolDraw2DSVG(200, 200)
    rdMolDraw2D.SetDarkMode(d)
    assert d.drawOptions().backgroundColour == (0.0, 0.0, 0.0, 1.0)
    o = rdMolDraw2D.MolDrawOptions()
    rdMolDraw2D.SetMonochromeMode(o, (1, 0, 0), (0, 0, 1))
    assert o.symbolColour == (1.0, 0.0, 0.0, 1.0)
    assert o.backgroundColour == (0.0, 0.0, 1.0, 1.0)
    rdMolDraw2D.SetACS1996Mode(o, 1.5)
    assert o.bondLineWidth == 0.6 and o.fixedFontSize == 10
    with pytest.raises(ValueError):
        rdMolDraw2D.SetACS1996Mode(o, 0.0)


def test_primitives_and_coordinates():
    d = rdMolDraw2D.MolDraw2DSVG(200, 100)
    d.SetScale(200, 100, Point2D(0, 0), Point2D(2, 1))
    d.SetColour((0, 0, 1))
    d.DrawLine(Point2D(0, 0), Point2D(2, 1))
    d.SetFillPolys(False)
    d.DrawRect(Point2D(10, 10), Point2D(20, 20), rawCoords=True)
    d.DrawEllipse(Point2D(10, 10), Point2D(30, 20), True)
    d.DrawArc(Point2D(50, 50), 10, 0, 90, True)
    d.DrawArrow(Point2D(10, 50), Point2D(80, 50), True, 0.1, 0.5, (1, 0, 0), True)
    d.DrawTriangle(Point2D(0, 0), Point2D(5, 0), Point2D(0, 5), True)
    d.DrawPolygon([Point2D(0, 0), Point2D(5, 0), Point2D(0, 5)], True)
    d.DrawWavyLine(Point2D(0, 0), Point2D(50, 0), (1, 0, 0), (1, 0, 0), 8, 0.5, True)
    d.DrawString("hello", Point2D(50, 50), True)
    d.DrawString("hello", Point2D(50, 50), 1, True)
    with pytest.raises(ValueError):
        d.DrawString("x", Point2D(1, 1), 5)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    assert "stroke:#0000FF" in svg and "<ellipse" in svg
    assert d.LineWidth() == 2.0
    p = d.GetDrawCoords(Point2D(1, 0.5))
    assert 0 < p.x < 200 and 0 < p.y < 100
    assert d.Offset().x == 0


def test_metadata_and_tags():
    m = Chem.MolFromSmiles("CO")
    from rdoxide.Chem import rdDepictor

    rdDepictor.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DSVG(200, 200)
    d.DrawMolecule(m)
    d.AddMoleculeMetadata(m)
    d.TagAtoms(m, events={"onclick": "pick"})
    d.FinishDrawing()
    svg = d.GetDrawingText()
    assert "<rdkit:atom idx=\"1\" atom-smiles=\"[CH3]\"" in svg
    assert "bond-smiles=\"-\"" in svg
    assert "class='atom-selector atom-1'" in svg and "onclick='pick(1);'" in svg


def test_helper_functions():
    m = Chem.MolFromSmiles("C[C@H](F)Cl")
    prepared = rdMolDraw2D.PrepareMolForDrawing(m)
    assert prepared.GetNumConformers() == 1 and m.GetNumConformers() == 0
    with pytest.raises(RuntimeError):
        rdMolDraw2D.PrepareMolForDrawing(None)
    svg = rdMolDraw2D.MolToSVG(m, 200, 200, highlightAtoms=[0])
    assert "<ellipse" in svg and svg.startswith("<?xml")
    d = rdMolDraw2D.MolDraw2DSVG(300, 300)
    rdMolDraw2D.PrepareAndDrawMolecule(d, m, legend="x")
    assert "<path" in d.GetDrawingText()
    acs = rdMolDraw2D.MolToACS1996SVG(m)
    assert "stroke-width:0.6px" in acs
    assert rdMolDraw2D.MeanBondLength(prepared) > 0


def test_native_text():
    m = Chem.MolFromSmiles("CCN")
    d = rdMolDraw2D.MolDraw2DSVG(250, 200, -1, -1, True)
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    assert "fill:#0000FF' >N</text>" in txt
    assert "font-family:sans-serif" in txt


def test_options_from_json():
    o = rdMolDraw2D.MolDrawOptions()
    rdMolDraw2D.UpdateMolDrawOptionsFromJSON(
        o,
        '{"addAtomIndices": 1, "padding": 0.2, "legendPosition": "Top", '
        '"highlightColour": [0, 1, 0], "atomLabels": {"3": "x"}, '
        '"atomColourPalette": "bw", "colourScheme": "dark", '
        '"drawingExtentsInclude": {"BONDS": true, "RADICALS": true}}',
    )
    assert o.addAtomIndices and o.padding == 0.2
    assert o.legendPosition == rdMolDraw2D.LegendPosition.Top
    assert o.highlightColour == (0.0, 1.0, 0.0, 1.0)
    assert dict(o.atomLabels.items()) == {3: "x"}
    assert o.drawingExtentsInclude == 2 | 32
    # the individual colours override the scheme; the palette was set after it
    assert o.backgroundColour == (0.0, 0.0, 0.0, 1.0)
    assert o.getAtomPalette() == {-1: (0.0, 0.0, 0.0, 1.0)}
    d = rdMolDraw2D.MolDraw2DSVG(200, 200)
    rdMolDraw2D.UpdateDrawerParamsFromJSON(d, '{"rotate": 90}')
    assert d.drawOptions().rotate == 90
    rdMolDraw2D.UpdateMolDrawOptionsFromJSON(o, "")
    with pytest.raises(Exception):
        rdMolDraw2D.UpdateMolDrawOptionsFromJSON(o, "{not json")


def test_set_own_draw_options():
    d = rdMolDraw2D.MolDraw2DSVG(100, 100)
    o = d.drawOptions()
    o.padding = 0.2
    d.SetDrawOptions(d.drawOptions())
    assert d.drawOptions().padding == 0.2


def test_contours():
    import numpy as np

    d = rdMolDraw2D.MolDraw2DSVG(200, 200)
    ps = rdMolDraw2D.ContourParams()
    ps.fillGrid = True
    assert ps.gridResolution == 0.15 and len(ps.colourMap) == 3
    ps.contourColour = (1, 0, 0)
    assert ps.contourColour == (1.0, 0.0, 0.0, 1.0)
    rdMolDraw2D.ContourAndDrawGaussians(
        d, [Point2D(0, 0), Point2D(1, 1)], [1.0, -0.5], [0.4, 0.4], nContours=4, params=ps
    )
    d.FinishDrawing()
    assert "stroke:#FF0000" in d.GetDrawingText()
    grid = np.zeros((5, 4))
    grid[2, 2] = 1
    d = rdMolDraw2D.MolDraw2DSVG(200, 200)
    rdMolDraw2D.ContourAndDrawGrid(d, grid, [0, 1, 2, 3, 4], [0, 1, 2, 3], nContours=2)
    assert "<path" in d.GetDrawingText()
    with pytest.raises(ValueError):
        rdMolDraw2D.ContourAndDrawGrid(d, grid, [0, 1, 2], [0, 1, 2, 3])
    with pytest.raises(ValueError):
        rdMolDraw2D.ContourAndDrawGrid(d, [[1.0]], [0], [0])
    with pytest.raises(ValueError):
        rdMolDraw2D.ContourAndDrawGaussians(d, [], [], [])


def test_multi_colour_highlights():
    m = Chem.MolFromSmiles("CCOC")
    d = rdMolDraw2D.MolDraw2DSVG(300, 300)
    red, green = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)
    d.DrawMoleculeWithHighlights(m, "legend", {0: [red], 1: [red, green]}, {0: [green]}, {}, {})
    d.FinishDrawing()
    svg = d.GetDrawingText()
    assert svg.count("<ellipse") == 1 and "#00FF00" in svg
    with pytest.raises(ValueError):
        d.DrawMoleculeWithHighlights(m, "", {9: [red]}, {}, {}, {})


def test_lasso_highlights():
    m = Chem.MolFromSmiles("CCOCC")
    d = rdMolDraw2D.MolDraw2DSVG(300, 300)
    d.drawOptions().multiColourHighlightStyle = rdMolDraw2D.MultiColourHighlightStyle.Lasso
    red, blue = (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)
    d.DrawMoleculeWithHighlights(m, "", {0: [red], 1: [red], 2: [red, blue], 3: [blue]}, {}, {}, {})
    d.FinishDrawing()
    svg = d.GetDrawingText()
    assert "stroke:#FF0000;stroke-width:2." in svg and "<ellipse" not in svg


def test_draw_reaction():
    from rdoxide.Chem import rdChemReactions

    rxn = rdChemReactions.ReactionFromSmarts("[C:1](=[O:2])O.[N:3]>CC(=O)C>[C:1](=[O:2])[N:3]")
    d = rdMolDraw2D.MolDraw2DSVG(800, 300)
    d.DrawReaction(rxn)
    d.FinishDrawing()
    plain = d.GetDrawingText()
    assert plain.count("<rect") == 1 and "class='bond-" in plain
    d = rdMolDraw2D.MolDraw2DSVG(800, 300)
    d.DrawReaction(rxn, highlightByReactant=True, highlightColorsReactants=[(1.0, 0.0, 0.0), (0.0, 0.0, 1.0)])
    d.FinishDrawing()
    assert d.GetDrawingText() != plain and "fill:#FF0000" in d.GetDrawingText()
    # a canvas that sizes itself
    d = rdMolDraw2D.MolDraw2DSVG(-1, -1)
    d.DrawReaction(rxn)
    assert d.Width() > 0 and d.Height() > 0
    # an empty reaction still has an arrow
    d = rdMolDraw2D.MolDraw2DSVG(300, 200)
    d.DrawReaction(rdChemReactions.ChemicalReaction())
    assert "<path" in d.GetDrawingText()


def test_comic_mode_and_unknown_stereo():
    mol = Chem.MolFromSmiles("CC(F)(Cl)Br")

    def svg(**opts):
        d = rdMolDraw2D.MolDraw2DSVG(250, 200)
        for k, v in opts.items():
            setattr(d.drawOptions(), k, v)
        d.DrawMolecule(mol)
        d.FinishDrawing()
        return d.GetDrawingText()

    plain = svg()
    comic = svg(comicMode=True)
    assert comic != plain and "class='bond-0" in comic
    assert svg(unspecifiedStereoIsUnknown=True) != plain


def test_prepare_for_drawing_kekulizes_canonically():
    # same Kekule form as MolOps::KekulizeIfPossible(mol, false) (canonical=true)
    from rdoxide import Chem
    from rdoxide.Chem.Draw import rdMolDraw2D

    m = rdMolDraw2D.PrepareMolForDrawing(Chem.MolFromSmiles("c1ccccc1"), kekulize=True)
    assert [b.GetBondType().name for b in m.GetBonds()] == ["SINGLE", "DOUBLE"] * 3


def test_svg_glyph_contours_end_at_curve_start():
    # a contour whose last segment is a curve back to its start gets no closing line
    from rdoxide import Chem
    from rdoxide.Chem.Draw import rdMolDraw2D

    d = rdMolDraw2D.MolDraw2DSVG(300, 250)
    d.DrawMolecule(Chem.MolFromSmiles("CCO"))
    d.FinishDrawing()
    svg = d.GetDrawingText()
    o = svg.split("<path class='atom-2'")[1].split("fill='#FF0000'")[0]
    lines = o.split("\n")
    assert not any(l.startswith("L ") for l in lines)
