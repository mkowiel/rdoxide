//! Tests of the SVG drawer.
#![allow(clippy::field_reassign_with_default)]

use std::collections::BTreeMap;

use rdoxide::geometry::Point2D;
use rdoxide::graphmol::depictor::{Compute2DCoordParameters, compute_2d_coords};
use rdoxide::graphmol::mol_draw2d::draw_mol::DrawMolHighlights;
use rdoxide::graphmol::mol_draw2d::drawer::{MolDraw2DSVG, MoleculeHighlights};
use rdoxide::graphmol::mol_draw2d::helpers::{DrawColour, OrientType};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use rdoxide::graphmol::{ROMol, RWMol};

fn smi(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn drawn(smiles: &str, w: i32, h: i32) -> (MolDraw2DSVG, String) {
    let mut d = MolDraw2DSVG::new(w, h, -1, -1).unwrap();
    d.draw_molecule(&smi(smiles), "", -1).unwrap();
    d.finish_drawing();
    let svg = d.get_drawing_text();
    (d, svg)
}

#[test]
fn simple_molecule() {
    let (_, svg) = drawn("CCO", 300, 250);
    assert!(svg.starts_with("<?xml version='1.0' encoding='iso-8859-1'?>\n<svg version='1.1'"));
    assert!(
        svg.contains(
            "width='300px' height='250px' viewBox='0 0 300 250'>\n<!-- END OF HEADER -->\n"
        )
    );
    assert!(svg.contains("<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300.0' height='250.0' x='0.0' y='0.0'> </rect>\n"));
    assert!(svg.ends_with("</svg>\n"));
    // two bonds, split in two for the colours of the C-O bond
    assert_eq!(svg.matches("class='bond-").count(), 3);
    assert!(svg.contains("class='bond-0 atom-0 atom-1' d="));
    assert!(svg.contains("fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px"));
    // the oxygen label is drawn as glyph paths of atom-2, in red
    assert!(svg.contains("class='atom-2' d='M"));
    assert!(svg.contains("fill='#FF0000'"));
}

#[test]
fn background_can_be_left_out() {
    let mut d = MolDraw2DSVG::new(200, 200, -1, -1).unwrap();
    d.draw_options_mut().clear_background = false;
    d.draw_molecule(&smi("CC"), "", -1).unwrap();
    assert!(!d.get_drawing_text().contains("<rect"));
}

#[test]
fn flexible_canvas_takes_size_from_molecule() {
    let mut d = MolDraw2DSVG::new(-1, -1, -1, -1).unwrap();
    d.draw_molecule(&smi("c1ccccc1CCO"), "", -1).unwrap();
    d.finish_drawing();
    assert!(d.width() > 0 && d.height() > 0);
    let svg = d.get_drawing_text();
    assert!(svg.contains(&format!(
        "width='{}px' height='{}px'",
        d.width(),
        d.height()
    )));
    assert_eq!(d.panel_width(), d.width());
}

#[test]
fn legend_and_highlights() {
    let m = smi("CCO");
    let mut plain = MolDraw2DSVG::new(300, 300, -1, -1).unwrap();
    plain.draw_molecule(&m, "", -1).unwrap();
    let mut legend = MolDraw2DSVG::new(300, 300, -1, -1).unwrap();
    legend.draw_molecule(&m, "ethanol", -1).unwrap();
    assert!(legend.get_drawing_text().len() > plain.get_drawing_text().len());

    let mut hl = MolDraw2DSVG::new(300, 300, -1, -1).unwrap();
    let mut colours = BTreeMap::new();
    colours.insert(1, DrawColour::new(1.0, 0.0, 0.0));
    hl.draw_molecule_atoms(&m, "", &[0, 1], Some(&colours), None, -1)
        .unwrap();
    let svg = hl.get_drawing_text();
    // the highlighted atoms get ellipses, the bond between them a line
    assert_eq!(svg.matches("<ellipse").count(), 2);
    assert!(svg.contains("fill:#FF0000;fill-rule:evenodd;"));
    assert_eq!(plain.get_drawing_text().matches("<ellipse").count(), 0);
}

#[test]
fn radicals_are_spots() {
    let (_, svg) = drawn("[CH2]C", 300, 300);
    // one spot, a filled polygon of atom-0
    let spots = svg
        .lines()
        .filter(|l| l.contains("class='atom-0'") && l.contains("fill-opacity:1"))
        .count();
    assert_eq!(spots, 1, "{svg}");
}

#[test]
fn several_molecules_in_a_grid() {
    let a = smi("CCO");
    let b = smi("c1ccccc1");
    let mols: Vec<Option<&ROMol>> = vec![Some(&a), None, Some(&b)];
    let mut d = MolDraw2DSVG::new(400, 200, 200, 200).unwrap();
    let legends = vec!["a".to_string(), "x".to_string(), "b".to_string()];
    d.draw_molecules(&mols, Some(&legends), None, None).unwrap();
    d.finish_drawing();
    let svg = d.get_drawing_text();
    // one background for the whole canvas
    assert_eq!(svg.matches("<rect").count(), 1);
    // a three-molecule list with a hole still has a size to match
    assert!(
        d.draw_molecules(&mols, Some(&legends[..2]), None, None)
            .is_err()
    );
    // the second panel's drawing is in the right half of the canvas
    let mut d1 = MolDraw2DSVG::new(400, 200, 200, 200).unwrap();
    d1.draw_molecules(&[Some(&b)], None, None, None).unwrap();
    let mut d2 = MolDraw2DSVG::new(400, 200, 200, 200).unwrap();
    d2.draw_molecules(&[Some(&a), Some(&b)], None, None, None)
        .unwrap();
    assert_ne!(d1.get_drawing_text(), d2.get_drawing_text());
    assert!(d2.set_active_mol_idx(1).is_ok());
    assert!(d2.set_active_mol_idx(2).is_err());
    let x0 = d2.get_draw_coords_atom(0).unwrap().x;
    assert!(x0 > 200.0, "{x0}");
}

#[test]
fn highlights_for_several() {
    let a = smi("CCO");
    let mut d = MolDraw2DSVG::new(200, 200, -1, -1).unwrap();
    let hl = [MoleculeHighlights {
        atoms: Some(vec![0, 1]),
        ..Default::default()
    }];
    d.draw_molecules(&[Some(&a)], None, Some(&hl), None)
        .unwrap();
    // continuous highlighting adds the bond between the atoms
    assert!(d.get_drawing_text().contains("<ellipse"));
}

#[test]
fn coordinates() {
    let (d, _) = drawn("CCO", 300, 300);
    let c = d.get_draw_coords_atom(2).unwrap();
    assert!(c.x > 0.0 && c.x < 300.0 && c.y > 0.0 && c.y < 300.0);
    let m = d.get_atom_coords(c).unwrap();
    let back = d.get_draw_coords(m).unwrap();
    assert!((back.x - c.x).abs() < 1e-6 && (back.y - c.y).abs() < 1e-6);
    // the atom coordinates are those of the drawing frame, with y inverted
    let a = d.get_atom_coords_atom(2).unwrap();
    assert_eq!(a, d.atom_coords().unwrap()[2]);
    assert!((a.y + m.y).abs() < 1e-6);
    assert_eq!(d.atom_coords().unwrap().len(), 3);
    assert_eq!(d.atom_syms().unwrap()[2].0, "OH");
    assert_eq!(d.atom_syms().unwrap().len(), 3);
    assert!(d.range().unwrap().x > 0.0);
    // sizes of strings in molecule units
    let (w, h) = d.get_string_size("OH").unwrap();
    assert!(w > 0.0 && h > 0.0);
    let (x0, y0, x1, y1) = d
        .get_string_extremes("OH", OrientType::W, Point2D::new(0.0, 0.0))
        .unwrap();
    assert!(x1 > x0 && y1 > y0);
    let (lw, lh) = d.get_label_size("NH2", OrientType::N).unwrap();
    assert!(lw > 0.0 && lh > 0.0);
}

#[test]
fn explicit_scale_and_primitives() {
    let mut d = MolDraw2DSVG::new(200, 100, -1, -1).unwrap();
    d.set_scale_for_range(
        200,
        100,
        Point2D::new(0.0, 0.0),
        Point2D::new(2.0, 1.0),
        None,
    )
    .unwrap();
    // the primitives need a molecule to have been drawn first only for the
    // active scale, the transformation is the global one
    d.set_colour(DrawColour::new(0.0, 0.0, 1.0));
    d.draw_line(Point2D::new(0.0, 0.0), Point2D::new(2.0, 1.0), false)
        .unwrap();
    d.draw_line(Point2D::new(10.0, 10.0), Point2D::new(20.0, 10.0), true)
        .unwrap();
    d.set_dash(vec![2.0, 2.0]);
    d.set_fill_polys(false);
    d.draw_rect(Point2D::new(1.0, 1.0), Point2D::new(5.0, 5.0), true)
        .unwrap();
    d.draw_ellipse(Point2D::new(1.0, 1.0), Point2D::new(5.0, 5.0), true)
        .unwrap();
    d.draw_wavy_line(
        Point2D::new(0.0, 0.0),
        Point2D::new(50.0, 0.0),
        DrawColour::new(1.0, 0.0, 0.0),
        DrawColour::new(1.0, 0.0, 0.0),
        8,
        0.5,
        true,
    )
    .unwrap();
    d.draw_arrow(
        Point2D::new(0.0, 0.0),
        Point2D::new(50.0, 0.0),
        true,
        0.1,
        0.5,
        DrawColour::new(0.0, 0.0, 0.0),
        true,
    )
    .unwrap();
    d.draw_plus(
        Point2D::new(50.0, 50.0),
        4,
        DrawColour::new(0.0, 0.0, 0.0),
        true,
    )
    .unwrap();
    d.draw_string("hi", Point2D::new(50.0, 50.0), true).unwrap();
    d.finish_drawing();
    let svg = d.get_drawing_text();
    assert!(svg.contains("stroke:#0000FF"));
    assert!(svg.contains("stroke-dasharray:2,2"));
    assert!(svg.contains("<ellipse cx='3.0' cy='3.0' rx='2.0' ry='2.0'"));
    assert!(svg.contains(" C"), "wavy line");
    // unfilled rectangles are not closed
    assert!(svg.contains("fill:none;stroke:#0000FF"));
    assert!((d.range().unwrap().x - 2.2).abs() < 1e-9);
}

#[test]
fn metadata_and_tags() {
    let mut m = smi("CO");
    compute_2d_coords(&mut m, &Compute2DCoordParameters::default()).unwrap();
    let mut d = MolDraw2DSVG::new(200, 200, -1, -1).unwrap();
    d.draw_molecule(&m, "", -1).unwrap();
    d.add_molecule_metadata(&m, -1).unwrap();
    let mut events = BTreeMap::new();
    events.insert("onclick".to_string(), "pick".to_string());
    d.tag_atoms(&m, 0.2, &events).unwrap();
    d.finish_drawing();
    let svg = d.get_drawing_text();
    assert!(svg.contains(
        "<metadata>\n<rdkit:mol xmlns:rdkit = \"http://www.rdkit.org/xml\" version=\"0.9\">\n"
    ));
    assert!(svg.contains("<rdkit:atom idx=\"1\" atom-smiles=\"[CH3]\" drawing-x=\""));
    assert!(svg.contains("<rdkit:atom idx=\"2\" atom-smiles=\"[OH]\""));
    assert!(svg.contains(
        "<rdkit:bond idx=\"1\" begin-atom-idx=\"1\" end-atom-idx=\"2\" bond-smiles=\"-\" />\n"
    ));
    assert!(svg.contains("</rdkit:mol></metadata>\n"));
    assert_eq!(
        svg.matches("class='bond-selector bond-0 atom-0 atom-1'")
            .count(),
        1
    );
    assert_eq!(svg.matches("class='atom-selector atom-").count(), 2);
    assert!(svg.contains(" onclick='pick(1);'"));
}

#[test]
fn classes_follow_the_active_class() {
    let mut d = MolDraw2DSVG::new(200, 200, -1, -1).unwrap();
    d.set_active_class("mine");
    d.draw_molecule(&smi("CO"), "", -1).unwrap();
    let svg = d.get_drawing_text();
    assert!(svg.contains("class='mine bond-0  atom-0 atom-1'") || svg.contains("class='mine"));
    assert!(svg.contains("class='mine atom-1'"));
}

#[test]
fn varied_molecules_draw() {
    for s in [
        "C[C@H](F)Cl",
        "c1ccccc1",
        "C/C=C/C",
        "[Na+].[Cl-]",
        "C1CC1.CCO",
        "C#N",
        "[13CH3][NH3+]",
        "OC(=O)[C@@H](N)C",
        "*C",
        "[O-][N+](=O)c1ccccc1",
    ] {
        let (_, svg) = drawn(s, 300, 300);
        assert!(svg.ends_with("</svg>\n"), "{s}");
        assert!(svg.contains("<path"), "{s}");
    }
}

#[test]
fn unsupported_options_are_errors() {
    let mut d = MolDraw2DSVG::new(100, 100, -1, -1).unwrap();
    d.draw_options_mut().comic_mode = true;
    assert!(d.draw_molecule(&smi("CC"), "", -1).is_err());
}

#[test]
fn mol_size() {
    let mut d = MolDraw2DSVG::new(-1, -1, -1, -1).unwrap();
    let (w, h) = d
        .get_mol_size(&smi("CCCCO"), "", DrawMolHighlights::default(), -1)
        .unwrap();
    assert!(w > 0 && h > 0);
}

#[test]
fn native_text_writes_text_elements() {
    let mut d = MolDraw2DSVG::with_text(
        250,
        200,
        -1,
        -1,
        rdoxide::graphmol::mol_draw2d::helpers::MolDrawOptions::default(),
        true,
    )
    .unwrap();
    d.draw_molecule(&smi("CCN"), "", -1).unwrap();
    d.finish_drawing();
    let svg = d.get_drawing_text();
    assert!(svg.contains("fill:#0000FF' >N</text>"), "{svg}");
    assert!(svg.contains("font-family:sans-serif;text-anchor:start;"));
    assert!(!svg.contains("class='atom-2' d='M"));
}

#[test]
fn options_from_json() {
    use rdoxide::graphmol::mol_draw2d::helpers::{LegendPosition, MolDrawOptions};
    use rdoxide::graphmol::mol_draw2d::options_json::update_mol_draw_options_from_json;
    let mut o = MolDrawOptions::default();
    update_mol_draw_options_from_json(
        &mut o,
        r#"{"addAtomIndices": true, "bondLineWidth": 3.5, "legendPosition": "Left",
            "symbolColour": [0.5, 0.25, 1], "fontFile": "BuiltinRobotoRegular",
            "atomColourPalette": {"7": [0, 0, 1]}}"#,
    )
    .unwrap();
    assert!(o.add_atom_indices);
    assert_eq!(o.bond_line_width, 3.5);
    assert_eq!(o.legend_position, LegendPosition::Left);
    assert_eq!(o.symbol_colour, DrawColour::new(0.5, 0.25, 1.0));
    assert_eq!(o.font_file, "BuiltinRobotoRegular");
    assert_eq!(o.atom_colour_palette[&7], DrawColour::new(0.0, 0.0, 1.0));
    // untouched options keep their values; bad JSON is an error
    assert_eq!(o.padding, 0.05);
    assert!(update_mol_draw_options_from_json(&mut o, "{oops").is_err());
    update_mol_draw_options_from_json(&mut o, "").unwrap();
}

#[test]
fn multi_colour_highlights() {
    let red = DrawColour::new(1.0, 0.0, 0.0);
    let green = DrawColour::new(0.0, 1.0, 0.0);
    let blue = DrawColour::new(0.0, 0.0, 1.0);
    let m = smi("CCOC");
    let draw = |fill: bool| {
        let mut d = MolDraw2DSVG::new(300, 300, -1, -1).unwrap();
        d.draw_options_mut().fill_highlights = fill;
        let mut atoms = BTreeMap::new();
        atoms.insert(0, vec![red]);
        atoms.insert(1, vec![red, green]);
        atoms.insert(2, vec![red, green, blue]);
        let mut bonds = BTreeMap::new();
        bonds.insert(0, vec![red]);
        bonds.insert(1, vec![green, blue]);
        d.draw_molecule_with_highlights(&m, "", atoms, bonds, BTreeMap::new(), BTreeMap::new(), -1)
            .unwrap();
        d.finish_drawing();
        d.get_drawing_text()
    };
    let filled = draw(true);
    // one ellipse for a single colour, arcs (polygons) for several
    assert_eq!(filled.matches("<ellipse").count(), 1);
    assert!(filled.contains("fill:#FF0000;"));
    assert!(filled.contains("fill:#00FF00;"));
    assert!(filled.contains("fill:#0000FF;"));
    let open = draw(false);
    assert_eq!(open.matches("<ellipse").count(), 1);
    assert!(open.contains("fill:none;fill-rule:evenodd;stroke:#00FF00"));
    assert_ne!(filled, open);
    // out of range indices are errors
    let mut d = MolDraw2DSVG::new(300, 300, -1, -1).unwrap();
    let mut bad = BTreeMap::new();
    bad.insert(17, vec![red]);
    assert!(
        d.draw_molecule_with_highlights(
            &m,
            "",
            bad.clone(),
            BTreeMap::new(),
            BTreeMap::new(),
            BTreeMap::new(),
            -1
        )
        .is_err()
    );
    assert!(
        d.draw_molecule_with_highlights(
            &m,
            "",
            BTreeMap::new(),
            bad,
            BTreeMap::new(),
            BTreeMap::new(),
            -1
        )
        .is_err()
    );
}
