//! Tests of the `DrawMol` core: coordinates, labels, scale, legend.
#![allow(clippy::field_reassign_with_default)]

use rdoxide::geometry::Point2D;
use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::mol_draw2d::draw_mol::{DrawMol, DrawMolHighlights};
use rdoxide::graphmol::mol_draw2d::draw_text::DrawTextFT;
use rdoxide::graphmol::mol_draw2d::helpers::{LegendPosition, MolDrawOptions};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn smi(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

fn text() -> DrawTextFT {
    DrawTextFT::new(40.0, 6.0, "").unwrap()
}

fn build(
    smiles: &str,
    w: i32,
    h: i32,
    legend: &str,
    opts: &MolDrawOptions,
) -> (DrawMol, DrawTextFT) {
    let mut t = text();
    let mut dm = DrawMol::new(
        &smi(smiles),
        legend,
        w,
        h,
        opts,
        &mut t,
        DrawMolHighlights::default(),
        true,
        -1,
        false,
    )
    .unwrap();
    dm.extract_all(1.0, &mut t).unwrap();
    dm.calc_mean_bond_length().unwrap();
    dm.calculate_scale(&mut t);
    (dm, t)
}

#[test]
fn coordinates_and_labels() {
    let (dm, _) = build("CCO", 300, 300, "", &MolDrawOptions::default());
    assert_eq!(dm.draw_mol.num_conformers(), 1);
    assert_eq!(dm.at_cds.len(), 3);
    // y is inverted relative to the conformer
    let c = dm.draw_mol.conformer(-1).unwrap().atom_pos(2).unwrap();
    assert_eq!(dm.at_cds[2], Point2D::new(c.x, -c.y));
    // only the oxygen is labelled
    assert!(dm.atom_labels[0].is_none() && dm.atom_labels[1].is_none());
    assert_eq!(dm.atom_labels[2].as_ref().unwrap().symbol, "OH");
    assert_eq!(dm.atomic_nums, vec![6, 6, 8]);
    assert!(dm.mean_bond_length > 0.5);
}

#[test]
fn rotation() {
    let mut opts = MolDrawOptions::default();
    let (plain, _) = build("CO", 300, 300, "", &opts);
    opts.rotate = 90.0;
    let (rot, _) = build("CO", 300, 300, "", &opts);
    // a quarter turn about the origin maps (x, y) to (-y, x)
    let (p, r) = (plain.at_cds[1], rot.at_cds[1]);
    assert!((r.x + p.y).abs() < 1e-9 && (r.y - p.x).abs() < 1e-9);
}

#[test]
fn scale_and_draw_coords() {
    // bonds are not extracted yet, so use a molecule where every atom has a
    // label to get extremes that cover all atoms
    let (mut dm, mut t) = build("NO", 300, 300, "", &MolDrawOptions::default());
    assert!(dm.scale > 1.0);
    assert!(dm.x_range > 1e-4 && dm.y_range > 1e-4);
    dm.finish_create_draw_objects(&mut t).unwrap();
    assert!(dm.drawing_initialised);
    // drawing coordinates of the atoms fall inside the canvas, and map back
    for i in 0..2 {
        let d = dm.atom_draw_coords(i).unwrap();
        assert!(
            d.x > 0.0 && d.x < 300.0 && d.y > 0.0 && d.y < 300.0,
            "{d:?}"
        );
        let back = dm.atom_coords_from_screen(d);
        let at = dm.atom_coords(i).unwrap();
        assert!((back.x - at.x).abs() < 1e-9 && (back.y - at.y).abs() < 1e-9);
    }
    assert!(dm.atom_draw_coords(7).is_err());
}

#[test]
fn flexi_canvas() {
    let (dm, _) = build("CCO", -1, -1, "", &MolDrawOptions::default());
    assert!(dm.width > 0 && dm.height > 0);
    assert_eq!(dm.mol_height, dm.draw_height);
    let mut o = MolDrawOptions::default();
    o.fixed_font_size = 12;
    let (dm, t) = build("CCO", 300, 300, "", &o);
    assert!((dm.font_scale - 12.0 / t.base_font_size()).abs() < 1e-12);
}

#[test]
fn legend_layout() {
    let (mut dm, mut t) = build("CCO", 300, 300, "ethanol", &MolDrawOptions::default());
    assert!(dm.legend_height > 0 && dm.mol_height + dm.legend_height == dm.draw_height);
    dm.finish_create_draw_objects(&mut t).unwrap();
    assert_eq!(dm.legends.len(), 1);
    // the legend sits at the bottom, centred
    let leg = &dm.legends[0];
    assert!(
        leg.pos.y > 150.0 && (leg.pos.x - 150.0).abs() < 1.0,
        "{:?}",
        leg.pos
    );

    let mut o = MolDrawOptions::default();
    o.legend_position = LegendPosition::Top;
    let (mut dm, mut t) = build("CCO", 300, 300, "a\nb", &o);
    dm.finish_create_draw_objects(&mut t).unwrap();
    assert_eq!(dm.legends.len(), 2);
    assert!(dm.legends[0].pos.y < dm.legends[1].pos.y && dm.legends[1].pos.y < 150.0);

    o.legend_position = LegendPosition::Left;
    o.legend_vertical_text = true;
    let (mut dm, mut t) = build("CCO", 300, 300, "abc", &o);
    assert!(dm.legend_width > 0 && dm.mol_width < dm.draw_width);
    dm.finish_create_draw_objects(&mut t).unwrap();
    // vertical text: one glyph per line
    assert_eq!(dm.legends.len(), 3);
}

#[test]
fn offsets_and_reset() {
    let (mut dm, mut t) = build("CCO", 300, 300, "", &MolDrawOptions::default());
    dm.finish_create_draw_objects(&mut t).unwrap();
    let before = dm.atom_labels[2].as_ref().unwrap().cds;
    dm.set_offsets(10.0, 20.0, &mut t).unwrap();
    let after = dm.atom_labels[2].as_ref().unwrap().cds;
    assert!((after.x - before.x - 10.0).abs() < 1e-9 && (after.y - before.y - 20.0).abs() < 1e-9);
    dm.reset_everything(&mut t);
    assert!(dm.at_cds.is_empty() && dm.atom_labels.is_empty());
    assert_eq!(dm.scale, 1.0);
}

#[test]
fn atom_indices_and_without_mol() {
    let mut o = MolDrawOptions::default();
    o.add_atom_indices = true;
    o.add_bond_indices = true;
    let (dm, _) = build("CO", 300, 300, "", &o);
    use rdoxide::rdgeneral::RDProps;
    assert_eq!(
        dm.draw_mol
            .atom(1)
            .get_prop("atomNote")
            .unwrap()
            .to_rd_string(),
        "1"
    );
    assert_eq!(
        dm.draw_mol
            .bond(0)
            .get_prop("bondNote")
            .unwrap()
            .to_rd_string(),
        "0"
    );
    let mut t = text();
    let dm = DrawMol::without_mol(
        200,
        100,
        &MolDrawOptions::default(),
        &mut t,
        0.0,
        4.0,
        1.0,
        3.0,
        25.0,
        25.0,
    );
    // y extremes are inverted and swapped
    assert_eq!((dm.y_min, dm.y_max), (-3.0, -1.0));
    assert_eq!(dm.x_range, 4.0);
    let mut o = MolDrawOptions::default();
    o.unspecified_stereo_is_unknown = true;
    let r = DrawMol::new(
        &smi("CO"),
        "",
        100,
        100,
        &o,
        &mut t,
        DrawMolHighlights::default(),
        true,
        -1,
        false,
    );
    assert!(r.is_err());
}

fn kinds(dm: &DrawMol) -> Vec<&'static str> {
    use rdoxide::graphmol::mol_draw2d::draw_shape::ShapeKind::*;
    dm.bonds
        .iter()
        .map(|b| match b.kind {
            Arrow { .. } => "arrow",
            Ellipse => "ellipse",
            SimpleLine { .. } => "line",
            PolyLine { .. } => "poly",
            SolidWedge { .. } => "solid",
            DashedWedge { .. } => "dashed",
            WavyLine { .. } => "wavy",
            Arc { .. } => "arc",
        })
        .collect()
}

#[test]
fn bond_shapes() {
    let o = MolDrawOptions::default();
    // a single bond between like atoms is one line; between unlike atoms two
    // half lines, one per colour
    let (dm, _) = build("CC", 300, 300, "", &o);
    assert_eq!(kinds(&dm), vec!["line"]);
    let (dm, _) = build("CN", 300, 300, "", &o);
    assert_eq!(kinds(&dm), vec!["line", "line"]);
    assert_ne!(dm.bonds[0].line_colour, dm.bonds[1].line_colour);
    // each half line stays on its own side of the midpoint
    let mid = (dm.bonds[0].points[1].x - dm.bonds[1].points[0].x).abs();
    assert!(mid < 1e-9);
    // a triple bond has three lines
    let (dm, _) = build("CC#CC", 300, 300, "", &o);
    assert_eq!(dm.bonds.iter().filter(|b| b.bond == 1).count(), 3);
    // a double bond has two
    let (dm, _) = build("CC=CC", 300, 300, "", &o);
    assert_eq!(dm.bonds.iter().filter(|b| b.bond == 1).count(), 2);
    let ls: Vec<_> = dm.bonds.iter().filter(|b| b.bond == 1).collect();
    // the two lines are parallel and offset by the double bond offset
    let d = |s: &rdoxide::graphmol::mol_draw2d::draw_shape::DrawShape| s.points[1] - s.points[0];
    let cross = d(ls[0]).x * d(ls[1]).y - d(ls[0]).y * d(ls[1]).x;
    assert!(cross.abs() < 1e-6);
}

#[test]
fn ring_bonds_and_aromatics() {
    let o = MolDrawOptions::default();
    let (dm, _) = build("c1ccccc1", 300, 300, "", &o);
    // kekulized with the aromatic flags kept: 3 double bonds of 2 lines and 3
    // single ones of 1, plus the little joins between the single lines
    let lines = dm
        .bonds
        .iter()
        .filter(|b| {
            matches!(
                b.kind,
                rdoxide::graphmol::mol_draw2d::draw_shape::ShapeKind::SimpleLine { .. }
            )
        })
        .count();
    assert_eq!(lines, 9);
}

#[test]
fn wedges_wavy_dative_zero() {
    let o = MolDrawOptions::default();
    let (dm, _) = build("C[C@H](N)C(=O)O", 300, 300, "", &o);
    let k = kinds(&dm);
    assert!(k.contains(&"solid") || k.contains(&"dashed"), "{k:?}");
    let mut m = smi("CC");
    m.bond_mut(0)
        .set_bond_dir(rdoxide::graphmol::bond::BondDir::Unknown);
    let mut t = text();
    let mut o2 = MolDrawOptions::default();
    o2.prepare_mols_before_drawing = false;
    let mut dm = DrawMol::new(
        &m,
        "",
        300,
        300,
        &o2,
        &mut t,
        DrawMolHighlights::default(),
        true,
        -1,
        false,
    )
    .unwrap();
    dm.extract_all(1.0, &mut t).unwrap();
    assert_eq!(kinds(&dm), vec!["wavy"]);

    let (dm, _) = build("N->[Cu]", 300, 300, "", &o);
    assert!(kinds(&dm).contains(&"arrow"), "{:?}", kinds(&dm));
    let (dm, _) = build("[Na+].[Cl-]", 300, 300, "", &o);
    assert!(dm.bonds.is_empty());
}

#[test]
fn query_bonds() {
    use rdoxide::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
    let q = mol_from_smarts("C~C", &SmartsParserParams::default())
        .unwrap()
        .unwrap();
    let mut t = text();
    let mut o = MolDrawOptions::default();
    o.prepare_mols_before_drawing = false;
    let mut dm = DrawMol::new(
        &q,
        "",
        300,
        300,
        &o,
        &mut t,
        DrawMolHighlights::default(),
        true,
        -1,
        false,
    )
    .unwrap();
    dm.extract_all(1.0, &mut t).unwrap();
    // the generic query bond is dotted, in the query colour
    assert_eq!(dm.bonds.len(), 1);
    assert_eq!(dm.bonds[0].line_colour, o.query_colour);
    let q = mol_from_smarts("CC", &SmartsParserParams::default())
        .unwrap()
        .unwrap();
    let mut dm = DrawMol::new(
        &q,
        "",
        300,
        300,
        &o,
        &mut t,
        DrawMolHighlights::default(),
        true,
        -1,
        false,
    )
    .unwrap();
    dm.extract_all(1.0, &mut t).unwrap();
    assert_eq!(dm.bonds.len(), 3);
}

fn build_hl(smiles: &str, opts: &MolDrawOptions, hl: DrawMolHighlights) -> (DrawMol, DrawTextFT) {
    let mut t = text();
    let mut dm = DrawMol::new(
        &smi(smiles),
        "",
        300,
        300,
        opts,
        &mut t,
        hl,
        true,
        -1,
        false,
    )
    .unwrap();
    dm.extract_all(1.0, &mut t).unwrap();
    dm.calculate_scale(&mut t);
    (dm, t)
}

#[test]
fn continuous_bond_and_atom_highlights() {
    use rdoxide::graphmol::mol_draw2d::draw_shape::ShapeKind;
    let o = MolDrawOptions::default();
    // an isolated highlighted bond is a rectangle; atoms get ellipses
    let hl = DrawMolHighlights {
        atoms: vec![0],
        bonds: vec![0],
        ..Default::default()
    };
    let (dm, _) = build_hl("CCC", &o, hl);
    let polys: Vec<_> = dm
        .highlights
        .iter()
        .filter(|h| matches!(h.kind, ShapeKind::PolyLine { .. }))
        .collect();
    let ells: Vec<_> = dm
        .highlights
        .iter()
        .filter(|h| matches!(h.kind, ShapeKind::Ellipse))
        .collect();
    assert_eq!(polys.len(), 1);
    assert_eq!(polys[0].points.len(), 4);
    assert_eq!(ells.len(), 1);
    assert_eq!(ells[0].line_colour, o.highlight_colour);
    // two adjacent highlighted bonds are mitred: 4 points each, sharing the
    // join, and the hull is convex (no butterfly)
    let hl = DrawMolHighlights {
        bonds: vec![0, 1],
        ..Default::default()
    };
    let (dm, _) = build_hl("CCC", &o, hl);
    assert_eq!(dm.highlights.len(), 2);
    // a branch point with three highlighted bonds gives triangular ends
    let hl = DrawMolHighlights {
        bonds: vec![0, 1, 2],
        ..Default::default()
    };
    let (dm, _) = build_hl("CC(C)C", &o, hl);
    assert_eq!(dm.highlights.len(), 3);
    assert!(dm.highlights.iter().any(|h| h.points.len() > 4));
}

#[test]
fn circle_highlights_radii_and_colours() {
    use rdoxide::graphmol::mol_draw2d::draw_shape::ShapeKind;
    let mut o = MolDrawOptions::default();
    o.continuous_highlight = false;
    let col = rdoxide::graphmol::mol_draw2d::helpers::DrawColour::new(0.0, 1.0, 0.0);
    let hl = DrawMolHighlights {
        atoms: vec![0, 1],
        atom_map: [(1, col)].into(),
        radii: [(0, 0.7)].into(),
        ..Default::default()
    };
    let (dm, _) = build_hl("CCO", &o, hl);
    assert_eq!(dm.highlights.len(), 2);
    assert!(matches!(dm.highlights[0].kind, ShapeKind::Ellipse));
    assert_eq!(dm.highlights[0].points[1], Point2D::new(0.7, 0.7));
    assert_eq!(dm.highlights[0].line_colour, o.highlight_colour);
    assert_eq!(
        dm.highlights[1].points[1],
        Point2D::new(o.highlight_radius, o.highlight_radius)
    );
    assert_eq!(dm.highlights[1].line_colour, col);
    // a highlighted atom with circles behind it is drawn in the colour of carbon
    assert_eq!(
        dm.get_colour(0).unwrap(),
        rdoxide::graphmol::mol_draw2d::helpers::DrawColour::new(0.0, 0.0, 0.0)
    );
}

#[test]
fn labelled_atom_highlight_surrounds_label() {
    let o = MolDrawOptions::default();
    let hl = DrawMolHighlights {
        atoms: vec![1],
        ..Default::default()
    };
    let (dm, _) = build_hl("CO", &o, hl);
    let e = &dm.highlights[0];
    let label = dm.atom_labels[1].as_ref().unwrap();
    // the ellipse is at least as big as the default radius and centred on the label
    assert!(e.points[1].x >= o.highlight_radius);
    assert!((e.points[0].x - label.cds.x).abs() < 1.0);
    let mut circles = o.clone();
    circles.atom_highlights_are_circles = true;
    let hl = DrawMolHighlights {
        atoms: vec![1],
        ..Default::default()
    };
    let (dm, _) = build_hl("CO", &circles, hl);
    assert_eq!(dm.highlights[0].points[1].x, dm.highlights[0].points[1].y);
}

#[test]
fn regions_and_attachments() {
    use rdoxide::graphmol::mol_draw2d::draw_shape::ShapeKind;
    let mut o = MolDrawOptions::default();
    o.atom_regions = vec![vec![0, 1, 2], vec![0]];
    let (dm, _) = build_hl("CCC", &o, DrawMolHighlights::default());
    // single-atom regions are ignored
    assert_eq!(dm.highlights.len(), 1);
    assert!(dm.highlights[0].fill && dm.highlights[0].points.len() == 4);

    let mut o = MolDrawOptions::default();
    o.dummies_are_attachments = true;
    let (dm, _) = build_hl("*CC", &o, DrawMolHighlights::default());
    let wavy = dm
        .bonds
        .iter()
        .filter(|b| matches!(b.kind, ShapeKind::WavyLine { .. }))
        .count();
    assert_eq!(wavy, 1);
    // the dummy has no label then
    assert!(dm.atom_labels[0].is_none());
}

#[test]
fn atom_and_bond_notes_avoid_the_drawing() {
    let mut o = MolDrawOptions::default();
    o.add_atom_indices = true;
    o.add_bond_indices = true;
    let (dm, _) = build_hl("CC(=O)N", &o, DrawMolHighlights::default());
    // one note per atom and per bond
    assert_eq!(dm.annotations.len(), 4 + 3);
    assert!(dm.annotations.iter().all(|a| a.font_scale > 0.0));
    // an atom note ends up away from the atom itself
    let first = &dm.annotations[0];
    assert!((first.pos - dm.at_cds[0]).length() > 0.2);
    // classes
    assert!(dm.annotations.iter().all(|a| a.class == "note"));
}

#[test]
fn cip_codes_and_stereo_groups() {
    use rdoxide::rdgeneral::RDProps;
    let mut o = MolDrawOptions::default();
    o.add_stereo_annotation = true;
    let mut t = text();
    let mut m = smi("C[C@H](N)C=CC");
    m.atom_mut(1).set_prop("_CIPCode", "S", true).unwrap();
    m.bond_mut(3)
        .set_stereo(rdoxide::graphmol::bond::BondStereo::StereoE)
        .unwrap();
    let mut o2 = o.clone();
    o2.prepare_mols_before_drawing = false;
    let mut dm = DrawMol::new(
        &m,
        "",
        300,
        300,
        &o2,
        &mut t,
        DrawMolHighlights::default(),
        true,
        -1,
        false,
    )
    .unwrap();
    dm.extract_all(1.0, &mut t).unwrap();
    let texts: Vec<_> = dm.annotations.iter().map(|a| a.text.as_str()).collect();
    assert!(
        texts.contains(&"(S)") && texts.contains(&"(E)"),
        "{texts:?}"
    );
    assert!(dm.annotations.iter().any(|a| a.class == "CIP_Code"));

    // enhanced stereo labels, and CIP codes hidden for OR/AND groups
    let mut m = smi("C[C@H](N)C[C@@H](O)C |o1:1,4|");
    m.atom_mut(1).set_prop("_CIPCode", "R", true).unwrap();
    let mut dm = DrawMol::new(
        &m,
        "",
        300,
        300,
        &o2,
        &mut t,
        DrawMolHighlights::default(),
        true,
        -1,
        false,
    )
    .unwrap();
    dm.extract_all(1.0, &mut t).unwrap();
    let texts: Vec<_> = dm.annotations.iter().map(|a| a.text.as_str()).collect();
    assert_eq!(
        texts.iter().filter(|t| **t == "or1").count(),
        2,
        "{texts:?}"
    );
    assert!(!texts.contains(&"(R)"));
    let mut shown = o2.clone();
    shown.show_all_cip_codes = true;
    let mut dm = DrawMol::new(
        &m,
        "",
        300,
        300,
        &shown,
        &mut t,
        DrawMolHighlights::default(),
        true,
        -1,
        false,
    )
    .unwrap();
    dm.extract_all(1.0, &mut t).unwrap();
    assert!(dm.annotations.iter().any(|a| a.text == "(R)"));
    // and none when the stereo group annotation is off
    let mut off = o2.clone();
    off.add_stereo_group_annotation = false;
    let mut dm = DrawMol::new(
        &m,
        "",
        300,
        300,
        &off,
        &mut t,
        DrawMolHighlights::default(),
        true,
        -1,
        false,
    )
    .unwrap();
    dm.extract_all(1.0, &mut t).unwrap();
    assert!(dm.annotations.iter().all(|a| a.class != "stereoGroup"));
}

#[test]
fn molecule_notes_and_radicals() {
    use rdoxide::rdgeneral::RDProps;
    let mut o = MolDrawOptions::default();
    o.prepare_mols_before_drawing = true;
    let mut t = text();
    let mut m = smi("CCO");
    m.set_prop("molNote", "my note", false).unwrap();
    let mut dm = DrawMol::new(
        &m,
        "",
        300,
        300,
        &o,
        &mut t,
        DrawMolHighlights::default(),
        true,
        -1,
        false,
    )
    .unwrap();
    dm.extract_all(1.0, &mut t).unwrap();
    dm.calculate_scale(&mut t);
    dm.finish_create_draw_objects(&mut t).unwrap();
    assert_eq!(dm.legends.len(), 1);
    assert_eq!(dm.legends[0].text, "my note");
    // in a corner of the canvas
    let p = dm.legends[0].pos;
    assert!(p.x >= 0.0 && p.x < 300.0 && p.y >= 0.0 && p.y <= 300.0);

    // the chiral flag gives ABS when asked for
    let mut m = smi("CCO");
    m.set_prop("_MolFileChiralFlag", 1, false).unwrap();
    let mut o2 = o.clone();
    o2.include_chiral_flag_label = true;
    let mut dm = DrawMol::new(
        &m,
        "",
        300,
        300,
        &o2,
        &mut t,
        DrawMolHighlights::default(),
        true,
        -1,
        false,
    )
    .unwrap();
    dm.extract_all(1.0, &mut t).unwrap();
    dm.calculate_scale(&mut t);
    dm.finish_create_draw_objects(&mut t).unwrap();
    assert_eq!(dm.legends[0].text, "ABS");

    // radicals get a rectangle and an orientation each
    let (dm, _) = build_hl("[CH2]C[O]", &o, DrawMolHighlights::default());
    assert_eq!(dm.radicals.len(), 2);
    let mut off = o.clone();
    off.include_radicals = false;
    let (dm, _) = build_hl("[CH2]C[O]", &off, DrawMolHighlights::default());
    assert!(dm.radicals.is_empty());
}

fn from_block(block: &str) -> RWMol {
    use rdoxide::graphmol::file_parsers::{MolFileParserParams, mol_from_mol_block};
    mol_from_mol_block(block, &MolFileParserParams::default())
        .unwrap()
        .unwrap()
}

fn run(m: &RWMol, o: &MolDrawOptions) -> (DrawMol, DrawTextFT) {
    let mut t = text();
    let mut dm = DrawMol::new(
        m,
        "",
        300,
        300,
        o,
        &mut t,
        DrawMolHighlights::default(),
        true,
        -1,
        false,
    )
    .unwrap();
    dm.extract_all(1.0, &mut t).unwrap();
    dm.calculate_scale(&mut t);
    dm.finish_create_draw_objects(&mut t).unwrap();
    (dm, t)
}

const SRU_BLOCK: &str = "
  test

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.5000    0.0000    0.0000 C   0  0
    3.0000    0.0000    0.0000 C   0  0
    4.5000    0.0000    0.0000 C   0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
M  STY  1   1 SRU
M  SAL   1  2   2   3
M  SBL   1  2   1   3
M  SDI   1  4    1.0000   -0.5000    1.0000    0.5000
M  SDI   1  4    3.5000   -0.5000    3.5000    0.5000
M  SCN  1   1 HT
M  SLB  1   1   1
M  SMT   1 n
M  END
";

#[test]
fn sgroup_brackets_and_labels() {
    let mut o = MolDrawOptions::default();
    o.prepare_mols_before_drawing = false;
    let (dm, _) = run(&from_block(SRU_BLOCK), &o);
    // two brackets of 4 points each
    let brackets: Vec<_> = dm
        .post_shapes
        .iter()
        .filter(|s| s.points.len() == 4)
        .collect();
    assert_eq!(brackets.len(), 2);
    let texts: Vec<_> = dm.annotations.iter().map(|a| a.text.as_str()).collect();
    assert!(texts.contains(&"n"), "{texts:?}");
    assert!(texts.contains(&"HT"), "{texts:?}");
    // the brackets are vertical lines to either side of the middle bonds
    let xs: Vec<f64> = brackets.iter().map(|b| b.points[1].x).collect();
    assert!(xs[0] < xs[1]);
    // no annotations when they aren't wanted
    let mut t = text();
    let mut dm = DrawMol::new(
        &from_block(SRU_BLOCK),
        "",
        300,
        300,
        &o,
        &mut t,
        DrawMolHighlights::default(),
        false,
        -1,
        false,
    )
    .unwrap();
    dm.extract_all(1.0, &mut t).unwrap();
    assert!(dm.annotations.is_empty());
    assert_eq!(dm.post_shapes.len(), 2);
}

#[test]
fn data_sgroup_labels() {
    use rdoxide::graphmol::mol_draw2d::draw_mol_sgroups::get_sgroup_data_labels;
    let block = "
  test

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.5000    0.0000    0.0000 O   0  0
  1  2  1  0
M  STY  1   1 DAT
M  SAL   1  1   2
M  SDT   1 FIELD                         F
M  SDD   1     0.0000    0.0000    DR    ALL  0       0
M  SED   1 some data
M  END
";
    let m = from_block(block);
    let labels = get_sgroup_data_labels(&m, 0.0).unwrap();
    assert_eq!(labels.len(), 1);
    assert_eq!(labels[0].text, "some data");
    assert_eq!(labels[0].atom_idx, 1);
    // no coordinates given relative to the atom: placed by the drawing code
    assert!(!labels[0].positioned);
    assert!(get_sgroup_data_labels(&smi("CC"), 0.0).unwrap().is_empty());
    let mut o = MolDrawOptions::default();
    o.prepare_mols_before_drawing = false;
    let (dm, _) = run(&m, &o);
    assert!(dm.annotations.iter().any(|a| a.text == "some data"));
}

#[test]
fn link_nodes() {
    use rdoxide::rdgeneral::RDProps;
    let mut m = smi("CCC");
    m.set_prop(
        rdoxide::rdgeneral::common_properties::molFileLinkNodes,
        "1 3 2 2 1 2 3",
        false,
    )
    .unwrap();
    let mut o = MolDrawOptions::default();
    o.prepare_mols_before_drawing = true;
    let (dm, _) = run(&m, &o);
    // a bracket across each of the two bonds, and the repeat label
    assert_eq!(dm.post_shapes.len(), 2);
    assert!(dm.annotations.iter().any(|a| a.text == "(1-3)"));
}

#[test]
fn variable_attachment_and_close_contacts() {
    let block = "
  test

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.5000    0.0000    0.0000 C   0  0
    1.5300    0.0000    0.0000 C   0  0
    3.0000    0.0000    0.0000 C   0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
M  END
";
    let mut o = MolDrawOptions::default();
    o.prepare_mols_before_drawing = false;
    let (dm, _) = run(&from_block(block), &o);
    // atoms 1 and 2 nearly coincide: a red box marks the first of them (the
    // other is skipped once flagged, as upstream)
    let boxes: Vec<_> = dm
        .post_shapes
        .iter()
        .filter(|s| s.points.len() == 5)
        .collect();
    assert_eq!(boxes.len(), 1);
    assert_eq!(boxes[0].atom1, 1);
    assert_eq!(
        boxes[0].line_colour,
        rdoxide::graphmol::mol_draw2d::helpers::DrawColour::new(1.0, 0.0, 0.0)
    );
    let mut off = o.clone();
    off.flag_close_contacts_dist = -1;
    let (dm, _) = run(&from_block(block), &off);
    assert!(dm.post_shapes.is_empty());

    // variable attachment points: a bond to a dummy with endpoints 2 and 3
    use rdoxide::rdgeneral::RDProps;
    let mut m = smi("*C(C)C");
    m.bond_mut(0)
        .set_prop("_MolFileBondEndPts", "(2 2 3)", false)
        .unwrap();
    m.bond_mut(0)
        .set_prop("_MolFileBondAttach", "ANY", false)
        .unwrap();
    let mut o = MolDrawOptions::default();
    o.prepare_mols_before_drawing = false;
    let mut t = text();
    let mut dm = DrawMol::new(
        &m,
        "",
        300,
        300,
        &o,
        &mut t,
        DrawMolHighlights::default(),
        true,
        -1,
        false,
    )
    .unwrap();
    dm.extract_all(1.0, &mut t).unwrap();
    // an ellipse for each end point and a fat line for the bond between them
    assert_eq!(dm.pre_shapes.len(), 3);
    // the dummy's * label is removed
    assert!(dm.atom_labels[0].is_none());
    // bad indices are rejected
    m.bond_mut(0)
        .set_prop("_MolFileBondEndPts", "(1 9)", false)
        .unwrap();
    let mut dm = DrawMol::new(
        &m,
        "",
        300,
        300,
        &o,
        &mut t,
        DrawMolHighlights::default(),
        true,
        -1,
        false,
    )
    .unwrap();
    assert!(dm.extract_all(1.0, &mut t).is_err());
}

#[test]
fn clashing_labels_are_reoriented() {
    // two unbonded atoms whose labels extend towards each other: N (E) and
    // O (W)
    let block = "
  test

  2  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0
    0.9000    0.5000    0.0000 O   0  0
M  END
";
    let mut o = MolDrawOptions::default();
    o.prepare_mols_before_drawing = false;
    let (dm, _) = run(&from_block(block), &o);
    let orig = [dm.atom_syms[0].1, dm.atom_syms[1].1];
    let now = [
        dm.atom_labels[0].as_ref().unwrap().orient,
        dm.atom_labels[1].as_ref().unwrap().orient,
    ];
    assert_ne!(orig, now);
}
