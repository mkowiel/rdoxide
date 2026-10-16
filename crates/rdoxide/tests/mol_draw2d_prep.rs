//! Tests of `MolDraw2DUtils` and the DrawMol helpers.
#![allow(clippy::field_reassign_with_default)]

use rdoxide::geometry::Point2D;
use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::bond::{BondDir, BondStereo};
use rdoxide::graphmol::mol_draw2d::draw_mol_utils::*;
use rdoxide::graphmol::mol_draw2d::helpers::{DrawColour, MolDrawOptions, OrientType};
use rdoxide::graphmol::mol_draw2d::utils::*;
use rdoxide::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use rdoxide::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

fn smi(s: &str) -> RWMol {
    mol_from_smiles(s, &SmilesParserParams::default())
        .unwrap()
        .unwrap()
}

/// Draw-frame coordinates (y inverted) for a prepared molecule.
fn draw_cds(m: &RWMol) -> Vec<Point2D> {
    let conf = m.conformer(-1).unwrap();
    (0..m.num_atoms())
        .map(|i| {
            let p = conf.atom_pos(i).unwrap();
            Point2D::new(p.x, -p.y)
        })
        .collect()
}

#[test]
fn prepare_generates_coords_and_wedges() {
    let mut m = smi("C[C@H](N)C(=O)O");
    assert_eq!(m.num_conformers(), 0);
    prepare_mol_for_drawing(&mut m, true, true, true, false, false).unwrap();
    assert_eq!(m.num_conformers(), 1);
    assert!(
        m.bonds()
            .iter()
            .any(|b| b.bond_dir() == BondDir::BeginWedge || b.bond_dir() == BondDir::BeginDash)
    );
    let mbl = mean_bond_length(&m, -1).unwrap();
    assert!(mbl > 0.5 && mbl < 3.0);
    // a molecule without bonds has mean length 0
    assert_eq!(mean_bond_length(&smi("[Na+]"), -1).unwrap(), 0.0);
}

#[test]
fn prepare_adds_wavy_bonds() {
    let mut m = smi("CC=CC");
    m.bond_mut(1).set_stereo(BondStereo::StereoAny).unwrap();
    prepare_mol_for_drawing(&mut m, true, false, true, false, true).unwrap();
    assert!(m.bonds().iter().any(|b| b.bond_dir() == BondDir::Unknown));
    assert_eq!(m.bond(1).stereo(), BondStereo::StereoNone);
}

#[test]
fn acs1996() {
    let mut o = MolDrawOptions::default();
    assert!(set_acs1996_options(&mut o, 0.0).is_err());
    set_acs1996_options(&mut o, 1.5).unwrap();
    assert!((o.fixed_bond_length - 9.6).abs() < 1e-12);
    assert_eq!(o.fixed_font_size, 10);
    assert!(o.font_file.ends_with("FreeSans.ttf") || o.font_file == "BuiltinRobotoRegular");
}

#[test]
fn atom_symbols_and_orientation() {
    let mut m = smi("C[NH3+].[2H]OC(=O)[13CH3].C[O-]");
    prepare_mol_for_drawing(&mut m, true, false, false, false, false).unwrap();
    let cds = draw_cds(&m);
    let opts = MolDrawOptions::default();
    let sym = |i: u32| {
        let o = get_atom_orientation(&m, i, &cds);
        get_atom_symbol(&m, i, &cds, o, &opts).unwrap()
    };
    // plain carbons have no label
    assert_eq!(sym(0), "");
    assert_eq!(sym(1), "N<sup>+</sup>H<sub>3</sub>");
    assert_eq!(sym(3), "O");
    assert_eq!(sym(2), "<sup>2</sup>H");
    assert_eq!(sym(6), "<sup>13</sup>C");
    assert_eq!(sym(8), "O<sup>-</sup>");
    let mut o2 = MolDrawOptions::default();
    o2.no_atom_labels = true;
    assert_eq!(
        get_atom_symbol(&m, 1, &cds, OrientType::E, &o2).unwrap(),
        ""
    );
    o2.no_atom_labels = false;
    o2.atom_labels.insert(1, "Foo".into());
    assert_eq!(
        get_atom_symbol(&m, 1, &cds, OrientType::E, &o2).unwrap(),
        "<lit>Foo</lit>"
    );
    o2.atom_labels.clear();
    o2.atom_label_deuterium_tritium = true;
    assert_eq!(
        get_atom_symbol(&m, 2, &cds, OrientType::E, &o2).unwrap(),
        "<lit>D</lit>"
    );
}

#[test]
fn orientation_of_terminal_atoms() {
    // a lone atom: halogen-like elements list the H first
    let m = smi("[Cl]");
    assert_eq!(
        get_atom_orientation(&m, 0, &[Point2D::new(0.0, 0.0)]),
        OrientType::W
    );
    let m = smi("[Na]");
    assert_eq!(
        get_atom_orientation(&m, 0, &[Point2D::new(0.0, 0.0)]),
        OrientType::E
    );
    // a bond going right makes the label of the other atom W, going left E
    let m = smi("CO");
    let cds = [Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)];
    assert_eq!(get_atom_orientation(&m, 0, &cds), OrientType::W);
    assert_eq!(get_atom_orientation(&m, 1, &cds), OrientType::E);
    // vertical bonds of terminal atoms give E rather than N/S
    let cds = [Point2D::new(0.0, 0.0), Point2D::new(0.0, 1.0)];
    assert_eq!(get_atom_orientation(&m, 1, &cds), OrientType::E);
}

#[test]
fn query_symbols() {
    let q = mol_from_smarts("[F,Cl,Br,I,At]C", &SmartsParserParams::default())
        .unwrap()
        .unwrap();
    let cds = [Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0)];
    let opts = MolDrawOptions::default();
    assert_eq!(
        get_atom_symbol(&q, 0, &cds, OrientType::E, &opts).unwrap(),
        "<lit>X</lit>"
    );
    let mut o = MolDrawOptions::default();
    o.use_complex_query_atom_symbols = false;
    assert_eq!(
        get_atom_symbol(&q, 0, &cds, OrientType::E, &o).unwrap(),
        "<lit>[F,Cl,Br,I,At]</lit>"
    );
    o.brackets_around_atom_lists = false;
    assert_eq!(
        get_atom_symbol(&q, 0, &cds, OrientType::E, &o).unwrap(),
        "<lit>F,Cl,Br,I,At</lit>"
    );
    assert_eq!(get_atom_list_text(q.atom(0)).unwrap(), "[F,Cl,Br,I,At]");
}

#[test]
fn colours_and_highlights() {
    let mut o = MolDrawOptions::default();
    assert_eq!(
        get_colour_by_atomic_num(7, &o),
        DrawColour::new(0.0, 0.0, 1.0)
    );
    // unknown elements use the -1 colour, hydrogens may use the 201 colour
    assert_eq!(
        get_colour_by_atomic_num(50, &o),
        DrawColour::new(0.0, 0.0, 0.0)
    );
    o.no_atom_labels = true;
    assert_eq!(
        get_colour_by_atomic_num(1, &o),
        DrawColour::new(0.68, 0.85, 0.90)
    );
    let o = MolDrawOptions::default();
    assert_eq!(get_highlight_bond_width(&o, 0, None), 2.0 * 8.0);
    let mut o2 = MolDrawOptions::default();
    o2.fill_highlights = false;
    assert_eq!(get_highlight_bond_width(&o2, 0, None), 2.0 * 4.0);
    let mult = std::collections::BTreeMap::from([(3, 5)]);
    assert_eq!(get_highlight_bond_width(&o, 3, Some(&mult)), 10.0);

    let m = smi("CCO");
    assert_eq!(get_bond_highlights_for_atoms(&m, &[0, 1, 2]), vec![0, 1]);
    let atom_map = std::collections::BTreeMap::from([(1, DrawColour::new(0.0, 1.0, 0.0))]);
    let col = get_highlight_bond_colour(&m, 0, &o, &[0], &Default::default(), &atom_map);
    assert_eq!(col, DrawColour::new(0.0, 1.0, 0.0));
    assert_eq!(
        get_highlight_bond_colour(&m, 1, &o, &[0], &Default::default(), &atom_map),
        DrawColour::new(0.0, 0.0, 0.0)
    );
}

#[test]
fn geometry() {
    let p = calc_perpendicular(Point2D::new(1.0, 0.0), Point2D::new(0.0, 0.0));
    assert_eq!((p.x, p.y), (-0.0, 1.0));
    let inner = calc_inner_perpendicular(
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 0.0),
        Point2D::new(1.0, 1.0),
    );
    assert!(inner.y > 0.0);
    let (a1, a2, a3) = (
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 0.0),
        Point2D::new(1.0, 0.0),
    );
    assert!(are_bonds_trans(a1, a2, a3, Point2D::new(2.0, 0.0)));
    assert!(!are_bonds_trans(a1, a2, a3, Point2D::new(1.0, 1.0)));
    assert!(
        are_bonds_parallel(
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(2.0, 0.0),
            Point2D::new(1.0, 0.0),
            1e-4
        )
        .unwrap()
    );
    let m = smi("CC(C)O");
    assert_eq!(other_neighbor(&m, 1, 0, 0), Some(2));
    assert_eq!(other_neighbor(&m, 1, 0, 1), Some(3));
    assert_eq!(other_neighbor(&m, 1, 0, 2), None);
    let mut p = Point2D::new(0.0, 0.0);
    adjust_bond_end_for_string(Point2D::new(10.0, 0.0), 0.0, &[], &mut p);
    assert_eq!((p.x, p.y), (0.0, 0.0));
}

#[test]
fn centering() {
    let mut m = smi("CCO");
    prepare_mol_for_drawing(&mut m, true, false, false, false, false).unwrap();
    center_mol_for_drawing(&mut m, -1).unwrap();
    let c = m.conformer(-1).unwrap();
    let n = m.num_atoms() as f64;
    let (mut x, mut y) = (0.0, 0.0);
    for i in 0..m.num_atoms() {
        let p = c.atom_pos(i).unwrap();
        x += p.x;
        y += p.y;
    }
    assert!((x / n).abs() < 1e-9 && (y / n).abs() < 1e-9);
    // the middle atom of an allene is linear
    let mut l = smi("CC=C=CC");
    prepare_mol_for_drawing(&mut l, true, false, false, false, false).unwrap();
    let cds = draw_cds(&l);
    assert!(is_linear_atom(&l, 2, &cds));
    assert!(!is_linear_atom(&l, 1, &cds));
}
