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
