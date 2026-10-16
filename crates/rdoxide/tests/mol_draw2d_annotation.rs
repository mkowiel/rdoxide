//! Tests of annotations and atom symbols.

use rdoxide::geometry::Point2D;
use rdoxide::graphmol::mol_draw2d::draw_annotation::*;
use rdoxide::graphmol::mol_draw2d::draw_text::DrawTextFT;
use rdoxide::graphmol::mol_draw2d::helpers::{DrawColour, OrientType, TextAlignType};

fn text() -> DrawTextFT {
    let mut t = DrawTextFT::new(40.0, 6.0, "").unwrap();
    t.set_font_size(12.0);
    t
}

#[test]
fn annotation_layout_scale_and_move() {
    let mut t = text();
    let col = DrawColour::new(0.0, 0.0, 0.0);
    let mut a = DrawAnnotation::new(
        "1.5",
        TextAlignType::Middle,
        "note",
        0.5,
        Point2D::new(10.0, 20.0),
        col,
        &mut t,
    )
    .unwrap();
    // the relative font scale was applied and the drawer's scale restored
    assert!((a.font_scale - 0.5 * t.font_scale()).abs() < 1e-12);
    let (w, h) = a.dimensions();
    assert!(w > 0.0 && h > 0.0);
    let mut ext = EMPTY_EXTREMES;
    a.find_extremes(&mut ext, 0.0);
    assert!(ext.0 < 10.0 && ext.1 > 10.0);

    a.translate(Point2D::new(5.0, 0.0));
    assert_eq!(a.pos, Point2D::new(15.0, 20.0));
    let (wm, hm) = a.dimensions();
    assert!((wm - w).abs() < 1e-9 && (hm - h).abs() < 1e-9);

    a.scale(Point2D::new(2.0, 2.0), &mut t).unwrap();
    let (w2, h2) = a.dimensions();
    assert!((w2 / w - 2.0).abs() < 1e-9 && (h2 / h - 2.0).abs() < 1e-9);
    assert_eq!(a.pos, Point2D::new(30.0, 40.0));

    let svg = a.draw_svg(&mut t, "bond-1").unwrap();
    assert!(svg.starts_with("<path class='bond-1 note' d='"), "{svg}");
}

#[test]
fn annotation_clash() {
    let mut t = text();
    let col = DrawColour::new(0.0, 0.0, 0.0);
    let a = DrawAnnotation::new(
        "X",
        TextAlignType::Middle,
        "",
        1.0,
        Point2D::new(0.0, 0.0),
        col,
        &mut t,
    )
    .unwrap();
    let probe = a.rects[0];
    assert!(a.does_rect_clash(&probe, 0.0));
    let mut far = probe;
    far.trans += Point2D::new(100.0, 100.0);
    assert!(!a.does_rect_clash(&far, 0.0));
}

#[test]
fn atom_symbol() {
    let t = text();
    let col = DrawColour::new(1.0, 0.0, 0.0);
    let mut s = AtomSymbol::new("NH", 3, OrientType::E, Point2D::new(50.0, 50.0), col, &t).unwrap();
    assert_eq!(s.rects().len(), 2);
    let mut ext = EMPTY_EXTREMES;
    s.find_extremes(&mut ext);
    assert!(ext.0 < 50.0 && ext.1 > 50.0);
    s.translate(Point2D::new(1.0, 1.0));
    assert_eq!(s.cds, Point2D::new(51.0, 51.0));
    let mut t2 = text();
    let svg = s.draw_svg(&mut t2, "").unwrap();
    assert!(svg.starts_with("<path class='atom-3' d='"));
    assert!(svg.contains("fill='#FF0000'"));
    // an empty symbol has no rects
    assert!(
        AtomSymbol::new("", 0, OrientType::C, Point2D::new(0.0, 0.0), col, &t)
            .unwrap()
            .rects()
            .is_empty()
    );
}

#[test]
fn colon_is_no_taller_than_its_neighbours() {
    let t = text();
    let col = DrawColour::new(0.0, 0.0, 0.0);
    let s = AtomSymbol::new("C:8", 0, OrientType::E, Point2D::new(0.0, 0.0), col, &t).unwrap();
    let r = s.rects();
    assert_eq!(r.len(), 3);
    assert_eq!(r[1].height, r[0].height.min(r[2].height));
    // markup is skipped when finding the colon
    let lit = AtomSymbol::new(
        "<lit>[C;X2:4]</lit>",
        0,
        OrientType::E,
        Point2D::new(0.0, 0.0),
        col,
        &t,
    )
    .unwrap();
    assert!(!lit.rects().is_empty());
}
