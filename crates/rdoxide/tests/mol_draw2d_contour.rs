//! Tests of contouring: CONREC and drawing contours.

use rdoxide::geometry::Point2D;
use rdoxide::graphmol::mol_draw2d::contour::{
    ContourParams, contour_and_draw_gaussians, contour_and_draw_grid,
};
use rdoxide::graphmol::mol_draw2d::drawer::MolDraw2DSVG;
use rdoxide::numerics::conrec::{connect_line_segments, contour};

/// A 3x3 grid with a peak of 1 in the middle.
fn peak() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut d = vec![0.0; 9];
    d[4] = 1.0;
    (d, vec![0.0, 1.0, 2.0], vec![0.0, 1.0, 2.0])
}

#[test]
fn conrec_contours_a_peak() {
    let (d, x, y) = peak();
    let segs = contour(&d, 0, 2, 0, 2, &x, &y, 1, &[0.5]).unwrap();
    assert!(!segs.is_empty());
    assert!(segs.iter().all(|s| s.iso_val == 0.5));
    // all points lie inside the grid, around the peak
    for s in &segs {
        for p in [s.p1, s.p2] {
            assert!((0.0..=2.0).contains(&p.x) && (0.0..=2.0).contains(&p.y));
            assert!((p.x - 1.0).abs() + (p.y - 1.0).abs() > 0.0);
        }
    }
    // the segments join up into one closed line
    let lines = connect_line_segments(&segs, 1000.0, 1e6);
    assert_eq!(lines.len(), 1);
    assert_eq!(lines[0].1, 0.5);
    assert_eq!(lines[0].0.len(), segs.len() + 1);
    assert_eq!(lines[0].0.first(), lines[0].0.last());
    // levels outside the data give nothing; bad input is an error
    assert!(
        contour(&d, 0, 2, 0, 2, &x, &y, 1, &[5.0])
            .unwrap()
            .is_empty()
    );
    assert!(contour(&d, 0, 0, 0, 2, &x, &y, 1, &[0.5]).is_err());
    assert!(contour(&d, 0, 2, 0, 2, &x, &y, 0, &[0.5]).is_err());
}

#[test]
fn contours_of_a_grid_are_drawn() {
    let (d, x, y) = peak();
    let mut drawer = MolDraw2DSVG::new(200, 200, -1, -1).unwrap();
    let mut levels = Vec::new();
    contour_and_draw_grid(
        &mut drawer,
        &d,
        &x,
        &y,
        3,
        &mut levels,
        &ContourParams::default(),
        None,
    )
    .unwrap();
    // levels were made from the data
    assert_eq!(levels, vec![0.0, 0.5, 1.0]);
    drawer.finish_drawing();
    let svg = drawer.get_drawing_text();
    assert!(svg.contains("<path"));
    assert!(svg.contains("stroke:#7F7F7F7F"));
    // too many contours for the levels given (the data are only scanned,
    // as upstream, when the levels are empty or the grid is filled)
    let mut few = vec![0.5];
    let fill = ContourParams {
        fill_grid: true,
        ..Default::default()
    };
    assert!(contour_and_draw_grid(&mut drawer, &d, &x, &y, 3, &mut few, &fill, None).is_err());
    // a flat grid draws nothing
    let flat = vec![1.0; 9];
    let mut lv = Vec::new();
    let before = drawer.get_drawing_text().len();
    contour_and_draw_grid(
        &mut drawer,
        &flat,
        &x,
        &y,
        3,
        &mut lv,
        &ContourParams::default(),
        None,
    )
    .unwrap();
    assert_eq!(drawer.get_drawing_text().len(), before);
}

#[test]
fn filled_grids_and_gaussians() {
    let mut drawer = MolDraw2DSVG::new(200, 200, -1, -1).unwrap();
    let params = ContourParams {
        fill_grid: true,
        dash_negative: true,
        ..Default::default()
    };
    let mut levels = Vec::new();
    contour_and_draw_gaussians(
        &mut drawer,
        &[Point2D::new(0.0, 0.0), Point2D::new(1.5, 0.5)],
        &[1.0, -1.0],
        &[0.4, 0.4],
        5,
        &mut levels,
        &params,
        None,
    )
    .unwrap();
    assert_eq!(levels.len(), 5);
    drawer.finish_drawing();
    let svg = drawer.get_drawing_text();
    // filled boxes and dashed negative contours
    assert!(svg.contains("fill-rule:evenodd;fill-opacity:0.5"));
    assert!(svg.contains("stroke-dasharray:2,2"));
    // the molecule coordinates are set by the contouring
    assert!(drawer.range().unwrap().x > 1.5);
    // mismatched sizes are an error
    assert!(
        contour_and_draw_gaussians(
            &mut drawer,
            &[Point2D::new(0.0, 0.0)],
            &[1.0, 2.0],
            &[0.4],
            5,
            &mut Vec::new(),
            &params,
            None
        )
        .is_err()
    );
}
