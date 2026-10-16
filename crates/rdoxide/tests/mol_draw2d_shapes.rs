//! Tests of the drawing shapes against a recording canvas.

use rdoxide::Result;
use rdoxide::geometry::Point2D;
use rdoxide::graphmol::mol_draw2d::draw_shape::*;
use rdoxide::graphmol::mol_draw2d::helpers::{DashPattern, DrawColour};
use rdoxide::graphmol::mol_draw2d::string_rect::StringRect;

#[derive(Default)]
struct Recorder {
    colour: DrawColour,
    width: f64,
    fill: bool,
    scale_bw: bool,
    split: bool,
    dash: DashPattern,
    log: Vec<String>,
    dashes_seen: Vec<DashPattern>,
}

impl ShapeCanvas for Recorder {
    fn colour(&self) -> DrawColour {
        self.colour
    }
    fn set_colour(&mut self, c: DrawColour) {
        self.colour = c;
    }
    fn line_width(&self) -> f64 {
        self.width
    }
    fn set_line_width(&mut self, w: f64) {
        self.width = w;
    }
    fn fill_polys(&self) -> bool {
        self.fill
    }
    fn set_fill_polys(&mut self, f: bool) {
        self.fill = f;
    }
    fn scale_bond_width(&self) -> bool {
        self.scale_bw
    }
    fn set_scale_bond_width(&mut self, s: bool) {
        self.scale_bw = s;
    }
    fn split_bonds(&self) -> bool {
        self.split
    }
    fn dash(&self) -> DashPattern {
        self.dash.clone()
    }
    fn set_dash(&mut self, d: DashPattern) {
        self.dash = d;
    }
    fn set_active_atm_idx(&mut self, _: i32, _: i32) {}
    fn set_active_bnd_idx(&mut self, _: i32) {}
    fn draw_arrow(
        &mut self,
        _: Point2D,
        _: Point2D,
        _: bool,
        _: f64,
        _: f64,
        _: DrawColour,
        _: bool,
    ) -> Result<()> {
        self.log.push("arrow".into());
        Ok(())
    }
    fn draw_ellipse(&mut self, a: Point2D, b: Point2D, _: bool) -> Result<()> {
        self.log
            .push(format!("ellipse {} {} {} {}", a.x, a.y, b.x, b.y));
        Ok(())
    }
    fn draw_line(
        &mut self,
        a: Point2D,
        b: Point2D,
        c1: DrawColour,
        _: DrawColour,
        _: bool,
    ) -> Result<()> {
        self.dashes_seen.push(self.dash.clone());
        self.log
            .push(format!("line {} {} {} {} r{}", a.x, a.y, b.x, b.y, c1.r));
        Ok(())
    }
    fn draw_polygon(&mut self, p: &[Point2D], _: bool) -> Result<()> {
        self.log.push(format!("polygon {}", p.len()));
        Ok(())
    }
    fn draw_triangle(&mut self, _: Point2D, _: Point2D, _: Point2D, _: bool) -> Result<()> {
        self.log.push("triangle".into());
        Ok(())
    }
    fn draw_wavy_line(
        &mut self,
        _: Point2D,
        _: Point2D,
        _: DrawColour,
        _: DrawColour,
        n: u32,
        _: f64,
        _: bool,
    ) -> Result<()> {
        self.log.push(format!("wavy {n}"));
        Ok(())
    }
    fn draw_arc(&mut self, _: Point2D, xr: f64, yr: f64, a1: f64, a2: f64, _: bool) -> Result<()> {
        self.log.push(format!("arc {xr} {yr} {a1} {a2}"));
        Ok(())
    }
}

fn pt(x: f64, y: f64) -> Point2D {
    Point2D::new(x, y)
}

const BLACK: DrawColour = DrawColour::new(0.0, 0.0, 0.0);
const RED: DrawColour = DrawColour::new(1.0, 0.0, 0.0);

#[test]
fn state_is_restored_and_ellipse_extremes() {
    let e = DrawShape::ellipse(vec![pt(5.0, 5.0), pt(2.0, 1.0)], 2.0, false, RED, true, 0).unwrap();
    let mut c = Recorder {
        width: 7.0,
        ..Default::default()
    };
    e.draw(&mut c).unwrap();
    assert_eq!(c.log, vec!["ellipse 3 4 7 6"]);
    assert_eq!((c.width, c.fill, c.colour), (7.0, false, BLACK));
    let mut ext = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    e.find_extremes(&mut ext);
    assert_eq!(ext, (3.0, 7.0, 4.0, 6.0));
    assert!(DrawShape::ellipse(vec![pt(0.0, 0.0)], 2.0, false, RED, true, 0).is_err());
}

#[test]
fn scale_and_move() {
    let mut l = DrawShape::simple_line(
        vec![pt(1.0, 1.0), pt(2.0, 3.0)],
        2.0,
        false,
        BLACK,
        (0, 1, 0),
        no_dash(),
    )
    .unwrap();
    l.scale(pt(2.0, 3.0));
    l.translate(pt(1.0, 1.0));
    assert_eq!(l.points, vec![pt(3.0, 4.0), pt(5.0, 10.0)]);
    let mut c = Recorder::default();
    l.draw(&mut c).unwrap();
    assert_eq!(c.log, vec!["line 3 4 5 10 r0"]);
}

#[test]
fn short_dashes_shrink_on_short_lines() {
    let draw = |len: f64| {
        let line = DrawShape::simple_line(
            vec![pt(0.0, 0.0), pt(len, 0.0)],
            2.0,
            false,
            BLACK,
            (0, 1, 0),
            short_dashes(),
        )
        .unwrap();
        let mut rec = Recorder::default();
        line.draw(&mut rec).unwrap();
        // restored to the original after drawing
        assert!(rec.dash.is_empty());
        rec.dashes_seen[0].clone()
    };
    assert_eq!(draw(5.0), vec![0.5, 2.0 / 3.0]);
    assert_eq!(draw(20.0), vec![1.0, 2.0 / 1.5]);
    assert_eq!(draw(40.0), short_dashes());
}

#[test]
fn dashed_wedge_lines() {
    let w = DrawShape::dashed_wedge(
        vec![pt(0.0, 0.0), pt(40.0, -5.0), pt(40.0, 5.0)],
        BLACK,
        RED,
        true,
        1.0,
        (0, 1, 0),
    )
    .unwrap();
    // one less dash: n_dashes - 1 lines, each a pair of points
    assert_eq!(w.points.len() % 2, 0);
    assert!(w.points.len() >= 6);
    let ShapeKind::DashedWedge { line_colours, .. } = &w.kind else {
        panic!()
    };
    assert_eq!(line_colours.len(), w.points.len() / 2);
    assert_eq!(line_colours[0], BLACK);
    assert_eq!(*line_colours.last().unwrap(), RED);
    let mut c = Recorder::default();
    w.draw(&mut c).unwrap();
    assert_eq!(c.log.len(), line_colours.len());
}

#[test]
fn solid_wedge_triangles() {
    let pts = vec![pt(0.0, 0.0), pt(10.0, -2.0), pt(10.0, 2.0)];
    let one =
        DrawShape::solid_wedge(pts.clone(), BLACK, BLACK, false, vec![], 2.0, (0, 1, 0)).unwrap();
    assert_eq!(one.points.len(), 3);
    let two =
        DrawShape::solid_wedge(pts.clone(), BLACK, RED, false, vec![], 2.0, (0, 1, 0)).unwrap();
    assert_eq!(two.points.len(), 9);
    let mut c = Recorder::default();
    two.draw(&mut c).unwrap();
    assert_eq!(c.log, vec!["triangle", "polygon 4"]);
    // three other bonds are trimmed to the two with the widest angle
    let y = DrawShape::solid_wedge(
        pts,
        BLACK,
        BLACK,
        false,
        vec![pt(1.0, 0.1), pt(0.0, 1.0), pt(0.0, -1.0)],
        2.0,
        (0, 1, 0),
    )
    .unwrap();
    assert_eq!(y.points.len(), 6);
}

#[test]
fn clashes_with_rectangles() {
    let rect = StringRect::new(pt(0.0, 0.0), pt(0.0, 0.0), 2.0, 2.0);
    let line = DrawShape::simple_line(
        vec![pt(-5.0, 0.0), pt(5.0, 0.0)],
        2.0,
        false,
        BLACK,
        (0, 1, 0),
        no_dash(),
    )
    .unwrap();
    assert!(line.does_rect_clash(&rect, 0.0));
    let far = DrawShape::simple_line(
        vec![pt(-5.0, 9.0), pt(5.0, 9.0)],
        2.0,
        false,
        BLACK,
        (0, 1, 0),
        no_dash(),
    )
    .unwrap();
    assert!(!far.does_rect_clash(&rect, 0.0));
    let arc = DrawShape::arc(
        vec![pt(0.0, 0.0), pt(10.0, 10.0)],
        0.0,
        90.0,
        2.0,
        false,
        BLACK,
        false,
        0,
    )
    .unwrap();
    assert!(arc.does_rect_clash(&rect, 0.0));
    assert!(
        DrawShape::arc(
            vec![pt(0.0, 0.0), pt(1.0, 1.0)],
            0.0,
            400.0,
            2.0,
            false,
            BLACK,
            false,
            0
        )
        .is_err()
    );
}

#[test]
fn arrow_and_wavy() {
    let (frac, angle) = DrawShape::default_arrow_head();
    let a = DrawShape::arrow(
        vec![pt(0.0, 0.0), pt(10.0, 0.0)],
        2.0,
        false,
        BLACK,
        false,
        (0, 1, 0),
        frac,
        angle,
    )
    .unwrap();
    let mut ext = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    a.find_extremes(&mut ext);
    // upstream passes the arrow tip as both ends of calcArrowHead, so the
    // stored barb ends coincide with the tip
    assert_eq!(ext, (0.0, 10.0, 0.0, 0.0));
    let w = DrawShape::wavy_line(
        vec![pt(0.0, 0.0), pt(10.0, 0.0)],
        2.0,
        false,
        BLACK,
        BLACK,
        0.5,
        (0, 1, 0),
    )
    .unwrap();
    let mut c = Recorder::default();
    w.draw(&mut c).unwrap();
    assert_eq!(c.log, vec!["wavy 30"]);
}
