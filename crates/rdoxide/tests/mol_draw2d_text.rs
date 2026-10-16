//! Tests of the drawing foundation: colours, options, string rectangles,
//! the TrueType reader and text layout.

use rdoxide::geometry::Point2D;
use rdoxide::graphmol::mol_draw2d::details::*;
use rdoxide::graphmol::mol_draw2d::draw_text::*;
use rdoxide::graphmol::mol_draw2d::font::*;
use rdoxide::graphmol::mol_draw2d::helpers::*;
use rdoxide::graphmol::mol_draw2d::string_rect::StringRect;

#[test]
fn colours_and_options() {
    let a = DrawColour::new(0.5, 0.25, 1.0);
    let b = DrawColour::with_alpha(0.5004, 0.25, 1.0, 0.1);
    assert_ne!(a, b);
    assert!(a.feq(&b, 0.001, true));
    assert!(!a.feq(&b, 0.001, false));
    assert_eq!((a * 2.0).r, 1.0);
    assert!(a.div_checked(0.0).is_err());
    assert_eq!(
        colour_to_svg(&DrawColour::new(1.0, 0.5, 0.0)).unwrap(),
        "#FF7F00"
    );
    assert_eq!(
        colour_to_svg(&DrawColour::with_alpha(0.0, 0.0, 0.0, 0.5)).unwrap(),
        "#0000007F"
    );
    assert!(colour_to_svg(&DrawColour::new(1.5, 0.0, 0.0)).is_err());

    let mut o = MolDrawOptions::default();
    assert_eq!(o.highlight_colour_palette.len(), 10);
    assert_eq!(o.atom_colour_palette[&7], DrawColour::new(0.0, 0.0, 1.0));
    assert_eq!(o.drawing_extents_include, draw_element::ALL);
    o.use_cdk_atom_palette();
    assert_eq!(o.atom_colour_palette[&5], DrawColour::new(1.0, 0.71, 0.71));
    o.set_dark_mode();
    assert_eq!(o.background_colour, DrawColour::new(0.0, 0.0, 0.0));
    o.set_monochrome_mode(
        DrawColour::new(1.0, 0.0, 0.0),
        DrawColour::new(0.0, 0.0, 1.0),
    );
    assert_eq!(o.atom_colour_palette.len(), 1);
}

#[test]
fn string_rect_geometry() {
    let r = StringRect::new(Point2D::new(1.0, 1.0), Point2D::new(1.0, 1.0), 2.0, 2.0);
    assert!(r.is_point_inside(Point2D::new(0.5, 0.5), 0.0));
    assert!(!r.is_point_inside(Point2D::new(3.0, 0.5), 0.0));
    let mut other = r;
    other.trans = Point2D::new(1.5, 0.0);
    assert!(r.does_it_intersect(&other, 0.0));
    other.trans = Point2D::new(5.0, 0.0);
    assert!(!r.does_it_intersect(&other, 0.0));
}

#[test]
fn geometry_helpers() {
    let ip = do_lines_intersect(
        Point2D::new(0.0, 0.0),
        Point2D::new(2.0, 2.0),
        Point2D::new(0.0, 2.0),
        Point2D::new(2.0, 0.0),
    )
    .unwrap();
    assert!((ip.x - 1.0).abs() < 1e-12 && (ip.y - 1.0).abs() < 1e-12);
    assert!(
        do_lines_intersect(
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(0.0, 1.0),
            Point2D::new(1.0, 1.0),
        )
        .is_none()
    );
    assert!(is_point_in_triangle(
        Point2D::new(0.2, 0.2),
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 0.0),
        Point2D::new(0.0, 1.0)
    ));
    let mut end = Point2D::new(0.0, 0.0);
    adjust_line_end_for_ellipse(
        Point2D::new(0.0, 0.0),
        1.0,
        1.0,
        Point2D::new(5.0, 0.0),
        &mut end,
    );
    assert!((end.x - 1.0).abs() < 1e-12 && end.y.abs() < 1e-12);
    let segs =
        get_wavy_line_segments(Point2D::new(0.0, 0.0), Point2D::new(3.0, 0.0), 3, 0.1).unwrap();
    assert_eq!(segs.len(), 4);
    assert_eq!(
        arc_points(Point2D::new(0.0, 0.0), Point2D::new(2.0, 2.0), 0.0, 90.0).len(),
        181
    );
    assert_eq!(format_double(2.25), "2.2");
}

#[test]
fn builtin_fonts_load() {
    for font in [
        TrueTypeFont::builtin_telex_regular(),
        TrueTypeFont::builtin_roboto_regular(),
    ] {
        assert!(font.units_per_em() > 0);
        // the glyph box of a letter is inside the em box and has an advance
        for c in b'A'..=b'z' {
            let g = font.load_char(u32::from(c)).unwrap();
            if !c.is_ascii_alphabetic() {
                continue;
            }
            let bb = g.outline.bbox();
            assert!(bb.x_max > bb.x_min && bb.y_max > bb.y_min, "char {c}");
            assert!(g.advance > 0);
            // the decomposition closes every contour on its start point
            let cmds = g.outline.decompose();
            let moves = cmds
                .iter()
                .filter(|c| matches!(c, PathCommand::MoveTo(..)))
                .count();
            assert_eq!(moves, g.outline.contour_ends.len());
        }
        // space has an advance but no outline
        let sp = font.load_char(u32::from(b' ')).unwrap();
        assert!(sp.advance > 0 && sp.outline.points.is_empty());
        assert_eq!(sp.outline.bbox(), BBox::default());
    }
}

#[test]
fn bbox_checks_curve_extrema() {
    // a quadratic whose control point is far above the end points: the curve
    // peaks halfway between, at 50 (y2 + (y1-y2)(y3-y2)/(y1-y2+y3-y2) = 100-50)
    let outline = Outline {
        points: vec![
            OutlinePoint {
                x: 0,
                y: 0,
                on_curve: true,
            },
            OutlinePoint {
                x: 50,
                y: 100,
                on_curve: false,
            },
            OutlinePoint {
                x: 100,
                y: 0,
                on_curve: true,
            },
        ],
        contour_ends: vec![2],
    };
    let bb = outline.bbox();
    assert_eq!((bb.x_min, bb.y_min, bb.x_max, bb.y_max), (0, 0, 100, 50));
}

#[test]
fn string_layout() {
    let dt = DrawTextFT::new(40.0, 6.0, "").unwrap();
    // pieces of atom labels
    assert_eq!(
        atom_label_to_pieces("CH<sub>3</sub>", OrientType::E),
        vec!["C".to_string(), "H<sub>3</sub>".to_string()]
    );
    assert_eq!(
        atom_label_to_pieces("NH<sub>2</sub>O", OrientType::N),
        vec![
            "N".to_string(),
            "H<sub>2</sub>".to_string(),
            "O".to_string()
        ]
    );
    assert_eq!(
        atom_label_to_pieces("<lit>ab</lit>", OrientType::N),
        vec!["ab".to_string()]
    );

    let (w1, h1) = dt.string_size("C").unwrap();
    let (w2, h2) = dt.string_size("CC").unwrap();
    assert!(w1 > 0.0 && h1 > 0.0);
    assert!(w2 > 1.5 * w1 && (h2 - h1).abs() < 1e-12);

    // sub- and superscripts are smaller and shifted
    let l = dt
        .string_rects(
            "N<sup>+</sup>H<sub>2</sub>",
            OrientType::E,
            true,
            TextAlignType::Start,
        )
        .unwrap();
    assert_eq!(l.draw_chars, b"N+H2");
    assert_eq!(
        l.draw_modes,
        vec![
            TextDrawType::Normal,
            TextDrawType::Superscript,
            TextDrawType::Normal,
            TextDrawType::Subscript
        ]
    );
    assert!(l.rects[1].rect_corr > 0.0 && l.rects[3].rect_corr < 0.0);
    assert!(
        dt.string_rects("", OrientType::E, false, TextAlignType::Middle)
            .is_err()
    );
}

#[test]
fn svg_for_text() {
    let mut dt = DrawTextFT::new(40.0, 6.0, "BuiltinRobotoRegular").unwrap();
    dt.set_font_size(0.6 * 20.0);
    dt.set_colour(DrawColour::new(1.0, 0.0, 0.0));
    let svg = dt
        .draw_string_svg(
            "O",
            Point2D::new(50.0, 50.0),
            TextAlignType::Middle,
            "atom-0",
        )
        .unwrap();
    assert!(svg.starts_with("<path class='atom-0' d='M "), "{svg}");
    assert!(svg.ends_with("' fill='#FF0000'/>\n"), "{svg}");
    // an O is two contours
    assert_eq!(svg.matches("M ").count(), 2);
    assert!(svg.contains("Q "));
    assert!(DrawTextFT::new(40.0, 6.0, "/no/such/font.ttf").is_err());
}
