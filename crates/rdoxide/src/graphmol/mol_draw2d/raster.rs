//! PNG output of the molecule drawer. `MolDraw2DCairo` of RDKit draws with Cairo; here the SVG that the drawer writes
//! (`<rect>`, `<path>` and `<ellipse>` elements: lines, polygons, arcs, wavy lines and the glyph outlines of the text)
//! is rendered with `tiny-skia` (pure Rust, anti-aliased fills and strokes) and encoded as PNG. Geometry, colours and
//! text layout are the drawer's; only the anti-aliasing of edge pixels differs from Cairo.

use tiny_skia::{
    FillRule, LineCap, LineJoin, Paint, PathBuilder, Pixmap, Stroke, StrokeDash, Transform,
};

use std::ops::{Deref, DerefMut};

use crate::error::{RDKitError, Result};

use super::drawer::MolDraw2DSVG;
use super::helpers::MolDrawOptions;
use crate::graphmol::ROMol;
use crate::graphmol::file_parsers::mol_file_writer::{MolWriterParams, mol_to_v3k_mol_block};
use crate::graphmol::file_parsers::png::{
    MOL_TAG, PKL_TAG, PngMetadataParams, SMILES_TAG, add_metadata_to_png, augment_tag_name,
};
use crate::graphmol::pickle::pickle_mol;
use crate::graphmol::smiles_parse::cx_write::mol_to_cx_smiles;

/// The attributes of a tag, `name='value'` or `name="value"`.
fn attributes(tag: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let b = tag.as_bytes();
    let mut i = 0;
    // skip the tag name
    while i < b.len() && !b[i].is_ascii_whitespace() {
        i += 1;
    }
    while i < b.len() {
        while i < b.len() && (b[i].is_ascii_whitespace() || b[i] == b'/') {
            i += 1;
        }
        let start = i;
        while i < b.len() && b[i] != b'=' && !b[i].is_ascii_whitespace() && b[i] != b'/' {
            i += 1;
        }
        let name = &tag[start..i];
        while i < b.len() && b[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= b.len() || b[i] != b'=' {
            continue;
        }
        i += 1;
        while i < b.len() && b[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= b.len() || (b[i] != b'\'' && b[i] != b'"') {
            continue;
        }
        let quote = b[i];
        i += 1;
        let vstart = i;
        while i < b.len() && b[i] != quote {
            i += 1;
        }
        out.push((name.to_string(), tag[vstart..i.min(b.len())].to_string()));
        i += 1;
    }
    out
}

fn attr<'a>(attrs: &'a [(String, String)], name: &str) -> Option<&'a str> {
    attrs
        .iter()
        .find(|(n, _)| n == name)
        .map(|(_, v)| v.as_str())
}

fn number(s: &str) -> f32 {
    s.trim()
        .trim_end_matches("px")
        .parse::<f32>()
        .unwrap_or(0.0)
}

/// `#RRGGBB` to rgb
fn colour(s: &str) -> Option<(u8, u8, u8)> {
    let s = s.trim();
    let hex = s.strip_prefix('#')?;
    if hex.len() != 6 {
        return None;
    }
    let v = u32::from_str_radix(hex, 16).ok()?;
    Some(((v >> 16) as u8, (v >> 8) as u8, v as u8))
}

#[derive(Debug, Default, Clone)]
struct Style {
    fill: Option<(u8, u8, u8)>,
    fill_opacity: f32,
    fill_evenodd: bool,
    stroke: Option<(u8, u8, u8)>,
    stroke_opacity: f32,
    stroke_width: f32,
    cap: Option<LineCap>,
    join: Option<LineJoin>,
    miter_limit: f32,
    dash: Vec<f32>,
}

impl Style {
    fn from_attributes(attrs: &[(String, String)]) -> Style {
        let mut st = Style {
            fill_opacity: 1.0,
            stroke_opacity: 1.0,
            stroke_width: 1.0,
            miter_limit: 4.0,
            ..Style::default()
        };
        // glyph outlines carry a plain fill attribute
        if let Some(f) = attr(attrs, "fill") {
            st.fill = colour(f);
        } else if attr(attrs, "style").is_none() {
            st.fill = Some((0, 0, 0));
        }
        if let Some(style) = attr(attrs, "style") {
            let mut opacity = 1.0f32;
            for item in style.split(';') {
                let Some((k, v)) = item.split_once(':') else {
                    continue;
                };
                let v = v.trim();
                match k.trim() {
                    "fill" => st.fill = colour(v),
                    "fill-opacity" => st.fill_opacity = number(v),
                    "fill-rule" => st.fill_evenodd = v == "evenodd",
                    "stroke" => st.stroke = colour(v),
                    "stroke-opacity" => st.stroke_opacity = number(v),
                    "stroke-width" => st.stroke_width = number(v),
                    "stroke-linecap" => {
                        st.cap = Some(match v {
                            "round" => LineCap::Round,
                            "square" => LineCap::Square,
                            _ => LineCap::Butt,
                        })
                    }
                    "stroke-linejoin" => {
                        st.join = Some(match v {
                            "round" => LineJoin::Round,
                            "bevel" => LineJoin::Bevel,
                            _ => LineJoin::Miter,
                        })
                    }
                    "stroke-miterlimit" => st.miter_limit = number(v),
                    "stroke-dasharray" => st.dash = v.split(',').map(number).collect(),
                    "opacity" => opacity = number(v),
                    _ => {}
                }
            }
            st.fill_opacity *= opacity;
            st.stroke_opacity *= opacity;
        }
        st
    }
}

/// The numbers and commands of a path `d` attribute.
fn path_tokens(d: &str) -> Vec<Result<char, f32>> {
    // Ok(command), Err(number): a Result used as a two-way enum
    let mut out = Vec::new();
    let b = d.as_bytes();
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        if c.is_ascii_alphabetic() && c != b'e' && c != b'E' {
            out.push(Ok(c as char));
            i += 1;
        } else if c == b'-' || c == b'+' || c == b'.' || c.is_ascii_digit() {
            let start = i;
            i += 1;
            while i < b.len() {
                let ch = b[i];
                if ch.is_ascii_digit() || ch == b'.' {
                    i += 1;
                } else if (ch == b'e' || ch == b'E') && i + 1 < b.len() {
                    i += 1;
                    if b[i] == b'-' || b[i] == b'+' {
                        i += 1;
                    }
                } else {
                    break;
                }
            }
            // a second '.' starts a new number
            let text = &d[start..i];
            if let Ok(v) = text.parse::<f32>() {
                out.push(Err(v));
            } else if let Some(pos) = text[1..].find('.') {
                let (a, rest) = text.split_at(pos + 1);
                out.push(Err(a.parse().unwrap_or(0.0)));
                out.push(Err(rest.parse().unwrap_or(0.0)));
            }
        } else {
            i += 1;
        }
    }
    out
}

fn build_path(d: &str) -> Option<tiny_skia::Path> {
    let tokens = path_tokens(d);
    let mut pb = PathBuilder::new();
    let (mut x, mut y) = (0.0f32, 0.0f32);
    let (mut sx, mut sy) = (0.0f32, 0.0f32);
    let mut i = 0;
    let mut cmd = 'M';
    let num = |i: &mut usize| -> Option<f32> {
        match tokens.get(*i) {
            Some(Err(v)) => {
                *i += 1;
                Some(*v)
            }
            _ => None,
        }
    };
    while i < tokens.len() {
        if let Ok(c) = tokens[i] {
            cmd = c;
            i += 1;
            if c == 'Z' || c == 'z' {
                pb.close();
                x = sx;
                y = sy;
                continue;
            }
        }
        let rel = cmd.is_ascii_lowercase();
        let (ox, oy) = if rel { (x, y) } else { (0.0, 0.0) };
        match cmd.to_ascii_uppercase() {
            'M' => {
                let (nx, ny) = (num(&mut i)? + ox, num(&mut i)? + oy);
                pb.move_to(nx, ny);
                (x, y, sx, sy) = (nx, ny, nx, ny);
                cmd = if rel { 'l' } else { 'L' };
            }
            'L' => {
                let (nx, ny) = (num(&mut i)? + ox, num(&mut i)? + oy);
                pb.line_to(nx, ny);
                (x, y) = (nx, ny);
            }
            'H' => {
                let nx = num(&mut i)? + ox;
                pb.line_to(nx, y);
                x = nx;
            }
            'V' => {
                let ny = num(&mut i)? + oy;
                pb.line_to(x, ny);
                y = ny;
            }
            'C' => {
                let (x1, y1) = (num(&mut i)? + ox, num(&mut i)? + oy);
                let (x2, y2) = (num(&mut i)? + ox, num(&mut i)? + oy);
                let (nx, ny) = (num(&mut i)? + ox, num(&mut i)? + oy);
                pb.cubic_to(x1, y1, x2, y2, nx, ny);
                (x, y) = (nx, ny);
            }
            'Q' => {
                let (x1, y1) = (num(&mut i)? + ox, num(&mut i)? + oy);
                let (nx, ny) = (num(&mut i)? + ox, num(&mut i)? + oy);
                pb.quad_to(x1, y1, nx, ny);
                (x, y) = (nx, ny);
            }
            _ => return None,
        }
    }
    pb.finish()
}

fn paint_for(rgb: (u8, u8, u8), opacity: f32) -> Paint<'static> {
    let mut paint = Paint::default();
    paint.set_color_rgba8(
        rgb.0,
        rgb.1,
        rgb.2,
        (opacity.clamp(0.0, 1.0) * 255.0).round() as u8,
    );
    paint.anti_alias = true;
    paint
}

fn draw_path(pixmap: &mut Pixmap, path: &tiny_skia::Path, st: &Style) {
    if let Some(fill) = st.fill {
        let rule = if st.fill_evenodd {
            FillRule::EvenOdd
        } else {
            FillRule::Winding
        };
        pixmap.fill_path(
            path,
            &paint_for(fill, st.fill_opacity),
            rule,
            Transform::identity(),
            None,
        );
    }
    if let Some(stroke_col) = st.stroke {
        if st.stroke_width > 0.0 {
            let mut stroke = Stroke {
                width: st.stroke_width,
                miter_limit: st.miter_limit,
                line_cap: st.cap.unwrap_or(LineCap::Butt),
                line_join: st.join.unwrap_or(LineJoin::Miter),
                ..Stroke::default()
            };
            if !st.dash.is_empty() {
                let mut dash = st.dash.clone();
                if dash.len() % 2 == 1 {
                    let again = dash.clone();
                    dash.extend(again);
                }
                stroke.dash = StrokeDash::new(dash, 0.0);
            }
            pixmap.stroke_path(
                path,
                &paint_for(stroke_col, st.stroke_opacity),
                &stroke,
                Transform::identity(),
                None,
            );
        }
    }
}

fn ellipse_path(cx: f32, cy: f32, rx: f32, ry: f32) -> Option<tiny_skia::Path> {
    let rect = tiny_skia::Rect::from_xywh(cx - rx, cy - ry, 2.0 * rx, 2.0 * ry)?;
    PathBuilder::from_oval(rect)
}

/// Renders the SVG of the drawer.
pub fn svg_to_pixmap(svg: &str) -> Result<Pixmap> {
    let mut pixmap: Option<Pixmap> = None;
    let mut rest = svg;
    while let Some(start) = rest.find('<') {
        rest = &rest[start..];
        if rest.starts_with("<!--") {
            rest = &rest[rest.find("-->").map_or(rest.len(), |e| e + 3)..];
            continue;
        }
        if rest.starts_with("<metadata") {
            rest = &rest[rest.find("</metadata>").map_or(rest.len(), |e| e + 11)..];
            continue;
        }
        if rest.starts_with("<text") {
            // text written as characters needs a font; the PNG drawer always uses glyph outlines
            rest = &rest[rest.find("</text>").map_or(rest.len(), |e| e + 7)..];
            continue;
        }
        let Some(end) = rest.find('>') else { break };
        let tag = &rest[1..end];
        rest = &rest[end + 1..];
        let name = tag
            .split(|c: char| c.is_whitespace() || c == '/')
            .next()
            .unwrap_or("");
        let attrs = attributes(tag);
        match name {
            "svg" => {
                let w = attr(&attrs, "width").map_or(0.0, number);
                let h = attr(&attrs, "height").map_or(0.0, number);
                if w < 1.0 || h < 1.0 {
                    return Err(RDKitError::precondition("the SVG has no size"));
                }
                pixmap = Pixmap::new(w as u32, h as u32);
                if pixmap.is_none() {
                    return Err(RDKitError::precondition(
                        "cannot make a canvas of this size",
                    ));
                }
            }
            "rect" => {
                let Some(pm) = pixmap.as_mut() else { continue };
                let st = Style::from_attributes(&attrs);
                let (x, y) = (
                    attr(&attrs, "x").map_or(0.0, number),
                    attr(&attrs, "y").map_or(0.0, number),
                );
                let (w, h) = (
                    attr(&attrs, "width").map_or(0.0, number),
                    attr(&attrs, "height").map_or(0.0, number),
                );
                if let Some(rect) = tiny_skia::Rect::from_xywh(x, y, w, h) {
                    draw_path(pm, &PathBuilder::from_rect(rect), &st);
                }
            }
            "path" => {
                let Some(pm) = pixmap.as_mut() else { continue };
                let st = Style::from_attributes(&attrs);
                if let Some(path) = attr(&attrs, "d").and_then(build_path) {
                    draw_path(pm, &path, &st);
                }
            }
            "ellipse" => {
                let Some(pm) = pixmap.as_mut() else { continue };
                let st = Style::from_attributes(&attrs);
                let g = |n: &str| attr(&attrs, n).map_or(0.0, number);
                if let Some(path) = ellipse_path(g("cx"), g("cy"), g("rx"), g("ry")) {
                    draw_path(pm, &path, &st);
                }
            }
            _ => {}
        }
    }
    pixmap.ok_or_else(|| RDKitError::precondition("the SVG has no <svg> element"))
}

/// The PNG of the SVG of the drawer.
pub fn svg_to_png(svg: &str) -> Result<Vec<u8>> {
    let pixmap = svg_to_pixmap(svg)?;
    pixmap
        .encode_png()
        .map_err(|e| RDKitError::precondition(format!("cannot encode the PNG: {e}")))
}

/// The molecule tags of the PNG metadata of a drawn molecule (`addMoleculeMetadata` of the Cairo drawer): the pickle,
/// the V3000 mol block if there are conformers, and the CXSMILES. Later molecules get the index as tag suffix.
fn metadata_entries(mol: &ROMol, conf_id: i32, idx: usize) -> Result<Vec<(String, Vec<u8>)>> {
    let params = PngMetadataParams::default();
    let suffix = if idx > 0 {
        idx.to_string()
    } else {
        String::new()
    };
    let mut out = vec![(
        augment_tag_name(&format!("{PKL_TAG}{suffix}")),
        pickle_mol(mol, params.property_flags)?,
    )];
    if mol.num_conformers() > 0 {
        let writer = MolWriterParams {
            include_stereo: true,
            kekulize: false,
            ..MolWriterParams::default()
        };
        let mb = mol_to_v3k_mol_block(mol, &writer, conf_id)?;
        out.push((
            augment_tag_name(&format!("{MOL_TAG}{suffix}")),
            mb.into_bytes(),
        ));
    }
    let smi = mol_to_cx_smiles(
        mol,
        &params.smiles_write_params,
        params.cx_smiles_flags,
        params.restore_bond_dirs,
    )?;
    out.push((
        augment_tag_name(&format!("{SMILES_TAG}{suffix}")),
        smi.into_bytes(),
    ));
    Ok(out)
}

/// The PNG of the SVG of a drawer with the metadata of the molecules it has drawn.
pub fn drawer_png(drawer: &MolDraw2DSVG) -> Result<Vec<u8>> {
    let png = svg_to_png(&drawer.get_drawing_text())?;
    if drawer.metadata_mols.is_empty() {
        return Ok(png);
    }
    let mut metadata = Vec::new();
    for (idx, (mol, conf_id)) in drawer.metadata_mols.iter().enumerate() {
        metadata.extend(metadata_entries(mol, *conf_id, idx)?);
    }
    add_metadata_to_png(&png, &metadata, true)
}

/// `MolDraw2DCairo`: draws like [`MolDraw2DSVG`] (all of its drawing methods are available through `Deref`) and returns
/// the drawing as PNG. The text is always drawn as glyph outlines, as Cairo's FreeType text drawer does.
#[derive(Debug)]
pub struct MolDraw2DCairo {
    inner: MolDraw2DSVG,
}

impl Deref for MolDraw2DCairo {
    type Target = MolDraw2DSVG;
    fn deref(&self) -> &MolDraw2DSVG {
        &self.inner
    }
}

impl DerefMut for MolDraw2DCairo {
    fn deref_mut(&mut self) -> &mut MolDraw2DSVG {
        &mut self.inner
    }
}

impl MolDraw2DCairo {
    /// A drawer for a canvas of the given size in pixels; a size of -1 follows from the molecule
    pub fn new(width: i32, height: i32, panel_width: i32, panel_height: i32) -> Result<Self> {
        Self::with_options(
            width,
            height,
            panel_width,
            panel_height,
            MolDrawOptions::default(),
        )
    }

    pub fn with_options(
        width: i32,
        height: i32,
        panel_width: i32,
        panel_height: i32,
        options: MolDrawOptions,
    ) -> Result<Self> {
        let mut inner =
            MolDraw2DSVG::with_text(width, height, panel_width, panel_height, options, false)?;
        inner.record_metadata = true;
        Ok(MolDraw2DCairo { inner })
    }

    /// `getDrawingText`: the PNG of what has been drawn
    pub fn get_drawing_text(&self) -> Result<Vec<u8>> {
        drawer_png(&self.inner)
    }

    /// `writeDrawingText`
    pub fn write_drawing_text(&self, path: &str) -> Result<()> {
        let png = self.get_drawing_text()?;
        std::fs::write(path, png)
            .map_err(|e| RDKitError::precondition(format!("cannot write {path}: {e}")))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

    const SVG: &str = "<?xml version='1.0'?>\n<svg width='40px' height='30px' viewBox='0 0 40 30'>\n<!-- END OF HEADER -->\n\
        <rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='40' height='30' x='0' y='0'> </rect>\n\
        <path class='bond-0 atom-0 atom-1' d='M 5,15 L 35,15' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:4px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n\
        <ellipse cx='20' cy='8' rx='4' ry='4' style='fill:#0000FF;fill-rule:evenodd;stroke:#0000FF;stroke-width:1px;stroke-opacity:1' />\n\
        <path class='atom-0' d='M 2 2 L 8 2 L 8 6 Z' fill='#00FF00'/>\n</svg>\n";

    fn px(pm: &Pixmap, x: u32, y: u32) -> (u8, u8, u8, u8) {
        let p = pm.pixel(x, y).unwrap().demultiply();
        (p.red(), p.green(), p.blue(), p.alpha())
    }

    #[test]
    fn renders_the_elements_of_the_drawer() {
        let pm = svg_to_pixmap(SVG).unwrap();
        assert_eq!((pm.width(), pm.height()), (40, 30));
        assert_eq!(px(&pm, 0, 29), (255, 255, 255, 255)); // background
        assert_eq!(px(&pm, 20, 15), (255, 0, 0, 255)); // line
        assert_eq!(px(&pm, 20, 8), (0, 0, 255, 255)); // ellipse
        assert_eq!(px(&pm, 6, 3), (0, 255, 0, 255)); // glyph-style path with a fill attribute
        assert_eq!(px(&pm, 20, 22), (255, 255, 255, 255)); // untouched
    }

    #[test]
    fn png_has_the_signature_and_decodes_to_the_same_size() {
        let png = svg_to_png(SVG).unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
        let pm = Pixmap::decode_png(&png).unwrap();
        assert_eq!((pm.width(), pm.height()), (40, 30));
        assert_eq!(px(&pm, 20, 15), (255, 0, 0, 255));
    }

    #[test]
    fn path_numbers_without_separators() {
        let t = path_tokens("M1.5-2.25L.5.25 3e1,4");
        let nums: Vec<f32> = t.iter().filter_map(|x| x.as_ref().err().copied()).collect();
        assert_eq!(nums, vec![1.5, -2.25, 0.5, 0.25, 30.0, 4.0]);
        assert!(build_path("M0 0 C1 1 2 2 3 3 Q4 4 5 5 H6 V7 Z").is_some());
    }

    #[test]
    fn draws_a_molecule_as_png() {
        let mol = mol_from_smiles("CCO", &SmilesParserParams::default())
            .unwrap()
            .unwrap();
        let mut d = MolDraw2DCairo::new(300, 250, -1, -1).unwrap();
        d.draw_molecule(&mol, "ethanol", -1).unwrap();
        d.finish_drawing();
        let png = d.get_drawing_text().unwrap();
        let pm = Pixmap::decode_png(&png).unwrap();
        assert_eq!((pm.width(), pm.height()), (300, 250));
        let pixels: Vec<(u8, u8, u8, u8)> = (0..250)
            .flat_map(|y| (0..300).map(move |x| (x, y)))
            .map(|(x, y)| px(&pm, x, y))
            .collect();
        let ink = pixels
            .iter()
            .filter(|p| **p != (255, 255, 255, 255))
            .count();
        assert!(ink > 300, "{ink}");
        // the oxygen label is red, the carbon bonds and the legend are black
        assert!(pixels.iter().any(|p| p.0 > 200 && p.1 < 60 && p.2 < 60));
        assert!(pixels.iter().any(|p| p.0 < 60 && p.1 < 60 && p.2 < 60));
    }

    #[test]
    fn errors() {
        assert!(svg_to_pixmap("<g/>").is_err());
        assert!(svg_to_pixmap("<svg width='0px' height='3px'>").is_err());
    }
}
