//! A small TrueType reader that reproduces what RDKit gets from FreeType when
//! it loads glyphs with `FT_LOAD_NO_SCALE | FT_LOAD_NO_BITMAP`: unhinted
//! outlines and advances in font units, the control box computed by
//! `FT_Outline_Get_BBox` and the callback sequence of `FT_Outline_Decompose`.
//!
//! Only `glyf` based (quadratic) fonts are supported, which covers the two
//! fonts built into RDKit (Telex and Roboto) and ordinary `.ttf` files.
//!
//! The fonts are in `crates/rdoxide/data/fonts/`; see LICENSES/README.md for
//! their licences.

use std::borrow::Cow;

use crate::{RDKitError, Result};

static TELEX_REGULAR: &[u8] = include_bytes!("../../../data/fonts/Telex-Regular.ttf");
static ROBOTO_REGULAR: &[u8] = include_bytes!("../../../data/fonts/Roboto-Regular.ttf");

/// An outline point in font units.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OutlinePoint {
    pub x: i32,
    pub y: i32,
    pub on_curve: bool,
}

/// A glyph outline: points plus the index of the last point of each contour.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Outline {
    pub points: Vec<OutlinePoint>,
    pub contour_ends: Vec<usize>,
}

/// Bounding box in font units (`FT_BBox`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct BBox {
    pub x_min: i32,
    pub y_min: i32,
    pub x_max: i32,
    pub y_max: i32,
}

/// One callback of `FT_Outline_Decompose`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PathCommand {
    MoveTo(i32, i32),
    LineTo(i32, i32),
    /// control point, end point
    ConicTo((i32, i32), (i32, i32)),
}

/// A loaded glyph.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Glyph {
    pub outline: Outline,
    /// Horizontal advance in font units.
    pub advance: i32,
}

fn font_error(msg: &str) -> RDKitError {
    RDKitError::Value(format!("invalid or unsupported TrueType font: {msg}"))
}

fn be_u16(d: &[u8], off: usize) -> Result<u16> {
    d.get(off..off + 2)
        .map(|b| u16::from_be_bytes([b[0], b[1]]))
        .ok_or_else(|| font_error("truncated data"))
}

fn be_i16(d: &[u8], off: usize) -> Result<i16> {
    be_u16(d, off).map(|v| v as i16)
}

fn be_u32(d: &[u8], off: usize) -> Result<u32> {
    d.get(off..off + 4)
        .map(|b| u32::from_be_bytes([b[0], b[1], b[2], b[3]]))
        .ok_or_else(|| font_error("truncated data"))
}

/// `FT_MulDiv`: `a * b / c` rounded to nearest, ties away from zero.
fn mul_div(a: i64, b: i64, c: i64) -> i64 {
    let sign = (a < 0) ^ (b < 0) ^ (c < 0);
    let (a, b, c) = (a.abs(), b.abs(), c.abs());
    let r = (a * b + c / 2) / c;
    if sign { -r } else { r }
}

/// `FT_MulFix`: `a * b / 65536` rounded to nearest, ties away from zero.
fn mul_fix(a: i64, b: i64) -> i64 {
    let sign = (a < 0) ^ (b < 0);
    let r = (a.abs() * b.abs() + 0x8000) >> 16;
    if sign { -r } else { r }
}

#[derive(Clone, Copy, Debug)]
enum CmapKind {
    Format0,
    Format4,
    Format6,
    Format12,
}

#[derive(Clone, Copy, Debug)]
struct Cmap {
    offset: usize,
    kind: CmapKind,
}

/// A parsed TrueType font.
#[derive(Clone, Debug)]
pub struct TrueTypeFont {
    data: Cow<'static, [u8]>,
    units_per_em: u16,
    num_glyphs: u16,
    num_h_metrics: u16,
    long_loca: bool,
    glyf: usize,
    loca: usize,
    hmtx: usize,
    cmap: Cmap,
}

const MAX_COMPOSITE_DEPTH: u32 = 8;

impl TrueTypeFont {
    /// The Telex Regular font, RDKit's default (`BuiltinTelexRegular`).
    pub fn builtin_telex_regular() -> Self {
        Self::from_static(TELEX_REGULAR).expect("embedded Telex font is valid")
    }

    /// The Roboto Regular font (`BuiltinRobotoRegular`).
    pub fn builtin_roboto_regular() -> Self {
        Self::from_static(ROBOTO_REGULAR).expect("embedded Roboto font is valid")
    }

    /// Reads a font from a file (the first face of a `.ttf`).
    pub fn from_file(path: &str) -> Result<Self> {
        let data = std::fs::read(path)
            .map_err(|_| RDKitError::Runtime(format!("Font file {path} not found.")))?;
        Self::parse(Cow::Owned(data))
    }

    fn from_static(data: &'static [u8]) -> Result<Self> {
        Self::parse(Cow::Borrowed(data))
    }

    fn parse(data: Cow<'static, [u8]>) -> Result<Self> {
        let d: &[u8] = &data;
        let tag = be_u32(d, 0)?;
        if tag != 0x0001_0000 && tag != 0x7472_7565 {
            return Err(font_error("not a TrueType font with glyf outlines"));
        }
        let num_tables = be_u16(d, 4)? as usize;
        let table = |name: &[u8; 4]| -> Result<usize> {
            for i in 0..num_tables {
                let rec = 12 + 16 * i;
                if d.get(rec..rec + 4) == Some(name) {
                    return Ok(be_u32(d, rec + 8)? as usize);
                }
            }
            Err(font_error(&format!(
                "missing '{}' table",
                String::from_utf8_lossy(name)
            )))
        };
        let head = table(b"head")?;
        let maxp = table(b"maxp")?;
        let hhea = table(b"hhea")?;
        let units_per_em = be_u16(d, head + 18)?;
        if units_per_em == 0 {
            return Err(font_error("zero units per em"));
        }
        let long_loca = be_i16(d, head + 50)? != 0;
        let num_glyphs = be_u16(d, maxp + 4)?;
        let num_h_metrics = be_u16(d, hhea + 34)?;
        let cmap = Self::find_cmap(d, table(b"cmap")?)?;
        Ok(TrueTypeFont {
            units_per_em,
            num_glyphs,
            num_h_metrics,
            long_loca,
            glyf: table(b"glyf")?,
            loca: table(b"loca")?,
            hmtx: table(b"hmtx")?,
            cmap,
            data,
        })
    }

    /// Picks a Unicode subtable the way FreeType's default charmap does:
    /// full repertoire (format 12) first, then the BMP one.
    fn find_cmap(d: &[u8], cmap: usize) -> Result<Cmap> {
        let n = be_u16(d, cmap + 2)? as usize;
        let mut best: Option<(u8, Cmap)> = None;
        for i in 0..n {
            let rec = cmap + 4 + 8 * i;
            let (platform, encoding) = (be_u16(d, rec)?, be_u16(d, rec + 2)?);
            let offset = cmap + be_u32(d, rec + 4)? as usize;
            let kind = match be_u16(d, offset)? {
                0 => CmapKind::Format0,
                4 => CmapKind::Format4,
                6 => CmapKind::Format6,
                12 => CmapKind::Format12,
                _ => continue,
            };
            let rank = match (platform, encoding) {
                (3, 10) | (0, 4) | (0, 6) => 3,
                (3, 1) | (0, _) => 2,
                _ => continue,
            };
            if best.is_none_or(|(r, _)| rank > r) {
                best = Some((rank, Cmap { offset, kind }));
            }
        }
        best.map(|(_, c)| c)
            .ok_or_else(|| font_error("no Unicode character map"))
    }

    pub fn units_per_em(&self) -> u16 {
        self.units_per_em
    }

    /// `1 / units_per_EM`.
    pub fn em_scale(&self) -> f64 {
        1.0 / f64::from(self.units_per_em)
    }

    /// Glyph index of a Unicode code point (0, `.notdef`, if unmapped).
    pub fn glyph_index(&self, code: u32) -> u32 {
        self.lookup(code).unwrap_or(0)
    }

    fn lookup(&self, code: u32) -> Result<u32> {
        let d: &[u8] = &self.data;
        let o = self.cmap.offset;
        match self.cmap.kind {
            CmapKind::Format0 => {
                if code < 256 {
                    Ok(u32::from(*d.get(o + 6 + code as usize).unwrap_or(&0)))
                } else {
                    Ok(0)
                }
            }
            CmapKind::Format6 => {
                let first = u32::from(be_u16(d, o + 6)?);
                let count = u32::from(be_u16(d, o + 8)?);
                if code >= first && code < first + count {
                    Ok(u32::from(be_u16(d, o + 10 + 2 * (code - first) as usize)?))
                } else {
                    Ok(0)
                }
            }
            CmapKind::Format4 => {
                if code > 0xffff {
                    return Ok(0);
                }
                let seg_x2 = be_u16(d, o + 6)? as usize;
                let ends = o + 14;
                let starts = ends + seg_x2 + 2;
                let deltas = starts + seg_x2;
                let ranges = deltas + seg_x2;
                for s in 0..seg_x2 / 2 {
                    let end = u32::from(be_u16(d, ends + 2 * s)?);
                    if code > end {
                        continue;
                    }
                    let start = u32::from(be_u16(d, starts + 2 * s)?);
                    if code < start {
                        return Ok(0);
                    }
                    let delta = u32::from(be_u16(d, deltas + 2 * s)?);
                    let range = be_u16(d, ranges + 2 * s)? as usize;
                    if range == 0 {
                        return Ok((code + delta) & 0xffff);
                    }
                    let at = ranges + 2 * s + range + 2 * (code - start) as usize;
                    let g = u32::from(be_u16(d, at)?);
                    return Ok(if g == 0 { 0 } else { (g + delta) & 0xffff });
                }
                Ok(0)
            }
            CmapKind::Format12 => {
                let n = be_u32(d, o + 12)? as usize;
                for g in 0..n {
                    let rec = o + 16 + 12 * g;
                    let (start, end) = (be_u32(d, rec)?, be_u32(d, rec + 4)?);
                    if code >= start && code <= end {
                        return Ok(be_u32(d, rec + 8)? + (code - start));
                    }
                }
                Ok(0)
            }
        }
    }

    /// Loads the glyph for a code point like `FT_Load_Char`: unmapped
    /// characters give the `.notdef` glyph.
    pub fn load_char(&self, code: u32) -> Result<Glyph> {
        self.load_glyph(self.glyph_index(code))
    }

    pub fn load_glyph(&self, index: u32) -> Result<Glyph> {
        if index >= u32::from(self.num_glyphs) {
            return Err(font_error("glyph index out of range"));
        }
        let outline = self.glyph_outline(index, 0)?;
        Ok(Glyph {
            outline,
            advance: self.advance(index)?,
        })
    }

    fn advance(&self, index: u32) -> Result<i32> {
        let n = u32::from(self.num_h_metrics);
        if n == 0 {
            return Ok(0);
        }
        let i = index.min(n - 1) as usize;
        Ok(i32::from(be_u16(&self.data, self.hmtx + 4 * i)?))
    }

    fn glyph_range(&self, index: u32) -> Result<Option<(usize, usize)>> {
        let d: &[u8] = &self.data;
        let i = index as usize;
        let (a, b) = if self.long_loca {
            (
                be_u32(d, self.loca + 4 * i)? as usize,
                be_u32(d, self.loca + 4 * i + 4)? as usize,
            )
        } else {
            (
                2 * be_u16(d, self.loca + 2 * i)? as usize,
                2 * be_u16(d, self.loca + 2 * i + 2)? as usize,
            )
        };
        Ok((b > a).then_some((self.glyf + a, self.glyf + b)))
    }

    fn glyph_outline(&self, index: u32, depth: u32) -> Result<Outline> {
        if depth > MAX_COMPOSITE_DEPTH {
            return Err(font_error("composite glyphs nested too deeply"));
        }
        let Some((start, _end)) = self.glyph_range(index)? else {
            return Ok(Outline::default());
        };
        let d: &[u8] = &self.data;
        let n_contours = be_i16(d, start)?;
        if n_contours >= 0 {
            self.simple_outline(start, n_contours as usize)
        } else {
            self.composite_outline(start + 10, depth)
        }
    }

    fn simple_outline(&self, start: usize, n_contours: usize) -> Result<Outline> {
        let d: &[u8] = &self.data;
        let mut p = start + 10;
        let mut contour_ends = Vec::with_capacity(n_contours);
        for i in 0..n_contours {
            contour_ends.push(be_u16(d, p + 2 * i)? as usize);
        }
        p += 2 * n_contours;
        let n_points = contour_ends.last().map_or(0, |&e| e + 1);
        let instr = be_u16(d, p)? as usize;
        p += 2 + instr;

        let mut flags = Vec::with_capacity(n_points);
        while flags.len() < n_points {
            let f = *d.get(p).ok_or_else(|| font_error("truncated glyph"))?;
            p += 1;
            flags.push(f);
            if f & 8 != 0 {
                let repeat = *d.get(p).ok_or_else(|| font_error("truncated glyph"))?;
                p += 1;
                for _ in 0..repeat {
                    flags.push(f);
                }
            }
        }
        flags.truncate(n_points);

        let mut read_coords = |short_bit: u8, same_bit: u8| -> Result<Vec<i32>> {
            let mut v = 0i32;
            let mut out = Vec::with_capacity(n_points);
            for &f in &flags {
                if f & short_bit != 0 {
                    let b = i32::from(*d.get(p).ok_or_else(|| font_error("truncated glyph"))?);
                    p += 1;
                    v += if f & same_bit != 0 { b } else { -b };
                } else if f & same_bit == 0 {
                    v += i32::from(be_i16(d, p)?);
                    p += 2;
                }
                out.push(v);
            }
            Ok(out)
        };
        let xs = read_coords(2, 16)?;
        let ys = read_coords(4, 32)?;
        let points = (0..n_points)
            .map(|i| OutlinePoint {
                x: xs[i],
                y: ys[i],
                on_curve: flags[i] & 1 != 0,
            })
            .collect();
        Ok(Outline {
            points,
            contour_ends,
        })
    }

    fn composite_outline(&self, mut p: usize, depth: u32) -> Result<Outline> {
        const ARG_WORDS: u16 = 0x1;
        const ARGS_ARE_XY: u16 = 0x2;
        const HAVE_SCALE: u16 = 0x8;
        const MORE: u16 = 0x20;
        const HAVE_XY_SCALE: u16 = 0x40;
        const HAVE_2X2: u16 = 0x80;
        const SCALED_OFFSET: u16 = 0x800;
        let d: &[u8] = &self.data;
        let mut result = Outline::default();
        loop {
            let flags = be_u16(d, p)?;
            let component = u32::from(be_u16(d, p + 2)?);
            p += 4;
            if flags & ARGS_ARE_XY == 0 {
                return Err(font_error("point-matched composite glyphs"));
            }
            if flags & SCALED_OFFSET != 0 {
                return Err(font_error("scaled composite glyph offsets"));
            }
            let (dx, dy) = if flags & ARG_WORDS != 0 {
                let r = (i32::from(be_i16(d, p)?), i32::from(be_i16(d, p + 2)?));
                p += 4;
                r
            } else {
                let r = (
                    i32::from(*d.get(p).ok_or_else(|| font_error("truncated"))? as i8),
                    i32::from(*d.get(p + 1).ok_or_else(|| font_error("truncated"))? as i8),
                );
                p += 2;
                r
            };
            // 2x2 matrix in 16.16: [xx xy; yx yy], x' = xx*x + yx*y
            let f2dot14 = |off: usize| -> Result<i64> { Ok(i64::from(be_i16(d, off)?) * 4) };
            let (mut xx, mut xy, mut yx, mut yy) = (0x10000i64, 0i64, 0i64, 0x10000i64);
            if flags & HAVE_SCALE != 0 {
                xx = f2dot14(p)?;
                yy = xx;
                p += 2;
            } else if flags & HAVE_XY_SCALE != 0 {
                xx = f2dot14(p)?;
                yy = f2dot14(p + 2)?;
                p += 4;
            } else if flags & HAVE_2X2 != 0 {
                xx = f2dot14(p)?;
                xy = f2dot14(p + 2)?;
                yx = f2dot14(p + 4)?;
                yy = f2dot14(p + 6)?;
                p += 8;
            }
            let mut part = self.glyph_outline(component, depth + 1)?;
            for pt in &mut part.points {
                let (x, y) = (i64::from(pt.x), i64::from(pt.y));
                pt.x = (mul_fix(x, xx) + mul_fix(y, yx)) as i32 + dx;
                pt.y = (mul_fix(x, xy) + mul_fix(y, yy)) as i32 + dy;
            }
            let base = result.points.len();
            result
                .contour_ends
                .extend(part.contour_ends.iter().map(|e| e + base));
            result.points.extend(part.points);
            if flags & MORE == 0 {
                return Ok(result);
            }
        }
    }
}

impl Outline {
    /// `FT_Outline_Get_BBox`: the exact bounding box of the outline, found by
    /// checking the extrema of the quadratic segments whose control point lies
    /// outside the box of the on-curve points.
    pub fn bbox(&self) -> BBox {
        let Some(first) = self.points.first() else {
            return BBox::default();
        };
        if self.contour_ends.is_empty() {
            return BBox::default();
        }
        let mut cbox = BBox {
            x_min: first.x,
            y_min: first.y,
            x_max: first.x,
            y_max: first.y,
        };
        let mut bbox = cbox;
        for pt in &self.points {
            cbox.x_min = cbox.x_min.min(pt.x);
            cbox.y_min = cbox.y_min.min(pt.y);
            cbox.x_max = cbox.x_max.max(pt.x);
            cbox.y_max = cbox.y_max.max(pt.y);
            if pt.on_curve {
                bbox.x_min = bbox.x_min.min(pt.x);
                bbox.y_min = bbox.y_min.min(pt.y);
                bbox.x_max = bbox.x_max.max(pt.x);
                bbox.y_max = bbox.y_max.max(pt.y);
            }
        }
        if cbox == bbox {
            return bbox;
        }
        let mut last = (0, 0);
        for cmd in self.decompose() {
            match cmd {
                PathCommand::MoveTo(x, y) => {
                    bbox.x_min = bbox.x_min.min(x);
                    bbox.y_min = bbox.y_min.min(y);
                    bbox.x_max = bbox.x_max.max(x);
                    bbox.y_max = bbox.y_max.max(y);
                    last = (x, y);
                }
                PathCommand::LineTo(x, y) => last = (x, y),
                PathCommand::ConicTo(c, to) => {
                    if c.0 < bbox.x_min || c.0 > bbox.x_max {
                        conic_check(last.0, c.0, to.0, &mut bbox.x_min, &mut bbox.x_max);
                    }
                    if c.1 < bbox.y_min || c.1 > bbox.y_max {
                        conic_check(last.1, c.1, to.1, &mut bbox.y_min, &mut bbox.y_max);
                    }
                    last = to;
                }
            }
        }
        bbox
    }

    /// `FT_Outline_Decompose`: the move/line/conic sequence for the outline.
    /// Every contour ends with a line back to its start.
    pub fn decompose(&self) -> Vec<PathCommand> {
        let mut cmds = Vec::new();
        let mut first = 0usize;
        for &last in &self.contour_ends {
            if last < first || last >= self.points.len() {
                break;
            }
            self.decompose_contour(first, last, &mut cmds);
            first = last + 1;
        }
        cmds
    }

    fn decompose_contour(&self, first: usize, last: usize, cmds: &mut Vec<PathCommand>) {
        let pts = &self.points[first..=last];
        let v = |p: &OutlinePoint| (p.x, p.y);
        let mid = |a: (i32, i32), b: (i32, i32)| ((a.0 + b.0) / 2, (a.1 + b.1) / 2);
        let mut limit = pts.len() - 1;
        let mut v_start = v(&pts[0]);
        let v_last = v(&pts[limit]);
        // index of the point before the first one to process
        let mut cur: isize = 0;
        if !pts[0].on_curve {
            // the first point is a control point
            if pts[limit].on_curve {
                v_start = v_last;
                limit -= 1;
            } else {
                v_start = mid(v_start, v_last);
            }
            cur = -1;
        }
        cmds.push(PathCommand::MoveTo(v_start.0, v_start.1));
        // a contour of a single control point has nothing to process
        let limit = limit as isize;
        'outer: while cur < limit {
            cur += 1;
            let pt = &pts[cur as usize];
            if pt.on_curve {
                cmds.push(PathCommand::LineTo(pt.x, pt.y));
                continue;
            }
            let mut control = v(pt);
            loop {
                if cur < limit {
                    cur += 1;
                    let next = &pts[cur as usize];
                    if next.on_curve {
                        cmds.push(PathCommand::ConicTo(control, v(next)));
                        continue 'outer;
                    }
                    let middle = mid(control, v(next));
                    cmds.push(PathCommand::ConicTo(control, middle));
                    control = v(next);
                } else {
                    cmds.push(PathCommand::ConicTo(control, v_start));
                    break 'outer;
                }
            }
        }
        cmds.push(PathCommand::LineTo(v_start.0, v_start.1));
    }
}

/// `BBox_Conic_Check`: grows `[min, max]` by the extremum of the quadratic
/// curve with the given end and control coordinates, if there is one.
fn conic_check(y1: i32, y2: i32, y3: i32, min: &mut i32, max: &mut i32) {
    let (y1, y2, y3) = (i64::from(y1), i64::from(y2), i64::from(y3));
    if y1 <= y3 && y2 == y1 {
        return;
    }
    if y1 < y3 {
        if y2 >= y1 && y2 <= y3 {
            return;
        }
    } else if y2 >= y3 && y2 <= y1 {
        return;
    }
    let (a, b) = (y1 - y2, y3 - y2);
    if a + b == 0 {
        return;
    }
    let ext = (y2 + mul_div(a, b, a + b)) as i32;
    if ext < *min {
        *min = ext;
    } else if ext > *max {
        *max = ext;
    }
}
