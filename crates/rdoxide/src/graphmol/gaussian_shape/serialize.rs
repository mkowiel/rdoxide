//! The Boost text-archive form of a [`ShapeInput`] (upstream
//! `ShapeInput::serialize`), as stored in SynthonSpace databases.
//!
//! Only the archive layout that upstream produces for `ShapeInput` is
//! supported. Class-information tokens (`0 0`) that Boost writes the first
//! time it meets a class are written and skipped where they occur in that
//! layout.

use super::shape_input::ShapeInput;
use crate::{RDKitError, Result};

const HEADER: &str = "22 serialization::archive 19";

/// `printf("%.17e")`
pub fn fmt_double(v: f64) -> String {
    let s = format!("{v:.17e}");
    let (mant, exp) = s.split_once('e').unwrap_or((&s, "0"));
    let (sign, digits) = match exp.strip_prefix('-') {
        Some(d) => ('-', d),
        None => ('+', exp),
    };
    format!("{mant}e{sign}{digits:0>2}")
}

fn push_doubles(out: &mut String, v: &[f64]) {
    for &x in v {
        out.push(' ');
        out.push_str(&fmt_double(x));
    }
}

fn push_bitset(out: &mut String, bits: &[bool]) {
    let blocks = bits.len().div_ceil(64);
    out.push_str(&format!(" {} {blocks} 0", bits.len()));
    for b in 0..blocks {
        let mut v = 0u64;
        for i in 0..64 {
            if bits.get(b * 64 + i).copied().unwrap_or(false) {
                v |= 1 << i;
            }
        }
        out.push_str(&format!(" {v}"));
    }
}

/// The tokens of a text archive.
#[derive(Debug)]
pub struct Tokens<'a> {
    toks: Vec<&'a str>,
    pos: usize,
}

fn bad(what: &str) -> RDKitError {
    RDKitError::Runtime(format!("Bad shape archive: {what}"))
}

impl<'a> Tokens<'a> {
    pub fn new(text: &'a str) -> Self {
        Tokens {
            toks: text.split_whitespace().collect(),
            pos: 0,
        }
    }
    pub fn token(&mut self) -> Result<&'a str> {
        let t = self.toks.get(self.pos).copied().ok_or_else(|| bad("eof"))?;
        self.pos += 1;
        Ok(t)
    }
    pub fn peek(&self) -> Option<&'a str> {
        self.toks.get(self.pos).copied()
    }
    pub fn usize(&mut self) -> Result<usize> {
        self.token()?.parse().map_err(|_| bad("integer"))
    }
    pub fn i64(&mut self) -> Result<i64> {
        self.token()?.parse().map_err(|_| bad("integer"))
    }
    pub fn f64(&mut self) -> Result<f64> {
        self.token()?.parse().map_err(|_| bad("double"))
    }
    pub fn expect_header(&mut self) -> Result<()> {
        for h in HEADER.split(' ') {
            if self.token()? != h {
                return Err(bad("header"));
            }
        }
        Ok(())
    }
    /// A class-information pair, `0 0`.
    pub fn class_info(&mut self) -> Result<()> {
        self.token()?;
        self.token()?;
        Ok(())
    }
    /// A string: its length then the text.
    pub fn string(&mut self) -> Result<String> {
        let len = self.usize()?;
        // the text may contain spaces, so it can't come from the tokens
        let first = self.token()?;
        if len == 0 {
            return Ok(String::new());
        }
        let mut s = first.to_string();
        while s.len() < len {
            s.push(' ');
            s.push_str(self.token()?);
        }
        Ok(s)
    }
    pub fn doubles(&mut self) -> Result<Vec<f64>> {
        let n = self.usize()?;
        self.token()?; // item version
        (0..n).map(|_| self.f64()).collect()
    }
    pub fn array<const N: usize>(&mut self) -> Result<[f64; N]> {
        let n = self.usize()?;
        if n != N {
            return Err(bad("array size"));
        }
        let mut a = [0.0; N];
        for x in &mut a {
            *x = self.f64()?;
        }
        Ok(a)
    }
    /// Skips the (variable) leading tokens of a `dynamic_bitset` up to its
    /// size, then reads it.
    fn bitset(&mut self, num_bits: usize) -> Result<Vec<bool>> {
        let blocks = num_bits.div_ceil(64);
        loop {
            let ok = self
                .toks
                .get(self.pos)
                .is_some_and(|t| t.parse() == Ok(num_bits))
                && self
                    .toks
                    .get(self.pos + 1)
                    .is_some_and(|t| t.parse() == Ok(blocks))
                && self.toks.get(self.pos + 2) == Some(&"0");
            if ok {
                break;
            }
            self.token()?;
        }
        self.pos += 3;
        let mut bits = vec![false; num_bits];
        for b in 0..blocks {
            let v: u64 = self.token()?.parse().map_err(|_| bad("bitset block"))?;
            for i in 0..64 {
                if b * 64 + i < num_bits {
                    bits[b * 64 + i] = v & (1 << i) != 0;
                }
            }
        }
        Ok(bits)
    }
}

impl ShapeInput {
    /// The Boost text archive of this shape input.
    pub fn to_text(&self) -> String {
        let mut out = String::from(HEADER);
        out.push_str(" 0 0");
        out.push_str(&format!(" {}", self.active_shape));
        out.push_str(&format!(" 0 0 {} 0", self.coords.len()));
        for c in &self.coords {
            out.push_str(&format!(" {} 0", c.len()));
            push_doubles(&mut out, c);
        }
        out.push_str(&format!(" {} 0", self.alphas.len()));
        push_doubles(&mut out, &self.alphas);
        out.push_str(&format!(" {} 0", self.types.len()));
        for t in &self.types {
            out.push_str(&format!(" {t}"));
        }
        out.push_str(&format!(" {} {}", self.num_atoms, self.num_feats));
        for v in [&self.self_overlap_shape_vols, &self.self_overlap_color_vols] {
            out.push_str(&format!(" {} 0", v.len()));
            push_doubles(&mut out, v);
        }
        out.push_str(&format!(" 0 0 {} 0 0 0", self.extreme_points.len()));
        for e in &self.extreme_points {
            out.push_str(" 6");
            for v in e {
                out.push_str(&format!(" {v}"));
            }
        }
        match &self.carbon_radii {
            None => out.push_str(" 0 0 -1"),
            Some(bits) => {
                out.push_str(" 0 0 0 1 0 0");
                push_bitset(&mut out, bits);
            }
        }
        out.push_str(&format!(" {} {}", self.smiles.len(), self.smiles));
        out.push_str(" 1 0\n0");
        push_bitset_tail(&mut out, &self.normalizeds);
        out.push_str("\n1");
        push_bitset_tail(&mut out, &self.normalization_oks);
        out.push_str(&format!(" 0 0 {} 0 0 0", self.canon_rots.len()));
        for r in &self.canon_rots {
            out.push_str(" 9");
            push_doubles(&mut out, r);
        }
        out.push_str(&format!(" 0 0 {} 0 0 0", self.canon_transs.len()));
        for r in &self.canon_transs {
            out.push_str(" 3");
            push_doubles(&mut out, r);
        }
        out.push_str(&format!(" {} 0", self.eigen_values.len()));
        for r in &self.eigen_values {
            out.push_str(" 3");
            push_doubles(&mut out, r);
        }
        out
    }

    /// Reads a shape input from its Boost text archive.
    pub fn from_text(text: &str) -> Result<Self> {
        let mut t = Tokens::new(text);
        t.expect_header()?;
        t.class_info()?;
        let mut s = ShapeInput {
            active_shape: t.usize()?,
            ..Default::default()
        };
        t.class_info()?;
        let n_shapes = t.usize()?;
        t.token()?;
        for _ in 0..n_shapes {
            s.coords.push(t.doubles()?);
        }
        s.alphas = t.doubles()?;
        let n_types = t.usize()?;
        t.token()?;
        for _ in 0..n_types {
            s.types.push(t.i64()? as i32);
        }
        s.num_atoms = t.usize()?;
        s.num_feats = t.usize()?;
        s.self_overlap_shape_vols = t.doubles()?;
        s.self_overlap_color_vols = t.doubles()?;
        t.class_info()?;
        let n = t.usize()?;
        t.token()?;
        t.class_info()?;
        for _ in 0..n {
            if t.usize()? != 6 {
                return Err(bad("extreme points"));
            }
            let mut e = [0usize; 6];
            for v in &mut e {
                *v = t.usize()?;
            }
            s.extreme_points.push(e);
        }
        t.class_info()?;
        if t.peek() == Some("-1") {
            t.token()?;
        } else {
            let n_bits = s.num_atoms + s.num_feats;
            s.carbon_radii = Some(t.bitset(n_bits.max(s.num_atoms))?);
        }
        s.smiles = t.string()?;
        s.normalizeds = t.bitset(n_shapes)?;
        s.normalization_oks = t.bitset(n_shapes)?;
        t.class_info()?;
        let n = t.usize()?;
        t.token()?;
        t.class_info()?;
        for _ in 0..n {
            s.canon_rots.push(t.array::<9>()?);
        }
        t.class_info()?;
        let n = t.usize()?;
        t.token()?;
        t.class_info()?;
        for _ in 0..n {
            s.canon_transs.push(t.array::<3>()?);
        }
        let n = t.usize()?;
        t.token()?;
        for _ in 0..n {
            s.eigen_values.push(t.array::<3>()?);
        }
        Ok(s)
    }
}

/// The size and blocks of a bitset, without its leading tokens.
fn push_bitset_tail(out: &mut String, bits: &[bool]) {
    push_bitset(out, bits);
}
