//! Port of the V2000 Molfile reader of the InChI library: the data model of
//! `mol_fmt.h`, `mol_fmt2.c` (field reader, coordinate normalization helpers)
//! and the header/counts/atoms/bonds/STEXT readers of `mol_fmt1.c`.
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust. Licensed under the MIT license
//! (see LICENSES/InChI-MIT.txt).
//!
//! Lines are byte vectors without the terminating zero; positions past the end
//! read as the C string terminator.

use super::ichi_io::{InchiIostream, inchi_fgets_lf};
use super::ichierr::add_error_message;
use super::ichiprt1::MolCoord;
use super::util::{
    dotify_non_printable_chars, extract_charges_and_radicals, get_atomic_mass, inchi_memicmp,
    lrtrim, mystrncpy, normalize_string, remove_one_lf,
};

pub const SD_FMT_END_OF_DATA: &[u8] = b"$$$$";
pub const MOL_FMT_INPLINELEN: i32 = 204;
pub const MOL_FMT_MAXLINELEN: usize = 200;
pub const MOL_FMT_MAX_VALUE_LEN: i32 = 32;
pub const RADICAL_DOUBLET: i8 = 2;
pub const INCHI_NAME: &[u8] = b"InChI";

/// the kinds of `MolfileReadField` numeric data
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum IntKind {
    /// `MOL_FMT_CHAR_INT_DATA`
    Char,
    /// `MOL_FMT_SHORT_INT_DATA`
    Short,
    /// `MOL_FMT_LONG_INT_DATA`
    Long,
}

#[derive(Clone, Debug, Default)]
pub struct MolFmtHeaderBlock {
    pub molname: Vec<u8>,
    pub line2: Vec<u8>,
    pub user_initls: Vec<u8>,
    pub prog_name: Vec<u8>,
    pub month: i8,
    pub day: i8,
    pub year: i8,
    pub hour: i8,
    pub minute: i8,
    pub dim_code: Vec<u8>,
    pub scaling_factor1: i16,
    pub scaling_factor2: f64,
    pub energy: f64,
    pub internal_regno: i64,
    pub comment: Vec<u8>,
}

#[derive(Clone, Debug, Default)]
pub struct MolFmtAtom {
    pub fx: f64,
    pub fy: f64,
    pub fz: f64,
    pub symbol: Vec<u8>,
    pub mass_difference: i8,
    pub charge: i8,
    pub radical: i8,
    pub stereo_parity: i8,
    pub valence: i8,
    pub my_n_imp_h: i8,
    pub display_tom: i8,
    pub atom_aliased_flag: i8,
}

#[derive(Clone, Debug, Default)]
pub struct MolFmtBond {
    pub atnum1: i16,
    pub atnum2: i16,
    pub bond_type: i8,
    pub bond_stereo: i8,
}

/// `INT_ARRAY` (the Rust vector replaces the growth bookkeeping)
pub type IntArray = Vec<i32>;
/// `NUM_LISTS`
pub type NumLists = Vec<Vec<i32>>;

#[derive(Clone, Debug, Default)]
pub struct MolFmtSgroup {
    pub id: i32,
    pub type_: i32,
    pub subtype: i32,
    pub conn: i32,
    pub label: i32,
    pub xbr1: [f64; 4],
    pub xbr2: [f64; 4],
    pub smt: Vec<u8>,
    pub alist: IntArray,
    pub blist: IntArray,
}

#[derive(Clone, Debug, Default)]
pub struct MolFmtV3000 {
    pub n_non_star_atoms: i32,
    pub n_star_atoms: i32,
    pub atom_index_orig: Vec<i32>,
    pub atom_index_fin: Vec<i32>,
    pub n_sgroups: i32,
    pub n_3d_constraints: i32,
    pub n_collections: i32,
    pub n_non_haptic_bonds: i32,
    pub n_haptic_bonds: i32,
    pub haptic_bonds: NumLists,
    pub n_steabs: i32,
    pub steabs: NumLists,
    pub n_sterel: i32,
    pub sterel: NumLists,
    pub n_sterac: i32,
    pub sterac: NumLists,
}

#[derive(Clone, Debug, Default)]
pub struct MolFmtCtab {
    pub n_atoms: i32,
    pub n_bonds: i32,
    pub chiral_flag: i8,
    pub n_stext_entries: i16,
    pub n_property_lines: i16,
    pub follow_inchi_1_treating_iso_mass: i16,
    pub version_string: Vec<u8>,
    pub atoms: Option<Vec<MolFmtAtom>>,
    pub bonds: Option<Vec<MolFmtBond>>,
    pub coords: Option<Vec<MolCoord>>,
    pub sgroups: Vec<MolFmtSgroup>,
    /// `MolFmtSgroups_Free` zeroes the container: nothing can be appended later
    pub sgroups_freed: bool,
    pub v3000: Option<Box<MolFmtV3000>>,
}

#[derive(Clone, Debug, Default)]
pub struct MolFmtData {
    pub hdr: MolFmtHeaderBlock,
    pub ctab: MolFmtCtab,
}

fn at(s: &[u8], i: usize) -> u8 {
    s.get(i).copied().unwrap_or(0)
}

fn is_space(c: u8) -> bool {
    matches!(c, b' ' | 0x09..=0x0d)
}

/// the part of `s` before the first zero byte (the C string)
fn c_str(s: &[u8]) -> &[u8] {
    &s[..s.iter().position(|&c| c == 0).unwrap_or(s.len())]
}

/// `strtol(s, &end, 10)` with `long` of 64 bits: (value, bytes consumed)
pub fn c_strtol10(s: &[u8]) -> (i64, usize) {
    let s = c_str(s);
    let mut i = 0;
    while i < s.len() && is_space(s[i]) {
        i += 1;
    }
    let mut neg = false;
    if i < s.len() && (s[i] == b'+' || s[i] == b'-') {
        neg = s[i] == b'-';
        i += 1;
    }
    let ds = i;
    let mut v: i128 = 0;
    while i < s.len() && s[i].is_ascii_digit() {
        v = (v * 10 + i128::from(s[i] - b'0')).min(i128::from(i64::MAX) + 1);
        i += 1;
    }
    if i == ds {
        return (0, 0);
    }
    let v = if neg { -v } else { v };
    (
        v.clamp(i128::from(i64::MIN), i128::from(i64::MAX)) as i64,
        i,
    )
}

/// `strtod(s, &end)`: decimal, hexadecimal, `inf`/`infinity` and `nan`
pub fn c_strtod(s: &[u8]) -> (f64, usize) {
    let s = c_str(s);
    let n = s.len();
    let mut i = 0;
    while i < n && is_space(s[i]) {
        i += 1;
    }
    let mut neg = false;
    if i < n && (s[i] == b'+' || s[i] == b'-') {
        neg = s[i] == b'-';
        i += 1;
    }
    let sign = |v: f64| if neg { -v } else { v };
    let lower = |k: usize| s.get(k).map(|c| c.to_ascii_lowercase());
    let word = |k: usize, w: &[u8]| (0..w.len()).all(|j| lower(k + j) == Some(w[j]));
    if word(i, b"inf") {
        let len = if word(i, b"infinity") { 8 } else { 3 };
        return (sign(f64::INFINITY), i + len);
    }
    if word(i, b"nan") {
        let mut k = i + 3;
        if at(s, k) == b'(' {
            let mut j = k + 1;
            while j < n && (s[j].is_ascii_alphanumeric() || s[j] == b'_') {
                j += 1;
            }
            if at(s, j) == b')' {
                k = j + 1;
            }
        }
        return (sign(f64::NAN), k);
    }
    if at(s, i) == b'0'
        && matches!(at(s, i + 1), b'x' | b'X')
        && (at(s, i + 2).is_ascii_hexdigit()
            || (at(s, i + 2) == b'.' && at(s, i + 3).is_ascii_hexdigit()))
    {
        let mut k = i + 2;
        let mut m = 0f64;
        let mut exp = 0i32;
        while at(s, k).is_ascii_hexdigit() {
            m = m * 16.0 + f64::from((at(s, k) as char).to_digit(16).unwrap());
            k += 1;
        }
        if at(s, k) == b'.' {
            k += 1;
            while at(s, k).is_ascii_hexdigit() {
                m = m * 16.0 + f64::from((at(s, k) as char).to_digit(16).unwrap());
                exp -= 4;
                k += 1;
            }
        }
        if matches!(at(s, k), b'p' | b'P') {
            let mut j = k + 1;
            let mut eneg = false;
            if matches!(at(s, j), b'+' | b'-') {
                eneg = at(s, j) == b'-';
                j += 1;
            }
            if at(s, j).is_ascii_digit() {
                let mut e: i32 = 0;
                while at(s, j).is_ascii_digit() {
                    e = e
                        .saturating_mul(10)
                        .saturating_add(i32::from(at(s, j) - b'0'));
                    j += 1;
                }
                exp = exp.saturating_add(if eneg { -e } else { e });
                k = j;
            }
        }
        return (sign(m * 2f64.powi(exp)), k);
    }
    let ds = i;
    let mut digits = 0;
    while at(s, i).is_ascii_digit() {
        i += 1;
        digits += 1;
    }
    if at(s, i) == b'.' {
        let mut k = i + 1;
        let mut f = 0;
        while at(s, k).is_ascii_digit() {
            k += 1;
            f += 1;
        }
        if digits + f > 0 {
            digits += f;
            i = k;
        }
    }
    if digits == 0 {
        return (0.0, 0);
    }
    if matches!(at(s, i), b'e' | b'E') {
        let mut k = i + 1;
        if matches!(at(s, k), b'+' | b'-') {
            k += 1;
        }
        let es = k;
        while at(s, k).is_ascii_digit() {
            k += 1;
        }
        if k > es {
            i = k;
        }
    }
    let txt = std::str::from_utf8(&s[ds..i]).unwrap_or("0");
    let v: f64 = txt.parse().unwrap_or(0.0);
    (sign(v), i)
}

/// `MolfileStrnread`: copies up to `len` bytes of `source` (stopping at the
/// zero) and finds where the trailing blanks start. Returns the copy, the
/// length of the copy and the index of the first trailing blank.
fn molfile_strnread(source: &[u8], len: i32) -> (Vec<u8>, usize, usize) {
    let len = len.max(0) as usize;
    let mut dest = Vec::new();
    for k in 0..len {
        let c = at(source, k);
        if c == 0 {
            break;
        }
        dest.push(c);
    }
    let n = dest.len();
    let mut i = n as i64 - 1;
    while i >= 0 && is_space(dest[i as usize]) {
        i -= 1;
    }
    (dest, n, (i + 1) as usize)
}

/// `MolfileReadField` with `MOL_FMT_STRING_DATA`: `data` gets the string;
/// returns its length
pub fn molfile_read_field_str(
    data: &mut Vec<u8>,
    field_len: i32,
    line: &[u8],
    pos: &mut usize,
) -> i32 {
    let mut i = 0usize;
    while (i as i32) < field_len && at(line, *pos + i) != 0 && is_space(at(line, *pos + i)) {
        i += 1;
    }
    let (mut d, len, q) =
        molfile_strnread(&line[(*pos + i).min(line.len())..], field_len - i as i32);
    d.truncate(q);
    *data = d;
    *pos += len + i;
    q as i32
}

/// `MolfileReadField` with the integer types: (value, return code)
pub fn molfile_read_field_int(
    kind: IntKind,
    field_len: i32,
    line: &[u8],
    pos: &mut usize,
) -> (i64, i32) {
    let mut ret = 1;
    let mut ldata: i64 = 0;
    let rest = &line[(*pos).min(line.len())..];
    if field_len > MOL_FMT_MAX_VALUE_LEN {
        ret = -1;
    } else if field_len <= 0 {
        let (v, used) = c_strtol10(rest);
        ldata = v;
        *pos += used;
        if used == 0 {
            ret = 0;
        }
    } else {
        let (s, len, q) = molfile_strnread(rest, field_len);
        *pos += len;
        if len == 0 || q == 0 {
            ret = 0;
        } else {
            let (v, used) = c_strtol10(&s[..q]);
            ldata = v;
            if used != q {
                ret = -1;
            }
        }
    }
    let value = match kind {
        IntKind::Char => {
            if (i64::from(i8::MIN)..=i64::from(i8::MAX)).contains(&ldata) {
                ldata
            } else {
                ret = -1;
                0
            }
        }
        IntKind::Short => {
            if (i64::from(i16::MIN)..=i64::from(i16::MAX)).contains(&ldata) {
                ldata
            } else {
                ret = -1;
                0
            }
        }
        IntKind::Long => {
            if i64::MIN < ldata && ldata < i64::MAX {
                ldata
            } else {
                ret = -1;
                0
            }
        }
    };
    (value, ret)
}

/// the value and the return code of the `double` / `float` field readers
fn read_field_real(field_len: i32, line: &[u8], pos: &mut usize) -> (f64, i32) {
    let mut ret = 1;
    let mut ddata = 0.0;
    let rest = &line[(*pos).min(line.len())..];
    if field_len > MOL_FMT_MAX_VALUE_LEN {
        ret = -1;
    } else if field_len <= 0 {
        let (v, used) = c_strtod(rest);
        ddata = v;
        *pos += used;
        if used == 0 {
            ret = 0;
        }
    } else {
        let (s, len, q) = molfile_strnread(rest, field_len);
        *pos += len;
        if len == 0 || q == 0 {
            ret = 0;
        } else {
            let (v, used) = c_strtod(&s[..q]);
            ddata = v;
            if used != q {
                ret = -1;
            }
        }
    }
    (ddata, ret)
}

/// `MolfileReadField` with `MOL_FMT_DOUBLE_DATA`
pub fn molfile_read_field_double(field_len: i32, line: &[u8], pos: &mut usize) -> (f64, i32) {
    let (d, ret) = read_field_real(field_len, line, pos);
    if !d.is_infinite() {
        (d, ret)
    } else {
        (0.0, -1)
    }
}

/// `MolfileReadField` with `MOL_FMT_FLOAT_DATA`
pub fn molfile_read_field_float(field_len: i32, line: &[u8], pos: &mut usize) -> (f32, i32) {
    let (d, ret) = read_field_real(field_len, line, pos);
    if d.abs() <= f64::from(f32::MIN_POSITIVE) {
        (0.0, ret)
    } else if d.abs() >= f64::from(f32::MAX) {
        (0.0, -1)
    } else {
        (d as f32, ret)
    }
}

/// `MolfileReadField` with `MOL_FMT_JUMP_TO_RIGHT`
pub fn molfile_read_field_jump(field_len: i32, line: &[u8], pos: &mut usize) -> i32 {
    let mut i = 0usize;
    while (i as i32) < field_len && at(line, *pos + i) != 0 {
        i += 1;
    }
    *pos += i;
    i as i32
}

/// `MolfileExtractStrucNum`: the number in a name line like "Structure #22"
pub fn molfile_extract_struc_num(hdr: &MolFmtHeaderBlock) -> i64 {
    const S_STRUCT: &[u8] = b"Structure #";
    let pad = |s: &[u8], n: usize| {
        let mut v = c_str(s).to_vec();
        v.resize(n.max(v.len()), 0);
        v
    };
    let name = pad(&hdr.molname, S_STRUCT.len());
    if inchi_memicmp(&name, S_STRUCT, S_STRUCT.len()) != 0 {
        return 0;
    }
    let rest = &name[S_STRUCT.len()..];
    let (num, used) = c_strtol10(rest);
    let line2 = pad(&hdr.line2, INCHI_NAME.len());
    let tail = &line2[INCHI_NAME.len()..];
    let has_sdfile = tail
        .windows(b"SDfile Output".len())
        .any(|w| w == b"SDfile Output");
    let end_ok = at(rest, used) == 0;
    if !end_ok || inchi_memicmp(&line2, INCHI_NAME, INCHI_NAME.len()) != 0 || !has_sdfile {
        return 0;
    }
    num
}

/// `MolfileHasNoChemStruc`
pub fn molfile_has_no_chem_struc(mfdata: Option<&MolFmtData>) -> bool {
    let Some(m) = mfdata else {
        return true;
    };
    if m.ctab.atoms.is_none() {
        return true;
    }
    if m.ctab.n_atoms <= 0 {
        return true;
    }
    m.ctab.n_bonds > 0 && m.ctab.bonds.is_none()
}

/// result of [`molfile_get_xyz_dim_and_norm_factors`]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct XyzNorm {
    pub num_dimensions: i32,
    pub x0: f64,
    pub y0: f64,
    pub z0: f64,
    pub xmin: f64,
    pub ymin: f64,
    pub zmin: f64,
    pub scaler: f64,
}

/// `MolfileGetXYZDimAndNormFactors`
pub fn molfile_get_xyz_dim_and_norm_factors(
    mfdata: &MolFmtData,
    find_norm_factors: bool,
    err: &mut i32,
    p_str_err: &mut String,
) -> XyzNorm {
    const MIN_C: f64 = 0.0;
    const MAX_C: f64 = 256.0;
    const MAX_AVE: f64 = 20.0;
    const MIN_AVE: f64 = 10.0;
    let fmax = |a: f64, b: f64| if a > b { a } else { b };
    let fmin = |a: f64, b: f64| if a < b { a } else { b };
    let (mut max_x, mut max_y, mut max_z) = (-1.0e32f64, -1.0e32f64, -1.0e32f64);
    let (mut min_x, mut min_y, mut min_z) = (1.0e32f64, 1.0e32f64, 1.0e32f64);
    let macheps = 1.0e-10;
    let small_coeff = 0.00001;
    let mut coeff = 1.0;
    let finish = |num_dimensions: i32, min: (f64, f64, f64), coeff: f64| XyzNorm {
        num_dimensions,
        x0: min.0,
        y0: min.1,
        z0: min.2,
        xmin: MIN_C,
        ymin: MIN_C,
        zmin: MIN_C,
        scaler: coeff,
    };
    if molfile_has_no_chem_struc(Some(mfdata)) {
        return finish(0, (min_x, min_y, min_z), coeff);
    }
    let atoms = mfdata.ctab.atoms.as_ref().unwrap();
    let num_atoms = mfdata.ctab.n_atoms;
    for a in atoms.iter().take(num_atoms as usize) {
        max_x = fmax(a.fx, max_x);
        min_x = fmin(a.fx, min_x);
        max_y = fmax(a.fy, max_y);
        min_y = fmin(a.fy, min_y);
        max_z = fmax(a.fz, max_z);
        min_z = fmin(a.fz, min_z);
    }
    let mut num_bonds = 0;
    let mut average_bond_length = 0.0;
    if let Some(bonds) = mfdata.ctab.bonds.as_ref() {
        for b in bonds.iter().take(mfdata.ctab.n_bonds.max(0) as usize) {
            let a1 = i32::from(b.atnum1) - 1;
            let a2 = i32::from(b.atnum2) - 1;
            if a1 < 0 || a1 >= num_atoms || a2 < 0 || a2 >= num_atoms || a1 == a2 {
                *err |= 1;
                treat_err(err, 0, p_str_err, "Bond to nonexistent atom");
                continue;
            }
            let (p1, p2) = (&atoms[a1 as usize], &atoms[a2 as usize]);
            let dx = p1.fx - p2.fx;
            let dy = p1.fy - p2.fy;
            let dz = p1.fz - p2.fz;
            average_bond_length += (dx * dx + dy * dy + dz * dz).sqrt();
            num_bonds += 1;
        }
    }
    let coef = |max: f64, min: f64| {
        if max - min <= small_coeff * (max.abs() + min.abs()) {
            0.0
        } else {
            (MAX_C - MIN_C) / (max - min)
        }
    };
    let x_coeff = coef(max_x, min_x);
    let y_coeff = coef(max_y, min_y);
    let z_coeff = coef(max_z, min_z);
    let num_dimensions = if (x_coeff > macheps || y_coeff > macheps) && z_coeff.abs() < macheps {
        2
    } else if z_coeff.abs() > macheps {
        3
    } else {
        0
    };
    if !find_norm_factors {
        return finish(num_dimensions, (min_x, min_y, min_z), coeff);
    }
    coeff = match num_dimensions {
        2 => {
            if x_coeff > macheps && y_coeff > macheps {
                fmin(x_coeff, y_coeff)
            } else if x_coeff > macheps {
                x_coeff
            } else if y_coeff > macheps {
                y_coeff
            } else {
                1.0
            }
        }
        3 => {
            if x_coeff > macheps && y_coeff > macheps {
                fmin(fmin(x_coeff, y_coeff), z_coeff)
            } else if x_coeff > macheps {
                fmin(x_coeff, z_coeff)
            } else if y_coeff > macheps {
                fmin(y_coeff, z_coeff)
            } else {
                z_coeff
            }
        }
        _ => 0.0,
    };
    if num_bonds > 0 {
        average_bond_length /= f64::from(num_bonds);
        if average_bond_length * coeff > MAX_AVE {
            coeff = MAX_AVE / average_bond_length;
        } else if average_bond_length * coeff < macheps {
            coeff = 1.0;
        } else if average_bond_length * coeff < MIN_AVE {
            coeff = MIN_AVE / average_bond_length;
        }
    }
    finish(num_dimensions, (min_x, min_y, min_z), coeff)
}

/// `TREAT_ERR`
fn treat_err(err: &mut i32, new_err: i32, p_str_err: &mut String, msg: &str) {
    if *err == 0 && new_err != 0 {
        *err = new_err;
    }
    add_error_message(p_str_err, msg);
}

/// the line `inchi_fgetsLf` reads (without the zero terminator), `None` at
/// the end of the input
fn get_line(inp: &mut InchiIostream) -> Option<Vec<u8>> {
    inchi_fgets_lf(MOL_FMT_INPLINELEN, inp).map(|mut l| {
        l.truncate(c_str(&l).len());
        l
    })
}

/// `MolfileReadHeaderLines`
pub fn molfile_read_header_lines(
    hdr: &mut MolFmtHeaderBlock,
    inp: &mut InchiIostream,
    _p_str_err: &mut String,
) -> i32 {
    let Some(mut line) = get_line(inp) else {
        return 1;
    };
    remove_one_lf(&mut line);
    let mut pos = 0;
    molfile_read_field_str(&mut hdr.molname, 200, &line, &mut pos);

    let Some(mut line) = get_line(inp) else {
        return 3;
    };
    remove_one_lf(&mut line);
    let mut pos = 0;
    molfile_read_field_str(&mut hdr.user_initls, 2, &line, &mut pos);
    molfile_read_field_str(&mut hdr.prog_name, 8, &line, &mut pos);
    for f in [
        &mut hdr.month,
        &mut hdr.day,
        &mut hdr.year,
        &mut hdr.hour,
        &mut hdr.minute,
    ] {
        *f = molfile_read_field_int(IntKind::Char, 2, &line, &mut pos).0 as i8;
    }
    molfile_read_field_str(&mut hdr.dim_code, 2, &line, &mut pos);
    hdr.scaling_factor1 = molfile_read_field_int(IntKind::Short, 2, &line, &mut pos).0 as i16;
    hdr.scaling_factor2 = molfile_read_field_double(10, &line, &mut pos).0;
    hdr.energy = molfile_read_field_double(12, &line, &mut pos).0;
    hdr.internal_regno = molfile_read_field_int(IntKind::Long, 6, &line, &mut pos).0;
    let mut pos = 0;
    molfile_read_field_str(&mut hdr.line2, 200, &line, &mut pos);

    let Some(mut line) = get_line(inp) else {
        return 7;
    };
    remove_one_lf(&mut line);
    let mut pos = 0;
    molfile_read_field_str(&mut hdr.comment, 80, &line, &mut pos);
    0
}

/// stores a `short` read by `MolfileReadField` into the low half of an `int`
/// (little-endian host: the upper half keeps its value)
fn store_short_in_int(dst: &mut i32, value: i64) {
    *dst = (*dst & !0xFFFF) | (value as u16 as i32);
}

/// `MolfileReadCountsLine`
pub fn molfile_read_counts_line(
    ctab: &mut MolFmtCtab,
    inp: &mut InchiIostream,
    p_str_err: &mut String,
) -> i32 {
    let mut err = 0;
    let Some(mut line) = get_line(inp) else {
        treat_err(&mut err, 1, p_str_err, "Cannot read counts line");
        return err;
    };
    remove_one_lf(&mut line);
    if at(&line, MOL_FMT_MAXLINELEN) != 0 {
        treat_err(&mut err, 0, p_str_err, "Too long counts line");
    }
    let mut pos = 0;
    let ok = (|| {
        let (v, r) = molfile_read_field_int(IntKind::Short, 3, &line, &mut pos);
        store_short_in_int(&mut ctab.n_atoms, v);
        if r < 0 {
            return false;
        }
        let (v, r) = molfile_read_field_int(IntKind::Short, 3, &line, &mut pos);
        store_short_in_int(&mut ctab.n_bonds, v);
        if r < 0 {
            return false;
        }
        molfile_read_field_jump(3, &line, &mut pos);
        molfile_read_field_jump(3, &line, &mut pos);
        let (v, r) = molfile_read_field_int(IntKind::Char, 3, &line, &mut pos);
        ctab.chiral_flag = v as i8;
        if r < 0 {
            return false;
        }
        let (v, r) = molfile_read_field_int(IntKind::Short, 3, &line, &mut pos);
        ctab.n_stext_entries = v as i16;
        if r < 0 {
            return false;
        }
        for _ in 0..4 {
            molfile_read_field_jump(3, &line, &mut pos);
        }
        let (v, r) = molfile_read_field_int(IntKind::Short, 3, &line, &mut pos);
        ctab.n_property_lines = v as i16;
        r >= 0
    })();
    if !ok {
        err = 3;
        treat_err(&mut err, 3, p_str_err, "Cannot interpret counts line:");
        dotify_non_printable_chars(&mut line);
        add_error_message(p_str_err, &String::from_utf8_lossy(&line));
        return err;
    }
    molfile_read_field_str(&mut ctab.version_string, 6, &line, &mut pos);
    ctab.v3000 = if ctab.version_string == b"V3000" {
        Some(Box::default())
    } else {
        None
    };
    ctab.sgroups = Vec::new();
    ctab.sgroups_freed = false;
    err
}

/// `MolfileReadAtomsBlock`
pub fn molfile_read_atoms_block(
    ctab: &mut MolFmtCtab,
    inp: &mut InchiIostream,
    mut err: i32,
    p_str_err: &mut String,
) -> i32 {
    const CHARGE_VAL: [i32; 8] = [0, 3, 2, 1, b'R' as i32, -1, -2, -3];
    for i in 0..ctab.n_atoms.max(0) as usize {
        let Some(mut line) = get_line(inp) else {
            if err == 0 {
                treat_err(&mut err, 2, p_str_err, "Cannot read atom block line");
            }
            break;
        };
        remove_one_lf(&mut line);
        if at(&line, MOL_FMT_MAXLINELEN) != 0 {
            treat_err(&mut err, 0, p_str_err, "Too long atom block line");
        }
        if err != 0 {
            if line == SD_FMT_END_OF_DATA {
                err = -err.abs();
                break;
            }
            continue;
        }
        if let Some(coords) = ctab.coords.as_mut() {
            mystrncpy(&mut coords[i], &line, 31);
        }
        let Some(atoms) = ctab.atoms.as_mut() else {
            continue;
        };
        let a = &mut atoms[i];
        let mut pos = 0;
        let first_half = (|| {
            let (v, r) = molfile_read_field_double(10, &line, &mut pos);
            a.fx = v;
            if r < 0 {
                return false;
            }
            let (v, r) = molfile_read_field_double(10, &line, &mut pos);
            a.fy = v;
            if r < 0 {
                return false;
            }
            let (v, r) = molfile_read_field_double(10, &line, &mut pos);
            a.fz = v;
            if r < 0 {
                return false;
            }
            molfile_read_field_jump(1, &line, &mut pos);
            if molfile_read_field_str(&mut a.symbol, 3, &line, &mut pos) == 0 {
                return false;
            }
            for (w, f) in [
                (2, &mut a.mass_difference),
                (3, &mut a.charge),
                (3, &mut a.stereo_parity),
            ] {
                let (v, r) = molfile_read_field_int(IntKind::Char, w, &line, &mut pos);
                *f = v as i8;
                if r < 0 {
                    return false;
                }
            }
            molfile_read_field_jump(3, &line, &mut pos);
            molfile_read_field_jump(3, &line, &mut pos);
            let (v, r) = molfile_read_field_int(IntKind::Char, 3, &line, &mut pos);
            a.valence = v as i8;
            r >= 0
        })();
        if !first_half {
            err = 4;
            treat_err(&mut err, 4, p_str_err, "Cannot interpret atom block line:");
            dotify_non_printable_chars(&mut line);
            add_error_message(p_str_err, &String::from_utf8_lossy(&line));
            if line == SD_FMT_END_OF_DATA {
                err = -err.abs();
                break;
            }
            continue;
        }
        if a.symbol.len() == 2 && a.symbol[1].is_ascii_uppercase() {
            a.symbol[1] = a.symbol[1].to_ascii_lowercase();
        }
        let chg = i32::from(a.charge);
        if !(0..CHARGE_VAL.len() as i32).contains(&chg) {
            a.charge = (4 - chg) as i8;
            a.radical = 0;
        } else if CHARGE_VAL[chg as usize] == i32::from(b'R') {
            a.charge = 0;
            a.radical = RADICAL_DOUBLET;
        } else {
            a.charge = CHARGE_VAL[chg as usize] as i8;
            a.radical = 0;
        }
        for _ in 0..6 {
            molfile_read_field_jump(3, &line, &mut pos);
        }
    }
    err
}

/// `MolfileReadBondsBlock`
pub fn molfile_read_bonds_block(
    ctab: &mut MolFmtCtab,
    inp: &mut InchiIostream,
    mut err: i32,
    p_str_err: &mut String,
) -> i32 {
    for i in 0..ctab.n_bonds.max(0) as usize {
        let Some(mut line) = get_line(inp) else {
            if err == 0 {
                treat_err(&mut err, 2, p_str_err, "Cannot read bond block line");
            }
            break;
        };
        remove_one_lf(&mut line);
        if at(&line, MOL_FMT_MAXLINELEN) != 0 {
            err = if err != 0 { err } else { 3 };
        }
        if err != 0 {
            if line == SD_FMT_END_OF_DATA {
                err = -err.abs();
                break;
            }
            continue;
        }
        let Some(bonds) = ctab.bonds.as_mut() else {
            continue;
        };
        let b = &mut bonds[i];
        let mut pos = 0;
        let ok = (|| {
            let (v, r) = molfile_read_field_int(IntKind::Short, 3, &line, &mut pos);
            b.atnum1 = v as i16;
            if r < 0 {
                return false;
            }
            let (v, r) = molfile_read_field_int(IntKind::Short, 3, &line, &mut pos);
            b.atnum2 = v as i16;
            if r < 0 {
                return false;
            }
            let (v, r) = molfile_read_field_int(IntKind::Char, 3, &line, &mut pos);
            b.bond_type = v as i8;
            if r < 0 {
                return false;
            }
            let (v, r) = molfile_read_field_int(IntKind::Char, 3, &line, &mut pos);
            b.bond_stereo = v as i8;
            if r < 0 {
                return false;
            }
            molfile_read_field_jump(3, &line, &mut pos);
            molfile_read_field_jump(3, &line, &mut pos);
            true
        })();
        if !ok {
            if err == 0 {
                treat_err(&mut err, 4, p_str_err, "Cannot interpret bond block line:");
                dotify_non_printable_chars(&mut line);
                add_error_message(p_str_err, &String::from_utf8_lossy(&line));
            }
            if line == SD_FMT_END_OF_DATA {
                err = -err.abs();
                break;
            }
        }
    }
    err
}

/// `MolfileReadSTextBlock`: passes by the STEXT entries
pub fn molfile_read_stext_block(
    ctab: &MolFmtCtab,
    inp: &mut InchiIostream,
    mut err: i32,
    p_str_err: &mut String,
) -> i32 {
    // `S_SHORT i` is compared with `2 * n_stext_entries` as an int
    let mut i: i16 = 0;
    while i32::from(i) < 2 * i32::from(ctab.n_stext_entries) {
        if get_line(inp).is_none() {
            if err == 0 {
                treat_err(&mut err, 2, p_str_err, "Cannot read STEXT block line");
            }
            break;
        }
        i = i.wrapping_add(1);
    }
    err
}

pub const MOL_FMT_M_STY_NON: i32 = 0;
pub const MOL_FMT_M_STY_SRU: i32 = 1;
pub const MOL_FMT_M_STY_MON: i32 = 2;
pub const MOL_FMT_M_STY_COP: i32 = 3;
pub const MOL_FMT_M_STY_MOD: i32 = 4;
pub const MOL_FMT_M_STY_CRO: i32 = 5;
pub const MOL_FMT_M_STY_MER: i32 = 6;
pub const MOL_FMT_M_SST_NON: i32 = 0;
pub const MOL_FMT_M_SST_ALT: i32 = 1;
pub const MOL_FMT_M_SST_RAN: i32 = 2;
pub const MOL_FMT_M_SST_BLK: i32 = 3;
pub const MOL_FMT_M_CONN_NON: i32 = 0;
pub const MOL_FMT_M_CONN_HT: i32 = 1;
pub const MOL_FMT_M_CONN_HH: i32 = 2;
pub const MOL_FMT_M_CONN_EU: i32 = 3;
/// `ZERO_ATW_DIFF`: marks the mass of the most abundant isotope
pub const ZERO_ATW_DIFF: i32 = 127;
pub const MAXVAL: i32 = 20;

/// `MolFmtSgroups_GetIndexBySgroupId`
pub fn molfmt_sgroups_get_index_by_sgroup_id(id: i32, sgroups: &[MolFmtSgroup]) -> Option<usize> {
    sgroups.iter().position(|g| g.id == id)
}

/// `MolFmtSgroups_Append`: false (C: -1) once the container was freed
fn sgroups_append(ctab: &mut MolFmtCtab, id: i32, type_: i32) -> bool {
    if ctab.sgroups_freed {
        return false;
    }
    ctab.sgroups.push(MolFmtSgroup {
        id,
        type_,
        ..Default::default()
    });
    true
}

/// `MolFmtSgroups_Free`
fn sgroups_free(ctab: &mut MolFmtCtab) {
    ctab.sgroups = Vec::new();
    ctab.sgroups_freed = true;
}

fn read_short(len: i32, line: &[u8], pos: &mut usize) -> (i16, i32) {
    let (v, r) = molfile_read_field_int(IntKind::Short, len, line, pos);
    (v as i16, r)
}

/// `MolfileReadSgroupOfPolymer`
pub fn molfile_read_sgroup_of_polymer(
    ctab: &mut MolFmtCtab,
    line: &[u8],
    sz_type: &[u8],
    pos: &mut usize,
    mut err: i32,
) -> i32 {
    let num_atoms = ctab.n_atoms as i16;
    let num_bonds = ctab.n_bonds as i16;
    let mut fail = false;
    let ty = sz_type;
    let multi = |pos: &mut usize| {
        let (n, r) = read_short(3, line, pos);
        (n, r > 0 && (1..=8).contains(&n))
    };
    'exit: {
        if ty == b"STY" && {
            let (n, ok) = multi(pos);
            if ok {
                for _ in 0..n {
                    let (num, r1) = read_short(0, line, pos);
                    fail = r1 < 0;
                    let mut stmp = Vec::new();
                    if !fail {
                        molfile_read_field_str(&mut stmp, 4, line, pos);
                        lrtrim(&mut stmp);
                        let id = i32::from(num);
                        let typ = match stmp.as_slice() {
                            b"SRU" => Some(MOL_FMT_M_STY_SRU),
                            b"MON" => Some(MOL_FMT_M_STY_MON),
                            b"COP" => Some(MOL_FMT_M_STY_COP),
                            b"MOD" => Some(MOL_FMT_M_STY_MOD),
                            b"CRO" => Some(MOL_FMT_M_STY_CRO),
                            b"MER" => Some(MOL_FMT_M_STY_MER),
                            _ => None,
                        };
                        match typ {
                            None => fail = true,
                            Some(typ) => {
                                let index = match molfmt_sgroups_get_index_by_sgroup_id(
                                    id,
                                    &ctab.sgroups,
                                ) {
                                    None => {
                                        if sgroups_append(ctab, id, typ) {
                                            Some(ctab.sgroups.len() - 1)
                                        } else {
                                            fail = true;
                                            None
                                        }
                                    }
                                    Some(i) => {
                                        ctab.sgroups[i].type_ = typ;
                                        Some(i)
                                    }
                                };
                                if let (false, Some(i)) = (fail, index) {
                                    ctab.sgroups[i].xbr1 = [-777777.777; 4];
                                    ctab.sgroups[i].xbr2 = [-777777.777; 4];
                                    ctab.sgroups[i].smt.clear();
                                }
                            }
                        }
                    }
                    if fail {
                        err = 5;
                        sgroups_free(ctab);
                        return err;
                    }
                }
            }
            ok
        } {
        } else if ty == b"SST" && {
            let (n, ok) = multi(pos);
            if ok {
                for _ in 0..n {
                    let (num, r1) = read_short(0, line, pos);
                    fail = r1 < 0;
                    let mut stmp = Vec::new();
                    if !fail {
                        molfile_read_field_str(&mut stmp, 4, line, pos);
                        let index =
                            molfmt_sgroups_get_index_by_sgroup_id(i32::from(num), &ctab.sgroups);
                        if index.is_none() {
                            fail = true;
                        }
                        if let (false, Some(i)) = (fail, index) {
                            ctab.sgroups[i].subtype = MOL_FMT_M_SST_NON;
                            lrtrim(&mut stmp);
                            match stmp.as_slice() {
                                b"ALT" => ctab.sgroups[i].subtype = MOL_FMT_M_SST_ALT,
                                b"RAN" => ctab.sgroups[i].subtype = MOL_FMT_M_SST_RAN,
                                b"BLO" | b"BLK" => ctab.sgroups[i].subtype = MOL_FMT_M_SST_BLK,
                                _ => {
                                    fail = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                if fail {
                    err = 6;
                    sgroups_free(ctab);
                    return err;
                }
            }
            ok
        } {
        } else if ty == b"SLB" && {
            let (n, ok) = multi(pos);
            if ok {
                for _ in 0..n {
                    let (num, r1) = read_short(0, line, pos);
                    let mut tmp = 0;
                    fail = r1 < 0 || {
                        let (t, r2) = read_short(0, line, pos);
                        tmp = t;
                        r2 < 0
                    };
                    if !fail {
                        match molfmt_sgroups_get_index_by_sgroup_id(i32::from(num), &ctab.sgroups) {
                            None => fail = true,
                            Some(i) => ctab.sgroups[i].label = i32::from(tmp),
                        }
                    }
                }
                if fail {
                    err = 7;
                    sgroups_free(ctab);
                    return err;
                }
            }
            ok
        } {
        } else if ty == b"SCN" && {
            let (n, ok) = multi(pos);
            if ok {
                for _ in 0..n {
                    let (num, r1) = read_short(0, line, pos);
                    fail = r1 < 0;
                    let mut stmp = Vec::new();
                    if !fail {
                        molfile_read_field_str(&mut stmp, 4, line, pos);
                        match molfmt_sgroups_get_index_by_sgroup_id(i32::from(num), &ctab.sgroups) {
                            None => fail = true,
                            Some(i) => {
                                ctab.sgroups[i].conn = MOL_FMT_M_CONN_NON;
                                lrtrim(&mut stmp);
                                match stmp.as_slice() {
                                    b"HT" => ctab.sgroups[i].conn = MOL_FMT_M_CONN_HT,
                                    b"HH" => ctab.sgroups[i].conn = MOL_FMT_M_CONN_HH,
                                    b"EU" => ctab.sgroups[i].conn = MOL_FMT_M_CONN_EU,
                                    _ => fail = true,
                                }
                            }
                        }
                    }
                    if fail {
                        err = 8;
                        sgroups_free(ctab);
                        return err;
                    }
                }
            }
            ok
        } {
        } else if ty == b"SAL" || ty == b"SBL" {
            let is_atoms = ty == b"SAL";
            let (sg_num, r1) = read_short(4, line, pos);
            let num_entries;
            if r1 > 0 && {
                let (n, r2) = read_short(3, line, pos);
                num_entries = n;
                r2 > 0
            } {
                let index = molfmt_sgroups_get_index_by_sgroup_id(i32::from(sg_num), &ctab.sgroups);
                if index.is_none() {
                    fail = true;
                }
                let mut items = Vec::new();
                if !fail {
                    let limit = if is_atoms { num_atoms } else { num_bonds };
                    for _ in 0..num_entries {
                        let (v, r) = read_short(0, line, pos);
                        if r < 0 || v <= 0 || v > limit {
                            fail = true;
                            break;
                        }
                        items.push(i32::from(v));
                    }
                }
                if !fail {
                    let g = &mut ctab.sgroups[index.unwrap()];
                    if is_atoms {
                        g.alist.extend(items);
                    } else {
                        g.blist.extend(items);
                    }
                }
            } else {
                fail = true;
            }
            if fail {
                err = if is_atoms { 9 } else { 10 };
                sgroups_free(ctab);
                return err;
            }
        } else if ty == b"SDI" {
            let (sg_num, r1) = read_short(4, line, pos);
            if r1 > 0 && {
                let (n, r2) = read_short(3, line, pos);
                if r2 > 0 {
                    let index =
                        molfmt_sgroups_get_index_by_sgroup_id(i32::from(sg_num), &ctab.sgroups);
                    if index.is_none() || n != 4 {
                        fail = true;
                    }
                    let mut x = [0f64; 4];
                    if !fail {
                        for xj in x.iter_mut().take(n as usize) {
                            let (v, r) = molfile_read_field_double(0, line, pos);
                            *xj = v;
                            if r < 0 {
                                fail = true;
                                break;
                            }
                        }
                    }
                    if !fail {
                        let g = &mut ctab.sgroups[index.unwrap()];
                        if (-g.xbr1[0].abs() + 777777.777).abs() < 1.0e-7 {
                            g.xbr1 = x;
                        } else {
                            g.xbr2 = x;
                        }
                    }
                }
                r2 > 0
            } {
            } else {
                fail = true;
            }
            if fail {
                err = 11;
                sgroups_free(ctab);
                return err;
            }
        } else if ty == b"SMT" {
            let mut index = None;
            let (sg_num, r1) = read_short(4, line, pos);
            let mut stmplong = Vec::new();
            if r1 > 0 && molfile_read_field_str(&mut stmplong, 80, line, pos) > 0 {
                index = molfmt_sgroups_get_index_by_sgroup_id(i32::from(sg_num), &ctab.sgroups);
            }
            match index {
                None => fail = true,
                Some(i) => {
                    lrtrim(&mut stmplong);
                    ctab.sgroups[i].smt = stmplong;
                }
            }
            if fail {
                err = 11;
                sgroups_free(ctab);
                return err;
            }
        }
        break 'exit;
    }
    let _ = err;
    0
}

/// `MolfileTreatPseudoElementAtoms`: returns the number of pseudo-atoms
pub fn molfile_treat_pseudo_element_atoms(
    ctab: &mut MolFmtCtab,
    pseudos_allowed: bool,
    err: &mut i32,
    p_str_err: &mut String,
) -> i32 {
    let mut nzz = 0;
    let Some(atoms) = ctab.atoms.as_mut() else {
        return 0;
    };
    for a in atoms.iter_mut().take(ctab.n_atoms.max(0) as usize) {
        if a.symbol == b"Zy" {
            treat_err(err, 76, p_str_err, "Invalid element(s):");
            treat_err(err, 76, p_str_err, "Zy");
        }
        let is_star = a.symbol == b"*";
        let is_zz = !is_star && a.symbol == b"Zz";
        if is_star || is_zz {
            nzz += 1;
            if !pseudos_allowed {
                treat_err(err, 76, p_str_err, "Invalid element(s):");
                treat_err(err, 76, p_str_err, &String::from_utf8_lossy(&a.symbol));
            } else if is_star {
                a.symbol = b"Zz".to_vec();
            }
        }
    }
    nzz
}

/// `MolfileReadPropBlock`
pub fn molfile_read_prop_block(
    ctab: &mut MolFmtCtab,
    hdr: &mut MolFmtHeaderBlock,
    inp: &mut InchiIostream,
    treat_polymers: bool,
    mut err: i32,
    p_str_err: &mut String,
    no_warnings: bool,
) -> i32 {
    const ALIAS: bool = true;
    let mut multi_line_mode = !ALIAS;
    let mut atom_number: i32 = 0;
    let mut skip_lines: i16 = 0;
    let num_atoms = ctab.n_atoms as i16;
    let mut charge_encountered = 0;
    let mut radical_encountered = 0;
    let mut isotope_encountered = 0;
    let mut polymer_occurred = false;
    let mut i: i16 = 0;
    'fin: {
        while !ctab.version_string.is_empty() || i < ctab.n_property_lines {
            i = i.wrapping_add(1);
            let Some(mut line) = get_line(inp) else {
                if err == 0 {
                    treat_err(&mut err, 2, p_str_err, "Cannot read properties block line");
                }
                break 'fin;
            };
            remove_one_lf(&mut line);
            if at(&line, MOL_FMT_MAXLINELEN) != 0 {
                treat_err(&mut err, 3, p_str_err, "Too long properties block line");
                continue;
            }
            if skip_lines > 0 {
                skip_lines -= 1;
                continue;
            }
            if multi_line_mode == ALIAS && atom_number != 0 {
                multi_line_mode = !ALIAS;
                let len = normalize_string(&mut line);
                if len <= 0 {
                    atom_number = 0;
                    continue;
                }
                if len < 6 {
                    let (n_rad, n_charge) = {
                        let (r, c) = extract_charges_and_radicals(&mut line);
                        (r, c)
                    };
                    let atom = &mut ctab.atoms.as_mut().unwrap()[(atom_number - 1) as usize];
                    atom.charge = n_charge as i8;
                    atom.radical = n_rad as i8;
                    if len == 1 && at(&line, 0) == b'D' {
                        line[0] = b'H';
                        atom.mass_difference = 1;
                    } else if len == 1 && at(&line, 0) == b'T' {
                        line[0] = b'H';
                        atom.mass_difference = 2;
                    } else {
                        atom.mass_difference = 0;
                    }
                    atom.symbol = if line.len() < 6 {
                        line.clone()
                    } else {
                        b"???".to_vec()
                    };
                    atom.atom_aliased_flag = atom.atom_aliased_flag.wrapping_add(1);
                }
                skip_lines = 0;
                atom_number = 0;
                continue;
            }
            let mut pos = 0;
            let mut char_m = Vec::new();
            let mut blank = Vec::new();
            let mut sz_type = Vec::new();
            if molfile_read_field_str(&mut char_m, 1, &line, &mut pos) != 1
                || molfile_read_field_str(&mut blank, 2, &line, &mut pos) != 0
                || molfile_read_field_str(&mut sz_type, 3, &line, &mut pos) <= 0
            {
                if line == SD_FMT_END_OF_DATA {
                    err = if err != 0 { -err.abs() } else { -4 };
                    break;
                }
                continue;
            }
            match char_m[0] {
                b'V' => {
                    skip_lines = 0;
                    continue;
                }
                b'G' => {
                    skip_lines = 1;
                    continue;
                }
                b'A' => {
                    atom_number = c_strtol10(&sz_type).0 as i32;
                    if ctab.atoms.is_some() && 0 < atom_number && atom_number <= ctab.n_atoms {
                        multi_line_mode = ALIAS;
                    } else {
                        atom_number = 0;
                        skip_lines = 1;
                    }
                    continue;
                }
                _ => {}
            }
            if char_m[0] == b'S' && sz_type == b"SKP" {
                let (v, r) = read_short(3, &line, &mut pos);
                skip_lines = if r <= 0 { 0 } else { v };
                continue;
            }
            if char_m[0] != b'M' {
                continue;
            }
            if sz_type == b"REG" {
                while at(&line, pos) == b' ' {
                    pos += 1;
                }
                let mut len = 0;
                while at(&line, pos + len) != 0 && at(&line, pos + len) != b' ' {
                    len += 1;
                }
                let len = len.min(MOL_FMT_MAX_VALUE_LEN as usize) as i32;
                hdr.internal_regno = molfile_read_field_int(IntKind::Long, len, &line, &mut pos).0;
                continue;
            }
            if sz_type == b"END" {
                if !ctab.version_string.is_empty() {
                    break;
                }
                continue;
            }
            if ctab.atoms.is_none() {
                continue;
            }
            let multi = |pos: &mut usize| {
                let (n, r) = read_short(3, &line, pos);
                (n, r > 0 && (1..=8).contains(&n))
            };
            if sz_type == b"CHG" {
                let (num_entries, ok) = multi(&mut pos);
                if ok {
                    if charge_encountered == 0 && radical_encountered == 0 {
                        charge_encountered = -1;
                    }
                    let mut atoms = [0i16; 8];
                    let mut charges = [0i16; 8];
                    let mut bad = false;
                    for j in 0..num_entries as usize {
                        let (a, r1) = read_short(0, &line, &mut pos);
                        atoms[j] = a;
                        let mut r2 = 0;
                        if r1 >= 0 {
                            let (c, r) = read_short(0, &line, &mut pos);
                            charges[j] = c;
                            r2 = r;
                        }
                        if r1 < 0
                            || r2 < 0
                            || atoms[j] <= 0
                            || atoms[j] > num_atoms
                            || !(-15..=15).contains(&charges[j])
                        {
                            bad = true;
                            break;
                        }
                    }
                    let list = ctab.atoms.as_mut().unwrap();
                    if bad {
                        treat_err(&mut err, 0, p_str_err, "Charge not recognized:");
                        dotify_non_printable_chars(&mut line);
                        add_error_message(p_str_err, &String::from_utf8_lossy(&line));
                        continue;
                    }
                    if charge_encountered == -1 {
                        for a in list.iter_mut().take(num_atoms.max(0) as usize) {
                            if a.atom_aliased_flag == 0 {
                                a.charge = 0;
                                a.radical = 0;
                            }
                        }
                        charge_encountered = 1;
                    }
                    for j in 0..num_entries as usize {
                        let a = &mut list[(atoms[j] - 1) as usize];
                        if a.atom_aliased_flag == 0 {
                            a.charge = charges[j] as i8;
                        }
                    }
                    continue;
                }
            }
            if sz_type == b"RAD" {
                let (num_entries, ok) = multi(&mut pos);
                if ok {
                    if charge_encountered == 0 && radical_encountered == 0 {
                        radical_encountered = -1;
                    }
                    let mut atoms = [0i16; 8];
                    let mut rads = [0i16; 8];
                    let mut bad = false;
                    for j in 0..num_entries as usize {
                        let (a, r1) = read_short(0, &line, &mut pos);
                        atoms[j] = a;
                        let mut r2 = 0;
                        if r1 >= 0 {
                            let (c, r) = read_short(0, &line, &mut pos);
                            rads[j] = c;
                            r2 = r;
                        }
                        if r1 < 0
                            || r2 < 0
                            || atoms[j] <= 0
                            || atoms[j] > num_atoms
                            || !(0..=3).contains(&rads[j])
                        {
                            bad = true;
                            break;
                        }
                    }
                    let list = ctab.atoms.as_mut().unwrap();
                    if bad {
                        treat_err(&mut err, 0, p_str_err, "Radical not recognized:");
                        dotify_non_printable_chars(&mut line);
                        add_error_message(p_str_err, &String::from_utf8_lossy(&line));
                        continue;
                    }
                    if radical_encountered == -1 {
                        for a in list.iter_mut().take(num_atoms.max(0) as usize) {
                            if a.atom_aliased_flag == 0 {
                                a.charge = 0;
                                a.radical = 0;
                            }
                        }
                        radical_encountered = 1;
                    }
                    for j in 0..num_entries as usize {
                        let a = &mut list[(atoms[j] - 1) as usize];
                        if a.atom_aliased_flag == 0 {
                            a.radical = rads[j] as i8;
                        }
                    }
                    continue;
                }
            }
            if sz_type == b"ISO" {
                let (num_entries, ok) = multi(&mut pos);
                if ok {
                    if isotope_encountered == 0 {
                        isotope_encountered = -1;
                    }
                    let mut atoms = [0i16; 8];
                    let mut iso_mass = [0i16; 8];
                    for j in 0..num_entries as usize {
                        let (a, r1) = read_short(0, &line, &mut pos);
                        atoms[j] = a;
                        let mut r2 = 0;
                        if r1 >= 0 {
                            let (c, r) = read_short(0, &line, &mut pos);
                            iso_mass[j] = c;
                            r2 = r;
                        }
                        if r1 < 0 || r2 < 0 || atoms[j] <= 0 || atoms[j] > num_atoms {
                            atoms[j] = -1;
                            treat_err(&mut err, 0, p_str_err, "Isotopic data not recognized:");
                            dotify_non_printable_chars(&mut line);
                            add_error_message(p_str_err, &String::from_utf8_lossy(&line));
                        }
                    }
                    let list = ctab.atoms.as_mut().unwrap();
                    if isotope_encountered == -1 {
                        for a in list.iter_mut().take(num_atoms.max(0) as usize) {
                            a.mass_difference = 0;
                        }
                        isotope_encountered = 1;
                    }
                    for j in 0..num_entries as usize {
                        if atoms[j] <= 0 {
                            continue;
                        }
                        let a = &mut list[(atoms[j] - 1) as usize];
                        if at(&a.symbol, 1) != 0
                            || (at(&a.symbol, 0) != b'D' && at(&a.symbol, 0) != b'T')
                        {
                            let atw = get_atomic_mass(&String::from_utf8_lossy(&a.symbol));
                            let atw_diff = i32::from(iso_mass[j]) - atw;
                            if atw != 0 && atw_diff.abs() < 20 {
                                a.mass_difference = (if atw_diff != 0 {
                                    atw_diff
                                } else {
                                    ZERO_ATW_DIFF
                                }) as i8;
                            }
                        }
                    }
                    continue;
                }
            }
            if matches!(
                sz_type.as_slice(),
                b"STY" | b"SST" | b"SLB" | b"SCN" | b"SAL" | b"SBL" | b"SDI" | b"SMT" | b"SBT"
            ) {
                if !treat_polymers {
                    polymer_occurred = true;
                    continue;
                }
                let result = molfile_read_sgroup_of_polymer(ctab, &line, &sz_type, &mut pos, err);
                if result != 0 {
                    treat_err(
                        &mut err,
                        result,
                        p_str_err,
                        "Could not interpret Molfile polymer data:",
                    );
                    dotify_non_printable_chars(&mut line);
                    add_error_message(p_str_err, &String::from_utf8_lossy(&line));
                    continue;
                }
            }
        }
    }
    if !treat_polymers && polymer_occurred && !no_warnings {
        add_error_message(p_str_err, "Ignore polymer data");
    }
    err
}

/// `MolfileReadDataLines` for V2000 input. `only` is the
/// `OnlyHeaderBlock`/`OnlyCTab` mode (nothing is allocated for atoms and bonds);
/// otherwise the data are returned. The V3000 reader is not ported yet.
/// Returns the data (read-all mode, no error) and the error code (negative
/// when the end of the data was met).
pub fn molfile_read_data_lines(
    inp: &mut InchiIostream,
    only: Option<&mut MolFmtData>,
    get_orig_coord: bool,
    treat_polymers: bool,
    p_str_err: &mut String,
    no_warnings: bool,
) -> (Option<MolFmtData>, i32) {
    let should_read_all = only.is_none();
    let mut local = MolFmtData::default();
    let m: &mut MolFmtData = match only {
        Some(o) => {
            *o = MolFmtData::default();
            o
        }
        None => &mut local,
    };
    let mut data_ended = false;
    let mut retcode;
    let mut prevcode;
    'fin: {
        retcode = molfile_read_header_lines(&mut m.hdr, inp, p_str_err);
        if retcode != 0 {
            retcode += 10;
            break 'fin;
        }
        retcode = molfile_read_counts_line(&mut m.ctab, inp, p_str_err);
        if retcode != 0 {
            retcode += 20;
            break 'fin;
        }
        assert!(m.ctab.v3000.is_none(), "V3000 molfiles are not ported yet");
        if should_read_all {
            let n_alloc = m.ctab.n_atoms.max(1) as usize;
            m.ctab.atoms = Some(vec![MolFmtAtom::default(); n_alloc]);
            if get_orig_coord {
                m.ctab.coords = Some(vec![[0u8; super::ichiprt1::MOL_COORD_LEN]; n_alloc]);
            }
        }
        retcode = molfile_read_atoms_block(&mut m.ctab, inp, retcode, p_str_err);
        if retcode != 0 {
            if retcode < 0 {
                retcode = -retcode;
                data_ended = true;
            }
            retcode += 30;
        }
        if should_read_all && retcode < 30 && !data_ended {
            m.ctab.bonds = Some(vec![MolFmtBond::default(); m.ctab.n_bonds.max(1) as usize]);
        }
        prevcode = retcode;
        if !data_ended {
            retcode = molfile_read_bonds_block(&mut m.ctab, inp, retcode, p_str_err);
            if retcode != 0 {
                if retcode < 0 {
                    retcode = -retcode;
                    data_ended = true;
                }
                retcode = if prevcode != 0 {
                    prevcode
                } else {
                    retcode + 40
                };
            }
        }
        prevcode = retcode;
        if !data_ended {
            retcode = molfile_read_stext_block(&m.ctab, inp, retcode, p_str_err);
            if retcode != 0 {
                retcode = if prevcode != 0 {
                    prevcode
                } else {
                    retcode + 50
                };
            }
        }
        prevcode = retcode;
        if !data_ended {
            let (hdr, ctab) = (&mut m.hdr, &mut m.ctab);
            retcode = molfile_read_prop_block(
                ctab,
                hdr,
                inp,
                treat_polymers,
                retcode,
                p_str_err,
                no_warnings,
            );
            if retcode != 0 {
                if retcode < 0 {
                    retcode = -retcode;
                    data_ended = true;
                }
                retcode = if prevcode != 0 {
                    prevcode
                } else {
                    retcode + 60
                };
            }
        }
        if let Some(atoms) = m.ctab.atoms.as_ref() {
            for a in atoms.iter().take(m.ctab.n_atoms.max(0) as usize) {
                if i32::from(a.valence) > MAXVAL {
                    retcode = 70 + 9;
                    treat_err(&mut retcode, 0, p_str_err, "Too large input atomic valence");
                    break;
                }
            }
            for a in atoms.iter().take(m.ctab.n_atoms.max(0) as usize) {
                let dmass = i32::from(a.mass_difference);
                if (a.symbol == b"H" && !matches!(dmass, 0 | 1 | 2 | 127))
                    || (a.symbol == b"D" && !(-1..=1).contains(&dmass))
                    || (a.symbol == b"T" && !(-2..=0).contains(&dmass))
                {
                    retcode = 70 + 8;
                    treat_err(
                        &mut retcode,
                        0,
                        p_str_err,
                        "Unacceptable isotope of hydrogen",
                    );
                    break;
                }
            }
        }
    }
    let err = if data_ended { -retcode } else { retcode };
    if should_read_all {
        (if retcode != 0 { None } else { Some(local) }, err)
    } else {
        (None, err)
    }
}
