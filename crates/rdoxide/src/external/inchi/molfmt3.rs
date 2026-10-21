//! Port of `mol_fmt3.c`: the V3000 Molfile reader of the InChI library
//! (CTAB head, atoms, bonds, haptic bonds, enhanced-stereo collections, tail).
//!
//! International Chemical Identifier (InChI), Version 1, Software version 1.07.
//! Copyright (c) 2024 IUPAC and InChI Trust. Licensed under the MIT license
//! (see LICENSES/InChI-MIT.txt).
//!
//! As in the C code, a "line pointer" `p` is an `Option<usize>` into the current
//! line (`None` is the NULL pointer) and `line` is the line buffer.

use super::ichi_io::{InchiIostream, inchi_ios_str_getc};
use super::ichierr::add_error_message;
use super::molfmt::{
    MolFmtCtab, SD_FMT_END_OF_DATA, ZERO_ATW_DIFF, c_strtod, c_strtol10, treat_err,
};
use super::util::{
    dotify_non_printable_chars, get_atomic_mass, mystrncpy, read_upto_delim, remove_one_lf,
};

pub const MOL_FMT_V3000_MAXFIELDLEN: i32 = 4096;
pub const MOL_FMT_V3000_STENON: i32 = -1;
pub const MOL_FMT_V3000_STEABS: i32 = 1;
pub const MOL_FMT_V3000_STEREL: i32 = 2;
pub const MOL_FMT_V3000_STERAC: i32 = 3;
const WS: &[u8] = b" \t\n\x0b\x0c\r";
const WS_DIGITS: &[u8] = b"1234567890 \t\n\x0b\x0c\r";
const WS_EQ: &[u8] = b"= \t\n\x0b\x0c\r";

/// the integer kinds of `MolfileV3000ReadField`
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum V3Int {
    Char,
    Short,
    Int,
    Long,
}

fn c_str(s: &[u8]) -> &[u8] {
    &s[..s.iter().position(|&c| c == 0).unwrap_or(s.len())]
}

/// the text from `p` on (what a C `strcmp(p, ..)` sees)
fn from_p(line: &[u8], p: Option<usize>) -> &[u8] {
    match p {
        Some(i) => c_str(&line[i.min(line.len())..]),
        None => b"",
    }
}

/// `MolfileV3000ReadField` for strings: returns the number of bytes read
/// (-1 on error, then `data` is unchanged)
pub fn v3000_read_field_str(line: &[u8], p: &mut Option<usize>, data: &mut Vec<u8>) -> i32 {
    let mut field = Vec::new();
    let nread = read_upto_delim(line, p, &mut field, MOL_FMT_V3000_MAXFIELDLEN, WS);
    if nread != 0 {
        if nread > 0 {
            *data = field[..nread as usize].to_vec();
        }
    } else {
        data.clear();
    }
    nread
}

/// `MolfileV3000ReadField` for the integer types: (value, bytes read)
pub fn v3000_read_field_int(kind: V3Int, line: &[u8], p: &mut Option<usize>) -> (i64, i32) {
    let mut field = Vec::new();
    let mut nread = read_upto_delim(line, p, &mut field, MOL_FMT_V3000_MAXFIELDLEN, WS);
    let (ldata, used) = c_strtol10(&field);
    if used == 0 {
        nread = 0;
    }
    let (lo, hi, strict) = match kind {
        V3Int::Long => (i64::MIN, i64::MAX, true),
        V3Int::Int => (i64::from(i32::MIN), i64::from(i32::MAX), false),
        V3Int::Char => (i64::from(i8::MIN), i64::from(i8::MAX), false),
        V3Int::Short => (i64::from(i16::MIN), i64::from(i16::MAX), false),
    };
    let ok = if strict {
        lo < ldata && ldata < hi
    } else {
        lo <= ldata && ldata <= hi
    };
    if ok { (ldata, nread) } else { (0, -1) }
}

/// `MolfileV3000ReadField` for doubles
pub fn v3000_read_field_double(line: &[u8], p: &mut Option<usize>) -> (f64, i32) {
    let mut field = Vec::new();
    let mut nread = read_upto_delim(line, p, &mut field, MOL_FMT_V3000_MAXFIELDLEN, WS);
    let (d, used) = c_strtod(&field);
    if used == 0 {
        nread = 0;
    }
    if d.is_infinite() {
        (0.0, -1)
    } else {
        (d, nread)
    }
}

/// `MolfileV3000ReadKeyword`
pub fn v3000_read_keyword(line: &[u8], p: &mut Option<usize>, key: &mut Vec<u8>) -> i32 {
    let mut field = Vec::new();
    let nread = read_upto_delim(line, p, &mut field, MOL_FMT_V3000_MAXFIELDLEN, WS_EQ);
    if nread != 0 {
        if nread > 0 {
            *key = field[..nread as usize].to_vec();
        }
        if let Some(i) = *p
            && line.get(i) == Some(&b'=')
        {
            *p = Some(i + 1);
        }
    } else {
        key.clear();
    }
    nread
}

/// `get_actual_atom_number`: the final atom number of the original `index`
/// (-1 for star atoms and unknown indices)
pub fn get_actual_atom_number(index: i32, n: i32, orig: &[i32], fin: &[i32]) -> i32 {
    for i in 0..n.max(0) as usize {
        if orig.get(i) == Some(&index) {
            return fin.get(i).copied().unwrap_or(0);
        }
    }
    -1
}

/// `MolfileV3000Init`
pub fn molfile_v3000_init(ctab: &mut MolFmtCtab) -> i32 {
    let n = ctab.n_atoms;
    let v = ctab.v3000.as_mut().expect("V3000 data");
    v.n_star_atoms = 0;
    v.n_non_star_atoms = 0;
    if n > 0 {
        v.atom_index_orig = vec![-1; n as usize];
        v.atom_index_fin = vec![-1; n as usize];
    } else {
        v.atom_index_orig = Vec::new();
        v.atom_index_fin = Vec::new();
    }
    v.n_haptic_bonds = 0;
    v.haptic_bonds = Vec::new();
    v.n_steabs = 0;
    v.steabs = Vec::new();
    v.n_sterel = 0;
    v.sterel = Vec::new();
    v.n_sterac = 0;
    v.sterac = Vec::new();
    0
}

/// `get_V3000_input_line_to_strbuf`: reads the next V3000 line (joining the
/// lines continued with '-', dropping the "M  V30 " prefixes) into `s`;
/// returns its length or -1
pub fn get_v3000_input_line_to_strbuf(s: &mut Vec<u8>, inp: &mut InchiIostream) -> i32 {
    const PREFIX: &[u8] = b"M  V30 ";
    s.clear();
    let mut old_used = 0;
    loop {
        // inchi_strbuf_addline( buf, inp, crlf2lf = 1, preserve_lf = 0 )
        let mut eof = false;
        loop {
            match inchi_ios_str_getc(inp) {
                None => {
                    eof = true;
                    break;
                }
                Some(c) => {
                    s.push(c);
                    if c == b'\n' {
                        break;
                    }
                }
            }
        }
        if !eof {
            if s.len() > 2 && s[s.len() - 2] == b'\r' {
                let k = s.len() - 2;
                s[k] = b'\n';
                s.pop();
            }
            s.pop();
        }
        if s.len() - old_used < 8 {
            return -1;
        }
        if &s[old_used..old_used + PREFIX.len()] != PREFIX {
            return -1;
        }
        s.drain(old_used..old_used + PREFIX.len());
        if s.last() != Some(&b'-') {
            break;
        }
        s.pop();
        old_used = s.len();
    }
    s.len() as i32
}

/// the `p` / `line` pair that the C functions fill with
/// `get_V3000_input_line_to_strbuf`
#[derive(Default)]
struct LineReader {
    buf: Vec<u8>,
    /// the C `line` pointer was assigned at least once
    have_line: bool,
    line: Vec<u8>,
    p: Option<usize>,
}

impl LineReader {
    fn next(&mut self, inp: &mut InchiIostream) {
        let nc = get_v3000_input_line_to_strbuf(&mut self.buf, inp);
        if nc < 1 {
            self.p = None;
        } else {
            self.p = Some(0);
            self.have_line = true;
        }
        self.line = if self.have_line {
            c_str(&self.buf).to_vec()
        } else {
            Vec::new()
        };
    }

    /// the text from `p` on, or "" for NULL
    fn p_text(&self) -> &[u8] {
        from_p(&self.line, self.p)
    }

    fn rm_lf(&mut self) {
        remove_one_lf(&mut self.line);
    }

    fn line_str(&self) -> String {
        String::from_utf8_lossy(c_str(&self.line)).into_owned()
    }

    fn dotify(&mut self) {
        dotify_non_printable_chars(&mut self.line);
    }
}

/// `MolfileV3000ReadCTABBeginAndCountsLine`
pub fn molfile_v3000_read_ctab_begin_and_counts_line(
    ctab: &mut MolFmtCtab,
    inp: &mut InchiIostream,
    p_str_err: &mut String,
) -> i32 {
    let mut err = 0;
    let mut rd = LineReader::default();
    rd.next(inp);
    if rd.p.is_none() || rd.p_text() != b"BEGIN CTAB" {
        treat_err(&mut err, 1, p_str_err, "Error: No V3000 CTab start marker");
        return err;
    }
    rd.rm_lf();
    ctab.n_atoms = -1;
    ctab.n_bonds = -1;
    ctab.chiral_flag = -1;
    ctab.n_stext_entries = -1;
    ctab.n_property_lines = 999;

    rd.next(inp);
    if rd.p.is_none() {
        treat_err(&mut err, 1, p_str_err, "Cannot read V3000 counts line");
        return err;
    }
    rd.rm_lf();
    let line = rd.line.clone();
    let mut p = rd.p;
    let mut field = Vec::new();
    v3000_read_field_str(&line, &mut p, &mut field);
    if field != b"COUNTS" {
        treat_err(&mut err, 1, p_str_err, "Cannot read V3000 counts line");
        return err;
    }
    let mut failed = 0;
    let v3 = ctab.v3000.as_mut().expect("V3000 data");
    let rd_int = |kind: V3Int, p: &mut Option<usize>| v3000_read_field_int(kind, &line, p);
    let (v, n) = rd_int(V3Int::Int, &mut p);
    ctab.n_atoms = v as i32;
    if n < 0 {
        failed = 2;
    } else {
        let (v, n) = rd_int(V3Int::Int, &mut p);
        ctab.n_bonds = v as i32;
        if n < 0 {
            failed = 1;
        } else {
            let (v, n) = rd_int(V3Int::Int, &mut p);
            v3.n_sgroups = v as i32;
            if n < 0 {
                failed = 1;
            } else {
                let (v, n) = rd_int(V3Int::Int, &mut p);
                v3.n_3d_constraints = v as i32;
                if n < 0 {
                    failed = 1;
                } else {
                    let (v, n) = rd_int(V3Int::Char, &mut p);
                    ctab.chiral_flag = v as i8;
                    if n < 0 {
                        failed = 1;
                    }
                }
            }
        }
    }
    if failed != 0 {
        err = 3;
        if failed == 2 {
            treat_err(
                &mut err,
                3,
                p_str_err,
                "Number of atoms too large. V3000 counts line:",
            );
        } else {
            treat_err(
                &mut err,
                3,
                p_str_err,
                "Cannot interpret V3000 counts line:",
            );
        }
        rd.dotify();
        add_error_message(p_str_err, &rd.line_str());
    }
    err
}

/// `MolfileV3000ReadSGroup`: skips the S-group block (the C loop does not
/// terminate without an "END SGROUP" line, so neither does this)
pub fn molfile_v3000_read_sgroup(inp: &mut InchiIostream) -> i32 {
    let mut rd = LineReader::default();
    loop {
        rd.next(inp);
        rd.rm_lf();
        if rd.p.is_some() && rd.p_text() == b"END SGROUP" {
            return 0;
        }
    }
}

/// `MolfileV3000Read3DBlock`
pub fn molfile_v3000_read_3d_block(
    inp: &mut InchiIostream,
    mut err: i32,
    p_str_err: &mut String,
) -> i32 {
    let mut rd = LineReader::default();
    rd.next(inp);
    rd.rm_lf();
    if rd.p.is_none() || rd.p_text() != b"END OBJ3D" {
        treat_err(&mut err, 1, p_str_err, "Error: No V3000 3DBlock end marker");
    }
    err
}

/// `MolfileV3000ReadStereoCollection`: the list `[-1, nnum, atoms...]`
pub fn molfile_v3000_read_stereo_collection(
    line: &[u8],
    p: &mut Option<usize>,
) -> Option<Vec<i32>> {
    let mut field = Vec::new();
    read_upto_delim(line, p, &mut field, MOL_FMT_V3000_MAXFIELDLEN, WS_DIGITS);
    if field != b"(" {
        return None;
    }
    read_upto_delim(line, p, &mut field, MOL_FMT_V3000_MAXFIELDLEN, WS);
    let (nnum, used) = c_strtol10(&field);
    if used == 0 || nnum < 0 {
        return None;
    }
    let nnum = nnum as i32;
    let mut list = vec![0i32; nnum as usize + 3];
    list[0] = -1;
    list[1] = nnum;
    for i in 2..nnum as usize + 2 {
        let (v, n) = v3000_read_field_int(V3Int::Int, line, p);
        list[i] = v as i32;
        if n < 0 {
            return None;
        }
    }
    Some(list)
}

/// `MolfileV3000ReadHapticBond`: the list `[-1, -1, nnum, atoms...]`
pub fn molfile_v3000_read_haptic_bond(line: &[u8], p: &mut Option<usize>) -> Option<Vec<i32>> {
    let mut field = Vec::new();
    read_upto_delim(line, p, &mut field, MOL_FMT_V3000_MAXFIELDLEN, WS_DIGITS);
    if field != b"(" {
        return None;
    }
    read_upto_delim(line, p, &mut field, MOL_FMT_V3000_MAXFIELDLEN, WS);
    let (nnum, used) = c_strtol10(&field);
    if used == 0 || nnum < 0 {
        return None;
    }
    let nnum = nnum as i32;
    let mut list = vec![0i32; nnum as usize + 3];
    list[0] = -1;
    list[1] = -1;
    list[2] = nnum;
    for i in 3..nnum as usize + 3 {
        let (v, n) = v3000_read_field_int(V3Int::Int, line, p);
        list[i] = v as i32;
        if n < 0 {
            return None;
        }
    }
    let nread = read_upto_delim(line, p, &mut field, MOL_FMT_V3000_MAXFIELDLEN, WS);
    if nread > 0 && field != b"ATTACH=ALL" {
        return None;
    }
    if nread < 0 {
        return None;
    }
    Some(list)
}

/// `MolfileV3000ReadCollections`
pub fn molfile_v3000_read_collections(
    ctab: &mut MolFmtCtab,
    inp: &mut InchiIostream,
    mut err: i32,
    p_str_err: &mut String,
) -> i32 {
    let mut rd = LineReader::default();
    let mut failed = false;
    let mut field: Vec<u8> = Vec::new();
    let mut n_coll: i32 = 0;
    rd.next(inp);
    rd.rm_lf();
    while rd.p.is_some() && rd.p_text() != b"END COLLECTION" {
        let line = rd.line.clone();
        let mut p = rd.p;
        let mut stereo_kind = MOL_FMT_V3000_STENON;
        let nread = read_upto_delim(&line, &mut p, &mut field, MOL_FMT_V3000_MAXFIELDLEN, b"/");
        if nread < 6 {
            failed = true;
            break;
        }
        if field != b"MDLV30" {
            failed = true;
            break;
        }
        read_upto_delim(
            &line,
            &mut p,
            &mut field,
            MOL_FMT_V3000_MAXFIELDLEN,
            WS_DIGITS,
        );
        if field == b"/STEABS" {
            n_coll = 1;
            stereo_kind = MOL_FMT_V3000_STEABS;
        } else if field == b"/STEREL" || field == b"/STERAC" {
            // the number is a `char` written into the low byte of an `int`
            let (v, n) = v3000_read_field_int(V3Int::Char, &line, &mut p);
            n_coll = (n_coll & !0xFF) | (v as i32 & 0xFF);
            if n < 0 {
                failed = true;
                break;
            }
            stereo_kind = if field == b"/STEREL" {
                MOL_FMT_V3000_STEREL
            } else {
                MOL_FMT_V3000_STERAC
            };
        }
        if stereo_kind != MOL_FMT_V3000_STENON {
            let len = v3000_read_keyword(&line, &mut p, &mut field);
            if len > 0 {
                if field == b"ATOMS" {
                    match molfile_v3000_read_stereo_collection(&line, &mut p) {
                        None => failed = true,
                        Some(mut list) => {
                            list[0] = n_coll;
                            let nnum = list[1];
                            let v3 = ctab.v3000.as_mut().unwrap();
                            for k in 2..nnum.max(2) as usize {
                                list[k] = get_actual_atom_number(
                                    list[k],
                                    v3.n_non_star_atoms + v3.n_star_atoms,
                                    &v3.atom_index_orig,
                                    &v3.atom_index_fin,
                                );
                            }
                            match stereo_kind {
                                MOL_FMT_V3000_STEABS => {
                                    v3.steabs.push(list);
                                    v3.n_steabs += 1;
                                }
                                MOL_FMT_V3000_STEREL => {
                                    v3.sterel.push(list);
                                    v3.n_sterel += 1;
                                }
                                _ => {
                                    v3.sterac.push(list);
                                    v3.n_sterac += 1;
                                }
                            }
                        }
                    }
                }
            } else {
                failed = true;
                break;
            }
        }
        rd.next(inp);
        rd.rm_lf();
    }
    if failed {
        rd.next(inp);
        rd.rm_lf();
    }
    if rd.p.is_none() {
        failed = true;
    }
    if failed {
        err = 7;
        treat_err(
            &mut err,
            7,
            p_str_err,
            "Cannot interpret V3000 collection line(s)",
        );
        rd.dotify();
        add_error_message(p_str_err, &rd.line_str());
        return err;
    }
    let v3 = ctab.v3000.as_ref().unwrap();
    if v3.n_steabs != 0 || v3.n_sterel != 0 || v3.n_sterac != 0 {
        add_error_message(p_str_err, "V3000 enhanced stereo read/stored but ignored");
    }
    err
}

/// `sprintf("%10g", v)`
fn fmt_g10(v: f64) -> String {
    let s = if v.is_nan() {
        if v.is_sign_negative() { "-nan" } else { "nan" }.to_string()
    } else if v.is_infinite() {
        if v < 0.0 { "-inf" } else { "inf" }.to_string()
    } else if v == 0.0 {
        if v.is_sign_negative() { "-0" } else { "0" }.to_string()
    } else {
        let e = format!("{:.5e}", v);
        let (mant, exp) = e.split_once('e').unwrap();
        let x: i32 = exp.parse().unwrap();
        if !(-4..6).contains(&x) {
            let mant = if mant.contains('.') {
                mant.trim_end_matches('0').trim_end_matches('.')
            } else {
                mant
            };
            format!("{mant}e{}{:02}", if x < 0 { '-' } else { '+' }, x.abs())
        } else {
            let t = format!("{:.*}", (5 - x) as usize, v);
            if t.contains('.') {
                t.trim_end_matches('0').trim_end_matches('.').to_string()
            } else {
                t
            }
        }
    };
    format!("{s:>10}")
}

/// `MolfileV3000ReadAtomsBlock`
pub fn molfile_v3000_read_atoms_block(
    ctab: &mut MolFmtCtab,
    inp: &mut InchiIostream,
    mut err: i32,
    p_str_err: &mut String,
) -> i32 {
    let mut rd = LineReader::default();
    rd.next(inp);
    if rd.p.is_none() || rd.p_text() != b"BEGIN ATOM" {
        treat_err(
            &mut err,
            1,
            p_str_err,
            "Error: No V3000 Atom block start marker",
        );
        return err;
    }
    rd.rm_lf();
    ctab.v3000.as_mut().unwrap().n_non_star_atoms = 0;
    let n_atoms = ctab.n_atoms;
    for i in 0..n_atoms.max(0) as usize {
        rd.next(inp);
        if rd.p.is_none() {
            if err == 0 {
                treat_err(&mut err, 2, p_str_err, "Cannot read V3000 atom block line");
            }
            break;
        }
        rd.rm_lf();
        if err != 0 {
            if c_str(&rd.line) == SD_FMT_END_OF_DATA {
                err = -err.abs();
                break;
            }
            continue;
        }
        if ctab.atoms.is_none() {
            continue;
        }
        let line = rd.line.clone();
        let mut p = rd.p;
        let mut symbol: Vec<u8> = Vec::new();
        let mut failed = false;
        let (index, n) = v3000_read_field_int(V3Int::Int, &line, &mut p);
        if n < 0 {
            failed = true;
        }
        if !failed && v3000_read_field_str(&line, &mut p, &mut symbol) < 0 {
            failed = true;
        }
        let mut xyz = [0f64; 3];
        for x in xyz.iter_mut() {
            if failed {
                break;
            }
            let (v, n) = v3000_read_field_double(&line, &mut p);
            *x = v;
            if n < 0 {
                failed = true;
            }
        }
        if !failed && v3000_read_field_int(V3Int::Int, &line, &mut p).1 < 0 {
            failed = true;
        }
        if failed {
            err = 4;
            treat_err(
                &mut err,
                4,
                p_str_err,
                "Cannot interpret V3000 atom block line:",
            );
            rd.dotify();
            add_error_message(p_str_err, &rd.line_str());
            if c_str(&rd.line) == SD_FMT_END_OF_DATA {
                err = -err.abs();
                break;
            }
            continue;
        }
        let index = index as i32;
        if let Some(coords) = ctab.coords.as_mut() {
            let s = format!("{}{}{}", fmt_g10(xyz[0]), fmt_g10(xyz[1]), fmt_g10(xyz[2]));
            let b = s.as_bytes();
            coords[i].fill(0);
            let n = b.len().min(coords[i].len() - 1);
            coords[i][..n].copy_from_slice(&b[..n]);
        }
        let v3 = ctab.v3000.as_mut().unwrap();
        if symbol == b"*" {
            v3.atom_index_orig[i] = index;
            v3.atom_index_fin[i] = -1;
            v3.n_star_atoms += 1;
            continue;
        }
        v3.n_non_star_atoms += 1;
        v3.atom_index_orig[i] = index;
        v3.atom_index_fin[i] = v3.n_non_star_atoms;
        let ii = (v3.n_non_star_atoms - 1) as usize;
        let atom = &mut ctab.atoms.as_mut().unwrap()[ii];
        let mut sym = vec![0u8; 6];
        let mut sbuf = symbol.clone();
        sbuf.push(0);
        mystrncpy(&mut sym, &sbuf, 6);
        sym.truncate(c_str(&sym).len());
        atom.symbol = sym;
        if atom.symbol.len() == 2 && atom.symbol[1].is_ascii_uppercase() {
            atom.symbol[1] = atom.symbol[1].to_ascii_lowercase();
        }
        atom.fx = xyz[0];
        atom.fy = xyz[1];
        atom.fz = xyz[2];
        let mut key = Vec::new();
        while p.is_some() && v3000_read_keyword(&line, &mut p, &mut key) > 0 {
            failed = false;
            let itmp = |p: &mut Option<usize>| v3000_read_field_int(V3Int::Int, &line, p);
            match key.as_slice() {
                b"CHG" => {
                    let (v, n) = v3000_read_field_int(V3Int::Char, &line, &mut p);
                    atom.charge = v as i8;
                    failed = n < 0;
                }
                b"RAD" => {
                    let (v, n) = v3000_read_field_int(V3Int::Char, &line, &mut p);
                    atom.radical = v as i8;
                    failed = n < 0;
                }
                b"CFG" => {
                    let (v, n) = v3000_read_field_int(V3Int::Char, &line, &mut p);
                    atom.stereo_parity = v as i8;
                    failed = n < 0;
                }
                b"MASS" => {
                    let (v, n) = v3000_read_field_int(V3Int::Short, &line, &mut p);
                    if n < 0 {
                        failed = true;
                        treat_err(&mut err, 0, p_str_err, "Isotopic data not recognized:");
                        add_error_message(p_str_err, &rd.line_str());
                    } else {
                        let atw = get_atomic_mass(&String::from_utf8_lossy(&atom.symbol));
                        let delta = v as i32 - atw;
                        atom.mass_difference =
                            (if delta != 0 { delta } else { ZERO_ATW_DIFF }) as i8;
                    }
                }
                b"VAL" => {
                    let (v, n) = itmp(&mut p);
                    if n < 0 {
                        failed = true;
                    } else {
                        atom.valence = if v == -1 { 15 } else { v as i8 };
                    }
                }
                b"HCOUNT" | b"STBOX" | b"INVRET" | b"EXACHG" | b"SUBST" | b"UNSAT" | b"RBCNT"
                | b"ATTCHPT" | b"SEQID" => {
                    itmp(&mut p);
                }
                b"RGROUPS" | b"ATTCHORD" | b"CLASS" => {
                    let mut s = Vec::new();
                    v3000_read_field_str(&line, &mut p, &mut s);
                }
                _ => {}
            }
            if failed {
                err = 4;
                treat_err(
                    &mut err,
                    4,
                    p_str_err,
                    "Cannot interpret V3000 atom block key-value pair",
                );
                rd.dotify();
                add_error_message(p_str_err, &rd.line_str());
                if c_str(&rd.line) == SD_FMT_END_OF_DATA {
                    err = -err.abs();
                    break;
                }
                continue;
            }
        }
    }
    let v3 = ctab.v3000.as_ref().unwrap();
    if v3.n_star_atoms != 0 {
        add_error_message(p_str_err, "V3000 star atoms ignored");
        ctab.n_atoms = v3.n_non_star_atoms;
    }
    rd.next(inp);
    if rd.p.is_none() || rd.p_text() != b"END ATOM" {
        treat_err(
            &mut err,
            1,
            p_str_err,
            "Error: No V3000 Atom block end marker",
        );
        return err;
    }
    rd.rm_lf();
    err
}

/// `MolfileV3000ReadBondsBlock`
pub fn molfile_v3000_read_bonds_block(
    ctab: &mut MolFmtCtab,
    inp: &mut InchiIostream,
    mut err: i32,
    p_str_err: &mut String,
) -> i32 {
    if ctab.n_bonds == 0 {
        return 0;
    }
    let mut rd = LineReader::default();
    rd.next(inp);
    if rd.p.is_none() || rd.p_text() != b"BEGIN BOND" {
        treat_err(
            &mut err,
            1,
            p_str_err,
            "Error: No V3000 Bond block start marker",
        );
        return err;
    }
    rd.rm_lf();
    {
        let v3 = ctab.v3000.as_mut().unwrap();
        v3.n_haptic_bonds = 0;
        v3.n_non_haptic_bonds = 0;
    }
    for _ in 0..ctab.n_bonds.max(0) {
        let mut is_haptic = false;
        rd.next(inp);
        if rd.p.is_none() {
            if err == 0 {
                treat_err(&mut err, 2, p_str_err, "Cannot read V3000 bond block line");
            }
            break;
        }
        rd.rm_lf();
        if err != 0 {
            if c_str(&rd.line) == SD_FMT_END_OF_DATA {
                err = -err.abs();
                break;
            }
            continue;
        }
        if ctab.bonds.is_none() {
            continue;
        }
        let line = rd.line.clone();
        let mut p = rd.p;
        let n_orig_at = {
            let v3 = ctab.v3000.as_ref().unwrap();
            v3.n_non_star_atoms + v3.n_star_atoms
        };
        let mut failed = false;
        let mut bond_type: i8 = 0;
        let mut stereo: i8 = 0;
        let mut atnum1: i16 = -1;
        let mut atnum2: i16 = -1;
        if v3000_read_field_int(V3Int::Int, &line, &mut p).1 < 0 {
            failed = true;
        } else {
            let (v, n) = v3000_read_field_int(V3Int::Char, &line, &mut p);
            bond_type = v as i8;
            if n < 0 {
                failed = true;
            } else {
                let (v, n) = v3000_read_field_int(V3Int::Short, &line, &mut p);
                atnum1 = v as i16;
                if n < 0 {
                    failed = true;
                } else {
                    let (v, n) = v3000_read_field_int(V3Int::Short, &line, &mut p);
                    atnum2 = v as i16;
                    if n < 0 {
                        failed = true;
                    }
                }
            }
        }
        {
            let v3 = ctab.v3000.as_ref().unwrap();
            atnum1 = get_actual_atom_number(
                i32::from(atnum1),
                n_orig_at,
                &v3.atom_index_orig,
                &v3.atom_index_fin,
            ) as i16;
            atnum2 = get_actual_atom_number(
                i32::from(atnum2),
                n_orig_at,
                &v3.atom_index_orig,
                &v3.atom_index_fin,
            ) as i16;
        }
        if atnum1 < 0 && atnum2 < 0 {
            failed = true;
        }
        if failed {
            if err == 0 {
                treat_err(
                    &mut err,
                    4,
                    p_str_err,
                    "Cannot interpret V3000 bond block line:",
                );
                rd.dotify();
                add_error_message(p_str_err, &rd.line_str());
            }
            if c_str(&rd.line) == SD_FMT_END_OF_DATA {
                err = -err.abs();
                break;
            }
        }
        let mut key = Vec::new();
        while p.is_some() && v3000_read_keyword(&line, &mut p, &mut key) > 0 {
            failed = false;
            match key.as_slice() {
                b"CFG" => {
                    let (v, n) = v3000_read_field_int(V3Int::Char, &line, &mut p);
                    stereo = v as i8;
                    if n < 0 {
                        failed = true;
                    } else if stereo == 2 {
                        stereo = 4;
                    } else if stereo == 3 {
                        stereo = 6;
                    }
                }
                b"TOPO" | b"RXCTR" | b"STBOX" => {
                    v3000_read_field_int(V3Int::Int, &line, &mut p);
                }
                b"ENDPTS" => match molfile_v3000_read_haptic_bond(&line, &mut p) {
                    None => failed = true,
                    Some(mut list) => {
                        let existent = if atnum1 >= 0 { atnum1 } else { atnum2 };
                        if existent < 0 {
                            failed = true;
                        } else {
                            let v3 = ctab.v3000.as_mut().unwrap();
                            let nnum = list[2];
                            list[1] = i32::from(existent);
                            for k in 3..nnum.max(3) as usize {
                                list[k] = get_actual_atom_number(
                                    list[k],
                                    n_orig_at,
                                    &v3.atom_index_orig,
                                    &v3.atom_index_fin,
                                );
                            }
                            v3.haptic_bonds.push(list);
                            is_haptic = true;
                        }
                    }
                },
                b"DISP" | b"ATTACH" => {
                    let mut s = Vec::new();
                    v3000_read_field_str(&line, &mut p, &mut s);
                }
                _ => {}
            }
            if failed {
                if err == 0 {
                    treat_err(
                        &mut err,
                        4,
                        p_str_err,
                        "Cannot interpret V3000 bond block line:",
                    );
                    rd.dotify();
                    add_error_message(p_str_err, &rd.line_str());
                }
                if c_str(&rd.line) == SD_FMT_END_OF_DATA {
                    err = -err.abs();
                    break;
                }
            }
        }
        let v3 = ctab.v3000.as_mut().unwrap();
        if is_haptic {
            let ii = v3.n_haptic_bonds as usize;
            v3.haptic_bonds[ii][0] = i32::from(bond_type);
            v3.n_haptic_bonds += 1;
            continue;
        }
        let ii = v3.n_non_haptic_bonds as usize;
        v3.n_non_haptic_bonds += 1;
        let b = &mut ctab.bonds.as_mut().unwrap()[ii];
        b.atnum1 = atnum1;
        b.atnum2 = atnum2;
        b.bond_type = bond_type;
        b.bond_stereo = stereo;
    }
    let v3 = ctab.v3000.as_ref().unwrap();
    if v3.n_haptic_bonds != 0 {
        add_error_message(p_str_err, "V3000 haptic bonds read/stored but ignored");
        ctab.n_bonds = v3.n_non_haptic_bonds;
    }
    rd.next(inp);
    if rd.p.is_none() || rd.p_text() != b"END BOND" {
        treat_err(
            &mut err,
            1,
            p_str_err,
            "Error: No V3000 Bond block end marker",
        );
        return err;
    }
    rd.rm_lf();
    err
}

/// `MolfileV3000ReadTailOfCTAB`
pub fn molfile_v3000_read_tail_of_ctab(
    ctab: &mut MolFmtCtab,
    inp: &mut InchiIostream,
    err: i32,
    p_str_err: &mut String,
) -> i32 {
    let mut retcode = err;
    let mut rd = LineReader::default();
    rd.next(inp);
    rd.rm_lf();
    if rd.p.is_some() && rd.p_text() == b"BEGIN SGROUP" {
        retcode = molfile_v3000_read_sgroup(inp);
        if retcode != 0 {
            return err;
        }
        rd.next(inp);
        rd.rm_lf();
    }
    if rd.p.is_some() && rd.p_text() == b"BEGIN OBJ3D" {
        retcode = molfile_v3000_read_3d_block(inp, retcode, p_str_err);
        if retcode != 0 {
            return err;
        }
        rd.next(inp);
        rd.rm_lf();
    }
    while rd.p.is_some() && rd.p_text() == b"LINKNODE" {
        rd.next(inp);
        rd.rm_lf();
    }
    while rd.p.is_some() && rd.p_text() == b"BEGIN COLLECTION" {
        retcode = molfile_v3000_read_collections(ctab, inp, retcode, p_str_err);
        if retcode != 0 {
            return err;
        }
        rd.next(inp);
        rd.rm_lf();
    }
    let mut err = err;
    if rd.p.is_none() || rd.p_text() != b"END CTAB" {
        treat_err(&mut err, 1, p_str_err, "Error: No V3000 CTAB end marker");
        return err;
    }
    rd.rm_lf();
    err
}
