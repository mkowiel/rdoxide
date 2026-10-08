//! V3000 CTAB reader (`ParseV3000CTAB` and helpers).
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/MolFileParser.cpp`
//! (2026.09.1). SGroup blocks are not ported yet and are reported as parse
//! errors.
//!
//! Copyright (C) 2002-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::mol_file_parser::{
    COMPLEX_QUERIES, GENERIC_GROUPS, calculate_3d_flag, complex_name_query, fpe,
    lookup_atomic_number, not_ported,
};
use super::utils::{
    LineStream, atof, lexical_cast_int, lossy, to_double_field, to_int, to_unsigned,
};
use crate::geometry::Point3D;
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::{Bond, BondDir, BondStereo, BondType};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::queries::{
    AtomFn, AtomQuery, CompositeQueryType, Query, QueryKind, make_atom_formal_charge_query,
    make_atom_implicit_h_count_query, make_atom_isotope_query, make_atom_null_query,
    make_atom_num_query, make_atom_ring_bond_count_query, make_atom_unsaturated_query,
    make_bond_is_in_ring_query, make_bond_null_query, make_double_or_aromatic_bond_query,
    make_single_or_aromatic_bond_query, make_single_or_double_bond_query,
};
use crate::graphmol::query_ops::{
    atom_expand_query, bond_expand_query, query_bond_from_bond, replace_atom_with_query_atom,
};
use crate::graphmol::romol::RWMol;
use crate::graphmol::stereo_group::{StereoGroup, StereoGroupType};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_warning};

/// `FileParserUtils::getV3000Line`: joins `-` continuation lines and strips
/// the `M  V30 ` prefix.
pub fn get_v3000_line(s: &mut LineStream<'_>, line: &mut u32) -> Result<Vec<u8>> {
    let mut res = Vec::new();
    *line += 1;
    let mut t = s.get_line();
    let bad = |line: u32| fpe(format!("Line {line} does not start with 'M  V30 '\n"));
    if t.len() < 7 || &t[..7] != b"M  V30 " {
        return Err(bad(*line));
    }
    while t.last() == Some(&b'-') {
        res.extend_from_slice(&t[7..t.len() - 1]);
        *line += 1;
        t = s.get_line();
        if t.len() < 7 || &t[..7] != b"M  V30 " {
            return Err(bad(*line));
        }
    }
    res.extend_from_slice(&t[7..]);
    Ok(res)
}

/// `FileParserUtils::strip(s)` with the default " \t\r\n"
fn strip(s: &[u8]) -> &[u8] {
    let ws = |c: &u8| matches!(c, b' ' | b'\t' | b'\r' | b'\n');
    let start = s.iter().position(|c| !ws(c)).unwrap_or(s.len());
    let end = s.iter().rposition(|c| !ws(c)).map_or(start, |e| e + 1);
    &s[start..end]
}

fn upper(s: &mut [u8]) {
    s.make_ascii_uppercase();
}

fn starts(t: &[u8], p: &[u8]) -> bool {
    t.starts_with(p)
}

/// `std::from_chars` for an unsigned value (0 on failure).
fn from_chars_u(t: &[u8]) -> u32 {
    let n = t.iter().take_while(|c| c.is_ascii_digit()).count();
    std::str::from_utf8(&t[..n])
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0)
}

/// `atoi`
fn atoi(t: &[u8]) -> i32 {
    let mut i = 0;
    while i < t.len() && t[i].is_ascii_whitespace() {
        i += 1;
    }
    let start = i;
    if i < t.len() && matches!(t[i], b'+' | b'-') {
        i += 1;
    }
    while i < t.len() && t[i].is_ascii_digit() {
        i += 1;
    }
    std::str::from_utf8(&t[start..i])
        .ok()
        .and_then(|s| s.parse::<i64>().ok())
        .map_or(0, |v| v as i32)
}

/// `tokenizeV3000Line`
pub fn tokenize_v3000_line(line: &[u8]) -> Vec<&[u8]> {
    let mut tokens = Vec::new();
    let mut in_quotes = false;
    let mut paren_depth = 0u32;
    let mut start = 0usize;
    let mut pos = 0usize;
    while pos < line.len() {
        let c = line[pos];
        if c == b' ' || c == b'\t' {
            if start == pos {
                start += 1;
                pos += 1;
            } else if !in_quotes && paren_depth == 0 {
                tokens.push(&line[start..pos]);
                pos += 1;
                start = pos;
            } else {
                pos += 1;
            }
        } else if c == b')' && paren_depth > 0 {
            paren_depth -= 1;
            pos += 1;
        } else if c == b'(' && !in_quotes {
            paren_depth += 1;
            pos += 1;
        } else if c == b'"' && paren_depth == 0 {
            if pos + 1 < line.len() && line[pos + 1] == b'"' {
                pos += 2;
            } else if in_quotes {
                tokens.push(&line[start + 1..pos]);
                pos += 1;
                start = pos;
                in_quotes = false;
            } else {
                pos += 1;
                in_quotes = true;
            }
        } else {
            pos += 1;
        }
    }
    if start != pos {
        tokens.push(&line[start..]);
    }
    tokens
}

/// `splitAssignToken`: `(PROP upper-cased, value)`.
fn split_assign_token(token: &[u8]) -> Option<(Vec<u8>, &[u8])> {
    let first = token.iter().position(|&c| c == b'=')?;
    let last = token.iter().rposition(|&c| c == b'=')?;
    if first != last {
        return None;
    }
    let mut prop = token[..first].to_vec();
    upper(&mut prop);
    Some((prop, &token[first + 1..]))
}

fn r_range(t: &[u8]) -> bool {
    t.first() == Some(&b'R') && t >= &b"R0"[..] && t <= &b"R99"[..]
}

/// `ParseV3000AtomSymbol`
fn parse_v3000_atom_symbol(token: &[u8], line: u32, strict: bool) -> Result<Atom> {
    let mut negate = false;
    let mut token = strip(token);
    if token.len() > 3 && token[..3].eq_ignore_ascii_case(b"NOT") {
        negate = true;
        token = strip(&token[3..]);
    }
    let t0 = token.first().copied().unwrap_or(0);
    if t0 == b'[' {
        if token.last() != Some(&b']') {
            return Err(fpe(format!(
                "Bad atom token '{}' on line: {line}",
                lossy(token)
            )));
        }
        let inner = &token[1..token.len() - 1];
        let mut res: Option<Atom> = None;
        for part in inner.split(|&c| c == b',') {
            let mut sym = strip(part).to_vec();
            if sym.is_empty() {
                continue;
            }
            if sym.len() == 2 && sym[1].is_ascii_uppercase() {
                sym[1] = sym[1].to_ascii_lowercase();
            }
            let at_num =
                crate::graphmol::PeriodicTable::get().get_atomic_number(&lossy(&sym))? as i32;
            match res.as_mut() {
                None => {
                    let mut a = Atom::new(at_num as u32);
                    a.set_query(Some(make_atom_num_query(at_num)));
                    res = Some(a);
                }
                Some(a) => {
                    atom_expand_query(a, make_atom_num_query(at_num), CompositeQueryType::Or, true)?
                }
            }
            res.as_mut().expect("atom").set_atomic_num(0);
        }
        // upstream dereferences a null atom for an empty list
        let mut res = res.ok_or_else(|| RDKitError::invariant("no atom built"))?;
        res.query.as_mut().expect("query").negate = negate;
        return Ok(res);
    }
    if negate {
        return Err(fpe(format!(
            "NOT tokens only supported for atom lists. line {line}"
        )));
    }
    let is_complex = COMPLEX_QUERIES.contains(&token);
    let mut res;
    if is_complex || token == b"R" || r_range(token) || token == b"R#" || token == b"*" {
        if is_complex || token == b"*" {
            res = Atom::new(0);
            res.set_query(Some(if token == b"*" {
                make_atom_null_query()
            } else {
                complex_name_query(token)
            }));
            res.set_no_implicit(true);
        } else {
            res = Atom::new(1);
            res.set_atomic_num(0);
        }
        if r_range(token) {
            let rnumber = lexical_cast_int(&token[1..]).unwrap_or(-1);
            if rnumber >= 0 {
                res.set_isotope(rnumber as u32);
            }
        }
        if t0 == b'R' {
            res.props.set_val(cp::dummyLabel, lossy(token));
        }
    } else if token == b"D" {
        res = Atom::new(1);
        res.set_isotope(2);
    } else if token == b"T" {
        res = Atom::new(1);
        res.set_isotope(3);
    } else if token == b"Pol" || token == b"Mod" {
        res = Atom::new(0);
        res.props.set_val(cp::dummyLabel, lossy(token));
    } else if GENERIC_GROUPS.contains(&token) {
        res = Atom::new(0);
        res.set_query(Some(make_atom_num_query(0)));
        res.props.set_val(cp::atomLabel, lossy(token));
    } else {
        res = Atom::new(0);
        lookup_atomic_number(&mut res, token, strict)?;
    }
    Ok(res)
}

fn and_q(mol: &mut RWMol, idx: u32, q: AtomQuery) -> Result<()> {
    atom_expand_query(mol.atom_mut(idx), q, CompositeQueryType::And, true)
}

/// `ParseV3000RGroups`
fn parse_v3000_rgroups(mol: &mut RWMol, idx: u32, text: &[u8], line: u32) -> Result<()> {
    if text.first() != Some(&b'(') || text.last() != Some(&b')') {
        return Err(fpe(format!(
            "Bad RGROUPS specification '{}' on line {line}. Missing parens.",
            lossy(text)
        )));
    }
    let resid = &text[1..text.len() - 1];
    let split: Vec<&[u8]> = resid.split(|&c| c == b' ').collect();
    let cast = |t: &[u8]| {
        super::utils::strip_spaces_and_cast_unsigned(t).map_err(|_| {
            fpe(format!(
                "Cannot convert '{}' to int on line{line}",
                lossy(t)
            ))
        })
    };
    let n_rs = cast(split[0])?;
    if split.len() < n_rs as usize + 1 {
        return Err(fpe(format!(
            "Bad RGROUPS specification '{}' on line {line}. Not enough values.",
            lossy(text)
        )));
    }
    for i in 0..n_rs as usize {
        let r_label = cast(split[i + 1])?;
        replace_atom_with_query_atom(mol, idx)?;
        let a = mol.atom_mut(idx);
        a.props.set_val(cp::_MolFileRLabel, r_label);
        a.props.set_val(cp::dummyLabel, format!("R{r_label}"));
        a.set_isotope(r_label);
        a.set_query(Some(make_atom_null_query()));
    }
    Ok(())
}

/// `ParseV3000AtomProps`
fn parse_v3000_atom_props(
    mol: &mut RWMol,
    idx: u32,
    tokens: &[&[u8]],
    line: u32,
    strict: bool,
) -> Result<()> {
    for &tok in tokens {
        let Some((prop, val)) = split_assign_token(tok) else {
            return Err(fpe(format!(
                "Invalid atom property: '{}' for atom {} on line {line}\n",
                lossy(tok),
                idx + 1
            )));
        };
        let ti = |v: &[u8]| to_int(v, true);
        match prop.as_slice() {
            b"CHG" => {
                let charge = ti(val)?;
                if !mol.atom(idx).has_query() {
                    mol.atom_mut(idx).set_formal_charge(charge);
                } else {
                    and_q(mol, idx, make_atom_formal_charge_query(charge))?;
                }
            }
            b"RAD" => {
                let n = match ti(val)? {
                    0 => None,
                    1 | 3 => Some(2),
                    2 => Some(1),
                    _ => {
                        return Err(fpe(format!(
                            "Unrecognized RAD value {} for atom {} on line {line}\n",
                            lossy(val),
                            idx + 1
                        )));
                    }
                };
                if let Some(n) = n {
                    mol.atom_mut(idx).set_num_radical_electrons(n);
                }
            }
            b"MASS" => {
                let v = match ti(val) {
                    Ok(v) => v,
                    Err(_) => match to_double_field(val, 0, val.len()) {
                        Ok(dv) => dv.floor() as i32,
                        Err(_) => -1,
                    },
                };
                if v < 0 {
                    return Err(fpe(format!(
                        "Bad value for MASS :{} for atom {} on line {line}\n",
                        lossy(val),
                        idx + 1
                    )));
                }
                if !mol.atom(idx).has_query() {
                    mol.atom_mut(idx).set_isotope(v as u32);
                } else {
                    and_q(mol, idx, make_atom_isotope_query(v))?;
                }
            }
            b"CFG" => {
                let cfg = ti(val)?;
                match cfg {
                    0 => {}
                    1..=3 => mol.atom_mut(idx).props.set_val(cp::molParity, cfg),
                    _ => {
                        return Err(fpe(format!(
                            "Unrecognized CFG value : {} for atom {} on line {line}\n",
                            lossy(val),
                            idx + 1
                        )));
                    }
                }
            }
            b"HCOUNT" => {
                if val != b"0" {
                    let mut hcount = ti(val)?;
                    replace_atom_with_query_atom(mol, idx)?;
                    if hcount == -1 {
                        hcount = 0;
                    }
                    let q = if hcount > 0 {
                        Query::with(
                            QueryKind::LessEqual,
                            hcount,
                            AtomFn::ImplicitHCount,
                            "less_AtomImplicitHCount",
                        )
                    } else {
                        make_atom_implicit_h_count_query(0)
                    };
                    and_q(mol, idx, q)?;
                }
            }
            b"UNSAT" => {
                if val == b"1" {
                    replace_atom_with_query_atom(mol, idx)?;
                    and_q(mol, idx, make_atom_unsaturated_query())?;
                }
            }
            b"RBCNT" => {
                if val != b"0" {
                    let mut rb = ti(val)?;
                    replace_atom_with_query_atom(mol, idx)?;
                    mol.atom_mut(idx).props.set_val(cp::molRingBondCount, rb);
                    if rb == -1 {
                        rb = 0;
                    } else if rb == -2 {
                        mol.props.set_val(cp::_NeedsQueryScan, 1i32);
                        rb = 0xDEADBEEFu32 as i32;
                    } else if rb > 4 {
                        rb = 4;
                    }
                    and_q(mol, idx, make_atom_ring_bond_count_query(rb))?;
                }
            }
            b"VAL" => {
                if val != b"0" {
                    let v = ti(val)?;
                    mol.atom_mut(idx).props.set_val(cp::molTotValence, v);
                }
            }
            b"RGROUPS" => parse_v3000_rgroups(mol, idx, val, line)?,
            b"STBOX" | b"SUBST" | b"EXACHG" | b"INVRET" | b"SEQID" => {
                if val != b"0" {
                    let v = ti(val)?;
                    let key = match prop.as_slice() {
                        b"STBOX" => cp::molStereoCare,
                        b"SUBST" => cp::molSubstCount,
                        b"EXACHG" => cp::molRxnExactChange,
                        b"INVRET" => cp::molInversionFlag,
                        _ => cp::molAtomSeqId,
                    };
                    mol.atom_mut(idx).props.set_val(key, v);
                }
            }
            b"ATTCHPT" => {
                if val != b"0" {
                    let v = ti(val)?;
                    if mol.atom(idx).has_prop(cp::molAttachPoint) {
                        let msg = format!(
                            "Multiple ATTCHPT values for atom {} on line {line}",
                            idx + 1
                        );
                        if strict {
                            return Err(fpe(msg));
                        }
                        rd_log_warning!("{msg}");
                    } else {
                        mol.atom_mut(idx).props.set_val(cp::molAttachPoint, v);
                    }
                }
            }
            b"ATTCHORD" => {
                if val.first() == Some(&b'(') {
                    return Err(not_ported("ATTCHORD templates", line));
                }
                let v = ti(val)?;
                mol.atom_mut(idx).props.set_val(cp::molAttachOrder, v);
            }
            b"CLASS" => mol
                .atom_mut(idx)
                .props
                .set_val(cp::molAtomClass, lossy(val)),
            b"SEQNAME" if !val.is_empty() => {
                mol.atom_mut(idx)
                    .props
                    .set_val(cp::molAtomSeqName, lossy(val));
            }
            _ => {}
        }
    }
    Ok(())
}

/// `ParseV3000AtomBlock`
fn parse_v3000_atom_block(
    s: &mut LineStream<'_>,
    line: &mut u32,
    n_atoms: u32,
    mol: &mut RWMol,
    conf: &mut Conformer,
    strict: bool,
    expect_macro_atoms: bool,
) -> Result<()> {
    let t = get_v3000_line(s, line)?;
    if !starts(&t, b"BEGIN ATOM") {
        return Err(fpe(format!("BEGIN ATOM line not found on line {line}")));
    }
    for _ in 0..n_atoms {
        let t = get_v3000_line(s, line)?;
        let tokens = tokenize_v3000_line(strip(&t));
        let bad = |line: u32| fpe(format!("Bad atom line : '{}' on line {line}", lossy(&t)));
        if tokens.len() < 2 {
            return Err(bad(*line));
        }
        let mol_idx = from_chars_u(tokens[0]);
        let is_macro = expect_macro_atoms
            && tokens[2..]
                .iter()
                .any(|t| split_assign_token(t).is_some_and(|(p, _)| p == b"CLASS"));
        let atom = if is_macro {
            let mut a = Atom::new(0);
            a.props.set_val(cp::dummyLabel, lossy(tokens[1]));
            a
        } else {
            parse_v3000_atom_symbol(tokens[1], *line, strict)?
        };
        if tokens.len() < 6 {
            return Err(bad(*line));
        }
        let pos = Point3D::new(atof(tokens[2]), atof(tokens[3]), atof(tokens[4]));
        let mut atom = atom;
        let map_num = atoi(tokens[5]);
        if map_num > 0 {
            atom.props.set_val(cp::molAtomMapNumber, map_num);
        }
        let aid = mol.add_atom(atom, false);
        parse_v3000_atom_props(mol, aid, &tokens[6..], *line, strict)?;
        mol.set_atom_bookmark(aid, mol_idx as i32);
        conf.set_atom_pos(aid, pos)?;
    }
    let t = get_v3000_line(s, line)?;
    if !starts(&t, b"END ATOM") {
        return Err(fpe(format!("END ATOM line not found on line {line}")));
    }
    Ok(())
}

/// `ParseV3000BondBlock`
fn parse_v3000_bond_block(
    s: &mut LineStream<'_>,
    line: &mut u32,
    n_bonds: u32,
    mol: &mut RWMol,
    chirality_possible: &mut bool,
) -> Result<()> {
    let t = get_v3000_line(s, line)?;
    if !starts(&t, b"BEGIN BOND") {
        return Err(fpe("BEGIN BOND line not found"));
    }
    for _ in 0..n_bonds {
        let t = get_v3000_line(s, line)?;
        let split = tokenize_v3000_line(strip(&t));
        if split.len() < 4 {
            return Err(fpe(format!("bond line {line} is too short")));
        }
        let bond_idx = from_chars_u(split[0]);
        let b_type = from_chars_u(split[1]);
        let a1 = from_chars_u(split[2]);
        let a2 = from_chars_u(split[3]);
        let mut bond = match b_type {
            1 => Bond::new(BondType::Single),
            2 => Bond::new(BondType::Double),
            3 => Bond::new(BondType::Triple),
            4 => {
                let mut b = Bond::new(BondType::Aromatic);
                b.set_is_aromatic(true);
                b
            }
            9 => Bond::new(BondType::Dative),
            10 => Bond::new(BondType::Hydrogen),
            0 => {
                rd_log_warning!(
                    "bond with order 0 found on line {line}. This is not part of the MDL specification."
                );
                Bond::new(BondType::Unspecified)
            }
            _ => {
                let mut b = Bond::new(BondType::Unspecified);
                match b_type {
                    8 => b.set_query(Some(make_bond_null_query())),
                    5..=7 => {
                        b.set_query(Some(match b_type {
                            5 => make_single_or_double_bond_query(),
                            6 => make_single_or_aromatic_bond_query(),
                            _ => make_double_or_aromatic_bond_query(),
                        }));
                        b.props.set_val(cp::_MolFileBondQuery, 1i32);
                    }
                    _ => {
                        b.set_query(Some(make_bond_null_query()));
                        rd_log_warning!(
                            "unrecognized query bond type, {b_type}, found on line {line}. Using an \"any\" query."
                        );
                    }
                }
                b
            }
        };
        bond.props.set_val(cp::_MolFileBondType, b_type);
        for &tok in &split[4..] {
            let Some((prop, val)) = split_assign_token(tok) else {
                return Err(fpe(format!(
                    "bad bond property '{}' on line {line}",
                    lossy(tok)
                )));
            };
            match prop.as_slice() {
                b"CFG" => {
                    let cfg = from_chars_u(val);
                    match cfg {
                        0 => {}
                        1 => {
                            bond.set_bond_dir(BondDir::BeginWedge);
                            *chirality_possible = true;
                        }
                        2 => {
                            if b_type == 1 {
                                bond.set_bond_dir(BondDir::Unknown);
                            } else if b_type == 2 {
                                bond.set_bond_dir(BondDir::EitherDouble);
                                bond.set_stereo(BondStereo::StereoAny)?;
                            }
                        }
                        3 => {
                            bond.set_bond_dir(BondDir::BeginDash);
                            *chirality_possible = true;
                        }
                        _ => {
                            return Err(fpe(format!(
                                "bad bond CFG {}' on line {line}",
                                lossy(val)
                            )));
                        }
                    }
                    bond.props.set_val(cp::_MolFileBondCfg, cfg);
                }
                b"TOPO" => {
                    if val != b"0" {
                        if !bond.has_query() {
                            bond = query_bond_from_bond(&bond);
                        }
                        let mut q = make_bond_is_in_ring_query();
                        if val == b"1" {
                        } else if val == b"2" {
                            q.negate = true;
                        } else {
                            return Err(fpe(format!(
                                "bad bond TOPO {}' on line {line}",
                                lossy(val)
                            )));
                        }
                        bond_expand_query(&mut bond, q, CompositeQueryType::And, true)?;
                    }
                }
                b"RXCTR" => {
                    let v = to_int(val, true)?;
                    bond.props.set_val(cp::molReactStatus, v);
                }
                b"STBOX" => bond.props.set_val(cp::molStereoCare, lossy(val)),
                b"ENDPTS" => bond.props.set_val(cp::_MolFileBondEndPts, lossy(val)),
                b"ATTACH" => bond.props.set_val(cp::_MolFileBondAttach, lossy(val)),
                _ => {}
            }
        }
        let b = mol.get_atom_with_bookmark(a1 as i32)?;
        let e = mol.get_atom_with_bookmark(a2 as i32)?;
        bond.set_begin_atom_idx_unchecked(b);
        bond.set_end_atom_idx_unchecked(e);
        let bidx = mol.add_bond_object(bond)? - 1;
        mol.set_bond_bookmark(bidx, bond_idx as i32);
        if !mol.bond(bidx).has_prop(cp::molStereoCare) {
            let care = |a: u32| {
                mol.atom(a)
                    .props
                    .get(cp::molStereoCare)
                    .and_then(RDValue::as_int)
            };
            if let (Some(c1), Some(c2)) = (care(b), care(e))
                && c1 == c2
            {
                mol.bond_mut(bidx).props.set_val(cp::molStereoCare, c1);
            }
        }
    }
    let t = get_v3000_line(s, line)?;
    if !starts(&t, b"END BOND") {
        return Err(fpe(format!("END BOND line not found at line {line}")));
    }
    Ok(())
}

/// Matches `MDLV30/STE(...)([0-9]*) +ATOMS=\(([0-9]+) +(.*)\) *` against
/// the whole line.
/// `(type, id, count, atoms)` captures of the stereo-group regex.
type StereoLabel<'a> = (&'a [u8], &'a [u8], &'a [u8], &'a [u8]);

fn match_stereo_label(t: &[u8]) -> Option<StereoLabel<'_>> {
    let rest = t.strip_prefix(b"MDLV30/STE")?;
    if rest.len() < 3 || rest[..3].contains(&b'\n') {
        return None;
    }
    let (kind, rest) = rest.split_at(3);
    let nd = rest.iter().take_while(|c| c.is_ascii_digit()).count();
    let (id, rest) = rest.split_at(nd);
    let ns = rest.iter().take_while(|&&c| c == b' ').count();
    if ns == 0 {
        return None;
    }
    let rest = rest[ns..].strip_prefix(b"ATOMS=(")?;
    let nd = rest.iter().take_while(|c| c.is_ascii_digit()).count();
    if nd == 0 {
        return None;
    }
    let (count, rest) = rest.split_at(nd);
    let ns = rest.iter().take_while(|&&c| c == b' ').count();
    if ns == 0 {
        return None;
    }
    let rest = &rest[ns..];
    // `(.*)\) *$`: the last ')' followed only by spaces
    let close = rest.iter().rposition(|&c| c != b' ')?;
    if rest[close] != b')' || rest[..close].contains(&b'\n') {
        return None;
    }
    Some((kind, id, count, &rest[..close]))
}

/// `parseEnhancedStereo`; returns the line after the collection block.
fn parse_enhanced_stereo(
    s: &mut LineStream<'_>,
    line: &mut u32,
    mol: &mut RWMol,
    strict: bool,
) -> Result<Vec<u8>> {
    let mut t = get_v3000_line(s, line)?;
    upper(&mut t);
    let mut groups = Vec::new();
    let mut abs_seen = 0;
    while !starts(&t, b"END") {
        if let Some((kind, id, count, atoms_txt)) = match_stereo_label(&t) {
            let (gtype, gid) = match kind {
                b"ABS" => {
                    if abs_seen == 1 {
                        let msg = format!("Seen a second ABS stereo group on line {line}\n");
                        if strict {
                            return Err(fpe(msg));
                        }
                        rd_log_warning!("{msg}");
                    }
                    abs_seen += 1;
                    (StereoGroupType::Absolute, 0)
                }
                b"REL" => (StereoGroupType::Or, to_unsigned(id, true)?),
                b"RAC" => (StereoGroupType::And, to_unsigned(id, true)?),
                _ => {
                    return Err(fpe(format!(
                        "Unrecognized stereogroup type : '{}' on line{line}",
                        lossy(&t)
                    )));
                }
            };
            let count = to_unsigned(count, true)?;
            let mut nums = atoms_txt
                .split(|c| c.is_ascii_whitespace())
                .filter(|x| !x.is_empty());
            let mut atoms: Vec<u32> = Vec::new();
            // `ss >> index`: once extraction fails, index is 0 from then on
            let mut failed = false;
            for _ in 0..count {
                let mut index = 0u32;
                if !failed {
                    match nums
                        .next()
                        .and_then(|x| std::str::from_utf8(x).ok().and_then(|x| x.parse().ok()))
                    {
                        Some(v) => index = v,
                        None => failed = true,
                    }
                }
                let a = index.wrapping_sub(1);
                mol.get_atom(a)?;
                if atoms.contains(&a) {
                    let msg = format!(
                        "Atom {index} appears more than once in stereo group specification on line {line}!"
                    );
                    if strict {
                        return Err(fpe(msg));
                    }
                    rd_log_warning!("{msg}");
                } else {
                    atoms.push(a);
                }
            }
            groups.push(StereoGroup::new(gtype, atoms, Vec::new(), gid)?);
        } else {
            rd_log_warning!(
                "Skipping unrecognized collection type at line {line}: {}",
                lossy(&t)
            );
        }
        t = get_v3000_line(s, line)?;
    }
    if !groups.is_empty() {
        mol.set_stereo_groups(groups)?;
    }
    get_v3000_line(s, line)
}

fn counts(t: &[u8]) -> Vec<&[u8]> {
    t.split(|&c| c == b' ' || c == b'\t')
        .filter(|x| !x.is_empty())
        .collect()
}

/// `FileParserUtils::ParseV3000CTAB`
#[allow(clippy::too_many_arguments)]
pub fn parse_v3000_ctab(
    s: &mut LineStream<'_>,
    line: &mut u32,
    mol: &mut RWMol,
    chirality_possible: &mut bool,
    strict: bool,
    expect_mend: bool,
    expect_macro_atoms: bool,
) -> Result<bool> {
    let mut t = get_v3000_line(s, line)?;
    upper(&mut t);
    if !starts(&t, b"BEGIN CTAB") {
        return Err(fpe(format!("BEGIN CTAB line not found on line {line}")));
    }
    let mut t = get_v3000_line(s, line)?;
    upper(&mut t);
    if t.len() < 8 || &t[..7] != b"COUNTS " {
        return Err(fpe(format!(
            "Bad counts line : '{}' on line {line}",
            lossy(&t)
        )));
    }
    let trimmed = super::mol_file_parser::trim(&t[7..]).to_vec();
    let split = counts(&trimmed);
    if split.len() < 2 {
        return Err(fpe(format!(
            "Bad counts line : '{}' on line {line}",
            lossy(&t)
        )));
    }
    let n_atoms = to_unsigned(split[0], true)?;
    let n_bonds = to_unsigned(split[1], true)?;
    let mut conf = Conformer::new(n_atoms as usize);
    let mut n_sgroups = if split.len() > 2 {
        to_unsigned(split[2], true)?
    } else {
        0
    };
    let n_3d = if split.len() > 3 {
        to_unsigned(split[3], true)?
    } else {
        0
    };
    let chiral_flag = if split.len() > 4 {
        to_unsigned(split[4], true)?
    } else {
        0
    };
    mol.props.set_val(cp::_MolFileChiralFlag, chiral_flag);

    if n_atoms > 0 {
        parse_v3000_atom_block(s, line, n_atoms, mol, &mut conf, strict, expect_macro_atoms)?;
    }
    if n_bonds > 0 {
        parse_v3000_bond_block(s, line, n_bonds, mol, chirality_possible)?;
    }

    let mut t = get_v3000_line(s, line)?;
    upper(&mut t);
    while t.len() > 8 && starts(&t, b"LINKNODE") {
        if t.len() > 9 {
            let mut existing = mol
                .props
                .get(cp::molFileLinkNodes)
                .and_then(|v| match v {
                    RDValue::String(s) => Some(s.clone()),
                    _ => None,
                })
                .map(|e| e + "|")
                .unwrap_or_default();
            existing += &lossy(&t[9..]);
            mol.props.set_val(cp::molFileLinkNodes, existing);
        }
        t = get_v3000_line(s, line)?;
    }

    let sgroup_found = false;
    let mut obj3d_found = false;
    upper(&mut t);
    while t.len() > 5 && starts(&t, b"BEGIN") {
        if starts(&t, b"BEGIN SGROUP") {
            if sgroup_found {
                return Err(fpe(format!(
                    "BEGIN SGROUP found more than once on line {line}"
                )));
            } else if n_sgroups == 0 {
                let msg = format!("BEGIN SGROUP  found but Sgroups NOT expected on line {line}");
                if strict {
                    return Err(fpe(msg));
                }
                rd_log_warning!("{msg}");
                n_sgroups = u32::MAX;
            }
            let _ = (n_sgroups, sgroup_found);
            return Err(not_ported("V3000 SGroup data", *line));
        } else if t.len() >= 16 && &t[6..16] == b"COLLECTION" {
            t = parse_enhanced_stereo(s, line, mol, strict)?;
            upper(&mut t);
        } else if starts(&t, b"BEGIN OBJ3D") {
            if obj3d_found {
                return Err(fpe(format!(
                    "BEGIN OBJ3D found more than once on line {line}"
                )));
            }
            if n_3d == 0 {
                let msg =
                    format!("BEGIN OBJ3D found but 3n3DConstraints NOT expected on line {line}");
                if strict {
                    return Err(fpe(msg));
                }
                rd_log_warning!("{msg}");
            }
            rd_log_warning!("3D constraint information in mol block ignored at line {line}");
            obj3d_found = true;
            for _ in 0..n_3d {
                get_v3000_line(s, line)?;
            }
            t = get_v3000_line(s, line)?;
            upper(&mut t);
            if !starts(&t, b"END OBJ3D") {
                let msg = format!("END OBJ3D line not found on line {line}");
                if strict {
                    return Err(fpe(msg));
                }
                rd_log_warning!("{msg}");
            }
            t = get_v3000_line(s, line)?;
            upper(&mut t);
        } else {
            rd_log_warning!("skipping block at line {line}: '{}'", lossy(&t));
            while !starts(&t, b"END") {
                t = get_v3000_line(s, line)?;
            }
            t = get_v3000_line(s, line)?;
            upper(&mut t);
        }
    }
    if n_sgroups != 0 && !sgroup_found {
        let msg = format!("BEGIN SGROUP line not found on line {line}");
        if strict {
            return Err(fpe(msg));
        }
        rd_log_warning!("{msg}");
    }
    if n_3d != 0 && !obj3d_found {
        let msg = format!("BEGIN OBJ3D line not found on line {line}");
        if strict {
            return Err(fpe(msg));
        }
        rd_log_warning!("{msg}");
    }
    upper(&mut t);
    if !starts(&t, b"END CTAB") {
        if strict {
            return Err(fpe("END CTAB line not found"));
        }
        rd_log_warning!("END CTAB line not found.");
    }
    let file_complete = if expect_mend {
        let t = s.get_line();
        *line += 1;
        t.starts_with(b"M  END")
    } else {
        true
    };
    let is_3d = calculate_3d_flag(mol, &conf, *chirality_possible);
    conf.set_3d(is_3d);
    mol.add_conformer(conf, true)?;
    Ok(file_complete)
}
