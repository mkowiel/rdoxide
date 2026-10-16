//! MDL Mol block reader (`MolFromMolBlock`, `MolFromMolDataStream`).
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/MolFileParser.cpp`
//! (2026.09.1). V2000 CTABs (SGroups included, see `sgroup_parsing`) and
//! V3000 blocks are ported; a few rarely used property lines
//! (MRV/APO/LIN) are reported as parse errors until they are ported.
//!
//! Copyright (C) 2002-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::sgroup_parsing::{V2000SGroups, is_sgroup_line, process_sgroups};
use super::utils::{
    LineStream, is_bad_cast, lexical_cast_int, lossy, strip_spaces, strip_spaces_and_cast_unsigned,
    substr, substr_from, to_double_field, to_int, to_unsigned,
};
use crate::geometry::Point3D;
use crate::graphmol::atom::{Atom, ValenceType};
use crate::graphmol::attachment_points::expand_attachment_points;
use crate::graphmol::bond::{Bond, BondDir, BondStereo, BondType};
use crate::graphmol::chirality::clear_single_bond_dir_flags;
use crate::graphmol::chirality::from_coords::{
    assign_chiral_types_from_3d, assign_chiral_types_from_bond_dirs, detect_bond_stereochemistry,
};
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::conformer::{Conformer, has_non_zero_z_coords};
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, remove_hs};
use crate::graphmol::mol_ops::sanitize_flags::{SANITIZE_ALL, SANITIZE_CLEANUP};
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::queries::{
    AtomFn, AtomQuery, CompositeQueryType, Query, QueryKind, complete_mol_queries,
    make_a_atom_query, make_ah_atom_query, make_atom_explicit_degree_query,
    make_atom_implicit_h_count_query, make_atom_null_query, make_atom_num_query,
    make_atom_ring_bond_count_query, make_atom_unsaturated_query, make_bond_is_in_ring_query,
    make_bond_null_query, make_double_or_aromatic_bond_query, make_m_atom_query,
    make_mh_atom_query, make_q_atom_query, make_qh_atom_query, make_recursive_structure_query,
    make_single_or_aromatic_bond_query, make_single_or_double_bond_query, make_x_atom_query,
    make_xh_atom_query,
};
use crate::graphmol::query_ops::{
    atom_expand_query, bond_expand_query, query_atom_from_atom, query_bond_from_bond,
    replace_atom_with_query_atom,
};
use crate::graphmol::romol::RWMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_error, rd_log_warning};

/// `v2::FileParsers::MolFileParserParams`
#[derive(Debug, Clone, PartialEq)]
pub struct MolFileParserParams {
    pub sanitize: bool,
    pub remove_hs: bool,
    pub strict_parsing: bool,
    pub expand_attachment_points: bool,
    pub parsing_scsr_mol: bool,
}

impl Default for MolFileParserParams {
    fn default() -> Self {
        MolFileParserParams {
            sanitize: true,
            remove_hs: true,
            strict_parsing: true,
            expand_attachment_points: false,
            parsing_scsr_mol: false,
        }
    }
}

pub(super) fn fpe(msg: impl Into<String>) -> RDKitError {
    RDKitError::FileParse(msg.into())
}

pub(super) fn not_ported(what: &str, line: u32) -> RDKitError {
    fpe(format!(
        "{what} on line {line} is not supported by rdoxide yet"
    ))
}

/// Maps a failed `lexical_cast` to the `FileParseException` upstream
/// throws; other errors (e.g. `std::out_of_range`) pass through.
pub(super) fn cast_err<T>(r: Result<T>, msg: impl FnOnce() -> String) -> Result<T> {
    r.map_err(|e| if is_bad_cast(&e) { fpe(msg()) } else { e })
}

fn atom_idx(mol: &RWMol, one_based: i64) -> Result<u32> {
    let idx = one_based - 1;
    let n = mol.num_atoms();
    if idx < 0 || idx >= n as i64 {
        // getAtomWithIdx(unsigned) -> URANGE_CHECK
        return Err(RDKitError::range("idx", (idx as u32) as i64, n));
    }
    Ok(idx as u32)
}

pub(super) const COMPLEX_QUERIES: [&[u8]; 8] = [b"A", b"AH", b"Q", b"QH", b"X", b"XH", b"M", b"MH"];

/// `convertComplexNameToQuery`
pub(super) fn complex_name_query(symb: &[u8]) -> AtomQuery {
    match symb {
        b"Q" => make_q_atom_query(),
        b"QH" => make_qh_atom_query(),
        b"A" => make_a_atom_query(),
        b"AH" => make_ah_atom_query(),
        b"X" => make_x_atom_query(),
        b"XH" => make_xh_atom_query(),
        b"M" => make_m_atom_query(),
        _ => make_mh_atom_query(),
    }
}

/// Keys of `GenericGroups::genericMatchers`.
pub(super) const GENERIC_GROUPS: [&[u8]; 68] = [
    b"Group",
    b"G",
    b"GroupH",
    b"GH",
    b"Group*",
    b"G*",
    b"GroupH*",
    b"GH*",
    b"Alkyl",
    b"ALK",
    b"AlkylH",
    b"ALH",
    b"Alkenyl",
    b"AEL",
    b"AlkenylH",
    b"AEH",
    b"Alkynyl",
    b"AYL",
    b"AlkynylH",
    b"AYH",
    b"Carbocyclic",
    b"CBC",
    b"CarbocyclicH",
    b"CBH",
    b"Carbocycloalkyl",
    b"CAL",
    b"CarbocycloalkylH",
    b"CAH",
    b"Carbocycloalkenyl",
    b"CEL",
    b"CarbocycloalkenylH",
    b"CEH",
    b"Carboaryl",
    b"ARY",
    b"CarboarylH",
    b"ARH",
    b"Cyclic",
    b"CYC",
    b"CyclicH",
    b"CYH",
    b"Acyclic",
    b"ACY",
    b"AcyclicH",
    b"ACH",
    b"Carboacyclic",
    b"ABC",
    b"CarboacyclicH",
    b"ABH",
    b"Heteroacyclic",
    b"AHC",
    b"HeteroacyclicH",
    b"AHH",
    b"Alkoxy",
    b"AOX",
    b"AlkoxyH",
    b"AOH",
    b"Heterocyclic",
    b"CHC",
    b"HeterocyclicH",
    b"CHH",
    b"Heteroaryl",
    b"HAR",
    b"HeteroarylH",
    b"HAH",
    b"NoCarbonRing",
    b"CXX",
    b"NoCarbonRingH",
    b"CXH",
];

pub(super) fn trim(s: &[u8]) -> &[u8] {
    let ws = |c: &u8| c.is_ascii_whitespace() || *c == 0x0b;
    let start = s.iter().position(|c| !ws(c)).unwrap_or(s.len());
    let end = s.iter().rposition(|c| !ws(c)).map_or(start, |e| e + 1);
    &s[start..end]
}

/// `lookupAtomicNumber`
pub(super) fn lookup_atomic_number(res: &mut Atom, symb: &[u8], strict: bool) -> Result<()> {
    let mut t = symb.to_vec();
    if t.len() == 2 && t[1].is_ascii_uppercase() {
        t[1] = t[1].to_ascii_lowercase();
    }
    match PeriodicTable::get().get_atomic_number(&lossy(&t)) {
        Ok(n) => res.set_atomic_num(n as i32),
        Err(e) => {
            if strict || symb.is_empty() {
                return Err(fpe(e.to_string()));
            }
            res.set_atomic_num(0);
            res.props.set_val(cp::dummyLabel, lossy(symb));
        }
    }
    Ok(())
}

/// An optional integer column: parsed only when the line is long enough
/// and the field differs from `skip`.
fn int_field(text: &[u8], pos: usize, len: usize, skip: &[u8], line: u32) -> Result<Option<i32>> {
    if text.len() < pos + len || &text[pos..pos + len] == skip {
        return Ok(None);
    }
    let f = &text[pos..pos + len];
    cast_err(to_int(f, true).map(Some), || {
        format!("Cannot convert '{}' to int on line {line}", lossy(f))
    })
}

/// `ParseMolFileAtomLine`
pub fn parse_mol_file_atom_line(text: &[u8], line: u32, strict: bool) -> Result<(Atom, Point3D)> {
    if (strict && text.len() < 34) || text.len() < 32 {
        return Err(fpe(format!(
            "Atom line too short: '{}' on line {line}",
            lossy(text)
        )));
    }
    let coords = || -> Result<Point3D> {
        Ok(Point3D::new(
            to_double_field(text, 0, 10)?,
            to_double_field(text, 10, 10)?,
            to_double_field(text, 20, 10)?,
        ))
    };
    let pos = cast_err(coords(), || {
        format!("Cannot process coordinates on line {line}")
    })?;
    let symb = trim(substr(text, 31, 3)?).to_vec();

    let mass_diff = int_field(text, 34, 2, b" 0", line)?.unwrap_or(0);
    let chg = int_field(text, 36, 3, b"  0", line)?.unwrap_or(0);
    let h_count = int_field(text, 42, 3, b"  0", line)?.unwrap_or(0);

    let mut res = Atom::new(0);
    let s0 = symb.first().copied().unwrap_or(0);
    let is_complex = COMPLEX_QUERIES.contains(&symb.as_slice());
    let in_r_range = s0 == b'R' && symb.as_slice() >= &b"R0"[..] && symb.as_slice() <= &b"R99"[..];
    if is_complex
        || symb == b"L"
        || symb == b"*"
        || symb == b"LP"
        || symb == b"R"
        || symb == b"R#"
        || in_r_range
    {
        if is_complex || symb == b"*" || symb == b"R" {
            let mut q = Atom::new(0);
            if symb == b"*" || symb == b"R" {
                q.set_query(Some(make_atom_null_query()));
            } else {
                q.set_query(Some(complex_name_query(&symb)));
            }
            res = q;
            res.set_no_implicit(true);
        } else {
            res.set_atomic_num(0);
        }
        if mass_diff == 0 && s0 == b'R' && symb.len() > 1 && in_r_range {
            let rnumber = lexical_cast_int(&symb[1..]).unwrap_or(-1);
            if rnumber >= 0 {
                res.set_isotope(rnumber as u32);
            }
        }
        if s0 == b'R' {
            res.props.set_val(cp::dummyLabel, lossy(&symb));
        }
    } else if symb == b"D" {
        res.set_atomic_num(1);
        res.set_isotope(2);
    } else if symb == b"T" {
        res.set_atomic_num(1);
        res.set_isotope(3);
    } else if symb == b"Pol" || symb == b"Mod" {
        res.set_atomic_num(0);
        res.props.set_val(cp::dummyLabel, lossy(&symb));
    } else if GENERIC_GROUPS.contains(&symb.as_slice()) {
        res = Atom::new(0);
        res.set_query(Some(make_atom_num_query(0)));
        res.props.set_val(cp::atomLabel, lossy(&symb));
    } else {
        lookup_atomic_number(&mut res, &symb, strict)?;
    }

    if chg != 0 {
        res.set_formal_charge(4 - chg);
    }

    if h_count >= 1 {
        if !res.has_query() {
            res = query_atom_from_atom(&res)?;
        }
        res.set_no_implicit(true);
        let q = if h_count > 1 {
            Query::with(
                QueryKind::LessEqual,
                h_count - 1,
                AtomFn::ImplicitHCount,
                "less_AtomImplicitHCount",
            )
        } else {
            make_atom_implicit_h_count_query(0)
        };
        atom_expand_query(&mut res, q, CompositeQueryType::And, true)?;
    }

    if mass_diff != 0 {
        let def_iso = PeriodicTable::get().get_most_common_isotope(res.atomic_num())?;
        let d_iso = def_iso + mass_diff;
        if d_iso < 0 {
            rd_log_warning!(
                " atom {} has a negative isotope offset. line:  {line}",
                res.idx()
            );
        }
        res.set_isotope(d_iso as u32);
    }

    if let Some(v) = int_field(text, 39, 3, b"  0", line)? {
        res.props.set_val(cp::molParity, v);
    }
    if let Some(v) = int_field(text, 45, 3, b"  0", line)? {
        res.props.set_val(cp::molStereoCare, v);
    }
    if let Some(v) = int_field(text, 48, 3, b"  0", line)?
        && v != 0
    {
        res.props.set_val(cp::molTotValence, v);
    }
    if let Some(v) = int_field(text, 54, 3, b"  0", line)?
        && v != 0
    {
        res.props.set_val(cp::molRxnRole, v);
    }
    if let Some(v) = int_field(text, 57, 3, b"  0", line)?
        && v != 0
    {
        res.props.set_val(cp::molRxnComponent, v);
    }
    if let Some(v) = int_field(text, 60, 3, b"  0", line)? {
        res.props.set_val(cp::molAtomMapNumber, v);
    }
    if let Some(v) = int_field(text, 63, 3, b"  0", line)? {
        res.props.set_val(cp::molInversionFlag, v);
    }
    if let Some(v) = int_field(text, 66, 3, b"  0", line)? {
        res.props.set_val(cp::molRxnExactChange, v);
    }
    Ok((res, pos))
}

/// `ParseMolFileBondLine`; returns the bond and its (0-based) end atoms.
pub fn parse_mol_file_bond_line(text: &[u8], line: u32) -> Result<(Bond, u32, u32)> {
    if text.len() < 9 {
        return Err(fpe(format!(
            "Bond line too short: '{}' on line {line}",
            lossy(text)
        )));
    }
    let mut spos = 0;
    let mut field = || -> Result<u32> {
        let f = substr(text, spos, 3)?;
        let v = cast_err(to_unsigned(f, true), || {
            format!("Cannot convert '{}' to int on line {line}", lossy(f))
        })?;
        spos += 3;
        Ok(v)
    };
    let idx1 = field()?.wrapping_sub(1);
    let idx2 = field()?.wrapping_sub(1);
    let b_type = field()?;

    let (mut res, btype) = match b_type {
        1 => (Bond::new(BondType::Single), BondType::Single),
        2 => (Bond::new(BondType::Double), BondType::Double),
        3 => (Bond::new(BondType::Triple), BondType::Triple),
        4 => (Bond::new(BondType::Aromatic), BondType::Aromatic),
        9 => (Bond::new(BondType::Dative), BondType::Dative),
        0 => {
            rd_log_warning!(
                "bond with order 0 found on line {line}. This is not part of the MDL specification."
            );
            (Bond::new(BondType::Unspecified), BondType::Unspecified)
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
            (b, BondType::Unspecified)
        }
    };
    res.set_bond_type(btype);
    res.props.set_val(cp::_MolFileBondType, b_type);

    if text.len() >= 12
        && &text[9..12] != b"  0"
        && let Ok(stereo) = to_unsigned(&text[9..12], true)
    {
        match stereo {
            0 => res.set_bond_dir(BondDir::None),
            1 => res.set_bond_dir(BondDir::BeginWedge),
            6 => res.set_bond_dir(BondDir::BeginDash),
            3 => {
                res.set_bond_dir(BondDir::EitherDouble);
                res.set_stereo(BondStereo::StereoAny)?;
            }
            4 => res.set_bond_dir(BondDir::Unknown),
            _ => {}
        }
        res.props.set_val(cp::_MolFileBondStereo, stereo);
    }
    if text.len() >= 18
        && &text[15..18] != b"  0"
        && let Ok(topology) = to_int(&text[15..18], true)
        && topology != 0
    {
        if !res.has_query() {
            res = query_bond_from_bond(&res);
        }
        let mut q = make_bond_is_in_ring_query();
        match topology {
            1 => {}
            2 => q.negate = true,
            _ => {
                return Err(fpe(format!(
                    "Unrecognized bond topology specifier: {topology} on line {line}"
                )));
            }
        }
        bond_expand_query(&mut res, q, CompositeQueryType::And, true)?;
    }
    if text.len() >= 21
        && &text[18..21] != b"  0"
        && let Ok(rs) = to_int(&text[18..21], true)
    {
        res.props.set_val(cp::molReactStatus, rs);
    }
    Ok((res, idx1, idx2))
}

/// `ParseOldAtomList`
fn parse_old_atom_list(mol: &mut RWMol, text: &[u8], line: u32) -> Result<()> {
    let f = substr(text, 0, 3)?;
    let idx = cast_err(strip_spaces_and_cast_unsigned(f), || {
        format!("Cannot convert '{}' to int on line {line}", lossy(f))
    })?
    .wrapping_sub(1);
    if idx >= mol.num_atoms() {
        return Err(RDKitError::range("idx", idx, mol.num_atoms()));
    }
    let mut a = query_atom_from_atom(mol.atom(idx))?;
    let mut q: AtomQuery = Query::new(QueryKind::Or);
    q.description = "AtomOr".into();
    match text.get(4).copied().unwrap_or(0) {
        b'T' => q.negate = true,
        b'F' => q.negate = false,
        c => {
            return Err(fpe(format!(
                "Unrecognized atom-list query modifier: '{}' on line {line}",
                c as char
            )));
        }
    }
    let f = substr(text, 9, 1).map_err(|_| {
        fpe(format!(
            "Cannot convert position 9 of '{}' to int on line {line}",
            lossy(text)
        ))
    })?;
    let n_queries = cast_err(to_int(f, true), || {
        format!("Cannot convert '{}' to int on line {line}", lossy(f))
    })?;
    if !(0..=5).contains(&n_queries) {
        return Err(RDKitError::range("nQueries", n_queries, 5));
    }
    for i in 0..n_queries as usize {
        let pos = 11 + i * 4;
        let f = substr(text, pos, 3).map_err(|_| {
            fpe(format!(
                "Cannot convert position {pos} of '{}' to int on line {line}",
                lossy(text)
            ))
        })?;
        let at_num = cast_err(to_int(f, true), || {
            format!("Cannot convert '{}' to int on line {line}", lossy(f))
        })?;
        if !(0..=200).contains(&at_num) {
            return Err(RDKitError::range("atNum", at_num, 200));
        }
        q.children.push(make_atom_num_query(at_num));
        if i == 0 {
            a.set_atomic_num(at_num);
        }
    }
    a.set_query(Some(q));
    a.props.set_val(cp::_MolFileAtomQuery, 1i32);
    mol.replace_atom(idx, &a, false)
}

/// The `nent` count of an `M  XXX` line.
fn count_field(text: &[u8], line: u32, unsigned_cast: bool) -> Result<i64> {
    let f = substr(text, 6, 3)?;
    let msg = || format!("Cannot convert '{}' to int on line {line}", lossy(f));
    if unsigned_cast {
        cast_err(strip_spaces_and_cast_unsigned(f), msg).map(i64::from)
    } else {
        cast_err(to_int(f, true), msg).map(i64::from)
    }
}

/// `ParseChargeLine` / `ParseRadicalLine`
fn parse_charge_or_radical_line(
    mol: &mut RWMol,
    text: &[u8],
    first_call: bool,
    line: u32,
    radical: bool,
) -> Result<()> {
    if first_call {
        for i in 0..mol.num_atoms() {
            mol.atom_mut(i).set_formal_charge(0);
        }
    }
    let nent = count_field(text, line, false)?;
    let mut spos = 9;
    for _ in 0..nent {
        let f = substr(text, spos, 4)?;
        let aid = cast_err(to_int(f, true), || {
            format!("Cannot convert '{}' to int on line {line}", lossy(f))
        })?;
        spos += 4;
        let f = substr(text, spos, 4)?;
        let val = cast_err(to_int(f, true), || {
            format!("Cannot convert '{}' to int on line {line}", lossy(f))
        })?;
        spos += 4;
        let idx = atom_idx(mol, aid as i64)?;
        if radical {
            let n = match val {
                0 => 0,
                1 | 3 => 2,
                2 => 1,
                _ => {
                    return Err(fpe(format!(
                        "Unrecognized radical value {val} for atom {} on line {line}\n",
                        aid - 1
                    )));
                }
            };
            mol.atom_mut(idx).set_num_radical_electrons(n);
        } else {
            mol.atom_mut(idx).set_formal_charge(val);
        }
    }
    Ok(())
}

/// `ParsePXALine`
fn parse_pxa_line(mol: &mut RWMol, text: &[u8], line: u32) -> Result<()> {
    let f = substr(text, 7, 3)?;
    let at = cast_err(strip_spaces_and_cast_unsigned(f), || {
        format!("Cannot convert '{}' to int on line {line}", lossy(f))
    })?;
    let idx = atom_idx(mol, at as i64)?;
    let v = lossy(substr_from(text, 10)?);
    mol.atom_mut(idx).props.set_val("_MolFile_PXA", v);
    Ok(())
}

/// Iterates the `(atom, value)` entries of ISO/SUB/UNS/RBC lines; the
/// value is `None` for a blank field.
fn atom_value_entries(text: &[u8], line: u32) -> Result<Vec<(u32, Option<i32>)>> {
    let nent = count_field(text, line, true)?;
    let mut res = Vec::new();
    let mut spos = 9;
    for _ in 0..nent {
        let f = substr(text, spos, 4)?;
        let aid = cast_err(strip_spaces_and_cast_unsigned(f), || {
            format!("Cannot convert '{}' to int on line {line}", lossy(f))
        })?;
        spos += 4;
        let mut v = None;
        if text.len() >= spos + 4 && &text[spos..spos + 4] != b"    " {
            let f = &text[spos..spos + 4];
            v = Some(cast_err(to_int(f, true), || {
                format!("Cannot convert '{}' to int on line {line}", lossy(f))
            })?);
        }
        spos += 4;
        res.push((aid, v));
    }
    Ok(res)
}

/// `ParseIsotopeLine`
fn parse_isotope_line(mol: &mut RWMol, text: &[u8], line: u32) -> Result<()> {
    let nent = count_field(text, line, true)?;
    let mut spos = 9;
    for _ in 0..nent {
        let f = substr(text, spos, 4)?;
        let aid = cast_err(strip_spaces_and_cast_unsigned(f), || {
            format!("Cannot convert '{}' to int on line {line}", lossy(f))
        })?;
        spos += 4;
        let idx = atom_idx(mol, aid as i64)?;
        if text.len() >= spos + 4 && &text[spos..spos + 4] != b"    " {
            let f = &text[spos..spos + 4];
            let iso = cast_err(to_int(f, true), || {
                format!("Cannot convert '{}' to int on line {line}", lossy(f))
            })?;
            if iso < 0 {
                crate::rd_log_error!(" atom {aid} has a negative isotope value. line:  {line}");
            } else {
                mol.atom_mut(idx).set_isotope(iso as u32);
            }
        }
        spos += 4;
    }
    Ok(())
}

fn and_expand(mol: &mut RWMol, idx: u32, q: AtomQuery) -> Result<()> {
    replace_atom_with_query_atom(mol, idx)?;
    atom_expand_query(mol.atom_mut(idx), q, CompositeQueryType::And, true)
}

/// `ParseSubstitutionCountLine`
fn parse_substitution_count_line(mol: &mut RWMol, text: &[u8], line: u32) -> Result<()> {
    for (aid, count) in atom_value_entries(text, line)? {
        let idx = atom_idx(mol, aid as i64)?;
        let count = count.unwrap_or(0);
        if count == 0 {
            continue;
        }
        let val = match count {
            -1 => 0,
            -2 => mol.atom_degree(idx) as i32,
            1..=5 => count,
            6 => {
                rd_log_warning!(
                    " atom degree query with value 6 found. This will not match degree >6. The MDL spec says it should.  line: {line}"
                );
                6
            }
            _ => {
                return Err(fpe(format!(
                    "Value {count} is not supported as a degree query. line: {line}"
                )));
            }
        };
        and_expand(mol, idx, make_atom_explicit_degree_query(val))?;
    }
    Ok(())
}

/// `ParseUnsaturationLine`
fn parse_unsaturation_line(mol: &mut RWMol, text: &[u8], line: u32) -> Result<()> {
    for (aid, count) in atom_value_entries(text, line)? {
        let idx = atom_idx(mol, aid as i64)?;
        match count.unwrap_or(0) {
            0 => {}
            1 => and_expand(mol, idx, make_atom_unsaturated_query())?,
            c => {
                return Err(fpe(format!(
                    "Value {c} is not supported as an unsaturation query (only 0 and 1 are allowed). line: {line}"
                )));
            }
        }
    }
    Ok(())
}

/// `ParseRingBondCountLine`
fn parse_ring_bond_count_line(mol: &mut RWMol, text: &[u8], line: u32) -> Result<()> {
    for (aid, count) in atom_value_entries(text, line)? {
        let idx = atom_idx(mol, aid as i64)?;
        let count = count.unwrap_or(0);
        if count == 0 {
            continue;
        }
        let q = match count {
            -1 => make_atom_ring_bond_count_query(0),
            -2 => {
                mol.props.set_val(cp::_NeedsQueryScan, 1i32);
                make_atom_ring_bond_count_query(0xDEADBEEFu32 as i32)
            }
            1..=3 => make_atom_ring_bond_count_query(count),
            4 => Query::with(
                QueryKind::LessEqual,
                4,
                AtomFn::RingBondCount,
                "AtomRingBondCount",
            ),
            _ => {
                return Err(fpe(format!(
                    "Value {count} is not supported as a ring-bond count query. line: {line}"
                )));
            }
        };
        and_expand(mol, idx, q)?;
    }
    Ok(())
}

/// `getAndQueries`: the children of a (possibly nested) AND query, or
/// nothing when the tree contains anything else.
fn get_and_queries<'a>(q: &'a AtomQuery, out: &mut Vec<&'a AtomQuery>) -> bool {
    let mut ok = true;
    let mut is_orig = true;
    for c in &q.children {
        if is_orig && q.description != "AtomAnd" {
            ok = false;
            break;
        }
        is_orig = false;
        ok = get_and_queries(c, out);
    }
    // upstream: q is reset to the last child's result; when that equals the
    // original (no children), q itself is collected.
    if q.children.is_empty() {
        out.push(q);
    }
    if !ok {
        out.clear();
    }
    ok
}

/// `ParseNewAtomList`
fn parse_new_atom_list(mol: &mut RWMol, text: &[u8], line: u32) -> Result<()> {
    if text.len() < 15 {
        return Err(fpe(format!("Atom list line too short: '{}'", lossy(text))));
    }
    let f = substr(text, 7, 3)?;
    let idx = cast_err(strip_spaces_and_cast_unsigned(f), || {
        format!("Cannot convert '{}' to int on line {line}", lossy(f))
    })?
    .wrapping_sub(1);
    if idx >= mol.num_atoms() {
        return Err(RDKitError::range("idx", idx, mol.num_atoms()));
    }
    let f = substr(text, 10, 3)?;
    let n_queries = cast_err(to_int(f, true), || {
        format!("Cannot convert '{}' to int on line {line}", lossy(f))
    })?;
    if n_queries == 0 {
        rd_log_warning!("Empty atom list: '{}' on line {line}.", lossy(text));
        return Ok(());
    }
    if n_queries < 0 {
        return Err(fpe(format!(
            "negative length atom list: '{}' on line {line}.\n",
            lossy(text)
        )));
    }
    let a_orig = mol.atom(idx).clone();
    let mut q_orig: Option<AtomQuery> = a_orig.query().cloned();
    let mut a: Option<Atom> = None;
    for i in 0..n_queries as usize {
        let pos = 16 + i * 4;
        if text.len() < pos + 4 {
            return Err(fpe(format!(
                "Atom list line too short: '{}' on line {line}",
                lossy(text)
            )));
        }
        let field = &text[pos..pos + 4];
        let Some(sp) = field.iter().position(|&c| c == b' ') else {
            return Err(RDKitError::Index("basic_string::erase".into()));
        };
        let at_num = PeriodicTable::get().get_atomic_number(&lossy(&field[..sp]))? as i32;
        match a.as_mut() {
            None => {
                let mut na = query_atom_from_atom(&a_orig)?;
                na.set_atomic_num(at_num);
                if q_orig.is_none() {
                    q_orig = na.query().cloned();
                }
                na.set_query(Some(make_atom_num_query(at_num)));
                a = Some(na);
            }
            Some(na) => {
                atom_expand_query(
                    na,
                    make_atom_num_query(at_num),
                    CompositeQueryType::Or,
                    true,
                )?;
                na.set_atomic_num(0);
            }
        }
    }
    let mut a = a.ok_or_else(|| RDKitError::invariant("no atom built"))?;
    if let Some(qo) = &q_orig {
        let mut qs = Vec::new();
        if get_and_queries(qo, &mut qs) {
            for q in qs {
                if q.description != "AtomAtomicNum" {
                    atom_expand_query(&mut a, q.clone(), CompositeQueryType::And, true)?;
                }
            }
        }
    }
    a.props.set_val(cp::_MolFileAtomQuery, 1i32);
    match text[14] {
        b'T' => a.query.as_mut().expect("query").negate = true,
        b'F' => a.query.as_mut().expect("query").negate = false,
        c => {
            return Err(fpe(format!(
                "Unrecognized atom-list query modifier: '{}' on line {line}",
                c as char
            )));
        }
    }
    mol.replace_atom(idx, &a, false)
}

/// `ParseRGroupLabels`
fn parse_rgroup_labels(mol: &mut RWMol, text: &[u8], line: u32) -> Result<()> {
    let n_labels = count_field(text, line, false)?;
    for i in 0..n_labels.max(0) as usize {
        let pos = 10 + i * 8;
        let f = substr(text, pos, 3)?;
        let at = cast_err(strip_spaces_and_cast_unsigned(f), || {
            format!("Cannot convert '{}' to int on line {line}", lossy(f))
        })?;
        let f = substr(text, pos + 4, 3)?;
        let r_label = cast_err(strip_spaces_and_cast_unsigned(f), || {
            format!("Cannot convert '{}' to int on line {line}", lossy(f))
        })?;
        let at_idx = at.wrapping_sub(1);
        if at_idx > mol.num_atoms() {
            return Err(fpe(format!(
                "Attempt to set R group label on nonexistent atom {at_idx} on line {line}"
            )));
        }
        let mut qatom = query_atom_from_atom(mol.get_atom(at_idx)?)?;
        qatom.props.set_val(cp::_MolFileRLabel, r_label);
        qatom.props.set_val(cp::dummyLabel, format!("R{r_label}"));
        if r_label > 0 && r_label < 999 {
            qatom.set_isotope(r_label);
        }
        qatom.set_query(Some(make_atom_null_query()));
        mol.replace_atom(at_idx, &qatom, false)?;
    }
    Ok(())
}

/// `ParseAtomAlias` / `ParseAtomValue`
fn parse_atom_alias_or_value(
    mol: &mut RWMol,
    text: &[u8],
    value: &[u8],
    key: &str,
    line: u32,
) -> Result<()> {
    let f = substr(text, 3, 3)?;
    let idx = cast_err(strip_spaces_and_cast_unsigned(f), || {
        format!("Cannot convert '{}' to int on line {line}", lossy(f))
    })?
    .wrapping_sub(1);
    if idx >= mol.num_atoms() {
        return Err(RDKitError::range("idx", idx, mol.num_atoms()));
    }
    mol.atom_mut(idx).props.set_val(key, lossy(value));
    Ok(())
}

/// `ParseZCHLine`, `ParseHYDLine` and `ParseZBOLine` (Alex Clark's ZBO
/// proposal, JCIM 51:3149-57 (2011)): `M  XXX nn` followed by (index,
/// value) pairs; `apply` gets the zero-based index and the value, which is
/// `None` for a blank field.
fn parse_zero_order_line(
    mol: &mut RWMol,
    text: &[u8],
    line: u32,
    code: &str,
    n_items: u32,
    apply: impl Fn(&mut RWMol, u32, Option<i32>) -> Result<()>,
) -> Result<()> {
    for (id, value) in atom_value_entries(text, line)? {
        if id == 0 || id > n_items {
            return Err(fpe(format!("Bad {code} specification on line {line}")));
        }
        apply(mol, id - 1, value)?;
    }
    Ok(())
}

/// `ParseZCHLine`
fn parse_zch_line(mol: &mut RWMol, text: &[u8], line: u32) -> Result<()> {
    let n = mol.num_atoms();
    parse_zero_order_line(mol, text, line, "ZCH", n, |mol, idx, v| {
        mol.atom_mut(idx).set_formal_charge(v.unwrap_or(0));
        Ok(())
    })
}

/// `ParseHYDLine`
fn parse_hyd_line(mol: &mut RWMol, text: &[u8], line: u32) -> Result<()> {
    let n = mol.num_atoms();
    parse_zero_order_line(mol, text, line, "HYD", n, |mol, idx, v| {
        if let Some(hs) = v.filter(|&hs| hs >= 0) {
            let atom = mol.atom_mut(idx);
            atom.props.set_val("_ZBO_H", true);
            atom.set_num_explicit_hs(hs as u32);
        }
        Ok(())
    })
}

/// `ParseZBOLine`
fn parse_zbo_line(mol: &mut RWMol, text: &[u8], line: u32) -> Result<()> {
    let n = mol.num_bonds();
    parse_zero_order_line(mol, text, line, "ZBO", n, |mol, idx, order| {
        let bond_type = match order.unwrap_or(0) {
            0 => Some(BondType::Zero),
            o => u8::try_from(o).ok().and_then(BondType::from_u8),
        };
        let t =
            bond_type.ok_or_else(|| fpe(format!("Bad ZBO bond order {order:?} on line {line}")))?;
        mol.bond_mut(idx).set_bond_type(t);
        Ok(())
    })
}

/// The 4 character field at `pos` of an `M  APO`/`M  LIN` line.
fn field4(text: &[u8], pos: usize) -> Result<&[u8]> {
    substr(text, pos, 4)
}

fn is_blank_field(f: &[u8]) -> bool {
    f == b"    "
}

fn cannot_convert_int(text: &[u8], pos: usize, line: u32) -> impl FnOnce() -> String {
    let field = lossy(field4(text, pos).unwrap_or_default());
    move || format!("Cannot convert '{field}' to int on line {line}")
}

/// The leading entry count of an `M  APO` / `M  LIN` line.
fn entry_count(text: &[u8], line: u32) -> Result<u32> {
    let field = substr(text, 6, 3)?;
    cast_err(strip_spaces_and_cast_unsigned(field), || {
        format!("Cannot convert '{}' to int on line {line}", lossy(field))
    })
}

/// An optional numeric field: blank (or beyond the end of the line) means
/// unset.
fn optional_field(text: &[u8], pos: usize, line: u32) -> Result<Option<u32>> {
    if text.len() < pos + 4 || is_blank_field(field4(text, pos)?) {
        return Ok(None);
    }
    let field = field4(text, pos)?;
    cast_err(
        strip_spaces_and_cast_unsigned(field),
        cannot_convert_int(text, pos, line),
    )
    .map(Some)
}

/// `ParseMarvinSmartsLine`
fn parse_marvin_smarts_line(mol: &mut RWMol, text: &[u8], line: u32) -> Result<()> {
    const ATOM_NUM_START: usize = 10;
    const SMARTS_START: usize = 15;
    // M  MRV SMA   1 [*;A]
    if !text.starts_with(b"M  MRV SMA") {
        return Ok(());
    }
    let idx_txt = substr(text, ATOM_NUM_START, SMARTS_START - ATOM_NUM_START)?;
    let idx = cast_err(strip_spaces_and_cast_unsigned(idx_txt), || {
        format!(
            "Cannot convert '{}' to an atom index on line {line}",
            lossy(idx_txt)
        )
    })?
    .wrapping_sub(1);
    if idx >= mol.num_atoms() {
        return Err(RDKitError::range("idx", idx as i64, mol.num_atoms()));
    }
    let sma = lossy(substr_from(text, SMARTS_START)?);
    mol.atom_mut(idx).props.set_val(cp::MRV_SMA, sma.clone());
    let query_mol = mol_from_smarts(&sma, &SmartsParserParams::default())
        .ok()
        .flatten()
        .ok_or_else(|| fpe(format!("Cannot parse smarts: '{sma}' on line {line}")))?;
    and_expand(mol, idx, make_recursive_structure_query(query_mol, 0))?;
    mol.atom_mut(idx).props.set_val(cp::_MolFileAtomQuery, 1i32);
    Ok(())
}

/// `ParseAttachPointLine`
fn parse_attach_point_line(mol: &mut RWMol, text: &[u8], line: u32, strict: bool) -> Result<()> {
    if !text.starts_with(b"M  APO") {
        return Err(RDKitError::precondition_at(
            "bad APO line",
            "text.substr(0, 6) == std::string(\"M  APO\")",
            "Code/GraphMol/FileParsers/MolFileParser.cpp",
            947,
        ));
    }
    let nent = entry_count(text, line)?;
    let mut spos = 9;
    for _ in 0..nent {
        let aid = cast_err(
            strip_spaces_and_cast_unsigned(field4(text, spos)?),
            cannot_convert_int(text, spos, line),
        )?;
        spos += 4;
        let mut val = 0i32;
        if text.len() >= spos + 4 && !is_blank_field(field4(text, spos)?) {
            val = cast_err(
                lexical_cast_int(strip_spaces(field4(text, spos)?)),
                cannot_convert_int(text, spos, line),
            )?;
        }
        if aid == 0 || aid > mol.num_atoms() {
            return Err(fpe(format!("Bad APO specification on line {line}")));
        }
        spos += 4;
        let idx = aid - 1;
        if !(0..=3).contains(&val) {
            return Err(fpe(format!(
                "Value {val} from APO specification on line {line} is invalid"
            )));
        }
        if val == 0 {
            continue;
        }
        // this is -1 in v3k mol blocks, so use that
        let val = if val == 3 { -1 } else { val };
        let atom = mol.atom_mut(idx);
        if atom.has_prop(cp::molAttachPoint) {
            let msg = format!(
                "Multiple ATTCHPT values for atom {} on line {line}",
                idx + 1
            );
            if strict {
                return Err(fpe(msg));
            }
            rd_log_warning!("{msg}");
        } else {
            atom.props.set_val(cp::molAttachPoint, val);
        }
    }
    Ok(())
}

/// `ParseLinkNodeLine`: the V2000 format is translated to the V3000 one.
fn parse_link_node_line(mol: &mut RWMol, text: &[u8], line: u32) -> Result<()> {
    if !text.starts_with(b"M  LIN") {
        return Err(RDKitError::precondition_at(
            "bad LIN line",
            "text.substr(0, 6) == std::string(\"M  LIN\")",
            "Code/GraphMol/FileParsers/MolFileParser.cpp",
            1020,
        ));
    }
    let nent = entry_count(text, line)?;
    let n_atoms = mol.num_atoms();
    let mut prop_val = String::new();
    let mut spos = 9;
    for _ in 0..nent {
        let aid = cast_err(
            strip_spaces_and_cast_unsigned(field4(text, spos)?),
            cannot_convert_int(text, spos, line),
        )?;
        if aid == 0 || aid > n_atoms {
            return Err(fpe(format!(
                "LIN specification has bad atom idx on line {line}"
            )));
        }
        spos += 4;
        let repeat_count = optional_field(text, spos, line)?.ok_or_else(|| {
            fpe(format!(
                "LIN specification missing repeat count on line {line}"
            ))
        })?;
        spos += 4;
        if repeat_count < 2 {
            return Err(fpe(format!(
                "LIN specification: repeat count must be >=2 on line {line}"
            )));
        }
        let subst_b = optional_field(text, spos, line)?.unwrap_or(0);
        spos += 4;
        let subst_c = optional_field(text, spos, line)?.unwrap_or(0);
        spos += 4;
        if subst_b == 0 || subst_b > n_atoms || subst_c > n_atoms {
            return Err(fpe(format!(
                "LIN specification has bad substituent idx on line {line}"
            )));
        }
        if !prop_val.is_empty() {
            prop_val.push('|');
        }
        if subst_c == 0 {
            // upstream feeds the substituent C to a format string without a
            // slot for it; `boost::format` throws on the surplus argument,
            // so LIN lines with a single substituent never parse.
            return Err(fpe(
                "boost::too_many_args: format-string referred to less arguments than were passed",
            ));
        }
        prop_val += &format!("1 {repeat_count} 2 {aid} {subst_b} {aid} {subst_c}");
        mol.props.set_val(cp::molFileLinkNodes, prop_val.clone());
    }
    Ok(())
}

fn starts(t: &[u8], p: &[u8]) -> bool {
    t.starts_with(p)
}

fn first(t: &[u8]) -> u8 {
    t.first().copied().unwrap_or(0)
}

/// `ParseMolBlockProperties`; returns whether `M  END` was found.
fn parse_mol_block_properties(
    s: &mut LineStream<'_>,
    line: &mut u32,
    mol: &mut RWMol,
    strict: bool,
) -> Result<bool> {
    let mut temp = s.get_line();
    *line += 1;
    if temp.is_empty() {
        if !strict {
            temp = s.get_line();
            *line += 1;
        } else {
            return Err(fpe(format!(
                "Problems encountered parsing Mol data, unexpected blank line found at line {line}"
            )));
        }
    } else if !matches!(temp[0], b'M' | b'A' | b'V' | b'G' | b'S') {
        parse_old_atom_list(mol, &temp, *line)?;
    }

    let mut sgroups = V2000SGroups::default();
    let mut first_charge_line = true;
    let line_beg = |t: &[u8]| t[..t.len().min(6)].to_vec();
    let mut lb = line_beg(&temp);
    while !s.eof() && !s.fail() && lb != b"M  END" && !starts(&temp, b"$$$$") {
        let l = *line;
        match first(&temp) {
            b'A' => {
                *line += 1;
                let next = s.get_line();
                parse_atom_alias_or_value(mol, &temp, &next, cp::molFileAlias, *line)?;
            }
            b'G' => {
                rd_log_warning!(" deprecated group abbreviation ignored on line {l}");
                *line += 1;
                let _ = s.get_line();
            }
            b'V' => {
                let v = substr_from(&temp, 7)?.to_vec();
                parse_atom_alias_or_value(mol, &temp, &v, cp::molFileValue, l)?;
            }
            _ => match lb.as_slice() {
                b"S  SKP" => {
                    let n = to_int(substr(&temp, 6, 3)?, true)
                        .map_err(|_| fpe(format!("bad SKP line {l}")))?;
                    if n < 0 {
                        return Err(fpe(format!("negative skip value {n} on line {l}")));
                    }
                    for _ in 0..n {
                        *line += 1;
                        let _ = s.get_line();
                    }
                }
                b"M  ALS" => parse_new_atom_list(mol, &temp, l)?,
                b"M  ISO" => parse_isotope_line(mol, &temp, l)?,
                b"M  RGP" => parse_rgroup_labels(mol, &temp, l)?,
                b"M  RBC" => parse_ring_bond_count_line(mol, &temp, l)?,
                b"M  SUB" => parse_substitution_count_line(mol, &temp, l)?,
                b"M  UNS" => parse_unsaturation_line(mol, &temp, l)?,
                b"M  CHG" => {
                    parse_charge_or_radical_line(mol, &temp, first_charge_line, l, false)?;
                    first_charge_line = false;
                }
                b"M  RAD" => {
                    parse_charge_or_radical_line(mol, &temp, first_charge_line, l, true)?;
                    first_charge_line = false;
                }
                b"M  PXA" => parse_pxa_line(mol, &temp, l)?,
                b"M  CRS" => {
                    return Err(fpe(format!(
                        "Unsupported SGroup subtype 'M  CRS' on line {l}"
                    )));
                }
                b if is_sgroup_line(b) => sgroups.parse_line(mol, &temp, l, strict)?,
                b"M  ZBO" => parse_zbo_line(mol, &temp, l)?,
                b"M  ZCH" => parse_zch_line(mol, &temp, l)?,
                b"M  HYD" => parse_hyd_line(mol, &temp, l)?,
                b"M  MRV" => parse_marvin_smarts_line(mol, &temp, l)?,
                b"M  APO" => parse_attach_point_line(mol, &temp, l, strict)?,
                b"M  LIN" => parse_link_node_line(mol, &temp, l)?,
                _ => {}
            },
        }
        *line += 1;
        temp = s.get_line();
        lb = line_beg(&temp);
    }
    let complete = first(&temp) == b'M' && starts(&temp, b"M  END");
    if complete {
        sgroups.finish(mol, strict)?;
    }
    Ok(complete)
}

/// `calculate3dFlag`
pub(super) fn calculate_3d_flag(
    mol: &mut RWMol,
    conf: &Conformer,
    chirality_possible: bool,
) -> bool {
    let marked_3d = mol
        .props
        .get(cp::_3DConf)
        .and_then(RDValue::as_int)
        .unwrap_or(0);
    mol.props.clear_val(cp::_3DConf);
    let nonzero_z = has_non_zero_z_coords(conf);
    if !nonzero_z && marked_3d == 1 {
        if chirality_possible {
            rd_log_warning!(
                "Warning: molecule is tagged as 3D, but all Z coords are zero and 2D stereo markers have been found, marking the mol as 2D."
            );
            return false;
        }
        return true;
    } else if marked_3d == 0 && nonzero_z {
        rd_log_warning!(
            "Warning: molecule is tagged as 2D, but at least one Z coordinate is not zero. Marking the mol as 3D."
        );
        return true;
    }
    nonzero_z
}

/// `ParseV2000CTAB` (atoms, bonds and properties)
fn parse_v2000_ctab(
    s: &mut LineStream<'_>,
    line: &mut u32,
    mol: &mut RWMol,
    chirality_possible: &mut bool,
    n_atoms: u32,
    n_bonds: u32,
    strict: bool,
) -> Result<bool> {
    let mut conf = Conformer::new(n_atoms as usize);
    if n_atoms == 0 {
        conf.set_3d(false);
    } else {
        // ParseMolBlockAtoms
        for i in 1..=n_atoms {
            *line += 1;
            let t = s.get_line();
            if s.eof() {
                return Err(fpe("EOF hit while reading atoms"));
            }
            let (atom, pos) = parse_mol_file_atom_line(&t, *line, strict)?;
            let aid = mol.add_atom(atom, false);
            conf.set_atom_pos(aid, pos)?;
            mol.set_atom_bookmark(aid, i as i32);
        }
    }
    // ParseMolBlockBonds
    for i in 1..=n_bonds {
        *line += 1;
        let t = s.get_line();
        if s.eof() {
            return Err(fpe("EOF hit while reading bonds"));
        }
        let (mut bond, a1, a2) = parse_mol_file_bond_line(&t, *line)?;
        if bond.bond_type() == BondType::Aromatic {
            bond.set_is_aromatic(true);
        }
        if !matches!(bond.bond_dir(), BondDir::None | BondDir::Unknown) {
            *chirality_possible = true;
        }
        if !bond.has_prop(cp::molStereoCare) {
            let care = |a: u32| {
                mol.get_atom(a)
                    .ok()
                    .and_then(|at| at.props.get(cp::molStereoCare).and_then(RDValue::as_int))
            };
            if let (Some(c1), Some(c2)) = (care(a1), care(a2))
                && c1 != 0
                && c2 != 0
            {
                bond.props.set_val(cp::molStereoCare, 1i32);
            }
        }
        bond.set_begin_atom_idx_unchecked(a1);
        bond.set_end_atom_idx_unchecked(a2);
        let bidx = mol.add_bond_object(bond)? - 1;
        mol.set_bond_bookmark(bidx, i as i32);
    }
    let is_3d = calculate_3d_flag(mol, &conf, *chirality_possible);
    conf.set_3d(is_3d);
    mol.add_conformer(conf, true)?;
    parse_mol_block_properties(s, line, mol, strict)
}

/// `ProcessMolProps`
fn process_mol_props(mol: &mut RWMol) -> Result<()> {
    for aidx in 0..mol.num_atoms() {
        let int_prop = |mol: &RWMol, k: &str| mol.atom(aidx).props.get(k).and_then(RDValue::as_int);
        if let Some(mut ival) = int_prop(mol, cp::molSubstCount)
            && ival != 0
        {
            replace_atom_with_query_atom(mol, aidx)?;
            let mut gt = false;
            if ival == -1 {
                ival = 0;
            } else if ival == -2 {
                ival = mol.atom_degree(aidx) as i32;
            } else if ival >= 6 {
                gt = true;
            }
            let q = if !gt {
                make_atom_explicit_degree_query(ival)
            } else {
                Query::with(
                    QueryKind::LessEqual,
                    ival,
                    AtomFn::ExplicitDegree,
                    "less_AtomExplicitDegree",
                )
            };
            atom_expand_query(mol.atom_mut(aidx), q, CompositeQueryType::And, true)?;
        }
        if let Some(ival) = int_prop(mol, cp::molTotValence)
            && ival != 0
            && !mol.atom(aidx).has_prop("_ZBO_H")
        {
            mol.atom_mut(aidx).set_no_implicit(true);
            if ival == 15 || ival == -1 {
                mol.atom_mut(aidx).set_num_explicit_hs(0);
            } else {
                let ev = mol.atom(aidx).valence(ValenceType::Explicit)? as i32;
                if ev > ival {
                    rd_log_warning!(
                        "atom {aidx} has specified valence ({ival}) smaller than the drawn valence {ev}."
                    );
                    mol.atom_mut(aidx).set_num_explicit_hs(0);
                } else {
                    mol.atom_mut(aidx).set_num_explicit_hs((ival - ev) as u32);
                }
            }
        }
        mol.atom_mut(aidx).props.clear_val(cp::molTotValence);
    }
    process_sgroups(mol)
}

/// `FileParserUtils::finishMolProcessing`
pub fn finish_mol_processing(
    res: &mut RWMol,
    chirality_possible: bool,
    params: &MolFileParserParams,
) -> Result<()> {
    res.clear_all_atom_bookmarks();
    res.clear_all_bond_bookmarks();
    if params.expand_attachment_points {
        expand_attachment_points(res, true, true)?;
    }
    for i in 0..res.num_atoms() {
        res.calc_explicit_valence(i, false)?;
    }
    process_mol_props(res)?;

    let (conf_id, is_3d) = {
        let c = res.conformer(-1)?;
        (c.id() as i32, c.is_3d())
    };
    if chirality_possible || is_3d {
        if !is_3d {
            assign_chiral_types_from_bond_dirs(res, conf_id, true)?;
        } else {
            res.update_property_cache(false)?;
            assign_chiral_types_from_3d(res, conf_id, true)?;
        }
    }
    clear_single_bond_dir_flags(res, false);

    if params.sanitize {
        if params.remove_hs {
            sanitize_mol(res, SANITIZE_CLEANUP).map_err(|(_, e)| e)?;
            detect_bond_stereochemistry(res, -1)?;
            remove_hs(res, &RemoveHsParameters::default(), true)?;
        } else {
            sanitize_mol(res, SANITIZE_ALL).map_err(|(_, e)| e)?;
            detect_bond_stereochemistry(res, -1)?;
        }
        assign_stereochemistry(res, true, true, true)?;
    } else {
        detect_bond_stereochemistry(res, -1)?;
    }

    if res.has_prop(cp::_NeedsQueryScan) {
        res.clear_prop(cp::_NeedsQueryScan);
        complete_mol_queries(res, 0xDEADBEEF)?;
    }
    Ok(())
}

/// `MolFromMolDataStream`: `Ok(None)` at the end of the stream.
pub fn mol_from_mol_data_stream(
    s: &mut LineStream<'_>,
    line: &mut u32,
    params: &MolFileParserParams,
) -> Result<Option<RWMol>> {
    let mut chirality_possible = false;
    *line += 1;
    let temp = s.get_line();
    if s.eof() {
        return Ok(None);
    }
    let mut res = RWMol::new();
    res.props.set_val(cp::_Name, lossy(&temp));

    *line += 1;
    let temp = s.get_line();
    res.props.set_val("_MolFileInfo", lossy(&temp));
    if temp.len() >= 22 && matches!(&temp[20..22], b"3d" | b"3D") {
        res.props.set_val(cp::_3DConf, 1i32);
    }
    *line += 1;
    let temp = s.get_line();
    res.props.set_val("_MolFileComments", lossy(&temp));

    *line += 1;
    let temp = s.get_line();
    if temp.len() < 6 {
        return Err(fpe(format!(
            "Counts line too short: '{}' on line{line}",
            lossy(&temp)
        )));
    }
    let f0 = substr(&temp, 0, 3)?;
    let n_atoms = cast_err(to_unsigned(f0, true), || {
        format!(
            "Cannot convert '{}' to unsigned int on line {line}",
            lossy(f0)
        )
    })?;
    let f1 = substr(&temp, 3, 3)?;
    let n_bonds = cast_err(to_unsigned(f1, true), || {
        format!(
            "Cannot convert '{}' to unsigned int on line {line}",
            lossy(f1)
        )
    })?;
    let mut chiral_flag = 0u32;
    // the remaining counts are optional; parse failures are ignored
    let _ = (|| -> Result<()> {
        if temp.len() >= 9 {
            to_unsigned(substr(&temp, 6, 3)?, true)?;
        }
        if temp.len() >= 15 {
            chiral_flag = to_unsigned(substr(&temp, 12, 3)?, true)?;
        }
        Ok(())
    })();

    let mut ctab_version = 2000;
    if temp.len() > 35 {
        if temp.len() < 39 || temp[34] != b'V' {
            let msg = format!("CTAB version string invalid at line {line}");
            if params.strict_parsing {
                return Err(fpe(msg));
            }
            rd_log_warning!("{msg}");
        } else if &temp[34..39] == b"V3000" {
            ctab_version = 3000;
        } else if &temp[34..39] != b"V2000" {
            let msg = format!(
                "Unsupported CTAB version: '{}' at line {line}",
                lossy(&temp[34..39])
            );
            if params.strict_parsing {
                return Err(fpe(msg));
            }
            rd_log_warning!("{msg}");
        } else if params.parsing_scsr_mol {
            return Err(fpe(format!("SCSR Mol files is not V3000 at line{line}")));
        }
    }
    res.props.set_val(cp::_MolFileChiralFlag, chiral_flag);

    let parsed = (|| -> Result<bool> {
        if ctab_version == 2000 {
            parse_v2000_ctab(
                s,
                line,
                &mut res,
                &mut chirality_possible,
                n_atoms,
                n_bonds,
                params.strict_parsing,
            )
        } else {
            if n_atoms != 0 || n_bonds != 0 {
                let msg = format!(
                    "V3000 mol blocks should have 0s in the initial counts line. (line: {line})"
                );
                if params.strict_parsing {
                    return Err(fpe(msg));
                }
                rd_log_warning!("{msg}");
            }
            let (expect_mend, expect_macro_atoms) = if params.parsing_scsr_mol {
                (false, true)
            } else {
                (true, false)
            };
            super::mol_file_v3000::parse_v3000_ctab(
                s,
                line,
                &mut res,
                &mut chirality_possible,
                params.strict_parsing,
                expect_mend,
                expect_macro_atoms,
            )
        }
    })();
    let file_complete = match parsed {
        Ok(complete) => complete,
        Err(RDKitError::UnhandledFeature(msg)) => {
            // an unsupported feature: skip the rest of the molecule
            rd_log_error!(" Unhandled CTAB feature: '{msg}'. Molecule skipped.");
            let mut temp = if s.eof() { Vec::new() } else { s.get_line() };
            *line += 1;
            while !s.eof()
                && !s.fail()
                && !temp.starts_with(b"M  END")
                && !temp.starts_with(b"$$$$")
            {
                temp = s.get_line();
                *line += 1;
            }
            if s.eof() && !temp.starts_with(b"M  END") && !temp.starts_with(b"$$$$") {
                return Err(fpe(format!(
                    "Problems encountered parsing Mol data, M  END missing around line {line}"
                )));
            }
            return Ok(None);
        }
        Err(e) => return Err(e),
    };
    if !file_complete {
        return Err(fpe(format!(
            "Problems encountered parsing Mol data, M  END missing around line {line}"
        )));
    }
    finish_mol_processing(&mut res, chirality_possible, params)?;
    Ok(Some(res))
}

/// `MolFromMolBlock`
pub fn mol_from_mol_block(mol_block: &str, params: &MolFileParserParams) -> Result<Option<RWMol>> {
    mol_from_mol_block_bytes(mol_block.as_bytes(), params)
}

/// `MolFromMolBlock` on raw bytes.
pub fn mol_from_mol_block_bytes(
    mol_block: &[u8],
    params: &MolFileParserParams,
) -> Result<Option<RWMol>> {
    let mut s = LineStream::new(mol_block);
    let mut line = 0;
    mol_from_mol_data_stream(&mut s, &mut line, params)
}

/// `MolFromMolFile`
pub fn mol_from_mol_file(file_name: &str, params: &MolFileParserParams) -> Result<Option<RWMol>> {
    let data = std::fs::read(file_name)
        .map_err(|_| RDKitError::BadFile(format!("Bad input file {file_name}")))?;
    mol_from_mol_block_bytes(&data, params)
}
