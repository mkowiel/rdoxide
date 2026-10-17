//! Tripos Mol2 file parser.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/Mol2FileParser.cpp`
//! (2026.09.1).
//!
//! Copyright (c) 2008-2023, Novartis Institutes for BioMedical Research Inc.
//! and other RDKit contributors
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::utils::{LineStream, lexical_cast_int, lexical_cast_unsigned, strip_spaces};
use crate::geometry::Point3D;
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::BondType;
use crate::graphmol::chirality::from_coords::{
    assign_chiral_types_from_3d, detect_bond_stereochemistry,
};
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::find_rings::find_sssr;
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, remove_hs};
use crate::graphmol::mol_ops::sanitize_flags::{
    SANITIZE_ALL, SANITIZE_CLEANUP, SANITIZE_CLEANUP_ORGANOMETALLICS,
};
use crate::graphmol::mol_ops::{clean_up, sanitize_mol};
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::queries::{CompositeQueryType, make_atom_null_query, make_atom_num_query};
use crate::graphmol::query_ops::atom_expand_query;
use crate::graphmol::romol::RWMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::{RDKitError, Result, rd_log_warning};

/// `Mol2ParserParams`
#[derive(Debug, Clone, Copy)]
pub struct Mol2ParserParams {
    pub sanitize: bool,
    pub remove_hs: bool,
    pub cleanup_substructures: bool,
}

impl Default for Mol2ParserParams {
    fn default() -> Self {
        Mol2ParserParams {
            sanitize: true,
            remove_hs: true,
            cleanup_substructures: true,
        }
    }
}

fn file_parse_error(message: impl Into<String>) -> RDKitError {
    RDKitError::FileParse(message.into())
}

/// `boost::char_separator<char>(" \t\n")`: empty tokens are dropped.
fn tokens(s: &str) -> Vec<&str> {
    s.split([' ', '\t', '\n'])
        .filter(|t| !t.is_empty())
        .collect()
}

fn line_string(stream: &mut LineStream<'_>) -> String {
    String::from_utf8_lossy(&stream.get_line()).into_owned()
}

fn cast_double(s: &str) -> Result<f64> {
    String::from_utf8_lossy(strip_spaces(s.as_bytes()))
        .parse::<f64>()
        .map_err(|_| RDKitError::Value("bad lexical cast".into()))
}

fn atom_prop(mol: &RWMol, idx: u32, name: &str) -> String {
    mol.atom(idx)
        .get_prop(name)
        .map(|v| v.to_rd_string())
        .unwrap_or_default()
}

fn mol_name(mol: &RWMol) -> String {
    mol.get_prop(cp::_Name)
        .map(|v| v.to_rd_string())
        .unwrap_or_default()
}

fn set_bond(mol: &mut RWMol, a: u32, b: u32, t: BondType, aromatic: bool) {
    if let Some(bi) = mol.bond_idx_between(a, b) {
        let bond = mol.bond_mut(bi);
        bond.set_bond_type(t);
        bond.set_is_aromatic(aromatic);
    }
}

fn bond_type_between(mol: &RWMol, a: u32, b: u32) -> Option<BondType> {
    mol.bond_idx_between(a, b)
        .map(|bi| mol.bond(bi).bond_type())
}

fn fix_nitro_substructure_and_charge(mol: &mut RWMol, at_idx: u32) {
    let mut n_o_dbl = 0;
    let mut to_mod = 0;
    let nbrs: Vec<u32> = mol.atom_neighbors(at_idx).collect();
    for n in nbrs {
        if mol.atom(n).atomic_num() == 8
            && bond_type_between(mol, at_idx, n) == Some(BondType::Double)
        {
            n_o_dbl += 1;
            to_mod = n;
        }
    }
    if n_o_dbl == 2 {
        set_bond(
            mol,
            at_idx,
            to_mod,
            BondType::Single,
            mol_bond_arom(mol, at_idx, to_mod),
        );
        mol.atom_mut(at_idx).set_formal_charge(1);
        mol.atom_mut(to_mod).set_formal_charge(-1);
    }
}

/// Aromatic flag of the bond between two atoms (unchanged by `setBondType`).
fn mol_bond_arom(mol: &RWMol, a: u32, b: u32) -> bool {
    mol.bond_idx_between(a, b)
        .is_some_and(|bi| mol.bond(bi).is_aromatic())
}

/// `readFormalChargesFromAttr`
fn read_formal_charges_from_attr(stream: &mut LineStream<'_>, mol: &mut RWMol) -> Result<()> {
    if stream.eof() {
        return Err(RDKitError::precondition("inStream is at eof"));
    }
    let mut temp = line_string(stream);
    if stream.eof() {
        return Err(file_parse_error("premature EOF in readFormalCharges"));
    }
    loop {
        let toks = tokens(&temp);
        let bad = || file_parse_error("Cannot process mol2 UnityAtomAttr.");
        if toks.len() < 2 {
            return Err(bad());
        }
        let atom_idx = lexical_cast_unsigned(toks[0].as_bytes()).map_err(|_| bad())?;
        let n_attr = lexical_cast_unsigned(toks[1].as_bytes()).map_err(|_| bad())?;
        for _ in 0..n_attr {
            let line = line_string(stream);
            if stream.eof() {
                return Err(file_parse_error("premature EOF in readFormalCharges"));
            }
            let t = tokens(&line);
            if t.first() == Some(&"AtomExpr") {
                let Some(expr) = t.get(1) else {
                    return Err(file_parse_error("Cannot process mol2 formal charge."));
                };
                if !expr.contains('=') {
                    let charge = lexical_cast_int(expr.as_bytes())
                        .map_err(|_| file_parse_error("Cannot process mol2 formal charge."))?;
                    if atom_idx == 0 || atom_idx > mol.num_atoms() {
                        return Err(RDKitError::precondition("Atom index out of range"));
                    }
                    mol.atom_mut(atom_idx - 1).set_formal_charge(charge);
                }
            }
        }
        if stream.eof() {
            break;
        }
        temp = line_string(stream);
        if temp.is_empty() || temp.starts_with('@') || temp.starts_with('#') {
            break;
        }
    }
    Ok(())
}

/// `guessFormalCharges`
fn guess_formal_charges(mol: &mut RWMol) -> Result<()> {
    let table = PeriodicTable::get();
    for idx in 0..mol.num_atoms() {
        let at = mol.atom(idx);
        if at.formal_charge() != 0 || at.symbol() == "C" || at.has_query() {
            continue;
        }
        let mut n_arom = 0;
        let mut accum = 0.0;
        for bnd in mol.atom_bonds(idx) {
            accum += bnd.valence_contrib(mol, idx);
            if bnd.bond_type() == BondType::Aromatic {
                n_arom += 1;
            }
        }
        if !mol.ring_info().is_sssr_or_better() {
            find_sssr(mol, false, false)?;
        }
        let at = mol.atom(idx);
        let tatt = atom_prop(mol, idx, cp::_TriposAtomType);
        if at.is_aromatic()
            && !tatt.contains("ar")
            && mol.ring_info().is_atom_in_ring_of_size(idx, 5)?
        {
            continue;
        }
        if n_arom == 3 && tatt == "N.ar" {
            rd_log_warning!(
                "{}: warning - aromatic N with 3 aromatic bonds - skipping charge guess for this atom",
                mol_name(mol)
            );
            continue;
        }
        let exp_val = (accum + 0.1).round() as i32;
        let valens = table.get_valence_list(at.atomic_num())?;
        let first = valens.first().copied().unwrap_or(0);
        let n_elec = table.get_n_outer_elecs(at.atomic_num())?;
        let mut assign = if n_elec >= 4 {
            exp_val - first
        } else {
            first - exp_val
        };
        if assign > 0 && n_elec >= 4 {
            for &vi in valens {
                assign = exp_val - vi;
                if vi <= exp_val && assign < 2 {
                    break;
                }
            }
        }
        if assign != 0 {
            let aromatic = at.is_aromatic();
            let is_n = at.symbol() == "N";
            let charge = if aromatic && assign.abs() > 1 {
                i32::from(assign > 0) - i32::from(assign < 0)
            } else {
                assign
            };
            mol.atom_mut(idx).set_formal_charge(charge);
            if assign == 2 && exp_val == 5 && is_n {
                fix_nitro_substructure_and_charge(mol, idx);
            }
        }
    }
    Ok(())
}

/// `chkNoHNeighbNOx`: the number of H neighbours of `at_idx`; `to_mod` is set
/// to `at_idx` when it is part of an N-oxide.
fn chk_no_h_neighb_nox(mol: &RWMol, at_idx: u32, to_mod: &mut i64) -> u32 {
    let mut n_h = 0;
    for n in mol.atom_neighbors(at_idx) {
        if mol.atom(n).atomic_num() == 1 {
            n_h += 1;
        } else if mol.atom(n).atomic_num() == 8 && mol.atom_degree(n) == 1 {
            *to_mod = i64::from(at_idx);
        }
    }
    n_h
}

/// The guanidinium / amidine part of `cleanUpMol2Substructures`.
fn clean_up_c_cat(mol: &mut RWMol, idx: u32, fixed: &mut [bool]) -> Result<bool> {
    fixed[idx as usize] = true;
    let nbrs: Vec<u32> = mol.atom_neighbors(idx).collect();
    let mut to_mod: i64 = -1;
    let n_n = nbrs
        .iter()
        .filter(|&&n| mol.atom(n).symbol() == "N")
        .count();
    if !(2..=3).contains(&n_n) {
        rd_log_warning!(
            "{}: Error - C.Cat with bad number of N neighbors.",
            mol_name(mol)
        );
        return Ok(false);
    }
    if n_n == 2 {
        let mut first = nbrs[0];
        let mut second = nbrs[0];
        let mut found_first = false;
        for &n in &nbrs {
            if mol.atom(n).symbol() == "N" {
                set_bond(mol, idx, n, BondType::Single, false);
                mol.atom_mut(n).set_is_aromatic(false);
                if fixed[n as usize] {
                    rd_log_warning!(
                        "{}: warning - charged amidine and isFixed atom.",
                        mol_name(mol)
                    );
                }
                fixed[n as usize] = true;
                if found_first {
                    second = n;
                } else {
                    first = n;
                    found_first = true;
                }
            }
        }
        let h1 = chk_no_h_neighb_nox(mol, first, &mut to_mod);
        let h2 = chk_no_h_neighb_nox(mol, second, &mut to_mod);
        if to_mod < 0 {
            if h1 != h2 {
                to_mod = i64::from(if h1 > h2 { first } else { second });
            } else {
                if !mol.ring_info().is_sssr_or_better() {
                    find_sssr(mol, false, false)?;
                }
                let r1 = mol.ring_info().num_atom_rings(first)?;
                let r2 = mol.ring_info().num_atom_rings(second)?;
                to_mod = i64::from(if r1 > r2 { first } else { second });
            }
        }
        let to_mod = to_mod as u32;
        set_bond(mol, idx, to_mod, BondType::Double, false);
        mol.atom_mut(to_mod).set_formal_charge(1);
        mol.atom_mut(idx).set_is_aromatic(false);
    } else {
        let mut lowest = 100;
        for &n in &nbrs {
            if !fixed[n as usize] {
                let mut hvy = 0;
                for nn in mol.atom_neighbors(n) {
                    if mol.atom(nn).atomic_num() > 1 {
                        hvy += if atom_prop(mol, nn, cp::_TriposAtomType) == "C.cat" {
                            2
                        } else {
                            1
                        };
                    }
                }
                if hvy < lowest {
                    to_mod = i64::from(n);
                    lowest = hvy;
                }
                set_bond(mol, idx, n, BondType::Single, false);
                mol.atom_mut(n).set_is_aromatic(false);
                fixed[n as usize] = true;
            } else {
                set_bond(mol, idx, n, BondType::Single, false);
            }
        }
        if to_mod < 0 {
            return Err(RDKitError::precondition("no nitrogen found to modify"));
        }
        let to_mod = to_mod as u32;
        set_bond(mol, idx, to_mod, BondType::Double, false);
        mol.atom_mut(to_mod).set_formal_charge(1);
        mol.atom_mut(idx).set_is_aromatic(false);
    }
    Ok(true)
}

/// `cleanUpMol2Substructures`: `false` when the molecule cannot be fixed.
fn clean_up_mol2_substructures(mol: &mut RWMol) -> Result<bool> {
    let mut fixed = vec![false; mol.num_atoms() as usize];
    for idx in 0..mol.num_atoms() {
        if fixed[idx as usize] {
            continue;
        }
        let tat = atom_prop(mol, idx, cp::_TriposAtomType);
        if tat == "N.4" {
            mol.atom_mut(idx).set_formal_charge(1);
        } else if tat == "O.co2" {
            if mol.atom_degree(idx) != 1 {
                rd_log_warning!("Warning - O.co2 with degree >1.");
                return Ok(false);
            }
            let nbr = mol.atom_neighbors(idx).next().unwrap_or(idx);
            let tatt = atom_prop(mol, nbr, cp::_TriposAtomType);
            if tatt == "P.3" {
                // phosphates: keep the first bond to O.co2 double, the rest single
                set_bond(mol, idx, nbr, BondType::Double, false);
                mol.atom_mut(idx).set_is_aromatic(false);
                fixed[idx as usize] = true;
                let onbrs: Vec<u32> = mol.atom_neighbors(nbr).collect();
                for o in onbrs {
                    if mol.atom(o).atomic_num() == 8
                        && !fixed[o as usize]
                        && atom_prop(mol, o, cp::_TriposAtomType) == "O.co2"
                    {
                        set_bond(mol, nbr, o, BondType::Single, false);
                        mol.atom_mut(o).set_formal_charge(-1);
                        mol.atom_mut(o).set_is_aromatic(false);
                        fixed[o as usize] = true;
                    }
                }
                mol.atom_mut(nbr).set_is_aromatic(false);
                fixed[nbr as usize] = true;
            } else if tatt == "C.2" || tatt == "S.o2" {
                if !fixed[nbr as usize] {
                    // the first occurrence is negatively charged with a single bond
                    set_bond(mol, idx, nbr, BondType::Single, false);
                    mol.atom_mut(idx).set_formal_charge(-1);
                    mol.atom_mut(idx).set_is_aromatic(false);
                    mol.atom_mut(nbr).set_is_aromatic(false);
                    fixed[idx as usize] = true;
                    fixed[nbr as usize] = true;
                } else {
                    set_bond(mol, idx, nbr, BondType::Double, false);
                    mol.atom_mut(idx).set_is_aromatic(false);
                    fixed[idx as usize] = true;
                }
            } else {
                rd_log_warning!(
                    "{}: warning - O.co2 with non C.2 or S.o2 neighbor.",
                    mol_name(mol)
                );
                return Ok(false);
            }
        } else if tat == "C.cat" && !clean_up_c_cat(mol, idx, &mut fixed)? {
            return Ok(false);
        }
    }
    Ok(true)
}

fn or_query(atom: &mut Atom, nums: &[i32]) -> Result<()> {
    for &n in nums {
        atom_expand_query(atom, make_atom_num_query(n), CompositeQueryType::Or, true)?;
    }
    Ok(())
}

/// `ParseMol2FileAtomLine`: `None` for lone pairs.
fn parse_atom_line(line: &str, pos: &mut Point3D) -> Result<Option<Atom>> {
    let toks = tokens(line);
    if toks.is_empty() {
        return Err(file_parse_error("no info in mol2 atom line"));
    }
    let premature = || file_parse_error("premature end of mol2 atom line");
    // skip TriposAtomId
    let name = *toks.get(1).ok_or_else(premature)?;
    let coord = |i: usize| -> Result<f64> {
        let t = toks.get(i).ok_or_else(premature)?;
        cast_double(t).map_err(|_| file_parse_error("Cannot process mol2 coordinates."))
    };
    pos.x = coord(2)?;
    pos.y = coord(3)?;
    pos.z = coord(4)?;
    let tat = *toks.get(5).ok_or_else(premature)?;
    let symb = tat.split('.').next().unwrap_or(tat);

    let mut res;
    match symb {
        "LP" => return Ok(None),
        "ANY" | "Du" => {
            res = Atom::new(0);
            res.set_query(Some(make_atom_null_query()));
        }
        "HEV" => {
            res = Atom::new(1);
            let mut q = make_atom_num_query(1);
            q.negate = true;
            res.set_query(Some(q));
        }
        "HET" => {
            res = Atom::new(7);
            res.set_query(Some(make_atom_num_query(7)));
            or_query(&mut res, &[8, 15, 16])?;
        }
        "HAL" => {
            res = Atom::new(9);
            res.set_query(Some(make_atom_num_query(9)));
            or_query(&mut res, &[17, 35, 53])?;
        }
        _ => {
            res = Atom::new(PeriodicTable::get().get_atomic_number(symb)?);
        }
    }
    res.set_prop("_TriposAtomName", name.to_string(), false)?;
    res.set_prop(cp::_TriposAtomType, tat.to_string(), false)?;
    // no implicit hydrogens for mol2 files
    res.set_no_implicit(true);

    // subst_id and subst_name are skipped; the partial charge follows
    if let Some(charge) = toks.get(8) {
        res.set_prop("_TriposPartialCharge", (*charge).to_string(), false)?;
    }
    Ok(Some(res))
}

/// `ParseMol2FileBondLine`: `None` for bonds that are skipped.
fn parse_bond_line(line: &str, corresp: &[Option<u32>]) -> Result<Option<(u32, u32, BondType)>> {
    let toks = tokens(line);
    let no_info = || file_parse_error("no info in mol2 bond line");
    if toks.is_empty() {
        return Err(no_info());
    }
    let cast = |i: usize| -> Result<u32> {
        let t = toks.get(i).ok_or_else(no_info)?;
        lexical_cast_unsigned(t.as_bytes())
            .map_err(|_| file_parse_error("Cannot process mol2 bonds."))
    };
    let idx1 = cast(1)?;
    let idx2 = cast(2)?;
    let tb = *toks.get(3).ok_or_else(no_info)?;
    if idx1 == 0 || idx2 == 0 {
        return Err(file_parse_error("Mol2 bond index starts at 1, not 0."));
    }
    let (i1, i2) = ((idx1 - 1) as usize, (idx2 - 1) as usize);
    if i1 >= corresp.len() || i2 >= corresp.len() {
        return Err(file_parse_error("index mismatch"));
    }
    let (Some(a), Some(b)) = (corresp[i1], corresp[i2]) else {
        return Ok(None);
    };
    let t = match tb {
        "1" | "am" => BondType::Single,
        "2" => BondType::Double,
        "3" => BondType::Triple,
        "ar" => BondType::Aromatic,
        "du" | "un" => BondType::Unspecified,
        _ => {
            rd_log_warning!("Warning - unsupported bond type: {tb} ignored!");
            return Ok(None);
        }
    };
    Ok(Some((a, b, t)))
}

fn parse_atom_block(
    stream: &mut LineStream<'_>,
    mol: &mut RWMol,
    n_atoms: u32,
    corresp: &mut [Option<u32>],
) -> Result<()> {
    if stream.eof() {
        return Err(RDKitError::precondition("inStream is at eof"));
    }
    let mut n_lp = 0;
    let mut has_h = false;
    let mut positions = Vec::with_capacity(n_atoms as usize);
    for slot in corresp.iter_mut() {
        let line = line_string(stream);
        if stream.eof() {
            return Err(file_parse_error("premature EOF"));
        }
        let mut pos = Point3D::default();
        match parse_atom_line(&line, &mut pos)? {
            Some(atom) => {
                if atom.symbol() == "H" {
                    has_h = true;
                }
                *slot = Some(mol.add_atom(atom, false));
                positions.push(pos);
            }
            None => n_lp += 1,
        }
    }
    // formal charge estimation needs explicit hydrogens
    if !has_h {
        rd_log_warning!(
            "{}: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.",
            mol_name(mol)
        );
    }
    let mut conf = Conformer::new((n_atoms - n_lp) as usize);
    for (i, p) in positions.into_iter().enumerate() {
        conf.set_atom_pos(i as u32, p)?;
    }
    mol.add_conformer(conf, true)?;
    Ok(())
}

fn parse_bond_block(
    stream: &mut LineStream<'_>,
    mol: &mut RWMol,
    n_bonds: u32,
    corresp: &[Option<u32>],
) -> Result<()> {
    if stream.eof() {
        return Err(RDKitError::precondition("inStream is at eof"));
    }
    for _ in 0..n_bonds {
        let line = line_string(stream);
        if stream.eof() {
            return Err(file_parse_error("premature EOF"));
        }
        if let Some((a, b, t)) = parse_bond_line(&line, corresp)? {
            // adding an aromatic bond also flags both of its atoms
            mol.add_bond(a, b, t)?;
        }
    }
    Ok(())
}

/// Sanitization and stereo perception of a freshly parsed molecule.
fn finish_mol(mol: &mut RWMol, params: &Mol2ParserParams) -> Result<()> {
    clean_up(mol)?;
    // the organometallics cleanup does not fit mol2 semantics
    let flags = SANITIZE_ALL ^ SANITIZE_CLEANUP_ORGANOMETALLICS;
    let sanitize = |mol: &mut RWMol, ops: u32| -> Result<()> {
        sanitize_mol(mol, ops).map_err(|(_, e)| {
            rd_log_warning!("sanitize {}: ", mol_name(mol));
            e
        })
    };
    if params.remove_hs {
        // bond stereo must be detected before the Hs go, so cleanup runs first
        sanitize(mol, SANITIZE_CLEANUP)?;
        detect_bond_stereochemistry(mol, -1)?;
        remove_hs(mol, &RemoveHsParameters::default(), false)?;
        sanitize(mol, flags)?;
    } else {
        sanitize(mol, flags)?;
        detect_bond_stereochemistry(mol, -1)?;
    }
    mol.update_property_cache(false)?;
    assign_stereochemistry(mol, true, true, false)
}

/// `MolFromMol2DataStream` / `MolFromMol2Block`; `Ok(None)` when the
/// substructure cleanup fails.
pub fn mol_from_mol2_block(block: &str, params: &Mol2ParserParams) -> Result<Option<RWMol>> {
    let mut stream = LineStream::new(block.as_bytes());
    let (mut mol_start, mut atom_start, mut bond_start, mut charge_start) = (0, 0, 0, 0);
    while !stream.eof() && !stream.fail() {
        let line = line_string(&mut stream);
        if stream.eof() {
            break;
        }
        if line.starts_with('@') {
            let first = tokens(&line).first().copied().unwrap_or("");
            let pos = stream.tellg().unwrap_or(0);
            match first {
                "@<TRIPOS>MOLECULE" => {
                    if mol_start == 0 {
                        mol_start = pos;
                    } else {
                        break;
                    }
                }
                "@<TRIPOS>ATOM" => atom_start = pos,
                "@<TRIPOS>BOND" => bond_start = pos,
                "@<TRIPOS>UNITY_ATOM_ATTR" => charge_start = pos,
                _ => {}
            }
        }
    }
    if mol_start == 0 {
        return Err(file_parse_error("No MOLECULE block found in Mol2 data"));
    }
    if atom_start == 0 {
        return Err(file_parse_error("No ATOM block found in Mol2 data"));
    }

    stream.seek(mol_start);
    let mut mol = RWMol::new();
    let name = line_string(&mut stream);
    mol.set_prop(cp::_Name, name.trim_end().to_string(), false)?;

    let counts = line_string(&mut stream);
    let toks = tokens(&counts);
    if toks.is_empty() {
        return Err(file_parse_error("Empty counts line"));
    }
    let cannot = |t: &str| file_parse_error(format!("Cannot convert {t} to unsigned int"));
    let n_atoms = lexical_cast_unsigned(toks[0].as_bytes()).map_err(|_| cannot(toks[0]))?;
    let n_bonds = match toks.get(1) {
        Some(t) => lexical_cast_unsigned(t.as_bytes()).map_err(|_| cannot(t))?,
        None => 0,
    };
    if n_atoms == 0 {
        return Err(file_parse_error("molecule has no atoms"));
    }
    line_string(&mut stream); // mol_type - ignored
    let charge_type = line_string(&mut stream);
    mol.set_prop("_TriposChargeType", charge_type.trim().to_string(), false)?;

    let mut corresp = vec![None; n_atoms as usize];
    stream.seek(atom_start);
    parse_atom_block(&mut stream, &mut mol, n_atoms, &mut corresp)?;
    if n_bonds != 0 {
        stream.seek(bond_start);
        parse_bond_block(&mut stream, &mut mol, n_bonds, &corresp)?;
    }

    if charge_start == 0 {
        if params.cleanup_substructures && !clean_up_mol2_substructures(&mut mol)? {
            return Ok(None);
        }
        // mol2 has no formal charges, so guess them from the valences
        guess_formal_charges(&mut mol)?;
    } else {
        stream.seek(charge_start);
        read_formal_charges_from_attr(&mut stream, &mut mol)?;
    }

    // chirality comes from the 3D coordinates and must be set before the Hs go
    assign_chiral_types_from_3d(&mut mol, -1, true)?;
    if params.sanitize {
        finish_mol(&mut mol, params)?;
    }
    Ok(Some(mol))
}

/// `MolFromMol2File`
pub fn mol_from_mol2_file(file_name: &str, params: &Mol2ParserParams) -> Result<Option<RWMol>> {
    let data = std::fs::read(file_name)
        .map_err(|_| RDKitError::BadFile(format!("Bad input file {file_name}")))?;
    if data.is_empty() {
        return Ok(None);
    }
    mol_from_mol2_block(&String::from_utf8_lossy(&data), params)
}
