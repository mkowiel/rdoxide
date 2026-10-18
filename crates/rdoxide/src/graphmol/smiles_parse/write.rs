//! SMILES writer: `MolToSmiles`, `MolFragmentToSmiles`,
//! `MolToRandomSmilesVect`, `GetAtomSmiles`, `GetBondSmiles`.
//!
//! Ported from RDKit `Code/GraphMol/SmilesParse/SmilesWrite.{h,cpp}`
//! (2026.09.1). The CXSMILES extensions are written by `cx_write`.
//!
//! Copyright (C) 2002-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;

use super::ops::check_chiral_permutation;
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::canon::{
    AtomColor, FragmentOptions, MolStackElem, canonicalize_enhanced_stereo, canonicalize_fragment,
};
use crate::graphmol::chirality::get_allow_nontetrahedral_chirality;
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::kekulize::{DEFAULT_MAX_BACKTRACKS, kekulize, kekulize_fragment};
use crate::graphmol::mol_ops::{get_mol_frags_as_mols, get_mol_frags_atoms, is_metal};
use crate::graphmol::new_canon::{
    FragmentRankOptions, RankOptions, rank_fragment_atoms, rank_mol_atoms,
};
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::random::with_global_generator;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// `SmilesWriteParams`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesWriteParams {
    pub do_isomeric_smiles: bool,
    pub do_kekule: bool,
    pub canonical: bool,
    pub clean_stereo: bool,
    pub all_bonds_explicit: bool,
    pub all_hs_explicit: bool,
    pub do_random: bool,
    pub rooted_at_atom: i32,
    pub include_dative_bonds: bool,
    pub ignore_atom_map_numbers: bool,
}

impl Default for SmilesWriteParams {
    fn default() -> Self {
        SmilesWriteParams {
            do_isomeric_smiles: true,
            do_kekule: false,
            canonical: true,
            clean_stereo: true,
            all_bonds_explicit: false,
            all_hs_explicit: false,
            do_random: false,
            rooted_at_atom: -1,
            include_dative_bonds: true,
            ignore_atom_map_numbers: false,
        }
    }
}

/// `SmilesWrite::inOrganicSubset`
pub fn in_organic_subset(atomic_number: i32) -> bool {
    matches!(
        atomic_number,
        0 | 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53
    )
}

fn string_prop(p: &impl RDProps, key: &str) -> Option<String> {
    match p.get_prop(key) {
        Some(RDValue::String(s)) => Some(s.clone()),
        Some(v) => Some(v.to_rd_string()),
        None => None,
    }
}

/// anonymous-namespace `getAtomChiralityInfo`
fn atom_chirality_info(mol: &ROMol, idx: u32) -> Result<String> {
    let atom = mol.atom(idx);
    let mut s = match atom.chiral_tag() {
        ChiralType::TetrahedralCW => "@@".to_string(),
        ChiralType::TetrahedralCCW => "@".to_string(),
        _ => String::new(),
    };
    if s.is_empty() && get_allow_nontetrahedral_chirality() {
        s = match atom.chiral_tag() {
            ChiralType::SquarePlanar => "@SP",
            ChiralType::TrigonalBipyramidal => "@TB",
            ChiralType::Octahedral => "@OH",
            _ => "",
        }
        .to_string();
        if !s.is_empty() {
            let perm = atom
                .get_prop(cp::_chiralPermutation)
                .and_then(RDValue::as_int);
            if let Some(p) = perm {
                if atom.chiral_tag() > ChiralType::Other
                    && !check_chiral_permutation(atom.chiral_tag(), p)
                {
                    return Err(RDKitError::Value("bad chirality spec".into()));
                }
                if p != 0 {
                    s += &p.to_string();
                }
            }
        }
    }
    Ok(s)
}

/// anonymous-namespace `atomNeedsBracket`
fn atom_needs_bracket(
    mol: &ROMol,
    idx: u32,
    at_string: &str,
    params: &SmilesWriteParams,
) -> Result<bool> {
    let atom = mol.atom(idx);
    let num = atom.atomic_num();
    if !in_organic_subset(num as i32) {
        return Ok(true);
    }
    if atom.formal_charge() != 0 {
        return Ok(true);
    }
    if params.do_isomeric_smiles && (atom.isotope() != 0 || !at_string.is_empty()) {
        return Ok(true);
    }
    if atom.has_prop(cp::molAtomMapNumber) {
        return Ok(true);
    }
    let default_vs = PeriodicTable::get().get_valence_list(num)?;
    let total_valence = atom.total_valence()? as i32;
    let non_standard = if atom.num_radical_electrons() != 0 {
        true
    } else if (num == 7 || num == 15) && atom.is_aromatic() && atom.num_explicit_hs() != 0 {
        true
    } else {
        total_valence != default_vs[0] && mol.total_num_hs(idx, false)? != 0
    };
    if non_standard {
        return Ok(true);
    }
    if atom.has_owning_mol() {
        for nbr in mol.atom_neighbors(idx) {
            if is_metal(mol.atom(nbr).atomic_num()) {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// `SmilesWrite::GetAtomSmiles(atom, params)`
pub fn get_atom_smiles(mol: &ROMol, idx: u32, params: &SmilesWriteParams) -> Result<String> {
    let atom = mol.atom(idx);
    let fc = atom.formal_charge();
    let num = atom.atomic_num();
    let isotope = atom.isotope();
    let custom = string_prop(atom, cp::smilesSymbol);
    let has_custom = custom.is_some();
    let mut symb = match custom {
        Some(s) => s,
        None => PeriodicTable::get().get_element_symbol(num)?.to_string(),
    };
    let mut at_string = String::new();
    if params.do_isomeric_smiles
        && atom.chiral_tag() != ChiralType::Unspecified
        && !atom.has_prop(cp::_brokenChirality)
    {
        at_string = atom_chirality_info(mol, idx)?;
    }
    let mut needs_bracket = true;
    if !has_custom && !params.all_hs_explicit {
        needs_bracket = atom_needs_bracket(mol, idx, &at_string, params)?;
    }
    let mut res = String::new();
    if needs_bracket {
        res.push('[');
    }
    if isotope != 0 && params.do_isomeric_smiles {
        res += &isotope.to_string();
    }
    if !params.do_kekule
        && atom.is_aromatic()
        && symb.as_bytes().first().is_some_and(u8::is_ascii_uppercase)
        && matches!(num, 5 | 6 | 7 | 8 | 14 | 15 | 16 | 33 | 34 | 52)
    {
        let mut b = symb.into_bytes();
        b[0] = b[0].to_ascii_lowercase();
        symb = String::from_utf8(b).unwrap_or_default();
    }
    res += &symb;
    res += &at_string;
    if needs_bracket {
        let tot_hs = mol.total_num_hs(idx, false)?;
        if tot_hs > 0 {
            res.push('H');
            if tot_hs > 1 {
                res += &tot_hs.to_string();
            }
        }
        if fc > 0 {
            res.push('+');
            if fc > 1 {
                res += &fc.to_string();
            }
        } else if fc < 0 {
            if fc < -1 {
                res += &fc.to_string();
            } else {
                res.push('-');
            }
        }
        if let Some(map) = atom.atom_map_num_prop() {
            res.push(':');
            res += &map.to_string();
        }
        res.push(']');
    }
    if let Some(label) = string_prop(atom, cp::_supplementalSmilesLabel) {
        res += &label;
    }
    Ok(res)
}

/// `SmilesWrite::GetBondSmiles(bond, params, atomToLeftIdx)`; clears the
/// bond's `_TraversalRingClosureBond` property like upstream.
pub fn get_bond_smiles(
    mol: &mut ROMol,
    bond: u32,
    params: &SmilesWriteParams,
    atom_to_left_idx: i32,
) -> Result<String> {
    let b = mol.bond(bond);
    let left = if atom_to_left_idx < 0 {
        b.begin_atom_idx() as i32
    } else {
        atom_to_left_idx
    };
    let bt = b.bond_type();
    let mut aromatic = false;
    if !params.do_kekule && matches!(bt, BondType::Single | BondType::Double | BondType::Aromatic) {
        let a1 = mol.atom(left as u32);
        let a2 = mol.atom(b.other_atom_idx(left as u32)?);
        if a1.is_aromatic() && a2.is_aromatic() && (a1.atomic_num() != 0 || a2.atomic_num() != 0) {
            aromatic = true;
        }
    }
    let dir = b.bond_dir();
    let b_arom = b.is_aromatic();
    let b_begin = b.begin_atom_idx();
    mol.bond_mut(bond).clear_prop(cp::_TraversalRingClosureBond);
    let dir_symbol = |explicit_default: &'static str, default_ok: bool| -> &'static str {
        match dir {
            BondDir::EndDownRight if params.all_bonds_explicit || params.do_isomeric_smiles => "\\",
            BondDir::EndUpRight if params.all_bonds_explicit || params.do_isomeric_smiles => "/",
            BondDir::EndDownRight | BondDir::EndUpRight => "",
            _ if default_ok => explicit_default,
            _ => "",
        }
    };
    let res = match bt {
        BondType::Single => {
            if dir != BondDir::None && dir != BondDir::Unknown {
                dir_symbol("-", params.all_bonds_explicit)
            } else if params.all_bonds_explicit || (aromatic && !b_arom) {
                "-"
            } else {
                ""
            }
        }
        BondType::Double => {
            if !aromatic || !b_arom || params.all_bonds_explicit {
                "="
            } else {
                ""
            }
        }
        BondType::Triple => "#",
        BondType::Quadruple => "$",
        BondType::Aromatic => {
            if dir != BondDir::None && dir != BondDir::Unknown {
                dir_symbol(":", params.all_bonds_explicit || !aromatic)
            } else if params.all_bonds_explicit || !aromatic {
                ":"
            } else {
                ""
            }
        }
        BondType::Dative => {
            if left >= 0 && b_begin == left as u32 {
                "->"
            } else {
                "<-"
            }
        }
        _ => "~",
    };
    Ok(res.to_string())
}

/// Atoms/bonds in play and custom symbols for `FragmentSmilesConstruct`.
#[derive(Default, Clone, Copy)]
struct FragmentSpec<'a> {
    atoms_in_play: Option<&'a [bool]>,
    bonds_in_play: Option<&'a [bool]>,
    atom_symbols: Option<&'a [String]>,
    bond_symbols: Option<&'a [String]>,
}

/// `SmilesWrite::FragmentSmilesConstruct`
#[allow(clippy::too_many_arguments)]
fn fragment_smiles_construct(
    mol: &mut ROMol,
    atom_idx: u32,
    colors: &mut [AtomColor],
    ranks: &[u32],
    params: &SmilesWriteParams,
    atom_ordering: &mut Vec<u32>,
    bond_ordering: &mut Vec<u32>,
    spec: FragmentSpec<'_>,
) -> Result<String> {
    if params.do_kekule {
        match (spec.atoms_in_play, spec.bonds_in_play) {
            (Some(a), Some(b)) => kekulize_fragment(mol, a, b, true, true, DEFAULT_MAX_BACKTRACKS)?,
            _ => kekulize(mol, true, true, DEFAULT_MAX_BACKTRACKS)?,
        }
    }
    let mut mol_stack = Vec::with_capacity((mol.num_atoms() + mol.num_bonds()) as usize);
    let mut ring_closure_map: BTreeMap<i32, i32> = BTreeMap::new();
    if !params.canonical {
        mol.set_prop(cp::_StereochemDone, 1i32, false)?;
    }
    let mut ring_closures_to_erase: Vec<i32> = Vec::new();
    if params.canonical && params.do_isomeric_smiles {
        canonicalize_enhanced_stereo(mol, Some(ranks))?;
    }
    canonicalize_fragment(
        mol,
        atom_idx,
        colors,
        ranks,
        &mut mol_stack,
        FragmentOptions {
            atoms_in_play: spec.atoms_in_play,
            bonds_in_play: spec.bonds_in_play,
            bond_symbols: spec.bond_symbols,
            do_isomeric_smiles: params.do_isomeric_smiles,
            do_random: params.do_random,
            do_chiral_inversions: true,
        },
    )?;
    let mut res = String::new();
    for elem in mol_stack {
        match elem {
            MolStackElem::Atom(a) => {
                for r in ring_closures_to_erase.drain(..) {
                    ring_closure_map.remove(&r);
                }
                match spec.atom_symbols {
                    None => res += &get_atom_smiles(mol, a, params)?,
                    Some(s) => res += &s[a as usize],
                }
                atom_ordering.push(a);
            }
            MolStackElem::Bond(b, left) => {
                match spec.bond_symbols {
                    None => res += &get_bond_smiles(mol, b, params, left)?,
                    Some(s) => res += &s[b as usize],
                }
                bond_ordering.push(b);
            }
            MolStackElem::Ring(ring_idx) => {
                let closure_val = if let Some(&v) = ring_closure_map.get(&ring_idx) {
                    ring_closures_to_erase.push(ring_idx);
                    v
                } else {
                    let mut v = 1;
                    while ring_closure_map.values().any(|&x| x == v) {
                        v += 1;
                    }
                    ring_closure_map.insert(ring_idx, v);
                    v
                };
                if closure_val < 10 {
                    res.push((b'0' + closure_val as u8) as char);
                } else if closure_val < 100 {
                    res.push('%');
                    res += &closure_val.to_string();
                } else {
                    res += "%(";
                    res += &closure_val.to_string();
                    res.push(')');
                }
            }
            MolStackElem::BranchOpen(_) => res.push('('),
            MolStackElem::BranchClose(_) => res.push(')'),
        }
    }
    Ok(res)
}

/// `MolToSmiles(mol, params)`: sets `_smilesAtomOutputOrder` /
/// `_smilesBondOutputOrder` on `mol`.
pub fn mol_to_smiles(mol: &mut ROMol, params: &SmilesWriteParams) -> Result<String> {
    mol_to_smiles_detail(mol, params, false, true)
}

/// `SmilesWrite::detail::MolToSmiles`. When writing CXSMILES
/// (`doing_cx_smiles`) the features that live in the extension (stereo
/// groups if `include_stereo_groups`, wiggly and unspecified double bonds,
/// dative and hydrogen bonds) are kept out of the SMILES proper.
pub(crate) fn mol_to_smiles_detail(
    mol: &mut ROMol,
    params: &SmilesWriteParams,
    doing_cx_smiles: bool,
    include_stereo_groups: bool,
) -> Result<String> {
    if mol.num_atoms() == 0 {
        return Ok(String::new());
    }
    if params.rooted_at_atom >= 0 && params.rooted_at_atom as u32 >= mol.num_atoms() {
        return Err(RDKitError::precondition(
            "rootedAtAtom must be less than the number of atoms",
        ));
    }
    let (mut mols, frags_atom_mapping) = get_mol_frags_as_mols(mol, false, false)?;
    let mut frags_rooted = Vec::new();
    let mut frags_bond_mapping: Vec<Vec<u32>> = Vec::new();
    let mut ats_present = vec![false; mol.num_atoms() as usize];
    for ats in &frags_atom_mapping {
        ats_present.iter_mut().for_each(|x| *x = false);
        for &a in ats {
            ats_present[a as usize] = true;
        }
        let mut rooted = -1;
        if params.rooted_at_atom >= 0 && ats_present[params.rooted_at_atom as usize] {
            let first = ats_present.iter().position(|&x| x).unwrap_or(0) as i32;
            rooted = params.rooted_at_atom - first;
        }
        frags_rooted.push(rooted);
        frags_bond_mapping.push(
            mol.bonds()
                .iter()
                .filter(|b| {
                    ats_present[b.begin_atom_idx() as usize]
                        && ats_present[b.end_atom_idx() as usize]
                })
                .map(|b| b.idx())
                .collect(),
        );
    }

    let mut vfragsmi = Vec::with_capacity(mols.len());
    let mut all_atom_ordering = Vec::new();
    let mut all_bond_ordering = Vec::new();
    for (frag_idx, tmol) in mols.iter_mut().enumerate() {
        let mut atom_map_nums = vec![0i32; tmol.num_atoms() as usize];
        for i in 0..tmol.num_atoms() {
            if params.canonical && params.ignore_atom_map_numbers {
                atom_map_nums[i as usize] = tmol.atom(i).atom_map_num();
                tmol.atom_mut(i).set_atom_map_num(0, true)?;
            }
            tmol.update_atom_property_cache(i, false)?;
        }
        if params.do_isomeric_smiles {
            tmol.set_prop(cp::_doIsoSmiles, 1i32, false)?;
            if !tmol.has_prop(cp::_StereochemDone) {
                assign_stereochemistry(tmol, params.clean_stereo, false, false)?;
            }
        }
        if !doing_cx_smiles || !include_stereo_groups {
            // otherwise they would be used in the canonicalization
            tmol.set_stereo_groups(Vec::new())?;
        }
        if !doing_cx_smiles {
            for bi in 0..tmol.num_bonds() {
                let b = tmol.bond_mut(bi);
                if matches!(b.bond_dir(), BondDir::Unknown | BondDir::EitherDouble) {
                    b.set_bond_dir(BondDir::None);
                }
                if b.stereo() == BondStereo::StereoAny {
                    b.set_stereo(BondStereo::StereoNone)?;
                }
            }
        }
        if doing_cx_smiles || !params.include_dative_bonds {
            // the extension has the coordinate bonds
            for bi in 0..tmol.num_bonds() {
                if tmol.bond(bi).bond_type() == BondType::Dative {
                    tmol.bond_mut(bi).set_bond_type(BondType::Single);
                    let begin = tmol.bond(bi).begin_atom_idx();
                    tmol.calc_explicit_valence(begin, false)?;
                }
            }
        }
        if doing_cx_smiles {
            // hydrogen bonds are single bonds here and listed in the extension
            for bi in 0..tmol.num_bonds() {
                if tmol.bond(bi).bond_type() == BondType::Hydrogen {
                    tmol.bond_mut(bi).set_bond_type(BondType::Single);
                }
            }
        }
        let mut rooted = frags_rooted[frag_idx];
        if params.do_random && rooted == -1 {
            let n = tmol.num_atoms();
            rooted = (with_global_generator(0, |g| g.next()) % n) as i32;
        }
        let n_atoms = tmol.num_atoms() as usize;
        let ranks: Vec<u32> = if params.canonical {
            let r = rank_mol_atoms(
                tmol,
                RankOptions {
                    break_ties: true,
                    include_chirality: params.do_isomeric_smiles,
                    include_isotopes: params.do_isomeric_smiles,
                    include_atom_maps: true,
                    include_chiral_presence: false,
                    include_stereo_groups: params.do_isomeric_smiles,
                    use_non_stereo_ranks: false,
                    ..RankOptions::default()
                },
            )?;
            if params.ignore_atom_map_numbers {
                for i in 0..tmol.num_atoms() {
                    tmol.atom_mut(i)
                        .set_atom_map_num(atom_map_nums[i as usize], true)?;
                }
            }
            r
        } else {
            (0..n_atoms as u32).collect()
        };
        let mut colors = vec![AtomColor::White; n_atoms];
        let next_atom_idx = if rooted >= 0 {
            rooted as u32
        } else {
            let mut next_rank = n_atoms as u32 + 1;
            let mut next = None;
            for i in 0..n_atoms {
                if colors[i] == AtomColor::White && ranks[i] < next_rank {
                    next_rank = ranks[i];
                    next = Some(i as u32);
                }
            }
            next.ok_or_else(|| RDKitError::invariant("no start atom found"))?
        };
        let mut atom_ordering = Vec::new();
        let mut bond_ordering = Vec::new();
        let smi = fragment_smiles_construct(
            tmol,
            next_atom_idx,
            &mut colors,
            &ranks,
            params,
            &mut atom_ordering,
            &mut bond_ordering,
            FragmentSpec::default(),
        )?;
        vfragsmi.push(smi);
        all_atom_ordering.push(
            atom_ordering
                .iter()
                .map(|&v| frags_atom_mapping[frag_idx][v as usize] as u32)
                .collect::<Vec<u32>>(),
        );
        all_bond_ordering.push(
            bond_ordering
                .iter()
                .map(|&v| frags_bond_mapping[frag_idx][v as usize])
                .collect::<Vec<u32>>(),
        );
    }

    let mut result = String::new();
    let mut flat_atoms: Vec<u32> = Vec::with_capacity(mol.num_atoms() as usize);
    let mut flat_bonds: Vec<u32> = Vec::with_capacity(mol.num_bonds() as usize);
    let mut tmp: Vec<(String, Vec<u32>, Vec<u32>)> = vfragsmi
        .into_iter()
        .zip(all_atom_ordering)
        .zip(all_bond_ordering)
        .map(|((s, a), b)| (s, a, b))
        .collect();
    if params.canonical {
        // std::tuple comparison: string (bytewise), then the vectors
        tmp.sort();
    }
    let n = tmp.len();
    for (i, (s, a, b)) in tmp.into_iter().enumerate() {
        result += &s;
        if i + 1 < n {
            result.push('.');
        }
        flat_atoms.extend(a);
        flat_bonds.extend(b);
    }
    mol.set_prop(cp::_smilesAtomOutputOrder, flat_atoms, true)?;
    mol.set_prop(cp::_smilesBondOutputOrder, flat_bonds, true)?;
    Ok(result)
}

/// `MolToRandomSmilesVect`
pub fn mol_to_random_smiles_vect(
    mol: &mut ROMol,
    num_smiles: u32,
    random_seed: u32,
    params: &SmilesWriteParams,
) -> Result<Vec<String>> {
    if random_seed > 0 {
        with_global_generator(random_seed as i32, |_| ());
    }
    let p = SmilesWriteParams {
        canonical: false,
        rooted_at_atom: -1,
        do_random: true,
        ..params.clone()
    };
    (0..num_smiles).map(|_| mol_to_smiles(mol, &p)).collect()
}

/// `MolFragmentToSmiles(mol, params, atomsToUse, bondsToUse, atomSymbols, bondSymbols)`
pub fn mol_fragment_to_smiles(
    mol: &mut ROMol,
    params: &SmilesWriteParams,
    atoms_to_use: &[i32],
    bonds_to_use: Option<&[i32]>,
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
) -> Result<String> {
    if atoms_to_use.is_empty() {
        return Err(RDKitError::precondition("no atoms provided"));
    }
    let na = mol.num_atoms() as usize;
    let nb = mol.num_bonds() as usize;
    if params.rooted_at_atom >= 0
        && (params.rooted_at_atom as usize >= na || !atoms_to_use.contains(&params.rooted_at_atom))
    {
        return Err(RDKitError::precondition(
            "rootedAtAtom not found in atomsToUse",
        ));
    }
    if atom_symbols.is_some_and(|s| s.len() < na) {
        return Err(RDKitError::precondition("bad atomSymbols vector"));
    }
    if bond_symbols.is_some_and(|s| s.len() < nb) {
        return Err(RDKitError::precondition("bad bondSymbols vector"));
    }
    if na == 0 {
        return Ok(String::new());
    }
    for &a in atoms_to_use {
        if a < 0 || a as usize >= na {
            return Err(RDKitError::range("idx", a, na));
        }
    }
    let mut rooted = params.rooted_at_atom;
    let mut tmol = ROMol::copy_from(mol, true, -1);
    if params.do_isomeric_smiles {
        tmol.set_prop(cp::_doIsoSmiles, 1i32, false)?;
    }
    let mut atoms_in_play = vec![false; na];
    for &a in atoms_to_use {
        atoms_in_play[a as usize] = true;
    }
    let mut bonds_in_play = vec![false; nb];
    if let Some(bonds) = bonds_to_use {
        for &b in bonds {
            if b < 0 || b as usize >= nb {
                return Err(RDKitError::range("idx", b, nb));
            }
            bonds_in_play[b as usize] = true;
        }
    } else {
        if params.rooted_at_atom >= 0 && get_mol_frags_atoms(mol).len() != 1 {
            return Err(RDKitError::precondition(
                "rootedAtAtom can only be used with molecules that have a single fragment",
            ));
        }
        for &a in atoms_to_use {
            for b in mol.atom_bonds(a as u32) {
                if atoms_in_play[b.other_atom_idx(a as u32)? as usize] {
                    bonds_in_play[b.idx() as usize] = true;
                }
            }
        }
    }
    if mol.ring_info().is_initialized() {
        let ri = tmol.ring_info_mut();
        ri.reset(true);
        ri.initialize(crate::graphmol::ring_info::FindRingType::OtherOrUnknown);
        let src = mol.ring_info();
        for (aring, bring) in src.atom_rings().iter().zip(src.bond_rings()) {
            if aring.iter().all(|&a| atoms_in_play[a as usize])
                && bring.iter().all(|&b| bonds_in_play[b as usize])
            {
                tmol.ring_info_mut()
                    .add_ring(aring.clone(), bring.clone())?;
            }
        }
    }
    if tmol.needs_update_property_cache() {
        for i in 0..tmol.num_atoms() {
            tmol.update_atom_property_cache(i, false)?;
        }
    }

    if params.do_isomeric_smiles {
        if !mol.has_prop(cp::_StereochemDone) {
            assign_stereochemistry(&mut tmol, true, false, false)?;
        } else {
            tmol.set_prop(cp::_StereochemDone, 1i32, false)?;
            for &a in atoms_to_use {
                if let Some(code) = string_prop(mol.atom(a as u32), cp::_CIPCode) {
                    tmol.atom_mut(a as u32)
                        .set_prop(cp::_CIPCode, code, false)?;
                }
            }
        }
        for bi in 0..tmol.num_bonds() {
            let b = tmol.bond(bi);
            if !bonds_in_play[bi as usize]
                || b.bond_type() != BondType::Double
                || b.stereo() == BondStereo::StereoNone
            {
                continue;
            }
            let sa = b.stereo_atoms().to_vec();
            if sa.len() != 2 {
                continue;
            }
            for (k, bond_atom) in [(0usize, b.begin_atom_idx()), (1, b.end_atom_idx())] {
                let stereo_atom = sa[k];
                if atoms_in_play[stereo_atom as usize] {
                    continue;
                }
                if tmol.atom_degree(bond_atom) > 2 {
                    let cur = tmol.bond(bi).stereo_atoms().to_vec();
                    let replacement = tmol
                        .atom_neighbors(bond_atom)
                        .find(|&n| n != stereo_atom as u32 && atoms_in_play[n as usize]);
                    if let Some(n) = replacement {
                        let new_atoms = if stereo_atom == cur[0] {
                            vec![n as i32, cur[1]]
                        } else {
                            vec![cur[0], n as i32]
                        };
                        let bm = tmol.bond_mut(bi);
                        *bm.stereo_atoms_mut() = new_atoms;
                        match bm.stereo() {
                            BondStereo::StereoZ | BondStereo::StereoCis => {
                                bm.set_stereo(BondStereo::StereoTrans)?
                            }
                            BondStereo::StereoE | BondStereo::StereoTrans => {
                                bm.set_stereo(BondStereo::StereoCis)?
                            }
                            _ => {}
                        }
                    } else {
                        tmol.bond_mut(bi).set_stereo(BondStereo::StereoNone)?;
                        break;
                    }
                } else {
                    tmol.bond_mut(bi).set_stereo(BondStereo::StereoNone)?;
                    break;
                }
            }
        }
    }
    let ranks: Vec<u32> = if params.canonical {
        rank_fragment_atoms(
            &mut tmol,
            &atoms_in_play,
            &bonds_in_play,
            &FragmentRankOptions {
                atom_symbols,
                bond_symbols,
                break_ties: true,
                include_chirality: params.do_isomeric_smiles,
                include_isotopes: params.do_isomeric_smiles,
                include_atom_maps: !params.ignore_atom_map_numbers,
                include_chiral_presence: false,
                include_ring_stereo: true,
                compare_codes: None,
            },
        )?
    } else {
        (0..tmol.num_atoms()).collect()
    };

    let mut colors = vec![AtomColor::Black; tmol.num_atoms() as usize];
    for &a in atoms_to_use {
        colors[a as usize] = AtomColor::White;
    }
    let mut res = String::new();
    let mut atom_ordering = Vec::new();
    let mut bond_ordering = Vec::new();
    loop {
        let next_atom_idx = if rooted >= 0 {
            let r = rooted as u32;
            rooted = -1;
            r
        } else {
            let mut next_rank = tmol.num_atoms() + 1;
            let mut next = None;
            for &i in atoms_to_use {
                if colors[i as usize] == AtomColor::White && ranks[i as usize] < next_rank {
                    next_rank = ranks[i as usize];
                    next = Some(i as u32);
                }
            }
            next.ok_or_else(|| RDKitError::invariant("no start atom found"))?
        };
        let sub = fragment_smiles_construct(
            &mut tmol,
            next_atom_idx,
            &mut colors,
            &ranks,
            params,
            &mut atom_ordering,
            &mut bond_ordering,
            FragmentSpec {
                atoms_in_play: Some(&atoms_in_play),
                bonds_in_play: Some(&bonds_in_play),
                atom_symbols,
                bond_symbols,
            },
        )?;
        res += &sub;
        if colors.contains(&AtomColor::White) {
            res.push('.');
        } else {
            break;
        }
    }
    mol.set_prop(cp::_smilesAtomOutputOrder, atom_ordering, true)?;
    mol.set_prop(cp::_smilesBondOutputOrder, bond_ordering, true)?;
    Ok(res)
}
