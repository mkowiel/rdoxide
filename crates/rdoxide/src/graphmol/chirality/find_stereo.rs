//! New stereo perception: `StereoInfo` and `Chirality::findPotentialStereo`.
//!
//! Ported from RDKit `Code/GraphMol/FindStereo.cpp` (2026.09.1).
//!
//! Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! The algorithm repeatedly ranks the atoms with a "symmetry class" refinement
//! in which every specified or possible stereo element carries a distinguishing
//! label. A stereo atom or bond is real when its controlling neighbours end up
//! with different ranks. Unlike upstream the result is not also stored in the
//! molecule's private `_potentialStereo` property (nothing reads it).

use super::legacy::{
    find_chiral_atom_special_cases, get_atom_nonzero_degree, is_atom_potential_tetrahedral_center,
};
use super::{MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO, get_allow_nontetrahedral_chirality};
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{Bond, BondDir, BondStereo, BondType};
use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, symmetrize_sssr};
use crate::graphmol::new_canon::{AtomCompareCode, FragmentRankOptions, rank_fragment_atoms};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::utils::count_swaps_to_interconvert;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// `Atom::NOATOM`: marks a missing (implicit H) controlling atom.
pub const NOATOM: u32 = u32::MAX;

/// `Chirality::StereoType`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(u8)]
pub enum StereoType {
    #[default]
    Unspecified = 0,
    AtomTetrahedral = 1,
    AtomSquarePlanar = 2,
    AtomTrigonalBipyramidal = 3,
    AtomOctahedral = 4,
    /// a single double bond and odd-numbered cumulenes
    BondDouble = 5,
    /// even-numbered cumulenes
    BondCumuleneEven = 6,
    BondAtropisomer = 7,
}

/// `Chirality::StereoDescriptor`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(u8)]
pub enum StereoDescriptor {
    #[default]
    None = 0,
    TetCw = 1,
    TetCcw = 2,
    BondCis = 3,
    BondTrans = 4,
    BondAtropCw = 5,
    BondAtropCcw = 6,
}

/// `Chirality::StereoSpecified`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(u8)]
pub enum StereoSpecified {
    /// no information provided
    #[default]
    Unspecified = 0,
    Specified = 1,
    /// deliberately marked as unknown
    Unknown = 2,
}

/// `Chirality::StereoInfo`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoInfo {
    pub stereo_type: StereoType,
    pub specified: StereoSpecified,
    /// the atom or bond index the stereo element is centered on
    pub centered_on: u32,
    pub descriptor: StereoDescriptor,
    /// for the non-tetrahedral cases
    pub permutation: u32,
    /// all atoms around the atom or bond; the order is important
    pub controlling_atoms: Vec<u32>,
}

impl Default for StereoInfo {
    fn default() -> Self {
        StereoInfo {
            stereo_type: StereoType::Unspecified,
            specified: StereoSpecified::Unspecified,
            centered_on: NOATOM,
            descriptor: StereoDescriptor::None,
            permutation: 0,
            controlling_atoms: Vec::new(),
        }
    }
}

fn under_construction<T>(what: &str) -> Result<T> {
    Err(RDKitError::Runtime(format!("Under construction: {what}")))
}

/// `translateEZLabelToCisTrans`
fn translate_ez_label_to_cis_trans(label: BondStereo) -> BondStereo {
    match label {
        BondStereo::StereoE => BondStereo::StereoTrans,
        BondStereo::StereoZ => BondStereo::StereoCis,
        other => other,
    }
}

fn int_flag(props: &impl RDProps, key: &str) -> bool {
    props
        .get_prop(key)
        .and_then(RDValue::as_int)
        .is_some_and(|v| v != 0)
}

// ------------------------------------------------------------ candidates

/// `detail::isAtomPotentialNontetrahedralCenter`
pub fn is_atom_potential_nontetrahedral_center(mol: &ROMol, idx: u32) -> Result<bool> {
    let atom = mol.atom(idx);
    let tnz_degree = get_atom_nonzero_degree(mol, idx) + mol.total_num_hs(idx, false)?;
    let anum = atom.atomic_num();
    if !(2..=6).contains(&tnz_degree) || (anum < 12 && anum != 4) {
        return Ok(false);
    }
    let tag = atom.chiral_tag();
    if (ChiralType::SquarePlanar..=ChiralType::Octahedral).contains(&tag) {
        return Ok(true);
    }
    // with at least four neighbors but nothing specified we can start to
    // imagine that it might be enhanced stereo
    Ok(tag == ChiralType::Unspecified && tnz_degree >= 4)
}

/// `detail::isAtomPotentialStereoAtom`
pub fn is_atom_potential_stereo_atom(
    mol: &ROMol,
    idx: u32,
    allow_nontetrahedral_stereo: bool,
) -> Result<bool> {
    Ok(is_atom_potential_tetrahedral_center(mol, idx)?
        || (allow_nontetrahedral_stereo && is_atom_potential_nontetrahedral_center(mol, idx)?))
}

/// `detail::isBondPotentialStereoBond`
pub fn is_bond_potential_stereo_bond(mol: &ROMol, idx: u32) -> Result<bool> {
    let bond = mol.bond(idx);
    if bond.bond_type() != BondType::Double {
        return Ok(false);
    }
    // each end needs at least 2 and at most 3 neighbours (implicit Hs
    // included); a ring bond must be in no ring smaller than 8
    let (beg, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
    let (beg_degree, end_degree) = (mol.total_degree(beg)?, mol.total_degree(end)?);
    if beg_degree > 1
        && beg_degree < 4
        && end_degree > 1
        && end_degree < 4
        && mol.total_num_hs(beg, true)? < 2
        && mol.total_num_hs(end, true)? < 2
    {
        let in_small_ring = mol.ring_info().bond_rings().iter().any(|ring| {
            ring.len() < MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO as usize
                && ring.contains(&(idx as i32))
        });
        return Ok(!in_small_ring);
    }
    Ok(false)
}

// ------------------------------------------------------------ stereo info

/// `detail::getStereoInfo(const Bond *)`
pub fn get_stereo_info_bond(mol: &ROMol, idx: u32) -> Result<StereoInfo> {
    let bond = mol.bond(idx);
    let (beg, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
    let (beg_degree, end_degree) = (mol.atom_degree(beg), mol.atom_degree(end));
    let mut sinfo = StereoInfo::default();
    if bond.bond_type() == BondType::Double {
        if beg_degree < 1 || end_degree < 1 || beg_degree > 3 || end_degree > 3 {
            return Err(RDKitError::Value(
                "invalid atom degree in getStereoInfo(bond)".into(),
            ));
        }
        sinfo.stereo_type = StereoType::BondDouble;
        sinfo.centered_on = idx;
        let mut seen_squiggle = false;
        for atom in [beg, end] {
            for nbr in mol.atom_bonds(atom) {
                if nbr.idx() != idx {
                    seen_squiggle |= nbr.bond_dir() == BondDir::Unknown;
                    sinfo.controlling_atoms.push(nbr.other_atom_idx(atom)?);
                }
            }
            for _ in mol.atom_degree(atom)..3 {
                sinfo.controlling_atoms.push(NOATOM);
            }
        }
        if !seen_squiggle {
            // a squiggle bond to an H sets _UnknownStereo on the atom
            seen_squiggle = int_flag(mol.atom(beg), cp::_UnknownStereo)
                || int_flag(mol.atom(end), cp::_UnknownStereo);
        }
        let mut stereo = bond.stereo();
        if stereo == BondStereo::StereoAny
            || bond.bond_dir() == BondDir::EitherDouble
            || seen_squiggle
        {
            sinfo.specified = StereoSpecified::Unknown;
        } else if stereo != BondStereo::StereoNone {
            if matches!(stereo, BondStereo::StereoE | BondStereo::StereoZ) {
                stereo = translate_ez_label_to_cis_trans(stereo);
            }
            sinfo.specified = StereoSpecified::Specified;
            let satoms = bond.stereo_atoms();
            if satoms.len() != 2 {
                return Err(RDKitError::Value(
                    "only can support 2 stereo neighbors".into(),
                ));
            }
            let ca = &sinfo.controlling_atoms;
            let first_at_begin = if satoms[0] == ca[0] as i32 {
                true
            } else if satoms[0] == ca[1] as i32 {
                false
            } else {
                return Err(RDKitError::Value(
                    "controlling atom mismatch at begin".into(),
                ));
            };
            let first_at_end = if satoms[1] == ca[2] as i32 {
                true
            } else if satoms[1] == ca[3] as i32 {
                false
            } else {
                return Err(RDKitError::Value("controlling atom mismatch at end".into()));
            };
            if first_at_begin ^ first_at_end {
                stereo = if stereo == BondStereo::StereoCis {
                    BondStereo::StereoTrans
                } else {
                    BondStereo::StereoCis
                };
            }
            sinfo.descriptor = match stereo {
                BondStereo::StereoCis => StereoDescriptor::BondCis,
                BondStereo::StereoTrans => StereoDescriptor::BondTrans,
                _ => return under_construction("unrecognized bond stereo type"),
            };
        }
    } else if bond.bond_type() == BondType::Single
        && matches!(
            bond.stereo(),
            BondStereo::StereoAtropCCW | BondStereo::StereoAtropCW
        )
    {
        if beg_degree < 2 || end_degree < 2 || beg_degree > 3 || end_degree > 3 {
            return Err(RDKitError::Value(
                "invalid atom degree in getStereoInfo(bond)".into(),
            ));
        }
        sinfo.stereo_type = StereoType::BondAtropisomer;
        sinfo.centered_on = idx;
        for (atom, degree) in [(beg, beg_degree), (end, end_degree)] {
            for nbr in mol.atom_bonds(atom) {
                if nbr.idx() != idx {
                    sinfo.controlling_atoms.push(nbr.other_atom_idx(atom)?);
                }
            }
            if degree == 2 {
                sinfo.controlling_atoms.push(NOATOM);
            }
        }
        sinfo.specified = StereoSpecified::Specified;
        sinfo.descriptor = match bond.stereo() {
            BondStereo::StereoAtropCW => StereoDescriptor::BondAtropCw,
            _ => StereoDescriptor::BondAtropCcw,
        };
    }
    Ok(sinfo)
}

/// `detail::getStereoInfo(const Atom *)`
pub fn get_stereo_info_atom(mol: &ROMol, idx: u32) -> Result<StereoInfo> {
    let atom = mol.atom(idx);
    let mut sinfo = StereoInfo {
        stereo_type: StereoType::AtomTetrahedral,
        centered_on: idx,
        ..StereoInfo::default()
    };
    let mut explicit_unknown = false;
    for bnd in mol.atom_bonds(idx) {
        if bnd.bond_dir() == BondDir::Unknown {
            explicit_unknown = true;
        } else if !explicit_unknown {
            explicit_unknown = int_flag(bnd, cp::_UnknownStereo);
        }
        sinfo.controlling_atoms.push(bnd.other_atom_idx(idx)?);
    }
    let orig_nbr_order = sinfo.controlling_atoms.clone();
    sinfo.controlling_atoms.sort_unstable();

    if explicit_unknown {
        sinfo.specified = StereoSpecified::Unknown;
        return Ok(sinfo);
    }
    let mut stereo = atom.chiral_tag();
    if matches!(
        stereo,
        ChiralType::TetrahedralCCW | ChiralType::TetrahedralCW
    ) {
        sinfo.specified = StereoSpecified::Specified;
        let n_swaps = count_swaps_to_interconvert(&orig_nbr_order, &sinfo.controlling_atoms)?;
        if n_swaps % 2 == 1 {
            stereo = if stereo == ChiralType::TetrahedralCCW {
                ChiralType::TetrahedralCW
            } else {
                ChiralType::TetrahedralCCW
            };
        }
        sinfo.descriptor = if stereo == ChiralType::TetrahedralCCW {
            StereoDescriptor::TetCcw
        } else {
            StereoDescriptor::TetCw
        };
    } else if get_allow_nontetrahedral_chirality()
        && is_atom_potential_nontetrahedral_center(mol, idx)?
    {
        if stereo == ChiralType::Unspecified {
            match mol.total_degree(idx)? {
                // don't assume non-tetrahedral chirality
                4 => stereo = ChiralType::Tetrahedral,
                5 => stereo = ChiralType::TrigonalBipyramidal,
                6 => stereo = ChiralType::Octahedral,
                _ => {}
            }
        }
        sinfo.descriptor = StereoDescriptor::None;
        match stereo {
            ChiralType::Tetrahedral => sinfo.stereo_type = StereoType::AtomTetrahedral,
            ChiralType::SquarePlanar => sinfo.stereo_type = StereoType::AtomSquarePlanar,
            ChiralType::TrigonalBipyramidal => {
                sinfo.stereo_type = StereoType::AtomTrigonalBipyramidal
            }
            ChiralType::Octahedral => sinfo.stereo_type = StereoType::AtomOctahedral,
            _ => {}
        }
        if let Some(permutation) = atom
            .props
            .get(cp::_chiralPermutation)
            .and_then(RDValue::as_uint)
        {
            sinfo.permutation = permutation;
            // a permutation of zero explicitly says the chirality is unknown
            sinfo.specified = if permutation == 0 {
                StereoSpecified::Unknown
            } else {
                StereoSpecified::Specified
            };
        }
    }
    Ok(sinfo)
}

/// `getBondSymbol`
pub fn get_bond_symbol(bond: &Bond) -> &'static str {
    if bond.is_aromatic() {
        return ":";
    }
    match bond.bond_type() {
        BondType::Single => "-",
        BondType::Double => "=",
        BondType::Triple => "#",
        BondType::Aromatic => ":",
        _ => "?",
    }
}

// ------------------------------------------------------------ the algorithm

/// The state a stereo-perception run threads through its rounds.
struct Perception {
    known_atoms: Vec<bool>,
    possible_atoms: Vec<bool>,
    fixed_atoms: Vec<bool>,
    atom_codes: Vec<AtomCompareCode>,
    known_bonds: Vec<bool>,
    possible_bonds: Vec<bool>,
    fixed_bonds: Vec<bool>,
    bond_symbols: Vec<String>,
    /// number of rings with possible ring stereo each atom is in
    ring_stereo_atoms: Vec<u32>,
    /// number of rings with possible ring stereo each bond is in
    ring_stereo_bonds: Vec<u32>,
}

/// `AtomCompareState`: the label an atom contributes to its compare code.
#[derive(Clone, Copy)]
#[repr(u32)]
enum AtomCompareState {
    None = 0,
    Unknown = 1,
    Possible = 2,
    Ccw = 3,
    Cw = 4,
    Stereo = 5,
}

fn set_atom_compare_state(code: &mut AtomCompareCode, state: AtomCompareState, discriminator: u32) {
    code.secondary = (u64::from(state as u32) << 32) | u64::from(discriminator);
}

fn make_atom_compare_code(mol: &ROMol, idx: u32, symbol_rank: u32) -> AtomCompareCode {
    let atom = mol.atom(idx);
    let atom_map_number = if atom.atomic_num() == 0 {
        atom.props
            .get(cp::molAtomMapNumber)
            .and_then(RDValue::as_int)
            .unwrap_or(0)
    } else {
        0
    };
    // pack isotope, symbol rank and biased formal charge in comparison order
    let charge = (atom.formal_charge().wrapping_add(128)) as u8;
    AtomCompareCode {
        atom_map_number,
        degree: mol.atom_degree(idx),
        primary: (u64::from(atom.isotope()) << 40)
            | (u64::from(symbol_rank) << 8)
            | u64::from(charge),
        secondary: 0,
    }
}

impl Perception {
    fn new(mol: &ROMol) -> Self {
        let (na, nb) = (mol.num_atoms() as usize, mol.num_bonds() as usize);
        Perception {
            known_atoms: vec![false; na],
            possible_atoms: vec![false; na],
            fixed_atoms: vec![false; na],
            atom_codes: vec![AtomCompareCode::default(); na],
            known_bonds: vec![false; nb],
            possible_bonds: vec![false; nb],
            fixed_bonds: vec![false; nb],
            bond_symbols: vec![String::new(); nb],
            ring_stereo_atoms: vec![0; na],
            ring_stereo_bonds: vec![0; nb],
        }
    }

    /// `initAtomInfo`
    fn init_atom_info(
        &mut self,
        mol: &mut ROMol,
        flag_possible: bool,
        clean_it: bool,
    ) -> Result<()> {
        // rank of each distinct atom symbol, ordered by the symbol text
        let symbols: Vec<String> = mol.atoms().iter().map(|a| a.symbol()).collect();
        let mut unique: Vec<&String> = symbols.iter().collect();
        unique.sort();
        unique.dedup();
        let symbol_rank = |s: &String| unique.binary_search(&s).map_or(0, |r| r as u32);

        let allow_nontetrahedral = get_allow_nontetrahedral_chirality();
        for aidx in 0..mol.num_atoms() {
            if mol.atom(aidx).needs_update_property_cache() {
                mol.update_atom_property_cache(aidx, false)?;
            }
            let i = aidx as usize;
            self.atom_codes[i] = make_atom_compare_code(mol, aidx, symbol_rank(&symbols[i]));
            if is_atom_potential_stereo_atom(mol, aidx, allow_nontetrahedral)? {
                let sinfo = get_stereo_info_atom(mol, aidx)?;
                match sinfo.specified {
                    StereoSpecified::Unknown => {
                        self.known_atoms[i] = true;
                        set_atom_compare_state(
                            &mut self.atom_codes[i],
                            AtomCompareState::Unknown,
                            aidx,
                        );
                    }
                    StereoSpecified::Specified => {
                        self.known_atoms[i] = true;
                        let state = match sinfo.descriptor {
                            StereoDescriptor::TetCcw => AtomCompareState::Ccw,
                            StereoDescriptor::TetCw => AtomCompareState::Cw,
                            _ => AtomCompareState::Stereo,
                        };
                        set_atom_compare_state(&mut self.atom_codes[i], state, 0);
                    }
                    StereoSpecified::Unspecified => {
                        if flag_possible {
                            self.possible_atoms[i] = true;
                            if !clean_it {
                                set_atom_compare_state(
                                    &mut self.atom_codes[i],
                                    AtomCompareState::Possible,
                                    aidx,
                                );
                            }
                        }
                    }
                }
            } else if clean_it {
                mol.atom_mut(aidx).set_chiral_tag(ChiralType::Unspecified);
            }
        }
        Ok(())
    }

    /// `initBondInfo`
    fn init_bond_info(
        &mut self,
        mol: &mut ROMol,
        flag_possible: bool,
        clean_it: bool,
    ) -> Result<()> {
        for bidx in 0..mol.num_bonds() {
            let i = bidx as usize;
            self.bond_symbols[i] = get_bond_symbol(mol.bond(bidx)).to_string();
            if is_bond_potential_stereo_bond(mol, bidx)? {
                let sinfo = get_stereo_info_bond(mol, bidx)?;
                match sinfo.specified {
                    StereoSpecified::Unknown => {
                        self.known_bonds[i] = true;
                        self.bond_symbols[i] += &format!("_{bidx}");
                    }
                    StereoSpecified::Specified => {
                        self.known_bonds[i] = true;
                        self.bond_symbols[i] += match sinfo.descriptor {
                            StereoDescriptor::BondCis => "_cis",
                            StereoDescriptor::BondTrans => "_trans",
                            _ => "_STEREO",
                        };
                    }
                    StereoSpecified::Unspecified => {
                        if flag_possible {
                            self.possible_bonds[i] = true;
                            if !clean_it {
                                self.bond_symbols[i] += &format!("_{bidx}");
                            }
                        }
                    }
                }
            } else {
                let current = mol.bond(bidx).stereo();
                match current {
                    BondStereo::StereoAtropCW => {
                        self.known_bonds[i] = true;
                        self.bond_symbols[i] += "_atropcw";
                    }
                    BondStereo::StereoAtropCCW => {
                        self.known_bonds[i] = true;
                        self.bond_symbols[i] += "_atropccw";
                    }
                    _ if clean_it => mol.bond_mut(bidx).set_stereo(BondStereo::StereoNone)?,
                    _ => {}
                }
            }
        }
        Ok(())
    }

    /// `flagRingStereo`: flags the atoms and bonds that can transmit ring
    /// stereo, either through opposite atoms of an even-sized ring or through
    /// a bond shared by fused rings. With `use_possible` unset only the known
    /// atoms/bonds count.
    fn flag_ring_stereo(&mut self, mol: &mut ROMol, use_possible: bool) -> Result<()> {
        let atom_rings: Vec<Vec<i32>> = mol.ring_info().atom_rings().to_vec();
        let bond_rings: Vec<Vec<i32>> = mol.ring_info().bond_rings().to_vec();
        let is_atom = |s: &Self, a: u32| {
            s.known_atoms[a as usize] || (use_possible && s.possible_atoms[a as usize])
        };
        let is_bond = |s: &Self, b: u32| {
            s.known_bonds[b as usize] || (use_possible && s.possible_bonds[b as usize])
        };

        let mut possible_in_ring = vec![false; mol.num_atoms() as usize];
        for (aring, bring) in atom_rings.iter().zip(&bond_rings) {
            let sz = aring.len();
            let half_size = sz / 2 + sz % 2;
            let mut n_here = 0u32;
            possible_in_ring.iter_mut().for_each(|x| *x = false);
            for ai in 0..sz {
                let aidx = aring[ai] as u32;
                if !is_atom(self, aidx) {
                    continue;
                }
                for ring_divisor in [2usize, 3] {
                    if sz % ring_divisor != 0 {
                        continue;
                    }
                    let increment = sz / ring_divisor;
                    // look at the atoms 1/2 (or 1/3, 2/3) of the way around
                    let (mut found_by_bond, mut found_by_atom) = (0usize, 0usize);
                    let mut step = increment;
                    while step < sz {
                        let other = aring[(ai + step) % sz] as u32;
                        let hit = mol
                            .atom_bond_indices(other)
                            .any(|b| is_bond(self, b) && !bring.contains(&(b as i32)));
                        if hit {
                            found_by_bond += 1;
                        }
                        if found_by_bond == 0 && is_atom(self, other) {
                            found_by_atom += 1;
                        }
                        step += increment;
                    }
                    if found_by_bond == ring_divisor - 1 || found_by_atom == ring_divisor - 1 {
                        n_here += 1 + found_by_bond as u32;
                        let mut step = 0;
                        while step < sz {
                            possible_in_ring[aring[(ai + step) % sz] as usize] = true;
                            step += increment;
                        }
                        if ring_divisor == 2 {
                            let other = aring[(ai + increment) % sz];
                            mol.atom_mut(aidx)
                                .set_prop(cp::_ringStereoOtherAtom, other, false)?;
                        }
                    }
                }
                // an atom in more than one ring: explore the common edge for
                // another potential stereo atom
                if mol.ring_info().num_atom_rings(aidx)? > 1 {
                    let mut previous = aidx;
                    for step in 1..=half_size {
                        let other = aring[(ai + step) % sz] as u32;
                        let Some(bnd) = mol.get_bond_between_atoms(previous, other)?.map(Bond::idx)
                        else {
                            break;
                        };
                        if mol.ring_info().num_bond_rings(bnd)? < 2 {
                            break; // the end of the common edge
                        }
                        if is_atom(self, other) {
                            n_here += 2;
                            possible_in_ring[aidx as usize] = true;
                            possible_in_ring[other as usize] = true;
                            mol.atom_mut(aidx).set_prop(
                                cp::_ringStereoOtherAtom,
                                other as i32,
                                false,
                            )?;
                            mol.atom_mut(other).set_prop(
                                cp::_ringStereoOtherAtom,
                                aidx as i32,
                                false,
                            )?;
                            break;
                        }
                        previous = other;
                    }
                }
            }
            // a ring with at least two stereo atoms transmits ring stereo
            if n_here > 1 {
                for &a in aring {
                    if possible_in_ring[a as usize] {
                        self.ring_stereo_atoms[a as usize] += 1;
                    }
                }
                for &b in bring {
                    self.ring_stereo_bonds[b as usize] += 1;
                }
            }
        }
        Ok(())
    }

    /// `updateAtoms`: returns true when another round is needed.
    fn update_atoms(
        &mut self,
        mol: &ROMol,
        aranks: &[u32],
        sinfos: &mut Vec<StereoInfo>,
    ) -> Result<bool> {
        let mut need_another_round = false;
        for aidx in 0..mol.num_atoms() {
            let i = aidx as usize;
            if !(self.known_atoms[i] || self.possible_atoms[i]) {
                continue;
            }
            let mut sinfo = get_stereo_info_atom(mol, aidx)?;
            if self.fixed_atoms[i] {
                sinfos.push(sinfo);
                continue;
            }
            let mut nbrs: Vec<u32> = Vec::with_capacity(sinfo.controlling_atoms.len());
            let mut have_a_dupe = false;
            if sinfo.stereo_type == StereoType::AtomTetrahedral {
                for &nbr_idx in &sinfo.controlling_atoms {
                    let rnk = aranks[nbr_idx as usize];
                    if !nbrs.contains(&rnk) {
                        nbrs.push(rnk);
                        continue;
                    }
                    // a duplicate rank: fine for a ring-stereo candidate when
                    // the bond to that neighbour is a ring-stereo bond
                    if self.ring_stereo_atoms[i] > 0 {
                        let bond = mol.get_bond_between_atoms(aidx, nbr_idx)?;
                        if bond.is_none_or(|b| self.ring_stereo_bonds[b.idx() as usize] == 0) {
                            have_a_dupe = true;
                            break;
                        }
                    } else {
                        have_a_dupe = true;
                        break;
                    }
                }
            }
            if !have_a_dupe {
                let mut acs = self.atom_codes[i];
                if !self.possible_atoms[i] {
                    let mut sorted = nbrs.clone();
                    sorted.sort_unstable();
                    // FIX (upstream): only works for tetrahedral stereo
                    if sinfo.stereo_type == StereoType::AtomTetrahedral {
                        let n_swaps = count_swaps_to_interconvert(&nbrs, &sorted)?;
                        if n_swaps % 2 == 1
                            && matches!(
                                sinfo.descriptor,
                                StereoDescriptor::TetCcw | StereoDescriptor::TetCw
                            )
                        {
                            sinfo.descriptor = if sinfo.descriptor == StereoDescriptor::TetCcw {
                                StereoDescriptor::TetCw
                            } else {
                                StereoDescriptor::TetCcw
                            };
                        }
                        match sinfo.descriptor {
                            StereoDescriptor::TetCw => {
                                set_atom_compare_state(&mut acs, AtomCompareState::Cw, 0)
                            }
                            StereoDescriptor::TetCcw => {
                                set_atom_compare_state(&mut acs, AtomCompareState::Ccw, 0)
                            }
                            _ => {}
                        }
                    }
                    self.fixed_atoms[i] = true;
                }
                if self.atom_codes[i] != acs {
                    self.atom_codes[i] = acs;
                    need_another_round = true;
                }
                sinfos.push(sinfo);
            } else {
                // only do another round if anything changes here
                need_another_round |= self.possible_atoms[i];
                self.possible_atoms[i] = false;
                set_atom_compare_state(&mut self.atom_codes[i], AtomCompareState::None, 0);
                if self.ring_stereo_atoms[i] > 0 {
                    self.ring_stereo_atoms[i] = 0;
                    need_another_round = true;
                    self.unflag_ring_stereo(mol);
                }
            }
        }
        Ok(need_another_round)
    }

    /// An atom lost its ring stereo: update the other atoms/bonds of the
    /// rings it is in (the ring-stereo bookkeeping of `updateAtoms`).
    fn unflag_ring_stereo(&mut self, mol: &ROMol) {
        let ri = mol.ring_info();
        for (aring, bring) in ri.atom_rings().iter().zip(ri.bond_rings()) {
            let mut n_here = 0;
            for &ra in aring {
                // un-fix the ring's atoms so that they are rechecked: they may
                // no longer be stereo after the current atom was declared not
                self.fixed_atoms[ra as usize] = false;
                n_here += u32::from(self.ring_stereo_atoms[ra as usize] > 0);
            }
            if n_here <= 1 {
                if n_here == 1 {
                    // the last ring-stereo atom cannot have ring stereo alone
                    if let Some(&ra) = aring
                        .iter()
                        .find(|&&ra| self.ring_stereo_atoms[ra as usize] > 0)
                    {
                        self.ring_stereo_atoms[ra as usize] -= 1;
                    }
                }
                for &rb in bring {
                    if self.ring_stereo_bonds[rb as usize] > 0 {
                        self.ring_stereo_bonds[rb as usize] -= 1;
                    }
                }
            }
        }
    }

    /// `updateBonds`: returns true when another round is needed.
    fn update_bonds(
        &mut self,
        mol: &ROMol,
        aranks: &[u32],
        sinfos: &mut Vec<StereoInfo>,
    ) -> Result<bool> {
        let mut need_another_round = false;
        for bidx in 0..mol.num_bonds() {
            let i = bidx as usize;
            if !(self.known_bonds[i] || self.possible_bonds[i]) {
                continue;
            }
            let mut sinfo = get_stereo_info_bond(mol, bidx)?;
            if sinfo.stereo_type == StereoType::Unspecified {
                continue; // neither a double bond nor an atropisomer bond
            }
            if sinfo.controlling_atoms.len() != 4 {
                return Err(RDKitError::invariant("bad controlling atoms size"));
            }
            let ca = &sinfo.controlling_atoms;
            if (ca[0] == NOATOM && ca[1] == NOATOM) || (ca[2] == NOATOM && ca[3] == NOATOM) {
                // no neighbours on one side: a single implicit H whose
                // position (cis or trans) cannot be known
                if sinfo.specified == StereoSpecified::Specified {
                    return Err(RDKitError::invariant(
                        "stereo bond without neighbors can only be unspecified",
                    ));
                }
                self.fixed_bonds[i] = true;
            }
            if self.fixed_bonds[i] {
                sinfos.push(sinfo);
                continue;
            }
            let mut have_a_dupe = false;
            let mut needs_swap = false;
            for (a, b) in [(0usize, 1usize), (2, 3)] {
                let (c1, c2) = (sinfo.controlling_atoms[a], sinfo.controlling_atoms[b]);
                if c1 == NOATOM || c2 == NOATOM {
                    continue;
                }
                if self.are_controlling_atoms_dupes(mol, bidx, c1, c2, aranks)? {
                    have_a_dupe = true;
                } else if aranks[c1 as usize] < aranks[c2 as usize] {
                    sinfo.controlling_atoms.swap(a, b);
                    needs_swap = !needs_swap;
                }
            }
            if !have_a_dupe {
                if needs_swap
                    && matches!(
                        sinfo.descriptor,
                        StereoDescriptor::BondCis | StereoDescriptor::BondTrans
                    )
                {
                    sinfo.descriptor = if sinfo.descriptor == StereoDescriptor::BondCis {
                        StereoDescriptor::BondTrans
                    } else {
                        StereoDescriptor::BondCis
                    };
                }
                let mut symbol = self.bond_symbols[i].clone();
                match sinfo.specified {
                    StereoSpecified::Specified => match sinfo.descriptor {
                        StereoDescriptor::BondCis => symbol += "_cis",
                        StereoDescriptor::BondTrans => symbol += "_trans",
                        _ => {}
                    },
                    StereoSpecified::Unknown => symbol += "_unk",
                    StereoSpecified::Unspecified => {}
                }
                if self.bond_symbols[i] != symbol {
                    self.bond_symbols[i] = symbol;
                    need_another_round = true;
                }
                if !self.possible_bonds[i] {
                    self.fixed_bonds[i] = true;
                }
                sinfos.push(sinfo);
            } else if self.possible_bonds[i] {
                self.possible_bonds[i] = false;
                self.bond_symbols[i] = get_bond_symbol(mol.bond(bidx)).to_string();
                need_another_round = true;
            }
        }
        Ok(need_another_round)
    }

    /// `areStereobondControllingAtomsDupes`: the two controlling atoms have
    /// the same rank and no tie-breaking stereo element sits opposite the bond
    /// in a common even-sized ring.
    fn are_controlling_atoms_dupes(
        &self,
        mol: &ROMol,
        bidx: u32,
        c1: u32,
        c2: u32,
        aranks: &[u32],
    ) -> Result<bool> {
        if aranks[c1 as usize] != aranks[c2 as usize] {
            return Ok(false);
        }
        let ri = mol.ring_info();
        let (m1, m2) = (ri.atom_members(c1)?, ri.atom_members(c2)?);
        let bond = mol.bond(bidx);
        let (mut it1, mut it2) = (0, 0);
        while it1 < m1.len() && it2 < m2.len() {
            if m1[it1] < m2[it2] {
                it1 += 1;
                continue;
            } else if m1[it1] > m2[it2] {
                it2 += 1;
                continue;
            }
            let ring = &ri.atom_rings()[m1[it1] as usize];
            it1 += 1;
            it2 += 1;
            if ring.len() % 2 == 1 {
                continue; // no atom directly across an odd ring
            }
            for end in [bond.begin_atom_idx(), bond.end_atom_idx()] {
                let Some(pos) = ring.iter().position(|&a| a == end as i32) else {
                    continue;
                };
                let opposite = ring[(pos + ring.len() / 2) % ring.len()] as u32;
                if self.possible_atoms[opposite as usize] || self.known_atoms[opposite as usize] {
                    return Ok(false);
                }
                // atropisomer test: the non-ring bond of the opposite atom
                if mol.atom_degree(opposite) == 3 {
                    for nbr in mol.atom_bonds(opposite) {
                        let out = nbr.other_atom_idx(opposite)?;
                        if !ring.contains(&(out as i32))
                            && (self.possible_bonds[nbr.idx() as usize]
                                || self.known_bonds[nbr.idx() as usize])
                        {
                            return Ok(false);
                        }
                    }
                }
            }
        }
        Ok(true)
    }

    /// Ranks the atoms with the current compare codes and bond symbols
    /// (`updateRanks`): every atom and bond takes part, ties are not broken.
    fn rank(&self, mol: &mut ROMol) -> Result<Vec<u32>> {
        let all_atoms = vec![true; mol.num_atoms() as usize];
        let all_bonds = vec![true; mol.num_bonds() as usize];
        rank_fragment_atoms(
            mol,
            &all_atoms,
            &all_bonds,
            &FragmentRankOptions {
                atom_symbols: None,
                bond_symbols: Some(&self.bond_symbols),
                break_ties: false,
                include_chirality: false,
                include_isotopes: true,
                include_atom_maps: true,
                include_chiral_presence: true,
                include_ring_stereo: false,
                compare_codes: Some(&self.atom_codes),
            },
        )
    }

    /// One pass of rank / update until nothing changes.
    fn iterate(&mut self, mol: &mut ROMol, sinfos: &mut Vec<StereoInfo>) -> Result<Vec<u32>> {
        loop {
            sinfos.clear();
            let aranks = self.rank(mol)?;
            let mut again = self.update_atoms(mol, &aranks, sinfos)?;
            again |= self.update_bonds(mol, &aranks, sinfos)?;
            if !again {
                return Ok(aranks);
            }
        }
    }

    /// `cleanMolStereo`: removes the stereo specifications of elements that
    /// cannot be stereo.
    fn clean_mol_stereo(&self, mol: &mut ROMol) -> Result<()> {
        for aidx in 0..mol.num_atoms() {
            let i = aidx as usize;
            if self.fixed_atoms[i] || !self.known_atoms[i] {
                continue;
            }
            match mol.atom(aidx).chiral_tag() {
                ChiralType::TetrahedralCCW | ChiralType::TetrahedralCW => {
                    mol.atom_mut(aidx).set_chiral_tag(ChiralType::Unspecified);
                    for b in mol.atom_bond_indices(aidx).collect::<Vec<_>>() {
                        if matches!(
                            mol.bond(b).bond_dir(),
                            BondDir::BeginDash | BondDir::BeginWedge
                        ) {
                            mol.bond_mut(b).set_bond_dir(BondDir::None);
                        }
                    }
                }
                ChiralType::Tetrahedral
                | ChiralType::SquarePlanar
                | ChiralType::TrigonalBipyramidal
                | ChiralType::Octahedral => {
                    mol.atom_mut(aidx)
                        .set_prop(cp::_chiralPermutation, 0i32, false)?;
                }
                _ => {}
            }
        }
        let mut removed_stereo = false;
        for bidx in 0..mol.num_bonds() {
            let i = bidx as usize;
            if !self.fixed_bonds[i] && self.known_bonds[i] {
                let bond = mol.bond_mut(bidx);
                bond.set_stereo(BondStereo::StereoNone)?;
                bond.set_bond_dir(BondDir::None);
                bond.stereo_atoms_mut().clear();
                removed_stereo = true;
            }
        }
        // remove any slash bond directions without a stereo neighbour
        if removed_stereo {
            for bidx in 0..mol.num_bonds() {
                let dir = mol.bond(bidx).bond_dir();
                if !matches!(dir, BondDir::EndDownRight | BondDir::EndUpRight) {
                    continue;
                }
                let (beg, end) = (
                    mol.bond(bidx).begin_atom_idx(),
                    mol.bond(bidx).end_atom_idx(),
                );
                let mut dir_ok = false;
                for atom in [beg, end] {
                    dir_ok = mol
                        .atom_bonds(atom)
                        .any(|n| n.idx() != bidx && n.stereo() != BondStereo::StereoNone);
                    if dir_ok {
                        break;
                    }
                    mol.bond_mut(bidx).set_bond_dir(BondDir::None);
                }
                let _ = dir_ok;
            }
        }
        Ok(())
    }
}

/// `runCleanup`: up to two passes of "canonicalization" identify the stereo
/// elements. With `clean_it` the first pass ignores possible stereo and
/// removes the specification from every atom/bond without unique neighbours;
/// with `flag_possible` a pass treats each unspecified possible stereo center
/// as different from all others, which finds every possible element.
fn run_cleanup(mol: &mut ROMol, flag_possible: bool, clean_it: bool) -> Result<Vec<StereoInfo>> {
    let mut state = Perception::new(mol);
    state.init_atom_info(mol, flag_possible, clean_it)?;
    state.init_bond_info(mol, flag_possible, clean_it)?;
    // the original sets of possible elements are needed for the second pass
    let orig_possible_atoms = state.possible_atoms.clone();
    let orig_possible_bonds = state.possible_bonds.clone();
    state.flag_ring_stereo(mol, !clean_it)?;

    let mut res = Vec::new();
    if mol.num_atoms() == 0 {
        return Ok(res);
    }
    let mut aranks = state.iterate(mol, &mut res)?;

    if clean_it {
        state.clean_mol_stereo(mol)?;
    }
    if flag_possible
        && (state.possible_atoms != orig_possible_atoms
            || state.possible_bonds != orig_possible_bonds)
    {
        // some possible elements were removed: flag every center/bond where
        // stereo was removed as possible and iterate again
        state.possible_atoms = orig_possible_atoms;
        for i in 0..mol.num_atoms() as usize {
            if !state.fixed_atoms[i] && state.known_atoms[i] {
                state.possible_atoms[i] = true;
                state.known_atoms[i] = false;
            }
            if state.possible_atoms[i] {
                set_atom_compare_state(
                    &mut state.atom_codes[i],
                    AtomCompareState::Possible,
                    i as u32,
                );
            }
        }
        state.possible_bonds = orig_possible_bonds;
        for i in 0..mol.num_bonds() as usize {
            if !state.fixed_bonds[i] && state.known_bonds[i] {
                state.possible_bonds[i] = true;
                state.known_bonds[i] = false;
            }
            if state.possible_bonds[i] {
                state.bond_symbols[i] += &format!("_{i}");
            }
        }
        state.flag_ring_stereo(mol, true)?;
        aranks = state.iterate(mol, &mut res)?;
    }

    let mut special_cases = vec![false; mol.num_atoms() as usize];
    find_chiral_atom_special_cases(mol, &mut special_cases, &aranks)?;
    for (i, &rank) in aranks.iter().enumerate() {
        mol.atom_mut(i as u32)
            .set_prop(cp::_ChiralAtomRank, rank, true)?;
    }
    Ok(res)
}

/// `Chirality::findPotentialStereo(mol, cleanIt, flagPossible)`: identifies
/// the potential (and specified) stereo atoms and bonds of a molecule.
/// With `clean_it` the chirality/stereo specifications are removed from the
/// atoms and bonds that cannot be stereo.
pub fn find_potential_stereo(
    mol: &mut ROMol,
    clean_it: bool,
    flag_possible: bool,
) -> Result<Vec<StereoInfo>> {
    if !mol.ring_info().is_symm_sssr() {
        symmetrize_sssr(mol, SymmetrizeSSSRAlgorithm::Default, true, false, false)?;
    }
    if mol.needs_update_property_cache() {
        mol.update_property_cache(false)?;
    }
    run_cleanup(mol, flag_possible, clean_it)
}

/// The non-modifying overload `findPotentialStereo(const ROMol &)`: the
/// molecule is analysed as it is, on a copy.
pub fn find_potential_stereo_of(mol: &ROMol) -> Result<Vec<StereoInfo>> {
    find_potential_stereo(&mut mol.clone(), false, true)
}
