//! Ported from RDKit `Code/GraphMol/Atom.{h,cpp}`.
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! Upstream's `Atom` keeps a back-pointer to its owning molecule; here atoms
//! live inside [`ROMol`] and the operations that need the molecule (degree,
//! valence, ...) take `(&ROMol, atom index)`.

use super::bond::BondType;
use super::periodic_table::PeriodicTable;
use super::romol::ROMol;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{Dict, RDProps, RDValue};
use crate::{RDKitError, Result, SanitizeError};

/// `Atom::HybridizationType`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
#[repr(u8)]
pub enum HybridizationType {
    #[default]
    Unspecified = 0,
    S,
    SP,
    SP2,
    SP3,
    SP2D,
    SP3D,
    SP3D2,
    Other,
}

impl HybridizationType {
    pub const ALL: [HybridizationType; 9] = [
        Self::Unspecified,
        Self::S,
        Self::SP,
        Self::SP2,
        Self::SP3,
        Self::SP2D,
        Self::SP3D,
        Self::SP3D2,
        Self::Other,
    ];
    /// Name of the Python enum value (`rdchem.HybridizationType.X`).
    pub fn name(self) -> &'static str {
        match self {
            Self::Unspecified => "UNSPECIFIED",
            Self::S => "S",
            Self::SP => "SP",
            Self::SP2 => "SP2",
            Self::SP3 => "SP3",
            Self::SP2D => "SP2D",
            Self::SP3D => "SP3D",
            Self::SP3D2 => "SP3D2",
            Self::Other => "OTHER",
        }
    }
    pub fn from_u8(v: u8) -> Option<Self> {
        Self::ALL.get(v as usize).copied()
    }
}

/// `Atom::ChiralType`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
#[repr(u8)]
pub enum ChiralType {
    #[default]
    Unspecified = 0,
    TetrahedralCW,
    TetrahedralCCW,
    Other,
    Tetrahedral,
    Allene,
    SquarePlanar,
    TrigonalBipyramidal,
    Octahedral,
}

impl ChiralType {
    pub const ALL: [ChiralType; 9] = [
        Self::Unspecified,
        Self::TetrahedralCW,
        Self::TetrahedralCCW,
        Self::Other,
        Self::Tetrahedral,
        Self::Allene,
        Self::SquarePlanar,
        Self::TrigonalBipyramidal,
        Self::Octahedral,
    ];
    pub fn name(self) -> &'static str {
        match self {
            Self::Unspecified => "CHI_UNSPECIFIED",
            Self::TetrahedralCW => "CHI_TETRAHEDRAL_CW",
            Self::TetrahedralCCW => "CHI_TETRAHEDRAL_CCW",
            Self::Other => "CHI_OTHER",
            Self::Tetrahedral => "CHI_TETRAHEDRAL",
            Self::Allene => "CHI_ALLENE",
            Self::SquarePlanar => "CHI_SQUAREPLANAR",
            Self::TrigonalBipyramidal => "CHI_TRIGONALBIPYRAMIDAL",
            Self::Octahedral => "CHI_OCTAHEDRAL",
        }
    }
    pub fn from_u8(v: u8) -> Option<Self> {
        Self::ALL.get(v as usize).copied()
    }
}

/// `Atom::ValenceType`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ValenceType {
    Implicit = 0,
    Explicit = 1,
}

impl ValenceType {
    pub const ALL: [ValenceType; 2] = [ValenceType::Implicit, ValenceType::Explicit];
    pub fn name(self) -> &'static str {
        match self {
            ValenceType::Implicit => "IMPLICIT",
            ValenceType::Explicit => "EXPLICIT",
        }
    }
    pub fn from_u8(v: u8) -> Option<Self> {
        Self::ALL.get(v as usize).copied()
    }
}

/// An atom (`RDKit::Atom`).
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub(crate) atomic_num: u8,
    pub(crate) formal_charge: i8,
    pub(crate) num_explicit_hs: u8,
    pub(crate) no_implicit: bool,
    pub(crate) is_aromatic: bool,
    pub(crate) num_radical_electrons: u8,
    pub(crate) chiral_tag: ChiralType,
    pub(crate) hybrid: HybridizationType,
    pub(crate) isotope: u16,
    pub(crate) implicit_valence: i8,
    pub(crate) explicit_valence: i8,
    pub(crate) idx: u32,
    pub(crate) flags: u64,
    /// whether the atom belongs to a molecule (`dp_mol != nullptr`)
    pub(crate) has_owner: bool,
    pub(crate) props: Dict,
    /// `QueryAtom::dp_query` (a `QueryAtom` is an atom with a query)
    pub(crate) query: Option<Box<crate::graphmol::queries::AtomQuery>>,
}

impl Default for Atom {
    fn default() -> Self {
        Atom::new(0)
    }
}

impl RDProps for Atom {
    fn dict(&self) -> &Dict {
        &self.props
    }
    fn dict_mut(&mut self) -> &mut Dict {
        &mut self.props
    }
}

impl Atom {
    /// `Atom(unsigned int num)`
    pub fn new(atomic_num: u32) -> Self {
        Atom {
            atomic_num: atomic_num as u8,
            formal_charge: 0,
            num_explicit_hs: 0,
            no_implicit: false,
            is_aromatic: false,
            num_radical_electrons: 0,
            chiral_tag: ChiralType::Unspecified,
            hybrid: HybridizationType::Unspecified,
            isotope: 0,
            implicit_valence: -1,
            explicit_valence: -1,
            idx: 0,
            flags: 0,
            has_owner: false,
            props: Dict::new(),
            query: None,
        }
    }

    /// `Atom(const std::string &symbol)`
    pub fn from_symbol(symbol: &str) -> Result<Self> {
        Ok(Atom::new(PeriodicTable::get().get_atomic_number(symbol)?))
    }

    /// `Atom::copy()` / copy-constructor: copies everything except ownership
    /// and index.
    pub fn copy_detached(&self) -> Atom {
        let mut a = self.clone();
        a.has_owner = false;
        a.idx = 0;
        a
    }

    pub fn atomic_num(&self) -> u32 {
        u32::from(self.atomic_num)
    }
    pub fn set_atomic_num(&mut self, n: i32) {
        self.atomic_num = n as u8;
    }
    pub fn formal_charge(&self) -> i32 {
        i32::from(self.formal_charge)
    }
    pub fn set_formal_charge(&mut self, c: i32) {
        self.formal_charge = c as i8;
    }
    pub fn num_explicit_hs(&self) -> u32 {
        u32::from(self.num_explicit_hs)
    }
    pub fn set_num_explicit_hs(&mut self, n: u32) {
        self.num_explicit_hs = n as u8;
    }
    pub fn no_implicit(&self) -> bool {
        self.no_implicit
    }
    pub fn set_no_implicit(&mut self, v: bool) {
        self.no_implicit = v;
    }
    pub fn is_aromatic(&self) -> bool {
        self.is_aromatic
    }
    pub fn set_is_aromatic(&mut self, v: bool) {
        self.is_aromatic = v;
    }
    pub fn num_radical_electrons(&self) -> u32 {
        u32::from(self.num_radical_electrons)
    }
    pub fn set_num_radical_electrons(&mut self, n: u32) {
        self.num_radical_electrons = n as u8;
    }
    pub fn chiral_tag(&self) -> ChiralType {
        self.chiral_tag
    }
    pub fn set_chiral_tag(&mut self, t: ChiralType) {
        self.chiral_tag = t;
    }
    pub fn hybridization(&self) -> HybridizationType {
        self.hybrid
    }
    pub fn set_hybridization(&mut self, h: HybridizationType) {
        self.hybrid = h;
    }
    pub fn isotope(&self) -> u32 {
        u32::from(self.isotope)
    }
    pub fn set_isotope(&mut self, i: u32) {
        self.isotope = i as u16;
    }
    pub fn idx(&self) -> u32 {
        self.idx
    }
    pub fn set_idx(&mut self, i: u32) {
        self.idx = i;
    }
    pub fn flags(&self) -> u64 {
        self.flags
    }
    pub fn set_flags(&mut self, f: u64) {
        self.flags = f;
    }
    pub fn has_owning_mol(&self) -> bool {
        self.has_owner
    }
    pub fn has_query(&self) -> bool {
        self.query.is_some()
    }
    /// `getQuery()`
    pub fn query(&self) -> Option<&crate::graphmol::queries::AtomQuery> {
        self.query.as_deref()
    }
    /// `setQuery()`
    pub fn set_query(&mut self, q: Option<crate::graphmol::queries::AtomQuery>) {
        self.query = q.map(Box::new);
    }
    /// `getQueryType()` (the type label of the query)
    pub fn query_type(&self) -> &str {
        self.query.as_deref().map_or("", |q| q.type_label.as_str())
    }

    /// `getSymbol()`: dummy atoms may carry a `dummyLabel`.
    pub fn symbol(&self) -> String {
        if self.atomic_num == 0 {
            if let Some(RDValue::String(s)) = self.props.get(cp::dummyLabel) {
                return s.clone();
            }
        }
        PeriodicTable::get()
            .element(self.atomic_num())
            .symbol
            .clone()
    }

    /// `getMass()`
    pub fn mass(&self) -> f64 {
        let pt = PeriodicTable::get();
        if self.isotope != 0 {
            let res = pt
                .get_mass_for_isotope(self.atomic_num(), self.isotope())
                .unwrap_or(0.0);
            if self.atomic_num != 0 && res == 0.0 {
                f64::from(self.isotope)
            } else {
                res
            }
        } else {
            pt.element(self.atomic_num()).mass
        }
    }

    pub fn needs_update_property_cache(&self) -> bool {
        !(self.explicit_valence >= 0 && (self.no_implicit || self.implicit_valence >= 0))
    }

    pub fn clear_property_cache(&mut self) {
        self.explicit_valence = -1;
        self.implicit_valence = -1;
    }

    /// Raw cached explicit valence (-1 when not computed).
    pub fn cached_explicit_valence(&self) -> i32 {
        i32::from(self.explicit_valence)
    }

    /// Raw cached implicit valence (-1 when not computed).
    pub fn cached_implicit_valence(&self) -> i32 {
        i32::from(self.implicit_valence)
    }

    /// `getValence(which)`
    pub fn valence(&self, which: ValenceType) -> Result<u32> {
        if !self.has_owner {
            return Ok(0);
        }
        if which == ValenceType::Explicit && self.explicit_valence <= -1 {
            return Err(RDKitError::precondition(
                "getValence(ValenceType::EXPLICIT) called without call to calcExplicitValence()",
            ));
        }
        if which == ValenceType::Implicit && !self.no_implicit && self.implicit_valence <= -1 {
            return Err(RDKitError::precondition(
                "getValence(ValenceType::IMPLICIT) called without call to calcImplicitValence()",
            ));
        }
        Ok(match which {
            ValenceType::Explicit => self.explicit_valence as u32,
            ValenceType::Implicit => {
                if self.no_implicit {
                    0
                } else {
                    self.implicit_valence as u32
                }
            }
        })
    }

    /// `getNumImplicitHs()`
    pub fn num_implicit_hs(&self) -> Result<u32> {
        if self.no_implicit {
            return Ok(0);
        }
        if self.implicit_valence <= -1 {
            return Err(RDKitError::precondition(
                "getNumImplicitHs() called without preceding call to calcImplicitValence()",
            ));
        }
        self.valence(ValenceType::Implicit)
    }

    /// `getTotalValence()`
    pub fn total_valence(&self) -> Result<u32> {
        Ok(self.valence(ValenceType::Explicit)? + self.valence(ValenceType::Implicit)?)
    }

    /// `Atom::Match`: does this (query-less) atom match `what`?
    pub fn matches(&self, what: &Atom) -> bool {
        let mut res = self.atomic_num == what.atomic_num;
        if res {
            if self.atomic_num == 0 {
                let (tgt, test) = (self.isotope, what.isotope);
                if tgt != 0 && test != 0 && tgt != test {
                    res = false;
                }
            } else if (self.formal_charge != 0 && self.formal_charge != what.formal_charge)
                || (self.isotope != 0 && self.isotope != what.isotope)
                || (self.num_radical_electrons != 0
                    && self.num_radical_electrons != what.num_radical_electrons)
            {
                res = false;
            }
        }
        res
    }

    /// `setAtomMapNum`
    pub fn set_atom_map_num(&mut self, mapno: i32, strict: bool) -> Result<()> {
        if strict && !(0..1000).contains(&mapno) {
            return Err(RDKitError::precondition(
                "atom map number out of range [0..1000], use strict=false to override",
            ));
        }
        if mapno != 0 {
            self.props.set_val(cp::molAtomMapNumber, mapno);
        } else {
            self.props.clear_val(cp::molAtomMapNumber);
        }
        Ok(())
    }

    /// `getAtomMapNum`
    pub fn atom_map_num(&self) -> i32 {
        self.props
            .get(cp::molAtomMapNumber)
            .and_then(RDValue::as_int)
            .unwrap_or(0)
    }

    /// `invertChirality()`
    pub fn invert_chirality(&mut self) -> bool {
        const OCTAHEDRAL_INVERT: [u8; 31] = [
            0, 2, 1, 16, 14, 15, 18, 17, 10, 11, 8, 9, 13, 12, 4, 5, 3, 7, 6, 24, 23, 22, 21, 20,
            19, 30, 29, 28, 27, 26, 25,
        ];
        const TBP_INVERT: [u8; 21] = [
            0, 2, 1, 4, 3, 6, 5, 8, 7, 11, 12, 9, 10, 14, 13, 20, 19, 18, 17, 16, 15,
        ];
        match self.chiral_tag {
            ChiralType::TetrahedralCW => {
                self.chiral_tag = ChiralType::TetrahedralCCW;
                true
            }
            ChiralType::TetrahedralCCW => {
                self.chiral_tag = ChiralType::TetrahedralCW;
                true
            }
            tag @ (ChiralType::Tetrahedral
            | ChiralType::TrigonalBipyramidal
            | ChiralType::Octahedral) => {
                let Some(perm) = self
                    .props
                    .get(cp::_chiralPermutation)
                    .and_then(RDValue::as_uint)
                else {
                    return false;
                };
                let new = match tag {
                    ChiralType::Tetrahedral => match perm {
                        1 => 2,
                        2 => 1,
                        _ => 0,
                    },
                    ChiralType::TrigonalBipyramidal => {
                        TBP_INVERT.get(perm as usize).map_or(0, |&p| u32::from(p))
                    }
                    _ => OCTAHEDRAL_INVERT
                        .get(perm as usize)
                        .map_or(0, |&p| u32::from(p)),
                };
                self.props.set_val(cp::_chiralPermutation, new);
                new != 0
            }
            _ => false,
        }
    }
}

/// `setAtomRLabel`
pub fn set_atom_rlabel(atom: &mut Atom, rlabel: i32) -> Result<()> {
    if !(0..100).contains(&rlabel) {
        return Err(RDKitError::precondition_at(
            "rlabel out of range for MDL files",
            "rlabel >= 0 && rlabel < 100",
            "Code/GraphMol/Atom.cpp",
            873,
        ));
    }
    if rlabel != 0 {
        atom.props.set_val(cp::_MolFileRLabel, rlabel as u32);
    } else {
        atom.clear_prop(cp::_MolFileRLabel);
    }
    Ok(())
}

/// `getAtomRLabel`
pub fn get_atom_rlabel(atom: &Atom) -> i32 {
    atom.props
        .get(cp::_MolFileRLabel)
        .and_then(RDValue::as_uint)
        .unwrap_or(0) as i32
}

fn set_or_clear_string(atom: &mut Atom, key: &str, value: &str) {
    if value.is_empty() {
        atom.clear_prop(key);
    } else {
        atom.props.set_val(key, value);
    }
}

fn get_string(atom: &Atom, key: &str) -> String {
    atom.props
        .get(key)
        .map(RDValue::to_rd_string)
        .unwrap_or_default()
}

pub fn set_atom_alias(atom: &mut Atom, alias: &str) {
    set_or_clear_string(atom, cp::molFileAlias, alias);
}
pub fn get_atom_alias(atom: &Atom) -> String {
    get_string(atom, cp::molFileAlias)
}
pub fn set_atom_value(atom: &mut Atom, value: &str) {
    set_or_clear_string(atom, cp::molFileValue, value);
}
pub fn get_atom_value(atom: &Atom) -> String {
    get_string(atom, cp::molFileValue)
}
pub fn set_supplemental_smiles_label(atom: &mut Atom, label: &str) {
    set_or_clear_string(atom, cp::_supplementalSmilesLabel, label);
}
pub fn get_supplemental_smiles_label(atom: &Atom) -> String {
    get_string(atom, cp::_supplementalSmilesLabel)
}

/// `isEarlyAtom`
pub fn is_early_atom(atomic_num: i32) -> bool {
    const EARLY: [u8; 42] = [
        3, 4, 5, 11, 12, 13, 19, 20, 21, 22, 30, 31, 32, 37, 38, 39, 40, 41, 48, 49, 50, 51, 55,
        56, 57, 58, 59, 60, 61, 72, 73, 80, 81, 82, 83, 87, 88, 89, 90, 91, 92, 93,
    ];
    // all of 104..=118 are early
    (0..119).contains(&atomic_num) && (atomic_num >= 104 || EARLY.contains(&(atomic_num as u8)))
}

/// `isAromaticAtom`: aromatic flag, or any aromatic bond.
pub fn is_aromatic_atom(mol: &ROMol, idx: u32) -> bool {
    let atom = mol.atom(idx);
    if atom.is_aromatic {
        return true;
    }
    if atom.has_owner {
        return mol
            .atom_bonds(idx)
            .any(|b| b.is_aromatic() || b.bond_type() == BondType::Aromatic);
    }
    false
}

/// `getEffectiveAtomicNum`
fn effective_atomic_num(atom: &Atom, check_value: bool) -> Result<u32> {
    let max = PeriodicTable::get().get_max_atomic_number() as i32;
    let eff = atom.atomic_num() as i32 - atom.formal_charge();
    if check_value && !(0..=max).contains(&eff) {
        return Err(RDKitError::Sanitize(SanitizeError::AtomValence {
            atom_idx: atom.idx,
            message: "Effective atomic number out of range".into(),
        }));
    }
    Ok(eff.clamp(0, max) as u32)
}

fn can_be_hypervalent(atom: &Atom, eff: u32) -> bool {
    let an = atom.atomic_num();
    (eff > 16 && (an == 15 || an == 16)) || (eff > 34 && (an == 33 || an == 34))
}

fn valence_error(atom: &Atom, message: String) -> RDKitError {
    crate::rd_log_error!("{message}");
    RDKitError::Sanitize(SanitizeError::AtomValence {
        atom_idx: atom.idx,
        message,
    })
}

/// `calculateExplicitValence` (returns -1 on failure when `strict` is false
/// and `check_it` is true).
pub fn calculate_explicit_valence(
    mol: &ROMol,
    idx: u32,
    strict: bool,
    check_it: bool,
) -> Result<i32> {
    let pt = PeriodicTable::get();
    let atom = mol.atom(idx);
    let mut accum: f64 = mol
        .atom_bonds(idx)
        .map(|b| b.valence_contrib(mol, idx))
        .sum();
    accum += f64::from(atom.num_explicit_hs);
    let ovalens = &pt.element(atom.atomic_num()).valences;
    let mut eff = atom.atomic_num();
    if ovalens.len() > 1 || ovalens[0] != -1 {
        eff = effective_atomic_num(atom, check_it)?;
    }
    let dv = pt.element(eff).valences[0];
    let valens = &pt.element(eff).valences;
    if accum > f64::from(dv) && is_aromatic_atom(mol, idx) {
        let mut pval = dv;
        for &val in valens {
            if val == -1 || f64::from(val) > accum {
                break;
            }
            pval = val;
        }
        if accum - f64::from(pval) <= 1.5 {
            accum = f64::from(pval);
        }
    }
    accum += 0.1;
    let res = accum.round() as i32;
    if strict || check_it {
        let mut max_valence = *valens.last().expect("valence list");
        let mut offset = 0;
        if can_be_hypervalent(atom, eff) {
            max_valence = *ovalens.last().expect("valence list");
            offset -= atom.formal_charge();
        }
        if atom.atomic_num == 1 && atom.formal_charge == -1 {
            max_valence = 2;
        }
        if max_valence >= 0
            && *ovalens.last().expect("valence list") >= 0
            && res + offset > max_valence
        {
            if strict {
                return Err(valence_error(
                    atom,
                    format!(
                        "Explicit valence for atom # {} {}, {}, is greater than permitted",
                        atom.idx,
                        pt.element(atom.atomic_num()).symbol,
                        res
                    ),
                ));
            }
            return Ok(-1);
        }
    }
    Ok(res)
}

/// `calculateImplicitValence`
pub fn calculate_implicit_valence(
    mol: &ROMol,
    idx: u32,
    strict: bool,
    check_it: bool,
) -> Result<i32> {
    let atom = mol.atom(idx);
    if atom.no_implicit {
        return Ok(0);
    }
    let pt = PeriodicTable::get();
    let mut explicit_valence = i32::from(atom.explicit_valence);
    if explicit_valence == -1 {
        explicit_valence = calculate_explicit_valence(mol, idx, strict, check_it)?;
    }
    let atomic_num = atom.atomic_num();
    if atomic_num == 0 {
        return Ok(0);
    }
    if mol.atom_bonds(idx).any(|b| b.has_complex_bond_type_query()) {
        return Ok(0);
    }
    let formal_charge = atom.formal_charge();
    if explicit_valence == 0 && atom.num_radical_electrons == 0 && atomic_num == 1 {
        return match formal_charge {
            1 | -1 => Ok(0),
            0 => Ok(1),
            _ => {
                if strict {
                    Err(valence_error(
                        atom,
                        format!("Unreasonable formal charge on atom # {}.", atom.idx),
                    ))
                } else if check_it {
                    Ok(-1)
                } else {
                    Ok(0)
                }
            }
        };
    }
    // NOTE: upstream uses the cached d_explicitValence here
    let mut explicit_plus_rad =
        i32::from(atom.explicit_valence) + i32::from(atom.num_radical_electrons);
    let ovalens = &pt.element(atomic_num).valences;
    let mut eff = atomic_num;
    if ovalens.len() > 1 || ovalens[0] != -1 {
        eff = effective_atomic_num(atom, check_it)?;
    }
    if eff == 0 {
        return Ok(0);
    }
    let dv = pt.element(eff).valences[0];
    if dv == -1 {
        return Ok(0);
    }
    if can_be_hypervalent(atom, eff) {
        eff = atomic_num;
        explicit_plus_rad -= formal_charge;
    }
    let valens = &pt.element(eff).valences;
    let mut res;
    if is_aromatic_atom(mol, idx) {
        if explicit_plus_rad <= dv {
            res = dv - explicit_plus_rad;
        } else {
            let satis = valens
                .iter()
                .take_while(|&&v| v > 0)
                .any(|&v| v == explicit_plus_rad);
            if !satis && (strict || check_it) {
                if strict {
                    return Err(valence_error(
                        atom,
                        format!(
                            "Explicit valence for aromatic atom # {} not equal to any accepted valence\n",
                            atom.idx
                        ),
                    ));
                }
                return Ok(-1);
            }
            res = 0;
        }
    } else {
        res = -1;
        for &tot in valens.iter().take_while(|&&v| v >= 0) {
            if explicit_plus_rad <= tot {
                res = tot - explicit_plus_rad;
                break;
            }
        }
        if res < 0 {
            if (strict || check_it)
                && *valens.last().expect("valences") != -1
                && *ovalens.last().expect("valences") > 0
            {
                if strict {
                    return Err(valence_error(
                        atom,
                        format!(
                            "Explicit valence for atom # {} {} greater than permitted",
                            atom.idx,
                            pt.element(atomic_num).symbol
                        ),
                    ));
                }
                return Ok(-1);
            }
            res = 0;
        }
    }
    Ok(res)
}

/// `Atom::hasValenceViolation`
pub fn has_valence_violation(mol: &ROMol, idx: u32) -> bool {
    let atom = mol.atom(idx);
    if atom.atomic_num == 0 || atom.has_query() || mol.atom_bonds(idx).any(|b| b.has_query()) {
        return false;
    }
    let Ok(eff) = effective_atomic_num(atom, true) else {
        return true;
    };
    let pt = PeriodicTable::get();
    if atom.atomic_num == 1 {
        if atom.formal_charge > 1 || atom.formal_charge < -1 {
            return true;
        }
    } else if atom.formal_charge() > atom.atomic_num() as i32
        || pt.element(atom.atomic_num()).row != pt.element(eff).row
    {
        return true;
    }
    matches!(
        calculate_explicit_valence(mol, idx, false, true),
        Ok(-1) | Err(_)
    ) || matches!(
        calculate_implicit_valence(mol, idx, false, true),
        Ok(-1) | Err(_)
    )
}

/// `numPiElectrons`
pub fn num_pi_electrons(mol: &ROMol, idx: u32) -> Result<u32> {
    let atom = mol.atom(idx);
    if atom.is_aromatic {
        return Ok(1);
    }
    if atom.hybrid == HybridizationType::SP3 {
        return Ok(0);
    }
    let val = atom.valence(ValenceType::Explicit)?;
    let mut physical = atom.num_explicit_hs();
    physical += mol
        .atom_bonds(idx)
        .filter(|b| b.valence_contrib(mol, idx) != 0.0)
        .count() as u32;
    if val < physical {
        return Err(RDKitError::invariant(
            "explicit valence exceeds atom degree",
        ));
    }
    Ok(val - physical)
}
