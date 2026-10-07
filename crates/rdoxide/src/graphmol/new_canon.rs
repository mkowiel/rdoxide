//! Canonical atom ranking — ported from RDKit `Code/GraphMol/new_canon.{h,cpp}`
//! and `Code/RDGeneral/hanoiSort.h` (Schneider, Sayle & Landrum,
//! J. Chem. Inf. Model. 2015, 55, 2111).
//!
//! Copyright (C) 2014-2025 Greg Landrum, Nadine Schneider and other RDKit
//! contributors. Licensed under the BSD 3-Clause license
//! (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::{BTreeSet, VecDeque};

use super::atom::ChiralType;
use super::bond::{BondStereo, BondType};
use super::find_rings::fast_find_rings;
use super::romol::ROMol;
use super::stereo_group::StereoGroupType;
use crate::Result;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::utils::count_swaps_to_interconvert;
use crate::rdgeneral::{RDProps, RDValue};

const ATNUM_CLASS_OFFSET: u32 = 10000;

/// `Canon::bondholder`
#[derive(Debug, Clone)]
pub struct BondHolder {
    pub bond_type: BondType,
    pub bond_stereo: u32,
    pub nbr_sym_class: u32,
    pub nbr_idx: u32,
    pub stype: BondStereo,
    /// canon-atom indices whose current `index` drives cis/trans comparison
    pub controlling_atoms: [Option<u32>; 4],
    pub symbol: Option<String>,
    pub bond_idx: u32,
}

impl BondHolder {
    fn new(
        bond_type: BondType,
        bond_stereo: u32,
        nbr_idx: u32,
        nbr_sym_class: u32,
        bond_idx: u32,
    ) -> Self {
        BondHolder {
            bond_type,
            bond_stereo,
            nbr_sym_class,
            nbr_idx,
            stype: BondStereo::StereoNone,
            controlling_atoms: [None; 4],
            symbol: None,
            bond_idx,
        }
    }

    fn flip_if_needed(st: &mut BondStereo, ctrl: &[Option<u32>; 4], atoms: &[CanonAtom]) {
        let idx = |o: Option<u32>| o.map(|a| atoms[a as usize].index);
        let mut flip = false;
        if let (Some(c1), Some(c0)) = (idx(ctrl[1]), idx(ctrl[0])) {
            if c1 > c0 {
                flip = !flip;
            }
        }
        if let (Some(c3), Some(c2)) = (idx(ctrl[3]), idx(ctrl[2])) {
            if c3 > c2 {
                flip = !flip;
            }
        }
        if flip {
            *st = match *st {
                BondStereo::StereoCis => BondStereo::StereoTrans,
                BondStereo::StereoTrans => BondStereo::StereoCis,
                s => s,
            };
        }
    }

    /// `bondholder::compareStereo`
    fn compare_stereo(&self, o: &BondHolder, atoms: &[CanonAtom]) -> i32 {
        let (mut st1, mut st2) = (self.stype, o.stype);
        use BondStereo::*;
        if st1 == StereoNone {
            return if st2 == StereoNone { 0 } else { -1 };
        }
        if st2 == StereoNone {
            return 1;
        }
        if st1 == StereoAny {
            return if st2 == StereoAny { 0 } else { -1 };
        }
        if st2 == StereoAny {
            return 1;
        }
        if matches!(st1, StereoE | StereoZ) && matches!(st2, StereoE | StereoZ) {
            return (st1 as i32 - st2 as i32).signum();
        }
        Self::flip_if_needed(&mut st1, &self.controlling_atoms, atoms);
        Self::flip_if_needed(&mut st2, &o.controlling_atoms, atoms);
        (st1 as i32 - st2 as i32).signum()
    }

    /// `bondholder::compare(x, y, div)`
    fn compare(x: &BondHolder, y: &BondHolder, atoms: &[CanonAtom], div: u32) -> i32 {
        if let (Some(a), Some(b)) = (&x.symbol, &y.symbol) {
            match a.cmp(b) {
                std::cmp::Ordering::Less => return -1,
                std::cmp::Ordering::Greater => return 1,
                std::cmp::Ordering::Equal => {}
            }
        }
        if x.bond_type != y.bond_type {
            return if x.bond_type < y.bond_type { -1 } else { 1 };
        }
        if x.bond_stereo != y.bond_stereo {
            return if x.bond_stereo < y.bond_stereo { -1 } else { 1 };
        }
        // unsigned difference converted to int, as upstream
        let scdiv = (x.nbr_sym_class / div).wrapping_sub(y.nbr_sym_class / div) as i32;
        if scdiv != 0 {
            return scdiv;
        }
        if x.bond_stereo != 0 && y.bond_stereo != 0 {
            let cs = x.compare_stereo(y, atoms);
            if cs != 0 {
                return cs;
            }
        }
        0
    }

    fn greater(x: &BondHolder, y: &BondHolder, atoms: &[CanonAtom]) -> bool {
        Self::compare(x, y, atoms, 1) > 0
    }
}

/// `Canon::canon_atom`
#[derive(Debug, Clone, Default)]
pub struct CanonAtom {
    pub atom: u32,
    pub index: i32,
    pub degree: u32,
    pub total_num_hs: u32,
    pub has_ring_nbr: bool,
    pub is_ring_stereo_atom: bool,
    pub which_stereo_group: u32,
    pub type_of_stereo_group: u8,
    pub nbr_ids: Vec<i32>,
    pub symbol: Option<String>,
    pub neighbor_num: Vec<i32>,
    pub revisited_neighbors: Vec<i32>,
    pub bonds: Vec<BondHolder>,
}

/// Stable descending sort of bond holders (`std::sort(..., greater)`; the
/// small sorts used by libstdc++/libc++ keep equal elements in order).
fn sort_bonds(bonds: &mut [BondHolder], atoms: &[CanonAtom]) {
    bonds.sort_by(|a, b| BondHolder::compare(b, a, atoms, 1).cmp(&0));
}

/// `updateAtomNeighborIndex`: refresh neighbor classes and re-sort (insertion sort).
fn update_atom_neighbor_index(atoms: &mut [CanonAtom], i: usize) {
    let mut nbrs = std::mem::take(&mut atoms[i].bonds);
    for nbr in &mut nbrs {
        nbr.nbr_sym_class = atoms[nbr.nbr_idx as usize].index as u32;
    }
    for k in 1..nbrs.len() {
        if !BondHolder::greater(&nbrs[k], &nbrs[k - 1], atoms) {
            continue;
        }
        // insertion of nbrs[k] into the sorted prefix
        let value = nbrs.remove(k);
        let mut j = k;
        while j > 0 && BondHolder::greater(&value, &nbrs[j - 1], atoms) {
            j -= 1;
        }
        nbrs.insert(j, value);
    }
    atoms[i].bonds = nbrs;
}

/// `updateAtomNeighborNumSwaps`
fn update_atom_neighbor_num_swaps(
    mol: &ROMol,
    atoms: &[CanonAtom],
    atom_idx: usize,
) -> Vec<(u32, u32)> {
    let mut result = Vec::new();
    let is_ring_atom = mol
        .ring_info()
        .num_atom_rings(atoms[atom_idx].atom)
        .unwrap_or(0)
        != 0;
    for nbr in &atoms[atom_idx].bonds {
        let nidx = nbr.nbr_idx as usize;
        let ntag = mol.atom(atoms[nidx].atom).chiral_tag();
        if is_ring_atom && ntag != ChiralType::Unspecified {
            let mut seen: Vec<i32> = Vec::new();
            let mut too_many = false;
            let mut refv: Vec<i32> = Vec::new();
            for k in 0..atoms[nidx].degree as usize {
                let nn = atoms[nidx].nbr_ids[k];
                refv.push(nn);
                if atom_idx as i32 != nn {
                    let ix = atoms[nn as usize].index;
                    if seen.contains(&ix) {
                        too_many = true;
                    } else {
                        seen.push(ix);
                    }
                }
            }
            let mut probe = vec![atom_idx as i32];
            for b in &atoms[nidx].bonds {
                if b.nbr_idx as usize != atom_idx {
                    probe.push(b.nbr_idx as i32);
                }
            }
            if too_many {
                result.push((nbr.nbr_sym_class, 0));
            } else {
                let n_swaps = count_swaps_to_interconvert(&refv, &probe).unwrap_or(0);
                match ntag {
                    ChiralType::TetrahedralCW => {
                        result.push((nbr.nbr_sym_class, if n_swaps % 2 == 1 { 2 } else { 1 }))
                    }
                    ChiralType::TetrahedralCCW => {
                        result.push((nbr.nbr_sym_class, if n_swaps % 2 == 1 { 1 } else { 2 }))
                    }
                    _ => {}
                }
            }
        } else {
            result.push((nbr.nbr_sym_class, 0));
        }
    }
    result.sort();
    result
}

fn compare_bond_lists(atoms: &[CanonAtom], i: usize, j: usize, div: u32) -> i32 {
    for (a, b) in atoms[i].bonds.iter().zip(&atoms[j].bonds) {
        let c = BondHolder::compare(a, b, atoms, div);
        if c != 0 {
            return c;
        }
    }
    0
}

fn cmp_len(a: usize, b: usize) -> i32 {
    match a.cmp(&b) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Greater => 1,
        std::cmp::Ordering::Equal => 0,
    }
}

/// The compare-functor interface used by the partition refinement.
trait Comparer {
    fn compare(&mut self, atoms: &mut [CanonAtom], i: i32, j: i32) -> i32;
}

fn in_play(sel: Option<&[bool]>, i: usize) -> bool {
    sel.is_none_or(|s| s[i])
}

// ------------------------------------------------------------------ functors

struct SpecialChirality<'a> {
    mol: &'a ROMol,
    atoms_in_play: Option<&'a [bool]>,
}

impl Comparer for SpecialChirality<'_> {
    fn compare(&mut self, atoms: &mut [CanonAtom], i: i32, j: i32) -> i32 {
        let (i, j) = (i as usize, j as usize);
        if let Some(s) = self.atoms_in_play {
            if !(s[i] || s[j]) {
                return 0;
            }
        }
        if in_play(self.atoms_in_play, i) {
            update_atom_neighbor_index(atoms, i);
        }
        if in_play(self.atoms_in_play, j) {
            update_atom_neighbor_index(atoms, j);
        }
        let c = compare_bond_lists(atoms, i, j, 1);
        if c != 0 {
            return c;
        }
        let si = if in_play(self.atoms_in_play, i) {
            update_atom_neighbor_num_swaps(self.mol, atoms, i)
        } else {
            Vec::new()
        };
        let sj = if in_play(self.atoms_in_play, j) {
            update_atom_neighbor_num_swaps(self.mol, atoms, j)
        } else {
            Vec::new()
        };
        for (a, b) in si.iter().zip(&sj) {
            let c = a.1.wrapping_sub(b.1) as i32;
            if c != 0 {
                return c;
            }
        }
        0
    }
}

struct SpecialSymmetry<'a> {
    atoms_in_play: Option<&'a [bool]>,
}

impl Comparer for SpecialSymmetry<'_> {
    fn compare(&mut self, atoms: &mut [CanonAtom], i: i32, j: i32) -> i32 {
        let (i, j) = (i as usize, j as usize);
        if let Some(s) = self.atoms_in_play {
            if !(s[i] || s[j]) {
                return 0;
            }
        }
        match atoms[i].neighbor_num.cmp(&atoms[j].neighbor_num) {
            std::cmp::Ordering::Less => return -1,
            std::cmp::Ordering::Greater => return 1,
            _ => {}
        }
        match atoms[i]
            .revisited_neighbors
            .cmp(&atoms[j].revisited_neighbors)
        {
            std::cmp::Ordering::Less => return -1,
            std::cmp::Ordering::Greater => return 1,
            _ => {}
        }
        if in_play(self.atoms_in_play, i) {
            update_atom_neighbor_index(atoms, i);
        }
        if in_play(self.atoms_in_play, j) {
            update_atom_neighbor_index(atoms, j);
        }
        let c = compare_bond_lists(atoms, i, j, 1);
        if c != 0 {
            return c;
        }
        cmp_len(atoms[i].bonds.len(), atoms[j].bonds.len())
    }
}

/// `getChiralRank`
fn get_chiral_rank(mol: &ROMol, atoms: &[CanonAtom], i: usize) -> u32 {
    let at = atoms[i].atom;
    let mut perm: Vec<i32> = Vec::new();
    for nbr in mol.atom_neighbors(at) {
        let rnk = atoms[nbr as usize].index;
        if perm.contains(&rnk) {
            break;
        }
        perm.push(rnk);
    }
    if perm.len() as u32 != mol.atom_degree(at) {
        return 0;
    }
    let ctag = mol.atom(at).chiral_tag();
    if matches!(ctag, ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW) {
        let mut sorted = perm.clone();
        sorted.sort();
        let nswaps = count_swaps_to_interconvert(&perm, &sorted).unwrap_or(0);
        let mut res = if ctag == ChiralType::TetrahedralCW {
            2
        } else {
            1
        };
        if nswaps % 2 == 1 {
            res = if res == 2 { 1 } else { 2 };
        }
        return res;
    }
    0
}

/// `Canon::AtomCompareFunctor`
pub struct AtomCompareFunctor<'a> {
    mol: &'a ROMol,
    atoms_in_play: Option<&'a [bool]>,
    pub use_nbrs: bool,
    pub use_isotopes: bool,
    pub use_chirality: bool,
    pub use_chirality_rings: bool,
    pub use_atom_maps: bool,
    pub use_non_stereo_ranks: bool,
    pub use_chiral_presence: bool,
    pub use_atom_maps_on_dummies: bool,
}

impl std::fmt::Debug for AtomCompareFunctor<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AtomCompareFunctor").finish_non_exhaustive()
    }
}

impl<'a> AtomCompareFunctor<'a> {
    pub fn new(mol: &'a ROMol, atoms_in_play: Option<&'a [bool]>) -> Self {
        AtomCompareFunctor {
            mol,
            atoms_in_play,
            use_nbrs: false,
            use_isotopes: true,
            use_chirality: true,
            use_chirality_rings: true,
            use_atom_maps: true,
            use_non_stereo_ranks: false,
            use_chiral_presence: true,
            use_atom_maps_on_dummies: true,
        }
    }

    fn atom_ring_nbr_code(&self, atoms: &[CanonAtom], i: usize) -> u32 {
        if !atoms[i].has_ring_nbr {
            return 0;
        }
        let mut code: u32 = 0;
        for k in 0..atoms[i].degree as usize {
            let n = atoms[i].nbr_ids[k] as usize;
            if atoms[n].is_ring_stereo_atom {
                code =
                    code.wrapping_add((atoms[n].index as u32).wrapping_mul(10000).wrapping_add(1));
            }
        }
        code
    }

    fn int_prop(&self, atom: u32, key: &str) -> i32 {
        self.mol
            .atom(atom)
            .get_prop(key)
            .and_then(RDValue::as_int)
            .unwrap_or(0)
    }

    fn basecomp(&self, atoms: &[CanonAtom], i: usize, j: usize) -> i32 {
        macro_rules! cmp {
            ($a:expr, $b:expr) => {{
                let (a, b) = ($a, $b);
                if a < b {
                    return -1;
                } else if a > b {
                    return 1;
                }
            }};
        }
        let mol = self.mol;
        let (ai, aj) = (mol.atom(atoms[i].atom), mol.atom(atoms[j].atom));
        cmp!(atoms[i].index as u32, atoms[j].index as u32);
        if self.use_non_stereo_ranks {
            cmp!(
                self.int_prop(atoms[i].atom, cp::_CanonicalRankingNumber),
                self.int_prop(atoms[j].atom, cp::_CanonicalRankingNumber)
            );
        }
        if self.use_atom_maps || self.use_atom_maps_on_dummies {
            let mi =
                if self.use_atom_maps || (self.use_atom_maps_on_dummies && ai.atomic_num() == 0) {
                    self.int_prop(atoms[i].atom, cp::molAtomMapNumber)
                } else {
                    0
                };
            let mj =
                if self.use_atom_maps || (self.use_atom_maps_on_dummies && aj.atomic_num() == 0) {
                    self.int_prop(atoms[j].atom, cp::molAtomMapNumber)
                } else {
                    0
                };
            cmp!(mi, mj);
        }
        cmp!(atoms[i].degree, atoms[j].degree);
        if let (Some(a), Some(b)) = (&atoms[i].symbol, &atoms[j].symbol) {
            return match a.cmp(b) {
                std::cmp::Ordering::Less => -1,
                std::cmp::Ordering::Greater => 1,
                std::cmp::Ordering::Equal => 0,
            };
        }
        cmp!(ai.atomic_num(), aj.atomic_num());
        if self.use_isotopes {
            cmp!(ai.isotope(), aj.isotope());
        }
        cmp!(atoms[i].total_num_hs, atoms[j].total_num_hs);
        // compared as unsigned ints upstream
        cmp!(ai.formal_charge() as u32, aj.formal_charge() as u32);
        if self.use_chiral_presence {
            cmp!(
                ai.chiral_tag() != ChiralType::Unspecified,
                aj.chiral_tag() != ChiralType::Unspecified
            );
        }
        if self.use_chirality {
            let (gi, gj) = (atoms[i].which_stereo_group, atoms[j].which_stereo_group);
            if gi != 0 || gj != 0 {
                if gi != 0 && gj == 0 {
                    return 1;
                } else if gj != 0 && gi == 0 {
                    return -1;
                } else {
                    cmp!(atoms[i].type_of_stereo_group, atoms[j].type_of_stereo_group);
                    if gi != gj {
                        let groups = mol.stereo_groups();
                        let sgi: BTreeSet<u32> = groups[gi as usize - 1]
                            .atoms()
                            .iter()
                            .map(|&a| atoms[a as usize].index as u32)
                            .collect();
                        let sgj: BTreeSet<u32> = groups[gj as usize - 1]
                            .atoms()
                            .iter()
                            .map(|&a| atoms[a as usize].index as u32)
                            .collect();
                        cmp!(sgi, sgj);
                    } else if atoms[i].type_of_stereo_group == StereoGroupType::Absolute as u8 {
                        cmp!(
                            get_chiral_rank(mol, atoms, i),
                            get_chiral_rank(mol, atoms, j)
                        );
                    }
                }
            } else {
                let ci = ai.chiral_tag() != ChiralType::Unspecified;
                let cj = aj.chiral_tag() != ChiralType::Unspecified;
                cmp!(ci, cj);
                if ci && cj {
                    cmp!(
                        get_chiral_rank(mol, atoms, i),
                        get_chiral_rank(mol, atoms, j)
                    );
                }
            }
        }
        if self.use_chirality_rings {
            cmp!(
                self.atom_ring_nbr_code(atoms, i),
                self.atom_ring_nbr_code(atoms, j)
            );
        }
        0
    }
}

impl Comparer for AtomCompareFunctor<'_> {
    fn compare(&mut self, atoms: &mut [CanonAtom], i: i32, j: i32) -> i32 {
        let (i, j) = (i as usize, j as usize);
        if let Some(s) = self.atoms_in_play {
            if !(s[i] || s[j]) {
                return 0;
            }
        }
        let v = self.basecomp(atoms, i, j);
        if v != 0 {
            return v;
        }
        if self.use_nbrs {
            if in_play(self.atoms_in_play, i) {
                update_atom_neighbor_index(atoms, i);
            }
            if in_play(self.atoms_in_play, j) {
                update_atom_neighbor_index(atoms, j);
            }
            let c = compare_bond_lists(atoms, i, j, 1);
            if c != 0 {
                return c;
            }
            return cmp_len(atoms[i].bonds.len(), atoms[j].bonds.len());
        }
        0
    }
}

/// `Canon::ChiralAtomCompareFunctor`
struct ChiralAtomCompareFunctor<'a> {
    mol: &'a ROMol,
    use_nbrs: bool,
}

impl ChiralAtomCompareFunctor<'_> {
    fn neighborhood(&self, atoms: &mut [CanonAtom], i: usize) {
        let mut nbrs = std::mem::take(&mut atoms[i].bonds);
        for b in &mut nbrs {
            if b.nbr_idx == ATNUM_CLASS_OFFSET {
                continue;
            }
            let n = &atoms[b.nbr_idx as usize];
            b.nbr_sym_class =
                self.mol.atom(n.atom).atomic_num() * ATNUM_CLASS_OFFSET + n.index as u32 + 1;
        }
        sort_bonds(&mut nbrs, atoms);
        atoms[i].bonds = nbrs;
    }

    fn basecomp(&self, atoms: &[CanonAtom], i: usize, j: usize) -> i32 {
        let (ai, aj) = (self.mol.atom(atoms[i].atom), self.mol.atom(atoms[j].atom));
        let pairs = [
            (atoms[i].index as u32, atoms[j].index as u32),
            (ai.atomic_num(), aj.atomic_num()),
            (ai.isotope(), aj.isotope()),
        ];
        for (a, b) in pairs {
            if a != b {
                return if a < b { -1 } else { 1 };
            }
        }
        let cip = |a: &super::atom::Atom| match a.get_prop(cp::_CIPCode) {
            Some(v) => {
                if v.to_rd_string() == "R" {
                    2
                } else {
                    1
                }
            }
            None => 0,
        };
        cip(ai).cmp(&cip(aj)) as i32
    }
}

impl Comparer for ChiralAtomCompareFunctor<'_> {
    fn compare(&mut self, atoms: &mut [CanonAtom], i: i32, j: i32) -> i32 {
        let (i, j) = (i as usize, j as usize);
        let v = self.basecomp(atoms, i, j);
        if v != 0 {
            return v;
        }
        if self.use_nbrs {
            self.neighborhood(atoms, i);
            self.neighborhood(atoms, j);
            let c = compare_bond_lists(atoms, i, j, ATNUM_CLASS_OFFSET);
            if c != 0 {
                return c;
            }
            let c = compare_bond_lists(atoms, i, j, 1);
            if c != 0 {
                return c;
            }
            return cmp_len(atoms[i].bonds.len(), atoms[j].bonds.len());
        }
        0
    }
}

trait UseNbrs {
    fn set_use_nbrs(&mut self, v: bool);
}
impl UseNbrs for AtomCompareFunctor<'_> {
    fn set_use_nbrs(&mut self, v: bool) {
        self.use_nbrs = v;
    }
}
impl UseNbrs for ChiralAtomCompareFunctor<'_> {
    fn set_use_nbrs(&mut self, v: bool) {
        self.use_nbrs = v;
    }
}

// ------------------------------------------------------------------ hanoi sort

/// `RDKit::detail::hanoi`: merge sort over `mem[base..base+nel]` using
/// `mem[temp..temp+nel]` as scratch, recording equivalence-class sizes in
/// `count`. Returns true when the result ended up in the scratch area.
fn hanoi(
    mem: &mut [i32],
    base: usize,
    nel: usize,
    temp: usize,
    count: &mut [i32],
    changed: &[i32],
    cmp: &mut dyn FnMut(i32, i32) -> i32,
) -> bool {
    if nel == 1 {
        count[mem[base] as usize] = 1;
        return false;
    }
    if nel == 2 {
        let (n1, n2) = (mem[base], mem[base + 1]);
        let stat = if changed[n1 as usize] != 0 || changed[n2 as usize] != 0 {
            cmp(n1, n2)
        } else {
            0
        };
        if stat == 0 {
            count[n1 as usize] = 2;
            count[n2 as usize] = 0;
        } else if stat < 0 {
            count[n1 as usize] = 1;
            count[n2 as usize] = 1;
        } else {
            count[n1 as usize] = 1;
            count[n2 as usize] = 1;
            mem[base] = n2;
            mem[base + 1] = n1;
        }
        return false;
    }
    let mut n1 = nel / 2;
    let mut n2 = nel - n1;
    let (b1, t1, b2, t2) = (base, temp, base + n1, temp + n1);
    let (mut s1, mut s2, mut ptr, result);
    if hanoi(mem, b1, n1, t1, count, changed, cmp) {
        s2 = if hanoi(mem, b2, n2, t2, count, changed, cmp) {
            t2
        } else {
            b2
        };
        result = false;
        ptr = base;
        s1 = t1;
    } else {
        s2 = if hanoi(mem, b2, n2, t2, count, changed, cmp) {
            t2
        } else {
            b2
        };
        result = true;
        ptr = temp;
        s1 = b1;
    }
    loop {
        let (a1, a2) = (mem[s1], mem[s2]);
        let stat = if changed[a1 as usize] != 0 || changed[a2 as usize] != 0 {
            cmp(a1, a2)
        } else {
            0
        };
        let len1 = count[a1 as usize] as usize;
        let len2 = count[a2 as usize] as usize;
        if stat == 0 {
            count[a1 as usize] = (len1 + len2) as i32;
            count[a2 as usize] = 0;
            mem.copy_within(s1..s1 + len1, ptr);
            ptr += len1;
            n1 -= len1;
            if n1 == 0 {
                if ptr != s2 {
                    mem.copy_within(s2..s2 + n2, ptr);
                }
                return result;
            }
            s1 += len1;
            mem.copy_within(s2..s2 + len2, ptr);
            ptr += len2;
            n2 -= len2;
            if n2 == 0 {
                mem.copy_within(s1..s1 + n1, ptr);
                return result;
            }
            s2 += len2;
        } else if stat < 0 && len1 > 0 {
            mem.copy_within(s1..s1 + len1, ptr);
            ptr += len1;
            n1 -= len1;
            if n1 == 0 {
                if ptr != s2 {
                    mem.copy_within(s2..s2 + n2, ptr);
                }
                return result;
            }
            s1 += len1;
        } else if stat > 0 && len2 > 0 {
            mem.copy_within(s2..s2 + len2, ptr);
            ptr += len2;
            n2 -= len2;
            if n2 == 0 {
                mem.copy_within(s1..s1 + n1, ptr);
                return result;
            }
            s2 += len2;
        } else {
            unreachable!("hanoi: bad partition counts");
        }
    }
}

// ------------------------------------------------------------------ partitions

struct Partitions {
    order: Vec<i32>,
    count: Vec<i32>,
    next: Vec<i32>,
    changed: Vec<i32>,
    touched: Vec<bool>,
    activeset: i32,
}

impl Partitions {
    fn new(n: usize) -> Self {
        Partitions {
            order: vec![0; n],
            count: vec![0; n],
            next: vec![0; n],
            changed: vec![1; n],
            touched: vec![false; n],
            activeset: -1,
        }
    }

    /// `CreateSinglePartition`
    fn create_single(&mut self, atoms: &mut [CanonAtom]) {
        let n = atoms.len();
        for i in 0..n {
            atoms[i].index = 0;
            self.order[i] = i as i32;
            self.count[i] = 0;
        }
        self.count[0] = n as i32;
    }

    /// `ActivatePartitions`
    fn activate(&mut self) {
        let n = self.order.len();
        self.activeset = -1;
        for x in self.next.iter_mut() {
            *x = -2;
        }
        let mut i = 0;
        loop {
            let j = self.order[i] as usize;
            if self.count[j] > 1 {
                self.next[j] = self.activeset;
                self.activeset = j as i32;
                i += self.count[j] as usize;
            } else {
                i += 1;
            }
            if i >= n {
                break;
            }
        }
        for i in 0..n {
            let j = self.order[i] as usize;
            self.changed[j] = 1;
        }
    }

    fn queue_touched(&mut self) {
        for ii in 0..self.order.len() {
            if self.touched[ii] {
                let partition = self.order[ii] as usize;
                if self.count[partition] > 1 && self.next[partition] == -2 {
                    self.next[partition] = self.activeset;
                    self.activeset = partition as i32;
                }
                self.touched[ii] = false;
            }
        }
    }

    /// `RefinePartitions`
    fn refine(&mut self, atoms: &mut [CanonAtom], compar: &mut dyn Comparer, mode: bool) {
        let mut symclass: i32 = 0;
        while self.activeset != -1 {
            let partition = self.activeset as usize;
            self.activeset = self.next[partition];
            self.next[partition] = -2;
            let len = self.count[partition] as usize;
            let offset = atoms[partition].index as usize;
            let mut mem: Vec<i32> = Vec::with_capacity(2 * len);
            mem.extend_from_slice(&self.order[offset..offset + len]);
            mem.resize(2 * len, 0);
            let in_temp = {
                let mut f = |a: i32, b: i32| compar.compare(atoms, a, b);
                hanoi(
                    &mut mem,
                    0,
                    len,
                    len,
                    &mut self.count,
                    &self.changed,
                    &mut f,
                )
            };
            let sorted = if in_temp { &mem[len..] } else { &mem[..len] };
            self.order[offset..offset + len].copy_from_slice(sorted);
            for k in 0..len {
                let a = self.order[offset + k] as usize;
                self.changed[a] = 0;
            }
            let first = self.order[offset] as usize;
            for i in self.count[first] as usize..len {
                let index = self.order[offset + i] as usize;
                if self.count[index] != 0 {
                    symclass = (offset + i) as i32;
                }
                atoms[index].index = symclass;
                for k in 0..atoms[index].degree as usize {
                    let nb = atoms[index].nbr_ids[k] as usize;
                    self.changed[nb] = 1;
                }
            }
            if mode {
                let first = self.order[offset] as usize;
                for i in self.count[first] as usize..len {
                    let index = self.order[offset + i] as usize;
                    for k in 0..atoms[index].degree as usize {
                        let nb = atoms[index].nbr_ids[k] as usize;
                        self.touched[atoms[nb].index as usize] = true;
                    }
                }
                self.queue_touched();
            }
        }
    }

    /// `BreakTies`
    fn break_ties(&mut self, atoms: &mut [CanonAtom], compar: &mut dyn Comparer, mode: bool) {
        let n = atoms.len() as isize;
        let mut i: isize = 0;
        while i < n {
            let partition = self.order[i as usize] as usize;
            let old_part = atoms[partition].index;
            while self.count[partition] > 1 {
                let len = self.count[partition];
                let offset = atoms[partition].index + len - 1;
                let index = self.order[offset as usize] as usize;
                atoms[index].index = offset;
                self.count[partition] = len - 1;
                self.count[index] = 1;
                if atoms[index].degree < 1 {
                    continue;
                }
                for k in 0..atoms[index].degree as usize {
                    let nb = atoms[index].nbr_ids[k] as usize;
                    self.touched[atoms[nb].index as usize] = true;
                    self.changed[nb] = 1;
                }
                self.queue_touched();
                self.refine(atoms, compar, mode);
            }
            if atoms[partition].index != old_part {
                i -= 1;
            }
            i += 1;
        }
    }
}

/// `compareRingAtomsConcerningNumNeighbors`
fn compare_ring_atoms_concerning_num_neighbors(atoms: &mut [CanonAtom], mol: &ROMol) {
    let n = atoms.len();
    let ri = mol.ring_info();
    let in_ring = |a: u32| ri.is_initialized() && ri.num_atom_rings(a).unwrap_or(0) >= 1;
    let mut visited = vec![false; n];
    let mut last_level = vec![false; n];
    let mut current_level = vec![false; n];
    let mut revisited = vec![0i32; n];
    let mut visited_idx: Vec<usize> = Vec::new();
    let mut last_level_idx: Vec<usize> = Vec::new();
    let mut modified: Vec<usize> = Vec::new();
    for idx in 0..n {
        if !in_ring(atoms[idx].atom) {
            continue;
        }
        let mut neighbors: VecDeque<usize> = VecDeque::from([idx]);
        let mut neighbor_num: Vec<i32> = Vec::new();
        let mut rev: Vec<i32> = Vec::new();
        for &i in &visited_idx {
            visited[i] = false;
        }
        visited_idx.clear();
        for &i in &last_level_idx {
            last_level[i] = false;
        }
        last_level_idx.clear();
        while !neighbors.is_empty() {
            let mut num_level = 0;
            let mut next_level: Vec<usize> = Vec::new();
            while let Some(nidx) = neighbors.pop_front() {
                if !in_ring(atoms[nidx].atom) {
                    continue;
                }
                last_level[nidx] = true;
                last_level_idx.push(nidx);
                visited[nidx] = true;
                visited_idx.push(nidx);
                for k in 0..atoms[nidx].degree as usize {
                    let iidx = atoms[nidx].nbr_ids[k] as usize;
                    if !visited[iidx] {
                        current_level[iidx] = true;
                        num_level += 1;
                        visited[iidx] = true;
                        visited_idx.push(iidx);
                        next_level.push(iidx);
                    }
                }
            }
            for &i in &next_level {
                for k in 0..atoms[i].degree as usize {
                    let j = atoms[i].nbr_ids[k] as usize;
                    if current_level[j] || last_level[j] {
                        if revisited[j] == 0 {
                            modified.push(j);
                        }
                        revisited[j] += 1;
                    }
                }
            }
            for &i in &last_level_idx {
                last_level[i] = false;
            }
            last_level_idx.clear();
            for &i in &next_level {
                last_level[i] = true;
                last_level_idx.push(i);
            }
            for &i in &next_level {
                current_level[i] = false;
            }
            let mut tmp: Vec<i32> = modified.iter().map(|&i| revisited[i]).collect();
            tmp.sort();
            tmp.push(-1);
            rev.extend(tmp);
            for &i in &modified {
                revisited[i] = 0;
            }
            modified.clear();
            neighbor_num.push(num_level);
            neighbor_num.push(-1);
            neighbors.extend(next_level);
        }
        atoms[idx].neighbor_num = neighbor_num;
        atoms[idx].revisited_neighbors = rev;
    }
}

/// `detail::rankWithFunctor`
fn rank_with_functor<F: Comparer + UseNbrs>(
    ftor: &mut F,
    mol: &ROMol,
    atoms: &mut [CanonAtom],
    break_ties: bool,
    use_special: bool,
    use_chirality: bool,
    include_ring_stereo: bool,
    atoms_in_play: Option<&[bool]>,
) -> Vec<i32> {
    let n = atoms.len();
    let mut p = Partitions::new(n);
    p.create_single(atoms);
    ftor.set_use_nbrs(true);
    p.activate();
    p.refine(atoms, ftor, true);
    let mut ties = p.count.contains(&0);
    if use_chirality && ties && include_ring_stereo {
        let mut sc = SpecialChirality { mol, atoms_in_play };
        p.activate();
        p.refine(atoms, &mut sc, true);
    }
    ties = false;
    let mut sym_ring_atoms = 0u32;
    let mut ring_atoms = 0u32;
    let mut branching = false;
    let ri = mol.ring_info();
    for i in 0..n {
        let oi = p.order[i] as u32;
        if ri.is_initialized() && ri.num_atom_rings(oi).unwrap_or(0) != 0 {
            if p.count[oi as usize] > 2 {
                sym_ring_atoms += p.count[oi as usize] as u32;
            }
            ring_atoms += 1;
            if ri.num_atom_rings(oi).unwrap_or(0) > 1 && p.count[oi as usize] > 1 {
                branching = true;
            }
        }
        if p.count[i] == 0 {
            ties = true;
        }
    }
    if use_special
        && ties
        && ring_atoms > 0
        && (sym_ring_atoms as f32) / (ring_atoms as f32) > 0.5
        && branching
    {
        let mut sf = SpecialSymmetry { atoms_in_play };
        compare_ring_atoms_concerning_num_neighbors(atoms, mol);
        p.activate();
        p.refine(atoms, &mut sf, true);
    }
    if break_ties {
        p.break_ties(atoms, ftor, true);
    }
    p.order
}

// ------------------------------------------------------------------ initialization

fn make_bond_holder(
    mol: &ROMol,
    bond_idx: u32,
    other_idx: u32,
    include_chirality: bool,
) -> BondHolder {
    let bond = mol.bond(bond_idx);
    let stereo = if include_chirality {
        bond.stereo() as u32
    } else {
        0
    };
    let bt = if bond.is_aromatic() {
        BondType::Aromatic
    } else {
        bond.bond_type()
    };
    let mut res = BondHolder::new(bt, stereo, other_idx, 0, bond_idx);
    if include_chirality {
        res.stype = bond.stereo();
        let (b, e) = (bond.begin_atom_idx(), bond.end_atom_idx());
        if matches!(res.stype, BondStereo::StereoCis | BondStereo::StereoTrans) {
            let sa = bond.stereo_atoms();
            res.controlling_atoms[0] = Some(sa[0] as u32);
            res.controlling_atoms[2] = Some(sa[1] as u32);
            if mol.atom_degree(b) > 2 {
                for nbr in mol.atom_neighbors(b) {
                    if nbr != e && nbr != sa[0] as u32 {
                        res.controlling_atoms[1] = Some(nbr);
                    }
                }
            }
            if mol.atom_degree(e) > 2 {
                for nbr in mol.atom_neighbors(e) {
                    if nbr != b && nbr != sa[1] as u32 {
                        res.controlling_atoms[3] = Some(nbr);
                    }
                }
            }
        }
        if matches!(
            res.stype,
            BondStereo::StereoAtropCCW | BondStereo::StereoAtropCW
        ) {
            // Atropisomers::getAtropisomerAtomsAndBonds
            let ends = [b, e];
            for (side, &end) in ends.iter().enumerate() {
                let mut nbr_bonds: Vec<u32> = mol
                    .atom_bond_indices(end)
                    .filter(|&x| x != bond_idx)
                    .collect();
                let other = |bi: u32| mol.bond(bi).other_atom_idx(end).expect("in bond");
                if nbr_bonds.len() == 2 && other(nbr_bonds[1]) < other(nbr_bonds[0]) {
                    nbr_bonds.swap(0, 1);
                }
                if let Some(&first) = nbr_bonds.first() {
                    res.controlling_atoms[side * 2] = Some(other(first));
                }
                if nbr_bonds.len() > 1 {
                    res.controlling_atoms[side * 2 + 1] = Some(other(nbr_bonds[1]));
                }
            }
        }
    }
    res
}

fn has_ring_nbr(mol: &ROMol, at: u32) -> bool {
    mol.atom_neighbors(at).any(|n| {
        let a = mol.atom(n);
        matches!(
            a.chiral_tag(),
            ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW
        ) && a.has_prop(cp::_ringStereoAtoms)
    })
}

fn advanced_init(mol: &ROMol, ca: &mut CanonAtom) {
    let a = mol.atom(ca.atom);
    ca.total_num_hs = mol.total_num_hs(ca.atom, false).unwrap_or(0);
    ca.is_ring_stereo_atom = matches!(
        a.chiral_tag(),
        ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW
    ) && a.has_prop(cp::_ringStereoAtoms);
    ca.has_ring_nbr = has_ring_nbr(mol, ca.atom);
}

/// `initCanonAtoms`
fn init_canon_atoms(
    mol: &ROMol,
    include_chirality: bool,
    include_stereo_groups: bool,
) -> Vec<CanonAtom> {
    let n = mol.num_atoms() as usize;
    let mut atoms: Vec<CanonAtom> = (0..n)
        .map(|i| {
            let at = i as u32;
            let mut ca = CanonAtom {
                atom: at,
                index: i as i32,
                degree: mol.atom_degree(at),
                nbr_ids: mol.atom_neighbors(at).map(|x| x as i32).collect(),
                ..Default::default()
            };
            advanced_init(mol, &mut ca);
            ca.bonds = mol
                .atom_bond_indices(at)
                .map(|b| {
                    make_bond_holder(
                        mol,
                        b,
                        mol.bond(b).other_atom_idx(at).expect("in bond"),
                        include_chirality,
                    )
                })
                .collect();
            ca
        })
        .collect();
    for i in 0..n {
        let mut bonds = std::mem::take(&mut atoms[i].bonds);
        sort_bonds(&mut bonds, &atoms);
        atoms[i].bonds = bonds;
    }
    if include_chirality && include_stereo_groups {
        for (k, sg) in mol.stereo_groups().iter().enumerate() {
            for &a in sg.atoms() {
                atoms[a as usize].which_stereo_group = k as u32 + 1;
                atoms[a as usize].type_of_stereo_group = sg.group_type() as u8;
            }
        }
    }
    atoms
}

/// Ensure ring info for ranking; returns true when it must be cleared later.
fn ensure_rings(mol: &mut ROMol) -> Result<bool> {
    if !mol.ring_info().is_find_fast_or_better() {
        fast_find_rings(mol)?;
        return Ok(true);
    }
    Ok(false)
}

/// Options of `Canon::rankMolAtoms`.
#[derive(Debug, Clone, Copy)]
pub struct RankOptions {
    pub break_ties: bool,
    pub include_chirality: bool,
    pub include_isotopes: bool,
    pub include_atom_maps: bool,
    pub include_chiral_presence: bool,
    pub include_stereo_groups: bool,
    pub use_non_stereo_ranks: bool,
    pub include_ring_stereo: bool,
}

impl Default for RankOptions {
    fn default() -> Self {
        RankOptions {
            break_ties: true,
            include_chirality: true,
            include_isotopes: true,
            include_atom_maps: true,
            include_chiral_presence: false,
            include_stereo_groups: true,
            use_non_stereo_ranks: false,
            include_ring_stereo: true,
        }
    }
}

/// `Canon::rankMolAtoms`
pub fn rank_mol_atoms(mol: &mut ROMol, opts: RankOptions) -> Result<Vec<u32>> {
    if mol.num_atoms() == 0 {
        return Ok(Vec::new());
    }
    let clear_rings = ensure_rings(mol)?;
    let res = {
        let m: &ROMol = mol;
        let mut atoms = init_canon_atoms(m, opts.include_chirality, opts.include_stereo_groups);
        let mut ftor = AtomCompareFunctor::new(m, None);
        ftor.use_isotopes = opts.include_isotopes;
        ftor.use_chirality = opts.include_chirality;
        ftor.use_chirality_rings = opts.include_chirality && opts.include_ring_stereo;
        ftor.use_atom_maps = opts.include_atom_maps;
        ftor.use_non_stereo_ranks = opts.use_non_stereo_ranks;
        ftor.use_chiral_presence = opts.include_chiral_presence;
        let order = rank_with_functor(
            &mut ftor,
            m,
            &mut atoms,
            opts.break_ties,
            true,
            opts.include_chirality,
            opts.include_ring_stereo,
            None,
        );
        let mut res = vec![0u32; atoms.len()];
        for &o in &order {
            res[o as usize] = atoms[o as usize].index as u32;
        }
        res
    };
    if clear_rings {
        mol.ring_info_mut().reset(true);
    }
    Ok(res)
}

/// Options of `Canon::rankFragmentAtoms`.
#[derive(Debug, Clone)]
pub struct FragmentRankOptions<'a> {
    pub atom_symbols: Option<&'a [String]>,
    pub bond_symbols: Option<&'a [String]>,
    pub break_ties: bool,
    pub include_chirality: bool,
    pub include_isotopes: bool,
    pub include_atom_maps: bool,
    pub include_chiral_presence: bool,
    pub include_ring_stereo: bool,
}

impl Default for FragmentRankOptions<'_> {
    fn default() -> Self {
        FragmentRankOptions {
            atom_symbols: None,
            bond_symbols: None,
            break_ties: true,
            include_chirality: true,
            include_isotopes: true,
            include_atom_maps: true,
            include_chiral_presence: false,
            include_ring_stereo: true,
        }
    }
}

/// `Canon::rankFragmentAtoms`
pub fn rank_fragment_atoms(
    mol: &mut ROMol,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    opts: &FragmentRankOptions<'_>,
) -> Result<Vec<u32>> {
    if mol.num_atoms() == 0 {
        return Ok(Vec::new());
    }
    let clear_rings = ensure_rings(mol)?;
    let res = {
        let m: &ROMol = mol;
        let n = m.num_atoms() as usize;
        // initFragmentCanonAtoms
        let mut atoms: Vec<CanonAtom> = (0..n)
            .map(|i| CanonAtom {
                atom: i as u32,
                index: i as i32,
                ..Default::default()
            })
            .collect();
        for i in 0..n {
            if atoms_in_play[i] {
                atoms[i].symbol = opts.atom_symbols.map(|s| s[i].clone());
                advanced_init(m, &mut atoms[i]);
            }
        }
        for b in m.bonds() {
            let (bi, be, ee) = (b.idx(), b.begin_atom_idx(), b.end_atom_idx());
            if !bonds_in_play[bi as usize]
                || !atoms_in_play[be as usize]
                || !atoms_in_play[ee as usize]
            {
                continue;
            }
            atoms[be as usize].nbr_ids.push(ee as i32);
            atoms[be as usize].degree += 1;
            atoms[ee as usize].nbr_ids.push(be as i32);
            atoms[ee as usize].degree += 1;
            let mut h1 = make_bond_holder(m, bi, ee, opts.include_chirality);
            let mut h2 = make_bond_holder(m, bi, be, opts.include_chirality);
            if let Some(bs) = opts.bond_symbols {
                h1.symbol = Some(bs[bi as usize].clone());
                h2.symbol = Some(bs[bi as usize].clone());
            }
            atoms[be as usize].bonds.push(h1);
            atoms[ee as usize].bonds.push(h2);
        }
        for i in 0..n {
            if !atoms_in_play[i] {
                continue;
            }
            atoms[i].total_num_hs += m.atom_degree(i as u32) - atoms[i].degree;
            let mut bonds = std::mem::take(&mut atoms[i].bonds);
            sort_bonds(&mut bonds, &atoms);
            atoms[i].bonds = bonds;
        }
        let mut ftor = AtomCompareFunctor::new(m, Some(atoms_in_play));
        ftor.use_isotopes = opts.include_isotopes;
        ftor.use_chirality = opts.include_chirality;
        ftor.use_atom_maps = opts.include_atom_maps;
        ftor.use_chirality_rings = opts.include_chirality;
        ftor.use_chiral_presence = opts.include_chiral_presence;
        let order = rank_with_functor(
            &mut ftor,
            m,
            &mut atoms,
            opts.break_ties,
            true,
            opts.include_chirality,
            opts.include_ring_stereo,
            Some(atoms_in_play),
        );
        let mut res = vec![0u32; n];
        for &o in &order {
            res[o as usize] = atoms[o as usize].index as u32;
        }
        res
    };
    if clear_rings {
        mol.ring_info_mut().reset(true);
    }
    Ok(res)
}

/// `Canon::chiralRankMolAtoms`
pub fn chiral_rank_mol_atoms(mol: &mut ROMol) -> Result<Vec<u32>> {
    if mol.num_atoms() == 0 {
        return Ok(Vec::new());
    }
    let clear_rings = ensure_rings(mol)?;
    let res = {
        let m: &ROMol = mol;
        let n = m.num_atoms() as usize;
        let mut atoms: Vec<CanonAtom> = (0..n)
            .map(|i| {
                let at = i as u32;
                CanonAtom {
                    atom: at,
                    index: i as i32,
                    degree: m.atom_degree(at),
                    nbr_ids: m.atom_neighbors(at).map(|x| x as i32).collect(),
                    ..Default::default()
                }
            })
            .collect();
        for i in 0..n {
            atoms[i].bonds = get_chiral_bonds(m, i as u32);
        }
        let mut ftor = ChiralAtomCompareFunctor {
            mol: m,
            use_nbrs: false,
        };
        let order = rank_with_functor(&mut ftor, m, &mut atoms, false, false, false, false, None);
        let mut res = vec![0u32; n];
        for &o in &order {
            res[o as usize] = atoms[o as usize].index as u32;
        }
        res
    };
    if clear_rings {
        mol.ring_info_mut().reset(true);
    }
    Ok(res)
}

/// `getChiralBonds`
fn get_chiral_bonds(mol: &ROMol, at: u32) -> Vec<BondHolder> {
    let mut nbrs: Vec<BondHolder> = Vec::new();
    let empty: Vec<CanonAtom> = Vec::new();
    for bi in mol.atom_bond_indices(at) {
        let bond = mol.bond(bi);
        let nbr_idx = bond.other_atom_idx(at).expect("in bond");
        let nbr = mol.atom(nbr_idx);
        let degree_nbr = mol.atom_degree(nbr_idx);
        let stereo = match bond.stereo() {
            BondStereo::StereoZ | BondStereo::StereoCis => 1,
            BondStereo::StereoE | BondStereo::StereoTrans => 2,
            _ => 0,
        };
        let n_reps = if bond.bond_type() == BondType::Double
            && nbr.atomic_num() == 15
            && (degree_nbr == 4 || degree_nbr == 3)
        {
            1
        } else {
            (2.0 * bond.bond_type_as_double().unwrap_or(0.0)).floor() as usize
        };
        let symclass = nbr.atomic_num() * ATNUM_CLASS_OFFSET + nbr_idx + 1;
        let bh = BondHolder::new(BondType::Single, stereo, nbr_idx, symclass, bi);
        // std::lower_bound with operator< (compare < 0)
        let pos = nbrs.partition_point(|x| BondHolder::compare(x, &bh, &empty, 1) < 0);
        for _ in 0..n_reps {
            nbrs.insert(pos, bh.clone());
        }
    }
    nbrs.reverse();
    if !mol.atom(at).needs_update_property_cache() {
        for _ in 0..mol.total_num_hs(at, false).unwrap_or(0) {
            for _ in 0..2 {
                nbrs.push(BondHolder::new(
                    BondType::Single,
                    0,
                    ATNUM_CLASS_OFFSET,
                    ATNUM_CLASS_OFFSET,
                    0,
                ));
            }
        }
    }
    nbrs
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphmol::mol_ops::{sanitize_flags::SANITIZE_ALL, sanitize_mol};
    use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};

    fn mol(s: &str) -> ROMol {
        let p = SmilesParserParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        };
        let mut m = mol_from_smiles(s, &p).unwrap().unwrap();
        sanitize_mol(&mut m, SANITIZE_ALL).unwrap();
        m
    }

    #[test]
    fn ranks_are_a_permutation_and_symmetric() {
        let mut m = mol("OCC(C)O");
        let r = rank_mol_atoms(&mut m, RankOptions::default()).unwrap();
        let mut s = r.clone();
        s.sort();
        assert_eq!(s, (0..5).collect::<Vec<u32>>());
        let mut m = mol("OCC(C)O");
        let nb = rank_mol_atoms(
            &mut m,
            RankOptions {
                break_ties: false,
                ..Default::default()
            },
        )
        .unwrap();
        // the two terminal oxygens are not equivalent here, but invariant ranks
        // must be reproducible
        assert_eq!(
            nb,
            rank_mol_atoms(
                &mut m,
                RankOptions {
                    break_ties: false,
                    ..Default::default()
                }
            )
            .unwrap()
        );
        let mut m = mol("c1ccccc1");
        let nb = rank_mol_atoms(
            &mut m,
            RankOptions {
                break_ties: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(nb.iter().all(|&x| x == nb[0]));
    }

    #[test]
    fn hanoi_sorts_with_counts() {
        let vals: [i32; 5] = [5, 3, 5, 1, 3];
        let mut mem: Vec<i32> = (0..5).collect();
        mem.resize(10, 0);
        let mut count = vec![0; 5];
        let changed = vec![1; 5];
        let mut cmp = |a: i32, b: i32| (vals[a as usize] - vals[b as usize]).signum();
        let in_temp = hanoi(&mut mem, 0, 5, 5, &mut count, &changed, &mut cmp);
        let order: Vec<i32> = if in_temp {
            mem[5..].to_vec()
        } else {
            mem[..5].to_vec()
        };
        let sorted_vals: Vec<i32> = order.iter().map(|&i| vals[i as usize]).collect();
        assert_eq!(sorted_vals, vec![1, 3, 3, 5, 5]);
        assert_eq!(count[order[0] as usize], 1);
        assert_eq!(count[order[1] as usize], 2);
        assert_eq!(count[order[3] as usize], 2);
    }
}
