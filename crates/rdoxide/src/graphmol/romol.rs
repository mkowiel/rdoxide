//! Ported from RDKit `Code/GraphMol/ROMol.{h,cpp}` and `RWMol.{h,cpp}`.
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! Upstream stores the graph in a `boost::adjacency_list<vecS, vecS,
//! undirectedS>`. Iteration orders that are observable from Python
//! (neighbors/bonds of an atom, the bond list) follow boost's semantics:
//! per-atom incidence lists are kept in insertion order, the bond list is in
//! insertion order, and a bond's index always equals its position.
//!
//! Upstream's `ROMol` and `RWMol` are a single type here; [`RWMol`] is an alias.

use std::collections::BTreeMap;

use super::atom::{self, Atom, ValenceType};
use super::bond::{Bond, BondStereo, BondType};
use super::conformer::Conformer;
use super::ring_info::RingInfo;
use super::stereo_group::{self, StereoGroup, StereoGroupType};
use super::substance_group::{
    SubstanceGroup, insert_substance_groups, remove_substance_groups_referencing_atom,
    remove_substance_groups_referencing_bond,
};
use crate::geometry::Point3D;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{Dict, RDProps, RDValue};
use crate::{RDKitError, Result};

/// Bookmark used by parsers for the most recently added atom.
pub const CI_RIGHTMOST_ATOM: i32 = -0xBAD_BEEF;
pub const CI_LEADING_BOND: i32 = -0xBAD_BEEF + 1;
pub const CI_ATOM_HOLDER: i32 = -0xDEA_DD06;

/// A molecule (`RDKit::ROMol` / `RDKit::RWMol`).
#[derive(Debug, Clone, Default)]
pub struct ROMol {
    pub(crate) atoms: Vec<Atom>,
    pub(crate) bonds: Vec<Bond>,
    /// per atom: (neighbor atom index, bond index), in insertion order
    pub(crate) adj: Vec<Vec<(u32, u32)>>,
    pub(crate) ring_info: RingInfo,
    pub(crate) conformers: Vec<Conformer>,
    pub(crate) props: Dict,
    pub(crate) stereo_groups: Vec<StereoGroup>,
    pub(crate) sgroups: Vec<SubstanceGroup>,
    pub(crate) atom_bookmarks: BTreeMap<i32, Vec<u32>>,
    pub(crate) bond_bookmarks: BTreeMap<i32, Vec<u32>>,
    pub(crate) del_atoms: Option<Vec<bool>>,
    pub(crate) del_bonds: Option<Vec<bool>>,
}

/// Read/write molecule; same type as [`ROMol`].
pub type RWMol = ROMol;

impl RDProps for ROMol {
    fn dict(&self) -> &Dict {
        &self.props
    }
    fn dict_mut(&mut self) -> &mut Dict {
        &mut self.props
    }
}

impl ROMol {
    pub fn new() -> Self {
        Self::default()
    }

    /// Copy constructor `ROMol(other, quickCopy, confId)`.
    pub fn copy_from(other: &ROMol, quick_copy: bool, conf_id: i32) -> ROMol {
        let mut m = ROMol {
            atoms: other.atoms.clone(),
            bonds: other.bonds.clone(),
            adj: other.adj.clone(),
            ring_info: other.ring_info.clone(),
            stereo_groups: other.stereo_groups.clone(),
            del_atoms: other.del_atoms.clone(),
            del_bonds: other.del_bonds.clone(),
            ..ROMol::default()
        };
        if !quick_copy {
            m.conformers = other
                .conformers
                .iter()
                .filter(|c| conf_id < 0 || c.id as i32 == conf_id)
                .cloned()
                .collect();
            m.sgroups = other.sgroups.clone();
            m.props = other.props.clone();
            m.atom_bookmarks = other.atom_bookmarks.clone();
            m.bond_bookmarks = other.bond_bookmarks.clone();
        } else {
            m.props.set_val(
                crate::rdgeneral::dict::COMPUTED_PROP_NAME,
                Vec::<String>::new(),
            );
        }
        m
    }

    // ------------------------------------------------------------------ atoms

    /// `getNumAtoms()` (explicit atoms only).
    pub fn num_atoms(&self) -> u32 {
        self.atoms.len() as u32
    }

    /// `getNumAtoms(onlyExplicit=false)`: includes implicit Hs.
    pub fn num_atoms_with_hs(&self) -> Result<u32> {
        let mut res = self.num_atoms();
        for i in 0..self.num_atoms() {
            res += self.total_num_hs(i, false)?;
        }
        Ok(res)
    }

    pub fn num_heavy_atoms(&self) -> u32 {
        self.atoms.iter().filter(|a| a.atomic_num > 1).count() as u32
    }

    /// Unchecked access; panics on a bad index (internal use).
    pub fn atom(&self, idx: u32) -> &Atom {
        &self.atoms[idx as usize]
    }
    pub fn atom_mut(&mut self, idx: u32) -> &mut Atom {
        &mut self.atoms[idx as usize]
    }

    /// `getAtomWithIdx` with upstream's range check.
    pub fn get_atom(&self, idx: u32) -> Result<&Atom> {
        self.check_atom(idx)?;
        Ok(&self.atoms[idx as usize])
    }
    pub fn get_atom_mut(&mut self, idx: u32) -> Result<&mut Atom> {
        self.check_atom(idx)?;
        Ok(&mut self.atoms[idx as usize])
    }

    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    fn check_atom(&self, idx: u32) -> Result<()> {
        if idx >= self.num_atoms() {
            return Err(RDKitError::range_at(
                "idx",
                idx,
                self.num_atoms(),
                "Code/GraphMol/ROMol.cpp",
                202,
            ));
        }
        Ok(())
    }

    /// `getAtomDegree`
    pub fn atom_degree(&self, idx: u32) -> u32 {
        self.adj[idx as usize].len() as u32
    }

    /// Neighbor atom indices in upstream order.
    pub fn atom_neighbors(&self, idx: u32) -> impl Iterator<Item = u32> + '_ {
        self.adj[idx as usize].iter().map(|&(n, _)| n)
    }

    /// Bond indices of an atom in upstream order.
    pub fn atom_bond_indices(&self, idx: u32) -> impl Iterator<Item = u32> + '_ {
        self.adj[idx as usize].iter().map(|&(_, b)| b)
    }

    /// Bonds of an atom in upstream order (`atomBonds`).
    pub fn atom_bonds(&self, idx: u32) -> impl Iterator<Item = &Bond> + '_ {
        self.adj[idx as usize]
            .iter()
            .map(move |&(_, b)| &self.bonds[b as usize])
    }

    // ------------------------------------------------------------------ bonds

    /// `getNumBonds()` (heavy, i.e. explicit, bonds only).
    pub fn num_bonds(&self) -> u32 {
        self.bonds.len() as u32
    }

    /// `getNumBonds(onlyHeavy=false)`: includes bonds to implicit Hs.
    pub fn num_bonds_with_hs(&self) -> Result<u32> {
        let mut res = self.num_bonds();
        for i in 0..self.num_atoms() {
            res += self.total_num_hs(i, false)?;
        }
        Ok(res)
    }

    pub fn bond(&self, idx: u32) -> &Bond {
        &self.bonds[idx as usize]
    }
    pub fn bond_mut(&mut self, idx: u32) -> &mut Bond {
        &mut self.bonds[idx as usize]
    }
    pub fn bonds(&self) -> &[Bond] {
        &self.bonds
    }

    /// `getBondWithIdx` with upstream's range check.
    pub fn get_bond(&self, idx: u32) -> Result<&Bond> {
        if idx >= self.num_bonds() {
            return Err(RDKitError::range_at(
                "idx",
                idx,
                self.num_bonds(),
                "Code/GraphMol/ROMol.cpp",
                322,
            ));
        }
        Ok(&self.bonds[idx as usize])
    }
    pub fn get_bond_mut(&mut self, idx: u32) -> Result<&mut Bond> {
        if idx >= self.num_bonds() {
            return Err(RDKitError::range_at(
                "idx",
                idx,
                self.num_bonds(),
                "Code/GraphMol/ROMol.cpp",
                322,
            ));
        }
        Ok(&mut self.bonds[idx as usize])
    }

    /// Index of the bond between two atoms (no range check).
    pub fn bond_idx_between(&self, a1: u32, a2: u32) -> Option<u32> {
        self.adj
            .get(a1 as usize)?
            .iter()
            .find(|&&(n, _)| n == a2)
            .map(|&(_, b)| b)
    }

    /// `getBondBetweenAtoms` (range-checked; `None` when not bonded).
    pub fn get_bond_between_atoms(&self, a1: u32, a2: u32) -> Result<Option<&Bond>> {
        self.check_atom(a1)?;
        self.check_atom(a2)?;
        Ok(self
            .bond_idx_between(a1, a2)
            .map(|b| &self.bonds[b as usize]))
    }

    /// `Bond::setBeginAtomIdx` for an owned bond (range-checked).
    pub fn set_bond_begin_atom_idx(&mut self, bond: u32, what: u32) -> Result<()> {
        if what >= self.num_atoms() {
            return Err(RDKitError::range("what", what, self.num_atoms()));
        }
        self.bonds[bond as usize].begin = what;
        Ok(())
    }

    /// `Bond::setEndAtomIdx` for an owned bond (range-checked).
    pub fn set_bond_end_atom_idx(&mut self, bond: u32, what: u32) -> Result<()> {
        if what >= self.num_atoms() {
            return Err(RDKitError::range("what", what, self.num_atoms()));
        }
        self.bonds[bond as usize].end = what;
        Ok(())
    }

    /// `Bond::setStereoAtoms`
    pub fn set_bond_stereo_atoms(&mut self, bond: u32, bgn: u32, end: u32) -> Result<()> {
        let b = &self.bonds[bond as usize];
        if self.bond_idx_between(b.begin, bgn).is_none() {
            return Err(RDKitError::precondition(
                "bgnIdx not connected to begin atom of bond",
            ));
        }
        if self.bond_idx_between(b.end, end).is_none() {
            return Err(RDKitError::precondition(
                "endIdx not connected to end atom of bond",
            ));
        }
        self.bonds[bond as usize].stereo_atoms = vec![bgn as i32, end as i32];
        Ok(())
    }

    // ------------------------------------------------------------------ editing

    /// `ROMol::addAtom(atom, updateLabel, takeOwnership)`; returns the new index.
    pub fn add_atom(&mut self, mut atom: Atom, update_label: bool) -> u32 {
        let which = self.num_atoms();
        atom.has_owner = true;
        atom.idx = which;
        self.atoms.push(atom);
        self.adj.push(Vec::new());
        if update_label {
            self.replace_atom_bookmark(which, CI_RIGHTMOST_ATOM);
        }
        for conf in &mut self.conformers {
            let _ = conf.set_atom_pos(which, Point3D::default());
        }
        which
    }

    /// `RWMol::addAtom(updateLabel)`: add a dummy atom.
    pub fn add_new_atom(&mut self, update_label: bool) -> u32 {
        self.add_atom(Atom::new(0), update_label)
    }

    /// `ROMol::addBond(Bond*)`: inserts a copy of `bond` using its begin/end
    /// indices; returns the new number of bonds.
    pub fn add_bond_object(&mut self, mut bond: Bond) -> Result<u32> {
        self.check_atom(bond.begin)?;
        self.check_atom(bond.end)?;
        if bond.begin == bond.end {
            return Err(RDKitError::precondition("attempt to add self-bond"));
        }
        if self.bond_idx_between(bond.begin, bond.end).is_some() {
            return Err(RDKitError::precondition("bond already exists"));
        }
        let idx = self.num_bonds();
        bond.has_owner = true;
        bond.idx = idx;
        self.adj[bond.begin as usize].push((bond.end, idx));
        self.adj[bond.end as usize].push((bond.begin, idx));
        self.bonds.push(bond);
        Ok(self.num_bonds())
    }

    /// `RWMol::addBond(idx1, idx2, type)`; returns the new number of bonds.
    pub fn add_bond(&mut self, a1: u32, a2: u32, bond_type: BondType) -> Result<u32> {
        self.check_atom(a1)?;
        self.check_atom(a2)?;
        if a1 == a2 {
            return Err(RDKitError::precondition("attempt to add self-bond"));
        }
        if self.bond_idx_between(a1, a2).is_some() {
            return Err(RDKitError::precondition("bond already exists"));
        }
        let mut b = Bond::new(bond_type);
        if bond_type == BondType::Aromatic {
            b.is_aromatic = true;
            self.atoms[a1 as usize].is_aromatic = true;
            self.atoms[a2 as usize].is_aromatic = true;
        }
        b.begin = a1;
        b.end = a2;
        let n = self.add_bond_object(b)?;
        self.atoms[a1 as usize].clear_property_cache();
        self.atoms[a2 as usize].clear_property_cache();
        if let Some(del_atoms) = &self.del_atoms {
            let deleted = |i: u32| del_atoms.get(i as usize).copied().unwrap_or(false);
            if deleted(a1) || deleted(a2) {
                let del_bonds = self.del_bonds.get_or_insert_with(Vec::new);
                if del_bonds.len() < n as usize {
                    del_bonds.resize(n as usize, false);
                }
                del_bonds[n as usize - 1] = true;
            }
        }
        Ok(n)
    }

    /// Remove bond `idx` from the graph, renumbering higher bond indices.
    fn erase_bond(&mut self, idx: u32) {
        let (b, e) = (self.bonds[idx as usize].begin, self.bonds[idx as usize].end);
        self.adj[b as usize].retain(|&(n, _)| n != e);
        self.adj[e as usize].retain(|&(n, _)| n != b);
        self.bonds.remove(idx as usize);
        for list in &mut self.adj {
            for entry in list.iter_mut() {
                if entry.1 > idx {
                    entry.1 -= 1;
                }
            }
        }
        for (i, bond) in self.bonds.iter_mut().enumerate() {
            bond.idx = i as u32;
        }
        renumber_after_removal(&mut self.bond_bookmarks, idx);
        for g in &mut self.stereo_groups {
            for b in &mut g.bonds {
                if *b > idx {
                    *b -= 1;
                }
            }
        }
    }

    /// Clear stereo on bonds that used the bond `idx` to define stereo atoms.
    fn clear_stereo_referencing_bond(&mut self, idx: u32) {
        let (bgn, end) = (self.bonds[idx as usize].begin, self.bonds[idx as usize].end);
        for (a0, a1) in [(bgn, end), (end, bgn)] {
            let others: Vec<u32> = self.atom_bond_indices(a0).filter(|&o| o != idx).collect();
            for o in others {
                let ob = &mut self.bonds[o as usize];
                if ob.stereo_atoms.contains(&(a1 as i32)) {
                    if matches!(ob.stereo, BondStereo::StereoCis | BondStereo::StereoTrans) {
                        ob.stereo = BondStereo::StereoNone;
                    }
                    ob.stereo_atoms.clear();
                }
            }
        }
    }

    fn clear_bond_bookmarks_for(&mut self, idx: u32) {
        self.bond_bookmarks.retain(|_, list| {
            if let Some(p) = list.iter().position(|&b| b == idx) {
                list.remove(p);
            }
            !list.is_empty()
        });
    }

    fn clear_atom_bookmarks_for(&mut self, idx: u32) {
        self.atom_bookmarks.retain(|_, list| {
            if let Some(p) = list.iter().position(|&a| a == idx) {
                list.remove(p);
            }
            !list.is_empty()
        });
    }

    /// `RWMol::removeBond(aid1, aid2)`
    pub fn remove_bond(&mut self, a1: u32, a2: u32) -> Result<()> {
        self.check_atom(a1)?;
        self.check_atom(a2)?;
        let Some(idx) = self.bond_idx_between(a1, a2) else {
            return Ok(());
        };
        if let Some(del) = &mut self.del_bonds {
            if del.len() < self.bonds.len() {
                del.resize(self.bonds.len(), false);
            }
            del[idx as usize] = true;
            return Ok(());
        }
        self.clear_bond_bookmarks_for(idx);
        self.clear_stereo_referencing_bond(idx);
        self.ring_info.reset(true);
        remove_substance_groups_referencing_bond(self, idx)?;
        stereo_group::remove_bond_from_groups(idx, &mut self.stereo_groups);
        self.erase_bond(idx);
        Ok(())
    }

    /// Update `_MolFileBondEndPts` (V3000 multi-attachment bonds) after atom
    /// `idx` is removed.
    fn fix_bond_end_pts(&mut self, idx: u32) {
        for bond in &mut self.bonds {
            let Some(RDValue::String(sprop)) = bond.props.get(cp::_MolFileBondEndPts) else {
                continue;
            };
            if !(sprop.starts_with('(') && sprop.ends_with(')')) {
                continue;
            }
            let inner = &sprop[1..sprop.len() - 1];
            let mut toks = inner.split(' ').filter(|t| !t.is_empty());
            let mut num_ats = toks
                .next()
                .and_then(|t| t.parse::<f64>().ok())
                .unwrap_or(0.0) as u32;
            let mut oats: Vec<u32> = toks
                .filter_map(|t| t.parse::<f64>().ok())
                .map(|v| v as u32)
                .collect();
            if let Some(p) = oats.iter().position(|&a| a == idx + 1) {
                oats.remove(p);
                num_ats = num_ats.saturating_sub(1);
            }
            if num_ats == 0 {
                bond.clear_prop(cp::_MolFileBondEndPts);
                bond.clear_prop(cp::_MolFileBondAttach);
            } else {
                let mut s = format!("({num_ats} ");
                for i in &mut oats {
                    if *i > idx + 1 {
                        *i -= 1;
                    }
                    s.push_str(&format!("{i} "));
                }
                s.pop();
                s.push(')');
                bond.props.set_val(cp::_MolFileBondEndPts, s);
            }
        }
    }

    /// `RWMol::removeAtom(idx)` (clears computed props).
    pub fn remove_atom(&mut self, idx: u32) -> Result<()> {
        self.remove_atom_ext(idx, true)
    }

    /// `RWMol::removeAtom(atom, clearProps)`
    pub fn remove_atom_ext(&mut self, idx: u32, clear_props: bool) -> Result<()> {
        self.check_atom(idx)?;
        let nbrs: Vec<u32> = self.atom_neighbors(idx).collect();
        for nbr in nbrs {
            self.remove_bond(idx, nbr)?;
        }
        if let Some(del) = &mut self.del_atoms {
            if del.len() < self.atoms.len() {
                del.resize(self.atoms.len(), false);
            }
            del[idx as usize] = true;
            return Ok(());
        }
        self.clear_atom_bookmarks_for(idx);
        for conf in &mut self.conformers {
            if (idx as usize) < conf.positions.len() {
                conf.positions.remove(idx as usize);
            }
        }
        self.fix_bond_end_pts(idx);
        for bond in &mut self.bonds {
            if bond.begin > idx {
                bond.begin -= 1;
            }
            if bond.end > idx {
                bond.end -= 1;
            }
            let mut clear = false;
            for sa in &mut bond.stereo_atoms {
                if *sa == idx as i32 {
                    clear = true;
                    break;
                } else if *sa > idx as i32 {
                    *sa -= 1;
                }
            }
            if clear {
                bond.stereo_atoms.clear();
            }
        }
        remove_substance_groups_referencing_atom(self, idx)?;
        stereo_group::remove_atom_from_groups(idx, &mut self.stereo_groups);
        if clear_props {
            self.clear_computed_props(true);
        }
        self.remove_vertex(idx);
        Ok(())
    }

    /// `boost::remove_vertex`: drop atom `idx` (no incident bonds left) and
    /// renumber everything referring to higher atom indices.
    fn remove_vertex(&mut self, idx: u32) {
        self.atoms.remove(idx as usize);
        self.adj.remove(idx as usize);
        for (i, a) in self.atoms.iter_mut().enumerate() {
            a.idx = i as u32;
        }
        for list in &mut self.adj {
            for e in list.iter_mut() {
                if e.0 > idx {
                    e.0 -= 1;
                }
            }
        }
        renumber_after_removal(&mut self.atom_bookmarks, idx);
        for g in &mut self.stereo_groups {
            for a in &mut g.atoms {
                if *a > idx {
                    *a -= 1;
                }
            }
        }
    }

    /// `RWMol::replaceAtom(idx, atom, updateLabel, preserveProps)`
    pub fn replace_atom(&mut self, idx: u32, atom: &Atom, preserve_props: bool) -> Result<()> {
        self.check_atom(idx)?;
        let mut new = atom.copy_detached();
        new.has_owner = true;
        new.idx = idx;
        if preserve_props {
            let old = self.atoms[idx as usize].props.clone();
            new.update_props(&old, true);
        }
        self.atoms[idx as usize] = new;
        Ok(())
    }

    /// `RWMol::replaceBond(idx, bond, preserveProps, keepSGroups)`
    pub fn replace_bond(
        &mut self,
        idx: u32,
        bond: &Bond,
        preserve_props: bool,
        keep_sgroups: bool,
    ) -> Result<()> {
        if idx >= self.num_bonds() {
            return Err(RDKitError::range("idx", idx, self.num_bonds()));
        }
        let old = &self.bonds[idx as usize];
        let mut new = bond.copy_detached();
        new.has_owner = true;
        new.idx = idx;
        new.begin = old.begin;
        new.end = old.end;
        let order_difference = new.bond_type.as_double()? - old.bond_type.as_double()?;
        if order_difference > 0.0 {
            for a in [new.begin, new.end] {
                let atom = &mut self.atoms[a as usize];
                if atom.num_explicit_hs > 0 {
                    let new_hs = (f64::from(atom.num_explicit_hs) - order_difference) as i32;
                    atom.num_explicit_hs = new_hs.max(0) as u8;
                }
            }
        }
        if preserve_props {
            let old_props = self.bonds[idx as usize].props.clone();
            new.update_props(&old_props, true);
        }
        self.bonds[idx as usize] = new;
        if !keep_sgroups {
            remove_substance_groups_referencing_bond(self, idx)?;
        }
        Ok(())
    }

    /// `RWMol::insertMol`
    pub fn insert_mol(&mut self, other: &ROMol) -> Result<()> {
        let orig_atoms = self.num_atoms();
        let orig_bonds = self.num_bonds();
        for oatom in &other.atoms {
            let mut new = oatom.copy_detached();
            if let Some(RDValue::VecInt(n_atoms)) = new.props.get(cp::_ringStereoAtoms).cloned() {
                let shifted: Vec<i32> = n_atoms
                    .into_iter()
                    .map(|v| {
                        if v < 0 {
                            v - orig_atoms as i32
                        } else {
                            v + orig_atoms as i32
                        }
                    })
                    .collect();
                new.set_prop(cp::_ringStereoAtoms, shifted, true)?;
            }
            if let Some(val) = oatom
                .props
                .get(cp::_ringStereoOtherAtom)
                .and_then(RDValue::as_uint)
            {
                new.set_prop(cp::_ringStereoOtherAtom, val + orig_atoms, true)?;
            }
            self.add_atom(new, false);
        }
        for obond in &other.bonds {
            let mut b = obond.copy_detached();
            b.begin += orig_atoms;
            b.end += orig_atoms;
            for v in &mut b.stereo_atoms {
                *v += orig_atoms as i32;
            }
            self.add_bond_object(b)?;
        }
        if !other.conformers.is_empty() && self.conformers.is_empty() {
            let n = self.num_atoms() as usize;
            for oconf in &other.conformers {
                let mut nconf = Conformer::new(n);
                nconf.is_3d = oconf.is_3d;
                nconf.id = oconf.id;
                for (i, p) in oconf.positions.iter().enumerate() {
                    nconf.set_atom_pos(i as u32 + orig_atoms, *p)?;
                }
                self.add_conformer(nconf, false)?;
            }
        } else if !self.conformers.is_empty() && other.conformers.len() == self.conformers.len() {
            for (conf, oconf) in self.conformers.iter_mut().zip(&other.conformers) {
                for (i, p) in oconf.positions.iter().enumerate() {
                    conf.set_atom_pos(i as u32 + orig_atoms, *p)?;
                }
            }
        }
        self.insert_stereo_groups(other, orig_atoms, orig_bonds)?;
        insert_substance_groups(self, other, orig_atoms, orig_bonds)
    }

    fn insert_stereo_groups(
        &mut self,
        other: &ROMol,
        orig_atoms: u32,
        orig_bonds: u32,
    ) -> Result<()> {
        if other.stereo_groups.is_empty() {
            return Ok(());
        }
        let mut abs_atoms = Vec::new();
        let mut abs_bonds = Vec::new();
        let mut new_groups = Vec::new();
        for sg in &self.stereo_groups {
            if sg.group_type == StereoGroupType::Absolute {
                abs_atoms.extend_from_slice(&sg.atoms);
                abs_bonds.extend_from_slice(&sg.bonds);
            } else {
                new_groups.push(sg.clone());
            }
        }
        for sg in &other.stereo_groups {
            let atoms: Vec<u32> = sg.atoms.iter().map(|a| a + orig_atoms).collect();
            let bonds: Vec<u32> = sg.bonds.iter().map(|b| b + orig_bonds).collect();
            if sg.group_type == StereoGroupType::Absolute {
                abs_atoms.extend(atoms);
                abs_bonds.extend(bonds);
            } else {
                let mut g = StereoGroup::new(sg.group_type, atoms, bonds, sg.read_id)?;
                g.write_id = 0;
                new_groups.push(g);
            }
        }
        if !abs_atoms.is_empty() || !abs_bonds.is_empty() {
            new_groups.push(StereoGroup::new(
                StereoGroupType::Absolute,
                abs_atoms,
                abs_bonds,
                0,
            )?);
        }
        self.set_stereo_groups(new_groups)
    }

    // ------------------------------------------------------------------ batch edits

    /// `RWMol::beginBatchEdit`
    pub fn begin_batch_edit(&mut self) -> Result<()> {
        if self.del_atoms.is_some() || self.del_bonds.is_some() {
            return Err(RDKitError::Value(
                "Attempt to re-enter batchEdit mode".into(),
            ));
        }
        self.del_atoms = Some(vec![false; self.atoms.len()]);
        self.del_bonds = Some(vec![false; self.bonds.len()]);
        Ok(())
    }

    /// `RWMol::rollbackBatchEdit`
    pub fn rollback_batch_edit(&mut self) {
        self.del_atoms = None;
        self.del_bonds = None;
    }

    /// `RWMol::commitBatchEdit`
    pub fn commit_batch_edit(&mut self) -> Result<()> {
        let (Some(del_atoms), Some(del_bonds)) = (self.del_atoms.take(), self.del_bonds.take())
        else {
            self.del_atoms = None;
            self.del_bonds = None;
            return Ok(());
        };
        if !del_atoms.iter().any(|&d| d) && !del_bonds.iter().any(|&d| d) {
            return Ok(());
        }
        // remove bonds from the highest index down
        for i in (0..del_bonds.len()).rev() {
            if !del_bonds[i] || i >= self.bonds.len() {
                continue;
            }
            let idx = i as u32;
            self.clear_bond_bookmarks_for(idx);
            self.clear_stereo_referencing_bond(idx);
            remove_substance_groups_referencing_bond(self, idx)?;
            stereo_group::remove_bond_from_groups(idx, &mut self.stereo_groups);
            self.erase_bond(idx);
        }
        // then atoms, again from the highest index down
        let n_before = self.atoms.len();
        let mut alive = vec![true; n_before];
        for i in (0..del_atoms.len()).rev() {
            if !del_atoms[i] || i >= self.atoms.len() {
                continue;
            }
            let idx = i as u32;
            // bonds to deleted atoms were removed above (removeAtom marks them)
            let nbrs: Vec<u32> = self.atom_neighbors(idx).collect();
            for n in nbrs {
                let b = self.bond_idx_between(idx, n).expect("bonded");
                self.erase_bond(b);
            }
            self.clear_atom_bookmarks_for(idx);
            self.fix_bond_end_pts(idx);
            remove_substance_groups_referencing_atom(self, idx)?;
            stereo_group::remove_atom_from_groups(idx, &mut self.stereo_groups);
            alive[i] = false;
            // stereo atoms referring to the deleted atom are dropped below
            self.remove_vertex_keep_bond_refs(idx);
        }
        // map old -> new atom indices
        let mut new_idx = vec![u32::MAX; n_before];
        let mut next = 0u32;
        for (i, a) in alive.iter().enumerate() {
            if *a {
                new_idx[i] = next;
                next += 1;
            }
        }
        for bond in &mut self.bonds {
            bond.begin = new_idx[bond.begin as usize];
            bond.end = new_idx[bond.end as usize];
            if !bond.stereo_atoms.is_empty() {
                bond.stereo_atoms = bond
                    .stereo_atoms
                    .iter()
                    .filter_map(|&sa| {
                        new_idx
                            .get(sa as usize)
                            .filter(|&&v| v != u32::MAX)
                            .map(|&v| v as i32)
                    })
                    .collect();
            }
        }
        for conf in &mut self.conformers {
            let pos: Vec<Point3D> = conf
                .positions
                .iter()
                .enumerate()
                .filter(|(i, _)| alive.get(*i).copied().unwrap_or(false))
                .map(|(_, p)| *p)
                .collect();
            if pos.len() != self.atoms.len() {
                return Err(RDKitError::invariant("Lost coordinates!"));
            }
            conf.positions = pos;
        }
        self.ring_info.reset(true);
        self.clear_computed_props(true);
        Ok(())
    }

    /// Vertex removal used by batch edits: bond begin/end indices are left in
    /// the *old* numbering and remapped afterwards (as upstream does).
    fn remove_vertex_keep_bond_refs(&mut self, idx: u32) {
        self.atoms.remove(idx as usize);
        self.adj.remove(idx as usize);
        for (i, a) in self.atoms.iter_mut().enumerate() {
            a.idx = i as u32;
        }
        for list in &mut self.adj {
            for e in list.iter_mut() {
                if e.0 > idx {
                    e.0 -= 1;
                }
            }
        }
        renumber_after_removal(&mut self.atom_bookmarks, idx);
        for g in &mut self.stereo_groups {
            for a in &mut g.atoms {
                if *a > idx {
                    *a -= 1;
                }
            }
        }
    }

    // ------------------------------------------------------------------ bookmarks

    pub fn set_atom_bookmark(&mut self, atom: u32, mark: i32) {
        self.atom_bookmarks.entry(mark).or_default().push(atom);
    }
    pub fn replace_atom_bookmark(&mut self, atom: u32, mark: i32) {
        self.atom_bookmarks.insert(mark, vec![atom]);
    }
    pub fn has_atom_bookmark(&self, mark: i32) -> bool {
        self.atom_bookmarks.contains_key(&mark)
    }
    pub fn get_atom_with_bookmark(&self, mark: i32) -> Result<u32> {
        self.atom_bookmarks
            .get(&mark)
            .and_then(|l| l.first().copied())
            .ok_or_else(|| RDKitError::precondition("atom bookmark not found"))
    }
    pub fn get_all_atoms_with_bookmark(&self, mark: i32) -> Result<&[u32]> {
        self.atom_bookmarks
            .get(&mark)
            .map(Vec::as_slice)
            .ok_or_else(|| RDKitError::precondition("atom bookmark not found"))
    }
    pub fn clear_atom_bookmark(&mut self, mark: i32) {
        self.atom_bookmarks.remove(&mark);
    }
    pub fn clear_atom_bookmark_for(&mut self, mark: i32, atom: u32) {
        if let Some(list) = self.atom_bookmarks.get_mut(&mark) {
            if let Some(p) = list.iter().position(|&a| a == atom) {
                list.remove(p);
            }
            if list.is_empty() {
                self.atom_bookmarks.remove(&mark);
            }
        }
    }
    pub fn clear_all_atom_bookmarks(&mut self) {
        self.atom_bookmarks.clear();
    }
    pub fn set_bond_bookmark(&mut self, bond: u32, mark: i32) {
        self.bond_bookmarks.entry(mark).or_default().push(bond);
    }
    pub fn has_bond_bookmark(&self, mark: i32) -> bool {
        self.bond_bookmarks.contains_key(&mark)
    }
    pub fn get_bond_with_bookmark(&self, mark: i32) -> Result<u32> {
        self.bond_bookmarks
            .get(&mark)
            .and_then(|l| l.first().copied())
            .ok_or_else(|| RDKitError::precondition("bond bookmark not found"))
    }
    /// `getUniqueBondWithBookmark`: the first bond with the bookmark.
    pub fn get_unique_bond_with_bookmark(&self, mark: i32) -> Result<u32> {
        self.bond_bookmarks
            .get(&mark)
            .and_then(|l| l.first().copied())
            .ok_or_else(|| RDKitError::precondition("bookmark not found"))
    }
    pub fn get_all_bonds_with_bookmark(&self, mark: i32) -> Result<&[u32]> {
        self.bond_bookmarks
            .get(&mark)
            .map(Vec::as_slice)
            .ok_or_else(|| RDKitError::precondition("bond bookmark not found"))
    }
    pub fn clear_bond_bookmark(&mut self, mark: i32) {
        self.bond_bookmarks.remove(&mark);
    }
    pub fn clear_all_bond_bookmarks(&mut self) {
        self.bond_bookmarks.clear();
    }

    // ------------------------------------------------------------------ conformers

    pub fn num_conformers(&self) -> u32 {
        self.conformers.len() as u32
    }
    pub fn conformers(&self) -> &[Conformer] {
        &self.conformers
    }
    pub fn conformers_mut(&mut self) -> &mut [Conformer] {
        &mut self.conformers
    }

    fn conformer_position(&self, id: i32) -> Result<usize> {
        if self.conformers.is_empty() {
            return Err(RDKitError::Conformer(
                "No conformations available on the molecule".into(),
            ));
        }
        if id < 0 {
            return Ok(0);
        }
        self.conformers
            .iter()
            .position(|c| c.id == id as u32)
            .ok_or_else(|| RDKitError::Conformer(format!("Can't find conformation with ID: {id}")))
    }

    /// `getConformer(id)`; `id < 0` returns the first one.
    pub fn conformer(&self, id: i32) -> Result<&Conformer> {
        let p = self.conformer_position(id)?;
        Ok(&self.conformers[p])
    }
    pub fn conformer_mut(&mut self, id: i32) -> Result<&mut Conformer> {
        let p = self.conformer_position(id)?;
        Ok(&mut self.conformers[p])
    }

    /// `addConformer(conf, assignId)`; returns the conformer id.
    pub fn add_conformer(&mut self, mut conf: Conformer, assign_id: bool) -> Result<u32> {
        if conf.num_atoms() != self.num_atoms() {
            return Err(RDKitError::precondition("Number of atom mismatch"));
        }
        if assign_id {
            let max_id = self
                .conformers
                .iter()
                .map(|c| c.id as i64)
                .max()
                .unwrap_or(-1);
            conf.id = (max_id + 1) as u32;
        }
        conf.has_owner = true;
        let id = conf.id;
        self.conformers.push(conf);
        Ok(id)
    }

    pub fn remove_conformer(&mut self, id: u32) {
        if let Some(p) = self.conformers.iter().position(|c| c.id == id) {
            self.conformers.remove(p);
        }
    }

    pub fn clear_conformers(&mut self) {
        self.conformers.clear();
    }

    // ------------------------------------------------------------------ misc

    pub fn ring_info(&self) -> &RingInfo {
        &self.ring_info
    }
    pub fn ring_info_mut(&mut self) -> &mut RingInfo {
        &mut self.ring_info
    }

    /// `getSubstanceGroups(mol)`
    pub fn substance_groups(&self) -> &[SubstanceGroup] {
        &self.sgroups
    }

    /// `getSubstanceGroups(mol)` for modification; the groups must keep
    /// referring to existing atoms and bonds.
    pub fn substance_groups_mut(&mut self) -> &mut Vec<SubstanceGroup> {
        &mut self.sgroups
    }

    pub fn stereo_groups(&self) -> &[StereoGroup] {
        &self.stereo_groups
    }

    /// `setStereoGroups`: multiple ABSOLUTE groups are merged into one.
    pub fn set_stereo_groups(&mut self, groups: Vec<StereoGroup>) -> Result<()> {
        let num_abs = groups
            .iter()
            .filter(|g| g.group_type == StereoGroupType::Absolute)
            .count();
        if num_abs <= 1 {
            self.stereo_groups = groups;
            return Ok(());
        }
        let mut abs_atoms: Vec<u32> = Vec::new();
        let mut abs_bonds: Vec<u32> = Vec::new();
        let mut new_groups = Vec::new();
        for g in groups {
            if g.group_type == StereoGroupType::Absolute {
                // upstream inserts at the front
                abs_atoms.splice(0..0, g.atoms.iter().copied());
                abs_bonds.splice(0..0, g.bonds.iter().copied());
            } else {
                new_groups.push(g);
            }
        }
        new_groups.push(StereoGroup::new(
            StereoGroupType::Absolute,
            abs_atoms,
            abs_bonds,
            0,
        )?);
        self.stereo_groups = new_groups;
        Ok(())
    }

    pub fn has_query(&self) -> bool {
        self.atoms.iter().any(Atom::has_query) || self.bonds.iter().any(Bond::has_query)
    }

    /// `setName` / `getName` (the `_Name` property)
    pub fn set_name(&mut self, name: &str) {
        self.props.set_val(cp::_Name, name);
    }
    pub fn name(&self) -> String {
        self.props
            .get(cp::_Name)
            .map(RDValue::to_rd_string)
            .unwrap_or_default()
    }

    /// `ROMol::clearComputedProps(includeRings)`
    pub fn clear_computed_props(&mut self, include_rings: bool) {
        if include_rings {
            self.ring_info.reset(true);
        }
        RDProps::clear_computed_props(self);
        for a in &mut self.atoms {
            a.clear_computed_props();
        }
        for b in &mut self.bonds {
            b.clear_computed_props();
        }
    }

    // ------------------------------------------------------------------ valence / property cache

    /// `Atom::calcExplicitValence(strict)`
    pub fn calc_explicit_valence(&mut self, idx: u32, strict: bool) -> Result<i32> {
        let v = atom::calculate_explicit_valence(self, idx, strict, false)?;
        self.atoms[idx as usize].explicit_valence = v as i8;
        Ok(v)
    }

    /// `Atom::calcImplicitValence(strict)`
    pub fn calc_implicit_valence(&mut self, idx: u32, strict: bool) -> Result<i32> {
        if self.atoms[idx as usize].explicit_valence == -1 {
            self.calc_explicit_valence(idx, strict)?;
        }
        let v = atom::calculate_implicit_valence(self, idx, strict, false)?;
        self.atoms[idx as usize].implicit_valence = v as i8;
        Ok(v)
    }

    /// `Atom::updatePropertyCache(strict)`
    pub fn update_atom_property_cache(&mut self, idx: u32, strict: bool) -> Result<()> {
        self.calc_explicit_valence(idx, strict)?;
        self.calc_implicit_valence(idx, strict)?;
        Ok(())
    }

    /// `ROMol::updatePropertyCache(strict)`
    pub fn update_property_cache(&mut self, strict: bool) -> Result<()> {
        for i in 0..self.num_atoms() {
            self.update_atom_property_cache(i, strict)?;
        }
        Ok(())
    }

    pub fn needs_update_property_cache(&self) -> bool {
        self.atoms.iter().any(Atom::needs_update_property_cache)
    }

    pub fn clear_property_cache(&mut self) {
        for a in &mut self.atoms {
            a.clear_property_cache();
        }
    }

    /// `Atom::getTotalNumHs(includeNeighbors)`
    pub fn total_num_hs(&self, idx: u32, include_neighbors: bool) -> Result<u32> {
        let a = &self.atoms[idx as usize];
        let mut res = a.num_explicit_hs() + a.num_implicit_hs()?;
        if include_neighbors {
            res += self
                .atom_neighbors(idx)
                .filter(|&n| self.atoms[n as usize].atomic_num == 1)
                .count() as u32;
        }
        Ok(res)
    }

    /// `Atom::getTotalDegree`
    pub fn total_degree(&self, idx: u32) -> Result<u32> {
        Ok(self.total_num_hs(idx, false)? + self.atom_degree(idx))
    }

    /// `Atom::getValence(which)` for an atom of this molecule.
    pub fn atom_valence(&self, idx: u32, which: ValenceType) -> Result<u32> {
        self.atoms[idx as usize].valence(which)
    }

    /// `Atom::getPerturbationOrder(probe)`
    pub fn perturbation_order(&self, idx: u32, probe: &[i32]) -> Result<i32> {
        let bonds: Vec<i32> = self.atom_bond_indices(idx).map(|b| b as i32).collect();
        Ok(crate::rdgeneral::utils::count_swaps_to_interconvert(probe, &bonds)? as i32)
    }
}

/// Shift bookmark indices after the element `removed` disappeared.
fn renumber_after_removal(marks: &mut BTreeMap<i32, Vec<u32>>, removed: u32) {
    for list in marks.values_mut() {
        for v in list.iter_mut() {
            if *v > removed {
                *v -= 1;
            }
        }
    }
}
