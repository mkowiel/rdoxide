//! The molecule building operations of the SLN parser.
//!
//! Ported from RDKit `Code/GraphMol/SLNParse/SLNParseOps.h` (2026.09.1).
//! Upstream keeps the ring-closure bookmarks and the partial bonds of the
//! molecules being built in their atom/bond bookmarks; here they live in
//! [`SlnMol`].
//!
//! Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;

use super::SlnError;
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::{Bond, BondType};
use crate::graphmol::romol::RWMol;
use crate::rdgeneral::{RDProps, common_properties as cp};

/// A molecule being built.
#[derive(Debug)]
pub struct SlnMol {
    pub mol: RWMol,
    /// The atom that new atoms are attached to (`getActiveAtom`).
    active: u32,
    /// Ring closure/atom ID -> atom
    atom_bookmarks: BTreeMap<i32, u32>,
    /// Ring closures without an opener yet: ID -> bonds (their end atom is
    /// the atom that closes the ring)
    bond_bookmarks: BTreeMap<i32, Vec<Bond>>,
}

impl SlnMol {
    /// The closures that never found their opener.
    pub fn open_bond_bookmarks(&self) -> impl Iterator<Item = i32> + '_ {
        self.bond_bookmarks.keys().copied()
    }

    /// `bookmarkAtomID`: remembers an atom with an ID.
    fn bookmark_atom_id(&mut self, atom_idx: u32) -> Result<(), SlnError> {
        let label = match self.mol.atom(atom_idx).get_prop(cp::_AtomID) {
            Some(v) => match v.as_uint() {
                Some(l) => l as i32,
                None => return Ok(()),
            },
            None => return Ok(()),
        };
        if self.atom_bookmarks.contains_key(&label) {
            return Err(SlnError::Sln(format!(
                "SLN Parser error: Atom ID {label} used a second time."
            )));
        }
        if self.bond_bookmarks.contains_key(&label) {
            return Err(SlnError::Sln(format!(
                "SLN Parser error: Atom ID {label} appears *after* its ring closure."
            )));
        }
        self.atom_bookmarks.insert(label, atom_idx);
        Ok(())
    }

    /// `addBondToMol`: adds a bond, being careful to handle aromaticity
    /// properly.
    fn add_bond(&mut self, bond: Bond) -> Result<(), SlnError> {
        let aromatic = bond.bond_type() == BondType::Aromatic;
        let (b, e) = (bond.begin_atom_idx(), bond.end_atom_idx());
        let n = self.mol.add_bond_object(bond)?;
        if aromatic {
            // SLN doesn't have aromatic atom types, aromaticity is a
            // property of the bonds themselves
            self.mol.bond_mut(n - 1).set_is_aromatic(true);
            self.mol.atom_mut(b).set_is_aromatic(true);
            self.mol.atom_mut(e).set_is_aromatic(true);
        }
        Ok(())
    }

    /// Adds the hydrogens that are set on the atom as atoms, "otherwise
    /// getting the numbering right is just too hard".
    fn add_explicit_hs(&mut self, atom_idx: u32) -> Result<(), SlnError> {
        for _ in 0..self.mol.atom(atom_idx).num_explicit_hs() {
            let h = self.mol.add_atom(Atom::new(1), false);
            self.mol.add_bond(atom_idx, h, BondType::Single)?;
        }
        self.mol.atom_mut(atom_idx).set_num_explicit_hs(0);
        Ok(())
    }

    /// `closeRingBond`: closes an indexed ring using the bond provided; the
    /// bond is formed from the atom with the ring ID to the active atom.
    fn close_ring_bond(
        &mut self,
        ring_idx: i32,
        mut bond: Bond,
        postpone_allowed: bool,
    ) -> Result<(), SlnError> {
        let Some(&opener) = self.atom_bookmarks.get(&ring_idx) else {
            if postpone_allowed {
                // save it for later
                bond.set_end_atom_idx_unchecked(self.active);
                self.bond_bookmarks.entry(ring_idx).or_default().push(bond);
                return Ok(());
            }
            return Err(SlnError::Sln(format!(
                "SLN Parser error: Ring closure {ring_idx} does not have a corresponding opener."
            )));
        };
        bond.set_begin_atom_idx_unchecked(opener);
        bond.set_end_atom_idx_unchecked(self.active);
        self.add_bond(bond)
    }
}

pub type MolList = Vec<Option<SlnMol>>;

fn mol_at(list: &mut MolList, idx: usize) -> &mut SlnMol {
    list[idx].as_mut().expect("null molecule")
}

/// `startMol`: initializes a molecule; returns its index in the list.
pub fn start_mol(
    list: &mut MolList,
    first_atom: Atom,
    doing_query: bool,
) -> Result<usize, SlnError> {
    let mut m = SlnMol {
        mol: RWMol::new(),
        active: 0,
        atom_bookmarks: BTreeMap::new(),
        bond_bookmarks: BTreeMap::new(),
    };
    let idx = m.mol.add_atom(first_atom, true);
    m.active = idx;
    m.bookmark_atom_id(idx)?;
    if !doing_query {
        m.add_explicit_hs(idx)?;
    }
    list.push(Some(m));
    Ok(list.len() - 1)
}

/// `addAtomToMol`: adds an atom bound to the active atom with the bond.
pub fn add_atom_to_mol(
    list: &mut MolList,
    idx: usize,
    atom: Atom,
    mut bond: Bond,
    doing_query: bool,
) -> Result<(), SlnError> {
    let mp = mol_at(list, idx);
    let atom_idx1 = mp.active;
    let atom_idx2 = mp.mol.add_atom(atom, true);
    mp.active = atom_idx2;
    mp.bookmark_atom_id(atom_idx2)?;
    bond.set_begin_atom_idx_unchecked(atom_idx1);
    bond.set_end_atom_idx_unchecked(atom_idx2);
    mp.add_bond(bond)?;
    if !doing_query {
        mp.add_explicit_hs(atom_idx2)?;
    }
    Ok(())
}

/// `addAtomToMol` with a single bond
pub fn add_atom_to_mol_single(
    list: &mut MolList,
    idx: usize,
    atom: Atom,
    doing_query: bool,
) -> Result<(), SlnError> {
    add_atom_to_mol(list, idx, atom, Bond::new(BondType::Single), doing_query)
}

/// `closeRingBond`
pub fn close_ring_bond(
    list: &mut MolList,
    mol_idx: usize,
    ring_idx: i32,
    bond: Bond,
) -> Result<(), SlnError> {
    mol_at(list, mol_idx).close_ring_bond(ring_idx, bond, true)
}

/// `addBranchToMol`: adds the branch to the molecule, connected with the
/// bond (the bond is dropped for ionic bonds).
pub fn add_branch_to_mol(
    list: &mut MolList,
    mol_idx: usize,
    branch_idx: usize,
    bond: Bond,
) -> Result<usize, SlnError> {
    let branch = list[branch_idx].take().expect("null branch");
    let mp = mol_at(list, mol_idx);
    let active_atom_idx = mp.active;
    let n_orig_atoms = mp.mol.num_atoms();

    // add the fragment's atoms and bonds to the molecule
    mp.mol.insert_mol(&branch.mol)?;

    // copy in the atom bookmarks of the branch
    for (&label, &atom) in &branch.atom_bookmarks {
        if label < 0 {
            continue;
        }
        if mp.atom_bookmarks.contains_key(&label) {
            return Err(SlnError::Sln(format!(
                "SLN Parser error: Atom ID {label} used a second time."
            )));
        } else if mp.bond_bookmarks.contains_key(&label) {
            return Err(SlnError::Sln(format!(
                "SLN Parser error: Atom ID {label} appears *after* its ring closure."
            )));
        }
        mp.atom_bookmarks.insert(label, atom + n_orig_atoms);
    }

    // close the rings of the bond bookmarks of the branch
    for (label, bonds) in branch.bond_bookmarks {
        for mut tgt in bonds {
            if label > 0 && mp.atom_bookmarks.contains_key(&label) {
                let tmp = mp.active;
                mp.active = tgt.end_atom_idx() + n_orig_atoms;
                mp.close_ring_bond(label, tgt, false)?;
                mp.active = tmp;
            } else {
                // no partner found yet, copy into this mol
                tgt.set_end_atom_idx_unchecked(tgt.end_atom_idx() + n_orig_atoms);
                mp.bond_bookmarks.entry(label).or_default().push(tgt);
            }
        }
    }

    // set the connecting bond
    if bond.bond_type() != BondType::Ionic {
        let mut bond = bond;
        bond.set_begin_atom_idx_unchecked(active_atom_idx);
        bond.set_end_atom_idx_unchecked(n_orig_atoms);
        mp.add_bond(bond)?;
    }

    if list.len() == branch_idx + 1 {
        list.pop();
    }
    Ok(mol_idx)
}

/// `addFragToMol`: adds the atoms and bonds of a fragment to the molecule
/// without a bond between them.
pub fn add_frag_to_mol(
    list: &mut MolList,
    mol_idx: usize,
    frag_idx: usize,
) -> Result<usize, SlnError> {
    add_branch_to_mol(list, mol_idx, frag_idx, Bond::new(BondType::Ionic))
}
