//! The molecule the labeler works on (upstream `CIPLabeler/CIPMol.{h,cpp}`).
//!
//! Upstream refers to the caller's molecule; here the labeler works on a copy,
//! so the digraphs can borrow it while the labels are written to the
//! original.

use super::mancude::{FractionalAtomicNum, calc_frac_atom_nums};
use crate::Result;
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::{Bond, BondType};
use crate::graphmol::find_rings::fast_find_rings;
use crate::graphmol::kekulize::kekulize;
use crate::graphmol::romol::ROMol;
use crate::{RDKitError, rd_log_warning};
use std::cell::RefCell;

pub struct CipMol {
    mol: ROMol,
    kekulized_bonds: RefCell<Option<Vec<BondType>>>,
    atomnums: RefCell<Option<Vec<FractionalAtomicNum>>>,
}

impl CipMol {
    pub fn new(mol: &ROMol) -> Result<Self> {
        let mut mol = mol.clone();
        if !mol.ring_info().is_find_fast_or_better() {
            fast_find_rings(&mut mol)?;
        }
        Ok(CipMol {
            mol,
            kekulized_bonds: RefCell::new(None),
            atomnums: RefCell::new(None),
        })
    }

    pub fn mol(&self) -> &ROMol {
        &self.mol
    }

    /// `getFractionalAtomicNum`
    pub fn fractional_atomic_num(&self, atom: u32) -> Result<FractionalAtomicNum> {
        if self.atomnums.borrow().is_none() {
            let nums = calc_frac_atom_nums(self)?;
            *self.atomnums.borrow_mut() = Some(nums);
        }
        Ok(self.atomnums.borrow().as_ref().unwrap()[atom as usize])
    }

    pub fn num_atoms(&self) -> u32 {
        self.mol.num_atoms()
    }

    pub fn num_bonds(&self) -> u32 {
        self.mol.num_bonds()
    }

    pub fn atom(&self, idx: u32) -> &Atom {
        self.mol.atom(idx)
    }

    pub fn bond(&self, idx: u32) -> &Bond {
        self.mol.bond(idx)
    }

    pub fn total_num_hs(&self, atom: u32) -> Result<u32> {
        self.mol.total_num_hs(atom, false)
    }

    /// Indices of the bonds of an atom.
    pub fn bonds_of(&self, atom: u32) -> impl Iterator<Item = u32> + '_ {
        self.mol.atom_bond_indices(atom)
    }

    pub fn neighbors_of(&self, atom: u32) -> impl Iterator<Item = u32> + '_ {
        self.mol.atom_neighbors(atom)
    }

    pub fn other_atom(&self, bond: u32, atom: u32) -> Result<u32> {
        self.mol.bond(bond).other_atom_idx(atom)
    }

    /// `isInRing`
    pub fn is_in_ring(&self, bond: u32) -> Result<bool> {
        Ok(self.mol.ring_info().num_bond_rings(bond)? != 0)
    }

    /// `getBondOrder`: of the kekulized molecule when it can be kekulized.
    pub fn bond_order(&self, bond: u32) -> Result<i32> {
        if self.kekulized_bonds.borrow().is_none() {
            let mut tmp = self.mol.clone();
            // a molecule that cannot be kekulized keeps its aromatic bonds
            let _ = kekulize(&mut tmp, true, false, 100);
            let types = tmp.bonds().iter().map(|b| b.bond_type()).collect();
            *self.kekulized_bonds.borrow_mut() = Some(types);
        }
        let bond_type = self.kekulized_bonds.borrow().as_ref().unwrap()[bond as usize];
        Ok(match bond_type {
            BondType::Zero
            | BondType::Hydrogen
            | BondType::Dative
            | BondType::DativeL
            | BondType::DativeR => 0,
            BondType::Single => 1,
            BondType::Aromatic => {
                rd_log_warning!("non kekulizable aromatic bond being treated as bond order 1\n");
                1
            }
            BondType::Double => 2,
            BondType::Triple => 3,
            BondType::Quadruple => 4,
            BondType::Quintuple => 5,
            BondType::Hextuple => 6,
            _ => {
                return Err(RDKitError::Runtime(
                    "Non integer-order bonds are not allowed.".into(),
                ));
            }
        })
    }
}
