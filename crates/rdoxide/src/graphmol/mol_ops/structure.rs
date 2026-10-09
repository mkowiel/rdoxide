//! Structure-level molecule operations of `MolOps`: `renumberAtoms`
//! (`Renumber.cpp`), `getFormalCharge` and `removeStereochemistry`.
//!
//! Copyright (C) 2013-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::ring_info::FindRingType;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::stereo_group::StereoGroup;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// Where each original atom lands: `rev[old] = new` for the order
/// `new_order[new] = old`.
fn reverse_order(new_order: &[u32], num_atoms: u32) -> Result<Vec<u32>> {
    let mut rev = vec![num_atoms; num_atoms as usize];
    for (new_idx, &old_idx) in new_order.iter().enumerate() {
        if old_idx >= num_atoms {
            return Err(RDKitError::Value("idx value exceeds numAtoms".into()));
        }
        if rev[old_idx as usize] != num_atoms {
            return Err(RDKitError::Value(
                "newOrder contains duplicate indices".into(),
            ));
        }
        rev[old_idx as usize] = new_idx as u32;
    }
    Ok(rev)
}

/// `MolOps::renumberAtoms(mol, newOrder)`: a copy of `mol` in which atom `i`
/// is the atom `new_order[i]` of the original. Molecule-level properties are
/// not carried over, as upstream.
pub fn renumber_atoms(mol: &ROMol, new_order: &[u32]) -> Result<RWMol> {
    let n_atoms = mol.num_atoms();
    if new_order.len() != n_atoms as usize {
        return Err(RDKitError::precondition("bad newOrder size"));
    }
    let rev = reverse_order(new_order, n_atoms)?;

    let mut res = RWMol::default();
    for &old_idx in new_order {
        let mut atom = mol.atom(old_idx).copy_detached();
        // properties that refer to atom numbers
        if let Some(RDValue::VecInt(atoms)) = atom.props.get(cp::_ringStereoAtoms).cloned() {
            let renumbered: Vec<i32> = atoms
                .into_iter()
                .map(|v| {
                    if v < 0 {
                        -(rev[(-v - 1) as usize] as i32 + 1)
                    } else {
                        rev[(v - 1) as usize] as i32 + 1
                    }
                })
                .collect();
            atom.set_prop(cp::_ringStereoAtoms, renumbered, true)?;
        }
        if let Some(other) = atom
            .props
            .get(cp::_ringStereoOtherAtom)
            .and_then(RDValue::as_uint)
        {
            atom.set_prop(cp::_ringStereoOtherAtom, rev[other as usize], true)?;
        }
        res.add_atom(atom, false);
    }

    for bond in mol.bonds() {
        let mut b = bond.copy_detached();
        b.set_begin_atom_idx_unchecked(rev[bond.begin_atom_idx() as usize]);
        b.set_end_atom_idx_unchecked(rev[bond.end_atom_idx() as usize]);
        for idx in b.stereo_atoms_mut() {
            *idx = rev[*idx as usize] as i32;
        }
        res.add_bond_object(b)?;
    }

    for conf in mol.conformers() {
        let mut c = Conformer::new(n_atoms as usize);
        for (i, &old_idx) in new_order.iter().enumerate() {
            c.set_atom_pos(i as u32, conf.positions()[old_idx as usize])?;
        }
        c.set_id(conf.id());
        c.is_3d = conf.is_3d();
        res.add_conformer(c, false)?;
    }

    let old_rings = mol.ring_info();
    if old_rings.is_initialized() {
        let rings = res.ring_info_mut();
        rings.reset(true);
        rings.initialize(FindRingType::OtherOrUnknown);
        for (atom_ring, bond_ring) in old_rings.atom_rings().iter().zip(old_rings.bond_rings()) {
            let renumbered = atom_ring.iter().map(|&a| rev[a as usize] as i32).collect();
            rings.add_ring(renumbered, bond_ring.clone())?;
        }
    }

    if !mol.stereo_groups().is_empty() {
        let groups = mol
            .stereo_groups()
            .iter()
            .map(|g| {
                StereoGroup::new(
                    g.group_type(),
                    g.atoms().iter().map(|&a| rev[a as usize]).collect(),
                    // bonds keep their order
                    g.bonds().to_vec(),
                    g.read_id(),
                )
            })
            .collect::<Result<Vec<_>>>()?;
        res.set_stereo_groups(groups)?;
    }
    Ok(res)
}

/// `MolOps::getFormalCharge`: the sum of the atoms' formal charges.
pub fn get_formal_charge(mol: &ROMol) -> i32 {
    mol.atoms().iter().map(|a| a.formal_charge()).sum()
}

/// `MolOps::removeStereochemistry`: clears chiral tags, CIP labels, double
/// bond stereo, bond directions and stereo groups.
pub fn remove_stereochemistry(mol: &mut ROMol) -> Result<()> {
    mol.clear_prop(cp::_StereochemDone);
    for idx in 0..mol.num_atoms() {
        let atom = mol.atom_mut(idx);
        atom.set_chiral_tag(ChiralType::Unspecified);
        atom.clear_prop(cp::_CIPCode);
        atom.clear_prop(cp::_CIPRank);
    }
    for idx in 0..mol.num_bonds() {
        let bond = mol.bond_mut(idx);
        match bond.bond_type() {
            BondType::Double => {
                bond.set_stereo(BondStereo::StereoNone)?;
                bond.stereo_atoms_mut().clear();
                bond.set_bond_dir(BondDir::None);
            }
            BondType::Single => bond.set_bond_dir(BondDir::None),
            _ => {}
        }
    }
    mol.set_stereo_groups(Vec::new())
}
