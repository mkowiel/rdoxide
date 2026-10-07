//! Atom/bond query predicates (in progress).
//!
//! Ported from RDKit `Code/GraphMol/QueryOps.{h,cpp}` (2026.09.1).
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::graphmol::romol::ROMol;

/// `queryIsAtomBridgehead`: at least three ring bonds, and every ring
/// containing the atom shares at least two bonds with another ring that also
/// contains the atom.
pub fn query_is_atom_bridgehead(mol: &ROMol, atom: u32) -> bool {
    if mol.atom_degree(atom) < 3 {
        return false;
    }
    let ri = mol.ring_info();
    if !ri.is_initialized() {
        return false;
    }
    let nb = mol.num_bonds() as usize;
    let mut atom_ring_bonds = vec![false; nb];
    let mut count = 0;
    for b in mol.atom_bond_indices(atom) {
        if ri.num_bond_rings(b).unwrap_or(0) > 0 {
            atom_ring_bonds[b as usize] = true;
            count += 1;
        }
    }
    if count < 3 {
        return false;
    }
    let bond_rings = ri.bond_rings();
    let mut bonds_in_ring_i = vec![false; nb];
    let mut rings_overlap = vec![false; bond_rings.len()];
    for i in 0..bond_rings.len() {
        bonds_in_ring_i.iter_mut().for_each(|x| *x = false);
        let mut atom_in_ring_i = false;
        for &bidx in &bond_rings[i] {
            bonds_in_ring_i[bidx as usize] = true;
            if atom_ring_bonds[bidx as usize] {
                atom_in_ring_i = true;
            }
        }
        if !atom_in_ring_i {
            continue;
        }
        for j in i + 1..bond_rings.len() {
            let mut overlap = 0;
            let mut atom_in_ring_j = false;
            for &bidx in &bond_rings[j] {
                if atom_ring_bonds[bidx as usize] {
                    atom_in_ring_j = true;
                }
                if bonds_in_ring_i[bidx as usize] {
                    overlap += 1;
                }
                if overlap >= 2 && atom_in_ring_j {
                    rings_overlap[i] = true;
                    rings_overlap[j] = true;
                    break;
                }
            }
        }
        if !rings_overlap[i] {
            return false;
        }
    }
    true
}
