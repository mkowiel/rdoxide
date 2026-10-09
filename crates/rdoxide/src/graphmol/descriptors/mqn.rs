//! Molecular Quantum Numbers (upstream `GraphMol/Descriptors/MQN.cpp`).
//!
//! 42 counts: atoms by element, ring membership and degree, bonds by order and
//! ring membership, polarity (acceptor/donor sites, charges), rotatable bonds
//! and ring sizes. As upstream, polarity uses RDKit's own donor/acceptor
//! definitions rather than the original paper's.
//!
//! Copyright (C) 2013 Greg Landrum.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::lipinski::{NumRotatableBondsOptions, calc_num_rotatable_bonds};
use crate::Result;
use crate::graphmol::bond::BondType;
use crate::graphmol::romol::ROMol;

pub const NUM_MQNS: usize = 42;

/// `calcMQNs(mol, force)` (the `force` flag is ignored upstream as well).
pub fn calc_mqns(mol: &ROMol) -> Result<Vec<u32>> {
    let mut res = vec![0u32; NUM_MQNS];
    let rings = mol.ring_info();

    // atom-centered counts
    for idx in 0..mol.num_atoms() {
        let atom = mol.atom(idx);
        let n_hs = mol.total_num_hs(idx, false)?;
        let n_rings = rings.num_atom_rings(idx)?;
        let in_ring = n_rings > 0;
        let charge = atom.formal_charge();
        match atom.atomic_num() {
            6 => res[0] += 1,
            9 => res[1] += 1,
            17 => res[2] += 1,
            35 => res[3] += 1,
            53 => res[4] += 1,
            16 => res[5] += 1,
            15 => res[6] += 1,
            7 => {
                res[if in_ring { 8 } else { 7 }] += 1;
                if mol.atom_degree(idx) != 4 {
                    res[19] += 1; // acceptor sites
                    res[20] += 1; // acceptor atoms
                }
                if n_hs > 0 {
                    res[21] += n_hs; // donor sites
                    res[22] += 1; // donor atoms
                }
            }
            8 => {
                res[if in_ring { 10 } else { 9 }] += 1;
                res[20] += 1; // acceptor atoms
                res[19] += if charge != -1 { 2 } else { 3 }; // acceptor sites
                if n_hs > 0 {
                    res[21] += n_hs;
                    res[22] += 1;
                }
            }
            _ => {}
        }
        if charge > 0 {
            res[24] += 1;
        } else if charge < 0 {
            res[23] += 1;
        }
        if atom.atomic_num() != 1 {
            // acyclic degree 1 to 4 count at 25..=28, cyclic degree 2 to 4 at 29..=31
            match (mol.atom_degree(idx), in_ring) {
                (1, _) => res[25] += 1,
                (d @ 2..=4, false) => res[24 + d as usize] += 1,
                (d @ 2..=4, true) => res[27 + d as usize] += 1,
                _ => {}
            }
            if n_rings >= 2 {
                res[40] += 1;
            }
        }
    }
    res[11] = mol.num_heavy_atoms();

    // bond counts
    let mut n_aromatic = 0u32;
    for bond in mol.bonds() {
        if bond.is_aromatic() {
            n_aromatic += 1;
        }
        let n_rings = rings.num_bond_rings(bond.idx())?;
        let offset = match bond.bond_type() {
            BondType::Single => Some(12),
            BondType::Double => Some(13),
            BondType::Triple => Some(14),
            _ => None,
        };
        if let Some(offset) = offset {
            res[if n_rings > 0 { offset + 3 } else { offset }] += 1;
        }
        if n_rings >= 2 {
            res[41] += 1;
        }
    }
    // rather than kekulizing, the aromatic bonds are split evenly between the
    // cyclic single and double bond bins, an odd one going to the single bonds
    res[15] += n_aromatic / 2 + n_aromatic % 2;
    res[16] += n_aromatic / 2;
    res[18] = calc_num_rotatable_bonds(mol, NumRotatableBondsOptions::Default)?;

    // ring size counts
    for ring in rings.atom_rings() {
        match ring.len() {
            n @ 0..=9 => res[n + 29] += 1,
            _ => res[39] += 1,
        }
    }
    Ok(res)
}
