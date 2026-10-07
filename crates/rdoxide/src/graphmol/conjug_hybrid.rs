//! Ported from RDKit `Code/GraphMol/ConjugHybrid.cpp`.
//!
//! Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::aromaticity::count_atom_elec;
use super::atom::{ChiralType, HybridizationType};
use super::bond::BondType;
use super::periodic_table::PeriodicTable;
use super::romol::ROMol;

#[derive(Debug, Clone, Copy)]
struct ConjAtomInfo {
    num_substituents: u32,
    is_candidate: bool,
}

fn is_dative(bt: BondType) -> bool {
    matches!(
        bt,
        BondType::Dative | BondType::DativeL | BondType::DativeR | BondType::DativeOne
    )
}

fn is_atom_conjug_cand(mol: &ROMol, at: u32) -> bool {
    let pt = PeriodicTable::get();
    let a = mol.atom(at);
    let vals = &pt.element(a.atomic_num()).valences;
    if a.formal_charge() == 0 && vals[0] >= 0 && a.total_valence().unwrap_or(0) > vals[0] as u32 {
        return false;
    }
    let nouter = pt.element(a.atomic_num()).n_outer_elecs;
    (a.atomic_num() <= 10
        || (nouter != 5 && nouter != 6)
        || (nouter == 6 && mol.total_degree(at).unwrap_or(0) < 2))
        && count_atom_elec(mol, at) > 0
}

fn mark_conj_atom_bonds(mol: &mut ROMol, at: u32, info: &[ConjAtomInfo]) {
    let ai = info[at as usize];
    if !ai.is_candidate || ai.num_substituents < 2 || ai.num_substituents > 3 {
        return;
    }
    let bonds: Vec<u32> = mol.atom_bond_indices(at).collect();
    for &b1 in &bonds {
        let o1 = mol.bond(b1).other_atom_idx(at).expect("in bond");
        if mol.bond(b1).valence_contrib(mol, at) < 1.5 || !info[o1 as usize].is_candidate {
            continue;
        }
        for &b2 in &bonds {
            if b1 == b2 {
                continue;
            }
            let o2 = mol.bond(b2).other_atom_idx(at).expect("in bond");
            let a2 = info[o2 as usize];
            if a2.num_substituents > 3 {
                continue;
            }
            if a2.is_candidate {
                mol.bond_mut(b1).set_is_conjugated(true);
                mol.bond_mut(b2).set_is_conjugated(true);
            }
        }
    }
}

fn num_bonds_plus_lone_pairs(mol: &ROMol, at: u32) -> i32 {
    let a = mol.atom(at);
    let mut deg = mol.total_degree(at).unwrap_or(0) as i32;
    for b in mol.atom_bonds(at) {
        if b.bond_type() == BondType::Zero || (is_dative(b.bond_type()) && at != b.end_atom_idx()) {
            deg -= 1;
        }
    }
    if a.atomic_num() <= 1 {
        return deg;
    }
    let nouter = PeriodicTable::get().element(a.atomic_num()).n_outer_elecs;
    let total_valence = a.total_valence().unwrap_or(0) as i32;
    let chg = a.formal_charge();
    let num_free = nouter - (total_valence + chg);
    if total_valence + nouter - chg < 8 {
        let nrad = a.num_radical_electrons() as i32;
        let lone_pairs = (num_free - nrad) / 2;
        deg + lone_pairs + nrad
    } else {
        deg + num_free / 2
    }
}

/// `MolOps::atomHasConjugatedBond`
pub fn atom_has_conjugated_bond(mol: &ROMol, at: u32) -> bool {
    mol.atom_bonds(at).any(|b| b.is_conjugated())
}

/// `MolOps::setConjugation`
pub fn set_conjugation(mol: &mut ROMol) {
    for i in 0..mol.num_bonds() {
        let arom = mol.bond(i).is_aromatic();
        mol.bond_mut(i).set_is_conjugated(arom);
    }
    let info: Vec<ConjAtomInfo> = (0..mol.num_atoms())
        .map(|i| {
            let c = is_atom_conjug_cand(mol, i);
            ConjAtomInfo {
                num_substituents: if c {
                    mol.atom_degree(i) + mol.total_num_hs(i, false).unwrap_or(0)
                } else {
                    0
                },
                is_candidate: c,
            }
        })
        .collect();
    for i in 0..mol.num_atoms() {
        mark_conj_atom_bonds(mol, i, &info);
    }
}

/// `MolOps::setHybridization`
pub fn set_hybridization(mol: &mut ROMol) {
    for i in 0..mol.num_atoms() {
        if mol.atom(i).atomic_num() == 0 {
            mol.atom_mut(i)
                .set_hybridization(HybridizationType::Unspecified);
            continue;
        }
        let td = mol.total_degree(i).unwrap_or(0);
        let fixed = match mol.atom(i).chiral_tag() {
            ChiralType::Tetrahedral | ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW
                if td == 4 =>
            {
                Some(HybridizationType::SP3)
            }
            ChiralType::SquarePlanar if (2..=4).contains(&td) => Some(HybridizationType::SP2D),
            ChiralType::TrigonalBipyramidal if (2..=5).contains(&td) => {
                Some(HybridizationType::SP3D)
            }
            ChiralType::Octahedral if (2..=6).contains(&td) => Some(HybridizationType::SP3D2),
            _ => None,
        };
        if let Some(h) = fixed {
            mol.atom_mut(i).set_hybridization(h);
            continue;
        }
        let norbs = if mol.atom(i).atomic_num() < 89 {
            num_bonds_plus_lone_pairs(mol, i)
        } else {
            td as i32
        };
        let h = match norbs {
            0 | 1 => HybridizationType::S,
            2 => HybridizationType::SP,
            3 => HybridizationType::SP2,
            4 => {
                if td > 3 || !atom_has_conjugated_bond(mol, i) {
                    HybridizationType::SP3
                } else {
                    HybridizationType::SP2
                }
            }
            5 => HybridizationType::SP3D,
            6 => HybridizationType::SP3D2,
            _ => HybridizationType::Unspecified,
        };
        mol.atom_mut(i).set_hybridization(h);
    }
}
