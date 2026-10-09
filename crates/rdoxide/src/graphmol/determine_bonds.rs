//! Bond perception from 3D coordinates (the xyz2mol algorithm).
//!
//! Ported from RDKit `Code/GraphMol/DetermineBonds/DetermineBonds.{h,cpp}` (2026.09.1),
//! which is based on xyz2mol (https://github.com/jensengroup/xyz2mol, MIT).
//! Extended Hueckel connectivity (YAeHMOP) is not available.
//!
//! Copyright (C) 2022 Rachel Walker, Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::HashMap;

use super::bond::BondType;
use super::chirality::bond_dirs::set_double_bond_neighbor_directions;
use super::chirality::from_coords::assign_stereochemistry_from_3d;
use super::file_parsers::proximity_bonds::{CTD_QUICKREMOVE_H_H_CONTACTS, connect_the_dots};
use super::matrices::get_3d_distance_mat;
use super::mol_ops::sanitize_flags::SANITIZE_ALL;
use super::mol_ops::{get_formal_charge, sanitize_mol};
use super::periodic_table::PeriodicTable;
use super::romol::RWMol;
use crate::rdgeneral::matching::max_cardinality_matching;
use crate::{RDKitError, Result};

type Matrix = Vec<Vec<u32>>;

/// Whether extended Hueckel connectivity (needs YAeHMOP) is available.
pub const HUECKEL_SUPPORT: bool = false;

fn max_iterations_exceeded() -> RDKitError {
    // `MaxFindBondOrdersItersExceeded` is a std::runtime_error
    RDKitError::Runtime("Max Iterations Exceeded in Determine Bond Orders".into())
}

/// The cartesian product of the candidate valences, enumerated with the last
/// atom varying fastest (`LazyCartesianProduct`).
struct ValenceProduct {
    lists: Vec<Vec<u32>>,
    position: Vec<usize>,
    at_end: bool,
}

impl ValenceProduct {
    fn new(lists: Vec<Vec<u32>>) -> Self {
        let at_end = lists.iter().any(Vec::is_empty);
        ValenceProduct {
            position: vec![0; lists.len()],
            lists,
            at_end,
        }
    }

    fn at_end(&self) -> bool {
        self.at_end
    }

    fn next_entry(&mut self) -> Vec<u32> {
        let entry = self
            .lists
            .iter()
            .zip(&self.position)
            .map(|(list, &p)| list[p])
            .collect();
        // odometer increment
        let mut i = self.lists.len();
        loop {
            if i == 0 {
                self.at_end = true;
                break;
            }
            i -= 1;
            self.position[i] += 1;
            if self.position[i] < self.lists[i].len() {
                break;
            }
            self.position[i] = 0;
        }
        entry
    }
}

fn possible_valences(atomic_num: i32, known: &HashMap<i32, Vec<u32>>) -> Result<Vec<u32>> {
    if let Some(v) = known.get(&atomic_num) {
        return Ok(v.clone());
    }
    let list = PeriodicTable::get().get_valence_list(atomic_num as u32)?;
    Ok(list
        .iter()
        .filter(|&&v| v >= 0)
        .map(|&v| v as u32)
        .collect())
}

fn valence_combinations(mol: &RWMol, mut iterations: usize) -> Result<ValenceProduct> {
    let known: HashMap<i32, Vec<u32>> = HashMap::from([
        (1, vec![1]),
        (5, vec![3, 4]),
        (6, vec![4]),
        (7, vec![3, 4]),
        (8, vec![2, 1, 3]),
        (9, vec![1]),
        (14, vec![4]),
        (15, vec![5, 3]),
        (16, vec![6, 3, 2, 1]),
        (17, vec![1]),
        (32, vec![4]),
        (35, vec![1]),
        (53, vec![1]),
    ]);
    let n = mol.num_atoms() as usize;
    let degree = |i: usize| mol.atom_bonds(i as u32).count() as u32;
    let mut current: Vec<Vec<u32>> = Vec::with_capacity(n);
    for i in 0..n {
        let atom = mol.atom(i as u32);
        let valences = possible_valences(atom.atomic_num() as i32 - atom.formal_charge(), &known)?;
        if valences.is_empty() {
            return Err(RDKitError::Value(format!(
                "Atom {i} with atomic number {} has no valences defined.",
                atom.atomic_num()
            )));
        }
        current.push(valences);
    }
    loop {
        let mut next: Vec<Vec<u32>> = vec![Vec::new(); n];
        for i in 0..n {
            if current[i].len() == 1 {
                next[i] = current[i].clone();
                continue;
            }
            let mut avail_bonds = 0i64;
            for bond in mol.atom_bonds(i as u32) {
                let nbr = bond.other_atom_idx(i as u32)? as usize;
                let max_pos = i64::from(*current[nbr].iter().max().expect("non-empty"));
                avail_bonds += max_pos - i64::from(degree(nbr));
            }
            let num_bonds = i64::from(degree(i));
            next[i] = current[i]
                .iter()
                .copied()
                .filter(|&v| num_bonds <= i64::from(v) && i64::from(v) <= num_bonds + avail_bonds)
                .collect();
            if next[i].is_empty() {
                return Err(RDKitError::Value(format!(
                    "Unable to determine valence of atom {i} with atomic number {} and {} bonds. Check the input molecules bonding.",
                    mol.atom(i as u32).atomic_num(),
                    degree(i)
                )));
            }
        }
        if next == current {
            break;
        }
        current = next;
        iterations -= 1;
        if iterations == 0 {
            break;
        }
    }
    Ok(ValenceProduct::new(current))
}

/// connectivity from the sum of the (scaled) covalent radii
fn connectivity_vdw(mol: &mut RWMol, cov_factor: f64) -> Result<()> {
    let n = mol.num_atoms() as usize;
    let dist = get_3d_distance_mat(mol, -1)?;
    let table = PeriodicTable::get();
    let rcov = (0..n)
        .map(|i| Ok(cov_factor * table.get_rcovalent(mol.atom(i as u32).atomic_num())?))
        .collect::<Result<Vec<f64>>>()?;
    for i in 0..n {
        for j in i + 1..n {
            if dist[i * n + j] <= rcov[i] + rcov[j] {
                mol.add_bond(i as u32, j as u32, BondType::Single)?;
            }
        }
    }
    Ok(())
}

/// `determineConnectivity`: replaces the bonds of `mol` by the ones perceived
/// from its (3D) default conformer.
pub fn determine_connectivity(
    mol: &mut RWMol,
    use_hueckel: bool,
    _charge: i32,
    cov_factor: f64,
    use_vdw: bool,
) -> Result<()> {
    if use_hueckel {
        return Err(RDKitError::Value(
            "The RDKit was not compiled with YAeHMOP support".into(),
        ));
    }
    mol.begin_batch_edit()?;
    let bonds: Vec<(u32, u32)> = mol
        .bonds()
        .iter()
        .map(|b| (b.begin_atom_idx(), b.end_atom_idx()))
        .collect();
    for (a, b) in bonds {
        mol.remove_bond(a, b)?;
    }
    mol.commit_batch_edit()?;
    for idx in 0..mol.num_atoms() {
        mol.atom_mut(idx).set_no_implicit(true);
    }
    if use_vdw {
        connectivity_vdw(mol, cov_factor)
    } else {
        connect_the_dots(mol, CTD_QUICKREMOVE_H_H_CONTACTS)
    }
}

fn unsaturated(order: &[u32], valency: &[u32]) -> Vec<usize> {
    (0..order.len())
        .filter(|&i| order[i] > valency[i])
        .collect()
}

fn unsaturated_pairs(con_mat: &Matrix, unsat: &[usize]) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for (i, &a) in unsat.iter().enumerate() {
        for &b in &unsat[i + 1..] {
            if con_mat[a][b] != 0 {
                pairs.push((a, b));
            }
        }
    }
    pairs
}

fn check_valency(order: &[u32], valency: &[u32]) -> bool {
    valency.iter().zip(order).all(|(v, o)| v <= o)
}

fn check_saturation(order: &[u32], valency: &[u32]) -> bool {
    order.iter().zip(valency).all(|(o, v)| o <= v)
}

fn atomic_charge(atomic_num: i32, valence: u32) -> Result<i32> {
    let valence = valence as i32;
    Ok(match (atomic_num, valence) {
        (1, _) => 1 - valence,
        (5, _) => 3 - valence,
        (15, 5) | (16, 6) => 0,
        _ => PeriodicTable::get().get_n_outer_elecs(atomic_num as u32)? - 8 + valence,
    })
}

fn check_charge(mol: &RWMol, valency: &[u32], charge: i32) -> Result<bool> {
    let mut mol_charge = 0;
    for i in 0..mol.num_atoms() {
        let atom = mol.atom(i);
        mol_charge += atomic_charge(atom.atomic_num() as i32, valency[i as usize])?;
        if atom.atomic_num() == 6 {
            let degree = mol.atom_bonds(i).count();
            if degree == 2 && valency[i as usize] == 2 {
                mol_charge += 1;
            } else if degree == 3 && mol_charge + 1 < charge {
                mol_charge += 2;
            }
        }
    }
    Ok(mol_charge == charge)
}

fn set_atomic_charges(mol: &mut RWMol, valency: &[u32], charge: i32) -> Result<()> {
    let mut mol_charge = 0;
    for i in 0..mol.num_atoms() {
        let atom = mol.atom(i);
        let mut atom_charge = atomic_charge(
            atom.atomic_num() as i32 - atom.formal_charge(),
            valency[i as usize],
        )?;
        mol_charge += atom_charge;
        if atom.atomic_num() == 6 {
            let degree = mol.atom_bonds(i).count();
            if degree == 2 && valency[i as usize] == 2 {
                mol_charge += 1;
                atom_charge = 0;
            } else if degree == 3 && mol_charge + 1 < charge {
                mol_charge += 2;
                atom_charge = 1;
            }
        }
        if atom_charge != 0 {
            mol.atom_mut(i).set_formal_charge(atom_charge);
        }
    }
    Ok(())
}

fn set_atomic_radicals(mol: &mut RWMol, valency: &[u32], charge: i32) -> Result<()> {
    for i in 0..mol.num_atoms() {
        let atom_charge = atomic_charge(mol.atom(i).atomic_num() as i32, valency[i as usize])?;
        if atom_charge != 0 {
            mol.atom_mut(i)
                .set_num_radical_electrons(charge.unsigned_abs());
        }
    }
    Ok(())
}

fn set_chirality(mol: &mut RWMol) -> Result<()> {
    sanitize_mol(mol, SANITIZE_ALL).map_err(|(_, e)| e)?;
    set_double_bond_neighbor_directions(mol, Some(0))?;
    assign_stereochemistry_from_3d(mol, -1, true)
}

fn add_bond_ordering(
    mol: &mut RWMol,
    ord_mat: &Matrix,
    valency: &[u32],
    allow_charged_fragments: bool,
    embed_chiral: bool,
    use_atom_map: bool,
    charge: i32,
) -> Result<()> {
    let n = mol.num_atoms() as usize;
    for i in 0..n {
        for j in i + 1..n {
            let bond_type = match ord_mat[i][j] {
                0 => continue,
                2 => BondType::Double,
                3 => BondType::Triple,
                _ => BondType::Single,
            };
            let idx = mol
                .get_bond_between_atoms(i as u32, j as u32)?
                .map(|b| b.idx())
                .ok_or_else(|| RDKitError::invariant("missing bond"))?;
            mol.bond_mut(idx).set_bond_type(bond_type);
        }
    }
    if allow_charged_fragments {
        set_atomic_charges(mol, valency, charge)?;
    } else {
        set_atomic_radicals(mol, valency, charge)?;
    }
    let final_charge = get_formal_charge(mol);
    if final_charge != charge {
        return Err(RDKitError::Value(format!(
            "Final molecular charge ({charge}) does not match input ({final_charge}); could not find valid bond ordering"
        )));
    }
    if use_atom_map {
        for i in 0..mol.num_atoms() {
            mol.atom_mut(i).set_atom_map_num(i as i32 + 1, false)?;
        }
    }
    if embed_chiral {
        set_chirality(mol)?;
    }
    Ok(())
}

/// Adds pair bonds between the unsaturated atoms until no more can be added.
fn saturate(con_mat: &Matrix, order: &[u32], ord_mat: &mut Matrix, valency: &mut [u32]) {
    let n = order.len();
    loop {
        let unsat = unsaturated(order, valency);
        let pairs = unsaturated_pairs(con_mat, &unsat);
        if pairs.is_empty() {
            return;
        }
        let mate = max_cardinality_matching(n, &pairs);
        let mut new_bonds = false;
        for (v, m) in mate.iter().enumerate() {
            if let Some(m) = *m
                && v < m
            {
                new_bonds = true;
                ord_mat[v][m] += 1;
                valency[v] += 1;
                valency[m] += 1;
            }
        }
        if !new_bonds {
            return;
        }
    }
}

/// `determineBondOrders`: assigns bond orders, charges (or radicals) to a
/// molecule with known connectivity. `max_iterations` of 0 means no limit.
pub fn determine_bond_orders(
    mol: &mut RWMol,
    charge: i32,
    allow_charged_fragments: bool,
    embed_chiral: bool,
    use_atom_map: bool,
    max_iterations: usize,
) -> Result<()> {
    let mut iterations = if max_iterations == 0 {
        usize::MAX
    } else {
        max_iterations
    };
    let n = mol.num_atoms() as usize;
    let mut con_mat: Matrix = vec![vec![0; n]; n];
    let mut orig_valency = vec![0u32; n];
    for bond in mol.bonds() {
        let (mut i, mut j) = (bond.begin_atom_idx() as usize, bond.end_atom_idx() as usize);
        if i > j {
            std::mem::swap(&mut i, &mut j);
        }
        con_mat[i][j] += 1;
        orig_valency[i] += 1;
        orig_valency[j] += 1;
    }
    let mut best = con_mat.clone();
    let mut best_valency = orig_valency.clone();
    let mut best_sum: u32 = orig_valency.iter().sum();

    let mut combos = valence_combinations(mol, 100)?;
    while !combos.at_end() {
        iterations -= 1;
        if iterations == 0 {
            return Err(max_iterations_exceeded());
        }
        let order = combos.next_entry();
        if unsaturated(&order, &orig_valency).is_empty() {
            if check_valency(&order, &orig_valency)
                && check_charge(mol, &orig_valency, charge)?
                && check_saturation(&order, &orig_valency)
            {
                return add_bond_ordering(
                    mol,
                    &con_mat,
                    &orig_valency,
                    allow_charged_fragments,
                    embed_chiral,
                    use_atom_map,
                    charge,
                );
            }
            continue;
        }
        let mut ord_mat = con_mat.clone();
        let mut valency = orig_valency.clone();
        saturate(&con_mat, &order, &mut ord_mat, &mut valency);
        if check_valency(&order, &valency) && check_charge(mol, &valency, charge)? {
            if check_saturation(&order, &valency) {
                return add_bond_ordering(
                    mol,
                    &ord_mat,
                    &valency,
                    allow_charged_fragments,
                    embed_chiral,
                    use_atom_map,
                    charge,
                );
            }
            let sum: u32 = valency.iter().sum();
            if sum > best_sum {
                best = ord_mat;
                best_sum = sum;
                best_valency = valency;
            }
        }
    }
    add_bond_ordering(
        mol,
        &best,
        &best_valency,
        allow_charged_fragments,
        embed_chiral,
        use_atom_map,
        charge,
    )
}

/// `determineBonds`: connectivity followed by bond orders.
#[allow(clippy::too_many_arguments)]
pub fn determine_bonds(
    mol: &mut RWMol,
    use_hueckel: bool,
    charge: i32,
    cov_factor: f64,
    allow_charged_fragments: bool,
    embed_chiral: bool,
    use_atom_map: bool,
    use_vdw: bool,
    max_iterations: usize,
) -> Result<()> {
    if mol.num_atoms() <= 1 {
        return Ok(());
    }
    determine_connectivity(mol, use_hueckel, charge, cov_factor, use_vdw)?;
    determine_bond_orders(
        mol,
        charge,
        allow_charged_fragments,
        embed_chiral,
        use_atom_map,
        max_iterations,
    )
}
