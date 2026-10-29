//! Extended Hueckel calculations: a port of the molecular part of YAeHMOP's `tightbind` library
//! (<https://github.com/greglandrum/yaehmop>, BSD-2-Clause) as used by RDKit's `External/YAeHMOP`.

mod diag;
pub mod overlap;
mod params;
mod params_data;

pub use params::{Atom, atom_parameters};

use crate::error::{RDKitError, Result};
use crate::graphmol::romol::ROMol;

/// Weighting constant of the Wolfsberg-Helmholz formula
const THE_CONST: f64 = 1.75;
/// Overlaps between atoms further apart than this (angstrom) are dropped
const RHO: f64 = 10.0;
const DEGEN_TOL: f64 = 0.000001;

/// The outcome of a calculation (`EHTTools::EHTResults`)
#[derive(Debug, Clone, Default)]
pub struct EhtResults {
    pub num_atoms: usize,
    pub num_orbitals: usize,
    pub num_electrons: usize,
    pub fermi_energy: f64,
    pub total_energy: f64,
    /// `n x n` (row-major, upper triangle only like the C library), present when requested
    pub overlap_matrix: Option<Vec<f64>>,
    pub hamiltonian_matrix: Option<Vec<f64>>,
    /// `num_atoms x num_orbitals`
    pub reduced_charge_matrix: Vec<f64>,
    /// lower triangle of the atom x atom matrix, row by row
    pub reduced_overlap_population_matrix: Vec<f64>,
    pub atomic_charges: Vec<f64>,
    pub orbital_energies: Vec<f64>,
}

fn full_hamiltonian(atoms: &[Atom], tab: &[i32], n: usize, overlap: &[f64]) -> Vec<f64> {
    let mut hamil = vec![0.0f64; n * n];
    let mut diagonal = vec![0.0f64; n];
    for (a, &t) in atoms.iter().zip(tab) {
        let t = t as usize;
        let mut o = t;
        let mut put = |count: usize, coul: f64| {
            for k in 0..count {
                diagonal[o + k] = coul;
                hamil[(o + k) * n + o + k] = coul;
            }
            o += count;
        };
        if a.ns != 0 {
            put(1, a.coul_s);
        }
        if a.np != 0 {
            put(3, a.coul_p);
        }
        if a.nd != 0 {
            put(5, a.coul_d);
        }
        if a.nf != 0 {
            put(7, a.coul_f);
        }
    }
    for i in 1..n {
        for j in 0..i {
            let temp = diagonal[i] + diagonal[j];
            let mut temp2 = (diagonal[i] - diagonal[j]) / temp;
            temp2 *= temp2;
            temp2 = 0.5 * (THE_CONST + temp2 + temp2 * temp2 * (1.0 - THE_CONST));
            hamil[j * n + i] = temp2 * temp * overlap[j * n + i];
        }
    }
    hamil
}

/// The occupation numbers of the orbitals, degenerate levels sharing their electrons (`calc_occupations`)
fn calc_occupations(num_electrons: f64, eigenvalues: &[f64]) -> Vec<f64> {
    let n = eigenvalues.len();
    let mut occ = vec![0.0f64; n];
    let mut left = num_electrons;
    let mut i = 0;
    while i < n && left > 0.0 {
        if left >= 2.0 {
            occ[i] = 2.0;
            left -= 2.0;
        } else {
            occ[i] = left;
            left = 0.0;
        }
        i += 1;
    }
    share_degenerate(&mut occ, i, |k| eigenvalues[k]);
    occ
}

/// Averages the occupations over the degenerate levels around the last occupied one
fn share_degenerate(occ: &mut [f64], filled: usize, energy: impl Fn(usize) -> f64) {
    let n = occ.len();
    if filled == 0 {
        return;
    }
    let last = filled - 1;
    let mut end = filled;
    while end < n && (energy(last) - energy(end)).abs() < DEGEN_TOL {
        end += 1;
    }
    let mut begin = last as isize - 1;
    while begin > -1 && (energy(last) - energy(begin as usize)).abs() < DEGEN_TOL {
        begin -= 1;
    }
    let begin = (begin + 1) as usize;
    let levels = end - begin;
    if levels > 1 {
        let total: f64 = occ[begin..end].iter().sum();
        let each = total / levels as f64;
        for o in &mut occ[begin..end] {
            *o = each;
        }
    }
}

/// Runs an extended Hueckel calculation for the given atoms (element symbols and coordinates in angstrom)
/// and total formal charge.
pub fn run_atoms(
    symbols: &[String],
    positions: &[[f64; 3]],
    charge: i32,
    keep_matrices: bool,
) -> Result<EhtResults> {
    let nat = symbols.len();
    let mut atoms = Vec::with_capacity(nat);
    for s in symbols {
        match atom_parameters(s) {
            Some(a) => atoms.push(a),
            None => {
                return Err(RDKitError::Value(format!(
                    "Can't find parameters for atom: {s}. Parameter acquisition failure"
                )));
            }
        }
    }
    let num_electrons =
        atoms.iter().map(|a| f64::from(a.num_valence)).sum::<f64>() - f64::from(charge);
    let (tab, n) = overlap::orbital_lookup_table(&atoms);
    if num_electrons < 0.0 {
        return Err(RDKitError::Value(
            "Negative number of electrons specified.".into(),
        ));
    }
    if num_electrons == 0.0 {
        return Err(RDKitError::Value(
            "You forgot to specify the number of electrons.".into(),
        ));
    }
    if num_electrons > 2.0 * n as f64 {
        return Err(RDKitError::Value(
            "num_electrons is greater than twice the num_orbs. This is probably wrong".into(),
        ));
    }

    let ov = overlap::overlap_matrix(&atoms, positions, &tab, n, RHO);
    let hamil = full_hamiltonian(&atoms, &tab, n, &ov);
    let (eigenvalues, vect) = diag::generalized_eigen(n, &hamil, &ov).ok_or_else(|| {
        RDKitError::Value("Problems in the diagonalization, try more overlaps.".into())
    })?;

    // occupations and the Mulliken analysis
    let occ = calc_occupations(num_electrons, &eigenvalues);
    let s = |i: usize, j: usize| if j > i { ov[i * n + j] } else { ov[j * n + i] };
    let mut op = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut sum = 0.0;
            for k in 0..n {
                sum += occ[k] * vect[k * n + i] * vect[k * n + j];
            }
            op[i * n + j] = if i != j {
                2.0 * sum * s(i, j)
            } else {
                sum * s(i, j)
            };
        }
    }
    let atom_orbs = |i: usize| -> (usize, usize) {
        let end = if i + 1 < nat { tab[i + 1] as usize } else { n };
        (tab[i] as usize, end)
    };
    let mut net = vec![0.0f64; nat];
    let mut rop = Vec::with_capacity(nat * (nat + 1) / 2);
    for i in 0..nat {
        let (ib, ie) = atom_orbs(i);
        for j in 0..i {
            let (jb, je) = atom_orbs(j);
            let mut t = 0.0;
            for k in ib..ie {
                for l in jb..je {
                    t += op[k * n + l];
                }
            }
            rop.push(t);
        }
        let mut t = 0.0;
        for j in ib..ie {
            t += op[j * n + j];
        }
        rop.push(t);
    }
    // charge matrix
    let mut chg = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut ao = s(j, j) * vect[i * n + j];
            for k in 0..n {
                if k != j {
                    ao += s(j, k) * vect[i * n + k];
                }
            }
            chg[i * n + j] = vect[i * n + j] * ao;
        }
    }
    let mut rchg = vec![0.0f64; nat * n];
    for i in 0..nat {
        let (b, e) = atom_orbs(i);
        for j in 0..n {
            let mut t = 0.0;
            for k in b..e {
                t += chg[j * n + k];
            }
            rchg[i * n + j] = t;
        }
    }

    // Fermi energy and total energy come from the average-properties code, which works on
    // single-precision orbital energies
    let energies: Vec<f32> = eigenvalues.iter().map(|&e| e as f32).collect();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        energies[a]
            .partial_cmp(&energies[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let en = |k: usize| f64::from(energies[order[k]]);
    let mut crystal_occ = vec![0.0f64; n];
    let mut left = num_electrons;
    let mut filled = 0;
    while filled < n && left > 0.0 {
        if left >= 2.0 {
            crystal_occ[filled] = 2.0;
            left -= 2.0;
        } else {
            crystal_occ[filled] = left;
            left = 0.0;
        }
        filled += 1;
    }
    share_degenerate(&mut crystal_occ, filled, en);
    let fermi = en(filled - 1);
    // the orbital occupations and, from them, the net charges are accumulated from the single-precision
    // charge matrix of the average-properties code, which is what ends up in the results
    let mut total = 0.0;
    let mut ao_occups = vec![0.0f64; n];
    for (k, &o) in crystal_occ.iter().enumerate() {
        if o > 0.0001 {
            total += en(k) * o * 1.0;
            let mo = order[k];
            for (j, acc) in ao_occups.iter_mut().enumerate() {
                *acc += f64::from(chg[mo * n + j] as f32) * o * 1.0;
            }
        } else {
            break;
        }
    }
    let filled_bands: f64 = crystal_occ[..filled].iter().map(|o| o / 2.0).sum();
    for x in &mut ao_occups {
        *x = *x * filled_bands / (1.0 * filled_bands * 1.0);
    }
    for (i, a) in atoms.iter().enumerate() {
        let (b, e) = atom_orbs(i);
        net[i] = f64::from(a.num_valence);
        for x in &ao_occups[b..e] {
            net[i] -= x;
        }
    }

    let (overlap_matrix, hamiltonian_matrix) = if keep_matrices {
        (
            Some(ov.clone()),
            Some(full_hamiltonian(&atoms, &tab, n, &ov)),
        )
    } else {
        (None, None)
    };
    Ok(EhtResults {
        num_atoms: nat,
        num_orbitals: n,
        num_electrons: num_electrons.round() as usize,
        fermi_energy: fermi,
        total_energy: total,
        overlap_matrix,
        hamiltonian_matrix,
        reduced_charge_matrix: rchg,
        reduced_overlap_population_matrix: rop,
        atomic_charges: net,
        orbital_energies: eigenvalues,
    })
}

/// `EHTTools::runMol`
pub fn run_mol(mol: &ROMol, conf_id: i32, keep_matrices: bool) -> Result<EhtResults> {
    let conf = mol.conformer(conf_id)?;
    let n = mol.num_atoms();
    let mut symbols = Vec::with_capacity(n as usize);
    let mut pos = Vec::with_capacity(n as usize);
    for i in 0..n {
        symbols.push(mol.atom(i).symbol());
        let p = conf.atom_pos(i)?;
        pos.push([p.x, p.y, p.z]);
    }
    run_atoms(
        &symbols,
        &pos,
        crate::graphmol::mol_ops::get_formal_charge(mol),
        keep_matrices,
    )
}
