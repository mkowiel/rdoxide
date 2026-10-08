//! Gasteiger-Marsili partial charges (upstream `GasteigerCharges.{h,cpp}`).
//!
//! Ref: J. Gasteiger, M. Marsili, "Iterative Equalization of Orbital
//! Electronegativity - A Rapid Access to Atomic Charges", Tetrahedron 36,
//! 3219 (1980). The floating-point expressions keep upstream's operation
//! order and fused multiply-adds (see [`crate::rdgeneral::fp`]).

use super::gasteiger_params::{DAMP, DAMP_SCALE, IONXH, get_params};
use crate::error::Result;
use crate::graphmol::atom::HybridizationType;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::dict::RDProps;
use crate::rdgeneral::fp::fmuladd;

/// Default number of iterations (`nIter = 12`).
pub const DEFAULT_ITERATIONS: i32 = 12;

/// `Gasteiger::splitChargeConjugated`: formal charges of atoms in a conjugated
/// system are shared equally between the atoms of the same element that are
/// two conjugated bonds away, so e.g. both nitrogens of a benzamidine start
/// the iteration with a charge of 0.5.
fn split_charge_conjugated(mol: &ROMol, charges: &mut [f64]) -> Result<()> {
    for at in mol.atoms() {
        let aix = at.idx();
        let mut formal = f64::from(at.formal_charge());
        if formal.abs() > f64::EPSILON && charges[aix as usize].abs() < f64::EPSILON {
            let mut marker = vec![aix];
            for bnd1 in mol.atom_bonds(aix).filter(|b| b.is_conjugated()) {
                let aax = bnd1.other_atom_idx(aix)?;
                for bnd2 in mol.atom_bonds(aax) {
                    if bnd1.idx() != bnd2.idx() && bnd2.is_conjugated() {
                        let yax = bnd2.other_atom_idx(aax)?;
                        let yat = mol.atom(yax);
                        if at.atomic_num() == yat.atomic_num() {
                            formal += f64::from(yat.formal_charge());
                            marker.push(yax);
                        }
                    }
                }
            }
            let share = formal / marker.len() as f64;
            for mci in marker {
                charges[mci as usize] = share;
            }
        }
    }
    Ok(())
}

/// The parameter set name ("mode") of an atom.
fn atom_mode(mol: &ROMol, idx: u32) -> &'static str {
    let atom = mol.atom(idx);
    match atom.hybridization() {
        HybridizationType::SP3 => "sp3",
        HybridizationType::SP2 => "sp2",
        HybridizationType::SP => "sp",
        _ if atom.atomic_num() == 1 => "*",
        // a sulfur without hybridization information is classified by its oxygens
        _ if atom.atomic_num() == 16 => {
            let oxygens = mol
                .atom_neighbors(idx)
                .filter(|&n| mol.atom(n).atomic_num() == 8)
                .count();
            match oxygens {
                2 => "so2",
                1 => "so",
                _ => "sp3",
            }
        }
        // unknown: no parameters, so the lookup fails (or yields zeros)
        _ => "",
    }
}

/// One charge-transfer term, `dx / ((sgn * (a - b)) + b)` with `sgn` 0 for
/// `dx < 0` and 1 otherwise.
fn charge_flow(dx: f64, ion_x_this: f64, ion_x_other: f64) -> f64 {
    let sgn = if dx < 0.0 { 0.0 } else { 1.0 };
    dx / fmuladd(sgn, ion_x_this - ion_x_other, ion_x_other)
}

/// `computeGasteigerCharges(mol, charges, nIter, throwOnParamFailure)`.
///
/// Returns the charge of every atom and stores it on the atoms as the
/// computed property `_GasteigerCharge`, together with the total charge of the
/// atom's implicit hydrogens as `_GasteigerHCharge`.
pub fn compute_gasteiger_charges(
    mol: &mut ROMol,
    n_iter: i32,
    throw_on_param_failure: bool,
) -> Result<Vec<f64>> {
    let table = PeriodicTable::get();
    let natms = mol.num_atoms() as usize;
    let mut charges = vec![0.0; natms];
    // total charge of the implicit hydrogens on each heavy atom
    let mut h_chrg = vec![0.0; natms];
    let mut ion_x = vec![0.0; natms];
    let mut energ = vec![0.0; natms];

    split_charge_conjugated(mol, &mut charges)?;

    let mut atm_ps = Vec::with_capacity(natms);
    for idx in 0..natms {
        let atom = mol.atom(idx as u32);
        let elem = table.get_element_symbol(atom.atomic_num())?;
        let ps = get_params(elem, atom_mode(mol, idx as u32), throw_on_param_failure)?;
        ion_x[idx] = if atom.atomic_num() == 1 {
            IONXH
        } else {
            ps[0] + ps[1] + ps[2]
        };
        atm_ps.push(ps);
    }
    // hydrogens that are not in the graph (implicit)
    let h_params = get_params("H", "*", throw_on_param_failure)?;
    let n_hs = (0..natms as u32)
        .map(|i| mol.total_num_hs(i, false))
        .collect::<Result<Vec<u32>>>()?;

    let mut damp = DAMP;
    for _ in 0..n_iter as u32 {
        for aix in 0..natms {
            let p = atm_ps[aix];
            energ[aix] = fmuladd(charges[aix], fmuladd(p[2], charges[aix], p[1]), p[0]);
        }
        for aix in 0..natms {
            let mut dq = 0.0;
            for nbr in mol.atom_neighbors(aix as u32) {
                let nbr = nbr as usize;
                dq += charge_flow(energ[nbr] - energ[aix], ion_x[aix], ion_x[nbr]);
            }
            // the implicit hydrogens' share; as they connect to nothing else
            // their charges are updated at the same time
            let ni_hs = f64::from(n_hs[aix]);
            if n_hs[aix] > 0 {
                let q_hs = h_chrg[aix] / ni_hs;
                let enr = fmuladd(q_hs, fmuladd(h_params[2], q_hs, h_params[1]), h_params[0]);
                let dq_h = charge_flow(enr - energ[aix], ion_x[aix], IONXH);
                dq = fmuladd(ni_hs, dq_h, dq);
                h_chrg[aix] = fmuladd(-(ni_hs * dq_h), damp, h_chrg[aix]);
            }
            charges[aix] = fmuladd(damp, dq, charges[aix]);
        }
        damp *= DAMP_SCALE;
    }

    for (aix, (&charge, &h_charge)) in charges.iter().zip(&h_chrg).enumerate() {
        let atom = mol.atom_mut(aix as u32);
        atom.set_prop(cp::_GasteigerCharge, charge, true)?;
        atom.set_prop(cp::_GasteigerHCharge, h_charge, true)?;
    }
    Ok(charges)
}
