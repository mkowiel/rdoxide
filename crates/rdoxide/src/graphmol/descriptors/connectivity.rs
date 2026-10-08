//! Hall-Kier alpha and Phi (upstream `GraphMol/Descriptors/ConnectivityDescriptors.cpp`).
//! The other connectivity descriptors (Chi, Kappa3, ...) are not ported yet.

use crate::Result;
use crate::graphmol::atom::{Atom, HybridizationType};
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::ROMol;
use crate::graphmol::subgraphs::find_all_paths_of_length_n;

pub const PHI_VERSION: &str = "1.0.0";

/// `detail::getAlpha`: the tabulated Hall-Kier alpha of an atom, `None` when the
/// element is not in the table.
fn get_alpha(atom: &Atom) -> Option<f64> {
    use HybridizationType::{SP, SP2};
    let hyb = atom.hybridization();
    Some(match atom.atomic_num() {
        1 => 0.0,
        6 => match hyb {
            SP => -0.22,
            SP2 => -0.13,
            _ => 0.00,
        },
        7 => match hyb {
            SP => -0.29,
            SP2 => -0.20,
            _ => -0.04,
        },
        8 => match hyb {
            SP2 => -0.20,
            _ => -0.04,
        },
        9 => -0.07,
        15 => match hyb {
            SP2 => 0.30,
            _ => 0.43,
        },
        16 => match hyb {
            SP2 => 0.22,
            _ => 0.35,
        },
        17 => 0.29,
        35 => 0.48,
        53 => 0.73,
        _ => return None,
    })
}

/// `calcHallKierAlpha`: the sum of the atomic alphas; elements outside the table
/// use `rA / rC - 1` with the covalent radii (`Rb0`) of the element and carbon.
pub fn calc_hall_kier_alpha(mol: &ROMol) -> Result<f64> {
    let tbl = PeriodicTable::get();
    let r_c = tbl.get_rb0(6)?;
    let mut alpha_sum = 0.0;
    for i in 0..mol.num_atoms() {
        let atom = mol.atom(i);
        let n = atom.atomic_num();
        if n == 0 {
            continue;
        }
        let alpha = match get_alpha(atom) {
            Some(a) => a,
            None => tbl.get_rb0(n)? / r_c - 1.0,
        };
        alpha_sum += alpha;
    }
    Ok(alpha_sum)
}

fn kappa1_helper(p1: f64, a: f64, alpha: f64) -> f64 {
    let denom = p1 + alpha;
    if denom != 0.0 {
        (a + alpha) * (a + alpha - 1.0) * (a + alpha - 1.0) / (denom * denom)
    } else {
        0.0
    }
}

fn kappa2_helper(p2: f64, a: f64, alpha: f64) -> f64 {
    let denom = (p2 + alpha) * (p2 + alpha);
    if denom != 0.0 {
        (a + alpha - 1.0) * (a + alpha - 2.0) * (a + alpha - 2.0) / denom
    } else {
        0.0
    }
}

/// `calcPhi`: kappa1 * kappa2 / A over the heavy atoms.
pub fn calc_phi(mol: &ROMol) -> Result<f64> {
    let heavy = mol.num_heavy_atoms();
    if heavy == 0 {
        return Ok(0.0);
    }
    let alpha = calc_hall_kier_alpha(mol)?;
    let a = f64::from(heavy);
    let kappa1 = kappa1_helper(f64::from(mol.num_bonds()), a, alpha);
    let p2 = find_all_paths_of_length_n(mol, 2, true, false, -1, false, None)?.len();
    let kappa2 = kappa2_helper(p2 as f64, a, alpha);
    Ok(kappa1 * kappa2 / a)
}
