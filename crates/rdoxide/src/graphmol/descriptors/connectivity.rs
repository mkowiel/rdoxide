//! Connectivity descriptors: Chi, Hall-Kier alpha, Kappa and Phi (upstream
//! `GraphMol/Descriptors/ConnectivityDescriptors.cpp`).
//!
//! From Rev. Comp. Chem. vol 2, 367-422 (1991). The per-atom delta values of
//! the Chi indices are cached on the molecule like upstream does.

use crate::Result;
use crate::error::RDKitError;
use crate::graphmol::atom::{Atom, HybridizationType};
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::ROMol;
use crate::graphmol::subgraphs::find_all_paths_of_length_n;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::dict::{RDProps, RDValue};
use crate::rdgeneral::fp::fmuladd;

pub const CHI_VERSION: &str = "1.2.0";
pub const HALL_KIER_ALPHA_VERSION: &str = "1.2.0";
pub const KAPPA_VERSION: &str = "1.1.0";
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
/// With `atom_contribs` (at least one slot per atom) the alpha of each atom is
/// stored too; dummy atoms do not contribute and leave their slot untouched.
pub fn calc_hall_kier_alpha(mol: &ROMol, mut atom_contribs: Option<&mut [f64]>) -> Result<f64> {
    if atom_contribs
        .as_ref()
        .is_some_and(|c| c.len() < mol.num_atoms() as usize)
    {
        return Err(RDKitError::precondition("bad atomContribs vector"));
    }
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
        if let Some(contribs) = atom_contribs.as_deref_mut() {
            contribs[i as usize] = alpha;
        }
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

fn kappa3_helper(p3: f64, a: i32, alpha: f64) -> f64 {
    let denom = (p3 + alpha) * (p3 + alpha);
    if denom == 0.0 {
        return 0.0;
    }
    let a_f = f64::from(a);
    if a % 2 != 0 {
        (a_f + alpha - 1.0) * (a_f + alpha - 3.0) * (a_f + alpha - 3.0) / denom
    } else {
        (a_f + alpha - 2.0) * (a_f + alpha - 3.0) * (a_f + alpha - 3.0) / denom
    }
}

/// `calcKappa1`
pub fn calc_kappa1(mol: &ROMol) -> Result<f64> {
    let alpha = calc_hall_kier_alpha(mol, None)?;
    Ok(kappa1_helper(
        f64::from(mol.num_bonds()),
        f64::from(mol.num_heavy_atoms()),
        alpha,
    ))
}

/// `calcKappa2`
pub fn calc_kappa2(mol: &ROMol) -> Result<f64> {
    let p2 = find_all_paths_of_length_n(mol, 2, true, false, -1, false, None)?.len();
    let alpha = calc_hall_kier_alpha(mol, None)?;
    Ok(kappa2_helper(
        p2 as f64,
        f64::from(mol.num_heavy_atoms()),
        alpha,
    ))
}

/// `calcKappa3`
pub fn calc_kappa3(mol: &ROMol) -> Result<f64> {
    let p3 = find_all_paths_of_length_n(mol, 3, true, false, -1, false, None)?.len();
    let a = i32::try_from(mol.num_heavy_atoms()).expect("atom count fits in i32");
    let alpha = calc_hall_kier_alpha(mol, None)?;
    Ok(kappa3_helper(p3 as f64, a, alpha))
}

/// `calcPhi`: kappa1 * kappa2 / A over the heavy atoms.
pub fn calc_phi(mol: &ROMol) -> Result<f64> {
    let heavy = mol.num_heavy_atoms();
    if heavy == 0 {
        return Ok(0.0);
    }
    let alpha = calc_hall_kier_alpha(mol, None)?;
    let a = f64::from(heavy);
    let kappa1 = kappa1_helper(f64::from(mol.num_bonds()), a, alpha);
    let p2 = find_all_paths_of_length_n(mol, 2, true, false, -1, false, None)?.len();
    let kappa2 = kappa2_helper(p2 as f64, a, alpha);
    Ok(kappa1 * kappa2 / a)
}

/// `1 / sqrt(v)`, 0 for 0: the delta of an atom whose valence term is `v`.
fn inverse_sqrt(v: f64) -> f64 {
    if v != 0.0 { 1.0 / v.sqrt() } else { v }
}

/// The per-atom values `compute` yields, cached on the molecule under `key`.
fn cached_atom_values(
    mol: &mut ROMol,
    key: &str,
    force: bool,
    compute: impl Fn(&ROMol, &Atom, u32) -> Result<f64>,
) -> Result<Vec<f64>> {
    if !force && let Some(RDValue::VecDouble(values)) = mol.get_prop(key) {
        return Ok(values.clone());
    }
    let values = (0..mol.num_atoms())
        .map(|i| compute(mol, mol.atom(i), i))
        .collect::<Result<Vec<f64>>>()?;
    mol.set_prop(key, values.clone(), true)?;
    Ok(values)
}

/// `detail::hkDeltas`: the Hall-Kier valence deltas `1/sqrt(delta_v)`. Upstream
/// computes `nOuterElecs - nHs` in unsigned arithmetic; that is kept (it only
/// matters for atoms with more hydrogens than outer electrons).
fn hk_deltas(mol: &mut ROMol, force: bool) -> Result<Vec<f64>> {
    let tbl = PeriodicTable::get();
    cached_atom_values(mol, cp::_connectivityHKDeltas, force, |mol, atom, idx| {
        let n = atom.atomic_num();
        let outer = tbl.get_n_outer_elecs(n)? as u32;
        let n_hs = mol.total_num_hs(idx, false)?;
        let delta = if n <= 1 {
            0.0
        } else if n <= 10 {
            f64::from(outer.wrapping_sub(n_hs))
        } else {
            f64::from(outer.wrapping_sub(n_hs)) / f64::from(n.wrapping_sub(outer).wrapping_sub(1))
        };
        Ok(inverse_sqrt(delta))
    })
}

/// `detail::nVals`: `1/sqrt(nOuterElecs - nHs)` per atom.
fn n_vals(mol: &mut ROMol, force: bool) -> Result<Vec<f64>> {
    let tbl = PeriodicTable::get();
    cached_atom_values(mol, cp::_connectivityNVals, force, |mol, atom, idx| {
        let outer = tbl.get_n_outer_elecs(atom.atomic_num())? as u32;
        let n_hs = mol.total_num_hs(idx, false)?;
        Ok(inverse_sqrt(f64::from(outer.wrapping_sub(n_hs))))
    })
}

/// The path-based Chi sum: for every atom path of `n + 1` atoms the product of
/// the values of its atoms (a ring closure counts its first atom once).
fn chi_over_paths(mol: &ROMol, values: &[f64], n: u32) -> Result<f64> {
    let mut res = 0.0;
    for path in find_all_paths_of_length_n(mol, n + 1, false, false, -1, false, None)? {
        let n = n as usize;
        debug_assert_eq!(path.len(), n + 1);
        let mut accum = 1.0;
        for &atom in &path[..n] {
            accum *= values[atom as usize];
        }
        // the last atom only counts if this is not a ring (github 463)
        if path[n] != path[0] {
            accum *= values[path[n] as usize];
        }
        res += accum;
    }
    Ok(res)
}

/// Sum over the bonds of the product of the values of their atoms (Chi1).
fn chi1(mol: &ROMol, values: &[f64]) -> f64 {
    let mut res = 0.0;
    for bond in mol.bonds() {
        // `res += a * b` is a fused multiply-add on the reference platform
        res = fmuladd(
            values[bond.begin_atom_idx() as usize],
            values[bond.end_atom_idx() as usize],
            res,
        );
    }
    res
}

/// `calcChiNv`
pub fn calc_chi_nv(mol: &mut ROMol, n: u32, force: bool) -> Result<f64> {
    let deltas = hk_deltas(mol, force)?;
    chi_over_paths(mol, &deltas, n)
}

/// `calcChiNn`: like [`calc_chi_nv`] with the number of outer electrons instead
/// of the valence, which matters beyond the first row.
pub fn calc_chi_nn(mol: &mut ROMol, n: u32, force: bool) -> Result<f64> {
    let values = n_vals(mol, force)?;
    chi_over_paths(mol, &values, n)
}

/// `calcChi0v`
pub fn calc_chi0v(mol: &mut ROMol, force: bool) -> Result<f64> {
    Ok(hk_deltas(mol, force)?.iter().fold(0.0, |acc, v| acc + v))
}

/// `calcChi1v`
pub fn calc_chi1v(mol: &mut ROMol, force: bool) -> Result<f64> {
    let deltas = hk_deltas(mol, force)?;
    Ok(chi1(mol, &deltas))
}

/// `calcChi0n`
pub fn calc_chi0n(mol: &mut ROMol, force: bool) -> Result<f64> {
    Ok(n_vals(mol, force)?.iter().fold(0.0, |acc, v| acc + v))
}

/// `calcChi1n`
pub fn calc_chi1n(mol: &mut ROMol, force: bool) -> Result<f64> {
    let values = n_vals(mol, force)?;
    Ok(chi1(mol, &values))
}
