//! Surface-area descriptors (upstream `GraphMol/Descriptors/MolSurf.{h,cpp}`):
//! Labute's approximate surface area, the topological polar surface area (TPSA)
//! and the VSA binnings built on them.
//!
//! Like upstream, the per-atom results are cached on the molecule
//! (`_labuteAtomContribs`, `_tpsaAtomContribs-<includeSandP>`, ...) and
//! served from there unless `force` is set. The floating-point expressions
//! keep upstream's operation order and fused multiply-adds (see
//! [`crate::rdgeneral::fp`]) so the sums are bit-identical.

use std::f64::consts::PI;

use crate::error::{RDKitError, Result};
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::BondType;
use crate::graphmol::descriptors::crippen::get_crippen_atom_contribs;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::dict::{RDProps, RDValue};
use crate::rdgeneral::fp::fmuladd;

pub const LABUTE_ASA_VERSION: &str = "1.0.2";
pub const TPSA_VERSION: &str = "2.0.0";

/// Per-atom Labute contributions with the contribution of the implicit Hs.
#[derive(Clone, Debug, PartialEq)]
pub struct LabuteContribs {
    pub atoms: Vec<f64>,
    pub h: f64,
    pub total: f64,
}

/// `Vi[i] += Rj*Rj - (Ri - dij)^2 / dij`, the overlap term shared by the bond
/// and the implicit-H loops of `getLabuteAtomContribs`.
fn overlap(r_other: f64, r_this: f64, dij: f64) -> f64 {
    let d = r_this - dij;
    fmuladd(r_other, r_other, -(d * d / dij))
}

/// `std::min(std::max(fabs(Ri - Rj), bij), Ri + Rj)`
fn bond_distance(ri: f64, rj: f64, bij: f64) -> f64 {
    let m = if (ri - rj).abs() < bij {
        bij
    } else {
        (ri - rj).abs()
    };
    if ri + rj < m { ri + rj } else { m }
}

/// `getLabuteAtomContribs`
pub fn get_labute_atom_contribs(
    mol: &mut ROMol,
    include_hs: bool,
    force: bool,
) -> Result<LabuteContribs> {
    if !force
        && let (
            Some(RDValue::VecDouble(atoms)),
            Some(RDValue::Double(h)),
            Some(RDValue::Double(total)),
        ) = (
            mol.get_prop(cp::_labuteAtomContribs),
            mol.get_prop(cp::_labuteAtomHContrib),
            mol.get_prop(cp::_labuteASA),
        )
    {
        return Ok(LabuteContribs {
            atoms: atoms.clone(),
            h: *h,
            total: *total,
        });
    }
    let table = PeriodicTable::get();
    let n = mol.num_atoms() as usize;
    let rads = mol
        .atoms()
        .iter()
        .map(|a| table.get_rb0(a.atomic_num()))
        .collect::<Result<Vec<f64>>>()?;
    let mut vi = vec![0.0; n];

    const BOND_SCALE_FACTS: [f64; 4] = [0.1, 0.0, 0.2, 0.3];
    for bond in mol.bonds() {
        let (begin, end) = (bond.begin_atom_idx() as usize, bond.end_atom_idx() as usize);
        let (ri, rj) = (rads[begin], rads[end]);
        let mut bij = ri + rj;
        if !bond.is_aromatic() {
            let bond_type = bond.bond_type() as usize;
            if bond_type < 4 {
                bij -= BOND_SCALE_FACTS[bond_type];
            }
        } else {
            bij -= BOND_SCALE_FACTS[0];
        }
        let dij = bond_distance(ri, rj, bij);
        vi[begin] += overlap(rj, ri, dij);
        vi[end] += overlap(ri, rj, dij);
    }

    let mut h_contrib = 0.0;
    if include_hs {
        let rj = table.get_rb0(1)?;
        for (i, &ri) in rads.iter().enumerate() {
            let dij = bond_distance(ri, rj, ri + rj);
            vi[i] += overlap(rj, ri, dij);
            h_contrib += overlap(ri, rj, dij);
        }
    }

    let mut total = 0.0;
    for (v, &ri) in vi.iter_mut().zip(&rads) {
        *v = PI * ri * fmuladd(4.0, ri, -*v);
        total += *v;
    }
    if include_hs && h_contrib.abs() > 1e-4 {
        let rj = table.get_rb0(1)?;
        h_contrib = PI * rj * fmuladd(4.0, rj, -h_contrib);
        total += h_contrib;
    }
    mol.set_prop(cp::_labuteAtomContribs, vi.clone(), true)?;
    mol.set_prop(cp::_labuteAtomHContrib, h_contrib, true)?;
    mol.set_prop(cp::_labuteASA, total, true)?;
    Ok(LabuteContribs {
        atoms: vi,
        h: h_contrib,
        total,
    })
}

/// `calcLabuteASA`
pub fn calc_labute_asa(mol: &mut ROMol, include_hs: bool, force: bool) -> Result<f64> {
    if !force && let Some(RDValue::Double(total)) = mol.get_prop(cp::_labuteASA) {
        return Ok(*total);
    }
    Ok(get_labute_atom_contribs(mol, include_hs, force)?.total)
}

/// Neighbour and bond-type counts of one atom, hydrogens excluded
/// (`nNbrs`, `nHs`, `nSing`, `nDoub`, `nTrip`, `nArom` of `getTPSAAtomContribs`).
#[derive(Clone, Copy, Debug, Default)]
struct TpsaEnv {
    nbrs: i32,
    hs: i32,
    sing: i32,
    doub: i32,
    trip: i32,
    arom: i32,
}

/// Table value for a nitrogen, or `None` when the environment is not listed.
fn tpsa_nitrogen(e: TpsaEnv, chg: i32, in3_ring: bool) -> Option<f64> {
    let TpsaEnv {
        hs,
        sing,
        doub,
        trip,
        arom,
        ..
    } = e;
    match e.nbrs {
        1 => {
            if hs == 0 && chg == 0 && trip == 1 {
                Some(23.79)
            } else if hs == 1 && chg == 0 && doub == 1 {
                Some(23.85)
            } else if hs == 2 && chg == 0 && sing == 1 {
                Some(26.02)
            } else if hs == 2 && chg == 1 && doub == 1 {
                Some(25.59)
            } else if hs == 3 && chg == 1 && sing == 1 {
                Some(27.64)
            } else {
                None
            }
        }
        2 => {
            if hs == 0 && chg == 0 && sing == 1 && doub == 1 {
                Some(12.36)
            } else if hs == 0 && chg == 0 && trip == 1 && doub == 1 {
                Some(13.60)
            } else if hs == 1 && chg == 0 && sing == 2 && in3_ring {
                Some(21.94)
            } else if hs == 1 && chg == 0 && sing == 2 && !in3_ring {
                Some(12.03)
            } else if hs == 0 && chg == 1 && trip == 1 && sing == 1 {
                Some(4.36)
            } else if hs == 1 && chg == 1 && doub == 1 && sing == 1 {
                Some(13.97)
            } else if hs == 2 && chg == 1 && sing == 2 {
                Some(16.61)
            } else if hs == 0 && chg == 0 && arom == 2 {
                Some(12.89)
            } else if hs == 1 && chg == 0 && arom == 2 {
                Some(15.79)
            } else if hs == 1 && chg == 1 && arom == 2 {
                Some(14.14)
            } else {
                None
            }
        }
        3 => {
            if hs == 0 && chg == 0 && sing == 3 && in3_ring {
                Some(3.01)
            } else if hs == 0 && chg == 0 && sing == 3 && !in3_ring {
                Some(3.24)
            } else if hs == 0 && chg == 0 && sing == 1 && doub == 2 {
                Some(11.68)
            } else if hs == 0 && chg == 1 && sing == 2 && doub == 1 {
                Some(3.01)
            } else if hs == 1 && chg == 1 && sing == 3 {
                Some(4.44)
            } else if hs == 0 && chg == 0 && arom == 3 {
                Some(4.41)
            } else if hs == 0 && chg == 0 && sing == 1 && arom == 2 {
                Some(4.93)
            } else if hs == 0 && chg == 0 && doub == 1 && arom == 2 {
                Some(8.39)
            } else if hs == 0 && chg == 1 && arom == 3 {
                Some(4.10)
            } else if hs == 0 && chg == 1 && sing == 1 && arom == 2 {
                Some(3.88)
            } else {
                None
            }
        }
        4 => (hs == 0 && sing == 4 && chg == 1).then_some(0.0),
        _ => None,
    }
}

/// Table value for an oxygen, or `None` when the environment is not listed.
fn tpsa_oxygen(e: TpsaEnv, chg: i32, in3_ring: bool) -> Option<f64> {
    let TpsaEnv {
        hs,
        sing,
        doub,
        arom,
        ..
    } = e;
    match e.nbrs {
        1 => {
            if hs == 0 && chg == 0 && doub == 1 {
                Some(17.07)
            } else if hs == 1 && chg == 0 && sing == 1 {
                Some(20.23)
            } else if hs == 0 && chg == -1 && sing == 1 {
                Some(23.06)
            } else {
                None
            }
        }
        2 => {
            if hs == 0 && chg == 0 && sing == 2 && in3_ring {
                Some(12.53)
            } else if hs == 0 && chg == 0 && sing == 2 && !in3_ring {
                Some(9.23)
            } else if hs == 0 && chg == 0 && arom == 2 {
                Some(13.14)
            } else {
                None
            }
        }
        _ => None,
    }
}

/// Phosphorus (only with `includeSandP`): unlisted environments contribute 0.
fn tpsa_phosphorus(e: TpsaEnv, chg: i32) -> f64 {
    let TpsaEnv { hs, sing, doub, .. } = e;
    match e.nbrs {
        2 if hs == 0 && chg == 0 && sing == 1 && doub == 1 => 34.14,
        3 if hs == 0 && chg == 0 && sing == 3 => 13.59,
        3 if hs == 1 && chg == 0 && sing == 2 && doub == 1 => 23.47,
        4 if hs == 0 && chg == 0 && sing == 3 && doub == 1 => 9.81,
        _ => 0.0,
    }
}

/// Sulfur (only with `includeSandP`): unlisted environments contribute 0.
fn tpsa_sulfur(e: TpsaEnv, chg: i32) -> f64 {
    let TpsaEnv {
        hs,
        sing,
        doub,
        arom,
        ..
    } = e;
    match e.nbrs {
        1 if hs == 0 && chg == 0 && doub == 1 => 32.09,
        1 if hs == 1 && chg == 0 && sing == 1 => 38.80,
        2 if hs == 0 && chg == 0 && sing == 2 => 25.30,
        2 if hs == 0 && chg == 0 && arom == 2 => 28.24,
        3 if hs == 0 && chg == 0 && arom == 2 && doub == 1 => 21.70,
        3 if hs == 0 && chg == 0 && sing == 2 && doub == 1 => 19.21,
        4 if hs == 0 && chg == 0 && sing == 2 && doub == 2 => 8.38,
        _ => 0.0,
    }
}

/// `max(0, base - nNbrs * nbrWeight + nHs * 1.5)`, the estimate used for
/// nitrogens and oxygens that are not in the table.
fn tpsa_fallback(base: f64, nbr_weight: f64, e: TpsaEnv) -> f64 {
    let tmp = fmuladd(
        f64::from(e.hs),
        1.5,
        fmuladd(-f64::from(e.nbrs), nbr_weight, base),
    );
    if tmp < 0.0 { 0.0 } else { tmp }
}

/// The cache keys upstream builds with `boost::format("%s-%s")` (a bool is
/// printed as 0 or 1).
fn tpsa_keys(include_s_and_p: bool) -> (String, String) {
    let flag = u8::from(include_s_and_p);
    (
        format!("{}-{flag}", cp::_tpsa),
        format!("{}-{flag}", cp::_tpsaAtomContribs),
    )
}

fn tpsa_contribution(
    mol: &ROMol,
    i: u32,
    atom: &Atom,
    mut env: TpsaEnv,
    include_s_and_p: bool,
) -> Result<Option<f64>> {
    let atomic_num = atom.atomic_num();
    if atomic_num != 7
        && atomic_num != 8
        && (!include_s_and_p || (atomic_num != 15 && atomic_num != 16))
    {
        return Ok(None);
    }
    env.hs += i32::try_from(mol.total_num_hs(i, false)?).expect("H count fits in i32");
    let chg = atom.formal_charge();
    let in3_ring = mol.ring_info().is_atom_in_ring_of_size(i, 3)?;
    env.nbrs += i32::try_from(mol.atom_degree(i)).expect("degree fits in i32");
    Ok(Some(match atomic_num {
        7 => tpsa_nitrogen(env, chg, in3_ring).unwrap_or_else(|| tpsa_fallback(30.5, 8.2, env)),
        8 => tpsa_oxygen(env, chg, in3_ring).unwrap_or_else(|| tpsa_fallback(28.5, 8.6, env)),
        15 => tpsa_phosphorus(env, chg),
        _ => tpsa_sulfur(env, chg),
    }))
}

/// `getTPSAAtomContribs`: the contributions per atom and their sum.
pub fn get_tpsa_atom_contribs(
    mol: &mut ROMol,
    force: bool,
    include_s_and_p: bool,
) -> Result<(Vec<f64>, f64)> {
    let (total_key, contribs_key) = tpsa_keys(include_s_and_p);
    if !force && mol.has_prop(&contribs_key) {
        if let (Some(RDValue::VecDouble(contribs)), Some(RDValue::Double(total))) =
            (mol.get_prop(&contribs_key), mol.get_prop(&total_key))
        {
            return Ok((contribs.clone(), *total));
        }
        return Err(RDKitError::Key(format!(
            "cached TPSA properties {contribs_key}/{total_key} have the wrong type"
        )));
    }
    let n = mol.num_atoms() as usize;
    let mut envs = vec![TpsaEnv::default(); n];
    for bond in mol.bonds() {
        let (begin, end) = (bond.begin_atom_idx() as usize, bond.end_atom_idx() as usize);
        if mol.atom(begin as u32).atomic_num() == 1 {
            envs[end].nbrs -= 1;
            envs[end].hs += 1;
        } else if mol.atom(end as u32).atomic_num() == 1 {
            envs[begin].nbrs -= 1;
            envs[begin].hs += 1;
        } else if bond.is_aromatic() {
            envs[begin].arom += 1;
            envs[end].arom += 1;
        } else {
            let counter: fn(&mut TpsaEnv) -> &mut i32 = match bond.bond_type() {
                BondType::Single => |e| &mut e.sing,
                BondType::Double => |e| &mut e.doub,
                BondType::Triple => |e| &mut e.trip,
                _ => continue,
            };
            *counter(&mut envs[begin]) += 1;
            *counter(&mut envs[end]) += 1;
        }
    }

    let mut contribs = vec![0.0; n];
    let mut total = 0.0;
    for (i, env) in envs.into_iter().enumerate() {
        let atom = mol.atom(i as u32);
        if let Some(tmp) = tpsa_contribution(mol, i as u32, atom, env, include_s_and_p)? {
            contribs[i] = tmp;
            total += tmp;
        }
    }
    mol.set_prop(&contribs_key, contribs.clone(), true)?;
    mol.set_prop(&total_key, total, true)?;
    Ok((contribs, total))
}

/// `calcTPSA`
pub fn calc_tpsa(mol: &mut ROMol, force: bool, include_s_and_p: bool) -> Result<f64> {
    let (total_key, _) = tpsa_keys(include_s_and_p);
    if !force && let Some(RDValue::Double(total)) = mol.get_prop(&total_key) {
        return Ok(*total);
    }
    Ok(get_tpsa_atom_contribs(mol, force, include_s_and_p)?.1)
}

/// `assignContribsToBins`: adds each contribution to the bin its property
/// value falls in. The bin is found with the same halving search as
/// `std::upper_bound`, so unsorted user bins behave like upstream too.
fn assign_contribs_to_bins(
    contribs: &[f64],
    bin_prop: &[f64],
    bins: &[f64],
    res: &mut [f64],
) -> Result<()> {
    if contribs.len() != bin_prop.len() {
        return Err(RDKitError::precondition("mismatched array sizes"));
    }
    if res.len() < bins.len() + 1 {
        return Err(RDKitError::precondition("mismatched array sizes"));
    }
    for (&contrib, &value) in contribs.iter().zip(bin_prop) {
        let (mut first, mut len) = (0usize, bins.len());
        while len > 0 {
            let half = len >> 1;
            let middle = first + half;
            if value < bins[middle] {
                len = half;
            } else {
                first = middle + 1;
                len -= half + 1;
            }
        }
        res[first] += contrib;
    }
    Ok(())
}

const SLOGP_BINS: [f64; 11] = [-0.4, -0.2, 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6];
const SMR_BINS: [f64; 9] = [1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 3.8, 4.0];

/// Bins the Labute contributions by `property`, with `default_bins` when no
/// bins are given.
fn vsa_by_property(
    mol: &mut ROMol,
    bins: Option<&[f64]>,
    default_bins: &[f64],
    force: bool,
    property: impl FnOnce(&mut ROMol) -> Result<Vec<f64>>,
) -> Result<Vec<f64>> {
    let bins = bins.unwrap_or(default_bins);
    let mut res = vec![0.0; bins.len() + 1];
    let vsa = get_labute_atom_contribs(mol, true, force)?.atoms;
    let prop = property(mol)?;
    assign_contribs_to_bins(&vsa, &prop, bins, &mut res)?;
    Ok(res)
}

/// `calcSlogP_VSA`: Labute contributions binned by the Crippen logP contributions.
pub fn calc_slogp_vsa(mol: &mut ROMol, bins: Option<&[f64]>, force: bool) -> Result<Vec<f64>> {
    vsa_by_property(mol, bins, &SLOGP_BINS, force, |m| {
        Ok(get_crippen_atom_contribs(m, force)?.logp)
    })
}

/// `calcSMR_VSA`: Labute contributions binned by the Crippen MR contributions.
pub fn calc_smr_vsa(mol: &mut ROMol, bins: Option<&[f64]>, force: bool) -> Result<Vec<f64>> {
    vsa_by_property(mol, bins, &SMR_BINS, force, |m| {
        Ok(get_crippen_atom_contribs(m, force)?.mr)
    })
}

/// `calcCustomProp_VSA`: Labute contributions binned by a double atom property
/// (1 for atoms that do not have it).
pub fn calc_custom_prop_vsa(
    mol: &mut ROMol,
    custom_prop_name: &str,
    bins: &[f64],
    force: bool,
) -> Result<Vec<f64>> {
    vsa_by_property(mol, Some(bins), &[], force, |m| {
        m.atoms()
            .iter()
            .map(|atom| match atom.get_prop(custom_prop_name) {
                None => Ok(1.0),
                // `getProp<double>` casts doubles and floats, and parses strings
                Some(RDValue::String(s)) => s.trim().parse::<f64>().map_err(|_| {
                    RDKitError::Runtime(
                        "bad lexical cast: source type value could not be interpreted as target"
                            .into(),
                    )
                }),
                Some(v) => v
                    .as_double()
                    .ok_or_else(|| RDKitError::Runtime("bad any cast".into())),
            })
            .collect()
    })
}
