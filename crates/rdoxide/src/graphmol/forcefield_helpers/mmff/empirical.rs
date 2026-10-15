//! Empirical rules for the MMFF parameters that are missing from the tables
//! (upstream `AtomTyper.cpp`, MMFF.V pages 625-628 and MMFF.IV page 609).

use super::props::MmffMolProperties;
use super::rings::{is_angle_in_ring_of_size_3_or_4, periodic_table_row, periodic_table_row_hl};
use crate::error::{RDKitError, Result};
use crate::forcefield::mmff::params::{Angle, Bond, Prop, Tor};
use crate::graphmol::bond::BondType;
use crate::graphmol::romol::ROMol;

const DEG2RAD: f64 = std::f64::consts::PI / 180.0;

impl MmffMolProperties {
    /// Bond stretching parameters from the covalent radii and Pauling
    /// electronegativities (equations 18 and 19, with Badger's rule in the
    /// Herschbach-Laurie form when no reference bond exists).
    pub fn bond_stretch_empirical_rule_params(&self, mol: &ROMol, bond: u32) -> Result<Bond> {
        self.require_valid()?;
        let params = self.params();
        let bond = mol.bond(bond);
        let z1 = mol.atom(bond.begin_atom_idx()).atomic_num();
        let z2 = mol.atom(bond.end_atom_idx()).atomic_num();
        let cov = |z: u32, n: u32| {
            params.cov_rad_pau_ele(z).ok_or_else(|| {
                RDKitError::precondition(format!(
                    "covalent radius/Pauling electronegativity parameters for atom {n} not found"
                ))
            })
        };
        let (cov1, cov2) = (cov(z1, 1)?, cov(z2, 2)?);
        let c = if z1 == 1 || z2 == 1 { 0.050 } else { 0.085 };
        let n = 1.4;
        // equation (18); the delta correction of upstream is zero
        let r0 = cov1.r0 + cov2.r0 - c * (cov1.chi - cov2.chi).abs().powf(n);
        let kb = if let Some(bndk) = params.bndk(z1, z2) {
            // equation (19)
            let coeff = bndk.r0 / r0;
            let coeff2 = coeff * coeff;
            bndk.kb * (coeff2 * coeff2 * coeff2)
        } else {
            let hl = params
                .herschbach_laurie(periodic_table_row_hl(z1), periodic_table_row_hl(z2))
                .ok_or_else(|| {
                    RDKitError::precondition("Herschbach-Laurie parameters not found")
                })?;
            10f64.powf(-(r0 - hl.a_ij) / hl.d_ij)
        };
        Ok(Bond { kb, r0 })
    }

    /// Torsion parameters from rules (a)-(h) of MMFF.IV for the central
    /// bond idx2-idx3.
    pub fn torsion_empirical_rule_params(&self, mol: &ROMol, idx2: u32, idx3: u32) -> Result<Tor> {
        self.require_valid()?;
        let params = self.params();
        let (j_type, k_type) = (self.atom_type(idx2), self.atom_type(idx3));
        let (j, k) = (self.prop(idx2), self.prop(idx3));
        let bond = mol
            .bond_idx_between(idx2, idx3)
            .map(|b| mol.bond(b))
            .ok_or_else(|| RDKitError::precondition("invalid bond"))?;
        let atomic_num = [mol.atom(idx2).atomic_num(), mol.atom(idx3).atomic_num()];
        let (mut u, mut v, mut w) = ([0.0f64; 2], [0.0f64; 2], [0.0f64; 2]);
        for i in 0..2 {
            let (ui, vi, wi) = match atomic_num[i] {
                6 => (2.0, 2.12, 0.0),
                7 => (2.0, 1.5, 0.0),
                8 => (2.0, 0.2, 2.0),
                14 => (1.25, 1.22, 0.0),
                15 => (1.25, 2.40, 0.0),
                16 => (1.25, 0.49, 8.0),
                _ => (0.0, 0.0, 0.0),
            };
            u[i] = ui;
            v[i] = vi;
            w[i] = wi;
        }
        let n_jk = f64::from((i32::from(j.crd) - 1) * (i32::from(k.crd) - 1));
        let mut tor = Tor {
            v1: 0.0,
            v2: 0.0,
            v3: 0.0,
        };
        let v2_of = |beta: f64, pi_jk: f64| beta * pi_jk * (u[0] * u[1]).sqrt();
        let v3_of = || (v[0] * v[1]).sqrt() / n_jk;
        // rules (e) and (f): an sp3 atom next to a planar or linear one
        let flat_neighbor = |p: &Prop| {
            (p.crd == 3 && (p.val == 4 || p.val == 34 || p.mltb != 0))
                || (p.crd == 2 && (p.val == 3 || p.mltb != 0))
        };
        if j.linh != 0 || k.linh != 0 {
            // rule (a): all terms stay 0
        } else if params.is_aromatic_type(j_type)
            && params.is_aromatic_type(k_type)
            && bond.is_aromatic()
        {
            // rule (b)
            let beta = if (j.val == 3 && k.val == 4) || (j.val == 4 && k.val == 3) {
                3.0
            } else {
                6.0
            };
            let pi_jk = if j.pilp == 0 && k.pilp == 0 { 0.5 } else { 0.3 };
            tor.v2 = v2_of(beta, pi_jk);
        } else if bond.bond_type() == BondType::Double {
            // rule (c)
            let pi_jk = if j.mltb == 2 && k.mltb == 2 { 1.0 } else { 0.4 };
            tor.v2 = v2_of(6.0, pi_jk);
        } else if j.crd == 4 && k.crd == 4 {
            // rule (d)
            tor.v3 = v3_of();
        } else if j.crd == 4 && k.crd != 4 {
            // rule (e)
            if !flat_neighbor(k) {
                tor.v3 = v3_of();
            }
        } else if k.crd == 4 && j.crd != 4 {
            // rule (f)
            if !flat_neighbor(j) {
                tor.v3 = v3_of();
            }
        } else if (bond.bond_type() == BondType::Single && j.mltb != 0 && k.mltb != 0)
            || (j.mltb != 0 && k.pilp != 0)
            || (j.pilp != 0 && k.mltb != 0)
        {
            // rule (g)
            let both_second_row =
                periodic_table_row(atomic_num[0]) == 2 && periodic_table_row(atomic_num[1]) == 2;
            let pi_from_row = if both_second_row { 0.3 } else { 0.15 };
            if j.pilp != 0 && k.pilp != 0 {
                // case (1): all terms stay 0
            } else if j.pilp != 0 && k.mltb != 0 {
                let pi_jk = if j.mltb == 1 { 0.5 } else { pi_from_row };
                tor.v2 = v2_of(6.0, pi_jk);
            } else if k.pilp != 0 && j.mltb != 0 {
                let pi_jk = if k.mltb == 1 { 0.5 } else { pi_from_row };
                tor.v2 = v2_of(6.0, pi_jk);
            } else if (j.mltb == 1 || k.mltb == 1) && (atomic_num[0] != 6 || atomic_num[1] != 6) {
                tor.v2 = v2_of(6.0, 0.4);
            } else {
                tor.v2 = v2_of(6.0, 0.15);
            }
        } else if matches!(atomic_num[0], 8 | 16) && matches!(atomic_num[1], 8 | 16) {
            // rule (h)
            tor.v2 = -(w[0] * w[1]).sqrt();
        } else {
            tor.v3 = v3_of();
        }
        Ok(tor)
    }
}

/// Angle bending parameters from the covalent radii (equation 20). The rest
/// value of the angle is taken from `old_params` when the table has one.
pub fn angle_bend_empirical_rule_params(
    mol: &ROMol,
    old_params: Option<&Angle>,
    central: &Prop,
    bond1: &Bond,
    bond2: &Bond,
    idx1: u32,
    idx2: u32,
    idx3: u32,
) -> Angle {
    let atomic_num = [idx1, idx2, idx3].map(|i| mol.atom(i).atomic_num());
    let ring_size = is_angle_in_ring_of_size_3_or_4(mol, idx1, idx2, idx3);
    let theta0 = match old_params {
        Some(old) => old.theta0,
        None => {
            let mut theta0 = 120.0;
            match central.crd {
                4 => theta0 = 109.45,
                2 => {
                    if atomic_num[1] == 8 {
                        theta0 = 105.0;
                    } else if central.linh == 1 {
                        theta0 = 180.0;
                    }
                }
                3 => {
                    if central.val == 3 && central.mltb == 0 {
                        theta0 = if atomic_num[1] == 7 { 107.0 } else { 92.0 };
                    }
                }
                _ => {}
            }
            match ring_size {
                3 => 60.0,
                4 => 90.0,
                _ => theta0,
            }
        }
    };
    let (mut z, mut c) = ([0.0f64; 3], [0.0f64; 3]);
    for i in 0..3 {
        // Table VI - MMFF.V, page 628
        let (zi, ci) = match atomic_num[i] {
            1 => (1.395, 0.0),
            6 => (2.494, 1.016),
            7 => (2.711, 1.113),
            8 => (3.045, 1.337),
            9 => (2.847, 0.0),
            14 => (2.350, 0.811),
            15 => (2.350, 1.068),
            16 => (2.980, 1.249),
            17 => (2.909, 1.078),
            35 => (3.017, 0.0),
            53 => (3.086, 0.0),
            _ => (0.0, 0.0),
        };
        z[i] = zi;
        c[i] = ci;
    }
    let mut beta = 1.75;
    let (r0_ij, r0_jk) = (bond1.r0, bond2.r0);
    let d = (r0_ij - r0_jk) * (r0_ij - r0_jk) / ((r0_ij + r0_jk) * (r0_ij + r0_jk));
    let theta0_rad = DEG2RAD * theta0;
    match ring_size {
        4 => beta *= 0.85,
        3 => beta *= 0.05,
        _ => {}
    }
    // equation (20)
    let ka =
        beta * z[0] * c[1] * z[2] / ((r0_ij + r0_jk) * theta0_rad * theta0_rad * (2.0 * d).exp());
    Angle { ka, theta0 }
}
