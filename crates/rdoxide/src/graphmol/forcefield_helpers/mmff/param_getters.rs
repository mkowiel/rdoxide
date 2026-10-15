//! Lookup of the MMFF parameters of the terms of a molecule, falling back to
//! the empirical rules (upstream `MMFFMolProperties::getMMFF*Params`).

use super::empirical::angle_bend_empirical_rule_params;
use super::props::{MmffMolProperties, stretch_bend_type};
use super::rings::periodic_table_row;
use crate::error::Result;
use crate::forcefield::mmff::params::{Angle, Bond, Oop, Stbn, Tor};
use crate::graphmol::romol::ROMol;

fn is_double_zero(x: f64) -> bool {
    x < 1.0e-10 && x > -1.0e-10
}

/// Parameters of a bond stretching term.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BondStretchParams {
    pub bond_type: u32,
    pub params: Bond,
}

/// Parameters of an angle bending term.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AngleBendParams {
    pub angle_type: u32,
    pub params: Angle,
}

/// Parameters of a stretch-bend term, with those of its bonds and angle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StretchBendParams {
    pub stretch_bend_type: u32,
    pub params: Stbn,
    pub bonds: [Bond; 2],
    pub angle: Angle,
}

/// Parameters of a torsion term.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TorsionParams {
    pub torsion_type: u32,
    pub params: Tor,
}

impl MmffMolProperties {
    /// Stretching parameters of the bond idx1-idx2 (`None` if there is no bond
    /// or the properties are invalid).
    pub fn bond_stretch_params(
        &self,
        mol: &ROMol,
        idx1: u32,
        idx2: u32,
    ) -> Result<Option<BondStretchParams>> {
        if !self.valid {
            return Ok(None);
        }
        let Some(bond) = mol.bond_idx_between(idx1, idx2) else {
            return Ok(None);
        };
        let bond_type = self.bond_type(mol, bond)?;
        let tabulated = self.params().bond(
            bond_type,
            u32::from(self.atom_type(idx1)),
            u32::from(self.atom_type(idx2)),
        );
        let params = match tabulated {
            Some(p) => *p,
            None => self.bond_stretch_empirical_rule_params(mol, bond)?,
        };
        Ok(Some(BondStretchParams { bond_type, params }))
    }

    /// Bending parameters of the angle idx1-idx2-idx3.
    pub fn angle_bend_params(
        &self,
        mol: &ROMol,
        idx1: u32,
        idx2: u32,
        idx3: u32,
    ) -> Result<Option<AngleBendParams>> {
        if !(self.valid
            && mol.bond_idx_between(idx1, idx2).is_some()
            && mol.bond_idx_between(idx2, idx3).is_some())
        {
            return Ok(None);
        }
        let angle_type = self.angle_type(mol, idx1, idx2, idx3)?;
        let types = [idx1, idx2, idx3].map(|i| u32::from(self.atom_type(i)));
        let tabulated = self
            .params()
            .angle(angle_type, types[0], types[1], types[2])
            .copied();
        let mut result = tabulated;
        // a missing force constant is derived from the bond stretching
        // parameters of the two bonds
        if tabulated.is_none_or(|a| is_double_zero(a.ka)) {
            if let (Some(b1), Some(b2)) = (
                self.bond_stretch_params(mol, idx1, idx2)?,
                self.bond_stretch_params(mol, idx2, idx3)?,
            ) {
                result = Some(angle_bend_empirical_rule_params(
                    mol,
                    tabulated.as_ref(),
                    self.prop(idx2),
                    &b1.params,
                    &b2.params,
                    idx1,
                    idx2,
                    idx3,
                ));
            }
        }
        Ok(result.map(|params| AngleBendParams { angle_type, params }))
    }

    /// Stretch-bend parameters of the angle idx1-idx2-idx3 (none for linear
    /// central atoms or when all coefficients vanish).
    pub fn stretch_bend_params(
        &self,
        mol: &ROMol,
        idx1: u32,
        idx2: u32,
        idx3: u32,
    ) -> Result<Option<StretchBendParams>> {
        if !self.valid || self.prop(idx2).linh != 0 {
            return Ok(None);
        }
        let types = [idx1, idx2, idx3].map(|i| u32::from(self.atom_type(i)));
        let (Some(b1), Some(b2)) = (
            self.bond_stretch_params(mol, idx1, idx2)?,
            self.bond_stretch_params(mol, idx2, idx3)?,
        ) else {
            return Ok(None);
        };
        let Some(angle) = self.angle_bend_params(mol, idx1, idx2, idx3)? else {
            return Ok(None);
        };
        let (bt1, bt2) = (b1.bond_type, b2.bond_type);
        let stretch_bend_type = stretch_bend_type(
            angle.angle_type,
            if types[0] <= types[2] { bt1 } else { bt2 },
            if types[0] < types[2] { bt2 } else { bt1 },
        );
        let params = self.params();
        let (mut swap, mut stbn) =
            params.stbn(stretch_bend_type, bt1, bt2, types[0], types[1], types[2]);
        if stbn.is_none() {
            let rows = [idx1, idx2, idx3].map(|i| periodic_table_row(mol.atom(i).atomic_num()));
            (swap, stbn) = params.dfsb(rows[0], rows[1], rows[2]);
        }
        let Some(stbn) = stbn else {
            return Ok(None);
        };
        if is_double_zero(stbn.kba_ijk) && is_double_zero(stbn.kba_kji) {
            return Ok(None);
        }
        let stbn = if swap {
            Stbn {
                kba_ijk: stbn.kba_kji,
                kba_kji: stbn.kba_ijk,
            }
        } else {
            *stbn
        };
        Ok(Some(StretchBendParams {
            stretch_bend_type,
            params: stbn,
            bonds: [b1.params, b2.params],
            angle: angle.params,
        }))
    }

    /// Torsion parameters of idx1-idx2-idx3-idx4 (none if all the terms
    /// vanish).
    pub fn torsion_params(
        &self,
        mol: &ROMol,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
    ) -> Result<Option<TorsionParams>> {
        if !(self.valid
            && mol.bond_idx_between(idx1, idx2).is_some()
            && mol.bond_idx_between(idx2, idx3).is_some()
            && mol.bond_idx_between(idx3, idx4).is_some())
        {
            return Ok(None);
        }
        let types = [idx1, idx2, idx3, idx4].map(|i| u32::from(self.atom_type(i)));
        let tor_types = self.torsion_type(mol, idx1, idx2, idx3, idx4)?;
        let (matched_type, tabulated) = self
            .params()
            .tor(tor_types, types[0], types[1], types[2], types[3]);
        let (torsion_type, params) = match tabulated {
            Some(p) => (
                if matched_type != 0 {
                    matched_type
                } else {
                    tor_types.0
                },
                *p,
            ),
            None => (
                tor_types.0,
                self.torsion_empirical_rule_params(mol, idx2, idx3)?,
            ),
        };
        if is_double_zero(params.v1) && is_double_zero(params.v2) && is_double_zero(params.v3) {
            return Ok(None);
        }
        Ok(Some(TorsionParams {
            torsion_type,
            params,
        }))
    }

    /// Out-of-plane bending parameters for the center idx2 with the
    /// outer atoms idx1, idx3, idx4 (none if the table has no row).
    pub fn oop_bend_params(
        &self,
        mol: &ROMol,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
    ) -> Option<Oop> {
        if !(self.valid
            && mol.bond_idx_between(idx1, idx2).is_some()
            && mol.bond_idx_between(idx2, idx3).is_some()
            && mol.bond_idx_between(idx2, idx4).is_some())
        {
            return None;
        }
        let types = [idx1, idx2, idx3, idx4].map(|i| u32::from(self.atom_type(i)));
        self.params()
            .oop(types[0], types[1], types[2], types[3])
            .copied()
    }
}
