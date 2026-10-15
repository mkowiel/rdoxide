//! MMFF formal and partial charges (upstream
//! `MMFFMolProperties::computeMMFFCharges`, after Halgren, MMFF.V,
//! J. Comput. Chem. 1996, 17, 616-641).

use super::props::MmffMolProperties;
use crate::error::Result;
use crate::graphmol::romol::ROMol;

/// Share of the charge on `n` terminal atoms of which one keeps the double
/// bond: `-(n - 1) / n`.
fn shared_charge(n: i32) -> f64 {
    -(f64::from(n - 1) / f64::from(n))
}

impl MmffMolProperties {
    /// Sets the formal and partial charges of all atoms.
    pub(super) fn compute_charges(&mut self, mol: &ROMol) -> Result<()> {
        self.require_valid()?;
        for idx in 0..mol.num_atoms() {
            let charge = self.formal_charge_of(mol, idx);
            self.atoms[idx as usize].formal_charge = charge;
        }
        let params = self.params();
        for idx in 0..mol.num_atoms() {
            let atom_type = u32::from(self.atom_type(idx));
            let mut q0 = self.formal_charge(idx);
            let m = f64::from(self.prop(idx).crd);
            let v = params
                .pbci(atom_type)
                .expect("MMFF atom types always have PBCI parameters")
                .fcadj;
            let (mut sum_formal, mut sum_partial) = (0.0, 0.0);
            if v.abs() < 1.0e-10 {
                // a negative formal charge on a neighbor influences ipso
                for nbr in mol.atom_neighbors(idx) {
                    let nbr_charge = self.formal_charge(nbr);
                    if nbr_charge < 0.0 {
                        q0 += nbr_charge / (2.0 * f64::from(mol.atom_degree(nbr)));
                    }
                }
            }
            // anionic divalent nitrogen with a positively charged neighbor
            if atom_type == 62 {
                for nbr in mol.atom_neighbors(idx) {
                    let nbr_charge = self.formal_charge(nbr);
                    if nbr_charge > 0.0 {
                        q0 -= nbr_charge / 2.0;
                    }
                }
            }
            for nbr in mol.atom_neighbors(idx) {
                let bond = mol
                    .bond_idx_between(idx, nbr)
                    .expect("neighbors are bonded");
                let nbr_type = u32::from(self.atom_type(nbr));
                let bond_type = self.bond_type(mol, bond)?;
                // the sign of the increment depends on the order of the
                // atom types in the table
                let (sign, bci) = params.chg(bond_type, atom_type, nbr_type);
                sum_partial += match bci {
                    Some(bci) => f64::from(sign) * bci,
                    None => {
                        params.pbci(atom_type).expect("PBCI").pbci
                            - params.pbci(nbr_type).expect("PBCI").pbci
                    }
                };
                sum_formal += self.formal_charge(nbr);
            }
            // equation 15 of MMFF.V
            self.atoms[idx as usize].partial_charge =
                (1.0 - m * v) * q0 + v * sum_formal + sum_partial;
        }
        Ok(())
    }

    /// The formal charge MMFF assigns to an atom from its type and
    /// surroundings.
    fn formal_charge_of(&self, mol: &ROMol, idx: u32) -> f64 {
        let atom_type = self.atom_type(idx);
        match atom_type {
            // O2CM, O3N, OP, O2P, ..., SM: the charge of a (thio)carboxylate,
            // nitrate, phosphate, sulfonate or perchlorate group is shared
            // among its terminal atoms
            32 | 72 => self.terminal_group_charge(mol, idx),
            // N5M
            76 => {
                let ring = mol
                    .ring_info()
                    .atom_rings()
                    .iter()
                    .find(|r| r.contains(&(idx as i32)));
                match ring {
                    Some(ring) => {
                        let n = ring
                            .iter()
                            .filter(|&&a| self.atom_type(a as u32) == 76)
                            .count();
                        if n > 0 { -(1.0 / n as f64) } else { 0.0 }
                    }
                    None => 0.0,
                }
            }
            // NCN+, NGD+, NIM+/N5A+/N5B+/N5+
            55 | 56 | 81 => self.conjugated_cation_charge(mol, idx),
            // diazonium: +1 on the secondary nitrogen
            61 => {
                if mol.atom_neighbors(idx).any(|n| self.atom_type(n) == 42) {
                    1.0
                } else {
                    0.0
                }
            }
            34 | 49 | 51 | 54 | 58 | 92 | 93 | 94 | 97 => 1.0,
            87 | 95 | 96 | 98 | 99 => 2.0,
            88 => 3.0,
            35 | 62 | 89 | 90 | 91 => -1.0,
            _ => 0.0,
        }
    }

    fn terminal_group_charge(&self, mol: &ROMol, idx: u32) -> f64 {
        for nbr in mol.atom_neighbors(idx) {
            let nbr_type = self.atom_type(nbr);
            let z = mol.atom(nbr).atomic_num();
            let mut n_sec_n = 0;
            let mut n_term_os = 0;
            for nbr2 in mol.atom_neighbors(nbr) {
                let a2 = mol.atom(nbr2);
                if a2.atomic_num() == 7 && mol.atom_degree(nbr2) == 2 && !a2.is_aromatic() {
                    n_sec_n += 1;
                }
                if matches!(a2.atomic_num(), 8 | 16) && mol.atom_degree(nbr2) == 1 {
                    n_term_os += 1;
                }
            }
            // a deprotonated sulfonamide: the nitrogen does not replace an
            // oxygen of a sulfone
            if z == 16 && n_term_os == 2 && n_sec_n == 1 {
                n_sec_n = 0;
            }
            if z == 6 && n_term_os != 0 {
                return if n_term_os == 1 {
                    -1.0
                } else {
                    shared_charge(n_term_os)
                };
            }
            if nbr_type == 45 && n_term_os == 3 {
                return -1.0 / 3.0;
            }
            if nbr_type == 25 && n_term_os != 0 {
                return if n_term_os == 1 {
                    0.0
                } else {
                    shared_charge(n_term_os)
                };
            }
            if nbr_type == 18 && n_term_os != 0 {
                return if n_sec_n + n_term_os == 2 {
                    0.0
                } else {
                    -(f64::from(n_sec_n + n_term_os - 2) / f64::from(n_term_os))
                };
            }
            if nbr_type == 73 && n_term_os != 0 {
                return if n_term_os == 1 {
                    0.0
                } else {
                    shared_charge(n_term_os)
                };
            }
            if nbr_type == 77 && n_term_os != 0 {
                return -(1.0 / f64::from(n_term_os));
            }
        }
        0.0
    }

    /// The charge of a nitrogen in a cationic conjugated system: the formal
    /// charge of all conjugated nitrogens of type 81, 55 or 56 (linked through
    /// carbons of type 57 or 80) divided equally among them.
    fn conjugated_cation_charge(&self, mol: &ROMol, idx: u32) -> f64 {
        let mut charge = f64::from(mol.atom(idx).formal_charge());
        let mut conjugated = vec![false; mol.num_atoms() as usize];
        conjugated[idx as usize] = true;
        let mut n_conj = 1;
        let mut old_n_conj = 0;
        while n_conj > old_n_conj {
            old_n_conj = n_conj;
            for i in 0..mol.num_atoms() {
                if !conjugated[i as usize] {
                    continue;
                }
                for nbr in mol.atom_neighbors(i) {
                    if !matches!(self.atom_type(nbr), 57 | 80) {
                        continue;
                    }
                    for nbr2 in mol.atom_neighbors(nbr) {
                        if !matches!(self.atom_type(nbr2), 55 | 56 | 81) {
                            continue;
                        }
                        if !conjugated[nbr2 as usize] {
                            conjugated[nbr2 as usize] = true;
                            charge += f64::from(mol.atom(nbr2).formal_charge());
                            n_conj += 1;
                        }
                    }
                }
            }
        }
        charge / f64::from(n_conj)
    }
}
