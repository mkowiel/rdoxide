//! Neutralization of molecules by removing protons from the most acidic
//! atoms.
//!
//! Ported from RDKit `Code/GraphMol/StructChecker/ReCharge.cpp` (2026.09.1).
//! The charge tables this needs have no loader upstream, so with the
//! default options nothing is ever changed; the upstream quirks (marked
//! below) are kept.
//!
//! Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::atom_symbol::atom_symbol_match;
use super::options::StructCheckerOptions;
use super::pattern::{aa_match, convert_bond_type, ligand_matches};
use super::stereo::{EVEN_PARITY, ODD_PARITY, atom_parity, cis_trans_perception};
use super::utilities::{get_mol_atom_points, setup_neighbourhood};
use super::{AABondType, PathEntry};
use crate::Result;
use crate::graphmol::romol::ROMol;

/// `TotalCharge`. Upstream sums the charges of the atoms but returns 0.
pub fn total_charge(mol: &ROMol) -> i32 {
    let _charge: i32 = mol.atoms().iter().map(|a| a.formal_charge()).sum();
    0
}

/// `AllCentersRefined`: are all the defined stereocentres unambiguously
/// numbered, i.e. do all their ligands get different numbers?
fn all_centers_refined(mol: &ROMol, numbering: &[u32]) -> Result<bool> {
    let neighbour_array = setup_neighbourhood(mol);
    for (i, nbp) in neighbour_array.iter().enumerate() {
        let element = mol.atom(i as u32).atomic_num();
        let n_ligands = nbp.bonds.len();
        if matches!(element, 6 | 16 | 7 | 8 | 15 | 14) && n_ligands > 2 && n_ligands <= 4 {
            let parity = atom_parity(mol, i as u32, nbp)?;
            if parity == EVEN_PARITY || parity == ODD_PARITY {
                let mut refnum: Vec<u32> =
                    nbp.atoms.iter().map(|&a| numbering[a as usize]).collect();
                refnum.sort_unstable();
                if refnum.windows(2).any(|w| w[0] == w[1]) {
                    return Ok(false);
                }
            }
        }
    }
    Ok(true)
}

#[derive(Debug, Clone, Copy, Default)]
struct AtomRank {
    atom: u32,
    rank: u32,
    n_ligands: u32,
    elneg: f64,
    rank_sum: u32,
}

/// Insertion sort with the swap condition `less(a[j], a[j-1])`.
fn insertion_sort(ranks: &mut [AtomRank], less: impl Fn(&AtomRank, &AtomRank) -> bool) {
    for i in 1..ranks.len() {
        for j in (1..=i).rev() {
            if less(&ranks[j], &ranks[j - 1]) {
                ranks.swap(j, j - 1);
            } else {
                break;
            }
        }
    }
}

fn by_atom(a: &AtomRank, b: &AtomRank) -> bool {
    a.atom < b.atom
}

/// Removes protons until the desired charge is reached.
#[derive(Debug)]
pub struct ChargeFix<'a> {
    options: &'a StructCheckerOptions,
    mol: &'a mut ROMol,
    bond_color: Vec<u32>,
    atom_color: Vec<u32>,
    atom_pka_value: Vec<f64>,
    atom_old_pka_value: Vec<f64>,
}

impl<'a> ChargeFix<'a> {
    pub fn new(options: &'a StructCheckerOptions, mol: &'a mut ROMol) -> Self {
        let mut c = ChargeFix {
            options,
            mol,
            bond_color: Vec::new(),
            atom_color: Vec::new(),
            atom_pka_value: Vec::new(),
            atom_old_pka_value: Vec::new(),
        };
        c.reset_colors();
        c.reset_values();
        c
    }

    fn reset_colors(&mut self) {
        self.bond_color = vec![0; self.mol.num_bonds() as usize];
        self.atom_color = vec![0; self.mol.num_atoms() as usize];
    }

    fn reset_values(&mut self) {
        self.atom_pka_value = vec![0.0; self.mol.num_atoms() as usize];
    }

    /// `rechargeMolecule`: removes hydrogens until the desired charge is
    /// reached. The positions are selected by "acidity" combined with a
    /// refinement algorithm. Returns whether the molecule could be
    /// neutralized, the number of deprotonations and of refinement cycles.
    pub fn recharge_molecule(&mut self) -> Result<(bool, u32, u32)> {
        let mut gap = 0.0;
        let pka_value = 0.0;
        let mut acidity_limit = 24.0;
        let mut numbering = vec![0u32; self.mol.num_atoms() as usize];
        let mut ndeprot = 0;
        let mut nrefine = 0;
        loop {
            while total_charge(self.mol) > self.options.desired_charge {
                self.set_pka_values()?;
                // (upstream passes the arguments in swapped order)
                let nacid = self.mark_most_acidic_atoms(&mut gap, &mut acidity_limit);
                if pka_value > acidity_limit {
                    return Ok((false, ndeprot, nrefine)); // not acidic enough
                } else if nacid == 1 || (nacid == total_charge(self.mol) && gap > 8.0) {
                    self.decrement_marked_charges();
                } else {
                    break;
                }
                ndeprot += 1;
            }

            if total_charge(self.mol) > self.options.desired_charge {
                nrefine += 1;
                let mut nacid = self.refine_acidic_atoms(&mut numbering)?;
                if nacid == 1 {
                    self.decrement_marked_charges();
                } else if all_centers_refined(self.mol, &numbering)? {
                    for i in 0..self.atom_color.len() {
                        // select one mark
                        if self.atom_color[i] != 0 && nacid != 1 {
                            nacid -= 1;
                            self.atom_color[i] = 0;
                        }
                    }
                    self.decrement_marked_charges();
                } else {
                    return Ok((false, ndeprot, nrefine));
                }
                ndeprot += 1;
            } else {
                self.reset_values();
                return Ok((true, ndeprot, nrefine));
            }
        }
    }

    /// `setpKaValues`: estimates the pKa values of the acidic atoms.
    /// Returns true if all went well.
    fn set_pka_values(&mut self) -> Result<bool> {
        let opts = self.options;
        let mut result = true;
        let na = self.mol.num_atoms() as usize;
        let neighbour_array = setup_neighbourhood(self.mol);
        self.reset_colors();
        self.atom_old_pka_value = self.atom_pka_value.clone();
        self.reset_values();
        let mol: &ROMol = self.mol;

        for i in 0..na {
            let ap = mol.atom(i as u32);
            // is the atom acidic?
            let acidic = opts
                .acidic_atoms
                .iter()
                .any(|aa| aa_match(mol, i as u32, aa, &[], &neighbour_array, opts.verbose));
            if !acidic {
                continue; // center not acidic -> next atom
            }
            self.atom_color[i] += 1;
            self.atom_pka_value[i] = 0.0;
            let mut pka = 0.0;

            // add the local charge increment
            if let Some(j) = opts
                .charge_inc_table
                .iter()
                .find(|j| atom_symbol_match(&ap.symbol(), &j.atom_symbol))
            {
                pka += f64::from(ap.formal_charge()) * j.local_inc;
            }
            // add the local atom acidity
            if let Some(j) = opts
                .atom_acidity
                .iter()
                .find(|j| atom_symbol_match(&ap.symbol(), &j.atom_symbol))
            {
                pka += j.local_inc;
            }

            let nbp = &neighbour_array[i];
            // all alpha neighbours
            for j in 0..nbp.atoms.len() {
                let aap = mol.atom(nbp.atoms[j]);
                let abp = mol.bond(nbp.bonds[j]);
                // the alpha path conductivity
                let ppa: Option<&PathEntry> = opts.alpha_path_table.iter().find(|k| {
                    atom_symbol_match(&ap.symbol(), &k.path.atom_symbol)
                        && k.path
                            .ligands
                            .first()
                            .is_some_and(|l| ligand_matches(aap, abp, l, false))
                });
                let Some(ppa) = ppa else {
                    result = false;
                    break;
                };

                if convert_bond_type(abp.bond_type()) != AABondType::Single {
                    // (upstream indexes the table with the neighbour number)
                    let mult_inc = opts.charge_inc_table.get(j).map_or(0.0, |e| e.mult_inc);
                    pka += f64::from(
                        convert_bond_type(abp.bond_type()) as i32 - AABondType::Single as i32,
                    ) * mult_inc;
                }

                // the alpha charge increment
                if aap.formal_charge() != 0 {
                    let Some(k) = opts
                        .charge_inc_table
                        .iter()
                        .find(|k| atom_symbol_match(&aap.symbol(), &k.atom_symbol))
                    else {
                        result = false;
                        break;
                    };
                    pka += f64::from(aap.formal_charge()) * k.alpha_inc;
                }

                // the alpha acidity increment
                let Some(k) = opts
                    .atom_acidity
                    .iter()
                    .find(|k| atom_symbol_match(&aap.symbol(), &k.atom_symbol))
                else {
                    result = false;
                    break;
                };
                pka += ppa.cond * k.alpha_inc;

                // all beta neighbours
                let nbph = &neighbour_array[nbp.atoms[j] as usize];
                for k in 0..nbph.atoms.len() {
                    if nbph.atoms[k] as usize == i {
                        continue; // no loop back
                    }
                    let bap = mol.atom(nbph.atoms[k]);
                    let bbp = mol.bond(nbph.bonds[k]);
                    // the beta conductivity
                    // (upstream dereferences a null pointer to match here)
                    let ppb = opts.beta_path_table.iter().find(|l| {
                        atom_symbol_match(&ap.symbol(), &l.path.atom_symbol)
                            && l.path
                                .ligands
                                .first()
                                .is_some_and(|x| ligand_matches(aap, abp, x, false))
                            && l.path
                                .ligands
                                .get(1)
                                .is_some_and(|x| ligand_matches(bap, bbp, x, false))
                    });
                    let Some(ppb) = ppb else {
                        result = false;
                        break;
                    };
                    // the beta acidity increment
                    let Some(l) = opts
                        .atom_acidity
                        .iter()
                        .find(|l| atom_symbol_match(&bap.symbol(), &l.atom_symbol))
                    else {
                        result = false;
                        break;
                    };
                    pka += ppb.cond * l.beta_inc;
                    // the beta charge increment
                    if bap.formal_charge() != 0 {
                        let Some(l) = opts
                            .charge_inc_table
                            .iter()
                            .find(|l| atom_symbol_match(&bap.symbol(), &l.atom_symbol))
                        else {
                            result = false;
                            break;
                        };
                        pka += f64::from(bap.formal_charge()) * l.beta_inc;
                    }
                }
            }
            // pKa = 7 + (pKa'-7)*beta + ((pKa' - 7)*alpha) ^ 3
            self.atom_pka_value[i] = 7.0
                + opts.beta * (pka - 7.0)
                + opts.alpha * (pka - 7.0) * opts.alpha * (pka - 7.0) * opts.alpha * (pka - 7.0);
        }
        Ok(result)
    }

    /// `markMostAcidicAtoms`: marks the acidic atoms with the smallest pKa,
    /// removes the other marks. `pka_value` is set to the minimum pKa and
    /// `gap` to the pKa difference to the next acidic set of atoms.
    fn mark_most_acidic_atoms(&mut self, pka_value: &mut f64, gap: &mut f64) -> i32 {
        let mut min_pka = 1000.0;
        let mut next_pka = 1000.0;
        let mut result = 0;
        let epsilon = 0.000001;
        for i in 0..self.atom_color.len() {
            if self.atom_color[i] != 0 && self.atom_pka_value[i] < min_pka {
                min_pka = self.atom_pka_value[i];
            }
        }
        for i in 0..self.atom_color.len() {
            if self.atom_color[i] != 0 && self.atom_pka_value[i] < min_pka + epsilon {
                result += 1;
                self.atom_color[i] = 1;
            } else {
                if self.atom_color[i] != 0 && self.atom_pka_value[i] < next_pka {
                    next_pka = self.atom_pka_value[i];
                }
                self.atom_color[i] = 0;
                self.atom_pka_value[i] = 0.0;
            }
        }
        *pka_value = min_pka;
        *gap = next_pka - min_pka;
        result
    }

    /// `decrementMarkedCharges`
    fn decrement_marked_charges(&mut self) {
        for i in 0..self.mol.num_atoms() {
            if self.atom_color[i as usize] != 0 {
                let atom = self.mol.atom_mut(i);
                atom.set_formal_charge(atom.formal_charge() - 1);
            }
        }
    }

    /// `refineAcidicAtoms`: refines the class of acidic atoms to a smaller
    /// one based on the electronegativity of the neighbouring atoms.
    /// Returns the size of this class (-1 if an element has no
    /// electronegativity).
    fn refine_acidic_atoms(&mut self, numbering: &mut [u32]) -> Result<i32> {
        let opts = self.options;
        let epsilon = 0.000001;
        let n = self.mol.num_atoms() as usize;
        let mut atom_ranks = vec![AtomRank::default(); n];

        // set the electronegativities
        for (i, ar) in atom_ranks.iter_mut().enumerate() {
            ar.atom = i as u32 + 1;
            ar.rank = 0;
            ar.rank_sum = 0;
            let atom = self.mol.atom(i as u32);
            match opts.elneg_table.get(&atom.atomic_num()) {
                Some(e) => {
                    ar.elneg = e + 3.0 * f64::from(atom.formal_charge()) - opts.elneg0;
                }
                None => return Ok(-1),
            }
        }
        // the preliminary ranking by decreasing electronegativity
        insertion_sort(&mut atom_ranks, |a, b| a.elneg > b.elneg + epsilon);
        if n > 0 {
            atom_ranks[0].rank = 0;
        }
        let mut j = 0;
        for i in 1..n {
            if atom_ranks[i].elneg < atom_ranks[i - 1].elneg - epsilon {
                j = i as u32;
            }
            atom_ranks[i].rank = j;
        }
        insertion_sort(&mut atom_ranks, by_atom); // resort by atom number

        // use unsaturation to split the ranks (rank_sum is misused here)
        for ar in &mut atom_ranks {
            ar.rank_sum = 0;
        }
        for b in 0..self.mol.num_bonds() {
            let bond = self.mol.bond(b);
            if convert_bond_type(bond.bond_type()) != AABondType::Single {
                atom_ranks[bond.begin_atom_idx() as usize].rank_sum += 1;
                atom_ranks[bond.end_atom_idx() as usize].rank_sum += 1;
            }
        }
        // sort by rank (ascending) and unsaturation (descending)
        insertion_sort(&mut atom_ranks, |a, b| {
            a.rank < b.rank || (a.rank == b.rank && a.rank_sum > b.rank_sum)
        });
        let mut j = 0;
        for i in 1..n {
            if atom_ranks[i].rank > atom_ranks[i - 1].rank
                || atom_ranks[i].rank_sum < atom_ranks[i - 1].rank_sum
            {
                j = i as u32;
            }
            atom_ranks[i].rank = j;
        }
        insertion_sort(&mut atom_ranks, by_atom); // restore the atom order

        for ar in &mut atom_ranks {
            ar.n_ligands = 0;
        }
        for b in 0..self.mol.num_bonds() {
            let bond = self.mol.bond(b);
            atom_ranks[bond.begin_atom_idx() as usize].n_ligands += 1;
            atom_ranks[bond.end_atom_idx() as usize].n_ligands += 1;
        }

        // refine the ranking using the rank sums of the neighbours
        let mut do_cis_trans = false;
        loop {
            for i in 0..n {
                numbering[i] = atom_ranks[i].rank;
            }
            for ar in &mut atom_ranks {
                ar.rank_sum = 0;
            }
            for b in 0..self.mol.num_bonds() {
                let bond = self.mol.bond(b);
                let (bi, ei) = (bond.begin_atom_idx() as usize, bond.end_atom_idx() as usize);
                atom_ranks[bi].rank_sum += atom_ranks[ei].rank;
                atom_ranks[ei].rank_sum += atom_ranks[bi].rank;
            }
            if do_cis_trans {
                let (atom_point, non_zero_z) = get_mol_atom_points(self.mol, false)?;
                if non_zero_z {
                    cis_trans_perception(self.mol, &atom_point, numbering, &mut self.bond_color);
                }
                for (i, &color) in self.bond_color.iter().enumerate() {
                    if color != 0 {
                        let bond = self.mol.bond(i as u32);
                        atom_ranks[bond.begin_atom_idx() as usize].rank_sum += color;
                        atom_ranks[bond.end_atom_idx() as usize].rank_sum += color;
                    }
                }
            }
            for ar in &mut atom_ranks {
                // use the average rank sum
                // (shifting the decimal point)
                if let Some(v) = (ar.rank_sum * 10).checked_div(ar.n_ligands) {
                    ar.rank_sum = v;
                }
            }
            // sort by rank and rank sum
            insertion_sort(&mut atom_ranks, |a, b| {
                a.rank < b.rank || (a.rank == b.rank && a.rank_sum < b.rank_sum)
            });
            // the new ranks
            let mut changed = false;
            let mut j = 0;
            for i in 1..n {
                if atom_ranks[i].rank > atom_ranks[i - 1].rank
                    || atom_ranks[i].rank_sum > atom_ranks[i - 1].rank_sum
                {
                    j = i as u32;
                }
                changed = changed || atom_ranks[i].rank != j;
                atom_ranks[i].rank = j;
            }
            insertion_sort(&mut atom_ranks, by_atom); // restore the atom order
            if !changed && !do_cis_trans {
                do_cis_trans = true;
                changed = true;
            } else {
                do_cis_trans = false;
            }
            if !changed {
                break;
            }
        }

        // the smallest rank of the coloured atoms
        let mut min_rank = n as u32;
        for i in 0..n {
            if self.atom_color[i] != 0 && atom_ranks[i].rank < min_rank {
                min_rank = atom_ranks[i].rank;
            }
        }
        let mut result = 0;
        for i in 0..n {
            if self.atom_color[i] != 0 && atom_ranks[i].rank == min_rank {
                // count the members of the minimum class
                result += 1;
                self.atom_pka_value[i] = f64::from(atom_ranks[i].rank) + 1.0;
            } else {
                // uncolour the non-minimum members
                self.atom_color[i] = 0;
                self.atom_pka_value[i] = 0.0;
            }
        }
        Ok(result)
    }
}
