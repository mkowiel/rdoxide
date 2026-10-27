//! `denormal.c`: turns AROMATIC bonds into single and double bonds (kekulization with hydrogen
//! ranges). Copyright 2001-2011 Novartis Pharma AG.

use super::graph::{BondSet, combine_rings, ring_list};
use super::perceive::{AROMATIC, DOUBLE, Neighbourhood, SINGLE, TRIPLE, setup_neighbourhood};
use super::reaccs::{DOUBLET, Molecule, NONE};
use super::set::BitSet;
use super::utilities::{cas_valence, implicit_hydrogens};

/// `ANY_BOND`, `MINUS_TAUTOMER`, `PLUS_TAUTOMER`
const ANY_BOND: i32 = 8;
const PLUS_TAUTOMER: i32 = 0x100;
const MINUS_TAUTOMER: i32 = 0x200;

const MULTIPLE_COMPOSITION: i32 = 1;
const CONSTRAINT_VIOLATION: i32 = 2;
const CANDIDATE: u8 = 1;
const TOUCHED: u8 = 2;
const VIOLATION: i32 = -1;

/// `struct atom_constraint_t`
#[derive(Debug, Clone, Copy, Default)]
pub struct AtomConstraint {
    pub dbl_min: i32,
    pub dbl_max: i32,
    pub is_oxigen: bool,
    pub open: u8,
}

/// `struct bond_constraint_t`
#[derive(Debug, Clone, Copy, Default)]
pub struct BondConstraint {
    pub can_be_single: bool,
    pub can_be_double: bool,
    pub scount: i32,
    pub open: u8,
}

/// `ColorNormalizedComponents`: atoms connected by AROMATIC bonds share a colour
fn color_normalized_components(mp: &mut Molecule) {
    for (i, a) in mp.atoms.iter_mut().enumerate() {
        a.color = i as i32 + 1;
    }
    loop {
        let mut nmerge = 0;
        for b in &mp.bonds {
            if b.bond_type != AROMATIC {
                continue;
            }
            let (i0, i1) = ((b.atoms[0] - 1) as usize, (b.atoms[1] - 1) as usize);
            let (c0, c1) = (mp.atoms[i0].color, mp.atoms[i1].color);
            if c0 != c1 {
                nmerge += 1;
                if c0 < c1 {
                    mp.atoms[i1].color = c0;
                } else {
                    mp.atoms[i0].color = c1;
                }
            }
        }
        if nmerge == 0 {
            break;
        }
    }
}

/// `CCTDenormalize`: the result has `MULTIPLE_COMPOSITION` / `CONSTRAINT_VIOLATION` bits
pub fn cct_denormalize(mp: &mut Molecule) -> i32 {
    if mp.bonds.is_empty() {
        return 0;
    }
    let nbp = setup_neighbourhood(mp).unwrap_or_default();
    let mut acp = vec![AtomConstraint::default(); mp.atoms.len()];
    let mut bcp = vec![BondConstraint::default(); mp.bonds.len()];
    if nbp.len() != mp.atoms.len() {
        return 0;
    }

    // 'AROMATIC' bonds between shortcuts become NONE bonds
    for i in 0..mp.bonds.len() {
        if mp.bonds[i].bond_type != ANY_BOND {
            continue;
        }
        for k in 0..2 {
            if !mp.atoms[(mp.bonds[i].atoms[k] - 1) as usize]
                .atext
                .is_empty()
            {
                mp.bonds[i].bond_type = NONE;
            }
        }
    }
    for (i, ap) in mp.atoms.iter().enumerate() {
        let mut valence = if ap.dummy1 == 0 {
            cas_valence(&ap.atom_symbol)
        } else {
            ap.dummy1
        };
        let (hco_min, hco_max);
        let (dmin, dmax) = (ap.dummy2 & 0xF, (ap.dummy2 >> 4) & 0xF);
        if dmin >= dmax {
            hco_min = dmin;
            hco_max = dmin;
        } else if valence > 4 && (ap.atom_symbol == "S" || ap.atom_symbol == "I") {
            hco_min = dmin;
            hco_max = dmin;
        } else {
            hco_min = dmin;
            hco_max = dmax;
        }
        let (mut hco_min, mut hco_max) = (hco_min, hco_max);
        let (mut nsingle, mut ndouble, mut ntriple, mut naromatic) = (0, 0, 0, 0);
        for &bj in &nbp[i].bonds {
            match mp.bonds[bj].bond_type {
                ANY_BOND => {}
                SINGLE => nsingle += 1,
                DOUBLE => ndouble += 1,
                TRIPLE => ntriple += 1,
                _ => naromatic += 1,
            }
        }
        let mut charge = ap.charge;
        if ap.dummy2 & PLUS_TAUTOMER != 0 {
            charge += 1;
        }
        if ap.dummy2 & MINUS_TAUTOMER != 0 {
            charge -= 1;
        }
        if valence == 0 {
            valence = (nsingle * 2 + ndouble * 4 + ntriple * 6 + naromatic * 3) / 2;
            valence += hco_max;
            valence += charge.abs();
        } else if ap.dummy1 == 0 && ap.dummy2 == 0 {
            hco_min = implicit_hydrogens(
                &ap.atom_symbol,
                nsingle,
                naromatic,
                ndouble,
                ntriple,
                ap.radical == DOUBLET,
                charge,
            );
            hco_max = hco_min;
            if nsingle == 0
                && ndouble == 0
                && ntriple == 0
                && naromatic == 2
                && hco_max == 0
                && ap.atom_symbol == "N"
            {
                hco_max = 1;
            }
        }
        let nbonds = nbp[i].atoms.len() as i32;
        acp[i].open = 0;
        acp[i].is_oxigen = ap.atom_symbol == "O";
        acp[i].dbl_min = valence - charge.abs() - hco_max - nbonds - 2 * ntriple;
        acp[i].dbl_max = valence - charge.abs() - hco_min - nbonds - 2 * ntriple;
        if acp[i].dbl_min < 0 {
            acp[i].dbl_max -= acp[i].dbl_min;
            acp[i].dbl_min = 0;
        }
    }
    for (i, bp) in mp.bonds.iter().enumerate() {
        bcp[i].open = 0;
        bcp[i].can_be_single = bp.bond_type == SINGLE || bp.bond_type == AROMATIC;
        let (a0, a1) = ((bp.atoms[0] - 1) as usize, (bp.atoms[1] - 1) as usize);
        bcp[i].can_be_double = if acp[a0].dbl_max == 0 || acp[a1].dbl_max == 0 {
            false
        } else {
            bp.bond_type == DOUBLE || bp.bond_type == AROMATIC
        };
        bcp[i].scount = nbp[a0].atoms.len() as i32 + nbp[a1].atoms.len() as i32 - 2;
    }

    let result = denormalize(mp, &mut acp, &mut bcp);

    for (i, ap) in mp.atoms.iter_mut().enumerate() {
        ap.dummy2 &= PLUS_TAUTOMER | MINUS_TAUTOMER;
        if ap.dummy2 != NONE {
            let total_order: i32 = nbp[i].bonds.iter().map(|&b| mp.bonds[b].bond_type).sum();
            if total_order == 4 && ap.dummy2 & PLUS_TAUTOMER != 0 && ap.atom_symbol == "N" {
                ap.charge += 1;
            }
            if total_order == 1 && ap.dummy2 & MINUS_TAUTOMER != 0 && ap.atom_symbol == "O" {
                ap.charge -= 1;
            }
        }
        ap.color = NONE;
    }
    result
}

/// `Refine`: propagates the constraints; returns the number of undecided bonds or `VIOLATION`
fn refine(
    acp: &mut [AtomConstraint],
    bcp: &mut [BondConstraint],
    mp: &Molecule,
    nbp: &[Neighbourhood],
) -> i32 {
    loop {
        let mut changed = false;
        for i in 0..mp.atoms.len() {
            if acp[i].open & CANDIDATE != 0 && acp[i].open & TOUCHED != 0 {
                acp[i].open = CANDIDATE;
                let (mut dbl_min, mut dbl_max) = (0, 0);
                for &bj in &nbp[i].bonds {
                    let b = bcp[bj];
                    if b.can_be_single && b.can_be_double {
                        dbl_max += 1;
                    } else if b.can_be_single {
                        // single only
                    } else if b.can_be_double {
                        dbl_min += 1;
                        dbl_max += 1;
                    } else {
                        return VIOLATION;
                    }
                }
                if dbl_min == acp[i].dbl_max {
                    for (&bj, &aj) in nbp[i].bonds.iter().zip(&nbp[i].atoms) {
                        if bcp[bj].can_be_single && bcp[bj].can_be_double {
                            changed = true;
                            bcp[bj].can_be_double = false;
                            bcp[bj].open |= TOUCHED;
                            acp[aj].open |= TOUCHED;
                        }
                    }
                } else if dbl_max == acp[i].dbl_min {
                    for (&bj, &aj) in nbp[i].bonds.iter().zip(&nbp[i].atoms) {
                        if bcp[bj].can_be_single && bcp[bj].can_be_double {
                            changed = true;
                            bcp[bj].can_be_single = false;
                            bcp[bj].open |= TOUCHED;
                            acp[aj].open |= TOUCHED;
                        }
                    }
                } else if dbl_max < acp[i].dbl_min || dbl_min > acp[i].dbl_max {
                    return VIOLATION;
                }
            }
        }
        for b in bcp.iter_mut() {
            if b.open & CANDIDATE != 0 && b.open & TOUCHED != 0 {
                b.open = CANDIDATE;
            }
        }
        if !changed {
            break;
        }
    }
    bcp.iter()
        .filter(|b| b.open & CANDIDATE != 0 && b.can_be_single == b.can_be_double)
        .count() as i32
}

/// `AddNewSolution`: a solution with a known composition replaces the old one
fn add_new_solution(new_sol: BondSet, mut old: Vec<BondSet>) -> Vec<BondSet> {
    if let Some(b) = old
        .iter_mut()
        .find(|b| b.cardinality == new_sol.cardinality)
    {
        b.bond_set = new_sol.bond_set;
    } else {
        old.insert(0, new_sol);
    }
    old
}

/// `RecDenormalize`
fn rec_denormalize(
    acp: &mut [AtomConstraint],
    bcp: &mut [BondConstraint],
    mp: &Molecule,
    nbp: &[Neighbourhood],
    prev: Vec<BondSet>,
) -> Vec<BondSet> {
    let nopen = refine(acp, bcp, mp, nbp);
    if nopen == VIOLATION {
        return prev;
    }
    if nopen == 0 {
        let mut sol = BondSet {
            cardinality: 0,
            bond_set: BitSet::new(mp.bonds.len()),
        };
        for (i, b) in bcp.iter().enumerate() {
            if b.open & CANDIDATE != 0 && b.can_be_double {
                sol.bond_set.put(i);
                sol.cardinality += 1;
            }
        }
        return add_new_solution(sol, prev);
    }
    let mut ntauto = -1;
    let mut ibranch: Option<usize> = None;
    for i in 0..mp.bonds.len() {
        if bcp[i].open & CANDIDATE != 0 && bcp[i].can_be_single && bcp[i].can_be_double {
            let a0 = (mp.bonds[i].atoms[0] - 1) as usize;
            let a1 = (mp.bonds[i].atoms[1] - 1) as usize;
            let mut itauto = i32::from(acp[a0].dbl_min < acp[a0].dbl_max);
            if acp[a1].dbl_min < acp[a1].dbl_max {
                itauto += 1;
            }
            if itauto > ntauto {
                ntauto = itauto;
                ibranch = Some(i);
            }
        }
    }
    let Some(ib) = ibranch else { return prev };
    let (a0, a1) = (
        (mp.bonds[ib].atoms[0] - 1) as usize,
        (mp.bonds[ib].atoms[1] - 1) as usize,
    );
    let mut acp_new = acp.to_vec();
    let mut bcp_new = bcp.to_vec();
    bcp_new[ib].can_be_single = false;
    acp_new[a0].open |= TOUCHED;
    acp_new[a1].open |= TOUCHED;
    let sol = rec_denormalize(&mut acp_new, &mut bcp_new, mp, nbp, prev);
    acp_new.copy_from_slice(acp);
    bcp_new.copy_from_slice(bcp);
    bcp_new[ib].can_be_double = false;
    acp_new[a0].open |= TOUCHED;
    acp_new[a1].open |= TOUCHED;
    rec_denormalize(&mut acp_new, &mut bcp_new, mp, nbp, sol)
}

/// `Denormalize`: converts the AROMATIC bonds one connected component at a time
pub fn denormalize(
    mp: &mut Molecule,
    acp: &mut [AtomConstraint],
    bcp: &mut [BondConstraint],
) -> i32 {
    let nbp = setup_neighbourhood(mp).unwrap_or_default();
    color_normalized_components(mp);
    // MoleculeRings: the rings are only needed by the (dummy) solution comparison
    if !mp.bonds.is_empty() {
        let graph: Vec<[u32; 2]> = mp
            .bonds
            .iter()
            .map(|b| [b.atoms[0] as u32, b.atoms[1] as u32])
            .collect();
        let mut rings = ring_list(&graph);
        combine_rings(&mut rings);
    }
    let mut result = 0;
    while let Some(first) = mp.bonds.iter().find(|b| b.bond_type == AROMATIC) {
        let col = mp.atoms[(first.atoms[0] - 1) as usize].color;
        for (i, ap) in mp.atoms.iter().enumerate() {
            acp[i].open = if ap.color == col {
                CANDIDATE | TOUCHED
            } else {
                0
            };
        }
        for (i, bp) in mp.bonds.iter().enumerate() {
            bcp[i].open =
                if bp.bond_type == AROMATIC && mp.atoms[(bp.atoms[0] - 1) as usize].color == col {
                    CANDIDATE
                } else {
                    0
                };
        }
        let solution = rec_denormalize(acp, bcp, mp, &nbp, Vec::new());
        if solution.is_empty() {
            result |= CONSTRAINT_VIOLATION;
            break;
        }
        if solution.len() > 1 {
            result |= MULTIPLE_COMPOSITION;
        }
        for i in 0..mp.bonds.len() {
            if bcp[i].open & CANDIDATE != 0 {
                mp.bonds[i].bond_type = if solution[0].bond_set.is_member(i) {
                    DOUBLE
                } else {
                    SINGLE
                };
            }
        }
    }
    result
}
