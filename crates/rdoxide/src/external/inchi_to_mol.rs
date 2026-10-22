//! The InChI -> molecule half of RDKit's InChI support
//! (`External/INCHI-API/inchi.cpp`: `InchiToMol` and its helpers
//! `assignBondDirs`, `findAlternatingBonds` and the `cleanUp` family that turns
//! the hypervalent structures the InChI library returns into the charge
//! separated forms RDKit accepts).
//!
//! Copyright (c) 2011-2025 Novartis Institutes for BioMedical Research Inc.
//! and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-ChemReactions-Novartis-BSD-3-Clause.txt, the same Novartis
//! license text).

use std::collections::{BTreeSet, VecDeque};

use crate::external::inchi::inchi_dll::{
    INCHI_BOND_STEREO_DOUBLE_EITHER, INCHI_BOND_STEREO_SINGLE_1DOWN,
    INCHI_BOND_STEREO_SINGLE_1EITHER, INCHI_BOND_STEREO_SINGLE_1UP, INCHI_BOND_STEREO_SINGLE_2DOWN,
    INCHI_BOND_STEREO_SINGLE_2UP, INCHI_BOND_TYPE_ALTERN, INCHI_BOND_TYPE_TRIPLE,
    ISOTOPIC_SHIFT_FLAG,
};
use crate::external::inchi::inchi_dll_c::{
    INCHI_RET_OKAY, INCHI_RET_WARNING, get_struct_from_inchi,
};
use crate::external::inchi::readinch::{
    INCHI_PARITY_EVEN, INCHI_PARITY_NONE, INCHI_PARITY_ODD, INCHI_PARITY_UNDEFINED,
    INCHI_PARITY_UNKNOWN, INCHI_STEREO_TYPE_ALLENE, INCHI_STEREO_TYPE_DOUBLE_BOND,
    INCHI_STEREO_TYPE_NONE, INCHI_STEREO_TYPE_TETRAHEDRAL,
};
use crate::external::inchi_api::ExtraInchiReturnValues;
use crate::graphmol::atom::{Atom, ChiralType};
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::chirality::legacy::{assign_atom_cip_ranks, assign_stereochemistry};
use crate::graphmol::mol_ops::{RemoveHsParameters, remove_hs, sanitize_mol};
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::{RDKitError, Result};

type IntPairVect = Vec<(u32, u32)>;

/// `assignBondDirs`: assign `ENDUPRIGHT` / `ENDDOWNRIGHT` to the bonds of the
/// double-bond stereo pairs; `false` when the constraints conflict.
fn assign_bond_dirs(mol: &mut RWMol, z_pairs: &IntPairVect, e_pairs: &IntPairVect) -> bool {
    let mut pending: BTreeSet<u32> = BTreeSet::new();
    for &(a, b) in z_pairs.iter().chain(e_pairs.iter()) {
        pending.insert(a);
        pending.insert(b);
    }
    let mut queue: VecDeque<(u32, BondDir)> = VecDeque::new();
    while !pending.is_empty() || !queue.is_empty() {
        if queue.is_empty() {
            let first = *pending.iter().next().unwrap();
            queue.push_back((first, BondDir::EndUpRight));
        } else {
            let (cur, dir) = queue.pop_front().unwrap();
            let existing = mol.bond(cur).bond_dir();
            if existing != BondDir::None {
                if existing != dir {
                    return false;
                }
            } else {
                mol.bond_mut(cur).set_bond_dir(dir);
                pending.remove(&cur);
                let other_dir = if dir == BondDir::EndUpRight {
                    BondDir::EndDownRight
                } else {
                    BondDir::EndUpRight
                };
                for k in 0..2 {
                    let (rules, d) = if k == 0 {
                        (z_pairs, dir)
                    } else {
                        (e_pairs, other_dir)
                    };
                    for &(first, second) in rules {
                        let other = if first == cur {
                            second
                        } else if second == cur {
                            first
                        } else {
                            continue;
                        };
                        if other == cur {
                            continue;
                        }
                        let od = mol.bond(other).bond_dir();
                        if od != BondDir::None {
                            if od != d {
                                return false;
                            }
                        } else {
                            queue.push_back((other, d));
                        }
                    }
                }
            }
        }
    }
    true
}

/// `findAlternatingBonds`: the shortest path of alternating bonds from
/// `current` to an atom of the wanted element and charge. The path (bond
/// indices) is returned in `path` with the first bond of the path at the
/// bottom and the last one on top (a `std::stack` in C++: `last()` is its
/// `top()`).
#[allow(clippy::too_many_arguments)]
fn find_alternating_bonds(
    mol: &ROMol,
    current: u32,
    atomic_num: u32,
    charge: i32,
    next_type: BondType,
    ending_type: BondType,
    current_len: u32,
    max_len: u32,
    last_bond: Option<u32>,
    path: &mut Vec<u32>,
    visited: &mut BTreeSet<u32>,
) -> Option<u32> {
    if last_bond.is_none() {
        visited.clear();
        path.clear();
    }
    visited.insert(current);

    if let Some(lb) = last_bond {
        let a = mol.atom(current);
        if a.atomic_num() == atomic_num
            && mol.bond(lb).bond_type() == ending_type
            && a.formal_charge() == charge
        {
            if path.is_empty() || path.len() > current_len as usize {
                path.clear();
                path.push(lb);
                return Some(current);
            }
            return None;
        }
    }

    if max_len <= current_len {
        return None;
    }

    let mut target = None;
    let nbrs: Vec<u32> = mol.atom_neighbors(current).collect();
    for nid in nbrs {
        if visited.contains(&nid) {
            continue;
        }
        let bidx = mol.bond_idx_between(current, nid).unwrap();
        let bt = mol.bond(bidx).bond_type();
        if bt == next_type {
            let nb = if next_type == BondType::Single {
                BondType::Double
            } else {
                BondType::Single
            };
            if let Some(t) = find_alternating_bonds(
                mol,
                nid,
                atomic_num,
                charge,
                nb,
                ending_type,
                current_len + 1,
                max_len,
                Some(bidx),
                path,
                visited,
            ) {
                target = Some(t);
            }
        } else if ending_type != BondType::Single
            && ending_type != BondType::Double
            && bt == ending_type
        {
            if let Some(t) = find_alternating_bonds(
                mol,
                nid,
                atomic_num,
                charge,
                BondType::Unspecified,
                ending_type,
                current_len + 1,
                0,
                Some(bidx),
                path,
                visited,
            ) {
                target = Some(t);
            }
        }
    }

    if target.is_some() {
        if let Some(lb) = last_bond {
            path.push(lb);
        }
        return target;
    }
    None
}

/// The common call shape of `findAlternatingBonds` from the clean-up helpers
/// (fresh path and visited set).
fn find_alt(
    mol: &ROMol,
    atom: u32,
    atomic_num: u32,
    charge: i32,
    next_type: BondType,
    ending_type: BondType,
    max_len: u32,
    path: &mut Vec<u32>,
) -> Option<u32> {
    let mut visited = BTreeSet::new();
    find_alternating_bonds(
        mol,
        atom,
        atomic_num,
        charge,
        next_type,
        ending_type,
        0,
        max_len,
        None,
        path,
        &mut visited,
    )
}

fn set_bt(mol: &mut RWMol, a: u32, b: u32, t: BondType) {
    let idx = mol.bond_idx_between(a, b).unwrap();
    mol.bond_mut(idx).set_bond_type(t);
}

fn flip_double_single(mol: &mut RWMol, path: &mut Vec<u32>) {
    while let Some(b) = path.pop() {
        let t = if mol.bond(b).bond_type() == BondType::Double {
            BondType::Single
        } else {
            BondType::Double
        };
        mol.bond_mut(b).set_bond_type(t);
    }
}

/// Build a query from `(atomic numbers, bonds)` as the helpers do.
fn make_query(atoms: &[u32], bonds: &[(u32, u32, BondType)]) -> Result<ROMol> {
    let mut q = ROMol::new();
    for &a in atoms {
        q.add_atom(Atom::new(a), false);
    }
    for &(a, b, t) in bonds {
        q.add_bond(a, b, t)?;
    }
    Ok(q)
}

/// The mapping query atom -> molecule atom when there is exactly one match.
fn unique_match(mol: &ROMol, q: &ROMol) -> Result<Option<Vec<u32>>> {
    let matches = substruct_match(mol, q, &SubstructMatchParameters::default())?;
    if matches.len() != 1 {
        return Ok(None);
    }
    let mut map = vec![0u32; q.num_atoms() as usize];
    for &(qi, mi) in &matches[0] {
        map[qi as usize] = mi as u32;
    }
    Ok(Some(map))
}

fn explicit_valence(mol: &mut RWMol, idx: u32) -> Result<i32> {
    mol.calc_explicit_valence(idx, false)
}

/// `_Valence4NCleanUp1`: clean C1=NN=[N-]=N1
fn valence4n_cleanup1(mol: &mut RWMol, atom: u32) -> Result<bool> {
    if mol.atom(atom).atomic_num() != 7
        || mol.atom(atom).formal_charge() != -1
        || explicit_valence(mol, atom)? != 4
    {
        return Ok(false);
    }
    mol.atom_mut(atom).set_atomic_num(50);
    mol.atom_mut(atom).set_formal_charge(0);
    let q = make_query(
        &[6, 7, 50, 7, 7],
        &[
            (0, 1, BondType::Single),
            (1, 2, BondType::Double),
            (2, 3, BondType::Double),
            (3, 4, BondType::Single),
            (4, 0, BondType::Double),
        ],
    )?;
    let Some(m) = unique_match(mol, &q)? else {
        mol.atom_mut(atom).set_atomic_num(7);
        mol.atom_mut(atom).set_formal_charge(-1);
        return Ok(false);
    };
    set_bt(mol, m[0], m[1], BondType::Double);
    set_bt(mol, m[1], m[2], BondType::Single);
    set_bt(mol, m[2], m[3], BondType::Single);
    set_bt(mol, m[3], m[4], BondType::Double);
    set_bt(mol, m[4], m[0], BondType::Single);
    mol.atom_mut(atom).set_atomic_num(7);
    mol.atom_mut(atom).set_formal_charge(-1);
    Ok(true)
}

/// `_Valence4NCleanUp2`: directly to a N via double bond
fn valence4n_cleanup2(mol: &mut RWMol, atom: u32) -> Result<bool> {
    let mut stack = Vec::new();
    let Some(target) = find_alt(
        mol,
        atom,
        7,
        0,
        BondType::Double,
        BondType::Double,
        1,
        &mut stack,
    ) else {
        return Ok(false);
    };
    let top = *stack.last().unwrap();
    mol.bond_mut(top).set_bond_type(BondType::Single);
    mol.atom_mut(atom).set_formal_charge(0);
    mol.atom_mut(target).set_formal_charge(-1);
    Ok(true)
}

/// `_Valence5NCleanUp1`: valence-5 N connected to a N+
fn valence5n_cleanup1(mol: &mut RWMol, atom: u32) -> Result<bool> {
    let mut stack = Vec::new();
    let Some(target) = find_alt(
        mol,
        atom,
        7,
        1,
        BondType::Double,
        BondType::Double,
        5,
        &mut stack,
    ) else {
        return Ok(false);
    };
    mol.atom_mut(target).set_formal_charge(0);
    explicit_valence(mol, target)?;
    flip_double_single(mol, &mut stack);
    mol.atom_mut(atom).set_formal_charge(1);
    Ok(true)
}

/// `_Valence5NCleanUp2`: N connected to N- through a triple then single bond
fn valence5n_cleanup2(mol: &mut RWMol, atom: u32) -> Result<bool> {
    let mut stack = Vec::new();
    let Some(target) = find_alt(
        mol,
        atom,
        7,
        -1,
        BondType::Triple,
        BondType::Single,
        2,
        &mut stack,
    ) else {
        return Ok(false);
    };
    let bond = stack.pop().unwrap();
    mol.bond_mut(bond).set_bond_type(BondType::Single);
    let (b, e) = (
        mol.bond(bond).begin_atom_idx(),
        mol.bond(bond).end_atom_idx(),
    );
    if b == atom {
        mol.atom_mut(e).set_formal_charge(-1);
    } else {
        mol.atom_mut(b).set_formal_charge(-1);
    }
    let top = *stack.last().unwrap();
    mol.bond_mut(top).set_bond_type(BondType::Double);
    mol.atom_mut(target).set_formal_charge(0);
    explicit_valence(mol, target)?;
    explicit_valence(mol, atom)?;
    Ok(true)
}

/// `_Valence5NCleanUp3`: directly to a N via double bond (github #1572)
fn valence5n_cleanup3(mol: &mut RWMol, atom: u32) -> Result<bool> {
    let mut stack = Vec::new();
    let Some(target) = find_alt(
        mol,
        atom,
        7,
        0,
        BondType::Double,
        BondType::Double,
        1,
        &mut stack,
    ) else {
        return Ok(false);
    };
    let mut stack2 = Vec::new();
    let target2 = find_alt(
        mol,
        atom,
        8,
        0,
        BondType::Double,
        BondType::Double,
        1,
        &mut stack2,
    );
    if target2.is_none() {
        mol.atom_mut(target).set_formal_charge(-1);
        explicit_valence(mol, target)?;
        let top = *stack.last().unwrap();
        mol.bond_mut(top).set_bond_type(BondType::Single);
        mol.atom_mut(atom).set_formal_charge(1);
        explicit_valence(mol, atom)?;
    }
    Ok(true)
}

/// `_Valence5NCleanUp4`: N connected to two Si- via double bonds
fn valence5n_cleanup4(mol: &mut RWMol, atom: u32) -> Result<bool> {
    let mut nbrs = Vec::new();
    let mut bonds = Vec::new();
    let mut n_si = 0;
    for nbr in mol.atom_neighbors(atom).collect::<Vec<_>>() {
        let bond = mol.bond_idx_between(nbr, atom).unwrap();
        if mol.atom(nbr).atomic_num() == 14
            && mol.atom(nbr).formal_charge() == -1
            && mol.bond(bond).bond_type() == BondType::Double
        {
            if n_si >= 2 {
                return Ok(false);
            }
            nbrs.push(nbr);
            bonds.push(bond);
            n_si += 1;
        }
    }
    if n_si != 2 {
        return Ok(false);
    }
    mol.atom_mut(nbrs[0]).set_formal_charge(0);
    mol.atom_mut(nbrs[1]).set_formal_charge(0);
    mol.bond_mut(bonds[0]).set_bond_type(BondType::Single);
    mol.bond_mut(bonds[1]).set_bond_type(BondType::Single);
    Ok(true)
}

/// `_Valence5NCleanUp5`: alternating bonds to O, S, F or Cl
fn valence5n_cleanup5(mol: &mut RWMol, atom: u32, atomic_num: u32) -> Result<bool> {
    assert!(
        matches!(atomic_num, 8 | 16 | 9 | 17),
        "this cleanup looks for O or S or Cl or F"
    );
    let mut stack_uncharged = Vec::new();
    let mut stack_charged = Vec::new();
    // the visited set is shared between the two searches (and only cleared
    // when a search starts without a previous bond, i.e. both times)
    let uncharged = find_alt(
        mol,
        atom,
        atomic_num,
        0,
        BondType::Double,
        BondType::Double,
        7,
        &mut stack_uncharged,
    );
    let charged = find_alt(
        mol,
        atom,
        atomic_num,
        1,
        BondType::Double,
        BondType::Double,
        7,
        &mut stack_charged,
    );
    if uncharged.is_none() && charged.is_none() {
        return Ok(false);
    }
    let use_charged = uncharged.is_none();
    if let (Some(_), Some(c)) = (uncharged, charged) {
        assert!(mol.atom(c).formal_charge() == 1, "expecting +1 charge");
        mol.atom_mut(c).set_formal_charge(0);
        mol.atom_mut(c).set_num_explicit_hs(0);
    }
    mol.atom_mut(atom).set_formal_charge(1);
    let stack = if use_charged {
        &mut stack_charged
    } else {
        &mut stack_uncharged
    };
    flip_double_single(mol, stack);
    match (uncharged, charged) {
        (Some(u), Some(_)) => mol.atom_mut(u).set_num_explicit_hs(1),
        (Some(u), None) => mol.atom_mut(u).set_formal_charge(-1),
        (None, Some(c)) => mol.atom_mut(c).set_formal_charge(0),
        (None, None) => {}
    }
    if let Some(c) = charged {
        explicit_valence(mol, c)?;
    }
    if let Some(u) = uncharged {
        explicit_valence(mol, u)?;
    }
    Ok(true)
}

fn n_val5_neutral(mol: &mut RWMol, atom: u32) -> Result<bool> {
    Ok(mol.atom(atom).atomic_num() == 7
        && mol.atom(atom).formal_charge() == 0
        && explicit_valence(mol, atom)? == 5)
}

/// `_Valence5NCleanUp6`: clean CN1=CCN=CC=1
fn valence5n_cleanup6(mol: &mut RWMol, atom: u32) -> Result<bool> {
    if !n_val5_neutral(mol, atom)? {
        return Ok(false);
    }
    mol.atom_mut(atom).set_atomic_num(50);
    let q = make_query(
        &[6, 6, 50, 6, 6, 7, 6],
        &[
            (0, 1, BondType::Single),
            (1, 2, BondType::Double),
            (2, 3, BondType::Double),
            (3, 4, BondType::Unspecified),
            (4, 5, BondType::Single),
            (5, 0, BondType::Double),
            (2, 6, BondType::Single),
        ],
    )?;
    let Some(m) = unique_match(mol, &q)? else {
        mol.atom_mut(atom).set_atomic_num(7);
        return Ok(false);
    };
    set_bt(mol, m[0], m[1], BondType::Double);
    set_bt(mol, m[1], m[2], BondType::Single);
    set_bt(mol, m[4], m[5], BondType::Double);
    set_bt(mol, m[5], m[0], BondType::Single);
    mol.atom_mut(atom).set_atomic_num(7);
    mol.atom_mut(atom).set_formal_charge(1);
    Ok(true)
}

/// `_Valence5NCleanUp7`: clean CN1=NCOCC=1 that is connected via alternating
/// bonds to O
fn valence5n_cleanup7(mol: &mut RWMol, atom: u32) -> Result<bool> {
    let mut stack = Vec::new();
    let Some(target) = find_alt(
        mol,
        atom,
        8,
        0,
        BondType::Double,
        BondType::Double,
        5,
        &mut stack,
    ) else {
        return Ok(false);
    };
    if !n_val5_neutral(mol, atom)? {
        return Ok(false);
    }
    mol.atom_mut(atom).set_atomic_num(50);
    let q = make_query(
        &[6, 6, 50, 7, 6, 8, 6],
        &[
            (0, 1, BondType::Unspecified),
            (1, 2, BondType::Double),
            (2, 3, BondType::Double),
            (3, 4, BondType::Single),
            (4, 5, BondType::Single),
            (5, 0, BondType::Single),
            (2, 6, BondType::Single),
        ],
    )?;
    let Some(m) = unique_match(mol, &q)? else {
        mol.atom_mut(atom).set_atomic_num(7);
        return Ok(false);
    };
    set_bt(mol, m[1], m[2], BondType::Single);
    flip_double_single(mol, &mut stack);
    mol.atom_mut(target).set_formal_charge(-1);
    mol.atom_mut(atom).set_atomic_num(7);
    Ok(true)
}

/// `_Valence5NCleanUp8`: clean [N]=C1N=CN=N1
fn valence5n_cleanup8(mol: &mut RWMol, atom: u32) -> Result<bool> {
    if !n_val5_neutral(mol, atom)? {
        return Ok(false);
    }
    mol.atom_mut(atom).set_atomic_num(50);
    let q = make_query(
        &[6, 7, 6, 7, 7, 50],
        &[
            (0, 1, BondType::Single),
            (1, 2, BondType::Double),
            (2, 3, BondType::Single),
            (3, 4, BondType::Double),
            (4, 0, BondType::Single),
            (5, 0, BondType::Double),
        ],
    )?;
    let Some(m) = unique_match(mol, &q)? else {
        mol.atom_mut(atom).set_atomic_num(7);
        return Ok(false);
    };
    set_bt(mol, m[1], m[2], BondType::Single);
    set_bt(mol, m[2], m[3], BondType::Double);
    set_bt(mol, m[3], m[4], BondType::Single);
    set_bt(mol, m[4], m[0], BondType::Double);
    set_bt(mol, m[5], m[0], BondType::Single);
    mol.atom_mut(m[1]).set_formal_charge(-1);
    mol.atom_mut(atom).set_atomic_num(7);
    mol.atom_mut(atom).set_formal_charge(1);
    Ok(true)
}

/// `_Valence5NCleanUp9`: clean [N]=C1C=CN=N1
fn valence5n_cleanup9(mol: &mut RWMol, atom: u32) -> Result<bool> {
    if !n_val5_neutral(mol, atom)? {
        return Ok(false);
    }
    mol.atom_mut(atom).set_atomic_num(50);
    let q = make_query(
        &[6, 7, 7, 6, 6, 50],
        &[
            (0, 1, BondType::Single),
            (1, 2, BondType::Double),
            (2, 3, BondType::Single),
            (3, 4, BondType::Double),
            (4, 0, BondType::Single),
            (5, 0, BondType::Double),
        ],
    )?;
    let Some(m) = unique_match(mol, &q)? else {
        mol.atom_mut(atom).set_atomic_num(7);
        return Ok(false);
    };
    set_bt(mol, m[0], m[1], BondType::Double);
    set_bt(mol, m[1], m[2], BondType::Single);
    set_bt(mol, m[5], m[0], BondType::Single);
    mol.atom_mut(m[2]).set_formal_charge(-1);
    mol.atom_mut(atom).set_atomic_num(7);
    mol.atom_mut(atom).set_formal_charge(1);
    Ok(true)
}

/// `_Valence5NCleanUpA`: N connected via alternating bonds to N=N
fn valence5n_cleanup_a(mol: &mut RWMol, atom: u32) -> Result<bool> {
    if !n_val5_neutral(mol, atom)? {
        return Ok(false);
    }
    let q = make_query(&[7, 7], &[(0, 1, BondType::Double)])?;
    let matches = substruct_match(mol, &q, &SubstructMatchParameters::default())?;
    if matches.is_empty() {
        return Ok(false);
    }
    let mut best: Vec<u32> = Vec::new();
    for m in &matches {
        let (a0, a1) = (m[0].1 as u32, m[1].1 as u32);
        if a0 == atom || a1 == atom {
            continue;
        }
        mol.atom_mut(a0).set_atomic_num(50);
        mol.atom_mut(a1).set_atomic_num(50);
        let mut stack = Vec::new();
        let target = find_alt(
            mol,
            atom,
            50,
            0,
            BondType::Double,
            BondType::Double,
            9,
            &mut stack,
        );
        if target.is_some() && (best.is_empty() || stack.len() < best.len()) {
            best = stack;
        }
        mol.atom_mut(a0).set_atomic_num(7);
        mol.atom_mut(a1).set_atomic_num(7);
    }
    if !best.is_empty() {
        while let Some(b) = best.pop() {
            let t = if mol.bond(b).bond_type() == BondType::Single {
                BondType::Double
            } else {
                BondType::Single
            };
            mol.bond_mut(b).set_bond_type(t);
        }
        mol.atom_mut(atom).set_formal_charge(1);
        explicit_valence(mol, atom)?;
        return Ok(true);
    }
    Ok(false)
}

/// `_Valence5NCleanUpB`: directly to a C via double bond; last resort
fn valence5n_cleanup_b(mol: &mut RWMol, atom: u32) -> Result<bool> {
    let mut stack = Vec::new();
    let Some(target) = find_alt(
        mol,
        atom,
        6,
        0,
        BondType::Double,
        BondType::Double,
        1,
        &mut stack,
    ) else {
        return Ok(false);
    };
    mol.atom_mut(target).set_formal_charge(-1);
    explicit_valence(mol, target)?;
    let top = *stack.last().unwrap();
    mol.bond_mut(top).set_bond_type(BondType::Single);
    mol.atom_mut(atom).set_formal_charge(1);
    explicit_valence(mol, atom)?;
    Ok(true)
}

fn s_val7_minus(mol: &mut RWMol, atom: u32) -> Result<bool> {
    Ok(mol.atom(atom).atomic_num() == 16
        && mol.atom(atom).formal_charge() == -1
        && explicit_valence(mol, atom)? == 7)
}

/// `_Valence7SCleanUp1`: C([S-](=O)(=O)=O) to C(S([O-])(=O)=O)
fn valence7s_cleanup1(mol: &mut RWMol, atom: u32) -> Result<bool> {
    if !s_val7_minus(mol, atom)? {
        return Ok(false);
    }
    let mut neighbors_c = 0;
    let mut neighbors_o = 0;
    let mut chosen: Option<u32> = None;
    for nbr in mol.atom_neighbors(atom).collect::<Vec<_>>() {
        let bt = mol
            .bond(mol.bond_idx_between(nbr, atom).unwrap())
            .bond_type();
        match mol.atom(nbr).atomic_num() {
            8 => {
                if bt != BondType::Double {
                    neighbors_o = 100;
                    break;
                }
                chosen = Some(nbr);
                neighbors_o += 1;
            }
            6 => {
                if bt != BondType::Single {
                    neighbors_c = 100;
                    break;
                }
                neighbors_c += 1;
            }
            _ => {
                neighbors_c = 100;
                break;
            }
        }
    }
    if let Some(o) = chosen
        && (neighbors_c == 1 || neighbors_o == 3)
    {
        set_bt(mol, o, atom, BondType::Single);
        mol.atom_mut(o).set_formal_charge(-1);
        mol.atom_mut(atom).set_formal_charge(0);
        explicit_valence(mol, o)?;
        explicit_valence(mol, atom)?;
        return Ok(true);
    }
    Ok(false)
}

/// `_Valence7SCleanUp2`: [S-]=CC#N
fn valence7s_cleanup2(mol: &mut RWMol, atom: u32) -> Result<bool> {
    if !s_val7_minus(mol, atom)? {
        return Ok(false);
    }
    let mut stack = Vec::new();
    let target = find_alt(
        mol,
        atom,
        7,
        0,
        BondType::Double,
        BondType::Triple,
        3,
        &mut stack,
    );
    if target.is_some() {
        while let Some(b) = stack.pop() {
            let t = match mol.bond(b).bond_type() {
                BondType::Single => BondType::Double,
                BondType::Double => BondType::Single,
                BondType::Triple => BondType::Double,
                other => other,
            };
            mol.bond_mut(b).set_bond_type(t);
        }
        mol.atom_mut(atom).set_formal_charge(0);
        explicit_valence(mol, atom)?;
        return Ok(true);
    }
    Ok(false)
}

/// `_Valence7SCleanUp3`: S- connected to a N via double bond
fn valence7s_cleanup3(mol: &mut RWMol, atom: u32) -> Result<bool> {
    if !s_val7_minus(mol, atom)? {
        return Ok(false);
    }
    let mut stack = Vec::new();
    let target = find_alt(
        mol,
        atom,
        7,
        0,
        BondType::Double,
        BondType::Double,
        1,
        &mut stack,
    );
    if let Some(t) = target {
        let top = *stack.last().unwrap();
        mol.bond_mut(top).set_bond_type(BondType::Single);
        mol.atom_mut(t).set_formal_charge(-1);
        mol.atom_mut(atom).set_formal_charge(0);
        explicit_valence(mol, atom)?;
        return Ok(true);
    }
    Ok(false)
}

/// `_Valence8SCleanUp1`: S- connected to a N via alternating bonds
fn valence8s_cleanup1(mol: &mut RWMol, atom: u32) -> Result<bool> {
    if !s_val7_minus(mol, atom)? {
        return Ok(false);
    }
    let mut stack = Vec::new();
    let Some(target) = find_alt(
        mol,
        atom,
        7,
        0,
        BondType::Double,
        BondType::Double,
        9,
        &mut stack,
    ) else {
        return Ok(false);
    };
    flip_double_single(mol, &mut stack);
    mol.atom_mut(target).set_formal_charge(-1);
    explicit_valence(mol, target)?;
    mol.atom_mut(target).set_num_explicit_hs(0);
    mol.atom_mut(atom).set_formal_charge(0);
    explicit_valence(mol, atom)?;
    Ok(true)
}

/// `_Valence8ClCleanUp1`: [Cl-](=O)(=O)(=O)(=O) to [Cl+3]([O-])([O-])([O-])[O-]
fn valence8cl_cleanup1(mol: &mut RWMol, atom: u32) -> Result<bool> {
    if explicit_valence(mol, atom)? != 8 || mol.atom(atom).formal_charge() != -1 {
        return Ok(false);
    }
    let nbrs: Vec<u32> = mol.atom_neighbors(atom).collect();
    if nbrs.iter().all(|&n| mol.atom(n).atomic_num() == 8) {
        mol.atom_mut(atom).set_formal_charge(3);
        for n in nbrs {
            let b = mol.bond_idx_between(atom, n).unwrap();
            if mol.bond(b).bond_type() == BondType::Double {
                mol.bond_mut(b).set_bond_type(BondType::Single);
                mol.atom_mut(n).set_formal_charge(-1);
                explicit_valence(mol, n)?;
            }
        }
        explicit_valence(mol, atom)?;
        return Ok(true);
    }
    Ok(false)
}

/// `_Valence5ClCleanUp1`: [Cl+][O-] to Cl=O
fn valence5cl_cleanup1(mol: &mut RWMol, atom: u32) -> Result<bool> {
    if explicit_valence(mol, atom)? != 6 || mol.atom(atom).formal_charge() != 1 {
        return Ok(false);
    }
    let mut stack = Vec::new();
    let Some(target) = find_alt(
        mol,
        atom,
        8,
        -1,
        BondType::Single,
        BondType::Single,
        1,
        &mut stack,
    ) else {
        return Ok(false);
    };
    let top = *stack.last().unwrap();
    mol.bond_mut(top).set_bond_type(BondType::Double);
    mol.atom_mut(atom).set_formal_charge(0);
    mol.atom_mut(target).set_formal_charge(0);
    explicit_valence(mol, atom)?;
    Ok(true)
}

/// `_Valence3ClCleanUp1`: Cl#S to ClS
fn valence3cl_cleanup1(mol: &mut RWMol, atom: u32) -> Result<bool> {
    if explicit_valence(mol, atom)? != 3 || mol.atom(atom).formal_charge() != 0 {
        return Ok(false);
    }
    let mut stack = Vec::new();
    let target = find_alt(
        mol,
        atom,
        16,
        0,
        BondType::Triple,
        BondType::Triple,
        1,
        &mut stack,
    );
    if target.is_none() {
        return Ok(false);
    }
    let top = *stack.last().unwrap();
    mol.bond_mut(top).set_bond_type(BondType::Single);
    explicit_valence(mol, atom)?;
    Ok(true)
}

/// `cleanUp`: make the structure returned by the InChI library acceptable to
/// RDKit.
#[allow(clippy::collapsible_match)]
fn clean_up(mol: &mut RWMol) -> Result<()> {
    for atom in 0..mol.num_atoms() {
        match mol.atom(atom).atomic_num() {
            7 => {
                if explicit_valence(mol, atom)? == 4 {
                    if valence4n_cleanup1(mol, atom)? {
                        continue;
                    }
                    if mol.atom(atom).formal_charge() == -1 && valence4n_cleanup2(mol, atom)? {
                        continue;
                    }
                    continue;
                }
                if mol.atom(atom).formal_charge() != 0 {
                    continue;
                }
                let arom_holder = mol.atom(atom).is_aromatic();
                mol.atom_mut(atom).set_is_aromatic(false);
                if explicit_valence(mol, atom)? == 5 {
                    let _ = valence5n_cleanup6(mol, atom)?
                        || valence5n_cleanup7(mol, atom)?
                        || valence5n_cleanup8(mol, atom)?
                        || valence5n_cleanup9(mol, atom)?
                        || valence5n_cleanup_a(mol, atom)?
                        || valence5n_cleanup1(mol, atom)?
                        || valence5n_cleanup2(mol, atom)?
                        || valence5n_cleanup3(mol, atom)?
                        || valence5n_cleanup4(mol, atom)?
                        || valence5n_cleanup5(mol, atom, 8)?
                        || valence5n_cleanup5(mol, atom, 16)?
                        || valence5n_cleanup5(mol, atom, 9)?
                        || valence5n_cleanup5(mol, atom, 17)?
                        || valence5n_cleanup_b(mol, atom)?;
                }
                if arom_holder {
                    mol.atom_mut(atom).set_is_aromatic(true);
                }
            }
            17 => {
                if explicit_valence(mol, atom)? == 8 && valence8cl_cleanup1(mol, atom)? {
                    continue;
                }
                if explicit_valence(mol, atom)? == 5 && valence5cl_cleanup1(mol, atom)? {
                    continue;
                }
                if explicit_valence(mol, atom)? == 3 && valence3cl_cleanup1(mol, atom)? {
                    continue;
                }
            }
            16 => {
                let v = explicit_valence(mol, atom)?;
                if v == 7 {
                    if valence7s_cleanup1(mol, atom)? {
                        continue;
                    }
                    if valence7s_cleanup2(mol, atom)? {
                        continue;
                    }
                    if valence7s_cleanup3(mol, atom)? {
                        continue;
                    }
                    valence8s_cleanup1(mol, atom)?;
                } else if v == 8 {
                    valence8s_cleanup1(mol, atom)?;
                }
            }
            35 => {
                if explicit_valence(mol, atom)? == 3
                    && mol.atom(atom).formal_charge() == 0
                    && mol.atom_degree(atom) == 1
                {
                    let nid = mol.atom_neighbors(atom).next().unwrap();
                    if mol.atom(nid).atomic_num() == 34 {
                        set_bt(mol, atom, nid, BondType::Single);
                    }
                }
            }
            _ => {}
        }
    }
    Ok(())
}

/// `InchiToMol`: a molecule from a (standard or non-standard) InChI string.
/// `Ok(None)` when the InChI library could not restore a structure.
pub fn inchi_to_mol(
    inchi: &str,
    rv: &mut ExtraInchiReturnValues,
    sanitize: bool,
    remove_hs_flag: bool,
) -> Result<Option<RWMol>> {
    let (retcode, out) = get_struct_from_inchi(Some(inchi), Some(""));
    rv.return_code = retcode;
    if let Some(m) = &out.sz_message {
        rv.message = m.clone();
    }
    if let Some(l) = &out.sz_log {
        rv.log = l.clone();
    }
    if retcode != INCHI_RET_OKAY && retcode != INCHI_RET_WARNING {
        return Ok(None);
    }

    let mut m = RWMol::new();
    let table = PeriodicTable::get();
    // (isotope, atom, repeat) of the isotopic hydrogens
    let mut isotopes: Vec<(u32, u32, u32)> = Vec::new();
    let n_atoms = out.atom.len() as u32;
    for (i, ia) in out.atom.iter().enumerate() {
        let name_len = ia
            .elname
            .iter()
            .position(|&c| c == 0)
            .unwrap_or(ia.elname.len());
        let name = std::str::from_utf8(&ia.elname[..name_len]).unwrap_or("");
        let atomic_number = table.get_atomic_number(name)?;
        let mut atom = Atom::new(atomic_number);
        let ref_weight = (atom.mass() + 0.5) as i32;
        let mut isotope = 0i32;
        if ia.isotopic_mass != 0 {
            isotope = i32::from(ia.isotopic_mass) - ISOTOPIC_SHIFT_FLAG;
        }
        if isotope != 0 {
            atom.set_isotope((isotope + ref_weight) as u32);
        }
        atom.set_formal_charge(i32::from(ia.charge));
        if ia.radical != 0 {
            if ia.radical != 3 && ia.radical != 2 {
                crate::rd_log_warning!(
                    "expect radical to be either 2 or 3 while getting {}. Ignore radical.\n",
                    ia.radical
                );
            } else {
                atom.set_num_radical_electrons(ia.radical as u32 - 1);
            }
        }
        atom.set_num_explicit_hs(ia.num_iso_h[0] as u32);
        if ia.num_iso_h[1] != 0 {
            isotopes.push((1, i as u32, ia.num_iso_h[1] as u32));
        } else if ia.num_iso_h[2] != 0 {
            isotopes.push((2, i as u32, ia.num_iso_h[2] as u32));
        } else if ia.num_iso_h[3] != 0 {
            isotopes.push((3, i as u32, ia.num_iso_h[3] as u32));
        }
        atom.set_no_implicit(true);
        m.add_atom(atom, false);
    }

    // adding bonds
    let mut register: BTreeSet<(u32, u32)> = BTreeSet::new();
    for (i, ia) in out.atom.iter().enumerate() {
        let i = i as u32;
        for b in 0..ia.num_bonds as usize {
            let nbr = ia.neighbor[b] as u32;
            if register.contains(&(i, nbr)) || register.contains(&(nbr, i)) {
                continue;
            }
            register.insert((i, nbr));
            let bt = ia.bond_type[b];
            let (btype, aromatic) = if (0..=INCHI_BOND_TYPE_TRIPLE).contains(&bt) {
                (BondType::from_u8(bt as u8).unwrap(), false)
            } else if bt == INCHI_BOND_TYPE_ALTERN {
                crate::rd_log_warning!(
                    "receive ALTERN bond type which should be avoided. This is treated as aromatic.\n"
                );
                (BondType::Aromatic, true)
            } else {
                crate::rd_log_error!("illegal bond type ({}) in InChI\n", bt as u8);
                return Ok(None);
            };
            let bidx = m.add_bond(i, nbr, btype)?;
            if aromatic {
                m.bond_mut(bidx).set_is_aromatic(true);
            }
            let s = ia.bond_stereo[b];
            if s == INCHI_BOND_STEREO_SINGLE_1UP || s == INCHI_BOND_STEREO_SINGLE_2DOWN {
                m.bond_mut(bidx).set_bond_dir(BondDir::BeginWedge);
            } else if s == INCHI_BOND_STEREO_SINGLE_1DOWN || s == INCHI_BOND_STEREO_SINGLE_2UP {
                m.bond_mut(bidx).set_bond_dir(BondDir::BeginDash);
            } else if s == INCHI_BOND_STEREO_SINGLE_1EITHER {
                m.bond_mut(bidx).set_bond_dir(BondDir::Unknown);
            } else if s == INCHI_BOND_STEREO_DOUBLE_EITHER {
                m.bond_mut(bidx).set_bond_dir(BondDir::EitherDouble);
            }
        }
    }

    // adding isotopes at the end
    for (isotope, aid, repeat) in isotopes {
        for _ in 0..repeat {
            let mut atom = Atom::new(1);
            atom.set_isotope(isotope);
            let j = m.add_atom(atom, false);
            m.add_bond(j, aid, BondType::Single)?;
        }
    }

    // basic topological structure is ready. calculate valence
    m.update_property_cache(false)?;

    // 0D stereo
    let mut e_pairs: IntPairVect = Vec::new();
    let mut z_pairs: IntPairVect = Vec::new();
    if !out.stereo_0d.is_empty() {
        let mut ranks: Vec<u32> = Vec::new();
        assign_atom_cip_ranks(&mut m, &mut ranks)?;
        for st in &out.stereo_0d {
            if i32::from(st.parity) == INCHI_PARITY_NONE {
                continue;
            }
            match st.type_ {
                INCHI_STEREO_TYPE_NONE => {}
                INCHI_STEREO_TYPE_DOUBLE_BOND => {
                    let left = st.neighbor[1] as u32;
                    let right = st.neighbor[2] as u32;
                    let original_left = st.neighbor[0] as i32;
                    let original_right = st.neighbor[3] as i32;
                    let Some(bond) = m.bond_idx_between(left, right) else {
                        crate::rd_log_warning!(
                            "Extended double-bond stereochemistry (e.g. C=C=C=C) ignored\n"
                        );
                        continue;
                    };
                    let find_nbrs = |m: &ROMol, reference: u32| -> (i32, i32) {
                        let mut nbr = -1i32;
                        let mut extra = -1i32;
                        let mut cip = -1i64;
                        for b in m.atom_bonds(reference) {
                            if b.bond_type() != BondType::Single
                                && b.bond_type() != BondType::Aromatic
                            {
                                continue;
                            }
                            let atom = b.other_atom_idx(reference).unwrap() as i32;
                            let c = i64::from(ranks[atom as usize]);
                            if c > cip {
                                if nbr >= 0 {
                                    extra = nbr;
                                }
                                nbr = atom;
                                cip = c;
                            } else {
                                extra = atom;
                            }
                        }
                        (nbr, extra)
                    };
                    let (left_nbr, extra_left) = find_nbrs(&m, left);
                    let (right_nbr, extra_right) = find_nbrs(&m, right);
                    if left_nbr < 0 || right_nbr < 0 {
                        crate::rd_log_warning!(
                            "Ignoring stereochemistry on double-bond without appropriate neighbors\n"
                        );
                        continue;
                    }
                    let switch_ez = (original_left == left_nbr && original_right != right_nbr)
                        || (original_left != left_nbr && original_right == right_nbr);
                    let mut parity = i32::from(st.parity);
                    if parity == INCHI_PARITY_ODD && switch_ez {
                        parity = INCHI_PARITY_EVEN;
                    } else if parity == INCHI_PARITY_EVEN && switch_ez {
                        parity = INCHI_PARITY_ODD;
                    }

                    let mut find_bond_pairs = |m: &ROMol, reference: u32, nbr: i32, extra: i32| {
                        let bond = m.bond_idx_between(reference, nbr as u32).unwrap();
                        if extra >= 0 {
                            let mut modifier = -1;
                            if m.bond(bond).begin_atom_idx() != reference {
                                modifier *= -1;
                            }
                            let extra_bond = m.bond_idx_between(reference, extra as u32).unwrap();
                            if m.bond(extra_bond).begin_atom_idx() != reference {
                                modifier *= -1;
                            }
                            if modifier == 1 {
                                z_pairs.push((bond, extra_bond));
                            } else {
                                e_pairs.push((bond, extra_bond));
                            }
                        }
                        bond
                    };
                    let left_bond = find_bond_pairs(&m, left, left_nbr, extra_left);
                    let right_bond = find_bond_pairs(&m, right, right_nbr, extra_right);

                    let mut modifier = -1;
                    if m.bond(left_bond).begin_atom_idx() != left {
                        modifier *= -1;
                    }
                    if m.bond(right_bond).begin_atom_idx() != right {
                        modifier *= -1;
                    }
                    if parity == INCHI_PARITY_ODD {
                        m.bond_mut(bond).set_stereo(BondStereo::StereoZ)?;
                        if modifier == 1 {
                            e_pairs.push((left_bond, right_bond));
                        } else {
                            z_pairs.push((left_bond, right_bond));
                        }
                    } else if parity == INCHI_PARITY_EVEN {
                        m.bond_mut(bond).set_stereo(BondStereo::StereoE)?;
                        if modifier == 1 {
                            z_pairs.push((left_bond, right_bond));
                        } else {
                            e_pairs.push((left_bond, right_bond));
                        }
                    } else if parity == INCHI_PARITY_NONE {
                        m.bond_mut(bond).set_stereo(BondStereo::StereoNone)?;
                    } else {
                        m.bond_mut(bond).set_stereo(BondStereo::StereoAny)?;
                    }
                    let sa = m.bond_mut(bond).stereo_atoms_mut();
                    sa.push(left_nbr);
                    sa.push(right_nbr);
                }
                INCHI_STEREO_TYPE_TETRAHEDRAL => {
                    let parity = i32::from(st.parity);
                    if parity == INCHI_PARITY_UNDEFINED || parity == INCHI_PARITY_UNKNOWN {
                        continue;
                    }
                    let c = st.central_atom as u32;
                    let mut n_swaps = 0;
                    let mut nid = 0usize;
                    if st.neighbor[0] == st.central_atom {
                        nid = 1;
                        if m.atom_degree(c) == 3 {
                            n_swaps = 1;
                        }
                    }
                    let mut neighbors: Vec<i32> = Vec::new();
                    while nid < 4 {
                        let end = st.neighbor[nid] as u32;
                        let b = m.bond_idx_between(c, end).unwrap();
                        neighbors.push(b as i32);
                        nid += 1;
                    }
                    n_swaps += m.perturbation_order(c, &neighbors)?;
                    if parity == INCHI_PARITY_ODD {
                        m.atom_mut(c).set_chiral_tag(ChiralType::TetrahedralCCW);
                    } else {
                        m.atom_mut(c).set_chiral_tag(ChiralType::TetrahedralCW);
                    }
                    if n_swaps % 2 != 0 {
                        m.atom_mut(c).invert_chirality();
                    }
                }
                INCHI_STEREO_TYPE_ALLENE => {
                    crate::rd_log_warning!(
                        "Allene-style stereochemistry is not supported yet and will be ignored.\n"
                    );
                }
                t => {
                    crate::rd_log_warning!("Unrecognized stereo0D type ({t}) is ignored!\n");
                }
            }
        }
        if !assign_bond_dirs(&mut m, &z_pairs, &e_pairs) {
            crate::rd_log_warning!("Cannot assign bond directions!\n");
        }
    }
    let _ = n_atoms;

    // clean up the molecule to be acceptable to RDKit
    clean_up(&mut m)?;
    if sanitize {
        if remove_hs_flag {
            remove_hs(&mut m, &RemoveHsParameters::default(), true)?;
        } else {
            sanitize_mol(
                &mut m,
                crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL,
            )
            .map_err(|(_, e): (u32, RDKitError)| e)?;
        }
    }
    // call assignStereochemistry just to be safe; otherwise, MolToSmiles may
    // overwrite E/Z and/or bond direction on double bonds.
    assign_stereochemistry(&mut m, true, true, false)?;
    Ok(Some(m))
}
