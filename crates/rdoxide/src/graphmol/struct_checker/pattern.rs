//! Augmented atom matching and transformation.
//!
//! Ported from RDKit `Code/GraphMol/StructChecker/Pattern.cpp` (2026.09.1).
//!
//! Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::atom_symbol::{atom_symbol_match, get_atomic_number};
use super::utilities::{Neighbourhood, log_neighbourhood, setup_neighbourhood};
use super::{AABondType, AATopology, ANY_CHARGE, AugmentedAtom, Ligand, RadicalType};
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::{Bond, BondType};
use crate::graphmol::queries::{CompositeQueryType, make_atom_num_query};
use crate::graphmol::query_ops::atom_expand_query;
use crate::graphmol::romol::ROMol;
use crate::{RDKitError, Result};

fn verbose_atom(mol: &ROMol, atom_idx: u32, neighbours: &[Neighbourhood]) {
    let mut nstr = String::new();
    for &j in &neighbours[atom_idx as usize].atoms {
        let a = mol.atom(j);
        nstr.push_str(&format!(", {}", a.symbol()));
        if a.formal_charge() != 0 {
            nstr.push_str(&format!("{:+}", a.formal_charge()));
        }
    }
    crate::rd_log_info!(
        "* CheckAtom idx={atom_idx}, {} neighbours={}{nstr}\n",
        mol.atom(atom_idx).symbol(),
        neighbours[atom_idx as usize].atoms.len()
    );
}

/// `convertBondType(AABondType)`
pub fn convert_aa_bond_type(bt: AABondType) -> BondType {
    match bt {
        AABondType::None => BondType::Unspecified,
        AABondType::Single => BondType::Single,
        AABondType::Double => BondType::Double,
        AABondType::Triple => BondType::Triple,
        AABondType::Aromatic => BondType::Aromatic,
        AABondType::SingleDouble => BondType::Unspecified,
        AABondType::SingleAromatic => BondType::OneAndAHalf,
        AABondType::DoubleAromatic => BondType::TwoAndAHalf,
        AABondType::Any | AABondType::AllBondTypes => BondType::Other,
    }
}

/// `convertBondType(BondType)`
pub fn convert_bond_type(bt: BondType) -> AABondType {
    match bt {
        BondType::Single => AABondType::Single,
        BondType::Double => AABondType::Double,
        BondType::Triple => AABondType::Triple,
        BondType::OneAndAHalf => AABondType::SingleAromatic,
        BondType::TwoAndAHalf => AABondType::DoubleAromatic,
        BondType::Aromatic => AABondType::Aromatic,
        _ => AABondType::None,
    }
}

/// `LigandMatches`: the ligand is compatible with the atom and the bond.
pub fn ligand_matches(a: &Atom, b: &Bond, l: &Ligand, use_charge: bool) -> bool {
    if l.bond_type != AABondType::Any && !is_bond_type_match(b, l.bond_type) {
        return false;
    }
    if l.radical != RadicalType::Any && a.num_radical_electrons() != l.radical as u32 {
        return false;
    }
    if (l.charge != ANY_CHARGE || use_charge)
        && l.charge != ANY_CHARGE
        && a.formal_charge() != l.charge
    {
        return false;
    }
    atom_symbol_match(&a.symbol(), &l.atom_symbol)
}

/// Replaces the atom by an atom of the given symbol, or by a query atom for
/// a comma separated list of symbols.
fn apply_atom_symbol_list(mol: &mut ROMol, symbol: &str, idx: u32) -> Result<()> {
    let number = |s: &str| {
        get_atomic_number(s).ok_or_else(|| RDKitError::Value(format!("unknown atom symbol '{s}'")))
    };
    if !symbol.contains(',') {
        mol.atom_mut(idx).set_atomic_num(number(symbol)? as i32);
        return Ok(());
    }
    let mut tokens = symbol.split(',').filter(|t| !t.is_empty());
    let Some(first) = tokens.next() else {
        return Ok(());
    };
    let first = number(first)?;
    mol.atom_mut(idx).set_atomic_num(first as i32);
    let mut a = mol.atom(idx).copy_detached();
    a.set_query(Some(make_atom_num_query(first as i32)));
    for tok in tokens {
        atom_expand_query(
            &mut a,
            make_atom_num_query(number(tok)? as i32),
            CompositeQueryType::Or,
            true,
        )?;
    }
    mol.replace_atom(idx, &a, false)
}

/// `TransformAugmentedAtoms`: replaces the atoms that match the first
/// augmented atom of a pair by the second; true if anything changed.
pub fn transform_augmented_atoms(
    mol: &mut ROMol,
    aapair: &[(AugmentedAtom, AugmentedAtom)],
    verbose: bool,
) -> Result<bool> {
    let mut transformed = false;
    let n_atoms = mol.num_atoms() as usize;
    // (neighbour atom, bond)
    let mut neighbours: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n_atoms];
    let mut atom_in_ring = vec![false; n_atoms];
    let mut bonds_to_remove: Vec<(u32, u32)> = Vec::new();
    let neighbours_d = setup_neighbourhood(mol);

    for i in 0..mol.num_bonds() {
        let bond = mol.bond(i);
        let a1 = bond.end_atom_idx();
        let a2 = bond.begin_atom_idx();
        neighbours[a1 as usize].push((a2, i));
        neighbours[a2 as usize].push((a1, i));
    }
    for ring in mol.ring_info().bond_rings() {
        for &j in ring {
            let bond = mol.bond(j as u32);
            atom_in_ring[bond.end_atom_idx() as usize] = true;
            atom_in_ring[bond.begin_atom_idx() as usize] = true;
        }
    }

    for (i, pair) in aapair.iter().enumerate() {
        let (aa1, aa2) = (&pair.0, &pair.1);
        for j in 0..n_atoms {
            let atom = mol.atom(j as u32);
            let nbrs = &neighbours[j];
            let matches = nbrs.len() == aa1.ligands.len()
                && (aa1.charge == ANY_CHARGE || atom.formal_charge() == aa1.charge)
                && (aa1.radical == RadicalType::Any
                    || atom.num_radical_electrons() == aa1.radical as u32)
                && (aa1.topology == AATopology::None
                    || (aa1.topology == AATopology::Ring && atom_in_ring[j])
                    || (aa1.topology == AATopology::Chain && !atom_in_ring[j]))
                && atom_symbol_match(&atom.symbol(), &aa1.atom_symbol);
            if !matches {
                continue;
            }
            let nl = aa1.ligands.len().max(aa2.ligands.len());
            let mut matched = 0usize;
            let mut visited = vec![false; nl];
            let mut match_a = vec![u32::MAX; nl];
            let mut match_b = vec![u32::MAX; nl];
            for &(nbr_atom_idx, nbr_bond_idx) in nbrs {
                let nbr_atom = mol.atom(nbr_atom_idx);
                let nbr_bond = mol.bond(nbr_bond_idx);
                for (l, ligand) in aa1.ligands.iter().enumerate() {
                    if visited[l] {
                        continue;
                    }
                    if (ligand.substitution_count == 0
                        || nbrs.len() == ligand.substitution_count as usize)
                        && (ligand.charge == ANY_CHARGE
                            || nbr_atom.formal_charge() == ligand.charge)
                        && (ligand.radical == RadicalType::Any
                            || nbr_atom.num_radical_electrons() == ligand.radical as u32)
                        && ((ligand.bond_type == AABondType::Any
                            || is_bond_type_match(nbr_bond, ligand.bond_type))
                            || (ligand.bond_type == AABondType::SingleDouble
                                && matches!(
                                    nbr_bond.bond_type(),
                                    BondType::Single | BondType::Double
                                )))
                        && atom_symbol_match(&nbr_atom.symbol(), &ligand.atom_symbol)
                    {
                        matched += 1;
                        match_a[l] = nbr_atom_idx;
                        match_b[l] = nbr_bond_idx;
                        visited[l] = true;
                        break;
                    }
                }
            }
            if matched != nbrs.len() {
                continue;
            }

            if verbose {
                crate::rd_log_info!(
                    "Replace {} with AAPair i={i} {} -> {}\n",
                    log_neighbourhood(mol, j as u32, &neighbours_d),
                    aa1.atom_symbol,
                    aa2.atom_symbol
                );
            }
            // change the central atom
            let jj = j as u32;
            if aa1.atom_symbol != aa2.atom_symbol {
                apply_atom_symbol_list(mol, &aa2.atom_symbol, jj)?;
            }
            let a = mol.atom_mut(jj);
            if a.formal_charge() != aa2.charge {
                a.set_formal_charge(aa2.charge);
            }
            if a.num_radical_electrons() != aa2.radical as u32 {
                a.set_num_radical_electrons(aa2.radical as u32);
            }
            // change the ligand atoms and their bonds with the central atom
            for (l, ligand) in aa2.ligands.iter().enumerate() {
                if match_a[l] == u32::MAX || l >= aa1.ligands.len() {
                    return Err(RDKitError::Index(format!(
                        "Atom index {} out of range",
                        match_a[l]
                    )));
                }
                if aa1.ligands[l].atom_symbol != ligand.atom_symbol {
                    apply_atom_symbol_list(mol, &ligand.atom_symbol, match_a[l])?;
                }
                let al = mol.atom_mut(match_a[l]);
                if al.formal_charge() != ligand.charge {
                    al.set_formal_charge(ligand.charge);
                }
                if al.num_radical_electrons() != ligand.radical as u32 {
                    al.set_num_radical_electrons(ligand.radical as u32);
                }
                if ligand.bond_type == AABondType::None {
                    // remove the bond
                    bonds_to_remove.push((jj, neighbours[j][l].0));
                } else if aa1.ligands[l].bond_type != ligand.bond_type {
                    mol.bond_mut(match_b[l])
                        .set_bond_type(convert_aa_bond_type(ligand.bond_type));
                }
            }
            transformed = true;
        }
    }
    for (a1, a2) in bonds_to_remove {
        mol.remove_bond(a1, a2)?;
    }
    Ok(transformed)
}

/// `isBondTypeMatch`
pub fn is_bond_type_match(b: &Bond, lbt: AABondType) -> bool {
    let bt = b.bond_type();
    let aromatic = b.is_aromatic();
    match lbt {
        AABondType::Any => true,
        AABondType::Aromatic => {
            aromatic
                || matches!(
                    bt,
                    BondType::Aromatic
                        | BondType::OneAndAHalf
                        | BondType::TwoAndAHalf
                        | BondType::ThreeAndAHalf
                )
        }
        AABondType::SingleAromatic => {
            bt == BondType::OneAndAHalf || (aromatic && bt == BondType::Single)
        }
        AABondType::DoubleAromatic => {
            bt == BondType::TwoAndAHalf || (aromatic && bt == BondType::Double)
        }
        AABondType::Single => {
            matches!(
                bt,
                BondType::Single | BondType::OneAndAHalf | BondType::Aromatic
            ) || aromatic
        }
        AABondType::Double => {
            matches!(
                bt,
                BondType::Double | BondType::TwoAndAHalf | BondType::Aromatic
            ) || aromatic
        }
        AABondType::Triple => matches!(
            bt,
            BondType::Triple | BondType::ThreeAndAHalf | BondType::Aromatic
        ),
        _ => false,
    }
}

/// `RecMatchNeigh`: is there a one to one assignment of the ligands to the
/// neighbours (`matched_pairs[i]` are the ligands that neighbour i fits)?
fn rec_match_neigh(matched_pairs: &[Vec<usize>], visited: &mut [bool], used: &mut [bool]) -> bool {
    for i in 0..matched_pairs.len() {
        if visited[i] {
            continue;
        }
        visited[i] = true;
        for &pair in &matched_pairs[i] {
            if used[pair] {
                continue;
            }
            used[pair] = true;
            if rec_match_neigh(matched_pairs, visited, used) {
                return true;
            }
            used[pair] = false;
        }
        visited[i] = false;
    }
    visited.iter().all(|&v| v) && used.iter().all(|&u| u)
}

/// `RecMatch`: do all the ligands of `aa` map to the neighbours of the atom?
pub fn rec_match(
    mol: &ROMol,
    atom_idx: u32,
    aa: &AugmentedAtom,
    nbp: &[Neighbourhood],
    verbose: bool,
) -> bool {
    let nbph = &nbp[atom_idx as usize];
    let mut matched_pairs: Vec<Vec<usize>> = vec![Vec::new(); nbph.atoms.len()];
    for i in 0..nbph.atoms.len() {
        let atom = mol.atom(nbph.atoms[i]);
        let bond = mol.bond(nbph.bonds[i]);
        let mut found = false;
        for (j, au_ligand) in aa.ligands.iter().enumerate() {
            if au_ligand.charge != ANY_CHARGE && atom.formal_charge() != au_ligand.charge {
                continue;
            }
            if au_ligand.radical != RadicalType::Any
                && atom.num_radical_electrons() != au_ligand.radical as u32
            {
                continue;
            }
            if au_ligand.substitution_count != 0
                && nbp[nbph.atoms[i] as usize].atoms.len() != au_ligand.substitution_count as usize
            {
                continue;
            }
            if !is_bond_type_match(bond, au_ligand.bond_type) {
                continue;
            }
            if !atom_symbol_match(&atom.symbol(), &au_ligand.atom_symbol) {
                continue;
            }
            found = true;
            matched_pairs[i].push(j);
        }
        if !found {
            return false;
        }
    }
    let mut visited = vec![false; nbph.atoms.len()];
    let mut used = vec![false; aa.ligands.len()];
    if rec_match_neigh(&matched_pairs, &mut visited, &mut used) {
        if verbose {
            crate::rd_log_info!("RecMatch ret TRUE {}\n", aa.short_name);
        }
        return true;
    }
    false
}

/// `AAMatch`: does atom `i` match the augmented atom? `atom_ring_status`
/// can be empty to ignore the topology.
pub fn aa_match(
    mol: &ROMol,
    i: u32,
    aa: &AugmentedAtom,
    atom_ring_status: &[u32],
    nbp: &[Neighbourhood],
    verbose: bool,
) -> bool {
    let atom = mol.atom(i);
    if nbp[i as usize].atoms.len() == aa.ligands.len()
        && (aa.charge == ANY_CHARGE || atom.formal_charge() == aa.charge)
        && (aa.radical == RadicalType::Any || atom.num_radical_electrons() == aa.radical as u32)
        && atom_symbol_match(&atom.symbol(), &aa.atom_symbol)
    {
        if !atom_ring_status.is_empty()
            && aa.topology == AATopology::Ring
            && atom_ring_status[i as usize] == 0
        {
            return false;
        }
        if !atom_ring_status.is_empty()
            && aa.topology == AATopology::Chain
            && atom_ring_status[i as usize] != 0
        {
            return false;
        }
        rec_match(mol, i, aa, nbp, verbose)
    } else {
        false
    }
}

/// `RingState`: how many ring bonds are attached to every atom. (Upstream
/// indexes the bond counts by the position in the ring, not the bond.)
fn ring_state(mol: &ROMol) -> (Vec<u32>, Vec<u32>) {
    let mut atom_status = vec![0u32; mol.num_atoms() as usize];
    let mut bond_status = vec![0u32; mol.num_bonds() as usize];
    if mol.num_bonds() == 0 {
        return (atom_status, bond_status);
    }
    for ring in mol.ring_info().bond_rings() {
        for i in 0..ring.len() {
            bond_status[i] += 1;
        }
    }
    for (i, &count) in bond_status.iter().enumerate() {
        if count > 0 {
            let bond = mol.bond(i as u32);
            atom_status[bond.begin_atom_idx() as usize] += 1;
            atom_status[bond.end_atom_idx() as usize] += 1;
        }
    }
    (atom_status, bond_status)
}

/// `CheckAtoms`: does every atom of the molecule match one of the augmented
/// atoms?
pub fn check_atoms(mol: &ROMol, good_atoms: &[AugmentedAtom], verbose: bool) -> bool {
    if good_atoms.is_empty() {
        return true;
    }
    let (atom_status, _) = ring_state(mol);
    let neighbours = setup_neighbourhood(mol);
    let mut nmatch = 0u32;
    for i in 0..mol.num_atoms() {
        let prevn = nmatch;
        if verbose {
            verbose_atom(mol, i, &neighbours);
        }
        let status = atom_status[i as usize];
        let mut j = 0;
        while j < good_atoms.len() {
            let good = &good_atoms[j];
            // a ring matches to ring only
            if good.topology == AATopology::Ring && status == 0 {
                if verbose && mol.atom(i).atomic_num() == 6 {
                    crate::rd_log_info!(
                        "UNMATCHED ring state RING of atom idx={i} {} status={status}\n",
                        mol.atom(i).symbol()
                    );
                }
                j += 1;
                continue;
            }
            if good.topology == AATopology::Chain && status != 0 {
                if verbose && mol.atom(i).atomic_num() == 6 {
                    crate::rd_log_info!(
                        "UNMATCHED ring state CHAIN of atom idx={i} {} status={status}\n",
                        mol.atom(i).symbol()
                    );
                }
                j += 1;
                continue;
            }
            if neighbours[i as usize].atoms.len() == good.ligands.len()
                && aa_match(mol, i, good, &atom_status, &neighbours, verbose)
            {
                if verbose {
                    crate::rd_log_info!("AAMatch idx={i} aa_idx={j} ret TRUE !\n");
                }
                nmatch += 1;
                break;
            }
            j += 1;
        }
        if verbose && nmatch == prevn {
            crate::rd_log_info!(
                "UNMATCHED atom idx={i} {} status={status} nmatch={nmatch}\n{}\n",
                mol.atom(i).symbol(),
                log_neighbourhood(mol, i, &neighbours)
            );
        }
        if verbose && j == good_atoms.len() {
            crate::rd_log_warning!("{}\n", log_neighbourhood(mol, i, &neighbours));
        }
    }
    nmatch == mol.num_atoms()
}
