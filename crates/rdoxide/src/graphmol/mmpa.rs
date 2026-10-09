//! Matched molecular pair fragmentation (Hussain and Rea).
//!
//! Ported from RDKit `Code/GraphMol/MMPA/MMPA.{h,cpp}` (2026.09.1).
//!
//! Copyright (C) 2015 Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::{BTreeMap, HashMap};

use super::atom::Atom;
use super::bond::BondType;
use super::chem_transforms::mol_fragmenter::fragment_on_bonds;
use super::find_rings::get_mol_frags;
use super::romol::ROMol;
use super::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use super::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use super::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// The default cut pattern: single, non-ring bonds from a neutral carbon.
pub const DEFAULT_PATTERN: &str = "[#6+0;!$(*=,#[!#6])]!@!=!#[*]";

/// pairs of (begin atom, end atom) indices
type BondVector = Vec<(u32, u32)>;

/// One fragmentation: the core (absent for single cuts) and the side chains.
pub type MmpaPair = (Option<ROMol>, ROMol);

/// `MMPA::detail::computeMorganCodeHash`: a cheap Morgan-like hash used to
/// spot duplicate fragmentations.
pub fn compute_morgan_code_hash(mol: &ROMol) -> u64 {
    let n_iterations = mol.num_bonds().min(5);
    let mut curr: Vec<u64> = mol
        .atoms()
        .iter()
        .map(|a| {
            let mut code = a.atomic_num();
            code |= a.isotope() << 8;
            let charge = a.formal_charge();
            code |= charge.unsigned_abs() << 16;
            if charge < 0 {
                code |= 1 << 29;
            }
            code |= u32::from(a.is_aromatic()) << 30;
            u64::from(code)
        })
        .collect();
    for _ in 0..n_iterations {
        let prev = curr.clone();
        for bond in mol.bonds() {
            let order = bond.bond_type() as u32;
            let (a1, a2) = (bond.begin_atom_idx() as usize, bond.end_atom_idx() as usize);
            // upstream truncates the codes to `unsigned` for the arithmetic
            let (v1, v2) = (prev[a1] as u32, prev[a2] as u32);
            let term = |v: u32| {
                u64::from(
                    v.wrapping_mul(v)
                        .wrapping_add(v.wrapping_add(23).wrapping_mul(order + 1721)),
                )
            };
            curr[a1] = curr[a1].wrapping_add(term(v2));
            curr[a2] = curr[a2].wrapping_add(term(v1));
        }
    }
    curr.iter().fold(0u64, |acc, &code| {
        acc.wrapping_add(code.wrapping_mul(code.wrapping_add(6849)).wrapping_add(29))
    })
}

/// Copies the bonds of `template` between atoms of the fragment; double bonds
/// are added in a second pass (`skip_double_bonds` selects the pass).
fn add_bonds_from_template(
    template: &ROMol,
    new_mol: &mut ROMol,
    new_atom_map: &HashMap<u32, u32>,
    skip_double_bonds: bool,
) -> Result<()> {
    for tb in template.bonds() {
        let is_double = tb.bond_type() == BondType::Double;
        if is_double == skip_double_bonds {
            continue;
        }
        let (Some(&a1), Some(&a2)) = (
            new_atom_map.get(&tb.begin_atom_idx()),
            new_atom_map.get(&tb.end_atom_idx()),
        ) else {
            continue;
        };
        let new_bond = new_mol.add_bond(a1, a2, tb.bond_type())? - 1;
        new_mol.bond_mut(new_bond).set_bond_dir(tb.bond_dir());
        if is_double && tb.stereo_atoms().len() == 2 {
            let mapped = |i: i32| new_atom_map.get(&(i as u32)).map(|&m| m as i32);
            if let (Some(s1), Some(s2)) =
                (mapped(tb.stereo_atoms()[0]), mapped(tb.stereo_atoms()[1]))
            {
                *new_mol.bond_mut(new_bond).stereo_atoms_mut() = vec![s1, s2];
            }
        }
    }
    Ok(())
}

/// Appends the atoms `atom_indices` of `template`, with their bonds, to `new_mol`.
fn extract_atoms(template: &ROMol, new_mol: &mut ROMol, atom_indices: &[i32]) -> Result<()> {
    let mut new_atom_map = HashMap::new();
    for &ai in atom_indices {
        let atom: Atom = template.atom(ai as u32).copy_detached();
        new_atom_map.insert(ai as u32, new_mol.add_atom(atom, true));
    }
    add_bonds_from_template(template, new_mol, &new_atom_map, true)?;
    add_bonds_from_template(template, new_mol, &new_atom_map, false)
}

fn atom_map_number(atom: &Atom) -> Option<i32> {
    atom.get_prop(cp::molAtomMapNumber)
        .and_then(RDValue::as_int)
}

/// Number of atoms of `frag` that carry an attachment label.
fn count_attachments(frag: &[i32], isotope_track: &BTreeMap<u32, u32>) -> usize {
    frag.iter()
        .filter(|&&ai| isotope_track.contains_key(&(ai as u32)))
        .count()
}

/// Whether a candidate matches (by size and Morgan hash) a stored pair.
fn is_duplicate(res: &[MmpaPair], core: &Option<ROMol>, side_chains: &ROMol) -> bool {
    res.iter().any(|(rcore, rside)| {
        let same_size =
            |a: &ROMol, b: &ROMol| a.num_atoms() == b.num_atoms() && a.num_bonds() == b.num_bonds();
        let cores_match = match (core, rcore) {
            (None, None) => true,
            (Some(c), Some(r)) => {
                same_size(c, r) && compute_morgan_code_hash(c) == compute_morgan_code_hash(r)
            }
            _ => false,
        };
        same_size(side_chains, rside)
            && cores_match
            && compute_morgan_code_hash(side_chains) == compute_morgan_code_hash(rside)
    })
}

/// Maps the dummy labels of the side chains to 1, 2, ... in canonical SMILES order.
fn canonical_atom_maps(side_chains: &ROMol) -> Result<BTreeMap<i32, i32>> {
    let mut tmp = side_chains.clone();
    let mut old_maps = vec![0i32; tmp.num_atoms() as usize];
    // the labels would take part in the canonicalization
    for idx in 0..tmp.num_atoms() {
        if tmp.atom(idx).atomic_num() == 0
            && let Some(label) = atom_map_number(tmp.atom(idx))
        {
            tmp.atom_mut(idx).clear_prop(cp::molAtomMapNumber);
            old_maps[idx as usize] = label;
        }
    }
    let params = SmilesWriteParams {
        do_isomeric_smiles: true,
        ..SmilesWriteParams::default()
    };
    mol_to_smiles(&mut tmp, &params)?;
    let order = match tmp.get_prop(cp::_smilesAtomOutputOrder) {
        Some(RDValue::VecUnsignedInt(v)) => v.clone(),
        Some(RDValue::VecInt(v)) => v.iter().map(|&i| i as u32).collect(),
        _ => return Err(RDKitError::Key(cp::_smilesAtomOutputOrder.into())),
    };
    let mut ranked: Vec<(usize, i32)> = order
        .iter()
        .enumerate()
        .filter(|&(_, &atom)| old_maps[atom as usize] > 0)
        .map(|(rank, &atom)| (rank, old_maps[atom as usize]))
        .collect();
    ranked.sort_unstable();
    let mut maps = BTreeMap::new();
    for (_, label) in ranked {
        let next = maps.len() as i32 + 1;
        maps.entry(label).or_insert(next);
    }
    Ok(maps)
}

fn relabel_dummies(mol: &mut ROMol, maps: &BTreeMap<i32, i32>) -> Result<()> {
    for idx in 0..mol.num_atoms() {
        if mol.atom(idx).atomic_num() == 0
            && let Some(label) = atom_map_number(mol.atom(idx))
        {
            let new_label = *maps
                .get(&label)
                .ok_or_else(|| RDKitError::Index("map::at".into()))?;
            mol.atom_mut(idx)
                .set_prop(cp::molAtomMapNumber, new_label, false)?;
        }
    }
    Ok(())
}

/// Cuts `bonds_selected` and, when the result is new, appends it to `res`.
fn add_result(
    res: &mut Vec<MmpaPair>,
    mol: &ROMol,
    bonds_selected: &BondVector,
    max_cuts: usize,
) -> Result<()> {
    let mut bond_indices = Vec::with_capacity(bonds_selected.len());
    let mut dummy_labels = Vec::with_capacity(bonds_selected.len());
    for (isotope, &(a1, a2)) in (1u32..).zip(bonds_selected) {
        let bond = mol
            .get_bond_between_atoms(a1, a2)?
            .ok_or_else(|| RDKitError::Value("bond not found".into()))?;
        bond_indices.push(bond.idx());
        dummy_labels.push((isotope, isotope));
    }
    let n_cuts = bond_indices.len() as u32;
    let mut fragmented =
        fragment_on_bonds(mol, &bond_indices, true, Some(&dummy_labels), None, None)?;
    let mut isotope_track = BTreeMap::new();
    for ai in mol.num_atoms()..fragmented.num_atoms() {
        let isotope = fragmented.atom(ai).isotope();
        if isotope == 0 {
            return Err(RDKitError::Value("isotope should be >0".into()));
        }
        let atom = fragmented.atom_mut(ai);
        atom.set_prop(cp::molAtomMapNumber, isotope as i32, false)?;
        atom.set_isotope(0);
        isotope_track.insert(ai, isotope);
    }

    let mut core: Option<ROMol> = None;
    let side_chains = if n_cuts == 1 {
        fragmented.clone()
    } else {
        let frags = get_mol_frags(&fragmented);
        // is there a fragment with maxCuts connection points?
        if n_cuts >= 3
            && !frags
                .iter()
                .any(|f| count_attachments(f, &isotope_track) >= max_cuts)
        {
            return Ok(());
        }
        let mut side_chains = ROMol::new();
        let mut i_core = None;
        let mut max_attachments = 0;
        for (i, frag) in frags.iter().enumerate() {
            let n_attachments = count_attachments(frag, &isotope_track);
            max_attachments = max_attachments.max(n_attachments);
            if n_attachments == 1 {
                extract_atoms(&fragmented, &mut side_chains, frag)?;
            } else if n_attachments >= max_attachments {
                // the core is the fragment with the most connection points
                i_core = Some(i);
            }
        }
        if let Some(i) = i_core {
            let mut c = ROMol::new();
            extract_atoms(&fragmented, &mut c, &frags[i])?;
            core = Some(c);
        }
        side_chains
    };
    let mut side_chains = side_chains;

    if is_duplicate(res, &core, &side_chains) {
        return Ok(());
    }
    // the first side chain gets label 1, the second 2, ...; the core follows
    let maps = canonical_atom_maps(&side_chains)?;
    if let Some(core) = core.as_mut() {
        relabel_dummies(core, &maps)?;
    }
    relabel_dummies(&mut side_chains, &maps)?;
    res.push((core, side_chains));
    Ok(())
}

fn process_cuts(
    start: usize,
    min_cuts: usize,
    max_cuts: usize,
    bonds_selected: &mut BondVector,
    matching_bonds: &[BondVector],
    mol: &ROMol,
    res: &mut Vec<MmpaPair>,
) -> Result<()> {
    if max_cuts < min_cuts {
        return Err(RDKitError::Value(
            "supplied maxCuts is less than minCuts".into(),
        ));
    }
    if min_cuts == 0 {
        return Err(RDKitError::Value("minCuts must be greater than 0".into()));
    }
    for x in start..matching_bonds.len() {
        let before = bonds_selected.len();
        bonds_selected.extend_from_slice(&matching_bonds[x]);
        if bonds_selected.len() >= min_cuts {
            add_result(res, mol, bonds_selected, max_cuts)?;
        }
        if bonds_selected.len() < max_cuts {
            process_cuts(
                x + 1,
                min_cuts,
                max_cuts,
                bonds_selected,
                matching_bonds,
                mol,
                res,
            )?;
        }
        bonds_selected.truncate(before);
    }
    Ok(())
}

/// `fragmentMol(mol, res, minCuts, maxCuts, maxCutBonds, pattern)`: cuts the
/// bonds matching `pattern` (one bond per match). `None` when the molecule is
/// not fragmented (no match or more than `max_cut_bonds` of them).
pub fn fragment_mol(
    mol: &ROMol,
    min_cuts: usize,
    max_cuts: usize,
    max_cut_bonds: usize,
    pattern: &str,
) -> Result<Option<Vec<MmpaPair>>> {
    let smarts = mol_from_smarts(pattern, &SmartsParserParams::default())?
        .ok_or_else(|| RDKitError::Value(format!("bad SMARTS: {pattern}")))?;
    let matches = substruct_match(mol, &smarts, &SubstructMatchParameters::default())?;
    if matches.is_empty() {
        return Ok(None);
    }
    // the pattern is assumed to describe a single bond
    let matching_bonds: Vec<BondVector> = matches
        .iter()
        .map(|m| vec![(m[0].1 as u32, m[1].1 as u32)])
        .collect();
    if matching_bonds.len() > max_cut_bonds {
        return Ok(None);
    }
    let mut res = Vec::new();
    process_cuts(
        0,
        min_cuts,
        max_cuts,
        &mut Vec::new(),
        &matching_bonds,
        mol,
        &mut res,
    )?;
    Ok(Some(res))
}

/// `fragmentMol(mol, res, bondsToCut, minCuts, maxCuts)`: cuts the given bonds.
pub fn fragment_mol_at_bonds(
    mol: &ROMol,
    bonds_to_cut: &[u32],
    min_cuts: usize,
    max_cuts: usize,
) -> Result<Vec<MmpaPair>> {
    let matching_bonds: Vec<BondVector> = bonds_to_cut
        .iter()
        .map(|&i| {
            let bond = mol.bond(i);
            vec![(bond.begin_atom_idx(), bond.end_atom_idx())]
        })
        .collect();
    let mut res = Vec::new();
    process_cuts(
        0,
        min_cuts,
        max_cuts,
        &mut Vec::new(),
        &matching_bonds,
        mol,
        &mut res,
    )?;
    Ok(res)
}
