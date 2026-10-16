//! `RascalMCES.cpp`: the search for the maximum clique of the modular product
//! of the line graphs of 2 molecules.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::{Arc, OnceLock};
use std::time::Instant;

use super::lap::lap_maximize;
use super::options::RascalOptions;
use super::partition_set::{BitMatrix, PartitionSet, SharedGraph};
use super::result::{RascalResult, ResultSettings, result_compare};
use crate::graphmol::new_canon::{RankOptions, rank_mol_atoms};
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::{RDKitError, Result, rd_log_error, rd_log_warning};

/// Degree and index of the atoms of each element, by descending degree.
pub type DegSeqs = BTreeMap<i32, Vec<(i32, i32)>>;

/// Everything gathered before the search proper starts.
struct StartPoint {
    /// The tier one and tier two similarities. The initial partition set is
    /// only constructed if they are both above the threshold.
    tier1_sim: f64,
    tier2_sim: f64,
    /// The lower bound on the maximum clique size.
    lower_bound: u32,
    /// A Delta-Y exchange requires extra treatment. They're rare, though.
    delta_y_poss: bool,
    /// Some parts require mol2 to be the larger molecule. This records if
    /// they have been swapped with respect to the input molecules.
    swapped: bool,
    /// Copies of the input molecules, swapped if necessary, so that mol1 is
    /// always the smaller molecule.
    mol1: ROMol,
    mol2: ROMol,
    /// The line graphs of the 2 molecules as adjacency matrices.
    adj_matrix1: Vec<Vec<i32>>,
    adj_matrix2: Vec<Vec<i32>>,
    /// Which bonds are symmetrically equivalent (-1: none).
    equiv_bonds1: Vec<i32>,
    equiv_bonds2: Vec<i32>,
    /// The initial partition set. None if the thresholds weren't met.
    part_set: Option<PartitionSet>,
}

fn sorted_degree_seqs(mol: &ROMol, deg_seqs: &mut DegSeqs) {
    for a in mol.atoms() {
        deg_seqs
            .entry(a.atomic_num() as i32)
            .or_default()
            .push((mol.atom_degree(a.idx()) as i32, a.idx() as i32));
    }
    for seq in deg_seqs.values_mut() {
        // (the order of equal degrees does not change the similarities)
        seq.sort_by_key(|a| std::cmp::Reverse(a.0));
    }
}

fn get_atom_labels(mol: &ROMol, opts: &RascalOptions) -> Vec<String> {
    let mut labels = vec![String::new(); mol.num_atoms() as usize];
    for a in mol.atoms() {
        let mut label = a.symbol();
        if opts.exact_connections_match {
            label += &format!("X{}", mol.atom_degree(a.idx()));
        }
        if !opts.ignore_atom_aromaticity && a.is_aromatic() {
            label = format!("a{label}");
        }
        labels[a.idx() as usize] = label;
    }
    labels
}

fn calc_cost(atomi_bls: &[u32], atomj_bls: &[u32]) -> i32 {
    let uniq: BTreeSet<u32> = atomi_bls.iter().copied().collect();
    uniq.iter()
        .map(|uai| {
            let ni = atomi_bls.iter().filter(|b| *b == uai).count() as i32;
            let nj = atomj_bls.iter().filter(|b| *b == uai).count() as i32;
            ni.min(nj)
        })
        .sum()
}

fn get_assignment_score(
    atom_degrees1: &[(i32, i32)],
    atom_degrees2: &[(i32, i32)],
    bond_labels1: &[u32],
    bond_labels2: &[u32],
    mol1: &ROMol,
    mol2: &ROMol,
) -> i32 {
    const BIG_SCORE: i32 = 9999;
    let mut costs = vec![vec![BIG_SCORE; atom_degrees2.len()]; atom_degrees1.len()];
    for (i, adi) in atom_degrees1.iter().enumerate() {
        let atomi_bls: Vec<u32> = mol1
            .atom_bonds(adi.1 as u32)
            .map(|b| bond_labels1[b.idx() as usize])
            .collect();
        for (j, adj) in atom_degrees2.iter().enumerate() {
            let atomj_bls: Vec<u32> = mol2
                .atom_bonds(adj.1 as u32)
                .map(|b| bond_labels2[b.idx() as usize])
                .collect();
            costs[i][j] = calc_cost(&atomi_bls, &atomj_bls);
        }
    }
    match lap_maximize(&costs) {
        // no solution for the LAP was possible.
        None => 0,
        Some((a, b)) => a.iter().zip(&b).map(|(&r, &c)| costs[r][c]).sum(),
    }
}

/// The tier 1 similarity estimate of the molecules (from the degree
/// sequences of their atoms, which are filled in).
pub fn tier1_sim(
    mol1: &ROMol,
    mol2: &ROMol,
    deg_seqs1: &mut DegSeqs,
    deg_seqs2: &mut DegSeqs,
) -> f64 {
    sorted_degree_seqs(mol1, deg_seqs1);
    sorted_degree_seqs(mol2, deg_seqs2);
    let mut vg1g2 = 0;
    let mut eg1g2 = 0;
    for (k, seq1) in deg_seqs1.iter() {
        if let Some(seq2) = deg_seqs2.get(k) {
            let n = seq1.len().min(seq2.len());
            vg1g2 += n as i32;
            for i in 0..n {
                eg1g2 += seq1[i].0.min(seq2[i].0);
            }
        }
    }
    eg1g2 /= 2;
    let num = f64::from((vg1g2 + eg1g2) * (vg1g2 + eg1g2));
    num / f64::from((mol1.num_atoms() + mol1.num_bonds()) * (mol2.num_atoms() + mol2.num_bonds()))
}

/// The tier 2 similarity estimate: as tier 1 but using the bond labels around
/// the atoms, via a linear assignment problem.
pub fn tier2_sim(
    mol1: &ROMol,
    mol2: &ROMol,
    deg_seqs1: &DegSeqs,
    deg_seqs2: &DegSeqs,
    bond_labels1: &[u32],
    bond_labels2: &[u32],
) -> f64 {
    let mut vg1g2 = 0;
    let mut eg1g2 = 0;
    for (k, seq1) in deg_seqs1 {
        if let Some(seq2) = deg_seqs2.get(k) {
            vg1g2 += seq1.len().min(seq2.len()) as i32;
            eg1g2 += get_assignment_score(seq1, seq2, bond_labels1, bond_labels2, mol1, mol2);
        }
    }
    eg1g2 /= 2;
    let num = f64::from((vg1g2 + eg1g2) * (vg1g2 + eg1g2));
    num / f64::from((mol1.num_atoms() + mol1.num_bonds()) * (mol2.num_atoms() + mol2.num_bonds()))
}

fn bond_label_strings(mol: &ROMol, opts: &RascalOptions) -> Vec<String> {
    let atom_labels = get_atom_labels(mol, opts);
    let mut labels = vec![String::new(); mol.num_bonds() as usize];
    for b in mol.bonds() {
        let bond_type = if opts.ignore_bond_orders {
            "0".to_string()
        } else {
            (b.bond_type() as u8).to_string()
        };
        // The atom labels need to be in consistent order irrespective of
        // input order.
        let mut lbl1 = atom_labels[b.begin_atom_idx() as usize].clone();
        let mut lbl2 = atom_labels[b.end_atom_idx() as usize].clone();
        if lbl1 < lbl2 {
            std::mem::swap(&mut lbl1, &mut lbl2);
        }
        labels[b.idx() as usize] = format!("{lbl1}{bond_type}{lbl2}");
    }
    labels
}

/// The bond labels of the 2 molecules as small sequential integers (that code
/// the atoms and the bond type).
pub fn get_bond_labels(mol1: &ROMol, mol2: &ROMol, opts: &RascalOptions) -> (Vec<u32>, Vec<u32>) {
    let tmp1 = bond_label_strings(mol1, opts);
    let tmp2 = bond_label_strings(mol2, opts);
    let all: BTreeSet<&String> = tmp1.iter().chain(&tmp2).collect();
    let recode = |labels: &[String]| -> Vec<u32> {
        labels
            .iter()
            .map(|l| {
                all.iter()
                    .position(|x| *x == l)
                    .expect("label is in the set") as u32
            })
            .collect()
    };
    (recode(&tmp1), recode(&tmp2))
}

/// The line graph of the molecule as an adjacency matrix, the entries being
/// the atomic number of the atom that the 2 bonds share.
fn make_line_graph(mol: &ROMol) -> Vec<Vec<i32>> {
    let n = mol.num_bonds() as usize;
    let mut adj = vec![vec![0; n]; n];
    for a in mol.atoms() {
        let bonds: Vec<usize> = mol.atom_bonds(a.idx()).map(|b| b.idx() as usize).collect();
        for &b1 in &bonds {
            for &b2 in &bonds {
                if b1 != b2 {
                    adj[b1][b2] = a.atomic_num() as i32;
                }
            }
        }
    }
    adj
}

/// The SMILES of the molecules' rings, in the order of the ring info, for
/// the check that rings match completely.
fn ring_smiles(mol: &ROMol) -> Result<Vec<String>> {
    let atom_rings = mol.ring_info().atom_rings().to_vec();
    let bond_rings = mol.ring_info().bond_rings().to_vec();
    let mut res = Vec::with_capacity(bond_rings.len());
    for i in 0..bond_rings.len() {
        let mut ring_mol = ROMol::copy_from(mol, false, -1);
        let mut in_ring = vec![false; mol.num_atoms() as usize];
        for &a in &atom_rings[i] {
            in_ring[a as usize] = true;
        }
        for &rb in &bond_rings[i] {
            ring_mol
                .bond_mut(rb as u32)
                .set_prop("ORIG_INDEX", rb, false)?;
        }
        ring_mol.begin_batch_edit()?;
        for idx in 0..ring_mol.num_bonds() {
            // (upstream looks for the property with a typo in its name, so
            // this applies to every bond)
            let (bgn, end) = {
                let b = ring_mol.bond(idx);
                (b.begin_atom_idx(), b.end_atom_idx())
            };
            if !in_ring[bgn as usize] {
                ring_mol.remove_atom(bgn)?;
            }
            if !in_ring[end as usize] {
                ring_mol.remove_atom(end)?;
            }
        }
        ring_mol.commit_batch_edit()?;
        res.push(mol_to_smiles(&mut ring_mol, &SmilesWriteParams::default())?);
    }
    Ok(res)
}

#[allow(clippy::too_many_arguments)]
fn check_rings(
    mol1: &ROMol,
    mol1_ring_smiles: &[String],
    mol1_bond_idx: usize,
    mol2: &ROMol,
    mol2_ring_smiles: &[String],
    mol2_bond_idx: usize,
    aromatic_rings_match_only: bool,
) -> Result<bool> {
    let mol1_bond = mol1.bond(mol1_bond_idx as u32);
    let mol2_bond = mol2.bond(mol2_bond_idx as u32);
    if aromatic_rings_match_only && (!mol1_bond.is_aromatic() || !mol2_bond.is_aromatic()) {
        return Ok(true);
    }
    if mol1_bond.bond_type() != mol2_bond.bond_type() {
        return Ok(false);
    }
    // If neither bond was in a ring, but they were marked aromatic, then the
    // two mols are fragments so it's ok to match these bonds.
    let mol1_bond_rings = mol1.ring_info().bond_rings();
    let mol2_bond_rings = mol2.ring_info().bond_rings();
    let in1 = mol1.ring_info().num_bond_rings(mol1_bond_idx as u32)? > 0;
    let in2 = mol2.ring_info().num_bond_rings(mol2_bond_idx as u32)? > 0;
    if !in1 && !in2 {
        return Ok(true);
    }
    for (i, r1) in mol1_bond_rings.iter().enumerate() {
        if !r1.contains(&(mol1_bond_idx as i32)) {
            continue;
        }
        for (j, r2) in mol2_bond_rings.iter().enumerate() {
            if !r2.contains(&(mol2_bond_idx as i32)) {
                continue;
            }
            if mol1_ring_smiles[i] == mol2_ring_smiles[j] {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn check_ring_matches_ring(
    mol1: &ROMol,
    mol1_bond_idx: usize,
    mol2: &ROMol,
    mol2_bond_idx: usize,
) -> Result<bool> {
    let in1 = mol1.ring_info().num_bond_rings(mol1_bond_idx as u32)? > 0;
    let in2 = mol2.ring_info().num_bond_rings(mol2_bond_idx as u32)? > 0;
    Ok(in1 == in2)
}

/// The pairs of bonds (vertices of the 2 line graphs) that are the same type.
fn build_pairs(
    mol1: &ROMol,
    vtx_labels1: &[u32],
    mol2: &ROMol,
    vtx_labels2: &[u32],
    opts: &RascalOptions,
) -> Result<Vec<(i32, i32)>> {
    let (mut mol1_ring_smiles, mut mol2_ring_smiles) = (Vec::new(), Vec::new());
    if opts.complete_aromatic_rings {
        // For these purposes, it is correct that n1cccc1 and [nH]1cccc1 match
        // - the former would be from an N-substituted pyrrole, the latter from
        // a plain one.
        let fix = |s: String| s.replace("[nH]", "n").replace("[pH]", "p");
        mol1_ring_smiles = ring_smiles(mol1)?.into_iter().map(fix).collect();
        mol2_ring_smiles = ring_smiles(mol2)?.into_iter().map(fix).collect();
    }
    let mut pairs = Vec::new();
    for (i, l1) in vtx_labels1.iter().enumerate() {
        for (j, l2) in vtx_labels2.iter().enumerate() {
            if l1 != l2 {
                continue;
            }
            // completeSmallestRings automatically implies completeAromaticRings
            // and ringMatchesRingsOnly
            if opts.complete_smallest_rings {
                if !check_rings(
                    mol1,
                    &mol1_ring_smiles,
                    i,
                    mol2,
                    &mol2_ring_smiles,
                    j,
                    false,
                )? {
                    continue;
                }
            } else if opts.complete_aromatic_rings
                && !check_rings(mol1, &mol1_ring_smiles, i, mol2, &mol2_ring_smiles, j, true)?
            {
                continue;
            }
            if !opts.complete_smallest_rings
                && opts.ring_matches_ring_only
                && !check_ring_matches_ring(mol1, i, mol2, j)?
            {
                continue;
            }
            pairs.push((i as i32, j as i32));
        }
    }
    Ok(pairs)
}

fn calc_dist_matrix(adj: &[Vec<i32>]) -> Vec<Vec<i32>> {
    let n = adj.len();
    let mut dist = vec![vec![n as i32 + 1; n]; n];
    for i in 0..n {
        dist[i][i] = 0;
        for j in 0..n {
            if i != j && adj[i][j] != 0 {
                dist[i][j] = 1;
            }
        }
    }
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                dist[i][j] = dist[i][j].min(dist[i][k] + dist[k][j]);
            }
        }
    }
    dist
}

/// The vertices (pairs of bonds) of a modular product and its edges.
type ModularProduct = (Vec<(i32, i32)>, BitMatrix);

/// The modular product (correspondence graph) of the 2 line graphs, with its
/// vertices. None if there is nothing in common, or too much.
#[allow(clippy::too_many_arguments)]
fn make_modular_product(
    mol1: &ROMol,
    adj1: &[Vec<i32>],
    vtx_labels1: &[u32],
    mol2: &ROMol,
    adj2: &[Vec<i32>],
    vtx_labels2: &[u32],
    opts: &RascalOptions,
) -> Result<Option<ModularProduct>> {
    let pairs = build_pairs(mol1, vtx_labels1, mol2, vtx_labels2, opts)?;
    if pairs.is_empty() {
        // There was nothing in common at all.
        return Ok(None);
    }
    if pairs.len() > opts.max_bond_match_pairs as usize {
        rd_log_error!(
            "Too many matching bond pairs ({}) so can't continue.",
            pairs.len()
        );
        return Ok(None);
    }
    let (dist1, dist2) = if opts.max_frag_separation > -1 {
        (calc_dist_matrix(adj1), calc_dist_matrix(adj2))
    } else {
        (Vec::new(), Vec::new())
    };
    let mut mod_prod = BitMatrix::new(pairs.len());
    for i in 0..pairs.len() - 1 {
        for j in i + 1..pairs.len() {
            let (pi, pj) = (pairs[i], pairs[j]);
            if pi.0 == pj.0 || pi.1 == pj.1 {
                continue;
            }
            let mut dists_ok = true;
            if opts.max_frag_separation != -1
                && (dist1[pi.0 as usize][pj.0 as usize] - dist2[pi.1 as usize][pj.1 as usize]).abs()
                    > opts.max_frag_separation
            {
                dists_ok = false;
            }
            if dists_ok && adj1[pi.0 as usize][pj.0 as usize] == adj2[pi.1 as usize][pj.1 as usize]
            {
                mod_prod.set(i, j);
                mod_prod.set(j, i);
            }
        }
    }
    Ok(Some((pairs, mod_prod)))
}

fn calc_lower_bound(mol1: &ROMol, mol2: &ROMol, sim_thresh: f64) -> u32 {
    let mut mol1_at_nos = BTreeSet::new();
    let mut max_at_no = 0;
    for a in mol1.atoms() {
        mol1_at_nos.insert(a.atomic_num() as i32);
        max_at_no = max_at_no.max(a.atomic_num() as i32);
    }
    let mut mol2_at_nos = vec![false; max_at_no as usize + 1];
    for a in mol2.atoms() {
        // since we're interested in the atoms that match in the 2 molecules,
        // it doesn't matter if mol2 has an atomic number higher than anything
        // in mol1 - that can't be a match.
        if (a.atomic_num() as i32) < max_at_no {
            mol2_at_nos[a.atomic_num() as usize] = true;
        }
    }
    let delta_vg1 = mol1_at_nos
        .iter()
        .filter(|&&n| !mol2_at_nos[n as usize])
        .count() as f64;
    let size = (mol1.num_atoms() + mol1.num_bonds()) * (mol2.num_atoms() + mol2.num_bonds());
    let mut lb = f64::from(size).sqrt();
    lb = lb * sim_thresh - f64::from(mol1.num_atoms()) + delta_vg1;
    if lb < 0.0 {
        lb = 0.0;
    }
    (lb as u32).max(1)
}

fn delta_y_in_clique(clique: &[u32], mol1: &ROMol, mol2: &ROMol, vtx_pairs: &[(i32, i32)]) -> bool {
    if clique.len() < 3 {
        // there must be 3 bonds for a delta-y exchange, obs.
        return false;
    }
    // Map the clique onto the 2 molecules, counting the degrees of the atoms
    // if they are involved in the clique. When sorted, they will be the same
    // if no delta-y exchange has occurred.
    let mut degs1 = vec![0; mol1.num_atoms() as usize];
    let mut degs2 = vec![0; mol2.num_atoms() as usize];
    for &mem in clique {
        let (p1, p2) = vtx_pairs[mem as usize];
        let b1 = mol1.bond(p1 as u32);
        degs1[b1.begin_atom_idx() as usize] += 1;
        degs1[b1.end_atom_idx() as usize] += 1;
        let b2 = mol2.bond(p2 as u32);
        degs2[b2.begin_atom_idx() as usize] += 1;
        degs2[b2.end_atom_idx() as usize] += 1;
    }
    degs1.retain(|&d| d != 0);
    degs1.sort_unstable();
    degs2.retain(|&d| d != 0);
    degs2.sort_unstable();
    degs1 != degs2
}

/// Whether the fragments of the clique in molecule 1 (the fragments should
/// match those in molecule 2) all have the minimum size.
fn clique_ok(clique: &[u32], opts: &RascalOptions, mol1: &ROMol, vtx_pairs: &[(i32, i32)]) -> bool {
    if opts.min_frag_size > 0 {
        // only need to do it for mol1, as the fragments should match. The
        // fragments are the connected pieces of the clique's bonds.
        let mut parent: BTreeMap<u32, u32> = BTreeMap::new();
        fn find(parent: &mut BTreeMap<u32, u32>, x: u32) -> u32 {
            let p = *parent.entry(x).or_insert(x);
            if p == x {
                return x;
            }
            let root = find(parent, p);
            parent.insert(x, root);
            root
        }
        for &mem in clique {
            let b = mol1.bond(vtx_pairs[mem as usize].0 as u32);
            let (a, c) = (
                find(&mut parent, b.begin_atom_idx()),
                find(&mut parent, b.end_atom_idx()),
            );
            if a != c {
                parent.insert(a, c);
            }
        }
        let atoms: Vec<u32> = parent.keys().copied().collect();
        let mut sizes: BTreeMap<u32, i32> = BTreeMap::new();
        for a in atoms {
            *sizes.entry(find(&mut parent, a)).or_default() += 1;
        }
        if sizes.values().any(|&s| s < opts.min_frag_size) {
            return false;
        }
    }
    true
}

#[allow(clippy::too_many_arguments)]
fn update_max_clique(
    clique: &[u32],
    delta_y_poss: bool,
    opts: &RascalOptions,
    mol1: &ROMol,
    mol2: &ROMol,
    vtx_pairs: &[(i32, i32)],
    max_cliques: &mut Vec<Vec<u32>>,
    lower_bound: &mut u32,
) {
    if !max_cliques.is_empty() && clique.len() < max_cliques[0].len() {
        return;
    }
    let did_delta_y = delta_y_poss && delta_y_in_clique(clique, mol1, mol2, vtx_pairs);
    if did_delta_y {
        return;
    }
    if max_cliques.is_empty() {
        if clique_ok(clique, opts, mol1, vtx_pairs) {
            max_cliques.push(clique.to_vec());
        }
    } else {
        let mut good_clique = false;
        let mut did_clique_ok = false;
        if clique.len() > max_cliques[0].len() {
            good_clique = clique_ok(clique, opts, mol1, vtx_pairs);
            did_clique_ok = true;
            if good_clique {
                max_cliques.clear();
            }
        }
        if !did_clique_ok {
            good_clique = clique_ok(clique, opts, mol1, vtx_pairs);
        }
        if good_clique
            && max_cliques.len() < opts.max_best_mcess as usize
            && (max_cliques.is_empty() || clique.len() == max_cliques[0].len())
        {
            max_cliques.push(clique.to_vec());
        }
    }
    if !max_cliques.is_empty() && max_cliques[0].len() as u32 > *lower_bound {
        *lower_bound = max_cliques[0].len() as u32;
    }
}

/// Checks the time every 100 steps; a timeout gives the best cliques so far.
fn check_timeout(
    start_time: &Instant,
    opts: &RascalOptions,
    clique: &[u32],
    max_cliques: &mut Vec<Vec<u32>>,
    num_steps: &mut u64,
) -> Option<u64> {
    *num_steps += 1;
    if *num_steps == 100 {
        // This clock is very convenient, but seems quite expensive. Calling it
        // every step added 10% to the runtime.
        let run_time = start_time.elapsed().as_secs();
        if run_time as i64 > i64::from(opts.timeout) {
            if max_cliques.is_empty() {
                max_cliques.push(clique.to_vec());
            } else {
                if clique.len() > max_cliques[0].len() {
                    max_cliques.clear();
                }
                if max_cliques.is_empty() || clique.len() >= max_cliques[0].len() {
                    max_cliques.push(clique.to_vec());
                }
            }
            return Some(run_time);
        }
        *num_steps = 0;
    }
    None
}

fn equivalent_root_already_done(
    root_vtx: u32,
    vtx_pairs: &[(i32, i32)],
    equiv_bonds1: &[i32],
    equiv_bonds2: &[i32],
    root_classes: &mut BTreeSet<(i32, i32)>,
) -> bool {
    let vp = vtx_pairs[root_vtx as usize];
    let new_classes = (equiv_bonds1[vp.0 as usize], equiv_bonds2[vp.1 as usize]);
    if new_classes.0 == -1 {
        return false;
    }
    if root_classes.contains(&new_classes) {
        return true;
    }
    root_classes.insert(new_classes);
    false
}

fn has_substruct_match(mol: &ROMol, query: &ROMol) -> Result<bool> {
    let ps = SubstructMatchParameters {
        max_matches: 1,
        ..Default::default()
    };
    Ok(!substruct_match(mol, query, &ps)?.is_empty())
}

fn smarts_query(smarts: &str) -> Result<ROMol> {
    mol_from_smarts(smarts, &SmartsParserParams::default())?
        .ok_or_else(|| RDKitError::Value(format!("bad SMARTS {smarts}")))
}

fn check_equivalents_allowed(mol: &ROMol) -> Result<bool> {
    const NOT_SMARTS: [&str; 4] = [
        "*~*",
        "*~*1~*~*~1",
        "*12~*~*~2~*~1",
        "*14~*(~*~2~3~4)~*~2~*~3~1",
    ];
    const NOT_STATS: [(u32, u32); 4] = [(2, 1), (4, 4), (4, 5), (5, 8)];
    static NOT_STRUCTS: OnceLock<Vec<ROMol>> = OnceLock::new();
    let structs = NOT_STRUCTS.get_or_init(|| {
        NOT_SMARTS
            .iter()
            .map(|s| smarts_query(s).expect("the SMARTS are valid"))
            .collect()
    });
    for (q, stats) in structs.iter().zip(NOT_STATS) {
        if mol.num_atoms() == stats.0 && mol.num_bonds() == stats.1 && has_substruct_match(mol, q)?
        {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The depth-first search of the partition sets for the largest cliques.
/// Returns the cliques found and, if the search timed out, the run time.
fn explore_partitions(
    starter: &StartPoint,
    start_time: &Instant,
    opts: &RascalOptions,
    max_cliques: &mut Vec<Vec<u32>>,
    lower_bound: &mut u32,
) -> Result<Option<u64>> {
    let vtx_pairs = starter
        .part_set
        .as_ref()
        .map(|p| p.graph_pairs())
        .expect("the search needs a partition set");
    let mut num_steps = 0u64;
    let mut parts: Vec<PartitionSet> = vec![starter.part_set.clone().expect("checked above")];
    let mut clique: Vec<u32> = Vec::new();
    let mut root_classes = BTreeSet::new();
    let can_do_equivs = opts.do_equiv_bond_pruning
        && check_equivalents_allowed(&starter.mol1)?
        && check_equivalents_allowed(&starter.mol2)?;
    while !parts.is_empty() {
        if opts.timeout != -1 {
            if let Some(run_time) =
                check_timeout(start_time, opts, &clique, max_cliques, &mut num_steps)
            {
                return Ok(Some(run_time));
            }
        }
        let part = parts.last().expect("not empty");
        let mut backtrack = false;
        if opts.all_best_mcess {
            if (clique.len() + part.num_parts()) < *lower_bound as usize {
                backtrack = true;
            }
        } else if clique.len() + part.num_parts() <= *lower_bound as usize {
            backtrack = true;
        }
        if !backtrack {
            let ub = clique.len() as i64 + i64::from(part.upper_bound());
            let go_deeper = if opts.all_best_mcess {
                ub >= i64::from(*lower_bound)
            } else {
                ub > i64::from(*lower_bound)
            };
            if go_deeper {
                if !part.is_empty() {
                    let mut next_part = part.clone();
                    clique.push(next_part.pop_last_vertex());
                    if clique.len() == 1
                        && can_do_equivs
                        && equivalent_root_already_done(
                            clique[0],
                            &vtx_pairs,
                            &starter.equiv_bonds1,
                            &starter.equiv_bonds2,
                            &mut root_classes,
                        )
                    {
                        clique.pop();
                        backtrack = true;
                    } else {
                        next_part.prune_vertices(*clique.last().expect("just pushed"));
                        update_max_clique(
                            &clique,
                            starter.delta_y_poss,
                            opts,
                            &starter.mol1,
                            &starter.mol2,
                            &vtx_pairs,
                            max_cliques,
                            lower_bound,
                        );
                        parts.push(next_part);
                    }
                } else {
                    backtrack = true;
                }
            } else {
                backtrack = true;
            }
        }
        if backtrack || parts.last().is_some_and(PartitionSet::is_empty) {
            while let Some(last) = parts.last_mut() {
                if last.is_empty() {
                    parts.pop();
                    clique.pop();
                } else {
                    last.pop_last_vertex();
                    if !last.is_empty() {
                        break;
                    }
                }
            }
        }
        if parts.is_empty() {
            break;
        }
    }
    Ok(None)
}

fn delta_y_exchange_possible(mol1: &ROMol, mol2: &ROMol) -> Result<bool> {
    // A Delta-y exchange is an incorrect match when a cyclopropyl ring (the
    // delta) is matched to a C(C)(C) group (the y) because they both have
    // isomorphic line graphs. This checks to see if that's something we need
    // to worry about for these molecules.
    static QUERIES: OnceLock<(ROMol, ROMol)> = OnceLock::new();
    let (delta, y) = QUERIES.get_or_init(|| {
        (
            smarts_query("C1CC1").expect("valid"),
            smarts_query("C(C)C").expect("valid"),
        )
    });
    Ok(
        (has_substruct_match(mol1, delta)? && has_substruct_match(mol2, y)?)
            || (has_substruct_match(mol2, delta)? && has_substruct_match(mol1, y)?),
    )
}

fn find_equivalent_bonds(mol: &ROMol) -> Result<Vec<i32>> {
    let mut equiv = vec![-1; mol.num_bonds() as usize];
    let mut tmp = ROMol::copy_from(mol, false, -1);
    let ranks = rank_mol_atoms(
        &mut tmp,
        RankOptions {
            break_ties: false,
            ..Default::default()
        },
    )?;
    let mut next_class = 0;
    for b1 in mol.bonds() {
        for b2 in mol.bonds() {
            if b1.idx() == b2.idx() {
                continue;
            }
            let r = |b: &crate::graphmol::bond::Bond| {
                (
                    ranks[b.begin_atom_idx() as usize],
                    ranks[b.end_atom_idx() as usize],
                )
            };
            let ((b1b, b1e), (b2b, b2e)) = (r(b1), r(b2));
            if (b1b == b2b && b1e == b2e) || (b1b == b2e && b1e == b2b) {
                let (i1, i2) = (b1.idx() as usize, b2.idx() as usize);
                if equiv[i1] == -1 && equiv[i2] == -1 {
                    equiv[i1] = next_class;
                    equiv[i2] = next_class;
                    next_class += 1;
                } else if equiv[i1] == -1 {
                    equiv[i1] = equiv[i2];
                } else if equiv[i2] == -1 {
                    equiv[i2] = equiv[i1];
                }
            }
        }
    }
    Ok(equiv)
}

/// Gives the atoms matching the SMARTS of each class of equivalent atoms an
/// atomic number of its own (110 up) so that they match each other.
fn assign_equivalent_atoms(mol: &mut ROMol, equivalent_atoms: &str) -> Result<()> {
    if equivalent_atoms.is_empty() {
        return Ok(());
    }
    let class_smarts: Vec<&str> = {
        // (boost::split with token_compress_on)
        let mut v: Vec<&str> = Vec::new();
        let mut prev_empty = false;
        for t in equivalent_atoms.split(' ') {
            if t.is_empty() {
                if prev_empty {
                    continue;
                }
                prev_empty = true;
            } else {
                prev_empty = false;
            }
            v.push(t);
        }
        v
    };
    if class_smarts.len() > 9 {
        return Err(RDKitError::Value(
            "Too many classes of equivalent atoms.  Maximum is 9.".into(),
        ));
    }
    let mut at_num = 110;
    for smt in class_smarts {
        if smt.is_empty() {
            continue;
        }
        let qmol = smarts_query(smt)?;
        let hits = substruct_match(mol, &qmol, &SubstructMatchParameters::default())?;
        for hv in &hits {
            for h in hv {
                mol.atom_mut(h.1 as u32).set_atomic_num(at_num);
            }
        }
        at_num += 1;
    }
    Ok(())
}

fn make_initial_partition_set(
    mol1: &ROMol,
    mol2: &ROMol,
    opts: &RascalOptions,
) -> Result<StartPoint> {
    let swapped = mol1.num_atoms() > mol2.num_atoms();
    let (m1, m2) = if swapped { (mol2, mol1) } else { (mol1, mol2) };
    let mut starter = StartPoint {
        tier1_sim: -1.0,
        tier2_sim: -1.0,
        lower_bound: 0,
        delta_y_poss: false,
        swapped,
        mol1: ROMol::copy_from(m1, false, -1),
        mol2: ROMol::copy_from(m2, false, -1),
        adj_matrix1: Vec::new(),
        adj_matrix2: Vec::new(),
        equiv_bonds1: Vec::new(),
        equiv_bonds2: Vec::new(),
        part_set: None,
    };
    assign_equivalent_atoms(&mut starter.mol1, &opts.equivalent_atoms)?;
    assign_equivalent_atoms(&mut starter.mol2, &opts.equivalent_atoms)?;

    let (mut deg_seqs1, mut deg_seqs2) = (DegSeqs::new(), DegSeqs::new());
    starter.tier1_sim = tier1_sim(&starter.mol1, &starter.mol2, &mut deg_seqs1, &mut deg_seqs2);
    if starter.tier1_sim < opts.similarity_threshold {
        return Ok(starter);
    }
    let (bond_labels1, bond_labels2) = get_bond_labels(&starter.mol1, &starter.mol2, opts);
    starter.tier2_sim = tier2_sim(
        &starter.mol1,
        &starter.mol2,
        &deg_seqs1,
        &deg_seqs2,
        &bond_labels1,
        &bond_labels2,
    );
    if starter.tier2_sim < opts.similarity_threshold {
        return Ok(starter);
    }

    // Get the line graphs for the two molecules as adjacency matrices.
    starter.adj_matrix1 = make_line_graph(&starter.mol1);
    starter.adj_matrix2 = make_line_graph(&starter.mol2);

    // pairs are vertices in the 2 line graphs that are the same type. The
    // modular product is the correspondence graph of the two line graphs.
    let Some((vtx_pairs, mod_prod)) = make_modular_product(
        &starter.mol1,
        &starter.adj_matrix1,
        &bond_labels1,
        &starter.mol2,
        &starter.adj_matrix2,
        &bond_labels2,
        opts,
    )?
    else {
        return Ok(starter);
    };
    starter.lower_bound = if opts.min_clique_size > 0 {
        opts.min_clique_size
    } else {
        calc_lower_bound(&starter.mol1, &starter.mol2, opts.similarity_threshold)
    };
    let graph = Arc::new(SharedGraph {
        mod_prod,
        vtx_pairs,
        vtx1_labels: bond_labels1,
        vtx2_labels: bond_labels2,
    });
    starter.part_set = Some(PartitionSet::new(graph, starter.lower_bound));

    starter.delta_y_poss = delta_y_exchange_possible(&starter.mol1, &starter.mol2)?;

    if opts.do_equiv_bond_pruning {
        // if equiv_bonds[i] and equiv_bonds[j] are equal, the bonds are
        // equivalent.
        starter.equiv_bonds1 = find_equivalent_bonds(&starter.mol1)?;
        starter.equiv_bonds2 = find_equivalent_bonds(&starter.mol2)?;
    } else {
        starter.equiv_bonds1 = vec![-1; starter.mol1.num_bonds() as usize];
        starter.equiv_bonds2 = vec![-1; starter.mol2.num_bonds() as usize];
    }
    Ok(starter)
}

fn find_mces(starter: &mut StartPoint, opts: &RascalOptions) -> Result<Vec<RascalResult>> {
    let start_time = Instant::now();
    let mut tmp_opts = opts.clone();
    if opts.single_largest_frag {
        tmp_opts.all_best_mcess = true;
    }
    let mut max_cliques: Vec<Vec<u32>> = Vec::new();
    let mut lower_bound = starter.lower_bound;
    let timed_out = match explore_partitions(
        starter,
        &start_time,
        &tmp_opts,
        &mut max_cliques,
        &mut lower_bound,
    )? {
        Some(run_time) => {
            rd_log_warning!(
                "Timed out after {} second{}.",
                run_time,
                if run_time == 1 { "" } else { "s" }
            );
            true
        }
        None => false,
    };
    starter.lower_bound = lower_bound;
    let vtx_pairs = starter
        .part_set
        .as_ref()
        .map(|p| p.graph_pairs())
        .unwrap_or_default();
    let settings = ResultSettings {
        timed_out,
        swapped: starter.swapped,
        tier1_sim: starter.tier1_sim,
        tier2_sim: starter.tier2_sim,
        ring_matches_ring_only: opts.ring_matches_ring_only,
        single_largest_frag: opts.single_largest_frag,
        min_frag_sep: opts.max_frag_separation,
        exact_connections_match: opts.exact_connections_match,
        equivalent_atoms: opts.equivalent_atoms.clone(),
        ignore_bond_orders: opts.ignore_bond_orders,
    };
    let mut results = Vec::new();
    for c in &max_cliques {
        results.push(RascalResult::new(
            &starter.mol1,
            &starter.mol2,
            &starter.adj_matrix1,
            &starter.adj_matrix2,
            c,
            &vtx_pairs,
            &settings,
        )?);
    }
    if opts.single_largest_frag {
        libcxx_sort_by(&mut results, |r1, r2| {
            let (a1, a2) = (r1.atom_matches(), r2.atom_matches());
            let (b1, b2) = (r1.bond_matches(), r2.bond_matches());
            if a1.len() == a2.len() {
                if b1.len() == b2.len() {
                    if a1 == a2 {
                        return b1 < b2;
                    }
                    return a1 < a2;
                }
                return b1.len() > b2.len();
            }
            a1.len() > a2.len()
        });

        // the singleLargestFrag method throws bits of solutions out, so there
        // may now be duplicates and results that are different sizes.
        results.dedup_by(|r1, r2| {
            r1.atom_matches() == r2.atom_matches() && r1.bond_matches() == r2.bond_matches()
        });
        if let Some(first) = results.first() {
            let (na, nb) = (first.atom_matches().len(), first.bond_matches().len());
            let mut keep = vec![true; results.len()];
            for (i, r) in results.iter().enumerate().skip(1) {
                if r.atom_matches().len() < na || r.bond_matches().len() < nb {
                    keep[i] = false;
                }
            }
            let mut k = keep.iter();
            results.retain(|_| *k.next().expect("one flag per result"));
        }
    } else {
        // If 2 cliques are the same size, this sort puts the one with the
        // smaller number of fragments first, which may have fewer atoms.
        let mut err = None;
        libcxx_sort_by(&mut results, |a, b| match result_compare(a, b) {
            Ok(v) => v,
            Err(e) => {
                err.get_or_insert(e);
                false
            }
        });
        if let Some(e) = err {
            return Err(e);
        }
    }
    Ok(results)
}

/// The MCESs of the 2 molecules: at most one unless `all_best_mcess` is set.
/// Empty if the similarity estimates are below the threshold (or, with
/// `return_empty_mces`, a result with just the estimates).
pub fn rascal_mces(mol1: &ROMol, mol2: &ROMol, opts: &RascalOptions) -> Result<Vec<RascalResult>> {
    let mut starter = make_initial_partition_set(mol1, mol2, opts)?;
    if starter.part_set.is_none() {
        if opts.return_empty_mces {
            return Ok(vec![RascalResult::empty(
                starter.tier1_sim,
                starter.tier2_sim,
            )]);
        }
        return Ok(Vec::new());
    }
    let mut results = find_mces(&mut starter, opts)?;
    if results.is_empty() && opts.return_empty_mces {
        return Ok(vec![RascalResult::empty(
            starter.tier1_sim,
            starter.tier2_sim,
        )]);
    }
    if !opts.all_best_mcess && results.len() > 1 {
        results.truncate(1);
    }
    Ok(results)
}
