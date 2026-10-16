//! `SynthonSpaceSearch_details.cpp`: helpers shared by the searchers.

use std::collections::{BTreeMap, VecDeque};
use std::sync::Arc;
use std::time::Instant;

use super::hit_set::SynthonSpaceHitSet;
use super::params::{FragSetUniquifyMode, MAX_CONNECTOR_NUM};
use crate::Result;
use crate::graphmol::chem_transforms::mol_fragmenter::fragment_on_bonds;
use crate::graphmol::chem_transforms::mol_zip::{MolzipLabel, MolzipParams, molzip_single};
use crate::graphmol::chirality::find_stereo::{StereoType, find_potential_stereo};
use crate::graphmol::find_rings::find_sssr;
use crate::graphmol::mol_ops::frags::{get_mol_frags_as_mols, get_mol_frags_atoms};
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::queries::{
    CompositeQueryType, make_a_atom_query, make_atom_isotope_query, make_atom_num_query,
};
use crate::graphmol::query_ops::atom_expand_query;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::rdgeneral::threads::par_map;

/// A set of connector numbers: bit `i` is set when connector `i` is present
/// (`boost::dynamic_bitset` of `MAX_CONNECTOR_NUM + 1` bits).
pub type ConnectorBits = u32;

pub fn smiles_of(mol: &ROMol) -> String {
    let mut m = mol.clone();
    mol_to_smiles(&mut m, &SmilesWriteParams::default()).unwrap_or_default()
}

/// `checkTimeOut`
pub fn check_time_out(end_time: Option<Instant>) -> bool {
    if let Some(t) = end_time
        && Instant::now() > t
    {
        crate::rd_log_warning!("Timed out.");
        return true;
    }
    false
}

fn prev_permutation<T: Ord>(v: &mut [T]) -> bool {
    if v.len() < 2 {
        return false;
    }
    let mut i = v.len() - 1;
    while i > 0 && v[i - 1] <= v[i] {
        i -= 1;
    }
    if i == 0 {
        v.reverse();
        return false;
    }
    let mut j = v.len() - 1;
    while v[j] >= v[i - 1] {
        j -= 1;
    }
    v.swap(i - 1, j);
    v[i..].reverse();
    true
}

fn next_permutation<T: Ord>(v: &mut [T]) -> bool {
    if v.len() < 2 {
        return false;
    }
    let mut i = v.len() - 1;
    while i > 0 && v[i - 1] >= v[i] {
        i -= 1;
    }
    if i == 0 {
        v.reverse();
        return false;
    }
    let mut j = v.len() - 1;
    while v[j] <= v[i - 1] {
        j -= 1;
    }
    v.swap(i - 1, j);
    v[i..].reverse();
    true
}

/// `combMFromN`: all the ways of choosing `m` of `n`, in upstream's order.
pub fn comb_m_from_n(m: usize, n: usize) -> Vec<Vec<usize>> {
    let mut all_n = vec![1u8; m];
    all_n.resize(n, 0);
    let mut combs = Vec::new();
    loop {
        combs.push((0..n).filter(|&i| all_n[i] != 0).collect::<Vec<_>>());
        if !prev_permutation(&mut all_n) {
            break;
        }
    }
    combs
}

/// `permMFromN`
pub fn perm_m_from_n(m: usize, n: usize) -> Vec<Vec<usize>> {
    let mut perms = Vec::new();
    for mut c in comb_m_from_n(m, n) {
        loop {
            perms.push(c.clone());
            if !next_permutation(&mut c) {
                break;
            }
        }
    }
    perms
}

fn flag_ring_bonds(mol: &ROMol) -> Result<Vec<bool>> {
    let mut ring_bonds = vec![false; mol.num_bonds() as usize];
    let mut local;
    let mol = if mol.ring_info().is_initialized() {
        mol
    } else {
        // Query molecules don't have the ring info generated on creation.
        local = mol.clone();
        find_sssr(&mut local, false, false)?;
        &local
    };
    for r in mol.ring_info().bond_rings() {
        for &b in r {
            ring_bonds[b as usize] = true;
        }
    }
    Ok(ring_bonds)
}

fn add_bonds_to_list(
    mol: &ROMol,
    atom: u32,
    ring_bonds: &mut [bool],
    done_atoms: &mut [bool],
    atoms: &mut VecDeque<u32>,
    ring_block: &mut Vec<u32>,
) {
    for nbond in mol.atom_bonds(atom) {
        let bi = nbond.idx() as usize;
        if ring_bonds[bi] {
            ring_bonds[bi] = false;
            ring_block.push(bi as u32);
            let other = nbond.other_atom_idx(atom).unwrap_or(atom);
            if !done_atoms[other as usize] {
                atoms.push_back(other);
                done_atoms[other as usize] = true;
            }
        }
    }
}

fn get_ring_blocks(mol: &ROMol, mut ring_bonds: Vec<bool>) -> Vec<Vec<u32>> {
    let mut ring_blocks: Vec<Vec<u32>> = Vec::new();
    while ring_bonds.iter().any(|&b| b) {
        for bond in mol.bonds() {
            let bi = bond.idx() as usize;
            if ring_bonds[bi] {
                let mut block = vec![bi as u32];
                ring_bonds[bi] = false;
                let mut done_atoms = vec![false; mol.num_atoms() as usize];
                let mut to_do = VecDeque::new();
                add_bonds_to_list(
                    mol,
                    bond.begin_atom_idx(),
                    &mut ring_bonds,
                    &mut done_atoms,
                    &mut to_do,
                    &mut block,
                );
                add_bonds_to_list(
                    mol,
                    bond.end_atom_idx(),
                    &mut ring_bonds,
                    &mut done_atoms,
                    &mut to_do,
                    &mut block,
                );
                while let Some(next) = to_do.pop_front() {
                    add_bonds_to_list(
                        mol,
                        next,
                        &mut ring_bonds,
                        &mut done_atoms,
                        &mut to_do,
                        &mut block,
                    );
                }
                block.sort_unstable();
                ring_blocks.push(block);
            }
        }
    }
    ring_blocks
}

fn bond_pair_fragments_block(
    mol: &ROMol,
    bondi: u32,
    bondj: u32,
    num_atoms: usize,
    adj: &mut [Vec<bool>],
) -> bool {
    let bi = mol.bond(bondi);
    let bj = mol.bond(bondj);
    let (bi0, bi1) = (bi.begin_atom_idx() as usize, bi.end_atom_idx() as usize);
    let (bj0, bj1) = (bj.begin_atom_idx() as usize, bj.end_atom_idx() as usize);
    adj[bi0][bi1] = false;
    adj[bi1][bi0] = false;
    adj[bj0][bj1] = false;
    adj[bj1][bj0] = false;

    let mut atoms = 1usize;
    let mut to_do = VecDeque::from([bi0]);
    let mut done = vec![false; adj.len()];
    done[bi0] = true;
    while let Some(next) = to_do.pop_front() {
        for i in 0..adj.len() {
            if adj[next][i] && !done[i] {
                done[i] = true;
                to_do.push_back(i);
                atoms += 1;
            }
        }
    }

    adj[bi0][bi1] = true;
    adj[bi1][bi0] = true;
    adj[bj0][bj1] = true;
    adj[bj1][bj0] = true;
    atoms < num_atoms
}

fn find_bond_pairs_that_fragment(
    mol: &ROMol,
    ring_bonds: &[bool],
    ring_blocks: &[Vec<u32>],
    ring_bond_pairs: &mut Vec<(u32, u32)>,
) {
    let n = mol.num_atoms() as usize;
    let mut adj = vec![vec![false; n]; n];
    for bond in mol.bonds() {
        if ring_bonds[bond.idx() as usize] {
            adj[bond.begin_atom_idx() as usize][bond.end_atom_idx() as usize] = true;
            adj[bond.end_atom_idx() as usize][bond.begin_atom_idx() as usize] = true;
        }
    }
    for block in ring_blocks {
        let mut ok = true;
        let mut block_atoms = vec![false; n];
        for &b in block {
            let bond = mol.bond(b);
            let (a0, a1) = (bond.begin_atom_idx() as usize, bond.end_atom_idx() as usize);
            block_atoms[a0] = true;
            block_atoms[a1] = true;
            if adj[a0].iter().filter(|&&x| x).count() > 2
                || adj[a1].iter().filter(|&&x| x).count() > 2
            {
                ok = false;
            }
        }
        if ok {
            for i in 0..block.len().saturating_sub(1) {
                for j in i + 1..block.len() {
                    ring_bond_pairs.push((block[i], block[j]));
                }
            }
        } else {
            let num_atoms = block_atoms.iter().filter(|&&x| x).count();
            for i in 0..block.len().saturating_sub(1) {
                for j in i + 1..block.len() {
                    if bond_pair_fragments_block(mol, block[i], block[j], num_atoms, &mut adj) {
                        ring_bond_pairs.push((block[i], block[j]));
                    }
                }
            }
        }
    }
}

type FragEntry = (String, Option<Arc<ROMol>>);

fn make_fragments_for_mol(
    mol: &ROMol,
    split_bonds: &[u32],
    dummy_labels: &[(u32, u32)],
    max_num_frags: usize,
    ring_bonds: &[bool],
) -> FragEntry {
    let mut num_ring_bonds = 0;
    let mut num_non_ring_bonds = 0;
    for &sb in split_bonds {
        if ring_bonds[sb as usize] {
            num_ring_bonds += 1;
        } else {
            num_non_ring_bonds += 1;
        }
    }
    if 1 + num_non_ring_bonds + num_ring_bonds / 2 > max_num_frags {
        return (String::new(), None);
    }
    match fragment_on_bonds(mol, split_bonds, true, Some(dummy_labels), None, None) {
        Ok(frag_mol) => (smiles_of(&frag_mol), Some(Arc::new(frag_mol))),
        Err(_) => (String::new(), None),
    }
}

fn build_split_bonds(
    bond_pairs: &[(u32, u32)],
    max_bond_splits: usize,
    split_bonds: &mut Vec<Vec<u32>>,
) {
    if bond_pairs.len() == 1 {
        split_bonds.push(vec![bond_pairs[0].0]);
        return;
    }
    for i in 1..max_bond_splits {
        for comb in comb_m_from_n(i, bond_pairs.len()) {
            let mut next: Vec<u32> = Vec::new();
            for c in comb {
                next.push(bond_pairs[c].0);
                next.push(bond_pairs[c].1);
            }
            next.sort_unstable();
            next.dedup();
            // Each split needs a connector number, so any split set that
            // would produce one more than the SynthonSpace is set up for is
            // a bust.
            if next.len() > MAX_CONNECTOR_NUM {
                continue;
            }
            split_bonds.push(next);
        }
    }
    split_bonds.sort();
    split_bonds.dedup();
}

fn uniquify_frags_by_smiles(tmp_frags: &mut Vec<FragEntry>, max_num_frag_sets: u64) {
    libcxx_sort_by(tmp_frags, |a, b| a.0 < b.0);
    tmp_frags.dedup_by(|b, a| a.0 == b.0);
    if tmp_frags.len() as u64 > max_num_frag_sets {
        tmp_frags.truncate(max_num_frag_sets as usize);
    }
}

fn uniquify_frags_by_atoms(tmp_frags: &mut Vec<FragEntry>, max_num_frag_sets: u64) {
    let mut keep: Vec<(String, usize)> = Vec::with_capacity(tmp_frags.len());
    for (i, (_, mol)) in tmp_frags.iter().enumerate() {
        let mol = mol.as_ref().expect("fragment");
        let mut mol_frags = get_mol_frags_atoms(mol);
        let mut at_strings: Vec<String> = mol_frags
            .iter_mut()
            .map(|mf| {
                mf.sort_unstable();
                mf.iter()
                    .map(|a| a.to_string())
                    .collect::<Vec<_>>()
                    .join("_")
            })
            .collect();
        at_strings.sort();
        keep.push((at_strings.join("."), i));
    }
    libcxx_sort_by(&mut keep, |a, b| a.0 < b.0);
    keep.dedup_by(|b, a| a.0 == b.0);
    let n = (keep.len() as u64).min(max_num_frag_sets) as usize;
    let new_frags: Vec<FragEntry> = keep[..n]
        .iter()
        .map(|(_, i)| tmp_frags[*i].clone())
        .collect();
    *tmp_frags = new_frags;
}

/// `splitMolecule`: all the ways of splitting the query into up to
/// `max_num_frags` fragments. Returns an empty vector if it timed out.
pub fn split_molecule(
    query: &ROMol,
    max_num_frags: u32,
    max_num_frag_sets: u64,
    end_time: Option<Instant>,
    num_threads: i32,
    uniquify_mode: FragSetUniquifyMode,
    timed_out: &mut bool,
) -> Result<Vec<Vec<Arc<ROMol>>>> {
    let max_num_frags = (max_num_frags.max(1)).min(query.num_bonds() + 1) as usize;

    let ring_bonds = flag_ring_bonds(query)?;
    let ring_blocks = get_ring_blocks(query, ring_bonds.clone());

    let mut bond_pairs: Vec<(u32, u32)> = Vec::new();
    find_bond_pairs_that_fragment(query, &ring_bonds, &ring_blocks, &mut bond_pairs);
    for b in query.bonds() {
        if !ring_bonds[b.idx() as usize]
            && query.atom(b.begin_atom_idx()).atomic_num() != 1
            && query.atom(b.end_atom_idx()).atomic_num() != 1
        {
            bond_pairs.push((b.idx(), b.idx()));
        }
    }

    let mut split_bonds: Vec<Vec<u32>> = Vec::new();
    if !bond_pairs.is_empty() {
        build_split_bonds(&bond_pairs, max_num_frags, &mut split_bonds);
    }
    let dummy_labels: Vec<(u32, u32)> = (1..=MAX_CONNECTOR_NUM as u32).map(|i| (i, i)).collect();

    let mut tmp_frags: Vec<FragEntry> = par_map(&split_bonds, num_threads, |sb| {
        make_fragments_for_mol(query, sb, &dummy_labels, max_num_frags, &ring_bonds)
    });
    *timed_out = check_time_out(end_time);
    if *timed_out {
        return Ok(Vec::new());
    }
    tmp_frags.retain(|f| f.1.is_some());
    match uniquify_mode {
        FragSetUniquifyMode::BySmiles => {
            uniquify_frags_by_smiles(&mut tmp_frags, max_num_frag_sets)
        }
        FragSetUniquifyMode::ByAtoms => uniquify_frags_by_atoms(&mut tmp_frags, max_num_frag_sets),
    }

    // Keep the molecule itself (i.e. 0 splits).
    let mut fragments: Vec<Vec<Arc<ROMol>>> = vec![vec![Arc::new(query.clone())]];
    let finals: Vec<Option<Vec<Arc<ROMol>>>> = par_map(&tmp_frags, num_threads, |(_, m)| {
        let m = m.as_ref()?;
        let (mol_frags, _) = get_mol_frags_as_mols(m, false, true).ok()?;
        if mol_frags.len() <= max_num_frags {
            Some(mol_frags.into_iter().map(Arc::new).collect())
        } else {
            None
        }
    });
    fragments.extend(finals.into_iter().flatten().filter(|fs| !fs.is_empty()));
    *timed_out = check_time_out(end_time);
    Ok(fragments)
}

/// `countConnections`
pub fn count_connections(mol: &ROMol) -> i32 {
    mol.atoms()
        .iter()
        .filter(|a| {
            a.atomic_num() == 0 && a.isotope() >= 1 && a.isotope() <= MAX_CONNECTOR_NUM as u32
        })
        .count() as i32
}

/// `getConnectorPatterns`
pub fn get_connector_patterns(mols: &[Arc<ROMol>]) -> Vec<ConnectorBits> {
    mols.iter()
        .map(|m| {
            let mut bits = 0;
            for a in m.atoms() {
                if a.atomic_num() == 0
                    && a.isotope() != 0
                    && a.isotope() <= MAX_CONNECTOR_NUM as u32
                {
                    bits |= 1 << a.isotope();
                }
            }
            bits
        })
        .collect()
}

/// `getConnectorPattern`
pub fn get_connector_pattern(frag_set: &[Arc<ROMol>]) -> ConnectorBits {
    get_connector_patterns(frag_set)
        .into_iter()
        .fold(0, |acc, b| acc | b)
}

/// For each permutation of connector numbers, for each fragment, the
/// (atom index, new isotope) pairs to apply.
pub type ConnectorPermutation = Vec<Vec<(u32, u32)>>;

/// `getConnectorPermutations`
pub fn get_connector_permutations(
    mol_frags: &[ROMol],
    frag_conns: ConnectorBits,
    reaction_conns: ConnectorBits,
) -> Vec<ConnectorPermutation> {
    let num_frag_conns = frag_conns.count_ones() as usize;
    let perms = perm_m_from_n(num_frag_conns, reaction_conns.count_ones() as usize);
    perms
        .iter()
        .map(|perm| {
            mol_frags
                .iter()
                .map(|f| {
                    let mut done = vec![false; f.num_atoms() as usize];
                    let mut res = Vec::new();
                    for atom in f.atoms() {
                        if atom.atomic_num() == 0 {
                            for (i, p) in perm.iter().enumerate() {
                                let idx = atom.idx() as usize;
                                if !done[idx] && atom.isotope() as usize == i + 1 {
                                    res.push((atom.idx(), *p as u32 + 1));
                                    done[idx] = true;
                                }
                            }
                        }
                    }
                    res
                })
                .collect()
        })
        .collect()
}

/// `getConnectorPermutations` on connector bit patterns.
pub fn get_connector_permutations_bits(
    frag_conn_patts: &[ConnectorBits],
    reaction_conns: ConnectorBits,
) -> Vec<Vec<ConnectorBits>> {
    let conns = frag_conn_patts.iter().fold(0, |acc, b| acc | b);
    let perms = perm_m_from_n(
        conns.count_ones() as usize,
        reaction_conns.count_ones() as usize,
    );
    perms
        .iter()
        .map(|perm| {
            frag_conn_patts
                .iter()
                .map(|fcp| {
                    let mut bs = 0;
                    for (i, p) in perm.iter().enumerate() {
                        if fcp & (1 << (i + 1)) != 0 {
                            bs |= 1 << (p + 1);
                        }
                    }
                    bs
                })
                .collect()
        })
        .collect()
}

/// `expandBitSet`
pub fn expand_bit_set(bit_sets: &mut [Vec<bool>]) {
    if bit_sets.iter().any(|bs| bs.iter().any(|&b| b)) {
        for bs in bit_sets.iter_mut() {
            if !bs.iter().any(|&b| b) {
                bs.iter_mut().for_each(|b| *b = true);
            }
        }
    }
}

/// `bitSetsToVectors`
pub fn bit_sets_to_vectors(bit_sets: &[Vec<bool>]) -> Vec<Vec<usize>> {
    bit_sets
        .iter()
        .map(|bs| {
            bs.iter()
                .enumerate()
                .filter_map(|(i, &b)| b.then_some(i))
                .collect()
        })
        .collect()
}

/// `removeQueryAtoms`: makes query atoms more generic so that fragments
/// match synthons that could be built into a match.
pub fn remove_query_atoms(mol: &mut ROMol) -> Result<bool> {
    let mut did_something = false;
    for idx in 0..mol.num_atoms() {
        let atom = mol.atom(idx);
        if (atom.atomic_num() != 0 || atom.isotope() == 0)
            && atom.has_query()
            && atom.query().is_some_and(|q| q.description != "AtomType")
        {
            let mut qat = crate::graphmol::atom::Atom::new(atom.atomic_num());
            if atom.atomic_num() != 0 {
                qat.set_query(Some(make_atom_num_query(atom.atomic_num() as i32)));
            } else {
                qat.set_query(Some(make_a_atom_query()));
            }
            mol.replace_atom(idx, &qat, false)?;
            did_something = true;
        }
    }
    Ok(did_something)
}

/// `buildConnRegion`: the parts of the molecule within 3 bonds of a
/// connector atom.
pub fn build_conn_region(mol: &ROMol) -> Result<Option<ROMol>> {
    let n = mol.num_atoms() as usize;
    let mut in_frag = vec![false; n];
    for a in mol.atoms() {
        if a.atomic_num() == 0 && a.isotope() != 0 {
            in_frag[a.idx() as usize] = true;
            for n1 in mol.atom_neighbors(a.idx()) {
                in_frag[n1 as usize] = true;
                for n2 in mol.atom_neighbors(n1) {
                    in_frag[n2 as usize] = true;
                    for n3 in mol.atom_neighbors(n2) {
                        in_frag[n3 as usize] = true;
                    }
                }
            }
        }
    }
    if !in_frag.iter().any(|&b| b) {
        return Ok(None);
    }
    let mut mol_cp = mol.clone();
    mol_cp.begin_batch_edit()?;
    for idx in 0..mol_cp.num_atoms() {
        if !in_frag[idx as usize] {
            mol_cp.remove_atom(idx)?;
        } else {
            let a = mol_cp.atom_mut(idx);
            if a.atomic_num() == 0 && a.isotope() != 0 {
                a.set_isotope(1);
                if a.has_query() {
                    atom_expand_query(a, make_atom_isotope_query(1), CompositeQueryType::Or, true)?;
                }
            }
        }
    }
    mol_cp.commit_batch_edit()?;
    Ok(Some(mol_cp))
}

/// `buildProductName`
pub fn build_product_name<S: AsRef<str>>(reaction_id: &str, frag_ids: &[S]) -> String {
    let mut name = frag_ids
        .iter()
        .map(|s| s.as_ref())
        .collect::<Vec<_>>()
        .join(";");
    name.push(';');
    name.push_str(reaction_id);
    name
}

/// `buildProductName(hitset, fragNums)`
pub fn build_product_name_hs(hitset: &SynthonSpaceHitSet<'_>, frag_nums: &[usize]) -> String {
    let names: Vec<&str> = frag_nums
        .iter()
        .enumerate()
        .map(|(i, &f)| hitset.synthons_to_use[i][f].0)
        .collect();
    build_product_name(hitset.reaction.id(), &names)
}

/// `buildProduct`: zips the synthons together.
pub fn build_product(synthons: &[&ROMol]) -> Option<ROMol> {
    let mz = MolzipParams {
        label: MolzipLabel::Isotope,
        align_coordinates: true,
        ..Default::default()
    };
    let mut prod = synthons[0].clone();
    for s in &synthons[1..] {
        prod.insert_mol(s).ok()?;
    }
    let res = (|| -> Result<Option<RWMol>> {
        let Some(mut zip) = molzip_single(&prod, &mz)? else {
            return Ok(None);
        };
        sanitize_mol(&mut zip, SANITIZE_ALL).map_err(|(_, e)| e)?;
        Ok(Some(zip))
    })();
    match res {
        Ok(m) => m,
        Err(e) => {
            println!("AWOOGA :: Failed to zip {} because \n{e}", smiles_of(&prod));
            None
        }
    }
}

/// `mapFragsBySmiles`: for each distinct fragment SMILES, the (pointers of
/// the) fragments with that SMILES. Also makes sure the fragments have ring
/// info.
pub fn map_frags_by_smiles(
    frag_sets: &mut [Vec<Arc<ROMol>>],
) -> Result<BTreeMap<String, Vec<(usize, usize)>>> {
    let mut map: BTreeMap<String, Vec<(usize, usize)>> = BTreeMap::new();
    for (i, frag_set) in frag_sets.iter_mut().enumerate() {
        for (j, frag) in frag_set.iter_mut().enumerate() {
            if !frag.ring_info().is_initialized() {
                let f = Arc::get_mut(frag).expect("unshared fragment");
                find_sssr(f, false, false)?;
            }
            map.entry(smiles_of(frag)).or_default().push((i, j));
        }
    }
    Ok(map)
}

/// `countChiralAtoms`: (number of chiral atoms, number excluding those with
/// 2 or more connectors attached).
pub fn count_chiral_atoms(mol: &mut ROMol) -> Result<(u32, u32)> {
    let sis = find_potential_stereo(mol, true, true)?;
    let (mut num, mut exc) = (0, 0);
    for si in &sis {
        if si.stereo_type != StereoType::AtomTetrahedral {
            continue;
        }
        num += 1;
        let num_dummies = mol
            .atom_neighbors(si.centered_on)
            .filter(|&n| {
                let a = mol.atom(n);
                a.atomic_num() == 0 && a.isotope() <= MAX_CONNECTOR_NUM as u32
            })
            .count();
        if num_dummies < 2 {
            exc += 1;
        }
    }
    Ok((num, exc))
}

/// `sortAndUniquifyToTry`: uniquify on hit name, keeping the order of the
/// same name from different hitsets.
pub fn sort_and_uniquify_to_try(to_try: &mut Vec<(&SynthonSpaceHitSet<'_>, Vec<usize>)>) {
    let mut tmp: Vec<(usize, String)> = to_try
        .iter()
        .enumerate()
        .map(|(i, (hs, nums))| (i, build_product_name_hs(hs, nums)))
        .collect();
    tmp.sort_by(|a, b| a.1.cmp(&b.1));
    tmp.dedup_by(|b, a| a.1 == b.1);
    let new_to_try = tmp.iter().map(|(i, _)| to_try[*i].clone()).collect();
    *to_try = new_to_try;
}

/// `Stepper`: steps through all combinations of lists of different sizes,
/// (0,0,0), (0,0,1), (0,1,0) etc.
#[derive(Debug, Clone)]
pub struct Stepper {
    pub curr_state: Vec<usize>,
    pub sizes: Vec<usize>,
}

impl Stepper {
    pub fn new(sizes: &[usize]) -> Self {
        Stepper {
            curr_state: vec![0; sizes.len()],
            sizes: sizes.to_vec(),
        }
    }
    pub fn at_end(&self) -> bool {
        self.curr_state[0] == self.sizes[0]
    }
    pub fn step(&mut self) {
        if self.at_end() {
            return;
        }
        let mut i = self.curr_state.len() as i64 - 1;
        while i >= 0 {
            let iu = i as usize;
            self.curr_state[iu] += 1;
            if self.at_end() {
                return;
            }
            if self.curr_state[iu] == self.sizes[iu] {
                self.curr_state[iu] = 0;
            } else {
                break;
            }
            i -= 1;
        }
    }
}
