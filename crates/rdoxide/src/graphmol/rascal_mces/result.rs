//! `RascalResult`: one MCES of 2 molecules.

use std::sync::Arc;
use std::sync::OnceLock;
use std::sync::atomic::{AtomicI32, Ordering};

use crate::Result;
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::{Bond, BondType};
use crate::graphmol::kekulize::kekulize_if_possible;
use crate::graphmol::matrices::get_distance_mat;
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, remove_hs};
use crate::graphmol::mol_ops::frags::{get_mol_frags_as_mols, get_mol_frags_mapping};
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::queries::{
    CompositeQueryType, make_atom_aliphatic_query, make_atom_aromatic_query,
    make_atom_explicit_degree_query, make_atom_in_ring_query, make_atom_num_query,
    make_bond_is_in_ring_query, make_bond_null_query, make_bond_order_equals_query,
};
use crate::graphmol::query_ops::{atom_expand_query, bond_expand_query};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::smarts_write::mol_to_smarts;
use crate::graphmol::smiles_parse::write::SmilesWriteParams;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;

/// A lazily computed integer; -1 means not yet computed.
#[derive(Debug)]
struct Cached(AtomicI32);

impl Cached {
    fn unset() -> Self {
        Cached(AtomicI32::new(-1))
    }
    fn get(&self) -> i32 {
        self.0.load(Ordering::Relaxed)
    }
    fn set(&self, v: i32) {
        self.0.store(v, Ordering::Relaxed);
    }
}

impl Clone for Cached {
    fn clone(&self) -> Self {
        Cached(AtomicI32::new(self.get()))
    }
}

/// One MCES of 2 molecules, as bond and atom matches between them.
#[derive(Debug, Clone)]
pub struct RascalResult {
    mol1: Option<Arc<ROMol>>,
    mol2: Option<Arc<ROMol>>,
    mces_mol: OnceLock<Arc<ROMol>>,
    bond_matches: Vec<(i32, i32)>,
    atom_matches: Vec<(i32, i32)>,
    smarts: OnceLock<String>,
    timed_out: bool,
    tier1_sim: f64,
    tier2_sim: f64,
    ring_matches_ring_only: bool,
    max_frag_sep: i32,
    exact_connections_match: bool,
    equivalent_atoms: String,
    ignore_bond_orders: bool,
    // These are used for sorting the results.
    num_frags: Cached,
    ring_non_ring_bond_score: Cached,
    atom_match_score: Cached,
    max_delta_atom_atom_dist: Cached,
    largest_frag_size: Cached,
}

/// The settings of a search that matter to its results.
#[derive(Debug, Clone)]
pub struct ResultSettings {
    pub timed_out: bool,
    pub swapped: bool,
    pub tier1_sim: f64,
    pub tier2_sim: f64,
    pub ring_matches_ring_only: bool,
    pub single_largest_frag: bool,
    pub min_frag_sep: i32,
    pub exact_connections_match: bool,
    pub equivalent_atoms: String,
    pub ignore_bond_orders: bool,
}

/// The atom (if any) that 2 bonds have in common.
fn common_atom_in_bonds(b1: &Bond, b2: &Bond) -> i32 {
    if b1.begin_atom_idx() == b2.begin_atom_idx() {
        b1.begin_atom_idx() as i32
    } else if b1.end_atom_idx() == b2.begin_atom_idx() {
        b1.end_atom_idx() as i32
    } else if b1.begin_atom_idx() == b2.end_atom_idx() {
        b1.begin_atom_idx() as i32
    } else if b1.end_atom_idx() == b2.end_atom_idx() {
        b1.end_atom_idx() as i32
    } else {
        -1
    }
}

fn orig_index(p: &impl RDProps) -> Option<i32> {
    p.get_prop("ORIG_INDEX").and_then(|v| v.as_int())
}

impl RascalResult {
    /// The result of the clique, given as indices into the pairs of bonds,
    /// found on (possibly swapped) copies of the molecules.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mol1: &ROMol,
        mol2: &ROMol,
        adj_matrix1: &[Vec<i32>],
        adj_matrix2: &[Vec<i32>],
        clique: &[u32],
        vtx_pairs: &[(i32, i32)],
        settings: &ResultSettings,
    ) -> Result<Self> {
        let (m1, m2, mol1_adj) = if settings.swapped {
            (mol2, mol1, adj_matrix2)
        } else {
            (mol1, mol2, adj_matrix1)
        };
        let mut res = RascalResult {
            mol1: Some(Arc::new(m1.clone())),
            mol2: Some(Arc::new(m2.clone())),
            mces_mol: OnceLock::new(),
            bond_matches: extract_clique(clique, vtx_pairs, settings.swapped),
            atom_matches: Vec::new(),
            smarts: OnceLock::new(),
            timed_out: settings.timed_out,
            tier1_sim: settings.tier1_sim,
            tier2_sim: settings.tier2_sim,
            ring_matches_ring_only: settings.ring_matches_ring_only,
            max_frag_sep: settings.min_frag_sep,
            exact_connections_match: settings.exact_connections_match,
            equivalent_atoms: settings.equivalent_atoms.clone(),
            ignore_bond_orders: settings.ignore_bond_orders,
            num_frags: Cached::unset(),
            ring_non_ring_bond_score: Cached::unset(),
            atom_match_score: Cached::unset(),
            max_delta_atom_atom_dist: Cached::unset(),
            largest_frag_size: Cached::unset(),
        };
        res.match_clique_atoms(mol1_adj)?;
        if res.max_frag_sep != -1 {
            res.apply_max_frag_sep()?;
        }
        if settings.single_largest_frag {
            res.largest_frag_only()?;
        }
        Ok(res)
    }

    /// A result with no MCES, only the estimates of the similarity.
    pub fn empty(tier1_sim: f64, tier2_sim: f64) -> Self {
        RascalResult {
            mol1: None,
            mol2: None,
            mces_mol: OnceLock::new(),
            bond_matches: Vec::new(),
            atom_matches: Vec::new(),
            smarts: OnceLock::new(),
            timed_out: false,
            tier1_sim,
            tier2_sim,
            ring_matches_ring_only: false,
            max_frag_sep: -1,
            exact_connections_match: false,
            equivalent_atoms: String::new(),
            ignore_bond_orders: false,
            num_frags: Cached::unset(),
            ring_non_ring_bond_score: Cached::unset(),
            atom_match_score: Cached::unset(),
            max_delta_atom_atom_dist: Cached::unset(),
            largest_frag_size: Cached::unset(),
        }
    }

    pub fn bond_matches(&self) -> &[(i32, i32)] {
        &self.bond_matches
    }

    pub fn atom_matches(&self) -> &[(i32, i32)] {
        &self.atom_matches
    }

    pub fn timed_out(&self) -> bool {
        self.timed_out
    }

    pub fn tier1_sim(&self) -> f64 {
        self.tier1_sim
    }

    pub fn tier2_sim(&self) -> f64 {
        self.tier2_sim
    }

    /// Cuts the result down to the single largest fragment. This is
    /// irrecoverably destructive.
    pub fn largest_frag_only(&mut self) -> Result<()> {
        self.largest_frags_only(1)
    }

    pub fn largest_frags_only(&mut self, num_frags: u32) -> Result<()> {
        let Some(mol1_frags) = self.make_mol_frags(1)? else {
            return Ok(());
        };
        let (mut frags, _) = get_mol_frags_as_mols(&mol1_frags, false, true)?;
        if num_frags < 1 || (frags.len() as u32) < num_frags {
            return Ok(());
        }
        let mut order: Vec<usize> = (0..frags.len()).collect();
        libcxx_sort_by(&mut order, |&a, &b| {
            let (f1, f2) = (&frags[a], &frags[b]);
            if f1.num_atoms() == f2.num_atoms() {
                f1.num_bonds() > f2.num_bonds()
            } else {
                f1.num_atoms() > f2.num_atoms()
            }
        });
        order.truncate(num_frags as usize);
        let keep: Vec<RWMol> = order
            .iter()
            .map(|&i| std::mem::replace(&mut frags[i], RWMol::new()))
            .collect();
        self.rebuild_from_frags(&keep);
        Ok(())
    }

    pub fn trim_small_frags(&mut self, min_frag_size: u32) -> Result<()> {
        let Some(mol1_frags) = self.make_mol_frags(1)? else {
            return Ok(());
        };
        let (mut frags, _) = get_mol_frags_as_mols(&mol1_frags, false, true)?;
        frags.retain(|f| f.num_atoms() >= min_frag_size);
        self.rebuild_from_frags(&frags);
        Ok(())
    }

    /// The Johnson similarity of the 2 molecules.
    pub fn similarity(&self) -> f64 {
        match (&self.mol1, &self.mol2) {
            (Some(m1), Some(m2)) => {
                johnson_similarity(&self.bond_matches, &self.atom_matches, m1, m2)
            }
            _ => 0.0,
        }
    }

    /// Assuming the frags are all part of the original MCES, cuts it down to
    /// what's in the frags.
    fn rebuild_from_frags(&mut self, frags: &[RWMol]) {
        // Force the re-creation of the SMARTS and other properties next time
        // they're needed.
        self.smarts = OnceLock::new();
        self.max_frag_sep = -1;
        self.ring_non_ring_bond_score.set(-1);
        self.max_delta_atom_atom_dist.set(-1);
        self.largest_frag_size.set(-1);

        let (n_atoms, n_bonds) = match (&self.mol1, &self.mol2) {
            (Some(m1), Some(m2)) => (
                m1.num_atoms().max(m2.num_atoms()),
                m1.num_bonds().max(m2.num_bonds()),
            ),
            _ => (0, 0),
        };
        let mut frag_atoms = vec![false; n_atoms as usize];
        let mut frag_bonds = vec![false; n_bonds as usize];
        for f in frags {
            for a in f.atoms() {
                if let Some(i) = orig_index(a) {
                    frag_atoms[i as usize] = true;
                }
            }
            for b in f.bonds() {
                if let Some(i) = orig_index(b) {
                    frag_bonds[i as usize] = true;
                }
            }
        }
        self.atom_matches.retain(|am| frag_atoms[am.0 as usize]);
        self.bond_matches.retain(|bm| frag_bonds[bm.0 as usize]);
        self.num_frags.set(frags.len() as i32);
        self.largest_frag_size
            .set(frags.first().map_or(0, |f| f.num_atoms() as i32));
    }

    /// The SMARTS of the MCES.
    pub fn smarts(&self) -> Result<String> {
        if self.mol1.is_none() || self.mol2.is_none() {
            return Ok(String::new());
        }
        if let Some(s) = self.smarts.get() {
            if !s.is_empty() {
                return Ok(s.clone());
            }
        }
        let s = self.create_smarts_string()?;
        let _ = self.smarts.set(s.clone());
        Ok(s)
    }

    fn create_smarts_string(&self) -> Result<String> {
        let (Some(mol1), Some(mol2)) = (&self.mol1, &self.mol2) else {
            return Ok(String::new());
        };
        let mut smarts_mol = RWMol::new();
        let mut atom_map = std::collections::BTreeMap::new();
        let (mol1_rings, mol2_rings) = (mol1.ring_info(), mol2.ring_info());
        for &(i1, i2) in &self.atom_matches {
            let mut a = Atom::new(0);
            let mol1_atom = mol1.atom(i1 as u32);
            a.set_query(Some(make_atom_num_query(mol1_atom.atomic_num() as i32)));
            let mol2_atom = mol2.atom(i2 as u32);
            if mol1_atom.atomic_num() != mol2_atom.atomic_num() {
                atom_expand_query(
                    &mut a,
                    make_atom_num_query(mol2_atom.atomic_num() as i32),
                    CompositeQueryType::Or,
                    true,
                )?;
            }
            let and = CompositeQueryType::And;
            if mol1_atom.is_aromatic() && mol2_atom.is_aromatic() {
                atom_expand_query(&mut a, make_atom_aromatic_query(), and, true)?;
            } else if !mol1_atom.is_aromatic() && !mol2_atom.is_aromatic() {
                atom_expand_query(&mut a, make_atom_aliphatic_query(), and, true)?;
            }
            if self.ring_matches_ring_only
                && !mol1_atom.is_aromatic()
                && !mol2_atom.is_aromatic()
                && mol1_rings.num_atom_rings(i1 as u32)? > 0
                && mol2_rings.num_atom_rings(i2 as u32)? > 0
            {
                atom_expand_query(&mut a, make_atom_in_ring_query(), and, true)?;
            }
            if self.exact_connections_match {
                atom_expand_query(
                    &mut a,
                    make_atom_explicit_degree_query(mol1.atom_degree(i1 as u32) as i32),
                    and,
                    true,
                )?;
            }
            let ai = smarts_mol.add_atom(a, true);
            atom_map.insert(i1, ai);
        }

        for &(b1, b2) in &self.bond_matches {
            let mut b = Bond::new(BondType::Unspecified);
            let mol1_bond = mol1.bond(b1 as u32);
            let mol2_bond = mol2.bond(b2 as u32);
            b.set_begin_atom_idx_unchecked(atom_map[&(mol1_bond.begin_atom_idx() as i32)]);
            b.set_end_atom_idx_unchecked(atom_map[&(mol1_bond.end_atom_idx() as i32)]);
            if self.ignore_bond_orders {
                b.set_query(Some(make_bond_null_query()));
            } else {
                b.set_query(Some(make_bond_order_equals_query(mol1_bond.bond_type())));
                if mol1_bond.bond_type() != mol2_bond.bond_type() {
                    bond_expand_query(
                        &mut b,
                        make_bond_order_equals_query(mol2_bond.bond_type()),
                        CompositeQueryType::Or,
                        true,
                    )?;
                }
            }
            if self.ring_matches_ring_only
                && !mol1_bond.is_aromatic()
                && !mol2_bond.is_aromatic()
                && mol1_rings.num_bond_rings(b1 as u32)? > 0
                && mol2_rings.num_bond_rings(b2 as u32)? > 0
            {
                bond_expand_query(
                    &mut b,
                    make_bond_is_in_ring_query(),
                    CompositeQueryType::And,
                    true,
                )?;
            }
            smarts_mol.add_bond_object(b)?;
        }
        let mut smt = mol_to_smarts(&mut smarts_mol, &SmilesWriteParams::default())?;
        clean_smarts(&mut smt, &self.equivalent_atoms);
        Ok(smt)
    }

    /// Matches the atoms of the cliques' bonds: the atoms with 2 matching
    /// bonds incident on them first, then the isolated bonds.
    fn match_clique_atoms(&mut self, mol1_adj: &[Vec<i32>]) -> Result<()> {
        if self.bond_matches.is_empty() {
            return Ok(());
        }
        let (mol1, mol2) = match (&self.mol1, &self.mol2) {
            (Some(a), Some(b)) => (Arc::clone(a), Arc::clone(b)),
            _ => return Ok(()),
        };
        let mut mol1_matches = vec![-1i32; mol1.num_atoms() as usize];
        // set the clique atoms to -2 in mol1_matches, to mark them as yet
        // undecided.
        for bm in &self.bond_matches {
            let bond1 = mol1.bond(bm.0 as u32);
            mol1_matches[bond1.begin_atom_idx() as usize] = -2;
            mol1_matches[bond1.end_atom_idx() as usize] = -2;
        }

        // First, use the line graphs to match atoms that have 2 matching
        // bonds incident on them.
        for i in 0..self.bond_matches.len() - 1 {
            let pair1 = self.bond_matches[i];
            let bond1_1 = mol1.bond(pair1.0 as u32);
            let bond2_1 = mol2.bond(pair1.1 as u32);
            for j in i + 1..self.bond_matches.len() {
                let pair2 = self.bond_matches[j];
                if mol1_adj[pair1.0 as usize][pair2.0 as usize] != 0 {
                    // the 2 bonds are incident on the same atom, so the 2
                    // atoms must match
                    let bond1_2 = mol1.bond(pair2.0 as u32);
                    let bond2_2 = mol2.bond(pair2.1 as u32);
                    let mol1_atom = common_atom_in_bonds(bond1_1, bond1_2);
                    let mol2_atom = common_atom_in_bonds(bond2_1, bond2_2);
                    if mol1_atom != -1 {
                        mol1_matches[mol1_atom as usize] = mol2_atom;
                        let omol1 = bond1_1.other_atom_idx(mol1_atom as u32)?;
                        let omol2 = bond2_1.other_atom_idx(mol2_atom as u32)?;
                        mol1_matches[omol1 as usize] = omol2 as i32;
                        let omol1 = bond1_2.other_atom_idx(mol1_atom as u32)?;
                        let omol2 = bond2_2.other_atom_idx(mol2_atom as u32)?;
                        mol1_matches[omol1 as usize] = omol2 as i32;
                    }
                }
            }
        }
        // if there are -2 entries in mol1_matches there's more to do.
        if mol1_matches.contains(&-2) {
            // Any -2 entries are down to isolated bonds, which are a bit
            // tricky.
            for pair1 in &self.bond_matches {
                let bond1_1 = mol1.bond(pair1.0 as u32);
                let (b1b, b1e) = (bond1_1.begin_atom_idx(), bond1_1.end_atom_idx());
                if mol1_matches[b1b as usize] == -2 && mol1_matches[b1e as usize] == -2 {
                    let bond2_1 = mol2.bond(pair1.1 as u32);
                    let (b2b, b2e) = (bond2_1.begin_atom_idx(), bond2_1.end_atom_idx());
                    let straight = if mol1.atom(b1b).atomic_num() != mol1.atom(b1e).atomic_num() {
                        // it's fairly straightforward:
                        mol1.atom(b1b).atomic_num() == mol2.atom(b2b).atomic_num()
                    } else if mol1.total_num_hs(b1b, false)? != mol1.total_num_hs(b1e, false)? {
                        // try it on number of hydrogens
                        mol1.total_num_hs(b1b, false)? > mol1.total_num_hs(b1e, false)?
                    } else {
                        // it probably doesn't matter
                        true
                    };
                    if straight {
                        mol1_matches[b1b as usize] = b2b as i32;
                        mol1_matches[b1e as usize] = b2e as i32;
                    } else {
                        mol1_matches[b1b as usize] = b2e as i32;
                        mol1_matches[b1e as usize] = b2b as i32;
                    }
                }
            }
        }
        for (i, &m) in mol1_matches.iter().enumerate() {
            if m >= 0 {
                self.atom_matches.push((i as i32, m));
            }
        }
        Ok(())
    }

    /// If the clique involves a fragment that is more than the maximum
    /// separation from any other frag in either molecule, discards the smaller
    /// frag.
    fn apply_max_frag_sep(&mut self) -> Result<()> {
        let (Some(mol1), Some(mol2)) = (self.mol1.clone(), self.mol2.clone()) else {
            return Ok(());
        };
        let Some(mol1_frags) = self.make_mol_frags(1)? else {
            return Ok(());
        };
        let (frags1, _) = get_mol_frags_as_mols(&mol1_frags, false, true)?;
        if frags1.len() < 2 {
            return Ok(());
        }
        let frag_frag_dist = |f1: &RWMol, f2: &RWMol, path: &[f64], num_atoms: usize| -> i32 {
            let mut min_dist = i32::MAX;
            for a1 in f1.atoms() {
                let at1 = orig_index(a1).unwrap_or(0) as usize;
                for a2 in f2.atoms() {
                    let at2 = orig_index(a2).unwrap_or(0) as usize;
                    let dist = path[at1 * num_atoms + at2].round_ties_even() as i32;
                    min_dist = min_dist.min(dist);
                }
            }
            min_dist
        };
        let Some(mol2_frags) = self.make_mol_frags(2)? else {
            return Ok(());
        };
        let (frags2, _) = get_mol_frags_as_mols(&mol2_frags, false, true)?;
        let mol1_dists = get_distance_mat(&mol1);
        let mol2_dists = get_distance_mat(&mol2);
        let mut frags1: Vec<Option<RWMol>> = frags1.into_iter().map(Some).collect();
        let mut frags2: Vec<Option<RWMol>> = frags2.into_iter().map(Some).collect();

        let mut deleted = false;
        for i in 0..frags1.len() - 1 {
            for j in i + 1..frags1.len() {
                let (Some(f1i), Some(f1j)) = (&frags1[i], &frags1[j]) else {
                    continue;
                };
                let (Some(f2i), Some(f2j)) = (
                    frags2.get(i).and_then(Option::as_ref),
                    frags2.get(j).and_then(Option::as_ref),
                ) else {
                    continue;
                };
                let d1 = frag_frag_dist(f1i, f1j, &mol1_dists, mol1.num_atoms() as usize);
                let d2 = frag_frag_dist(f2i, f2j, &mol2_dists, mol2.num_atoms() as usize);
                if d1 > self.max_frag_sep || d2 > self.max_frag_sep {
                    deleted = true;
                    let drop = if f1i.num_atoms() < f1j.num_atoms() {
                        i
                    } else {
                        j
                    };
                    frags1[drop] = None;
                    if drop < frags2.len() {
                        frags2[drop] = None;
                    }
                }
            }
        }

        if deleted {
            // rebuild the bond and atom matches
            let mut new_bond_matches = Vec::new();
            let mut new_atom_matches = Vec::new();
            for frag in frags1.iter().flatten() {
                for b in frag.bonds() {
                    let b_idx = orig_index(b).unwrap_or(-1);
                    if let Some(bm) = self.bond_matches.iter().find(|bm| bm.0 == b_idx) {
                        new_bond_matches.push(*bm);
                    }
                }
            }
            for frag in frags1.iter().flatten() {
                for a in frag.atoms() {
                    let a_idx = orig_index(a).unwrap_or(-1);
                    if let Some(am) = self.atom_matches.iter().find(|am| am.0 == a_idx) {
                        new_atom_matches.push(*am);
                    }
                }
            }
            self.bond_matches = new_bond_matches;
            self.atom_matches = new_atom_matches;
        }
        Ok(())
    }

    /// The fragments of molecule 1 or 2 that are in the MCES, with the
    /// original atom and bond indices in the `ORIG_INDEX` properties.
    fn make_mol_frags(&self, mol_num: u32) -> Result<Option<RWMol>> {
        let the_mol = match mol_num {
            1 => &self.mol1,
            2 => &self.mol2,
            _ => return Ok(None),
        };
        let Some(the_mol) = the_mol else {
            return Ok(None);
        };
        let mut mol_frags = ROMol::copy_from(the_mol, false, -1);
        let mut a_in_clique = vec![false; the_mol.num_atoms() as usize];
        for am in &self.atom_matches {
            a_in_clique[(if mol_num == 1 { am.0 } else { am.1 }) as usize] = true;
        }
        let mut b_in_clique = vec![false; the_mol.num_bonds() as usize];
        for bm in &self.bond_matches {
            b_in_clique[(if mol_num == 1 { bm.0 } else { bm.1 }) as usize] = true;
        }
        mol_frags.begin_batch_edit()?;
        for idx in 0..mol_frags.num_atoms() {
            if !a_in_clique[idx as usize] {
                mol_frags.remove_atom(idx)?;
            } else {
                mol_frags
                    .atom_mut(idx)
                    .set_prop("ORIG_INDEX", idx as i32, false)?;
            }
        }
        for idx in 0..mol_frags.num_bonds() {
            let (bgn, end) = {
                let b = mol_frags.bond(idx);
                (b.begin_atom_idx(), b.end_atom_idx())
            };
            if !b_in_clique[idx as usize] {
                mol_frags.remove_bond(bgn, end)?;
            } else {
                mol_frags
                    .bond_mut(idx)
                    .set_prop("ORIG_INDEX", idx as i32, false)?;
            }
        }
        mol_frags.commit_batch_edit()?;
        Ok(Some(mol_frags))
    }

    fn calc_ring_non_ring_score(&self) -> Result<i32> {
        let (Some(mol1), Some(mol2)) = (&self.mol1, &self.mol2) else {
            return Ok(0);
        };
        let mut score = 0;
        for bm in &self.bond_matches {
            let nbr1 = mol1.ring_info().num_bond_rings(bm.0 as u32)?;
            let nbr2 = mol2.ring_info().num_bond_rings(bm.1 as u32)?;
            if (nbr1 > 0) != (nbr2 > 0) {
                score += 1;
            }
        }
        Ok(score)
    }

    fn calc_atom_match_score(&self) -> Result<i32> {
        let (Some(mol1), Some(mol2)) = (&self.mol1, &self.mol2) else {
            return Ok(0);
        };
        let mut score = 0;
        for am in &self.atom_matches {
            let h1 = mol1.total_num_hs(am.0 as u32, false)? as i32;
            let h2 = mol2.total_num_hs(am.1 as u32, false)? as i32;
            score += (h1 - h2).abs();
        }
        Ok(score)
    }

    fn calc_max_delta_atom_atom_dist_score(&self) -> i32 {
        let (Some(mol1), Some(mol2)) = (&self.mol1, &self.mol2) else {
            return 0;
        };
        // Possibly this could be improved, to be the total of the minimum
        // distances between each fragment.
        if self.atom_matches.is_empty() {
            return 0;
        }
        let mol1_dists = get_distance_mat(mol1);
        let mol2_dists = get_distance_mat(mol2);
        let dist = |i1: i32, i2: i32, dists: &[f64], n: u32| -> i32 {
            dists[(i1 as usize) * n as usize + i2 as usize].round_ties_even() as i32
        };
        let mut score = 0;
        for i in 0..self.atom_matches.len() - 1 {
            for j in i + 1..self.atom_matches.len() {
                let d1 = dist(
                    self.atom_matches[i].0,
                    self.atom_matches[j].0,
                    &mol1_dists,
                    mol1.num_atoms(),
                );
                let d2 = dist(
                    self.atom_matches[i].1,
                    self.atom_matches[j].1,
                    &mol2_dists,
                    mol2.num_atoms(),
                );
                score = score.max((d1 - d2).abs());
            }
        }
        score
    }

    fn calc_largest_frag_size(&self) -> Result<i32> {
        if self.mol1.is_none() || self.mol2.is_none() {
            return Ok(0);
        }
        let Some(mol1_frags) = self.make_mol_frags(1)? else {
            return Ok(0);
        };
        let (n_frags, mapping) = get_mol_frags_mapping(&mol1_frags);
        let size = |i: i32| mapping.iter().filter(|&&m| m == i).count() as i32;
        let mut lfs = size(0);
        for i in 1..n_frags as i32 {
            lfs = lfs.max(size(i));
        }
        Ok(lfs)
    }

    /// The number of contiguous fragments in the MCES.
    pub fn num_frags(&self) -> Result<i32> {
        if self.mol1.is_none() || self.mol2.is_none() {
            return Ok(0);
        }
        if self.num_frags.get() == -1 {
            if let Some(mol1_frags) = self.make_mol_frags(1)? {
                self.num_frags
                    .set(get_mol_frags_mapping(&mol1_frags).0 as i32);
            }
        }
        Ok(self.num_frags.get())
    }

    /// How many bonds in the clique don't match cyclic/non-cyclic, i.e. count
    /// as a match in the MCES but are ring bonds in one of the molecules and
    /// not in the other.
    pub fn ring_non_ring_bond_score(&self) -> Result<i32> {
        if self.mol1.is_none() || self.mol2.is_none() {
            return Ok(0);
        }
        if self.ring_non_ring_bond_score.get() == -1 {
            self.ring_non_ring_bond_score
                .set(self.calc_ring_non_ring_score()?);
        }
        Ok(self.ring_non_ring_bond_score.get())
    }

    /// A score for how well the atoms in the clique from mol1 match the atoms
    /// of mol2: the difference in H count for matching atoms, summed. So that,
    /// for example, an OH in mol1 that could match an OH or OMe matches the OH
    /// for preference.
    pub fn atom_match_score(&self) -> Result<i32> {
        if self.mol1.is_none() || self.mol2.is_none() {
            return Ok(0);
        }
        if self.atom_match_score.get() == -1 {
            self.atom_match_score.set(self.calc_atom_match_score()?);
        }
        Ok(self.atom_match_score.get())
    }

    /// The maximum difference in through-bond distance for pairs of matching
    /// atoms in the 2 molecules.
    pub fn max_delta_atom_atom_dist(&self) -> i32 {
        if self.mol1.is_none() || self.mol2.is_none() {
            return 0;
        }
        if self.max_delta_atom_atom_dist.get() == -1 {
            self.max_delta_atom_atom_dist
                .set(self.calc_max_delta_atom_atom_dist_score());
        }
        self.max_delta_atom_atom_dist.get()
    }

    /// The number of atoms in the largest contiguous fragment of the MCES.
    pub fn largest_frag_size(&self) -> Result<u32> {
        if self.mol1.is_none() || self.mol2.is_none() {
            return Ok(0);
        }
        if self.largest_frag_size.get() == -1 {
            self.largest_frag_size.set(self.calc_largest_frag_size()?);
        }
        Ok(self.largest_frag_size.get() as u32)
    }

    /// The MCES as a molecule (a sanitized piece of molecule 1).
    pub fn mces_mol(&self) -> Result<Option<Arc<ROMol>>> {
        if let Some(m) = self.mces_mol.get() {
            return Ok(Some(Arc::clone(m)));
        }
        let Some(mol1) = &self.mol1 else {
            return Ok(None);
        };
        let mut mol1_bonds = vec![false; mol1.num_bonds() as usize];
        for bm in &self.bond_matches {
            mol1_bonds[bm.0 as usize] = true;
        }
        let mut mol1_atoms = vec![false; mol1.num_atoms() as usize];
        for am in &self.atom_matches {
            mol1_atoms[am.0 as usize] = true;
        }
        let mut tmp = ROMol::copy_from(mol1, false, -1);
        kekulize_if_possible(&mut tmp, true, true, 100)?;
        tmp.begin_batch_edit()?;
        for idx in 0..tmp.num_bonds() {
            if mol1_bonds[idx as usize] {
                continue;
            }
            let (bo, bgn, end) = {
                let b = tmp.bond(idx);
                (b.bond_type(), b.begin_atom_idx(), b.end_atom_idx())
            };
            // the bond order, as a number of Hs, is upstream's enum value
            let bo = bo as u32;
            for at in [bgn, end] {
                let a = tmp.atom(at);
                if a.no_implicit() || (a.is_aromatic() && a.atomic_num() != 6) {
                    let n = a.num_explicit_hs() + bo;
                    tmp.atom_mut(at).set_num_explicit_hs(n);
                }
            }
            tmp.remove_bond(bgn, end)?;
        }
        for idx in 0..tmp.num_atoms() {
            if !mol1_atoms[idx as usize] {
                tmp.remove_atom(idx)?;
            }
        }
        tmp.commit_batch_edit()?;
        remove_hs(&mut tmp, &RemoveHsParameters::default(), true)?;
        sanitize_mol(&mut tmp, SANITIZE_ALL).map_err(|(_, e)| e)?;
        let arc = Arc::new(tmp);
        let _ = self.mces_mol.set(Arc::clone(&arc));
        Ok(Some(arc))
    }
}

/// The bond matches of a clique, sorted.
fn extract_clique(clique: &[u32], vtx_pairs: &[(i32, i32)], swapped: bool) -> Vec<(i32, i32)> {
    let mut matches: Vec<(i32, i32)> = clique
        .iter()
        .map(|&m| {
            let p = vtx_pairs[m as usize];
            if swapped { (p.1, p.0) } else { p }
        })
        .collect();
    matches.sort_unstable();
    matches
}

/// `details::resultCompare`: whether `res1` ranks before `res2` (more bonds,
/// fewer fragments, bigger largest fragment, ...).
pub fn result_compare(res1: &RascalResult, res2: &RascalResult) -> Result<bool> {
    use std::cmp::Ordering::{Equal, Greater, Less};
    let by = |o: std::cmp::Ordering| -> Option<bool> {
        match o {
            Less => Some(true),
            Greater => Some(false),
            Equal => None,
        }
    };
    if let Some(r) = by(res2.bond_matches.len().cmp(&res1.bond_matches.len())) {
        return Ok(r);
    }
    if let Some(r) = by(res1.num_frags()?.cmp(&res2.num_frags()?)) {
        return Ok(r);
    }
    if let Some(r) = by(res2.largest_frag_size()?.cmp(&res1.largest_frag_size()?)) {
        return Ok(r);
    }
    if let Some(r) = by(res1
        .ring_non_ring_bond_score()?
        .cmp(&res2.ring_non_ring_bond_score()?))
    {
        return Ok(r);
    }
    if let Some(r) = by(res1.atom_match_score()?.cmp(&res2.atom_match_score()?)) {
        return Ok(r);
    }
    if let Some(r) = by(res1
        .max_delta_atom_atom_dist()
        .cmp(&res2.max_delta_atom_atom_dist()))
    {
        return Ok(r);
    }
    Ok(res1.smarts()? < res2.smarts()?)
}

/// The Johnson similarity of 2 molecules from their matches.
pub fn johnson_similarity(
    bond_matches: &[(i32, i32)],
    atom_matches: &[(i32, i32)],
    mol1: &ROMol,
    mol2: &ROMol,
) -> f64 {
    let n = (bond_matches.len() + atom_matches.len()) as f64;
    let num = n * n;
    let denom =
        f64::from((mol1.num_atoms() + mol1.num_bonds()) * (mol2.num_atoms() + mol2.num_bonds()));
    num / denom
}

/// Replaces every non-overlapping occurrence found by `find` (which gives the
/// length of the match at a position and the offset of the dash to drop).
fn drop_dashes(s: &str, find: impl Fn(&[u8], usize) -> Option<(usize, usize)>) -> String {
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        match find(b, i) {
            Some((len, dash)) => {
                out.extend_from_slice(&b[i..i + dash]);
                out.extend_from_slice(&b[i + dash + 1..i + len]);
                i += len;
            }
            None => {
                out.push(b[i]);
                i += 1;
            }
        }
    }
    String::from_utf8(out).expect("only ASCII bytes were dropped")
}

fn skip_digits(b: &[u8], mut i: usize) -> usize {
    while i < b.len() && (b'1'..=b'9').contains(&b[i]) {
        i += 1;
    }
    i
}

/// The replacements applied by `cleanSmarts` to make the SMARTS read like
/// SMILES.
const ATOM_REPLACEMENTS: &[(&str, &str)] = &[
    ("[#6&A]", "C"),
    ("[#6&A&R]", "[C&R]"),
    ("[#6&a]", "c"),
    ("[#7&A]", "N"),
    ("[#7&A&R]", "[N&R]"),
    ("[#7&a]", "n"),
    ("[#8&A]", "O"),
    ("[#8&A&R]", "[O&R]"),
    ("[#8&a]", "o"),
    ("[#9&A]", "F"),
    ("[#16&A]", "S"),
    ("[#16&a]", "s"),
    ("[#17&A]", "Cl"),
    ("[#35&A]", "Br"),
    ("[#53&A]", "I"),
];

/// Cleans up the SMARTS written for the MCES, and turns the placeholder atoms
/// of the equivalent atom classes (atomic numbers 110 up) back into the
/// original SMARTS.
pub fn clean_smarts(smarts: &mut String, equivalent_atoms: &str) {
    let upper = |c: u8| c.is_ascii_uppercase();
    let arom = |c: u8| b"cnops".contains(&c);
    // Sometimes it needs more than 1 pass through
    let mut start = String::new();
    while start != *smarts {
        start = smarts.clone();
        for (patt, repl) in ATOM_REPLACEMENTS {
            *smarts = smarts.replace(patt, repl);
        }
        // ([A-Z])-([cnops])
        *smarts = drop_dashes(smarts, |b, i| {
            (i + 2 < b.len() && upper(b[i]) && b[i + 1] == b'-' && arom(b[i + 2])).then_some((3, 1))
        });
        // ([cnops][1-9]*)-([A-Z])
        *smarts = drop_dashes(smarts, |b, i| {
            if !arom(b[i]) {
                return None;
            }
            let j = skip_digits(b, i + 1);
            (j + 1 < b.len() && b[j] == b'-' && upper(b[j + 1])).then_some((j + 2 - i, j - i))
        });
        // ([A-Z][1-9]*)-([A-Z])
        *smarts = drop_dashes(smarts, |b, i| {
            if !upper(b[i]) {
                return None;
            }
            let j = skip_digits(b, i + 1);
            (j + 1 < b.len() && b[j] == b'-' && upper(b[j + 1])).then_some((j + 2 - i, j - i))
        });
        // ([A-Z])-([1-9])
        *smarts = drop_dashes(smarts, |b, i| {
            (i + 2 < b.len()
                && upper(b[i])
                && b[i + 1] == b'-'
                && (b'1'..=b'9').contains(&b[i + 2]))
            .then_some((3, 1))
        });
    }

    // Convert the equivalent atoms from weird atomic numbers to the original
    // SMARTS pattern
    for (at_num, smt) in (110..).zip(equivalent_atoms.split(' ')) {
        // The SMARTS come out with &A or &a after the atomic number depending
        // on the aromaticity of the underlying atom but the original SMARTS
        // for the equivalent atom should take that into account. For example
        // [*] needs to match both aromatic and aliphatic atoms. Include the
        // case of no &[Aa] for good measure.
        for suffix in ["&A]", "&a]", "]"] {
            *smarts = smarts.replace(&format!("[#{at_num}{suffix}"), smt);
        }
        // If the ringMatchesRingOnly option has been used, there may be
        // [#110&A&R] or [#110&a&R] or [#110&R]. In these cases it needs to
        // end ;R] but again without the &A or &a. The SMARTS has to match to
        // a single atom, so must end with a ]. e.g. [O,S] or
        // [c,$(o1cccc1),$(n1cccc1)]
        let replace_with = format!("{};R", &smt[..smt.len().saturating_sub(1)]);
        let prefix = format!("[#{at_num}");
        let mut out = String::new();
        let mut rest = smarts.as_str();
        while let Some(pos) = rest.find(&prefix) {
            let mut j = pos + prefix.len();
            let b = rest.as_bytes();
            while j + 1 < b.len() && b[j] == b'&' && (b[j + 1] == b'A' || b[j + 1] == b'a') {
                j += 2;
            }
            if j + 1 < b.len() && b[j] == b'&' && b[j + 1] == b'R' {
                out.push_str(&rest[..pos]);
                out.push_str(&replace_with);
                rest = &rest[j + 2..];
            } else {
                out.push_str(&rest[..pos + prefix.len()]);
                rest = &rest[pos + prefix.len()..];
            }
        }
        out.push_str(rest);
        *smarts = out;
    }
}
