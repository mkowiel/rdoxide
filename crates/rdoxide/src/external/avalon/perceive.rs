//! `perceive.c`: classification of the bonds and atoms of a molecule. Copyright 2010 Novartis
//! Institutes for BioMedical Research Inc.

use super::graph::{BondSet, combine_rings, proper_ring_pairs, ring_list};
use super::reaccs::{Molecule, SymbolList};
use super::set::BitSet;
use super::utilities::atom_symbol_match;

/// `RING`, `CHAIN` of `mdl.h`
pub const RING: i32 = 1;
pub const CHAIN: i32 = 2;

/// `PerceiveRingBonds`: sets the bond topography and flags the ring sizes in `dummy`
pub fn perceive_ring_bonds(mp: &mut Molecule) {
    let graph: Vec<[u32; 2]> = mp
        .bonds
        .iter()
        .map(|b| [b.atoms[0] as u32, b.atoms[1] as u32])
        .collect();
    for b in &mut mp.bonds {
        b.topography = CHAIN;
        b.dummy = 0;
    }
    let mut rings = ring_list(&graph);
    combine_rings(&mut rings);
    for ring in &rings {
        for (i, b) in mp.bonds.iter_mut().enumerate() {
            if ring.bond_set.is_member(i) {
                b.topography = RING;
                b.dummy |= 1 << ring.cardinality.min(8);
            }
        }
    }
}

/// bond types of `mdl.h`
pub const SINGLE: i32 = 1;
pub const DOUBLE: i32 = 2;
pub const TRIPLE: i32 = 3;
pub const AROMATIC: i32 = 4;

/// `neighbourhood_t`: the ligand atoms and bonds (indices, 0 based) of one atom
#[derive(Debug, Clone, Default)]
pub struct Neighbourhood {
    pub atoms: Vec<usize>,
    pub bonds: Vec<usize>,
}

/// `MAXNEIGHBOURS`
const MAXNEIGHBOURS: usize = 20;

/// `SetupNeighbourhood` for all atoms; `None` if an atom has too many ligands
pub fn setup_neighbourhood(mp: &Molecule) -> Option<Vec<Neighbourhood>> {
    let mut nb = vec![Neighbourhood::default(); mp.atoms.len()];
    for (i, b) in mp.bonds.iter().enumerate() {
        let (a0, a1) = ((b.atoms[0] - 1) as usize, (b.atoms[1] - 1) as usize);
        nb[a0].atoms.push(a1);
        nb[a0].bonds.push(i);
        if nb[a0].atoms.len() > MAXNEIGHBOURS {
            return None;
        }
        nb[a1].atoms.push(a0);
        nb[a1].bonds.push(i);
        if nb[a1].atoms.len() > MAXNEIGHBOURS {
            return None;
        }
    }
    Some(nb)
}

/// `RingState`: number of basis rings per bond and of ring bonds per atom
pub fn ring_state(mp: &Molecule) -> (Vec<i32>, Vec<i32>) {
    let mut atom_status = vec![0; mp.atoms.len()];
    let mut bond_status = vec![0; mp.bonds.len()];
    if mp.bonds.is_empty() {
        return (atom_status, bond_status);
    }
    let graph: Vec<[u32; 2]> = mp
        .bonds
        .iter()
        .map(|b| [b.atoms[0] as u32, b.atoms[1] as u32])
        .collect();
    let mut rings = ring_list(&graph);
    combine_rings(&mut rings);
    for r in &rings {
        for (i, st) in bond_status.iter_mut().enumerate() {
            if r.bond_set.is_member(i) {
                *st += 1;
            }
        }
    }
    for (i, b) in mp.bonds.iter().enumerate() {
        if bond_status[i] > 0 {
            atom_status[(b.atoms[0] - 1) as usize] += 1;
            atom_status[(b.atoms[1] - 1) as usize] += 1;
        }
    }
    (atom_status, bond_status)
}

struct Walk<'a> {
    nbp: &'a [Neighbourhood],
    max_size: usize,
    touched_atoms: Vec<bool>,
    touched_bonds: Vec<bool>,
}

/// `MarkRecursive`
fn mark_recursive(
    mp: &mut Molecule,
    w: &mut Walk,
    start_index: usize,
    path_length: usize,
    current_index: usize,
) {
    for i in 0..w.nbp[current_index].atoms.len() {
        let ai = w.nbp[current_index].atoms[i];
        if ai < start_index {
            continue;
        }
        if ai == start_index {
            if path_length < 3 {
                continue;
            }
            for (j, a) in mp.atoms.iter_mut().enumerate() {
                if w.touched_atoms[j] {
                    a.rsize_flags |= 1 << path_length;
                }
            }
            for (j, b) in mp.bonds.iter_mut().enumerate() {
                if w.touched_bonds[j] {
                    b.rsize_flags |= 1 << path_length;
                }
            }
            continue;
        }
        if w.touched_atoms[ai] {
            continue;
        }
        if path_length + 1 > w.max_size {
            continue;
        }
        if mp.atoms[ai].rsize_flags == 0 {
            continue;
        }
        let bi = w.nbp[current_index].bonds[i];
        if mp.bonds[bi].rsize_flags == 0 {
            continue;
        }
        w.touched_atoms[ai] = true;
        w.touched_bonds[bi] = true;
        mark_recursive(mp, w, start_index, path_length + 1, ai);
        w.touched_atoms[ai] = false;
        w.touched_bonds[bi] = false;
    }
}

/// `SetRingSizeFlags`: ring sizes up to `max_size` as flag bits in `rsize_flags` (bit 0: ring member)
pub fn set_ring_size_flags(mp: &mut Molecule, max_size: usize, nbp: &[Neighbourhood]) {
    let (atom_status, bond_status) = ring_state(mp);
    for (a, st) in mp.atoms.iter_mut().zip(&atom_status) {
        a.rsize_flags = u32::from(*st > 0);
    }
    for (b, st) in mp.bonds.iter_mut().zip(&bond_status) {
        b.rsize_flags = u32::from(*st > 0);
    }
    let mut w = Walk {
        nbp,
        max_size,
        touched_atoms: vec![false; mp.atoms.len()],
        touched_bonds: vec![false; mp.bonds.len()],
    };
    for i in 0..mp.atoms.len() {
        if mp.atoms[i].rsize_flags == 0 {
            continue;
        }
        w.touched_atoms[i] = true;
        mark_recursive(mp, &mut w, i, 1, i);
        w.touched_atoms[i] = false;
    }
}

fn graph_of(mp: &Molecule) -> Vec<[u32; 2]> {
    mp.bonds
        .iter()
        .map(|b| [b.atoms[0] as u32, b.atoms[1] as u32])
        .collect()
}

/// the C idiom that moves the nodes of `pairs` one by one in front of `ring_list`
fn prepend_all(ring_list: &mut Vec<BondSet>, pairs: Vec<BondSet>) {
    for p in pairs {
        ring_list.insert(0, p);
    }
}

/// `PerceiveAromaticBonds`: ring bonds of rings of sp2 atoms of Hueckel size become aromatic
pub fn perceive_aromatic_bonds(mp: &mut Molecule) {
    let graph = graph_of(mp);
    let mut bond_is_in_ring = vec![false; mp.bonds.len()];
    let mut rings = ring_list(&graph);
    combine_rings(&mut rings);
    for r in &rings {
        for (i, f) in bond_is_in_ring.iter_mut().enumerate() {
            if r.bond_set.is_member(i) {
                *f = true;
            }
        }
    }
    let pairs = proper_ring_pairs(&rings, mp.atoms.len(), &graph);
    prepend_all(&mut rings, pairs);

    loop {
        let mut changed = false;
        for r in &rings {
            let members: Vec<usize> = (0..mp.bonds.len())
                .filter(|&i| r.bond_set.is_member(i))
                .collect();
            let mut sp_count = vec![0; mp.atoms.len() + 1];
            let (mut ndouble, mut nsingle) = (0, 0);
            for &i in &members {
                if mp.bonds[i].bond_type == SINGLE {
                    nsingle += 1;
                }
            }
            let mut is_cumulene = false;
            for &i in &members {
                let b = &mp.bonds[i];
                if b.bond_type == DOUBLE {
                    ndouble += 1;
                    for k in 0..2 {
                        let a = b.atoms[k] as usize;
                        sp_count[a] += 1;
                        if sp_count[a] > 1 {
                            is_cumulene = true;
                        }
                    }
                } else if b.bond_type == TRIPLE {
                    is_cumulene = true;
                }
            }
            for &i in &members {
                let b = &mp.bonds[i];
                if b.bond_type == AROMATIC {
                    for k in 0..2 {
                        let a = b.atoms[k] as usize;
                        if sp_count[a] == 0 {
                            sp_count[a] += 1;
                        }
                    }
                }
            }
            let mut aromatic = !is_cumulene && r.cardinality >= 2 && (r.cardinality - 2) % 4 == 0;
            for &i in &members {
                let b = &mp.bonds[i];
                if sp_count[b.atoms[0] as usize] != 1 || sp_count[b.atoms[1] as usize] != 1 {
                    aromatic = false;
                }
            }
            if aromatic && (ndouble > 0 || nsingle > 0) {
                for &i in &members {
                    if bond_is_in_ring[i] && mp.bonds[i].bond_type != AROMATIC {
                        changed = true;
                        mp.bonds[i].bond_type = AROMATIC;
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }
}

/// `getList`
pub fn get_list(mp: &Molecule, index: usize) -> Option<&SymbolList> {
    mp.symbol_lists.iter().find(|l| l.atom == index as i32 + 1)
}

/// `PerceiveDYAromaticity`: Daylight style aromaticity; remembers hydrogens of aromatic hetero
/// atoms in `query_h_count`
pub fn perceive_dy_aromaticity(mp: &mut Molecule, nbp: &[Neighbourhood]) {
    let graph0 = graph_of(mp);
    let mut bond_is_in_ring = vec![false; mp.bonds.len()];
    let rings = ring_list(&graph0);
    if rings.is_empty() {
        return;
    }
    let mut atom_is_in_ring = vec![false; mp.atoms.len()];
    let mut candidate: Vec<bool> = mp.atoms.iter().map(|a| a.atom_symbol != "C").collect();
    for r in &rings {
        for (i, f) in bond_is_in_ring.iter_mut().enumerate() {
            if r.bond_set.is_member(i) {
                *f = true;
            }
        }
    }
    for b in &mp.bonds {
        if b.bond_type > SINGLE && b.bond_type != TRIPLE {
            candidate[(b.atoms[0] - 1) as usize] = true;
            candidate[(b.atoms[1] - 1) as usize] = true;
        }
    }
    let mut graph = graph0;
    for i in 0..mp.bonds.len() {
        if !bond_is_in_ring[i] {
            graph[i] = [0, 0];
        } else {
            atom_is_in_ring[graph[i][0] as usize - 1] = true;
            atom_is_in_ring[graph[i][1] as usize - 1] = true;
            if !candidate[graph[i][0] as usize - 1] {
                graph[i][0] = 0;
            }
            if !candidate[graph[i][1] as usize - 1] {
                graph[i][1] = 0;
            }
        }
    }
    let mut ring_list_v = ring_list(&graph);
    if ring_list_v.is_empty() {
        return;
    }
    combine_rings(&mut ring_list_v);
    let pairs = proper_ring_pairs(&ring_list_v, mp.atoms.len(), &graph);
    prepend_all(&mut ring_list_v, pairs);

    // fused pairs that could be aromatic as well; they are tested before the plain rings
    let mut set = BitSet::new(ring_list_v[0].bond_set.max_member());
    let mut fused: Vec<BondSet> = Vec::new();
    for i in 0..ring_list_v.len() {
        for j in i + 1..ring_list_v.len() {
            if ring_list_v[i]
                .bond_set
                .intersection_is_empty(&ring_list_v[j].bond_set)
            {
                continue;
            }
            set.copy_from(&ring_list_v[i].bond_set);
            set.xor_with(&ring_list_v[j].bond_set);
            if set.cardinality() + 2 == ring_list_v[i].cardinality + ring_list_v[j].cardinality {
                fused.insert(
                    0,
                    BondSet {
                        cardinality: set.cardinality(),
                        bond_set: set.clone(),
                    },
                );
            }
        }
    }
    fused.extend(ring_list_v);
    let candidates = fused;

    loop {
        let mut changed = false;
        for r in &candidates {
            let mut npi = 0;
            let mut conjugated = true;
            for i in 0..mp.atoms.len() {
                if !atom_is_in_ring[i] {
                    continue;
                }
                let (mut in_ring_double, mut in_ring_aromatic) = (0, 0);
                let mut exo_pull = false;
                let mut is_in_ring = false;
                let mut local_pi = 0;
                for (&aj, &bj) in nbp[i].atoms.iter().zip(&nbp[i].bonds) {
                    let bt = mp.bonds[bj].bond_type;
                    if bond_is_in_ring[bj] && r.bond_set.is_member(bj) {
                        is_in_ring = true;
                        if bt == AROMATIC {
                            in_ring_aromatic += 1;
                        }
                        if bt == DOUBLE {
                            in_ring_double += 1;
                        }
                    } else if bt == DOUBLE
                        && mp.atoms[i].atom_symbol == "C"
                        && !bond_is_in_ring[bj]
                        && atom_symbol_match(&mp.atoms[aj].atom_symbol, "O,S,P,N,L")
                    {
                        exo_pull = true;
                    }
                }
                if !is_in_ring {
                    continue;
                }
                let sym = mp.atoms[i].atom_symbol.clone();
                let charge = mp.atoms[i].charge;
                let lst_ok =
                    get_list(mp, i).is_some_and(|l| l.logic != 0 && !l.string.contains('C'));
                if (in_ring_aromatic >= 1 || in_ring_double == 1)
                    && (atom_symbol_match(&sym, "C,N,A,*") || sym == "L")
                {
                    local_pi = 1;
                } else if in_ring_aromatic == 0
                    && in_ring_double == 0
                    && charge == 0
                    && (atom_symbol_match(&sym, "N,S,O")
                        || (atom_symbol_match(&sym, "L") && lst_ok))
                {
                    local_pi = 2;
                    if nbp[i].atoms.len() == 2 && sym == "N" {
                        mp.atoms[i].query_h_count = 1 + 1;
                    }
                } else if in_ring_aromatic == 0
                    && in_ring_double == 0
                    && charge == 0
                    && exo_pull
                    && sym == "C"
                {
                    local_pi = 0;
                } else {
                    conjugated = false;
                }
                if charge < 0 && atom_symbol_match(&sym, "C,N") {
                    conjugated = false;
                }
                npi += local_pi;
            }
            if !conjugated || npi % 4 != 2 {
                continue;
            }
            for i in 0..mp.bonds.len() {
                if bond_is_in_ring[i]
                    && r.bond_set.is_member(i)
                    && mp.bonds[i].bond_type != AROMATIC
                {
                    mp.bonds[i].bond_type = AROMATIC;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }

    let mut arom = vec![false; mp.atoms.len()];
    for b in &mp.bonds {
        if b.bond_type == AROMATIC {
            arom[(b.atoms[0] - 1) as usize] = true;
            arom[(b.atoms[1] - 1) as usize] = true;
        }
    }
    for (i, b) in mp.bonds.iter_mut().enumerate() {
        if arom[(b.atoms[0] - 1) as usize]
            && arom[(b.atoms[1] - 1) as usize]
            && bond_is_in_ring[i]
            && b.bond_type == SINGLE
        {
            b.bond_type = AROMATIC;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::avalon::denormal::cct_denormalize;
    use crate::external::avalon::reaccsio::mol_str_to_mol;

    /// every perception routine must leave the same molecule as the Avalon C code
    #[test]
    fn perception_matches_c() {
        let text = include_str!("../../../../../tests/data/avalon_rings.jsonl");
        let mut n = 0;
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().expect("input");
            let mode = v["mode"].as_str().expect("mode");
            let want = v["output"].as_str().expect("output");
            let mut m = mol_str_to_mol(input).expect("molecule");
            let nbp = setup_neighbourhood(&m).expect("neighbourhood");
            let mut got = String::new();
            match mode {
                "ring" => perceive_ring_bonds(&mut m),
                "arom" => perceive_aromatic_bonds(&mut m),
                "dy" => perceive_dy_aromaticity(&mut m, &nbp),
                "rsize" => set_ring_size_flags(&mut m, 8, &nbp),
                "den" => got += &format!("r {}\n", cct_denormalize(&mut m)),
                _ => unreachable!(),
            }
            for a in &m.atoms {
                got += &format!(
                    "a {} {} {} {} {}\n",
                    a.query_h_count, a.rsize_flags, a.charge, a.dummy2, a.color
                );
            }
            for b in &m.bonds {
                got += &format!(
                    "{} {} {} {}\n",
                    b.bond_type, b.topography, b.dummy, b.rsize_flags
                );
            }
            assert_eq!(got, want, "{mode}\n{input}");
            n += 1;
        }
        assert!(n > 1300);
    }
}
