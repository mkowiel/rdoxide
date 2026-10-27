//! `perceive.c`: classification of the bonds and atoms of a molecule. Copyright 2010 Novartis
//! Institutes for BioMedical Research Inc.

use super::graph::{combine_rings, ring_list};
use super::reaccs::Molecule;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::avalon::reaccsio::mol_str_to_mol;

    /// ring flags must match the Avalon C code on fused, bridged and caged ring systems
    #[test]
    fn ring_bonds_match_c() {
        let text = include_str!("../../../../../tests/data/avalon_rings.jsonl");
        let mut n = 0;
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().expect("input");
            let want = v["output"].as_str().expect("output");
            let mut m = mol_str_to_mol(input).expect("molecule");
            perceive_ring_bonds(&mut m);
            let got: String = m
                .bonds
                .iter()
                .map(|b| format!("{} {}\n", b.topography, b.dummy))
                .collect();
            assert_eq!(got, want, "{input}");
            n += 1;
        }
        assert!(n > 60);
    }
}
