//! `SynthonSpaceHitSet.cpp`: the synthons of one SynthonSet that might make
//! hits for a fragment set.

use std::sync::Arc;

use super::synthon::Synthon;
use super::synthon_set::SynthonSet;
use crate::datastructs::bit_vect::ExplicitBitVect;
use crate::graphmol::romol::ROMol;

/// The extra data of the hit sets of the fingerprint searcher.
pub struct FpHitData<'a> {
    pub synthon_fps: Vec<Vec<&'a ExplicitBitVect>>,
    pub add_fp: Option<&'a ExplicitBitVect>,
    pub subtract_fp: Option<&'a ExplicitBitVect>,
}

pub struct SynthonSpaceHitSet<'a> {
    pub reaction: &'a SynthonSet,
    /// `synthons_to_use[i]` is the selection from `reaction.synthons()[i]`:
    /// (synthon id, synthon).
    pub synthons_to_use: Vec<Vec<(&'a str, &'a Synthon)>>,
    pub num_hits: usize,
    /// The fragments this hitset is derived from.
    pub frags: Vec<Arc<ROMol>>,
    pub fp: Option<FpHitData<'a>>,
}

impl<'a> SynthonSpaceHitSet<'a> {
    pub fn new(
        reaction: &'a SynthonSet,
        pool: &'a [Synthon],
        stu: &[Vec<usize>],
        frag_set: &[Arc<ROMol>],
    ) -> Self {
        let synthons = reaction.synthons();
        let synthons_to_use: Vec<Vec<(&'a str, &'a Synthon)>> = stu
            .iter()
            .enumerate()
            .map(|(i, sel)| {
                sel.iter()
                    .map(|&j| {
                        let (id, idx) = &synthons[i][j];
                        (id.as_str(), &pool[*idx])
                    })
                    .collect()
            })
            .collect();
        let num_hits = stu.iter().fold(1usize, |acc, s| acc * s.len());
        SynthonSpaceHitSet {
            reaction,
            synthons_to_use,
            num_hits,
            frags: frag_set.to_vec(),
            fp: None,
        }
    }

    /// A hit set that also holds the synthon fingerprints.
    pub fn new_fp(
        reaction: &'a SynthonSet,
        pool: &'a [Synthon],
        stu: &[Vec<usize>],
        frag_set: &[Arc<ROMol>],
    ) -> Self {
        let mut hs = Self::new(reaction, pool, stu, frag_set);
        let synthons = reaction.synthons();
        let synthon_fps = stu
            .iter()
            .enumerate()
            .map(|(i, sel)| {
                sel.iter()
                    .filter_map(|&j| pool[synthons[i][j].1].fp())
                    .collect()
            })
            .collect();
        hs.fp = Some(FpHitData {
            synthon_fps,
            add_fp: reaction.add_fp(),
            subtract_fp: reaction.subtract_fp(),
        });
        hs
    }
}
