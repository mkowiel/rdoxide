//! Types shared by the bounds-matrix builder (upstream
//! `BoundsMatrixBuilderDetails.h`).

use std::f64::consts::PI;

use crate::distgeom::{TorsionCandidates, TorsionRange, TorsionValues};
use crate::error::{RDKitError, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TorsionType {
    Cis,
    Trans,
    Flexible,
    Custom,
    CisTrans,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Type14 {
    InChain,
    InRing,
    TwoInSameRing,
    TwoInDiffRing,
    ShareRingBond,
    MacrocycleTwoInSameRing,
    MacrocycleAllInSameRing,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TorsionValue {
    pub kind: TorsionType,
    pub value: Option<f64>,
    pub extra_dist: Option<f64>,
    pub is_forced: bool,
}

impl TorsionValue {
    pub fn of(kind: TorsionType) -> Self {
        TorsionValue {
            kind,
            value: None,
            extra_dist: None,
            is_forced: false,
        }
    }

    pub fn forced(kind: TorsionType) -> Self {
        TorsionValue {
            is_forced: true,
            ..Self::of(kind)
        }
    }
}

/// The torsion range allowed for a bond in a ring of the given size.
pub fn ring_torsion(r_size: usize) -> TorsionRange {
    let torsion = match r_size {
        4 | 5 => PI * 45.0 / 180.0,
        6 => PI * 60.0 / 180.0,
        7 => PI * 90.0 / 180.0,
        8 => PI * 100.0 / 180.0,
        _ => PI,
    };
    TorsionRange::new(-torsion, torsion).expect("ring torsion bounds are ordered")
}

/// A 1-4 path with its torsion type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Path14Configuration {
    pub bid1: u32,
    pub bid2: u32,
    pub bid3: u32,
    pub aid1: u32,
    pub aid2: u32,
    pub aid3: u32,
    pub aid4: u32,
    pub value: TorsionValue,
    pub type14: Type14,
    pub r_size: usize,
}

impl Path14Configuration {
    pub fn to_torsion_range(&self) -> TorsionCandidates {
        match self.value.kind {
            TorsionType::Cis => TorsionValues::new(&[0.0]).into(),
            TorsionType::Trans => TorsionValues::new(&[PI]).into(),
            TorsionType::Flexible => ring_torsion(self.r_size).into(),
            TorsionType::CisTrans => TorsionValues::new(&[0.0, PI]).into(),
            TorsionType::Custom => {
                TorsionValues::new(&[self.value.value.expect("custom torsion value")]).into()
            }
        }
    }
}

pub type Path14Vect = Vec<Path14Configuration>;

/// An id for the unordered pair (id1, id2), with both ids below `n`.
pub fn unified_id2(id1: u32, id2: u32, n: u32) -> usize {
    let (lo, hi) = if id1 < id2 { (id1, id2) } else { (id2, id1) };
    lo as usize * n as usize + hi as usize
}

/// An id for (id1, id2, id3) independent of the order of id1 and id3, with all
/// ids below `n`.
pub fn unified_id3(id1: u32, id2: u32, id3: u32, n: u32) -> usize {
    let n = n as usize;
    if id1 < id3 {
        id1 as usize * n * n + id2 as usize * n + id3 as usize
    } else {
        id3 as usize * n * n + id2 as usize * n + id1 as usize
    }
}

/// Converts an id made by [`unified_id2`] or [`unified_id3`] back into the
/// original ids (`N` of them).
pub fn unified_id_to_bond_ids<const N: usize>(mut id: usize, n: u32) -> [u32; N] {
    let mut ids = [0u32; N];
    for slot in ids.iter_mut().rev() {
        *slot = (id % n as usize) as u32;
        id /= n as usize;
    }
    ids
}

/// A range of distances between the atoms aid1 and aid4. Starts out invalid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bounds {
    pub lower: f64,
    pub upper: f64,
    pub aid1: u32,
    pub aid4: u32,
}

impl Default for Bounds {
    fn default() -> Self {
        Bounds {
            lower: 1.0,
            upper: -1.0,
            aid1: 0,
            aid4: 0,
        }
    }
}

impl Bounds {
    pub fn new(lower: f64, upper: f64, aid1: u32, aid4: u32) -> Self {
        Bounds {
            lower,
            upper,
            aid1,
            aid4,
        }
    }

    pub fn valid(&self) -> bool {
        self.lower <= self.upper
    }
}

/// Merges candidate bounds into one: the union of all maximal non-empty
/// intersections is approximated by its first lower and last upper bound.
pub fn merge_bounds(mut bounds: Vec<Bounds>) -> Result<Bounds> {
    if bounds.is_empty() {
        return Err(RDKitError::precondition(
            "Cannot merge empty list of bounds",
        ));
    }
    bounds.sort_by(|a, b| a.lower.total_cmp(&b.lower));

    // sweep over the bounds sorted by their lower value, tracking the active
    // non-empty intersection and the largest upper bound seen so far (the
    // active intersection can end earlier and would lose potential overlaps)
    let mut current = bounds[0];
    let mut component_upper = current.upper;
    let mut result_lower: Option<f64> = None;
    for b in &bounds[1..] {
        if b.lower <= current.upper {
            current.lower = b.lower;
            current.upper = current.upper.min(b.upper);
        } else {
            // the active intersection is maximal
            result_lower.get_or_insert(current.lower);
            current.lower = b.lower;
            current.upper = if b.lower <= component_upper {
                component_upper.min(b.upper)
            } else {
                b.upper
            };
        }
        component_upper = component_upper.max(b.upper);
    }
    Ok(Bounds {
        lower: result_lower.unwrap_or(current.lower),
        upper: current.upper,
        aid1: bounds[0].aid1,
        aid4: bounds[0].aid4,
    })
}

/// A growable set of bits.
#[derive(Debug, Clone)]
pub(super) struct BitSet {
    words: Vec<u64>,
}

impl BitSet {
    pub fn new(len: usize) -> Self {
        BitSet {
            words: vec![0; len.div_ceil(64)],
        }
    }

    pub fn set(&mut self, i: usize) {
        self.words[i / 64] |= 1 << (i % 64);
    }

    pub fn set_to(&mut self, i: usize, value: bool) {
        if value {
            self.set(i);
        } else {
            self.words[i / 64] &= !(1 << (i % 64));
        }
    }

    pub fn get(&self, i: usize) -> bool {
        self.words[i / 64] >> (i % 64) & 1 == 1
    }
}

/// A symmetric matrix of values for pairs of bonds.
#[derive(Debug, Clone)]
pub(super) struct SymmGrid<T: Copy> {
    n: usize,
    data: Vec<T>,
}

impl<T: Copy> SymmGrid<T> {
    pub fn new(n: usize, fill: T) -> Self {
        SymmGrid {
            n,
            data: vec![fill; n * (n + 1) / 2],
        }
    }

    pub fn num_rows(&self) -> usize {
        self.n
    }

    fn index(i: usize, j: usize) -> usize {
        let (hi, lo) = if i >= j { (i, j) } else { (j, i) };
        hi * (hi + 1) / 2 + lo
    }

    pub fn get(&self, i: u32, j: u32) -> T {
        self.data[Self::index(i as usize, j as usize)]
    }

    pub fn set(&mut self, i: u32, j: u32, v: T) {
        self.data[Self::index(i as usize, j as usize)] = v;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unified_ids_round_trip() {
        assert_eq!(unified_id2(3, 7, 10), unified_id2(7, 3, 10));
        assert_eq!(
            unified_id_to_bond_ids::<2>(unified_id2(7, 3, 10), 10),
            [3, 7]
        );
        assert_eq!(unified_id3(1, 5, 9, 10), unified_id3(9, 5, 1, 10));
        assert_eq!(
            unified_id_to_bond_ids::<3>(unified_id3(9, 5, 1, 10), 10),
            [1, 5, 9]
        );
    }

    #[test]
    fn merge_keeps_first_lower_and_last_upper() {
        let b = |l, u| Bounds::new(l, u, 0, 3);
        // overlapping: the intersection
        let m = merge_bounds(vec![b(2.0, 4.0), b(3.0, 5.0)]).unwrap();
        assert_eq!((m.lower, m.upper), (3.0, 4.0));
        // disjoint: first lower of the first intersection, last upper
        let m = merge_bounds(vec![b(1.0, 2.0), b(3.0, 4.0)]).unwrap();
        assert_eq!((m.lower, m.upper), (1.0, 4.0));
        assert!(merge_bounds(Vec::new()).is_err());
    }

    #[test]
    fn bit_set_set_and_clear() {
        let mut bits = BitSet::new(130);
        bits.set(129);
        bits.set_to(3, true);
        bits.set_to(3, false);
        assert!(bits.get(129) && !bits.get(3) && !bits.get(128));
    }
}
