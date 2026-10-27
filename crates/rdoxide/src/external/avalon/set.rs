//! `set.c`: bit sets with a fixed allocated size. Copyright 2010 Novartis Institutes for BioMedical
//! Research Inc.

/// `bit_set_t`: members `0..=max_member`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BitSet {
    bits: Vec<bool>,
}

impl BitSet {
    /// `NewSet`
    pub fn new(max_member: usize) -> Self {
        BitSet {
            bits: vec![false; max_member + 1],
        }
    }

    /// `MaxMember`: the allocated size
    pub fn max_member(&self) -> usize {
        self.bits.len() - 1
    }

    /// `IsMember`
    pub fn is_member(&self, m: usize) -> bool {
        self.bits.get(m).copied().unwrap_or(false)
    }

    /// `NextMember`
    pub fn next_member(&self, start: usize) -> Option<usize> {
        (start..self.bits.len()).find(|&i| self.bits[i])
    }

    /// `PutMember`: grows the set if needed
    pub fn put(&mut self, m: usize) {
        if m >= self.bits.len() {
            self.bits.resize(m + 1, false);
        }
        self.bits[m] = true;
    }

    /// `ClearSet`
    pub fn clear(&mut self) {
        self.bits.fill(false);
    }

    /// `CopySet`: the destination must not be smaller than the source; extra words are kept
    pub fn copy_from(&mut self, src: &BitSet) {
        assert!(
            self.bits.len() >= src.bits.len(),
            "destination set size < source set size"
        );
        self.bits[..src.bits.len()].copy_from_slice(&src.bits);
    }

    /// `SetExclusiveUnion`
    pub fn xor_with(&mut self, other: &BitSet) {
        assert!(self.bits.len() >= other.bits.len());
        for (a, b) in self.bits.iter_mut().zip(&other.bits) {
            *a ^= *b;
        }
    }

    /// `SetIntersection`
    pub fn intersect_with(&mut self, other: &BitSet) {
        assert!(self.bits.len() >= other.bits.len());
        for (i, a) in self.bits.iter_mut().enumerate() {
            *a &= other.bits.get(i).copied().unwrap_or(false);
        }
    }

    /// `IntersectionIsEmpty`
    pub fn intersection_is_empty(&self, other: &BitSet) -> bool {
        !self.bits.iter().zip(&other.bits).any(|(a, b)| *a && *b)
    }

    /// `Cardinality`
    pub fn cardinality(&self) -> usize {
        self.bits.iter().filter(|b| **b).count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_operations() {
        let mut a = BitSet::new(20);
        a.put(3);
        a.put(17);
        let mut b = BitSet::new(20);
        b.put(17);
        assert!(!a.intersection_is_empty(&b));
        a.xor_with(&b);
        assert_eq!(a.cardinality(), 1);
        assert_eq!(a.next_member(0), Some(3));
        assert!(a.intersection_is_empty(&b));
        a.put(40);
        assert_eq!(a.max_member(), 40);
    }
}
