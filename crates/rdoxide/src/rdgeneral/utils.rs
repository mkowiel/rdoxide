//! Ported from RDKit `Code/RDGeneral/utils.h` (BSD-3-Clause, RDKit contributors).

use crate::{RDKitError, Result};

/// `countSwapsToInterconvert(ref, probe)`: number of swaps performed while
/// rearranging `probe` (a copy) to match `reference`, scanning left to right.
pub fn count_swaps_to_interconvert<T: PartialEq + Clone>(
    reference: &[T],
    probe: &[T],
) -> Result<u32> {
    if reference.len() != probe.len() {
        return Err(RDKitError::precondition("size mismatch"));
    }
    let mut probe = probe.to_vec();
    let mut n_swaps = 0;
    for i in 0..reference.len() {
        if probe[i] != reference[i] {
            let j = (i..probe.len())
                .find(|&j| probe[j] == reference[i])
                .ok_or_else(|| RDKitError::invariant("could not find probe element"))?;
            probe.swap(i, j);
            n_swaps += 1;
        }
    }
    Ok(n_swaps)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn swaps() {
        // examples from Atom::getPerturbationOrder docs: our order [0,1,2,3]
        let r = [0, 1, 2, 3];
        assert_eq!(count_swaps_to_interconvert(&[1, 0, 2, 3], &r).unwrap(), 1);
        assert_eq!(count_swaps_to_interconvert(&[1, 2, 3, 0], &r).unwrap(), 3);
        assert_eq!(count_swaps_to_interconvert(&[1, 2, 0, 3], &r).unwrap(), 2);
        assert!(count_swaps_to_interconvert(&[1, 2], &r).is_err());
    }
}
