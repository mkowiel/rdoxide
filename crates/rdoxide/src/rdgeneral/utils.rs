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

/// `Union(rings, res)` (RDGeneral/types.cpp): ordered union of ring members.
pub fn union_rings(rings: &[Vec<i32>]) -> Vec<i32> {
    let mut res: Vec<i32> = Vec::new();
    for ring in rings {
        for &x in ring {
            if !res.contains(&x) {
                res.push(x);
            }
        }
    }
    res
}

/// `Intersect(r1, r2, res)`: elements of `r1` also in `r2`, in `r1` order.
pub fn intersect(r1: &[i32], r2: &[i32]) -> Vec<i32> {
    r1.iter().filter(|x| r2.contains(x)).copied().collect()
}

/// `nextCombination(comb, tot)`: advance to the next k-combination of
/// `0..tot`; returns the index of the first changed position or -1 when done.
pub fn next_combination(comb: &mut [i32], tot: i32) -> i32 {
    let nelem = comb.len() as i32;
    let mut celem = nelem - 1;
    while comb[celem as usize] == tot - nelem + celem {
        celem -= 1;
        if celem < 0 {
            return -1;
        }
    }
    comb[celem as usize] += 1;
    for i in (celem as usize + 1)..comb.len() {
        comb[i] = comb[i - 1] + 1;
    }
    celem
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn combinations() {
        let mut c = vec![0, 1];
        let mut seen = vec![c.clone()];
        while next_combination(&mut c, 4) != -1 {
            seen.push(c.clone());
        }
        assert_eq!(
            seen,
            vec![
                vec![0, 1],
                vec![0, 2],
                vec![0, 3],
                vec![1, 2],
                vec![1, 3],
                vec![2, 3]
            ]
        );
        assert_eq!(
            union_rings(&[vec![1, 2, 3], vec![3, 4, 1]]),
            vec![1, 2, 3, 4]
        );
        assert_eq!(intersect(&[1, 2, 3], &[3, 1]), vec![1, 3]);
    }

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
