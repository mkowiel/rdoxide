//! The MaxMin diversity picker (upstream `SimDivPickers/MaxMinPicker.h`).
//!
//! Ashton, M. et al., Quant. Struct.-Act. Relat., 21 (2002), 598-604. A
//! subset is grown one item at a time, each time adding the pool item whose
//! smallest distance to the items picked so far is largest. Distances are
//! evaluated lazily: an item that cannot beat the best candidate of the
//! current round is not compared with the remaining picks.

use super::{check_dist_mat, check_pool, dist_from_ltm, pool_index};
use crate::error::Result;
use crate::rdgeneral::random::Mt19937;

/// Result of a MaxMin pick.
#[derive(Debug, Clone, PartialEq)]
pub struct MaxMinPicks {
    /// Picked items, in pick order (the first picks come first).
    pub picks: Vec<i32>,
    /// The distance of the last pick to its nearest earlier pick (`-1` if
    /// nothing had to be picked).
    pub threshold: f64,
}

/// Per pool item bookkeeping of the lazy evaluation.
#[derive(Debug, Clone, Copy, Default)]
struct Candidate {
    /// An upper bound of the distance to the nearest pick compared so far.
    dist_bound: f64,
    /// How many of the picks `dist_bound` takes into account.
    n_compared: usize,
}

/// A random seed for `seed < 0` (upstream uses `std::random_device`).
fn random_seed() -> u32 {
    use std::hash::{BuildHasher, Hasher};
    std::collections::hash_map::RandomState::new()
        .build_hasher()
        .finish() as u32
}

/// `MaxMinPicker::lazyPick`: picks `pick_size` of `pool_size` items.
///
/// * `dist` - the distance between two pool items.
/// * `first_picks` - items to start from; without them the first item is
///   chosen at random (`seed >= 0` makes this reproducible).
/// * `threshold` - with a value `>= 0`, picking stops as soon as the best
///   remaining distance is not larger than it.
pub fn lazy_pick(
    mut dist: impl FnMut(u32, u32) -> Result<f64>,
    pool_size: u32,
    pick_size: u32,
    first_picks: &[i32],
    seed: i32,
    threshold: f64,
) -> Result<MaxMinPicks> {
    check_pool(pool_size, pick_size)?;
    let mut candidates = vec![Candidate::default(); pool_size as usize];
    let mut picked = vec![false; pool_size as usize];
    let mut picks: Vec<i32> = Vec::with_capacity(pick_size as usize);
    let mut n_picked = 0usize;

    if first_picks.is_empty() {
        let mut rng = Mt19937::new(if seed >= 0 {
            seed as u32
        } else {
            random_seed()
        });
        let first = rng.uniform_int(0, pool_size as i32 - 1);
        picks.push(first);
        picked[first as usize] = true;
        n_picked = 1;
    } else {
        for &p in first_picks {
            let p = pool_index(p, pool_size)?;
            picks.push(p as i32);
            picked[p as usize] = true;
            n_picked += 1;
        }
    }
    let done = |res: Vec<i32>, threshold| MaxMinPicks {
        picks: res,
        threshold,
    };
    if n_picked >= pick_size as usize {
        return Ok(done(picks, -1.0));
    }

    // the items still to be picked, in index order
    let mut pool: Vec<u32> = (0..pool_size).filter(|&i| !picked[i as usize]).collect();
    if pool.is_empty() {
        return Ok(done(picks, -1.0));
    }

    // first pass: compare everything with the first pick
    let first = picks[0] as u32;
    for &i in &pool {
        candidates[i as usize] = Candidate {
            dist_bound: dist(i, first)?,
            n_compared: 1,
        };
    }

    let mut last_distance = -1.0;
    while n_picked < pick_size as usize && !pool.is_empty() {
        let mut max_of_min = -1.0;
        let mut best: Option<usize> = None; // position in `pool`
        for (pos, &item) in pool.iter().enumerate() {
            let cand = &mut candidates[item as usize];
            let mut min_to_i = cand.dist_bound;
            if min_to_i > max_of_min {
                let mut compared = cand.n_compared;
                while compared < n_picked {
                    let other = picks[compared] as u32;
                    debug_assert_ne!(item, other, "pool index != pick index");
                    let d = dist(item, other)?;
                    compared += 1;
                    if d <= min_to_i {
                        min_to_i = d;
                        if min_to_i <= max_of_min {
                            break;
                        }
                    }
                }
                cand.dist_bound = min_to_i;
                cand.n_compared = compared;
                if min_to_i > max_of_min {
                    max_of_min = min_to_i;
                    best = Some(pos);
                }
            }
        }
        let Some(best) = best else { break };
        if max_of_min <= threshold && threshold >= 0.0 {
            break;
        }
        last_distance = max_of_min;
        picks.push(pool.remove(best) as i32);
        n_picked += 1;
    }
    Ok(done(picks, last_distance))
}

/// `MaxMinPicker::pick`: MaxMin picking from a packed distance matrix (see
/// [`dist_from_ltm`]).
pub fn pick(
    dist_mat: &[f64],
    pool_size: u32,
    pick_size: u32,
    first_picks: &[i32],
    seed: i32,
) -> Result<Vec<i32>> {
    check_pool(pool_size, pick_size)?;
    check_dist_mat(dist_mat, pool_size)?;
    let res = lazy_pick(
        |i, j| Ok(dist_from_ltm(dist_mat, i, j)),
        pool_size,
        pick_size,
        first_picks,
        seed,
        -1.0,
    )?;
    Ok(res.picks)
}
