//! Roger Sayle's Leader picker (upstream `SimDivPickers/LeaderPicker.h`),
//! analogous to sphere exclusion: the first remaining item is picked and
//! every item within `threshold` of it is dropped from the pool.
//!
//! With several threads the distances of one round are computed in parallel;
//! the result does not depend on the thread count.

use super::{check_pool, pool_index};
use crate::error::Result;
use crate::rdgeneral::threads::{num_threads_to_use, par_map};

/// Smallest pool for which a round is split over threads.
const PARALLEL_MIN_POOL: usize = 16 * 1024;

/// The pool items farther than `threshold` from `query`.
fn keep_far<F>(items: &[u32], query: u32, threshold: f64, dist: &F) -> Result<Vec<u32>>
where
    F: Fn(u32, u32) -> Result<f64> + Sync,
{
    let mut kept = Vec::with_capacity(items.len());
    for &item in items {
        if dist(query, item)? > threshold {
            kept.push(item);
        }
    }
    Ok(kept)
}

/// Drops from `pool` the items within `threshold` of `query`.
fn exclude_near<F>(
    pool: &mut Vec<u32>,
    query: u32,
    threshold: f64,
    dist: &F,
    threads: usize,
) -> Result<()>
where
    F: Fn(u32, u32) -> Result<f64> + Sync,
{
    if threads <= 1 || pool.len() < PARALLEL_MIN_POOL {
        *pool = keep_far(pool, query, threshold, dist)?;
        return Ok(());
    }
    let chunks: Vec<&[u32]> = pool.chunks(pool.len().div_ceil(threads)).collect();
    let parts = par_map(&chunks, threads as i32, |chunk| {
        keep_far(chunk, query, threshold, dist)
    });
    let mut kept = Vec::with_capacity(pool.len());
    for part in parts {
        kept.extend(part?);
    }
    *pool = kept;
    Ok(())
}

/// `LeaderPicker::lazyPick`
///
/// * `dist` - the distance between two pool items.
/// * `pick_size` - the number of items to pick (`0`: as many as there are).
/// * `first_picks` - items that are picked first (and exclude their
///   neighbors from the pool).
/// * `threshold` - items at this distance or less from a pick are excluded.
/// * `num_threads` - threads to use; zero or negative means all hardware
///   threads minus `|num_threads|`.
pub fn lazy_pick<F>(
    dist: F,
    pool_size: u32,
    pick_size: u32,
    first_picks: &[i32],
    threshold: f64,
    num_threads: i32,
) -> Result<Vec<i32>>
where
    F: Fn(u32, u32) -> Result<f64> + Sync,
{
    check_pool(pool_size, pick_size)?;
    let pick_size = if pick_size == 0 { pool_size } else { pick_size } as usize;
    let threads = num_threads_to_use(num_threads);
    let mut pool: Vec<u32> = (0..pool_size).collect();
    let mut picks = Vec::new();

    for &p in first_picks {
        let p = pool_index(p, pool_size)?;
        picks.push(p as i32);
        exclude_near(&mut pool, p, threshold, &dist, threads)?;
    }
    while picks.len() < pick_size && !pool.is_empty() {
        // the pick leaves the pool whatever its distance to itself
        let query = pool.remove(0);
        picks.push(query as i32);
        exclude_near(&mut pool, query, threshold, &dist, threads)?;
    }
    Ok(picks)
}
