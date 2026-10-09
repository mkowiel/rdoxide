//! Diversity pickers (upstream `Code/SimDivPickers`): MaxMin and Leader.
//!
//! Both pickers work on a pool of `pool_size` items identified by their
//! index, and only need the distance between two items, which a caller
//! supplies as a function (or as a packed distance matrix, see
//! [`dist_from_ltm`]).

pub mod leader;
pub mod max_min;

use crate::error::{RDKitError, Result};

/// `getDistFromLTM`: the distance between items `i` and `j` from the lower
/// triangle of a distance matrix stored row by row without the diagonal.
///
/// `dist_mat` must hold `pool_size * (pool_size - 1) / 2` values, which
/// [`check_dist_mat`] verifies once up front.
pub fn dist_from_ltm(dist_mat: &[f64], i: u32, j: u32) -> f64 {
    let (hi, lo) = match i.cmp(&j) {
        std::cmp::Ordering::Equal => return 0.0,
        std::cmp::Ordering::Greater => (i as usize, j as usize),
        std::cmp::Ordering::Less => (j as usize, i as usize),
    };
    dist_mat[hi * (hi - 1) / 2 + lo]
}

/// Checks that `dist_mat` is large enough for a pool of `pool_size` items.
pub fn check_dist_mat(dist_mat: &[f64], pool_size: u32) -> Result<()> {
    let n = pool_size as usize;
    if dist_mat.len() < n * n.saturating_sub(1) / 2 {
        return Err(RDKitError::Value(
            "the distance matrix is too small for the pool".into(),
        ));
    }
    Ok(())
}

/// The argument checks shared by the pickers.
fn check_pool(pool_size: u32, pick_size: u32) -> Result<()> {
    if pool_size == 0 {
        return Err(RDKitError::Value("empty pool to pick from".into()));
    }
    if pool_size < pick_size {
        return Err(RDKitError::Value(
            "pickSize cannot be larger than the poolSize".into(),
        ));
    }
    Ok(())
}

/// A user-supplied pick index as a pool index.
fn pool_index(pick: i32, pool_size: u32) -> Result<u32> {
    u32::try_from(pick)
        .ok()
        .filter(|&p| p < pool_size)
        .ok_or_else(|| RDKitError::Value("pick index was larger than the poolSize".into()))
}
