//! Diversity picking by hierarchical clustering (upstream
//! `SimDivPickers/HierarchicalClusterPicker.{h,cpp}`): the pool is clustered
//! (Murtagh's algorithm) until `pick_size` clusters are left and the most
//! central item of each cluster is picked.

use super::{check_dist_mat, dist_from_ltm};
use crate::error::{RDKitError, Result};
use crate::ml::cluster::murtagh::{Linkage, hierarchical_cluster};

/// `HierarchicalClusterPicker::ClusterMethod`; the values are the
/// agglomeration criteria of [`Linkage`].
pub type ClusterMethod = Linkage;

/// Clusters the pool with `dist_mat` and returns the clusters together with
/// the matrix as the clustering left it: the algorithm overwrites the
/// dissimilarities it merges.
fn cluster_in_place(
    dist_mat: &[f64],
    pool_size: u32,
    pick_size: u32,
    method: ClusterMethod,
) -> Result<(Vec<Vec<i32>>, Vec<f64>)> {
    if pool_size < pick_size {
        return Err(RDKitError::precondition(
            "pickSize cannot be larger than the poolSize",
        ));
    }
    if pick_size == 0 {
        return Err(RDKitError::Value("pickSize must be at least 1".into()));
    }
    check_dist_mat(dist_mat, pool_size)?;
    let n = pool_size as usize;
    let mut diss = dist_mat[..n * (n - 1) / 2].to_vec();
    let history = hierarchical_cluster(n, method, &mut diss)?;

    // every item starts as a cluster; a step moves the items of cluster
    // `ib` into cluster `ia`
    let mut clusters: Vec<Vec<i32>> = (0..pool_size as i32).map(|i| vec![i]).collect();
    let mut removed = vec![false; n];
    for step in 0..n - pick_size as usize {
        let (kept, gone) = (history.ia[step] as usize - 1, history.ib[step] as usize - 1);
        let moved = std::mem::take(&mut clusters[gone]);
        clusters[kept].extend(moved);
        removed[gone] = true;
    }
    let clusters = clusters
        .into_iter()
        .zip(removed)
        .filter_map(|(c, gone)| (!gone).then_some(c))
        .collect();
    Ok((clusters, diss))
}

/// `HierarchicalClusterPicker::cluster`: merges the `pool_size` items until
/// `pick_size` clusters are left and returns them (the items of each cluster
/// in the order they were merged).
pub fn cluster(
    dist_mat: &[f64],
    pool_size: u32,
    pick_size: u32,
    method: ClusterMethod,
) -> Result<Vec<Vec<i32>>> {
    Ok(cluster_in_place(dist_mat, pool_size, pick_size, method)?.0)
}

/// `HierarchicalClusterPicker::pick`: the item of each cluster with the
/// smallest sum of squared distances to the other items of the cluster.
///
/// Like upstream, the distances are those of the matrix *after* the
/// clustering overwrote it (upstream clusters in place and then reads the
/// same buffer), so the representatives are not always the true medoids.
pub fn pick(
    dist_mat: &[f64],
    pool_size: u32,
    pick_size: u32,
    method: ClusterMethod,
) -> Result<Vec<i32>> {
    let (clusters, clobbered) = cluster_in_place(dist_mat, pool_size, pick_size, method)?;
    if clusters.len() != pick_size as usize {
        return Err(RDKitError::invariant("wrong number of clusters"));
    }
    Ok(clusters
        .iter()
        .map(|members| {
            let sum_d2 = |candidate: usize| -> f64 {
                members
                    .iter()
                    .enumerate()
                    .filter(|&(pos, _)| pos != candidate)
                    .map(|(_, &other)| {
                        let d = dist_from_ltm(&clobbered, members[candidate] as u32, other as u32);
                        d * d
                    })
                    .sum()
            };
            let mut best = (members[0], f64::MAX);
            for (pos, &candidate) in members.iter().enumerate() {
                let d2 = sum_d2(pos);
                if d2 < best.1 {
                    best = (candidate, d2);
                }
            }
            best.0
        })
        .collect())
}
