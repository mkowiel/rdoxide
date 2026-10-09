//! Murtagh's hierarchical agglomerative clustering (upstream
//! `ML/Cluster/Murtagh/{hc.c,hcdriver.c,Clustering.cpp}`, a translation of
//! F. Murtagh's `HC` Fortran program).
//!
//! The nearest-neighbour-list algorithm works on the lower triangle of a
//! dissimilarity matrix stored row by row without the diagonal. Seven
//! agglomeration criteria are supported, see [`Linkage`]. The cluster ids in
//! the result are 1-based, as in the original program.

use crate::error::{RDKitError, Result};
use crate::rdgeneral::threads::par_map;

/// Agglomeration criterion (`iopt` of `HC`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Linkage {
    /// Ward's minimum variance method.
    Ward = 1,
    /// Single link.
    Slink = 2,
    /// Complete link.
    Clink = 3,
    /// Average link (group average, UPGMA).
    Upgma = 4,
    /// McQuitty's method.
    McQuitty = 5,
    /// Gower's median method.
    Gower = 6,
    /// Centroid method.
    Centroid = 7,
}

impl Linkage {
    /// The criterion for `HC`'s numeric option, if there is one.
    pub fn from_option(option: i64) -> Option<Self> {
        Some(match option {
            1 => Linkage::Ward,
            2 => Linkage::Slink,
            3 => Linkage::Clink,
            4 => Linkage::Upgma,
            5 => Linkage::McQuitty,
            6 => Linkage::Gower,
            7 => Linkage::Centroid,
            _ => return None,
        })
    }
}

/// The agglomeration history of a clustering of `n` items. Each of the first
/// `n - 1` steps joins clusters `ia[step]` and `ib[step]` (1-based item ids;
/// a merged cluster keeps the smaller id) at dissimilarity `crit[step]`. The
/// last entry of every vector is unused (zero).
#[derive(Debug, Clone, PartialEq)]
pub struct Dendrogram {
    pub ia: Vec<i64>,
    pub ib: Vec<i64>,
    pub crit: Vec<f64>,
}

/// `IOFFSET`: position in the packed matrix of the 1-based pair `(i, j)`.
fn ioffset(i: usize, j: usize) -> usize {
    let (lo, hi) = if j > i { (i, j) } else { (j, i) };
    (hi - 1) * (hi - 2) / 2 + lo
}

const INF: f64 = 1e20;

/// `HC`: clusters `n` items from the packed dissimilarity matrix `diss`,
/// which is overwritten (it holds `n * (n - 1) / 2` values).
///
/// Variable names follow the original program; the Fortran `SAVE`d
/// variables `im`, `jm` and `jj` keep their value between steps, which the
/// algorithm relies on in degenerate cases (e.g. ties at infinity).
pub fn hierarchical_cluster(n: usize, linkage: Linkage, diss: &mut [f64]) -> Result<Dendrogram> {
    if diss.len() < n * n.saturating_sub(1) / 2 {
        return Err(RDKitError::Value(
            "the dissimilarity matrix is too small".into(),
        ));
    }
    let mut res = Dendrogram {
        ia: vec![0; n],
        ib: vec![0; n],
        crit: vec![0.0; n],
    };
    if n < 2 {
        return Ok(res);
    }
    // 1-based work arrays, like the original
    let mut membr = vec![1.0f64; n + 1];
    let mut flag = vec![true; n + 1];
    let mut nn = vec![0usize; n + 1];
    let mut disnn = vec![0.0f64; n + 1];
    // 1-based view of the matrix
    let d = |ind: usize| ind - 1;
    let (mut im, mut jm, mut jj) = (0usize, 0usize, 0usize);

    if linkage == Linkage::Ward {
        // (merging criteria are variances rather than distances)
        diss[..n * (n - 1) / 2].iter_mut().for_each(|x| *x /= 2.0);
    }

    // first create the list of nearest neighbours
    for i in 1..n {
        let mut dmin = INF;
        for j in i + 1..=n {
            let v = diss[d(ioffset(i, j))];
            if v < dmin {
                dmin = v;
                jm = j;
            }
        }
        nn[i] = jm;
        disnn[i] = dmin;
    }

    let mut ncl = n;
    while ncl > 1 {
        // the least dissimilarity according to the list of neighbours
        let mut dmin = INF;
        for i in 1..n {
            if flag[i] && disnn[i] < dmin {
                dmin = disnn[i];
                im = i;
                jm = nn[i];
            }
        }
        ncl -= 1;

        // agglomerate
        let (i2, j2) = (im.min(jm), im.max(jm));
        res.ia[n - ncl - 1] = i2 as i64;
        res.ib[n - ncl - 1] = j2 as i64;
        res.crit[n - ncl - 1] = dmin;

        // update the dissimilarities from the new cluster
        flag[j2] = false;
        dmin = INF;
        for k in 1..=n {
            if !flag[k] || k == i2 {
                continue;
            }
            let x = membr[i2] + membr[j2] + membr[k];
            let ind1 = ioffset(i2, k);
            let ind2 = ioffset(j2, k);
            let xx = diss[d(ioffset(i2, j2))];
            let (d1, d2) = (diss[d(ind1)], diss[d(ind2)]);
            let (mi, mj, mk) = (membr[i2], membr[j2], membr[k]);
            diss[d(ind1)] = match linkage {
                Linkage::Ward => ((mi + mk) * d1 + (mj + mk) * d2 - mk * xx) / x,
                Linkage::Slink => d1.min(d2),
                Linkage::Clink => d1.max(d2),
                Linkage::Upgma => (mi * d1 + mj * d2) / (mi + mj),
                Linkage::McQuitty => d1 * 0.5 + d2 * 0.5,
                Linkage::Gower => d1 * 0.5 + d2 * 0.5 - xx * 0.25,
                Linkage::Centroid => (mi * d1 + mj * d2 - mi * mj * xx / (mi + mj)) / (mi + mj),
            };
            if i2 < k && diss[d(ind1)] < dmin {
                dmin = diss[d(ind1)];
                jj = k;
            }
        }
        membr[i2] += membr[j2];
        disnn[i2] = dmin;
        nn[i2] = jj;

        // update the list of neighbours where it is needed
        for i in 1..n {
            if flag[i] && (nn[i] == i2 || nn[i] == j2) {
                let mut dmin = INF;
                for j in i + 1..=n {
                    let v = diss[d(ioffset(i, j))];
                    if flag[j] && v < dmin {
                        dmin = v;
                        jj = j;
                    }
                }
                nn[i] = jj;
                disnn[i] = dmin;
            }
        }
    }
    Ok(res)
}

/// `MurtaghDistCluster`: clusters `n` items given their packed distances.
pub fn cluster_distances(n: usize, linkage: Linkage, dists: &[f64]) -> Result<Dendrogram> {
    let len = n * n.saturating_sub(1) / 2;
    let mut diss = dists
        .get(..len)
        .ok_or_else(|| RDKitError::Value("the dissimilarity matrix is too small".into()))?
        .to_vec();
    hierarchical_cluster(n, linkage, &mut diss)
}

/// Squared Euclidean distances between the rows of `points` (`n` rows of
/// `dim` values), packed as the lower triangle row by row. Rows are computed
/// on up to `num_threads` threads; the values do not depend on the count.
pub fn squared_distances(points: &[f64], n: usize, dim: usize, num_threads: i32) -> Vec<f64> {
    let rows: Vec<usize> = (1..n).collect();
    let per_row = par_map(&rows, num_threads, |&i| {
        let a = &points[i * dim..(i + 1) * dim];
        (0..i)
            .map(|j| {
                let b = &points[j * dim..(j + 1) * dim];
                a.iter()
                    .zip(b)
                    .fold(0.0, |acc, (x, y)| acc + (x - y) * (x - y))
            })
            .collect::<Vec<f64>>()
    });
    per_row.into_iter().flatten().collect()
}

/// `MurtaghCluster`: clusters `n` points of dimension `dim` (row-major
/// `points`) by their squared Euclidean distances.
pub fn cluster_points(
    points: &[f64],
    n: usize,
    dim: usize,
    linkage: Linkage,
    num_threads: i32,
) -> Result<Dendrogram> {
    if points.len() < n * dim {
        return Err(RDKitError::Value("the data array is too small".into()));
    }
    let mut diss = squared_distances(points, n, dim, num_threads);
    hierarchical_cluster(n, linkage, &mut diss)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Five points on a line: 0, 1, 3, 7, 8.
    fn line_distances() -> Vec<f64> {
        let x = [0.0, 1.0, 3.0, 7.0, 8.0f64];
        (1..5)
            .flat_map(|i| (0..i).map(move |j| (x[i] - x[j]).abs()))
            .collect()
    }

    #[test]
    fn single_link_joins_nearest_first() {
        let d = cluster_distances(5, Linkage::Slink, &line_distances()).unwrap();
        assert_eq!(d.ia, [1, 4, 1, 1, 0]);
        assert_eq!(d.ib, [2, 5, 3, 4, 0]);
        assert_eq!(d.crit, [1.0, 1.0, 2.0, 4.0, 0.0]);
    }

    #[test]
    fn complete_link_criteria_are_the_largest_distances() {
        let d = cluster_distances(5, Linkage::Clink, &line_distances()).unwrap();
        assert_eq!(d.crit, [1.0, 1.0, 3.0, 8.0, 0.0]);
    }

    #[test]
    fn short_inputs_are_rejected_or_trivial() {
        assert!(cluster_distances(5, Linkage::Ward, &[1.0]).is_err());
        assert_eq!(cluster_distances(1, Linkage::Ward, &[]).unwrap().ia, [0]);
    }

    #[test]
    fn distances_do_not_depend_on_the_thread_count() {
        let points: Vec<f64> = (0..300).map(|i| f64::from(i * 37 % 101)).collect();
        let serial = squared_distances(&points, 100, 3, 1);
        assert_eq!(serial.len(), 100 * 99 / 2);
        assert_eq!(serial, squared_distances(&points, 100, 3, 4));
    }
}
