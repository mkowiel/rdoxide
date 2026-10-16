//! The shortest augmenting path algorithm for the rectangular linear
//! assignment problem, after Crouse, "On implementing 2D rectangular
//! assignment algorithms", IEEE Trans. Aerosp. Electron. Syst. 52(4):1679
//! (2016), as implemented by PM Larsen in SciPy (`lap_a_la_scipy.cpp`).
//!
//! Copyright (c) 2001-2002 Enthought, Inc. 2003-2023, SciPy Developers.
//! All rights reserved. Licensed under the BSD 3-Clause license (see
//! LICENSES/SciPy-BSD-3-Clause.txt).

const UNSET: usize = usize::MAX;

/// The augmenting path for row `i`: the sink column and the minimum value, or
/// `None` for an infeasible cost matrix.
#[allow(clippy::too_many_arguments)]
fn augmenting_path(
    nc: usize,
    cost: &[i32],
    u: &[f64],
    v: &[f64],
    path: &mut [usize],
    row4col: &[usize],
    shortest_path_costs: &mut [f64],
    mut i: usize,
    sr: &mut [bool],
    sc: &mut [bool],
    remaining: &mut [usize],
) -> Option<(usize, f64)> {
    let mut min_val = 0.0;

    // Crouse's pseudocode uses set complements to keep track of remaining
    // nodes. Here we use a vector, as it is more efficient.
    let mut num_remaining = nc;
    for (it, r) in remaining.iter_mut().enumerate().take(nc) {
        // Filling this up in reverse order ensures that the solution of a
        // constant cost matrix is the identity matrix.
        *r = nc - it - 1;
    }
    sr.fill(false);
    sc.fill(false);
    shortest_path_costs.fill(f64::MAX);

    // find shortest augmenting path
    let mut sink = None;
    while sink.is_none() {
        let mut index = UNSET;
        let mut lowest = f64::MAX;
        sr[i] = true;

        for it in 0..num_remaining {
            let j = remaining[it];
            let r = min_val + f64::from(cost[i * nc + j]) - u[i] - v[j];
            if r < shortest_path_costs[j] {
                path[j] = i;
                shortest_path_costs[j] = r;
            }
            // When multiple nodes have the minimum cost, we select one which
            // gives us a new sink node. This is particularly important for
            // integer cost matrices with small co-efficients.
            if shortest_path_costs[j] < lowest
                || (shortest_path_costs[j] == lowest && row4col[j] == UNSET)
            {
                lowest = shortest_path_costs[j];
                index = it;
            }
        }

        min_val = lowest;
        if min_val == f64::MAX {
            // infeasible cost matrix
            return None;
        }
        let j = remaining[index];
        if row4col[j] == UNSET {
            sink = Some(j);
        } else {
            i = row4col[j];
        }
        sc[j] = true;
        num_remaining -= 1;
        remaining[index] = remaining[num_remaining];
    }
    sink.map(|s| (s, min_val))
}

/// The row/column assignment maximising the summed cost: `Some((a, b))` with
/// `a[k]`, `b[k]` the matched row and column of the `k`th assignment, there
/// being `min(rows, columns)` of them. `None` if no solution is possible.
pub fn lap_maximize(costs: &[Vec<i32>]) -> Option<(Vec<usize>, Vec<usize>)> {
    if costs.is_empty() || costs[0].is_empty() {
        return Some((Vec::new(), Vec::new()));
    }
    let mut nr = costs.len();
    let mut nc = costs[0].len();
    let transpose = nc < nr;
    let mut cost = vec![0; nc * nr];
    // for maximization, take -ve of costs.
    for (i, row) in costs.iter().enumerate() {
        for (j, &c) in row.iter().enumerate() {
            if transpose {
                cost[j * nr + i] = -c;
            } else {
                cost[i * nc + j] = -c;
            }
        }
    }
    if transpose {
        std::mem::swap(&mut nc, &mut nr);
    }
    let mut u = vec![0.0; nr];
    let mut v = vec![0.0; nc];
    let mut shortest_path_costs = vec![0.0; nc];
    let mut path = vec![UNSET; nc];
    let mut col4row = vec![UNSET; nr];
    let mut row4col = vec![UNSET; nc];
    let mut sr = vec![false; nr];
    let mut sc = vec![false; nc];
    let mut remaining = vec![0; nc];

    // iteratively build the solution
    for cur_row in 0..nr {
        let (sink, min_val) = augmenting_path(
            nc,
            &cost,
            &u,
            &v,
            &mut path,
            &row4col,
            &mut shortest_path_costs,
            cur_row,
            &mut sr,
            &mut sc,
            &mut remaining,
        )?;

        // update dual variables
        u[cur_row] += min_val;
        for i in 0..nr {
            if sr[i] && i != cur_row {
                u[i] += min_val - shortest_path_costs[col4row[i]];
            }
        }
        for j in 0..nc {
            if sc[j] {
                v[j] -= min_val - shortest_path_costs[j];
            }
        }

        // augment previous solution
        let mut j = sink;
        loop {
            let i = path[j];
            row4col[j] = i;
            std::mem::swap(&mut col4row[i], &mut j);
            if i == cur_row {
                break;
            }
        }
    }

    if transpose {
        let mut order: Vec<usize> = (0..col4row.len()).collect();
        order.sort_by_key(|&k| col4row[k]);
        // (the sort keys are distinct, so the order is unambiguous)
        let a = order.iter().map(|&k| col4row[k]).collect();
        Some((a, order))
    } else {
        Some(((0..nr).collect(), col4row))
    }
}
