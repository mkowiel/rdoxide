//! Linear assignment (Jonker–Volgenant) and the distance histograms that feed
//! the O3A cost matrix (`LAP` and `MolHistogram` of upstream
//! `O3AAlignMolecules.cpp`).

/// Cost of pairing an atom with a dummy.
pub const O3_DUMMY_COST: i32 = 100_000;
/// Cost of a constrained pair, so that the assignment always picks it.
pub const O3_LARGE_NEGATIVE_WEIGHT: i32 = -1_000_000_000;
/// Number of distance bins of a [`MolHistogram`].
pub const O3_MAX_H_BINS: usize = 20;

/// Per heavy atom histogram of the topological distances to all atoms.
#[derive(Debug, Clone)]
pub struct MolHistogram {
    rows: Vec<[i32; O3_MAX_H_BINS]>,
}

impl MolHistogram {
    /// `atomic_nums` are those of all atoms, `dmat` is the row-major
    /// `n x n` distance matrix.
    pub fn new(atomic_nums: &[u32], dmat: &[f64]) -> Self {
        let n = atomic_nums.len();
        let mut rows = Vec::new();
        for i in 0..n {
            if atomic_nums[i] == 1 {
                continue;
            }
            let mut row = [0i32; O3_MAX_H_BINS];
            for j in 0..n {
                let dist = dmat[i * n + j].floor() as u32 as usize;
                if dist < O3_MAX_H_BINS {
                    row[dist] += 1;
                }
            }
            rows.push(row);
        }
        Self { rows }
    }

    pub fn get(&self, y: usize, x: usize) -> i32 {
        self.rows[y][x]
    }
}

/// Square linear assignment problem solved with the Jonker–Volgenant
/// shortest augmenting path algorithm (Computing 1987, 38, 325–340).
#[derive(Debug, Clone)]
pub struct Lap {
    dim: usize,
    pub row_sol: Vec<i32>,
    col_sol: Vec<i32>,
    free: Vec<i32>,
    col_list: Vec<i32>,
    matches: Vec<i32>,
    d: Vec<i32>,
    v: Vec<i32>,
    pred: Vec<i32>,
    cost: Vec<i32>,
}

impl Lap {
    pub fn new(dim: usize) -> Self {
        Self {
            dim,
            row_sol: vec![0; dim],
            col_sol: vec![0; dim],
            free: vec![0; dim],
            col_list: vec![0; dim],
            matches: vec![0; dim],
            d: vec![0; dim],
            v: vec![0; dim],
            pred: vec![0; dim],
            cost: vec![0; dim * dim],
        }
    }

    pub fn cost(&self, i: usize, j: usize) -> i32 {
        self.cost[i * self.dim + j]
    }

    pub fn set_cost(&mut self, i: usize, j: usize, c: i32) {
        self.cost[i * self.dim + j] = c;
    }

    pub fn row_sol(&self, i: usize) -> usize {
        self.row_sol[i] as usize
    }

    /// `LAP::computeMinCostPath`. The scan loops start at the value `up` has on
    /// entry while `up` itself keeps growing, exactly as in the C original.
    #[allow(clippy::mut_range_bound)]
    pub fn compute_min_cost_path(&mut self, dim: usize) {
        let n = self.dim;
        let c = |cost: &Vec<i32>, i: usize, j: usize| cost[i * n + j];
        let mut num_free = 0usize;
        let mut j2 = 0usize;
        let mut last: isize = 0;
        let mut end_of_path = 0usize;

        for m in self.matches.iter_mut().take(dim) {
            *m = 0;
        }
        // column reduction, in reverse order
        for j in (0..dim).rev() {
            let mut min = c(&self.cost, 0, j);
            let mut imin = 0;
            for i in 1..dim {
                if c(&self.cost, i, j) < min {
                    min = c(&self.cost, i, j);
                    imin = i;
                }
            }
            self.v[j] = min;
            self.matches[imin] += 1;
            if self.matches[imin] == 1 {
                self.row_sol[imin] = j as i32;
                self.col_sol[j] = imin as i32;
            } else {
                self.col_sol[j] = -1;
            }
        }
        // reduction transfer
        for i in 0..dim {
            if self.matches[i] == 0 {
                self.free[num_free] = i as i32;
                num_free += 1;
            } else if self.matches[i] == 1 {
                let j1 = self.row_sol[i] as usize;
                let mut min = O3_DUMMY_COST;
                for j in 0..dim {
                    if j != j1 && c(&self.cost, i, j) - self.v[j] < min {
                        min = c(&self.cost, i, j) - self.v[j];
                    }
                }
                self.v[j1] -= min;
            }
        }
        // augmenting row reduction, twice
        for _ in 0..2 {
            let mut k = 0usize;
            let prv_num_free = num_free;
            num_free = 0;
            while k < prv_num_free {
                let i = self.free[k] as usize;
                k += 1;
                let mut u_min = c(&self.cost, i, 0) - self.v[0];
                let mut j1 = 0usize;
                let mut u_sub_min = O3_DUMMY_COST;
                for j in 1..dim {
                    let h = c(&self.cost, i, j) - self.v[j];
                    if h < u_sub_min {
                        if h >= u_min {
                            u_sub_min = h;
                            j2 = j;
                        } else {
                            u_sub_min = u_min;
                            u_min = h;
                            j2 = j1;
                            j1 = j;
                        }
                    }
                }
                let mut i0 = self.col_sol[j1];
                if u_min < u_sub_min {
                    self.v[j1] -= u_sub_min - u_min;
                } else if i0 >= 0 {
                    j1 = j2;
                    i0 = self.col_sol[j2];
                }
                self.row_sol[i] = j1 as i32;
                self.col_sol[j1] = i as i32;
                if i0 >= 0 {
                    if u_min < u_sub_min {
                        k -= 1;
                        self.free[k] = i0;
                    } else {
                        self.free[num_free] = i0;
                        num_free += 1;
                    }
                }
            }
        }
        // augment the solution for each free row
        for f in 0..num_free {
            let free_row = self.free[f] as usize;
            for j in 0..dim {
                self.d[j] = c(&self.cost, free_row, j) - self.v[j];
                self.pred[j] = free_row as i32;
                self.col_list[j] = j as i32;
            }
            let mut low = 0usize;
            let mut up = 0usize;
            let mut unassigned_found = false;
            let mut min = 0i32;
            while !unassigned_found {
                if up == low {
                    last = low as isize - 1;
                    min = self.d[self.col_list[up] as usize];
                    up += 1;
                    for k in up..dim {
                        let j = self.col_list[k] as usize;
                        let h = self.d[j];
                        if h <= min {
                            if h < min {
                                up = low;
                                min = h;
                            }
                            self.col_list[k] = self.col_list[up];
                            self.col_list[up] = j as i32;
                            up += 1;
                        }
                    }
                    for k in low..up {
                        if self.col_sol[self.col_list[k] as usize] < 0 {
                            end_of_path = self.col_list[k] as usize;
                            unassigned_found = true;
                            break;
                        }
                    }
                }
                if !unassigned_found {
                    let j1 = self.col_list[low] as usize;
                    low += 1;
                    let i = self.col_sol[j1] as usize;
                    let h = c(&self.cost, i, j1) - self.v[j1] - min;
                    for k in up..dim {
                        let j = self.col_list[k] as usize;
                        let v2 = c(&self.cost, i, j) - self.v[j] - h;
                        if v2 < self.d[j] {
                            self.pred[j] = i as i32;
                            if v2 == min {
                                if self.col_sol[j] < 0 {
                                    end_of_path = j;
                                    unassigned_found = true;
                                    break;
                                }
                                self.col_list[k] = self.col_list[up];
                                self.col_list[up] = j as i32;
                                up += 1;
                            }
                            self.d[j] = v2;
                        }
                    }
                }
            }
            // update column prices
            let mut k = 0isize;
            while k <= last {
                let j1 = self.col_list[k as usize] as usize;
                self.v[j1] += self.d[j1] - min;
                k += 1;
            }
            // reset the assignments along the alternating path
            loop {
                let i = self.pred[end_of_path] as usize;
                self.col_sol[end_of_path] = i as i32;
                let j1 = end_of_path;
                end_of_path = self.row_sol[i] as usize;
                self.row_sol[i] = j1 as i32;
                if i == free_row {
                    break;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solves_a_small_assignment() {
        let costs = [[4, 1, 3], [2, 0, 5], [3, 2, 2]];
        let mut lap = Lap::new(3);
        for (i, row) in costs.iter().enumerate() {
            for (j, &c) in row.iter().enumerate() {
                lap.set_cost(i, j, c);
            }
        }
        lap.compute_min_cost_path(3);
        let total: i32 = (0..3).map(|i| lap.cost(i, lap.row_sol(i))).sum();
        assert_eq!(total, 5); // 1 + 2 + 2
    }
}
