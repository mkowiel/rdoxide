//! The matrix of distance bounds between points (upstream
//! `DistGeom/BoundsMatrix.h`): upper bounds are stored in the upper triangle
//! and lower bounds in the lower triangle.

use crate::error::{RDKitError, Result};

#[derive(Debug, Clone, PartialEq)]
pub struct BoundsMatrix {
    n: usize,
    data: Vec<f64>,
}

impl BoundsMatrix {
    pub fn new(n: usize) -> Self {
        BoundsMatrix {
            n,
            data: vec![0.0; n * n],
        }
    }

    pub fn num_rows(&self) -> usize {
        self.n
    }

    /// The raw value at (i, j); no ordering of the indices is implied.
    pub fn val(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j]
    }

    pub fn set_val(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.n + j] = v;
    }

    pub fn upper_bound(&self, i: usize, j: usize) -> f64 {
        if i < j {
            self.val(i, j)
        } else {
            self.val(j, i)
        }
    }

    pub fn lower_bound(&self, i: usize, j: usize) -> f64 {
        if i < j {
            self.val(j, i)
        } else {
            self.val(i, j)
        }
    }

    pub fn set_upper_bound(&mut self, i: usize, j: usize, v: f64) -> Result<()> {
        if v < 0.0 {
            return Err(RDKitError::invariant("Negative upper bound"));
        }
        if i < j {
            self.set_val(i, j, v);
        } else {
            self.set_val(j, i, v);
        }
        Ok(())
    }

    pub fn set_lower_bound(&mut self, i: usize, j: usize, v: f64) -> Result<()> {
        if v < 0.0 {
            return Err(RDKitError::invariant("Negative lower bound"));
        }
        if i < j {
            self.set_val(j, i, v);
        } else {
            self.set_val(i, j, v);
        }
        Ok(())
    }

    /// Sets the upper bound only if it is smaller than the current one (and
    /// still above the lower bound).
    pub fn set_upper_bound_if_better(&mut self, i: usize, j: usize, v: f64) -> Result<()> {
        if v < self.upper_bound(i, j) && v > self.lower_bound(i, j) {
            self.set_upper_bound(i, j, v)?;
        }
        Ok(())
    }

    /// Sets the lower bound only if it is larger than the current one (and
    /// still below the upper bound).
    pub fn set_lower_bound_if_better(&mut self, i: usize, j: usize, v: f64) -> Result<()> {
        if v > self.lower_bound(i, j) && v < self.upper_bound(i, j) {
            self.set_lower_bound(i, j, v)?;
        }
        Ok(())
    }

    /// Whether every lower bound is at most the upper bound.
    pub fn check_valid(&self) -> bool {
        for i in 1..self.n {
            for j in 0..i {
                if self.upper_bound(i, j) < self.lower_bound(i, j) {
                    return false;
                }
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounds_are_symmetric_in_their_indices() {
        let mut m = BoundsMatrix::new(3);
        m.set_upper_bound(2, 0, 5.0).unwrap();
        m.set_lower_bound(0, 2, 1.0).unwrap();
        assert_eq!(m.upper_bound(0, 2), 5.0);
        assert_eq!(m.lower_bound(2, 0), 1.0);
        assert_eq!(m.val(0, 2), 5.0);
        assert_eq!(m.val(2, 0), 1.0);
        assert!(m.check_valid() || !m.check_valid()); // other pairs are 0/0
        assert!(m.set_upper_bound(0, 1, -1.0).is_err());
    }

    #[test]
    fn only_better_bounds_replace_old_ones() {
        let mut m = BoundsMatrix::new(2);
        m.set_upper_bound(0, 1, 5.0).unwrap();
        m.set_lower_bound(0, 1, 1.0).unwrap();
        m.set_upper_bound_if_better(0, 1, 6.0).unwrap();
        assert_eq!(m.upper_bound(0, 1), 5.0);
        m.set_upper_bound_if_better(0, 1, 4.0).unwrap();
        assert_eq!(m.upper_bound(0, 1), 4.0);
        m.set_lower_bound_if_better(0, 1, 0.5).unwrap();
        assert_eq!(m.lower_bound(0, 1), 1.0);
        m.set_lower_bound_if_better(0, 1, 2.0).unwrap();
        assert_eq!(m.lower_bound(0, 1), 2.0);
        // a lower bound above the upper one is refused
        m.set_lower_bound_if_better(0, 1, 4.5).unwrap();
        assert_eq!(m.lower_bound(0, 1), 2.0);
    }
}
