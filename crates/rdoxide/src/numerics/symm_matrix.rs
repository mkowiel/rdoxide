//! A symmetric matrix stored as its lower triangle (upstream
//! `Numerics/SymmMatrix.h`).

#[derive(Debug, Clone, PartialEq)]
pub struct SymmMatrix {
    n: usize,
    data: Vec<f64>,
}

impl SymmMatrix {
    pub fn new(n: usize) -> Self {
        Self::filled(n, 0.0)
    }

    pub fn filled(n: usize, value: f64) -> Self {
        SymmMatrix {
            n,
            data: vec![value; n * (n + 1) / 2],
        }
    }

    pub fn num_rows(&self) -> usize {
        self.n
    }

    fn index(i: usize, j: usize) -> usize {
        if i >= j {
            i * (i + 1) / 2 + j
        } else {
            j * (j + 1) / 2 + i
        }
    }

    pub fn val(&self, i: usize, j: usize) -> f64 {
        self.data[Self::index(i, j)]
    }

    pub fn set_val(&mut self, i: usize, j: usize, v: f64) {
        self.data[Self::index(i, j)] = v;
    }

    /// The lower triangle, row by row.
    pub fn data(&self) -> &[f64] {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut [f64] {
        &mut self.data
    }

    /// `y = A x`
    pub fn multiply(&self, x: &[f64], y: &mut [f64]) {
        assert!(x.len() == self.n && y.len() == self.n, "size mismatch");
        for i in 0..self.n {
            let mut acc = 0.0;
            let row = i * (i + 1) / 2;
            for j in 0..=i {
                acc += self.data[row + j] * x[j];
            }
            let mut id = row + i;
            for j in i + 1..self.n {
                id += j;
                acc += self.data[id] * x[j];
            }
            y[i] = acc;
        }
    }
}
