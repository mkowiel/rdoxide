//! `freesasa/src/nb.c`: neighbor lists of atoms whose (probe-enlarged) spheres overlap, built with cell lists. The order of
//! the neighbors follows the cell traversal exactly because the Lee-Richards sums depend on it in the last bits.

/// `nb_list`: for every atom its neighbors with the x, y and xy distances to them
#[derive(Debug, Clone, Default)]
pub struct NbList {
    pub nb: Vec<Vec<usize>>,
    pub xyd: Vec<Vec<f64>>,
    pub xd: Vec<Vec<f64>>,
    pub yd: Vec<Vec<f64>>,
}

struct Cells {
    nx: usize,
    ny: usize,
    nz: usize,
    d: f64,
    min: [f64; 3],
    atoms: Vec<Vec<usize>>,
}

impl Cells {
    /// `cell_list_bounds` and `fill_cells`
    fn new(cell_size: f64, xyz: &[f64]) -> Cells {
        let n = xyz.len() / 3;
        let d = cell_size;
        let (mut lo, mut hi) = ([xyz[0], xyz[1], xyz[2]], [xyz[0], xyz[1], xyz[2]]);
        for i in 1..n {
            for k in 0..3 {
                lo[k] = lo[k].min(xyz[3 * i + k]);
                hi[k] = hi[k].max(xyz[3 * i + k]);
            }
        }
        let min = [lo[0] - d / 2., lo[1] - d / 2., lo[2] - d / 2.];
        let max = [hi[0] + d / 2., hi[1] + d / 2., hi[2] + d / 2.];
        let cnt = |k: usize| ((max[k] - min[k]) / d).ceil() as usize;
        let (nx, ny, nz) = (cnt(0), cnt(1), cnt(2));
        let mut c = Cells {
            nx,
            ny,
            nz,
            d,
            min,
            atoms: vec![Vec::new(); nx * ny * nz],
        };
        for i in 0..n {
            let ix = ((xyz[3 * i] - min[0]) / d) as usize;
            let iy = ((xyz[3 * i + 1] - min[1]) / d) as usize;
            let iz = ((xyz[3 * i + 2] - min[2]) / d) as usize;
            let idx = c.index(ix, iy, iz);
            c.atoms[idx].push(i);
        }
        c
    }

    fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix + self.nx * (iy + self.ny * iz)
    }

    /// `fill_nb`: the cell itself and its forward neighbors, in upstream's loop order
    fn forward_neighbors(&self, ix: usize, iy: usize, iz: usize) -> Vec<usize> {
        let xmin = ix.saturating_sub(1);
        let xmax = if ix < self.nx - 1 { ix + 1 } else { ix };
        let ymin = iy.saturating_sub(1);
        let ymax = if iy < self.ny - 1 { iy + 1 } else { iy };
        let zmin = iz.saturating_sub(1);
        let zmax = if iz < self.nz - 1 { iz + 1 } else { iz };
        let mut out = Vec::new();
        for i in xmin..=xmax {
            for j in ymin..=ymax {
                for k in zmin..=zmax {
                    // the scalar product with (1,1,1) is non-negative: no double counting of cell pairs
                    if i as i64 - ix as i64 + j as i64 - iy as i64 + k as i64 - iz as i64 >= 0 {
                        out.push(self.index(i, j, k));
                    }
                }
            }
        }
        out
    }
}

impl NbList {
    /// `freesasa_nb_new`
    pub fn new(xyz: &[f64], radii: &[f64]) -> NbList {
        let n = xyz.len() / 3;
        let mut nb = NbList {
            nb: vec![Vec::new(); n],
            xyd: vec![Vec::new(); n],
            xd: vec![Vec::new(); n],
            yd: vec![Vec::new(); n],
        };
        let cell_size = 2.0 * radii.iter().fold(0.0f64, |m, &r| m.max(r));
        let cells = Cells::new(cell_size, xyz);
        let _ = cells.d;
        for ic in 0..cells.atoms.len() {
            let iz = ic / (cells.nx * cells.ny);
            let iy = (ic / cells.nx) % cells.ny;
            let ix = ic % cells.nx;
            for jc in cells.forward_neighbors(ix, iy, iz) {
                nb.cell_pair(xyz, radii, &cells.atoms[ic], &cells.atoms[jc], ic == jc);
            }
        }
        let _ = cells.min;
        nb
    }

    /// `nb_calc_cell_pair`
    fn cell_pair(&mut self, v: &[f64], radii: &[f64], ci: &[usize], cj: &[usize], same: bool) {
        for (i, &ia) in ci.iter().enumerate() {
            let ri = radii[ia];
            let (xi, yi, zi) = (v[ia * 3], v[ia * 3 + 1], v[ia * 3 + 2]);
            let start = if same { i + 1 } else { 0 };
            for &ja in &cj[start..] {
                let rj = radii[ja];
                let cut2 = (ri + rj) * (ri + rj);
                let (dx, dy, dz) = (v[ja * 3] - xi, v[ja * 3 + 1] - yi, v[ja * 3 + 2] - zi);
                if dx * dx + dy * dy + dz * dz < cut2 {
                    self.add_pair(ia, ja, dx, dy);
                }
            }
        }
    }

    /// `nb_add_pair`
    fn add_pair(&mut self, i: usize, j: usize, dx: f64, dy: f64) {
        let d = (dx * dx + dy * dy).sqrt();
        self.nb[i].push(j);
        self.nb[j].push(i);
        self.xyd[i].push(d);
        self.xyd[j].push(d);
        self.xd[i].push(dx);
        self.xd[j].push(-dx);
        self.yd[i].push(dy);
        self.yd[j].push(-dy);
    }
}
