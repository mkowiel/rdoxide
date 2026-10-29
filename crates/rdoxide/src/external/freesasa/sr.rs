//! `freesasa/src/sasa_sr.c`: Shrake & Rupley with test points on a golden-section spiral.

use super::nb::NbList;
use std::f64::consts::PI;

/// `test_points`: `n` points on the unit sphere
fn test_points(n: usize) -> Vec<f64> {
    let dlong = PI * (3.0 - 5.0f64.sqrt());
    let dz = 2.0 / n as f64;
    let (mut longitude, mut z) = (0.0f64, 1.0 - dz / 2.0);
    let mut tp = Vec::with_capacity(3 * n);
    for _ in 0..n {
        let r = (1.0 - z * z).sqrt();
        tp.push(longitude.cos() * r);
        tp.push(longitude.sin() * r);
        tp.push(z);
        z -= dz;
        longitude += dlong;
    }
    tp
}

/// `freesasa_shrake_rupley`: the area of every atom (single threaded: the sums per atom do not depend on the thread count)
pub fn shrake_rupley(xyz: &[f64], radii: &[f64], probe_radius: f64, n_points: usize) -> Vec<f64> {
    let n_atoms = xyz.len() / 3;
    let srp = test_points(n_points);
    let r: Vec<f64> = radii.iter().map(|&r| r + probe_radius).collect();
    let r2: Vec<f64> = r.iter().map(|&r| r * r).collect();
    let nb = NbList::new(xyz, &r);
    let mut sasa = vec![0.0; n_atoms];
    for i in 0..n_atoms {
        let ri = r[i];
        let vi = [xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]];
        let nbi = &nb.nb[i];
        let mut n_surface = 0usize;
        let mut current_nb = 0usize;
        for j in 0..n_points {
            let p = [
                srp[3 * j] * ri + vi[0],
                srp[3 * j + 1] * ri + vi[1],
                srp[3 * j + 2] * ri + vi[2],
            ];
            // an isolated atom is entirely exposed; upstream reads an uninitialized neighbor there
            if nbi.is_empty() {
                n_surface += 1;
                continue;
            }
            let buried_by = |a: usize| {
                let (dx, dy, dz) = (
                    p[0] - xyz[a * 3],
                    p[1] - xyz[a * 3 + 1],
                    p[2] - xyz[a * 3 + 2],
                );
                dx * dx + dy * dy + dz * dz <= r2[a]
            };
            // start with the neighbor that buried the previous point
            if !buried_by(nbi[current_nb]) {
                match nbi.iter().position(|&a| buried_by(a)) {
                    Some(k) => current_nb = k,
                    None => n_surface += 1,
                }
            }
        }
        sasa[i] = (4.0 * PI * ri * ri * n_surface as f64) / n_points as f64;
    }
    sasa
}
