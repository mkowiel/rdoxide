//! `freesasa/src/sasa_lr.c`: Lee & Richards' slicing algorithm.

use super::nb::NbList;
use std::f64::consts::PI;

const TWOPI: f64 = 2.0 * PI;

/// `sort_arcs`: insertion sort of the (start, end) pairs by start
fn sort_arcs(arc: &mut [f64], n: usize) {
    for i in 1..n {
        let (t0, t1) = (arc[2 * i], arc[2 * i + 1]);
        let mut j = i;
        while j > 0 && arc[2 * (j - 1)] > t0 {
            arc[2 * j] = arc[2 * (j - 1)];
            arc[2 * j + 1] = arc[2 * (j - 1) + 1];
            j -= 1;
        }
        arc[2 * j] = t0;
        arc[2 * j + 1] = t1;
    }
}

/// `exposed_arc_length`
fn exposed_arc_length(arc: &mut [f64], n: usize) -> f64 {
    if n == 0 {
        return TWOPI;
    }
    sort_arcs(arc, n);
    let mut sum = arc[0];
    let mut sup = arc[1];
    let mut i2 = 2;
    while i2 < 2 * n {
        if sup < arc[i2] {
            sum += arc[i2] - sup;
        }
        let tmp = arc[i2 + 1];
        if tmp > sup {
            sup = tmp;
        }
        i2 += 2;
    }
    sum + TWOPI - sup
}

/// `freesasa_lee_richards`
pub fn lee_richards(xyz: &[f64], radii: &[f64], probe_radius: f64, n_slices: usize) -> Vec<f64> {
    let n_atoms = xyz.len() / 3;
    let r: Vec<f64> = radii.iter().map(|&r| r + probe_radius).collect();
    let adj = NbList::new(xyz, &r);
    let max_nni = adj.nb.iter().map(Vec::len).max().unwrap_or(0);
    let mut arc = vec![0.0f64; 4 * max_nni];
    let mut sasa = vec![0.0; n_atoms];
    for i in 0..n_atoms {
        sasa[i] = atom_area(xyz, &r, &adj, i, n_slices, &mut arc);
    }
    sasa
}

/// `atom_area`
fn atom_area(v: &[f64], r: &[f64], adj: &NbList, i: usize, ns: usize, arc: &mut [f64]) -> f64 {
    let nbi = &adj.nb[i];
    let (xydi, xdi, ydi) = (&adj.xyd[i], &adj.xd[i], &adj.yd[i]);
    let (zi, ri) = (v[3 * i + 2], r[i]);
    let z_nb: Vec<f64> = nbi.iter().map(|&a| v[3 * a + 2]).collect();
    let r_nb: Vec<f64> = nbi.iter().map(|&a| r[a]).collect();
    let mut sasa = 0.0;
    let delta = 2.0 * ri / ns as f64;
    let mut z = zi - ri - 0.5 * delta;
    for _ in 0..ns {
        z += delta;
        let di = (zi - z).abs();
        let ri_prime2 = ri * ri - di * di;
        if ri_prime2 < 0.0 {
            continue; // round-off errors
        }
        let ri_prime = ri_prime2.sqrt();
        if ri_prime <= 0.0 {
            continue;
        }
        let mut n_arcs = 0usize;
        let mut is_buried = false;
        for j in 0..nbi.len() {
            let dj = (z_nb[j] - z).abs();
            let rj = r_nb[j];
            if dj < rj {
                let rj_prime2 = rj * rj - dj * dj;
                let rj_prime = rj_prime2.sqrt();
                let dij = xydi[j];
                if dij >= ri_prime + rj_prime {
                    continue; // atoms aren't in contact
                }
                if dij + ri_prime < rj_prime {
                    is_buried = true; // circle i is completely inside j
                    break;
                }
                if dij + rj_prime < ri_prime {
                    continue; // circle j is completely inside i
                }
                // arc of circle i intersected by circle j
                let alpha = ((ri_prime2 + dij * dij - rj_prime2) / (2.0 * ri_prime * dij)).acos();
                // position of the mid-point of the intersection along circle i
                let beta = ydi[j].atan2(xdi[j]) + PI;
                let mut inf = beta - alpha;
                let mut sup = beta + alpha;
                if inf < 0.0 {
                    inf += TWOPI;
                }
                if sup > 2.0 * PI {
                    sup -= TWOPI;
                }
                let narc2 = 2 * n_arcs;
                if sup < inf {
                    // the arc passes 2*PI: split in two
                    arc[narc2] = 0.0;
                    arc[narc2 + 1] = sup;
                    arc[narc2 + 2] = inf;
                    arc[narc2 + 3] = TWOPI;
                    n_arcs += 2;
                } else {
                    arc[narc2] = inf;
                    arc[narc2 + 1] = sup;
                    n_arcs += 1;
                }
            }
        }
        if !is_buried {
            sasa += delta * ri * exposed_arc_length(arc, n_arcs);
        }
    }
    sasa
}
