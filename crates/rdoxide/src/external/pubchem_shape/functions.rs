//! `shape_functions1.cpp`, `shape_functions2.cpp` and `shape_math_functions.hpp` of pubchem-align3d (NCBI, public domain):
//! Gaussian overlaps and their gradients, quaternion and rotation-translation matrix algebra, the steric frame.
//! Coordinates are `f32` and all arithmetic is done in `f64` exactly where the C++ code does it.

use std::collections::{BTreeMap, BTreeSet};

use std::f64::consts::PI;

pub type ColorMap = BTreeMap<u32, Vec<u32>>;

const P: f64 = 2.7;
const P2: f64 = P * P;

/// `kappa = pi / lambda^(2/3)` with `lambda = 4 pi / (3 p)`
fn kappa() -> f64 {
    let lambda = (4.0 * PI) / (3.0 * P);
    PI / lambda.powf(2.0 / 3.0)
}

/// `default_dCutOff` squared: the cut-off RDKit's wrapper switches on (`setUseCutOff(true)`)
pub const D2_CUTOFF: f64 = 4.5 * 4.5;

/// `expApproxDouble`: (1 + x/1024)^1024
pub fn exp_approx(x: f64) -> f64 {
    const C: f64 = 0.0009765625;
    let mut y = 1.0 + C * x;
    for _ in 0..10 {
        y *= y;
    }
    y
}

/// `getAlpha`
pub fn get_alpha(r: f64) -> f64 {
    if r > 0.0 { kappa() / (r * r) } else { 0.0 }
}

/// `getA_ak`
pub fn get_a_ak(alpha1: f64, alpha2: f64) -> f64 {
    let sum = alpha1 + alpha2;
    if sum > 0.0 {
        alpha1 * alpha2 / sum
    } else {
        0.0
    }
}

/// `getOverlap0`
pub fn get_overlap0(alpha1: f64, alpha2: f64) -> f64 {
    let sum = alpha1 + alpha2;
    let a_pi = if sum > 0.0 { PI / sum } else { 0.0 };
    P2 * a_pi * a_pi.sqrt()
}

/// `setAlpha`
pub fn set_alpha(rad: &[f64]) -> Vec<f64> {
    rad.iter().map(|&r| get_alpha(r)).collect()
}

/// `getJointColorTypeSet`
pub fn get_joint_color_type_set(ref_types: &[u32], fit_types: &[u32]) -> BTreeSet<u32> {
    let ref_set: BTreeSet<u32> = ref_types.iter().copied().filter(|&t| t > 0).collect();
    fit_types
        .iter()
        .copied()
        .filter(|t| *t > 0 && ref_set.contains(t))
        .collect()
}

/// `getVolumeAtomIndexVector`: the atoms without a colour type
pub fn get_volume_atom_index_vector(types: &[u32]) -> Vec<u32> {
    types
        .iter()
        .enumerate()
        .filter(|(_, t)| **t == 0)
        .map(|(k, _)| k as u32)
        .collect()
}

/// `getColorAtomType2IndexVectorMap`
pub fn get_color_atom_type_to_index_vector_map(types: &[u32]) -> ColorMap {
    let mut m = ColorMap::new();
    for (k, &t) in types.iter().enumerate() {
        if t > 0 {
            m.entry(t).or_default().push(k as u32);
        }
    }
    m
}

/// `restrictColorAtomType2IndexVectorMap`
pub fn restrict_color_map(map: &mut ColorMap, joint: &BTreeSet<u32>) {
    map.retain(|t, _| joint.contains(t));
}

/// the point `m * (x, y, z)` (`m` is 3x3 rotation row-major, then translation)
fn transform(m: &[f64; 12], x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    (
        m[9] + m[0] * x + m[1] * y + m[2] * z,
        m[10] + m[3] * x + m[4] * y + m[5] * z,
        m[11] + m[6] * x + m[7] * y + m[8] * z,
    )
}

/// `ComputeOverlapAndAnalyticGradient`: the shape overlap and, when `grad` is given, its gradient
pub fn overlap_and_gradient(
    m: Option<&[f64; 12]>,
    ref_coord: &[f32],
    alpha_ref: &[f64],
    ref_vol: &[u32],
    fit_coord: &[f32],
    alpha_fit: &[f64],
    fit_vol: &[u32],
    mut grad: Option<&mut [f64; 7]>,
) -> f64 {
    let mut o = 0.0;
    if let Some(g) = grad.as_deref_mut() {
        *g = [0.0; 7];
    }
    for &i in fit_vol {
        let i = i as usize;
        let (ox, oy, oz) = (
            f64::from(fit_coord[i * 3]),
            f64::from(fit_coord[i * 3 + 1]),
            f64::from(fit_coord[i * 3 + 2]),
        );
        let (x, y, z) = match m {
            None => (ox, oy, oz),
            Some(m) => transform(m, ox, oy, oz),
        };
        for &j in ref_vol {
            let j = j as usize;
            let mut x_diff = f64::from(ref_coord[j * 3]) - x;
            let mut y_diff = f64::from(ref_coord[j * 3 + 1]) - y;
            let mut z_diff = f64::from(ref_coord[j * 3 + 2]) - z;
            let d2 = x_diff * x_diff + y_diff * y_diff + z_diff * z_diff;
            if d2 < D2_CUTOFF {
                let alpha1 = alpha_fit[i];
                let alpha2 = alpha_ref[j];
                let a_ak = get_a_ak(alpha1, alpha2);
                let a_oe = get_overlap0(alpha1, alpha2);
                let exp_minus = exp_approx(-d2 * a_ak);
                o += a_oe * exp_minus;
                if let Some(g) = grad.as_deref_mut() {
                    let derivative0 = -2.0 * a_oe * a_ak;
                    x_diff *= derivative0 * exp_minus;
                    y_diff *= derivative0 * exp_minus;
                    z_diff *= derivative0 * exp_minus;
                    let (xp2, yp2, zp2) = (x + x, y + y, z + z);
                    g[0] += x_diff * xp2;
                    g[2] += x_diff * zp2;
                    g[3] -= x_diff * yp2;
                    g[4] += x_diff;
                    g[0] += y_diff * yp2;
                    g[1] -= y_diff * zp2;
                    g[3] += y_diff * xp2;
                    g[5] += y_diff;
                    g[0] += z_diff * zp2;
                    g[1] += z_diff * yp2;
                    g[2] -= z_diff * xp2;
                    g[6] += z_diff;
                }
            }
        }
    }
    o
}

/// `ComputeShapeOverlap`
pub fn shape_overlap(
    ref_coord: &[f32],
    alpha_ref: &[f64],
    ref_vol: &[u32],
    fit_coord: &[f32],
    alpha_fit: &[f64],
    fit_vol: &[u32],
) -> f64 {
    overlap_and_gradient(
        None, ref_coord, alpha_ref, ref_vol, fit_coord, alpha_fit, fit_vol, None,
    )
}

/// `ComputeFeatureOverlap`
pub fn feature_overlap(
    ref_coord: &[f32],
    alpha_ref: &[f64],
    ref_map: &ColorMap,
    fit_coord: &[f32],
    alpha_fit: &[f64],
    fit_map: &ColorMap,
    m: Option<&[f64; 12]>,
) -> f64 {
    let mut overlap = 0.0;
    for (atom_type, ref_idx) in ref_map {
        let Some(fit_idx) = fit_map.get(atom_type) else {
            continue;
        };
        for &i in fit_idx {
            let i = i as usize;
            let (ox, oy, oz) = (
                f64::from(fit_coord[i * 3]),
                f64::from(fit_coord[i * 3 + 1]),
                f64::from(fit_coord[i * 3 + 2]),
            );
            let (x, y, z) = match m {
                None => (ox, oy, oz),
                Some(m) => transform(m, ox, oy, oz),
            };
            for &j in ref_idx {
                let j = j as usize;
                let x_diff = x - f64::from(ref_coord[j * 3]);
                let y_diff = y - f64::from(ref_coord[j * 3 + 1]);
                let z_diff = z - f64::from(ref_coord[j * 3 + 2]);
                let d2 = x_diff * x_diff + y_diff * y_diff + z_diff * z_diff;
                if d2 < D2_CUTOFF {
                    let alpha1 = alpha_fit[i];
                    let alpha2 = alpha_ref[j];
                    let a_ak = get_a_ak(alpha1, alpha2);
                    let a_oe = get_overlap0(alpha1, alpha2);
                    overlap += a_oe * exp_approx(-d2 * a_ak);
                }
            }
        }
    }
    overlap
}

/// `VQuaternionTransArray2RotTransMatrix` (double)
pub fn quaternion_trans_to_matrix(q: &[f64; 7], normalize: bool) -> [f64; 12] {
    let (mut q00, mut q01, mut q02, mut q03) = (q[0], q[1], q[2], q[3]);
    if normalize {
        let s = 1.0 / (q00 * q00 + q01 * q01 + q02 * q02 + q03 * q03).sqrt();
        q00 *= s;
        q01 *= s;
        q02 *= s;
        q03 *= s;
    }
    let (sqw, sqx, sqy, sqz) = (q00 * q00, q01 * q01, q02 * q02, q03 * q03);
    let (qxqy, qwqz) = (q01 * q02, q00 * q03);
    let (qxqz, qwqy) = (q01 * q03, q00 * q02);
    let (qyqz, qwqx) = (q02 * q03, q00 * q01);
    [
        sqw + sqx - sqy - sqz,
        2.0 * (qxqy - qwqz),
        2.0 * (qxqz + qwqy),
        2.0 * (qxqy + qwqz),
        sqw - sqx + sqy - sqz,
        2.0 * (qyqz - qwqx),
        2.0 * (qxqz - qwqy),
        2.0 * (qyqz + qwqx),
        sqw - sqx - sqy + sqz,
        q[4],
        q[5],
        q[6],
    ]
}

/// `VCombineRotTransMatrixMatrix`
pub fn combine_matrices(t1: &[f64; 12], t2: &[f64; 12]) -> [f64; 12] {
    [
        t1[0] * t2[0] + t1[1] * t2[3] + t1[2] * t2[6],
        t1[0] * t2[1] + t1[1] * t2[4] + t1[2] * t2[7],
        t1[0] * t2[2] + t1[1] * t2[5] + t1[2] * t2[8],
        t1[3] * t2[0] + t1[4] * t2[3] + t1[5] * t2[6],
        t1[3] * t2[1] + t1[4] * t2[4] + t1[5] * t2[7],
        t1[3] * t2[2] + t1[4] * t2[5] + t1[5] * t2[8],
        t1[6] * t2[0] + t1[7] * t2[3] + t1[8] * t2[6],
        t1[6] * t2[1] + t1[7] * t2[4] + t1[8] * t2[7],
        t1[6] * t2[2] + t1[7] * t2[5] + t1[8] * t2[8],
        t1[0] * t2[9] + t1[1] * t2[10] + t1[2] * t2[11] + t1[9],
        t1[3] * t2[9] + t1[4] * t2[10] + t1[5] * t2[11] + t1[10],
        t1[6] * t2[9] + t1[7] * t2[10] + t1[8] * t2[11] + t1[11],
    ]
}

/// `VRotTransMatrix2QuaternionTransArray`
pub fn matrix_to_quaternion_trans(m: &[f64; 12]) -> [f64; 7] {
    let trace = 1.0 + m[0] + m[4] + m[8];
    let (qw, qx, qy, qz);
    if trace > 0.00001 {
        let s = 2.0 * trace.sqrt();
        qw = 0.25 * s;
        let si = 1.0 / s;
        qx = (m[7] - m[5]) * si;
        qy = (m[2] - m[6]) * si;
        qz = (m[3] - m[1]) * si;
    } else if m[0] > m[4] && m[0] > m[8] {
        let s = 2.0 * (1.0 + m[0] - m[4] - m[8]).sqrt();
        let si = 1.0 / s;
        qw = (m[7] - m[5]) * si;
        qx = 0.25 * s;
        qy = (m[3] + m[1]) * si;
        qz = (m[2] + m[6]) * si;
    } else if m[4] > m[8] {
        let s = 2.0 * (1.0 + m[4] - m[0] - m[8]).sqrt();
        let si = 1.0 / s;
        qw = (m[2] - m[6]) * si;
        qx = (m[3] + m[1]) * si;
        qy = 0.25 * s;
        qz = (m[7] + m[5]) * si;
    } else {
        let s = 2.0 * (1.0 + m[8] - m[0] - m[4]).sqrt();
        let si = 1.0 / s;
        qw = (m[3] - m[1]) * si;
        qx = (m[2] + m[6]) * si;
        qy = (m[7] + m[5]) * si;
        qz = 0.25 * s;
    }
    let s = 1.0 / (qw * qw + qx * qx + qy * qy + qz * qz).sqrt();
    [qw * s, qx * s, qy * s, qz * s, m[9], m[10], m[11]]
}

/// `VApplyRotTransMatrix(c, d, n, m)`: the transformed copy of `d` (`f32` matrix, `f64` arithmetic)
pub fn apply_rot_trans_matrix(d: &[f32], m: &[f32; 12]) -> Vec<f32> {
    let m: Vec<f64> = m.iter().map(|&v| f64::from(v)).collect();
    let mut c = vec![0f32; d.len()];
    for i in (0..d.len()).step_by(3) {
        let (x, y, z) = (f64::from(d[i]), f64::from(d[i + 1]), f64::from(d[i + 2]));
        c[i] = (m[9] + m[0] * x + m[1] * y + m[2] * z) as f32;
        c[i + 1] = (m[10] + m[3] * x + m[4] * y + m[5] * z) as f32;
        c[i + 2] = (m[11] + m[6] * x + m[7] * y + m[8] * z) as f32;
    }
    c
}

/// `Get3dStericCenterOfGravity` for `f32`
pub fn steric_center(coord: &[f32]) -> [f32; 3] {
    let n = coord.len() / 3;
    if n == 0 {
        return [0.0; 3];
    }
    let (mut x, mut y, mut z) = (coord[0], coord[1], coord[2]);
    for i in 1..n {
        x += coord[i * 3];
        y += coord[1 + i * 3];
        z += coord[2 + i * 3];
    }
    let nf = n as f32;
    [x / nf, y / nf, z / nf]
}

/// `VStericFrameRotation`: the rotation into the inertial frame (first 9 entries) and the eigenvalues
pub fn steric_frame_rotation(coord: &[f32]) -> ([f32; 12], [f32; 3]) {
    let n = coord.len() / 3;
    let mut it = [[0.0f64; 3]; 3];
    let (mut sx, mut sy, mut sz) = (
        f64::from(coord[0]),
        f64::from(coord[1]),
        f64::from(coord[2]),
    );
    let (mut xx, mut yy, mut zz) = (sx * sx, sy * sy, sz * sz);
    it[0][0] = yy + zz;
    it[1][1] = xx + zz;
    it[2][2] = xx + yy;
    it[0][1] = sx * sy;
    it[0][2] = sx * sz;
    it[1][2] = sy * sz;
    for i in 1..n {
        sx = f64::from(coord[i * 3]);
        sy = f64::from(coord[1 + i * 3]);
        sz = f64::from(coord[2 + i * 3]);
        xx = sx * sx;
        yy = sy * sy;
        zz = sz * sz;
        it[0][0] += yy + zz;
        it[1][1] += xx + zz;
        it[2][2] += xx + yy;
        it[0][1] += sx * sy;
        it[0][2] += sx * sz;
        it[1][2] += sy * sz;
    }
    // the eigenvectors are accumulated in float, as the C++ code stores them in the float matrix
    let mut m = [0f32; 12];
    m[0] = 1.0;
    m[4] = 1.0;
    m[8] = 1.0;
    let mut e = [it[0][0], it[1][1], it[2][2]];
    let max_iter = 50;
    let conv = 0.000001;
    let thresh_val = 1.0 / (5.0 * 3.0 * 3.0);
    if n > 1 {
        let mut iter = 0;
        while iter < max_iter {
            iter += 1;
            let sum = it[0][1].abs() + it[0][2].abs() + it[1][2].abs();
            if sum < conv {
                break;
            }
            let thresh = if iter < 4 { sum * thresh_val } else { 0.0 };
            for p in 0..2 {
                for q in p + 1..3 {
                    let g = 100.0 * it[p][q].abs();
                    if iter > 4 && g < conv {
                        it[p][q] = 0.0;
                    } else if thresh < it[p][q].abs() {
                        let t = if g < conv {
                            it[p][q] / (e[q] - e[p])
                        } else {
                            let theta = (e[q] - e[p]) / (2.0 * it[p][q]);
                            if theta < 0.0 {
                                -1.0 / ((1.0 + theta * theta).sqrt() - theta)
                            } else {
                                1.0 / ((1.0 + theta * theta).sqrt() + theta)
                            }
                        };
                        let c = 1.0 / (1.0 + t * t).sqrt();
                        let s = t * c;
                        let tau = s / (1.0 + c);
                        let tz = t * it[p][q];
                        it[p][q] = 0.0;
                        for r in 0..p {
                            let (g, h) = (it[r][p], it[r][q]);
                            it[r][p] = g - s * (h + g * tau);
                            it[r][q] = h + s * (g - h * tau);
                        }
                        for r in p + 1..q {
                            let (g, h) = (it[p][r], it[r][q]);
                            it[p][r] = g - s * (h + g * tau);
                            it[r][q] = h + s * (g - h * tau);
                        }
                        for r in q + 1..3 {
                            let (g, h) = (it[p][r], it[q][r]);
                            it[p][r] = g - s * (h + g * tau);
                            it[q][r] = h + s * (g - h * tau);
                        }
                        e[p] -= tz;
                        e[q] += tz;
                        for r in 0..3 {
                            let (g, h) = (f64::from(m[r * 3 + p]), f64::from(m[r * 3 + q]));
                            m[r * 3 + p] = (g - s * (h + g * tau)) as f32;
                            m[r * 3 + q] = (h + s * (g - h * tau)) as f32;
                        }
                    }
                }
            }
        }
        // sort eigenvalues and eigenvectors
        let swap_rows = |m: &mut [f32; 12], a: usize, b: usize| {
            for k in 0..3 {
                let tmp = -m[a * 3 + k];
                m[a * 3 + k] = m[b * 3 + k];
                m[b * 3 + k] = tmp;
            }
        };
        let mut diff = e[0] - e[1];
        if diff > 0.0001 {
            e.swap(0, 1);
            swap_rows(&mut m, 0, 1);
        }
        diff = e[1] - e[2];
        if diff > 0.0001 {
            e.swap(1, 2);
            swap_rows(&mut m, 1, 2);
            diff = e[0] - e[1];
            if diff > 0.0001 {
                e.swap(0, 1);
                swap_rows(&mut m, 0, 1);
            }
        }
    }
    (m, [e[0] as f32, e[1] as f32, e[2] as f32])
}

/// `VTransform2StericFrame`: rotates and centres `coord`; returns the matrix used and the eigenvalues
pub fn transform_to_steric_frame(coord: &mut Vec<f32>) -> ([f32; 12], [f32; 3]) {
    let (mut m, ev) = steric_frame_rotation(coord);
    *coord = apply_rot_trans_matrix(coord, &m);
    let c = steric_center(coord);
    m[9] = -c[0];
    m[10] = -c[1];
    m[11] = -c[2];
    for i in (0..coord.len()).step_by(3) {
        coord[i] += m[9];
        coord[i + 1] += m[10];
        coord[i + 2] += m[11];
    }
    (m, ev)
}

/// `CalculateQrat`
pub fn calculate_qrat(ev: &[f32; 3], threshold: f32) -> u32 {
    let mut d = [
        ev[1] + ev[2] - ev[0],
        ev[0] + ev[2] - ev[1],
        ev[0] + ev[1] - ev[2],
    ];
    d.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let mut qrat = 1000;
    if d[1] > 0.0 {
        let u_rqyx = u32::from(threshold < (d[1] / d[0]));
        let u_rqzy = u32::from(threshold < (d[2] / d[1]));
        qrat = u_rqyx + u_rqzy;
    }
    qrat
}

/// `TransformCoordToStericFrameAndCalculateQrat`
pub fn qrat_of(coord: &[f32], threshold: f32) -> u32 {
    let mut c = coord.to_vec();
    let (_, ev) = transform_to_steric_frame(&mut c);
    calculate_qrat(&ev, threshold)
}
