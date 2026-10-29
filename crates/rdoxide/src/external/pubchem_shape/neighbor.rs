//! `shape_neighbor.cpp` of pubchem-align3d (NCBI, public domain): `Neighbor_Conformers`, the two-phase optimization of the
//! overlay of a fit shape on a reference shape starting from up to 16 poses.

use super::consts::{MATRIX, PMATRIX, PQMATRIX, QMATRIX};
use super::functions::{
    ColorMap, combine_matrices, feature_overlap, matrix_to_quaternion_trans, overlap_and_gradient,
    qrat_of, quaternion_trans_to_matrix,
};

/// the data of one shape the optimization needs
#[derive(Debug)]
pub struct ShapeData<'a> {
    pub coord: &'a [f32],
    pub alpha: &'a [f64],
    pub volume_atoms: &'a [u32],
    pub color_map: &'a ColorMap,
    pub sov: f64,
    pub sof: f64,
}

/// the result: the rotation/translation matrix and the shape and colour Tanimotos
#[derive(Debug)]
pub struct NeighborResult {
    pub matrix: [f32; 12],
    pub st: f64,
    pub ct: f64,
}

const MAX_QUATERNION_STEP: f64 = 0.075;
const MAX_TRANSLATION_STEP: f64 = 0.500;
const MIN_QUATERNION_STEP: f64 = 0.0002;
const MIN_TRANSLATION_STEP: f64 = 0.0020;
const MIN_VOLUME_OVERLAP: f64 = 0.0020;

fn signbit(x: f64) -> bool {
    x.is_sign_negative()
}

/// the C++ code copies `7 * sizeof(float)` bytes where it means 7 doubles: the first three doubles and the low half of
/// the fourth (little endian)
fn copy_seven_floats(dst: &mut [f64; 7], src: &[f64; 7]) {
    dst[..3].copy_from_slice(&src[..3]);
    dst[3] = f64::from_bits(
        (dst[3].to_bits() & 0xFFFF_FFFF_0000_0000) | (src[3].to_bits() & 0x0000_0000_FFFF_FFFF),
    );
}

struct Ctx<'a> {
    r: &'a ShapeData<'a>,
    f: &'a ShapeData<'a>,
}

impl Ctx<'_> {
    fn overlap(&self, m: &[f64; 12], g: Option<&mut [f64; 7]>) -> f64 {
        overlap_and_gradient(
            Some(m),
            self.r.coord,
            self.r.alpha,
            self.r.volume_atoms,
            self.f.coord,
            self.f.alpha,
            self.f.volume_atoms,
            g,
        )
    }
}

/// `Neighbor_Conformers`
pub fn neighbor_conformers(
    r: &ShapeData,
    f: &ShapeData,
    consider_color_atoms: bool,
    max_preiters: u32,
    max_postiters: u32,
    opt_param: f64,
) -> NeighborResult {
    let ctx = Ctx { r, f };
    let st_sum_sov = r.sov + f.sov;
    let ct_sum_sov = r.sof + f.sof;
    let _ = (ctx.r.sov, ctx.f.sov);

    let qrat_threshold = 0.7225f32;
    let ref_qrat = qrat_of(&r.coord[..r.alpha.len() * 3], qrat_threshold);
    let fit_qrat = qrat_of(&f.coord[..f.alpha.len() * 3], qrat_threshold);
    let is_qrat_nonzero = ref_qrat > 0 || fit_qrat > 0;
    let npose = if is_qrat_nonzero { 16 } else { 4 };

    let mut nbr_st = 0.0f64;
    let mut nbr_ct = 0.0f64;
    let mut nbr_comb_prelim = 0.0f64;
    let mut nbr_comb = 0.0f64;
    let mut best_qom = [0.0f64; 7];

    let mut best_qom_pose = [[0.0f64; 7]; 16];
    let mut best_mom_pose = [[0.0f64; 12]; 16];

    // the Tanimoto-style objective of a pose
    let objective = |st_tot: f64, mc_q: &[f64; 7]| -> (f64, f64, f64, f64) {
        let mut fst = 0.0;
        if st_sum_sov > 0.0 && st_sum_sov > st_tot {
            fst = st_tot / (st_sum_sov - st_tot);
        }
        let mut fct = 0.0;
        if consider_color_atoms && ct_sum_sov > 0.0 {
            let mc = quaternion_trans_to_matrix(mc_q, false);
            let ct_tot = feature_overlap(
                r.coord,
                r.alpha,
                r.color_map,
                f.coord,
                f.alpha,
                f.color_map,
                Some(&mc),
            );
            if ct_sum_sov > ct_tot {
                fct = ct_tot / (ct_sum_sov - ct_tot);
            }
        }
        let comb = if consider_color_atoms {
            opt_param * (fst - fct) + fct
        } else {
            fst
        };
        (fst, fct, comb, 0.0)
    };

    for ibest in 0..npose {
        if is_qrat_nonzero {
            best_qom_pose[ibest] = QMATRIX[ibest];
            let mut mm = [0.0f64; 12];
            mm.copy_from_slice(&MATRIX[ibest]);
            best_mom_pose[ibest] = mm;
        } else {
            let mut m_loc = 0;
            let mut best_f = -1.0f64;
            for m_loc_ in 0..7 {
                let qom: [f64; 7] = {
                    let mut q = [0.0; 7];
                    q.copy_from_slice(&PQMATRIX[ibest][m_loc_][..7]);
                    q
                };
                let mom = &PMATRIX[ibest][m_loc_];
                let st_tot = ctx.overlap(mom, None);
                let (_, _, f_, _) = objective(st_tot, &qom);
                if f_ > best_f {
                    m_loc = m_loc_;
                    best_f = f_;
                }
            }
            best_qom_pose[ibest].copy_from_slice(&PQMATRIX[ibest][m_loc][..7]);
            best_mom_pose[ibest] = PMATRIX[ibest][m_loc];
        }
    }

    let mut old_quattrans_pose = [[0.0f64; 7]; 16];
    let mut old_gradient_pose = [[0.0f64; 7]; 16];
    let mut qstep_size_pose = [0.0f64; 16];
    let mut tstep_size_pose = [0.0f64; 16];
    let mut allowed_pose = [false; 16];
    let mut st_tot_pose = [0.0f64; 16];
    let mut fcomb_best_pose = [0.0f64; 16];

    for ibest in 0..npose {
        allowed_pose[ibest] = true;
        fcomb_best_pose[ibest] = 0.0;
        qstep_size_pose[ibest] = -0.0010;
        tstep_size_pose[ibest] = -0.0100;
        let mut g = [0.0; 7];
        st_tot_pose[ibest] = ctx.overlap(&best_mom_pose[ibest], Some(&mut g));
        old_gradient_pose[ibest] = g;
    }

    for icycle in 0..2 {
        for ibest in 0..npose {
            let mut g = old_gradient_pose[ibest];
            let mut qstep_size = qstep_size_pose[ibest];
            let mut tstep_size = tstep_size_pose[ibest];
            let mut old_gradient = old_gradient_pose[ibest];
            let mut old_quattrans = old_quattrans_pose[ibest];
            let mut st_tot = st_tot_pose[ibest];
            let mut qom = best_qom_pose[ibest];

            if icycle == 1 && fcomb_best_pose[ibest] <= 0.7 * nbr_comb_prelim {
                allowed_pose[ibest] = false;
            }

            if allowed_pose[ibest] {
                let (start_iter, end_iter) = if icycle == 0 {
                    (0, max_preiters)
                } else {
                    (max_preiters, max_preiters + max_postiters)
                };
                for iter in start_iter..end_iter {
                    let mut step = [0.0f64; 7];
                    let mut converged = false;
                    let size = |i: usize, q: f64, t: f64| if i <= 3 { q } else { t };
                    if iter == 0 {
                        for i in 1..=3 {
                            step[i] = qstep_size * g[i];
                        }
                        for i in 4..=6 {
                            step[i] = tstep_size * g[i];
                        }
                    } else {
                        for i in 1..=6 {
                            let ss = size(i, qstep_size, tstep_size);
                            if signbit(g[i]) != signbit(old_gradient[i]) {
                                step[i] = (((qom[i] * old_gradient[i].abs())
                                    + (old_quattrans[i] * g[i].abs()))
                                    / (old_gradient[i].abs() + g[i].abs()))
                                    - qom[i];
                                let new_step = ss * g[i];
                                if step[i].abs() > new_step.abs() {
                                    step[i] *= (new_step / step[i]).abs();
                                }
                            } else {
                                step[i] = ss * g[i];
                            }
                        }
                    }
                    let old_overlap = st_tot;
                    old_quattrans = qom;
                    old_gradient = g;

                    let mut line_iter = 0u32;
                    while !converged {
                        let mtstep = step[4].abs().max(step[5].abs()).max(step[6].abs());
                        if mtstep > MAX_TRANSLATION_STEP {
                            let tscale = MAX_TRANSLATION_STEP / mtstep;
                            for i in 4..=6 {
                                if step[i].abs() > MAX_TRANSLATION_STEP {
                                    step[i] *= tscale;
                                }
                            }
                        }
                        let mqstep = step[1].abs().max(step[2].abs()).max(step[3].abs());
                        if mqstep > MAX_QUATERNION_STEP {
                            let qscale = MAX_QUATERNION_STEP / mqstep;
                            for i in 1..=3 {
                                if step[i].abs() > MAX_QUATERNION_STEP {
                                    step[i] *= qscale;
                                }
                            }
                        }
                        if MIN_QUATERNION_STEP > mqstep && MIN_TRANSLATION_STEP > mtstep {
                            st_tot = 0.0;
                            converged = true;
                            break;
                        }

                        let quat_squared =
                            step[1] * step[1] + step[2] * step[2] + step[3] * step[3];
                        step[0] = (1.0 - quat_squared).sqrt();
                        let tm1 = quaternion_trans_to_matrix(&qom, false);
                        let tm2 = quaternion_trans_to_matrix(&step, false);
                        let mc = combine_matrices(&tm2, &tm1);
                        let qc = matrix_to_quaternion_trans(&mc);
                        st_tot = ctx.overlap(&mc, Some(&mut g));

                        if st_tot > old_overlap {
                            qom = qc;
                            break;
                        }
                        if line_iter == 2 {
                            qstep_size *= 0.1;
                            tstep_size *= 0.1;
                            for i in 1..=3 {
                                step[i] = qstep_size * old_gradient[i];
                            }
                            for i in 4..=6 {
                                step[i] = tstep_size * old_gradient[i];
                            }
                            qstep_size *= 5.0;
                            tstep_size *= 5.0;
                        } else if line_iter > 3 {
                            converged = true;
                            qom = qc;
                            break;
                        } else {
                            for i in 1..=6 {
                                if signbit(g[i]) != signbit(old_gradient[i]) {
                                    step[i] = (((qc[i] * old_gradient[i].abs())
                                        + (old_quattrans[i] * g[i].abs()))
                                        / (old_gradient[i].abs() + g[i].abs()))
                                        - qc[i];
                                    let new_step = qstep_size * g[i];
                                    if step[i].abs() > new_step.abs() {
                                        step[i] *= (new_step / step[i]).abs();
                                    }
                                } else if 1.0 > g[i].abs() {
                                    step[i] = qstep_size * g[i];
                                } else if g[i].abs() > old_gradient[i].abs() {
                                    step[i] += qstep_size * g[i];
                                } else {
                                    let mut delta = g[i] * (step[i] / (old_gradient[i] - g[i]));
                                    let (frac, floor, scale) = if i <= 3 {
                                        (0.1, 0.001, 0.0005)
                                    } else {
                                        (0.5, 0.01, 0.005)
                                    };
                                    if delta.abs() > (step[i] * frac).abs() && delta.abs() > floor {
                                        delta *= scale / delta.abs();
                                    }
                                    step[i] += delta;
                                }
                            }
                        }
                        qom = old_quattrans;
                        line_iter += 1;
                    }
                    if converged || MIN_VOLUME_OVERLAP > (st_tot - old_overlap) {
                        if old_overlap > st_tot {
                            st_tot = old_overlap;
                            copy_seven_floats(&mut qom, &old_quattrans);
                        }
                        break;
                    }
                }

                let (fst, fct, fcomb, _) = objective(st_tot, &qom);
                if fcomb_best_pose[ibest] < fcomb {
                    fcomb_best_pose[ibest] = fcomb;
                }
                if icycle == 0 {
                    if nbr_comb_prelim < fcomb {
                        nbr_comb_prelim = fcomb;
                    }
                } else if nbr_comb < fcomb {
                    best_qom = qom;
                    nbr_st = fst;
                    nbr_ct = fct;
                    nbr_comb = fcomb;
                }
            }

            old_gradient_pose[ibest] = old_gradient;
            old_quattrans_pose[ibest] = old_quattrans;
            qstep_size_pose[ibest] = qstep_size;
            tstep_size_pose[ibest] = tstep_size;
            st_tot_pose[ibest] = st_tot;
            best_qom_pose[ibest] = qom;
            let _ = g;
        }
    }

    let m = quaternion_trans_to_matrix(&best_qom, false);
    let mut matrix = [0f32; 12];
    for (o, v) in matrix.iter_mut().zip(m) {
        *o = v as f32;
    }
    NeighborResult {
        matrix,
        st: nbr_st,
        ct: nbr_ct,
    }
}
