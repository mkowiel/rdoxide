//! `GaussianShape.h`: overlaying and scoring shapes and molecules.

use super::alignment::{QuatTrans, SingleConformerAlignment};
use super::options::{OptimMode, ShapeOverlayOptions, StartMode};
use super::shape_input::{ShapeInput, ShapeInputOptions};
use crate::Result;
use crate::geometry::{Point3D, Transform3D};
use crate::graphmol::mol_transforms::transform_conformer;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::rdgeneral::threads::par_map;

/// Moves the fit to its initial centroid and principal axes, then into the
/// frame of the reference.
fn compute_final_transform(
    in_ref_trans: &[f64; 3],
    in_ref_rot: &[f64; 9],
    in_fit_trans: &[f64; 3],
    in_fit_rot: &[f64; 9],
    ov_xform: &Transform3D,
) -> Transform3D {
    // Move to fitShape's initial centroid and principal axes
    let mut transform0 = Transform3D::identity();
    transform0.set_translation(&Point3D::new(
        in_fit_trans[0],
        in_fit_trans[1],
        in_fit_trans[2],
    ));

    let mut transform1 = Transform3D::identity();
    for i in 0..3 {
        for j in 0..3 {
            transform1.set_val(i, j, in_fit_rot[3 * i + j]);
        }
    }

    // Rotate by the inverse of the ref shape's canonical rotation and
    // translate by the negative of its canonical translation.
    let mut to_ref_ref_frame = Transform3D::identity();
    for i in 0..3 {
        for j in 0..3 {
            to_ref_ref_frame.set_val(i, j, in_ref_rot[3 * j + i]);
        }
        to_ref_ref_frame.set_val(i, 3, -in_ref_trans[i]);
    }

    to_ref_ref_frame
        .compose(ov_xform)
        .compose(&transform1)
        .compose(&transform0)
}

/// The 16 rotations of `ROTATE_0`, `ROTATE_180` and `ROTATE_45`.
fn plain_quats() -> [[f64; 4]; 16] {
    let s = std::f64::consts::FRAC_PI_4.sin();
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [s, -s, 0.0, 0.0],
        [s, s, 0.0, 0.0],
        [0.0, 0.0, -s, s],
        [0.0, 0.0, s, s],
        [s, 0.0, 0.0, -s],
        [0.0, s, s, 0.0],
        [s, 0.0, 0.0, s],
        [0.0, -s, s, 0.0],
        [s, 0.0, s, 0.0],
        [0.0, s, 0.0, s],
        [0.0, -s, 0.0, s],
        [s, 0.0, -s, 0.0],
    ]
}

/// The 28 rotations of `ROTATE_180_WIGGLE`: 4 groups of 7.
fn wiggle_quats() -> [[f64; 4]; 28] {
    let (q1, q) = (0.977659114061, 0.210196709523); // 0.215 (un-normalized)
    [
        [1.0, 0.0, 0.0, 0.0], //  0   X,  Y,  Z
        [q1, q, 0.0, 0.0],
        [q1, -q, 0.0, 0.0],
        [q1, 0.0, q, 0.0],
        [q1, 0.0, -q, 0.0],
        [q1, 0.0, 0.0, q],
        [q1, 0.0, 0.0, -q],
        [0.0, 1.0, 0.0, 0.0], //  1   X, -Y, -Z
        [q, q1, 0.0, 0.0],
        [q, -q1, 0.0, 0.0],
        [0.0, q1, q, 0.0],
        [0.0, q1, -q, 0.0],
        [0.0, q1, 0.0, q],
        [0.0, q1, 0.0, -q],
        [0.0, 0.0, 0.0, 1.0], //  2  -X, -Y,  Z
        [q, 0.0, 0.0, q1],
        [q, 0.0, 0.0, -q1],
        [0.0, q, 0.0, q1],
        [0.0, -q, 0.0, q1],
        [0.0, 0.0, q, q1],
        [0.0, 0.0, -q, q1],
        [0.0, 0.0, 1.0, 0.0], //  3  -X,  Y, -Z
        [q, 0.0, q1, 0.0],
        [q, 0.0, -q1, 0.0],
        [0.0, q, q1, 0.0],
        [0.0, -q, q1, 0.0],
        [0.0, 0.0, q1, q],
        [0.0, 0.0, q1, -q],
    ]
}

fn quat_trans(q: &[f64; 4], disp: &Point3D) -> QuatTrans {
    [q[0], q[1], q[2], q[3], disp.x, disp.y, disp.z]
}

fn use_color(opts: &ShapeOverlayOptions) -> bool {
    opts.optim_mode != OptimMode::ShapeOnly
}

fn initial_rotation_plain(
    index: usize,
    ref_shape: &ShapeInput,
    fit_shape: &ShapeInput,
    ref_disp: &Point3D,
    opts: &ShapeOverlayOptions,
) -> ([f64; 4], f64) {
    let quats = plain_quats();
    let mut sca = SingleConformerAlignment::new(
        ref_shape.view(),
        fit_shape.view(),
        quat_trans(&quats[index], ref_disp),
        opts,
    );
    let scores = sca.calc_scores(use_color(opts));
    (quats[index], scores[0])
}

fn initial_rotation_wiggle(
    index: usize,
    ref_shape: &ShapeInput,
    fit_shape: &ShapeInput,
    ref_disp: &Point3D,
    opts: &ShapeOverlayOptions,
) -> ([f64; 4], f64) {
    let quats = wiggle_quats();
    let start = index * 7;
    let mut best_quat = 0;
    let mut best_score = 0.0;
    let mut sca = SingleConformerAlignment::new(
        ref_shape.view(),
        fit_shape.view(),
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        opts,
    );
    for i in start..start + 7 {
        sca.set_quat_trans(quat_trans(&quats[i], ref_disp));
        let scores = sca.calc_scores(use_color(opts));
        if scores[0] > best_score {
            best_score = scores[0];
            best_quat = i;
        }
    }
    (quats[best_quat], best_score)
}

/// The translation that moves the fit's extreme point `index - 1` onto the
/// reference's (none for index 0).
fn initial_translation(
    index: usize,
    ref_shape: &ShapeInput,
    ref_extremes: &[usize; 6],
    fit_shape: &ShapeInput,
    fit_extremes: &[usize; 6],
) -> Point3D {
    if !(1..=6).contains(&index) {
        return Point3D::new(0.0, 0.0, 0.0);
    }
    let disp = |shape: &ShapeInput, ext: &[usize; 6]| {
        let c = &shape.coords()[ext[index - 1] * 3..];
        Point3D::new(c[0], c[1], c[2])
    };
    disp(fit_shape, fit_extremes) - disp(ref_shape, ref_extremes)
}

fn calculate_qrat(eigen_values: &[f64; 3]) -> u32 {
    let mut d = [
        eigen_values[1] + eigen_values[2] - eigen_values[0],
        eigen_values[0] + eigen_values[2] - eigen_values[1],
        eigen_values[0] + eigen_values[1] - eigen_values[2],
    ];
    d.sort_by(|a, b| b.total_cmp(a));
    const QRAT_THRESHOLD: f64 = 0.7225; // 0.85*0.85
    let mut qrat = 1000;
    if d[1] > 0.0 {
        let rqyx = u32::from(QRAT_THRESHOLD < d[1] / d[0]);
        let rqzy = u32::from(QRAT_THRESHOLD < d[2] / d[1]);
        qrat = rqyx + rqzy;
    }
    qrat
}

fn decide_start_mode_from_eigen_values(
    ref_shape: &ShapeInput,
    fit_shape: &ShapeInput,
) -> Result<StartMode> {
    // The PubChem code uses the moments of inertia for this, rather than the
    // canonical transformation.
    let rqrat = calculate_qrat(&ref_shape.calc_moments_of_inertia(true)?);
    let fqrat = calculate_qrat(&fit_shape.calc_moments_of_inertia(true)?);
    Ok(if rqrat > 0 || fqrat > 0 {
        StartMode::Rotate45
    } else {
        StartMode::Rotate180Wiggle
    })
}

/// The overlay proper, of 2 shapes already in the frame required.
fn align_shape_impl(
    ref_shape: &mut ShapeInput,
    fit_shape: &mut ShapeInput,
    opts: &ShapeOverlayOptions,
) -> Result<([f64; 3], Transform3D)> {
    let mut start_mode = opts.start_mode;
    if start_mode == StartMode::ALaPubchem {
        start_mode = decide_start_mode_from_eigen_values(ref_shape, fit_shape)?;
    }
    let final_rot_index = match start_mode {
        StartMode::Rotate0 | StartMode::Rotate0Fragment => 1,
        StartMode::Rotate180 | StartMode::Rotate180Fragment | StartMode::Rotate180Wiggle => 4,
        StartMode::Rotate45 | StartMode::Rotate45Fragment => 16,
        StartMode::ALaPubchem => 1,
    };
    let final_trans_index = match start_mode {
        StartMode::Rotate0Fragment | StartMode::Rotate45Fragment | StartMode::Rotate180Fragment => {
            7
        }
        _ => 1,
    };
    let ref_extremes = ref_shape.calc_extremes();
    let fit_extremes = fit_shape.calc_extremes();
    let (ref_shape, fit_shape) = (&*ref_shape, &*fit_shape);

    let mut best_score = [0.0; 3];
    let mut best_total = -1.0;
    let mut best_xform = Transform3D::identity();

    // Get together the start transformations.
    let mut aligners = Vec::with_capacity(final_trans_index * final_rot_index);
    let mut best_score_for_start: Vec<(f64, usize)> = Vec::new();
    let mut k = 0;
    for j in 0..final_trans_index {
        let ref_disp = initial_translation(j, ref_shape, &ref_extremes, fit_shape, &fit_extremes);
        for i in 0..final_rot_index {
            let (quat, score) = if start_mode == StartMode::Rotate180Wiggle {
                initial_rotation_wiggle(i, ref_shape, fit_shape, &ref_disp, opts)
            } else {
                initial_rotation_plain(i, ref_shape, fit_shape, &ref_disp, opts)
            };
            aligners.push(SingleConformerAlignment::new(
                ref_shape.view(),
                fit_shape.view(),
                quat_trans(&quat, &ref_disp),
                opts,
            ));
            best_score_for_start.push((score, k));
            k += 1;
        }
    }

    // Do it in 2 cycles, a quick optimisation first, followed by an
    // additional longer one for those that look like they're going to win.
    for cycle in 0..2 {
        libcxx_sort_by(&mut best_score_for_start, |p1, p2| p1.0 > p2.0);
        let mut next = Vec::with_capacity(best_score_for_start.len());
        for &(bssf, m) in &best_score_for_start {
            if cycle == 1 && bssf < 0.7 * best_score[0] {
                continue;
            }
            let (out, _) = aligners[m].do_overlay(cycle);
            next.push((out[0], m));
            if out[0] > best_total {
                best_total = out[0];
                best_score = [out[0], out[1], out[2]];
                best_xform = aligners[m].final_transform();
            }
        }
        best_score_for_start = next;
    }
    Ok((best_score, best_xform))
}

/// Aligns `fit_shape` onto `ref_shape`, moving the fit shape. Returns the
/// combination score of the shape Tversky value and the colour Tversky value
/// (zero if colours are not used) and the individual values, and the
/// transformation that aligns the fit onto the reference.
pub fn align_shape(
    ref_shape: &ShapeInput,
    fit_shape: &mut ShapeInput,
    opts: &ShapeOverlayOptions,
) -> Result<([f64; 3], Transform3D)> {
    // The shapes aren't necessarily normalized (it's not done on creation, for
    // example) but they might need to be.
    let mut working_ref = ref_shape.clone();
    let mut working_fit = fit_shape.clone();
    let in_ref_trans = working_ref.calc_canonical_translation()?;
    let in_ref_rot = working_ref.calc_canonical_rotation()?;
    let in_fit_trans = working_fit.calc_canonical_translation()?;
    let in_fit_rot = working_fit.calc_canonical_rotation()?;
    // If we're not normalizing, translate both shapes so that the fit is at
    // the origin, so the rotations work.
    let mut move_to_origin = Transform3D::identity();
    let mut move_from_origin = Transform3D::identity();
    if opts.normalize {
        if !working_ref.is_normalized() {
            working_ref.normalize_coords()?;
        }
        if !working_fit.is_normalized() {
            working_fit.normalize_coords()?;
        }
    } else {
        let t = working_fit.calc_canonical_translation()?;
        move_to_origin.set_translation(&Point3D::new(t[0], t[1], t[2]));
        move_from_origin.set_translation(&Point3D::new(-t[0], -t[1], -t[2]));
        working_fit.transform_coords(&move_to_origin);
        working_ref.transform_coords(&move_to_origin);
    }

    let (scores, mut best_xform) = align_shape_impl(&mut working_ref, &mut working_fit, opts)?;
    best_xform = if opts.normalize {
        compute_final_transform(
            &in_ref_trans,
            &in_ref_rot,
            &in_fit_trans,
            &in_fit_rot,
            &best_xform,
        )
    } else {
        // Shove it back again.
        move_from_origin
            .compose(&best_xform)
            .compose(&move_to_origin)
    };
    fit_shape.transform_coords(&best_xform);
    Ok((scores, best_xform))
}

/// Aligns the conformer `fit_conf_id` of `fit` onto the reference shape,
/// moving the conformer.
pub fn align_molecule_to_shape(
    ref_shape: &ShapeInput,
    fit: &mut ROMol,
    fit_opts: &ShapeInputOptions,
    overlay_opts: &ShapeOverlayOptions,
    fit_conf_id: i32,
) -> Result<([f64; 3], Transform3D)> {
    let mut fit_shape = ShapeInput::new(fit, fit_conf_id, fit_opts, overlay_opts)?;
    let (scores, xform) = align_shape(ref_shape, &mut fit_shape, overlay_opts)?;
    transform_conformer(fit.conformer_mut(fit_conf_id)?, &xform);
    Ok((scores, xform))
}

/// Aligns the conformer `fit_conf_id` of `fit` onto the conformer
/// `ref_conf_id` of `reference`, moving the fit conformer.
#[allow(clippy::too_many_arguments)]
pub fn align_molecule(
    reference: &ROMol,
    fit: &mut ROMol,
    ref_opts: &ShapeInputOptions,
    fit_opts: &ShapeInputOptions,
    overlay_opts: &ShapeOverlayOptions,
    ref_conf_id: i32,
    fit_conf_id: i32,
) -> Result<([f64; 3], Transform3D)> {
    let ref_shape = ShapeInput::new(reference, ref_conf_id, ref_opts, overlay_opts)?;
    align_molecule_to_shape(&ref_shape, fit, fit_opts, overlay_opts, fit_conf_id)
}

/// The result of `score_molecule_all_conformers`.
#[derive(Debug, Clone)]
pub struct AllConformerScores {
    /// The (shape) number of the reference in the best scoring overlay.
    pub ref_conf_id: i32,
    /// The (shape) number of the fit in the best scoring overlay.
    pub fit_conf_id: i32,
    /// `comb_scores[i][j]` is the score of the jth fit conformer onto the ith
    /// reference conformer.
    pub comb_scores: Vec<Vec<f64>>,
    /// The transformation of the best scoring overlay.
    pub xform: Transform3D,
}

/// Overlays all the conformers of `fit` on all those of `reference`.
pub fn score_molecule_all_conformers(
    reference: &ROMol,
    fit: &ROMol,
    ref_opts: &ShapeInputOptions,
    fit_opts: &ShapeInputOptions,
    overlay_opts: &ShapeOverlayOptions,
) -> Result<AllConformerScores> {
    // Pruning the shapes wastes time and obviously removes the correspondence
    // between conformers and shapes.
    let mut ref_cp = ref_opts.clone();
    ref_cp.shape_prune_threshold = -1.0;
    ref_cp.sort_shapes = false;
    let mut fit_cp = fit_opts.clone();
    fit_cp.shape_prune_threshold = -1.0;
    fit_cp.sort_shapes = false;
    let ref_shape = ShapeInput::new(reference, -1, &ref_cp, overlay_opts)?;
    let fit_shape = ShapeInput::new(fit, -1, &fit_cp, overlay_opts)?;
    let (n_ref, n_fit) = (ref_shape.num_shapes(), fit_shape.num_shapes());
    // The overlays are independent, so they run on several threads; the best
    // is then picked in the serial order, so the result does not depend on the
    // number of threads.
    let pairs: Vec<(usize, usize)> = (0..n_ref)
        .flat_map(|i| (0..n_fit).map(move |j| (i, j)))
        .collect();
    let results = par_map(&pairs, 0, |&(i, j)| -> Result<([f64; 3], Transform3D)> {
        let mut r = ref_shape.clone();
        let mut f = fit_shape.clone();
        r.set_active_shape(i)?;
        f.set_active_shape(j)?;
        align_shape(&r, &mut f, overlay_opts)
    });
    let mut comb_scores = vec![vec![0.0; n_fit]; n_ref];
    let mut best_score = -1.0;
    let mut res = AllConformerScores {
        ref_conf_id: 0,
        fit_conf_id: 0,
        comb_scores: Vec::new(),
        xform: Transform3D::identity(),
    };
    for (&(i, j), result) in pairs.iter().zip(results) {
        let (scores, this_xform) = result?;
        comb_scores[i][j] = scores[0];
        if scores[0] > best_score {
            best_score = scores[0];
            res.ref_conf_id = i as i32;
            res.fit_conf_id = j as i32;
            res.xform = this_xform;
        }
    }
    res.comb_scores = comb_scores;
    Ok(res)
}

/// Scores the shapes as they are, without moving them. The scores are the
/// combination, shape Tversky and colour Tversky; the volumes are those of
/// the shape and colour overlaps.
pub fn score_shape(
    ref_shape: &ShapeInput,
    fit_shape: &ShapeInput,
    opts: &ShapeOverlayOptions,
) -> ([f64; 3], (f64, f64)) {
    let mut sca = SingleConformerAlignment::new(
        ref_shape.view(),
        fit_shape.view(),
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        opts,
    );
    let scores = sca.calc_scores_for(ref_shape.coords(), fit_shape.coords(), use_color(opts));
    ([scores[0], scores[1], scores[2]], (scores[3], scores[4]))
}

/// Scores the conformer `fit_conf_id` of `fit` against the reference shape.
pub fn score_molecule_to_shape(
    ref_shape: &ShapeInput,
    fit: &ROMol,
    fit_opts: &ShapeInputOptions,
    overlay_opts: &ShapeOverlayOptions,
    fit_conf_id: i32,
) -> Result<([f64; 3], (f64, f64))> {
    let fit_shape = ShapeInput::new(fit, fit_conf_id, fit_opts, overlay_opts)?;
    Ok(score_shape(ref_shape, &fit_shape, overlay_opts))
}

/// Scores the conformers of 2 molecules as they are placed.
#[allow(clippy::too_many_arguments)]
pub fn score_molecule(
    reference: &ROMol,
    fit: &ROMol,
    ref_opts: &ShapeInputOptions,
    fit_opts: &ShapeInputOptions,
    overlay_opts: &ShapeOverlayOptions,
    ref_conf_id: i32,
    fit_conf_id: i32,
) -> Result<([f64; 3], (f64, f64))> {
    let mut tmp_opts = overlay_opts.clone();
    tmp_opts.normalize = false;
    tmp_opts.start_mode = StartMode::Rotate0;
    let ref_shape = ShapeInput::new(reference, ref_conf_id, ref_opts, &tmp_opts)?;
    let fit_shape = ShapeInput::new(fit, fit_conf_id, fit_opts, &tmp_opts)?;
    Ok(score_shape(&ref_shape, &fit_shape, &tmp_opts))
}
