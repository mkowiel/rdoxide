//! `SingleConformerAlignment`: the overlay of one conformer of a fit shape
//! onto one conformer of a reference shape, optimising a rotation (as a
//! quaternion) and translation with the step rules of the PubChem code.

use std::f64::consts::PI;

use super::options::{OptimMode, ShapeOverlayOptions};
use super::shape_input::{KAPPA, apply_transform_to_shape_into, translate_shape_into};
use crate::geometry::{Point3D, Transform3D};

/// The pieces of a shape that the alignment needs: coordinates for the atoms
/// followed by the colour features, `alphas` and `types` with an entry for
/// each of them.
#[derive(Debug, Clone, Copy)]
pub struct ShapeView<'a> {
    pub coords: &'a [f64],
    pub alphas: &'a [f64],
    pub types: &'a [i32],
    pub carbon_radii: Option<&'a [bool]>,
    pub num_shape: usize,
    pub num_color: usize,
    pub shape_vol: f64,
    pub color_vol: f64,
}

/// Scores are `[combo, shape Tversky, colour Tversky, shape overlap volume,
/// colour overlap volume]`.
pub type Scores = [f64; 5];

/// The quaternion (first 4) and translation (last 3) of a transformation.
pub type QuatTrans = [f64; 7];

#[derive(Debug)]
pub struct SingleConformerAlignment<'a> {
    ref_coords: Vec<f64>,
    ref_temp: Vec<f64>,
    refv: ShapeView<'a>,
    fit_coords: Vec<f64>,
    fit_temp: Vec<f64>,
    fitv: ShapeView<'a>,
    init_quat_trans: QuatTrans,
    optim_mode: OptimMode,
    sim_alpha: f64,
    sim_beta: f64,
    mixing_param: f64,
    use_cutoff: bool,
    dist_cutoff2: f64,
    shape_convergence_criterion: f64,
    max_its: u32,
    /// the quaternion/translation as the optimisation proceeds
    quat_trans: QuatTrans,
    /// The step sizes of the quaternion and translation during the
    /// optimisation, from the PubChem code.
    q_step_size: f64,
    t_step_size: f64,
    /// Scratch space for the gradients dr/dQ of the fit molecule.
    grad_converters: Vec<f64>,
}

/// `cartToQuatGrads`
fn cart_to_quat_grads(quat: &[f64], mol: &[f64], num_pts: usize, conv: &mut [f64], offset: usize) {
    let (q, r, s, u) = (quat[0], quat[1], quat[2], quat[3]);
    let coef = 1.0 / (q * q + r * r + s * s + u * u);
    for i in 0..num_pts {
        let (x, y, z) = (mol[3 * i], mol[3 * i + 1], mol[3 * i + 2]);
        let dx_dq = coef * 2.0 * (q * x + u * y - s * z);
        let dx_dr = coef * 2.0 * (r * x + s * y + u * z);
        let dy_dr = coef * 2.0 * (s * x - r * y + q * z);
        let dx_du = coef * 2.0 * (-u * x + q * y + r * z);
        let dz_ds = dx_dq;
        let dy_du = -dx_dq;
        let dy_ds = dx_dr;
        let dz_du = dx_dr;
        let dx_ds = -dy_dr;
        let dz_dq = dy_dr;
        let dy_dq = dx_du;
        let dz_dr = -dx_du;
        let k = 12 * (offset + i);
        conv[k..k + 12].copy_from_slice(&[
            dx_dq, dy_dq, dz_dq, dx_dr, dy_dr, dz_dr, dx_ds, dy_ds, dz_ds, dx_du, dy_du, dz_du,
        ]);
    }
}

/// Adds the gradient of one pair term to `gradients`.
fn add_grads(gradients: &mut [f64; 7], conv: &[f64], k: usize, r: f64, d: [f64; 3]) {
    let [dx, dy, dz] = d;
    gradients[1] += r * (dx * conv[k + 3] + dy * conv[k + 4] + dz * conv[k + 5]);
    gradients[2] += r * (dx * conv[k + 6] + dy * conv[k + 7] + dz * conv[k + 8]);
    gradients[3] += r * (dx * conv[k + 9] + dy * conv[k + 10] + dz * conv[k + 11]);
    gradients[4] += r * dx;
    gradients[5] += r * dy;
    gradients[6] += r * dz;
}

/// The gradient request of the volume functions: the quaternion/translation
/// being evaluated and the array the gradients are added to.
pub type GradRequest<'g> = Option<(&'g QuatTrans, &'g mut [f64; 7])>;

/// `calcVolAndGrads` for the atoms: the overlap volume of two sets of
/// Gaussians (negative alphas are switched off), optionally with the
/// gradients with respect to the quaternion and translation.
#[allow(clippy::too_many_arguments)]
pub fn calc_shape_vol_and_grads(
    refc: &[f64],
    ref_alphas: &[f64],
    num_ref: usize,
    ref_carbon: Option<&[bool]>,
    fit: &[f64],
    fit_alphas: &[f64],
    num_fit: usize,
    fit_carbon: Option<&[bool]>,
    grad_converters: &mut [f64],
    use_cutoff: bool,
    dist_cutoff2: f64,
    mut grads: GradRequest<'_>,
) -> f64 {
    if let Some((quat, _)) = &grads {
        cart_to_quat_grads(&quat[..4], fit, num_fit, grad_converters, 0);
    }
    let carbon_a = KAPPA / (1.7 * 1.7);
    let carbon_bit = 8.0 * (PI / (2.0 * carbon_a)).powf(1.5);
    let mut vol = 0.0;
    // If either of the carbon radii flags isn't supplied, treat them both as
    // being all carbon: there isn't enough information to do otherwise.
    let carbons = ref_carbon.zip(fit_carbon);
    for i in 0..num_ref {
        let ai = ref_alphas[i];
        if ai < 0.0 {
            continue;
        }
        for j in 0..num_fit {
            let dx = refc[3 * i] - fit[3 * j];
            let dy = refc[3 * i + 1] - fit[3 * j + 1];
            let dz = refc[3 * i + 2] - fit[3 * j + 2];
            let d2 = dx * dx + dy * dy + dz * dz;
            if use_cutoff && d2 > dist_cutoff2 {
                continue;
            }
            let aj = fit_alphas[j];
            if aj < 0.0 {
                continue;
            }
            let mult = -(ai * aj) / (ai + aj);
            let kij = (mult * d2).exp();
            let vij = match carbons {
                Some((rc, fc)) if !(rc[i] && fc[j]) => {
                    let pi_ai_aj = PI / (ai + aj);
                    8.0 * kij * pi_ai_aj * pi_ai_aj.sqrt()
                }
                _ => kij * carbon_bit,
            };
            vol += vij;
            if let Some((_, g)) = grads.as_mut() {
                add_grads(g, grad_converters, 12 * j, 2.0 * vij * mult, [dx, dy, dz]);
            }
        }
    }
    vol
}

/// `calcVolAndGrads` for the colour features: only features of the same type
/// overlap.
///
/// Like upstream, the gradient converters are read from the start of the
/// scratch space, i.e. those of the first fit atoms rather than the ones
/// computed for the features (a quirk that is kept for compatibility).
#[allow(clippy::too_many_arguments)]
pub fn calc_color_vol_and_grads(
    refc: &[f64],
    ref_alphas: &[f64],
    num_ref: usize,
    ref_types: &[i32],
    fit: &[f64],
    fit_alphas: &[f64],
    num_fit: usize,
    fit_types: &[i32],
    num_fit_shape: usize,
    grad_converters: &mut [f64],
    use_cutoff: bool,
    dist_cutoff2: f64,
    mut grads: GradRequest<'_>,
) -> f64 {
    if let Some((quat, _)) = &grads {
        cart_to_quat_grads(&quat[..4], fit, num_fit, grad_converters, num_fit_shape);
    }
    let mut vol = 0.0;
    for i in 0..num_ref {
        let ai = ref_alphas[i];
        for j in 0..num_fit {
            if ref_types[i] != fit_types[j] {
                continue;
            }
            let dx = refc[3 * i] - fit[3 * j];
            let dy = refc[3 * i + 1] - fit[3 * j + 1];
            let dz = refc[3 * i + 2] - fit[3 * j + 2];
            let d2 = dx * dx + dy * dy + dz * dz;
            if use_cutoff && d2 > dist_cutoff2 {
                continue;
            }
            let aj = fit_alphas[j];
            let mult = -(ai * aj) / (ai + aj);
            let kij = (mult * d2).exp();
            let pi_ai_aj = PI / (ai + aj);
            let vij = 8.0 * kij * pi_ai_aj * pi_ai_aj.sqrt();
            vol += vij;
            if let Some((_, g)) = grads.as_mut() {
                add_grads(g, grad_converters, 12 * j, 2.0 * vij * mult, [dx, dy, dz]);
            }
        }
    }
    vol
}

fn one_step(grad: f64, step_size: f64, quat_trans: f64, old_grad: f64, old_quat_trans: f64) -> f64 {
    if grad.is_sign_negative() != old_grad.is_sign_negative() {
        let mut step = (quat_trans * old_grad.abs() + old_quat_trans * grad.abs())
            / (old_grad.abs() + grad.abs() + grad.abs())
            - quat_trans;
        let new_step = step_size * grad;
        if step.abs() > new_step.abs() {
            // This is definitely what the PubChem code says! It reads as
            // keeping the sign of step, but the value of newStep.
            step *= (new_step / step).abs();
        }
        step
    } else {
        step_size * grad
    }
}

fn calc_step(
    grad: &QuatTrans,
    q_step: f64,
    t_step: f64,
    old_grad: &QuatTrans,
    quat_trans: &QuatTrans,
    old_quat_trans: &QuatTrans,
    iter: u32,
) -> QuatTrans {
    let mut step = [0.0; 7];
    for i in 1..7 {
        let size = if i < 4 { q_step } else { t_step };
        step[i] = if iter == 0 {
            // 1st iteration, use default step sizes
            size * grad[i]
        } else {
            one_step(grad[i], size, quat_trans[i], old_grad[i], old_quat_trans[i])
        };
    }
    step
}

/// `constrainStep` for the 3 components in `step`; returns the largest
/// component before constraining.
fn constrain_step(max_step: f64, step: &mut [f64], check_size: bool) -> f64 {
    let m_step = step[0].abs().max(step[1].abs()).max(step[2].abs());
    if m_step > max_step {
        let scale = max_step / m_step;
        for s in step[..3].iter_mut() {
            if s.abs() > max_step {
                *s *= scale;
            }
        }
    }
    if check_size {
        let quat_squared = step[0] * step[0] + step[1] * step[1] + step[2] * step[2];
        if quat_squared > 1.0 {
            let scale = 1.0 / (2.0 * quat_squared);
            for s in step[..3].iter_mut() {
                *s *= scale;
            }
        }
    }
    m_step
}

fn combine_quat_trans(q1: &QuatTrans, q2: &QuatTrans) -> QuatTrans {
    // Multiply the quaternions, which are assumed to be normalised, and add
    // the translations.
    [
        q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3],
        q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2],
        q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1],
        q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0],
        q1[4] + q2[4],
        q1[5] + q2[5],
        q1[6] + q2[6],
    ]
}

fn one_reduce_step(
    grad: f64,
    old_grad: f64,
    quat_trans: f64,
    old_quat_trans: f64,
    step_size: f64,
    mut step: f64,
) -> f64 {
    if grad.is_sign_negative() != old_grad.is_sign_negative() {
        step = (quat_trans * old_grad.abs() + old_quat_trans * grad.abs())
            / (old_grad.abs() + grad.abs())
            - quat_trans;
        let new_step = step_size * grad;
        if step.abs() > new_step.abs() {
            step *= (new_step / step).abs();
        }
    } else if grad.abs() <= 1.0 {
        step = step_size * grad;
    } else if grad.abs() > old_grad.abs() {
        // Going wrong way relative to other components?
        step += step_size * grad;
    } else {
        let mut delta = grad * (step / (old_grad - grad));
        if delta.abs() > (step * 0.1).abs() && delta.abs() > 0.001 {
            delta *= 0.0005 / delta.abs();
        }
        step += delta;
    }
    step
}

#[allow(clippy::too_many_arguments)]
fn reduce_step(
    grad: &QuatTrans,
    old_grad: &QuatTrans,
    quat_trans: &QuatTrans,
    old_quat_trans: &QuatTrans,
    line_iter: u32,
    step: &mut QuatTrans,
    q_step_size: &mut f64,
    t_step_size: &mut f64,
) {
    if line_iter == 2 {
        *q_step_size *= 0.1;
        *t_step_size *= 0.1;
        for i in 1..7 {
            let size = if i < 4 { *q_step_size } else { *t_step_size };
            step[i] = size * old_grad[i];
        }
        *q_step_size *= 5.0;
        *t_step_size *= 5.0;
    } else {
        for i in 1..7 {
            // The original PubChem code used the quaternion step size for the
            // translations too, which looks like a cut-and-paste error.
            let size = if i < 4 { *q_step_size } else { *t_step_size };
            step[i] = one_reduce_step(
                grad[i],
                old_grad[i],
                quat_trans[i],
                old_quat_trans[i],
                size,
                step[i],
            );
        }
    }
}

impl<'a> SingleConformerAlignment<'a> {
    /// Sets up the overlay of `fit` onto `reference`; `init_quat_trans` is the
    /// starting rotation and translation.
    pub fn new(
        reference: ShapeView<'a>,
        fit: ShapeView<'a>,
        init_quat_trans: QuatTrans,
        opts: &ShapeOverlayOptions,
    ) -> Self {
        // Move the reference by the initial translation, leaving fit at the
        // origin where the rotations work properly. Apply the initial rotation
        // to the fit.
        let mut ref_coords = reference.coords.to_vec();
        translate_shape_into_place(
            &mut ref_coords,
            Point3D::new(init_quat_trans[4], init_quat_trans[5], init_quat_trans[6]),
        );
        let mut xform = Transform3D::identity();
        xform.set_rotation_from_quaternion(&quat_of(&init_quat_trans));
        let mut fit_coords = fit.coords.to_vec();
        super::shape_input::apply_transform_to_shape(&mut fit_coords, &xform);
        SingleConformerAlignment {
            ref_temp: vec![0.0; ref_coords.len()],
            fit_temp: vec![0.0; fit_coords.len()],
            ref_coords,
            refv: reference,
            fit_coords,
            fitv: fit,
            init_quat_trans,
            optim_mode: opts.optim_mode,
            sim_alpha: opts.sim_alpha,
            sim_beta: opts.sim_beta,
            mixing_param: opts.opt_param,
            use_cutoff: opts.use_dist_cutoff,
            dist_cutoff2: opts.dist_cutoff * opts.dist_cutoff,
            shape_convergence_criterion: opts.shape_convergence_criterion,
            max_its: opts.n_steps,
            quat_trans: [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            q_step_size: -0.001,
            t_step_size: -0.01,
            grad_converters: vec![0.0; 12 * (fit.num_shape + fit.num_color)],
        }
    }

    pub fn set_quat_trans(&mut self, quat_trans: QuatTrans) {
        self.quat_trans = quat_trans;
    }

    /// The final transformation: the initial transformation combined with
    /// the optimised one.
    pub fn final_transform(&self) -> Transform3D {
        let mut tmp = Transform3D::identity();
        tmp.set_rotation_from_quaternion(&quat_of(&self.quat_trans));
        tmp.set_translation(&Point3D::new(
            self.quat_trans[4],
            self.quat_trans[5],
            self.quat_trans[6],
        ));
        let mut reverse_initial_trans = Transform3D::identity();
        reverse_initial_trans.set_translation(&Point3D::new(
            -self.init_quat_trans[4],
            -self.init_quat_trans[5],
            -self.init_quat_trans[6],
        ));
        let mut initial_rot = Transform3D::identity();
        initial_rot.set_rotation_from_quaternion(&quat_of(&self.init_quat_trans));
        reverse_initial_trans.compose(&tmp).compose(&initial_rot)
    }

    /// The scores of the shapes at `refc` and `fit` (assumed to be already in
    /// the correct configurations). The colour score is computed if
    /// `include_color` is set whatever the optimisation mode, as it is still
    /// wanted for `ShapeOnly` optimisations.
    pub fn calc_scores_for(&mut self, refc: &[f64], fit: &[f64], include_color: bool) -> Scores {
        let (r, f) = (self.refv, self.fitv);
        let shape_vol = calc_shape_vol_and_grads(
            refc,
            r.alphas,
            r.num_shape,
            r.carbon_radii,
            fit,
            f.alphas,
            f.num_shape,
            f.carbon_radii,
            &mut self.grad_converters,
            self.use_cutoff,
            self.dist_cutoff2,
            None,
        );
        let mut color_vol = 0.0;
        if r.num_color > 0
            && f.num_color > 0
            && (self.optim_mode == OptimMode::ShapePlusColor || include_color)
        {
            color_vol = calc_color_vol_and_grads(
                &refc[3 * r.num_shape..],
                &r.alphas[r.num_shape..],
                r.num_color,
                &r.types[r.num_shape..],
                &fit[3 * f.num_shape..],
                &f.alphas[f.num_shape..],
                f.num_color,
                &f.types[f.num_shape..],
                f.num_shape,
                &mut self.grad_converters,
                self.use_cutoff,
                self.dist_cutoff2,
                None,
            );
        }
        self.scores_from_volumes(shape_vol, color_vol, include_color)
    }

    /// The scores after applying the current quaternion and translation.
    pub fn calc_scores(&mut self, include_color: bool) -> Scores {
        self.apply_quat_trans(&self.quat_trans.clone());
        let (refc, fit) = (
            std::mem::take(&mut self.ref_temp),
            std::mem::take(&mut self.fit_temp),
        );
        let scores = self.calc_scores_for(&refc, &fit, include_color);
        self.ref_temp = refc;
        self.fit_temp = fit;
        scores
    }

    /// The scores from the given overlap volumes. The colour score is only
    /// calculated if the colour volumes are non-zero.
    pub fn scores_from_volumes(
        &self,
        shape_ov_vol: f64,
        color_ov_vol: f64,
        include_color: bool,
    ) -> Scores {
        let (r, f) = (&self.refv, &self.fitv);
        let mut scores = [0.0; 5];
        scores[3] = shape_ov_vol;
        scores[4] = color_ov_vol;
        scores[1] = shape_ov_vol
            / (self.sim_alpha * (r.shape_vol - shape_ov_vol)
                + self.sim_beta * (f.shape_vol - shape_ov_vol)
                + shape_ov_vol);
        if r.num_color > 0
            && f.num_color > 0
            && r.color_vol > 0.0
            && f.color_vol > 0.0
            && include_color
        {
            scores[2] = color_ov_vol
                / (self.sim_alpha * (r.color_vol - color_ov_vol)
                    + self.sim_beta * (f.color_vol - color_ov_vol)
                    + color_ov_vol);
            scores[0] = scores[1] * (1.0 - self.mixing_param) + scores[2] * self.mixing_param;
        } else {
            scores[0] = scores[1];
        }
        scores
    }

    /// Puts the reference translated by the negative of the translation, and
    /// the fit rotated by the quaternion, into the temporary coordinates.
    pub fn apply_quat_trans(&mut self, quat_trans: &QuatTrans) {
        let n_ref = self.refv.num_shape + self.refv.num_color;
        translate_shape_into(
            &self.ref_coords,
            &mut self.ref_temp,
            n_ref,
            Point3D::new(-quat_trans[4], -quat_trans[5], -quat_trans[6]),
        );
        let mut transform = Transform3D::identity();
        transform.set_rotation_from_quaternion(&quat_of(quat_trans));
        let n_fit = self.fitv.num_shape + self.fitv.num_color;
        apply_transform_to_shape_into(&self.fit_coords, &mut self.fit_temp, n_fit, &transform);
    }

    /// The overlap volumes (shape and colour) and their gradients with
    /// respect to the quaternion and translation.
    pub fn calc_volume_and_gradients(&mut self, quat_trans: &QuatTrans) -> (f64, f64, QuatTrans) {
        self.apply_quat_trans(quat_trans);
        let (r, f) = (self.refv, self.fitv);
        let mut gradients = [0.0; 7];
        let shape_vol = calc_shape_vol_and_grads(
            &self.ref_temp,
            r.alphas,
            r.num_shape,
            r.carbon_radii,
            &self.fit_temp,
            f.alphas,
            f.num_shape,
            f.carbon_radii,
            &mut self.grad_converters,
            self.use_cutoff,
            self.dist_cutoff2,
            Some((quat_trans, &mut gradients)),
        );
        if self.optim_mode != OptimMode::ShapePlusColor {
            return (shape_vol, 0.0, gradients);
        }
        let mut color_grads = [0.0; 7];
        let color_vol = calc_color_vol_and_grads(
            &self.ref_temp[3 * r.num_shape..],
            &r.alphas[r.num_shape..],
            r.num_color,
            &r.types[r.num_shape..],
            &self.fit_temp[3 * f.num_shape..],
            &f.alphas[f.num_shape..],
            f.num_color,
            &f.types[f.num_shape..],
            f.num_shape,
            &mut self.grad_converters,
            self.use_cutoff,
            self.dist_cutoff2,
            Some((quat_trans, &mut color_grads)),
        );
        // The colour gradients are normally dwarfed by the shape gradients, so
        // normalize them and then mix by the same rule as the final score.
        let norm = |g: &QuatTrans| g[1..].iter().fold(0.0, |acc, v| acc + v * v).sqrt();
        let ratio = norm(&gradients) / norm(&color_grads);
        for i in 1..7 {
            gradients[i] = gradients[i] * (1.0 - self.mixing_param)
                + color_grads[i] * ratio * self.mixing_param;
        }
        (shape_vol, color_vol, gradients)
    }

    /// Runs the optimisation (a short first cycle, a longer second one) and
    /// returns the scores and the transformation found, in the layout of
    /// upstream's `scores` array: combo, shape and colour Tversky, shape and
    /// colour overlap volumes, the fit and ref volumes, quaternion, and
    /// translation. Also whether it finished within the iteration limit.
    pub fn do_overlay(&mut self, cycle: u32) -> (Scores, bool) {
        let max_iters = if cycle == 0 {
            10
        } else {
            self.max_its.saturating_sub(10)
        };
        let finished = self.optimise(max_iters);

        // Get the final coords for fit into fit_temp, and compute the scores
        let mut xform = Transform3D::identity();
        xform.set_rotation_from_quaternion(&quat_of(&self.quat_trans));
        xform.set_translation(&Point3D::new(
            self.quat_trans[4],
            self.quat_trans[5],
            self.quat_trans[6],
        ));
        let n_fit = self.fitv.num_shape + self.fitv.num_color;
        apply_transform_to_shape_into(&self.fit_coords, &mut self.fit_temp, n_fit, &xform);
        let fit = std::mem::take(&mut self.fit_temp);
        let refc = std::mem::take(&mut self.ref_coords);
        let scores = self.calc_scores_for(&refc, &fit, true);
        self.fit_temp = fit;
        self.ref_coords = refc;
        (scores, finished)
    }

    /// Finds the quaternion and translation that maximise the overlap
    /// appropriate to the optimisation mode. Returns false if it didn't finish
    /// within the allowed number of iterations.
    pub fn optimise(&mut self, max_iters: u32) -> bool {
        const MAX_QUATERNION_STEP: f64 = 0.075;
        const MAX_TRANSLATION_STEP: f64 = 0.500;
        const MIN_QUATERNION_STEP: f64 = 0.0002;
        const MIN_TRANSLATION_STEP: f64 = 0.0020;

        let mut old_grad = [0.0; 7];
        let mut old_quat_trans = [0.0; 7];
        let mut finished = false;
        for iter in 0..max_iters {
            let (shape_vol, color_vol, grad) =
                self.calc_volume_and_gradients(&self.quat_trans.clone());
            let mut grad = grad;
            // Note that the combo score will have a zero colour score so will
            // be half the shape score unless we're optimising with colour
            // gradients.
            let mut combo_score = self.scores_from_volumes(shape_vol, color_vol, true)[0];
            let mut step = calc_step(
                &grad,
                self.q_step_size,
                self.t_step_size,
                &old_grad,
                &self.quat_trans,
                &old_quat_trans,
                iter,
            );

            // In case we have to backtrack
            let old_combo_score = combo_score;
            old_quat_trans = self.quat_trans;
            old_grad = grad;

            // What the PubChem code calls "Line search (sort of) loop"
            let mut converged = false;
            let mut line_iter = 0;
            while !converged {
                // Check that the absolute max step size does not go beyond
                // some reasonable size
                let mq_step = constrain_step(MAX_QUATERNION_STEP, &mut step[1..4], true);
                let mt_step = constrain_step(MAX_TRANSLATION_STEP, &mut step[4..7], false);
                if mq_step <= MIN_QUATERNION_STEP && mt_step <= MIN_TRANSLATION_STEP {
                    converged = true;
                    combo_score = 0.0; // Make sure we return to the old one.
                    break;
                }
                // Calculate the 0th component of the quaternion. Obviously it
                // relies on the other 3 components being small.
                let quat_squared = step[1] * step[1] + step[2] * step[2] + step[3] * step[3];
                step[0] = (1.0 - quat_squared).sqrt();
                // Update the quaternion with the step, multiplying them.
                let new_quat_trans = combine_quat_trans(&self.quat_trans, &step);
                let (sv, cv, new_grad) = self.calc_volume_and_gradients(&new_quat_trans);
                grad = new_grad;
                combo_score = self.scores_from_volumes(sv, cv, true)[0];
                // if we made a good step, keep the quaternion and we're done
                if combo_score > old_combo_score {
                    self.quat_trans = new_quat_trans;
                    break;
                }
                if line_iter > 2 {
                    converged = true;
                    self.quat_trans = new_quat_trans;
                    break;
                }
                // It got worse, so reduce the step.
                reduce_step(
                    &grad,
                    &old_grad,
                    &new_quat_trans,
                    &old_quat_trans,
                    line_iter,
                    &mut step,
                    &mut self.q_step_size,
                    &mut self.t_step_size,
                );
                self.quat_trans = old_quat_trans;
                line_iter += 1;
            }
            // Did it converge?
            if converged || self.shape_convergence_criterion > combo_score - old_combo_score {
                if old_combo_score > combo_score {
                    // The previous step was better, so keep it.
                    self.quat_trans = old_quat_trans;
                }
                finished = true;
                break;
            }
        }
        finished
    }
}

fn quat_of(q: &QuatTrans) -> [f64; 4] {
    [q[0], q[1], q[2], q[3]]
}

fn translate_shape_into_place(shape: &mut [f64], t: Point3D) {
    for p in shape.chunks_exact_mut(3) {
        p[0] += t.x;
        p[1] += t.y;
        p[2] += t.z;
    }
}
