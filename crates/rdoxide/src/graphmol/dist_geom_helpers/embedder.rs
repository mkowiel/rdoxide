//! Embedding molecules in 3D with distance geometry and ETKDG (upstream
//! `Embedder.cpp`).
//!
//! Differences from upstream: a negative random seed draws a seed from the
//! global random number generator for each conformer instead of sharing that
//! generator between threads; there is no Ctrl-C handling; the failure
//! counters are collected in thread-safe counters and copied to the
//! parameters at the end.

use std::sync::Mutex;
use std::time::{Duration, Instant};

use super::bounds_builder::{TopolBoundsOptions, TopolBoundsOutput, set_topol_bounds};
use super::bounds_builder::{collect_bonds_and_angles, init_bounds_mat};
use super::embed_checks::{
    AMPLITUDE_FOURTH_COORDINATE_IC, ERROR_TOL, MAX_MINIMIZED_E_PER_ATOM,
    MAX_MINIMIZED_E_PER_ATOM_AIO, TETRAHEDRAL_CENTERINVOLUME_TOL, bounds_fulfilled,
    center_in_volume, check_chiral_centers, check_tetrahedral_centers, clash_check,
    double_bond_geometry_checks, point,
};
use super::embed_params::{
    EmbedFailureCause, EmbedParameters, InitialEmbeddingMode, NUM_EMBED_FAILURE_CAUSES,
};
use super::embed_setup::{
    CoordMap, DoubleBondEnd, StereoDoubleBond, find_chiral_sets, find_double_bonds, init_etkdg,
    setup_initial_bounds_matrix,
};
use super::internal_coordinates::InternalCoordinates;
use super::z_matrix_builder::{correct_chiral_centers, set_molecule_dfs};
use crate::distgeom::{
    BoundsMatrix, ChiralSet, ZMatrix, add_torsion_terms, compute_initial_coords,
    compute_random_coords, compute_z_matrix_coords, construct_3d_force_field,
    construct_3d_improper_force_field, construct_all_in_one_force_field, construct_force_field,
    construct_plain_3d_force_field, improper_term_energies, pick_random_dist_mat,
};
use crate::error::{RDKitError, Result};
use crate::forcefield::ForceField;
use crate::geometry::Point3D;
use crate::geometry::point::compute_dihedral_angle;
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::forcefield_helpers::crystal_ff::CrystalFFDetails;
use crate::graphmol::matrices::get_distance_mat;
use crate::graphmol::mol_align::symmetrize_terminal_atoms;
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, needs_hs, remove_hs};
use crate::graphmol::mol_ops::frags::{get_mol_frags_as_mols, get_mol_frags_mapping};
use crate::graphmol::romol::ROMol;
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::numerics::alignment::align_points;
use crate::numerics::symm_matrix::SymmMatrix;
use crate::rd_log_warning;
use crate::rdgeneral::random::{MinStdRand, with_global_generator};
use crate::rdgeneral::threads::{num_threads_to_use, par_map};

/// everything about the molecule (fragment) that the embedding needs
struct EmbedArgs {
    four_d: bool,
    mmat: BoundsMatrix,
    zmat: ZMatrix,
    chiral_centers: Vec<ChiralSet>,
    tetrahedral_carbons: Vec<ChiralSet>,
    double_bond_ends: Vec<DoubleBondEnd>,
    stereo_double_bonds: Vec<StereoDoubleBond>,
    details: CrystalFFDetails,
    /// number of heavy atoms
    hac: usize,
}

impl EmbedArgs {
    fn dim(&self) -> usize {
        if self.four_d { 4 } else { 3 }
    }
}

/// Counts of the failure causes, shared between threads.
struct Failures {
    counts: Option<Mutex<Vec<u32>>>,
}

impl Failures {
    fn new(track: bool) -> Self {
        Failures {
            counts: track.then(|| Mutex::new(vec![0; NUM_EMBED_FAILURE_CAUSES])),
        }
    }

    fn record(&self, cause: EmbedFailureCause) {
        if let Some(counts) = &self.counts {
            counts.lock().unwrap_or_else(|e| e.into_inner())[cause as usize] += 1;
        }
    }
}

fn timed_out(end_time: Option<Instant>) -> bool {
    end_time.is_some_and(|t| Instant::now() > t)
}

fn flat_to_points(positions: &[f64], dim: usize) -> Vec<Point3D> {
    positions
        .chunks_exact(dim)
        .map(|c| Point3D::new(c[0], c[1], c[2]))
        .collect()
}

/// Copies the positions of a force field back into the flat array.
fn copy_back(field: &ForceField, positions: &mut [f64]) {
    positions.copy_from_slice(field.positions());
}

fn fixed_points_from(params: &EmbedParameters) -> Vec<usize> {
    match (&params.coord_map, params.use_random_coords) {
        (Some(map), true) => map.keys().map(|&k| k as usize).collect(),
        _ => Vec::new(),
    }
}

/// Distance matrix (and random) coordinates for a start of the embedding.
fn generate_initial_coords(
    positions: &mut [f64],
    eargs: &EmbedArgs,
    params: &EmbedParameters,
    dist_mat: &mut SymmMatrix,
    rng: &mut MinStdRand,
) -> Result<bool> {
    let dim = eargs.dim();
    match params.initial_embedding_mode {
        InitialEmbeddingMode::DgEmbedding => {
            pick_random_dist_mat(&eargs.mmat, dist_mat, rng)?;
            compute_initial_coords(
                dist_mat,
                positions,
                dim,
                rng,
                params.rand_neg_eig,
                params.num_zero_fail as usize,
            )
        }
        InitialEmbeddingMode::InternalCoordinateEmbedding => {
            let mut points = flat_to_points(positions, dim);
            compute_z_matrix_coords(&eargs.zmat, &mut points, rng)?;
            for (chunk, p) in positions.chunks_exact_mut(dim).zip(&points) {
                chunk[..3].copy_from_slice(&[p.x, p.y, p.z]);
            }
            if eargs.four_d
                && (!check_tetrahedral_centers(
                    &eargs.tetrahedral_carbons,
                    positions,
                    dim,
                    TETRAHEDRAL_CENTERINVOLUME_TOL,
                ) || !check_chiral_centers(&eargs.chiral_centers, positions, dim))
            {
                // only use a non-zero fourth dimension if at least one chiral
                // center is incorrect
                for chunk in positions.chunks_exact_mut(dim) {
                    chunk[3] = AMPLITUDE_FOURTH_COORDINATE_IC * (rng.next_f64() - 0.5);
                }
            }
            Ok(true)
        }
        InitialEmbeddingMode::RandomCoordinateEmbedding => {
            let box_size = if params.box_size_mult > 0.0 {
                5.0 * params.box_size_mult
            } else {
                -params.box_size_mult
            };
            let got = compute_random_coords(positions, box_size, rng)?;
            if let (true, Some(map)) = (params.use_random_coords, &params.coord_map) {
                for (&idx, p) in map {
                    let chunk = &mut positions[dim * idx as usize..dim * (idx as usize + 1)];
                    chunk[..3].copy_from_slice(&[p.x, p.y, p.z]);
                    // zero out any higher dimensional components
                    chunk[3..].fill(0.0);
                }
            }
            Ok(got)
        }
    }
}

/// The first minimization with the distance and chirality terms.
fn first_minimization(
    positions: &mut [f64],
    eargs: &EmbedArgs,
    params: &EmbedParameters,
    basin_thresh: f64,
) -> Result<bool> {
    let dim = eargs.dim();
    let n = positions.len() / dim;
    let fixed = fixed_points_from(params);
    let mut fixed_pts = vec![false; n];
    for &f in &fixed {
        fixed_pts[f] = true;
    }
    let mut field = construct_force_field(
        &eargs.mmat,
        positions,
        dim,
        &eargs.chiral_centers,
        eargs.details.force_consts.distance,
        eargs.details.force_consts.chiral,
        0.1,
        None,
        basin_thresh,
        Some(&fixed_pts),
    )?;
    field.fixed_points_mut().extend(fixed);
    if field.calc_energy() > ERROR_TOL {
        let mut need_more = 1;
        while need_more != 0 {
            need_more = field.minimize(400, params.optimizer_force_tol, 1e-6)?;
        }
    }
    let local_e = field.calc_energy();
    copy_back(&field, positions);
    // check that the energy is not too high (github #971)
    Ok(local_e / (n as f64) < MAX_MINIMIZED_E_PER_ATOM)
}

/// Minimizes again without the chiral constraints but with a higher weight on
/// the fourth dimension.
fn minimize_fourth_dimension(
    positions: &mut [f64],
    eargs: &EmbedArgs,
    params: &EmbedParameters,
    basin_thresh: f64,
    end_time: Option<Instant>,
) -> Result<bool> {
    let dim = eargs.dim();
    let mut field = construct_force_field(
        &eargs.mmat,
        positions,
        dim,
        &eargs.chiral_centers,
        eargs.details.force_consts.distance,
        0.2,
        eargs.details.force_consts.fourth_dim,
        None,
        basin_thresh,
        None,
    )?;
    field.fixed_points_mut().extend(fixed_points_from(params));
    if field.calc_energy() > ERROR_TOL {
        let mut need_more = 1;
        while need_more != 0 {
            if timed_out(end_time) {
                return Ok(false);
            }
            need_more = field.minimize(200, params.optimizer_force_tol, 1e-6)?;
        }
    }
    copy_back(&field, positions);
    Ok(true)
}

/// The minimization with the experimental torsion preferences; false if the
/// planarity of the sp2 centers is not good enough.
fn minimize_with_exp_torsions(
    positions: &mut [f64],
    eargs: &EmbedArgs,
    params: &EmbedParameters,
) -> Result<bool> {
    let dim = eargs.dim();
    let n = positions.len() / dim;
    // to 3D positions
    let mut positions3d = Vec::with_capacity(3 * n);
    for chunk in positions.chunks_exact(dim) {
        positions3d.extend_from_slice(&chunk[..3]);
    }
    let details = &eargs.details;
    let mut field = if params.use_basic_knowledge {
        // ETKDG or KDG
        construct_3d_force_field(&eargs.mmat, &positions3d, details, params.cpci.as_deref())?
    } else {
        // plain ETDG
        construct_plain_3d_force_field(&eargs.mmat, &positions3d, details)?
    };
    let fixed = fixed_points_from(params);
    field.fixed_points_mut().extend(fixed.iter().copied());
    if field.calc_energy() > ERROR_TOL {
        field.minimize(300, params.optimizer_force_tol, 1e-6)?;
    }
    positions3d.copy_from_slice(field.positions());

    let mut planar = true;
    if params.use_basic_knowledge {
        // a force field with only the impropers: is the energy low enough?
        let mut field2 = construct_3d_improper_force_field(
            &eargs.mmat,
            &positions3d,
            &details.improper_atoms,
            &details.angles,
        )?;
        field2.fixed_points_mut().extend(fixed);
        let planarity_tolerance = 0.7;
        let n_centers = details.improper_atoms.len();
        if field2.calc_energy() > n_centers as f64 * planarity_tolerance {
            planar = false;
        }
    }
    // overwrite the positions
    for (chunk, p) in positions
        .chunks_exact_mut(dim)
        .zip(positions3d.chunks_exact(3))
    {
        chunk[..3].copy_from_slice(p);
    }
    Ok(planar)
}

/// The all-in-one minimization: distance terms first, then the torsion terms.
fn minimize_all_in_one(
    positions: &mut [f64],
    eargs: &EmbedArgs,
    params: &EmbedParameters,
    end_time: Option<Instant>,
) -> Result<bool> {
    let dim = eargs.dim();
    let n = positions.len() / dim;
    let fixed = fixed_points_from(params);
    let mut fixed_pts = vec![false; n];
    for &f in &fixed {
        fixed_pts[f] = true;
    }
    let mut field = construct_all_in_one_force_field(
        &eargs.mmat,
        positions,
        dim,
        &eargs.details,
        &eargs.chiral_centers,
        params.cpci.as_deref(),
        None,
        Some(&fixed_pts),
    )?;
    field.fixed_points_mut().extend(fixed);
    if field.calc_energy() > ERROR_TOL {
        field.minimize(eargs.hac as u32, params.optimizer_force_tol, 1e-6)?;
    }
    add_torsion_terms(
        &mut field,
        &eargs.details,
        params.use_basic_knowledge,
        params.use_exp_torsion_angle_prefs,
    )?;
    if field.calc_energy() > ERROR_TOL {
        let mut iter = 0;
        let mut need_more = 1;
        while need_more != 0 && iter < 8 {
            if timed_out(end_time) {
                return Ok(false);
            }
            need_more = field.minimize(100, params.optimizer_force_tol, 1e-6)?;
            iter += 1;
        }
    }
    copy_back(&field, positions);
    Ok(field.calc_energy() / (n as f64) < MAX_MINIMIZED_E_PER_ATOM_AIO)
}

/// The checks of the K terms (linear-angle constraints and inversions) on the
/// embedded 3D coordinates.
///
/// (Upstream also compares the total energy with the number of centers before
/// anything is summed, which can never fail; that comparison is left out.)
fn check_k_terms(positions: &[f64], eargs: &EmbedArgs) -> Result<bool> {
    let dim = eargs.dim();
    let mut positions3d = Vec::with_capacity(positions.len() / dim * 3);
    for chunk in positions.chunks_exact(dim) {
        positions3d.extend_from_slice(&chunk[..3]);
    }
    let (angle_energies, inversion_energies) = improper_term_energies(
        &eargs.mmat,
        &positions3d,
        &eargs.details.improper_atoms,
        &eargs.details.angles,
    )?;
    Ok(!angle_energies.iter().any(|&e| e > 0.05) && !inversion_energies.iter().any(|&e| e > 5.0))
}

/// Double bonds with stereo have their substituents on the right side.
fn double_bond_stereo_checks(positions: &[f64], eargs: &EmbedArgs, failures: &Failures) -> bool {
    let dim = eargs.dim();
    for (atoms, sign) in &eargs.stereo_double_bonds {
        let p = atoms.map(|a| point(positions, dim, a as usize));
        // the dihedral is between 0 and 180; subtracting 90 gives positive
        // values closer to trans and negative ones closer to cis, and the
        // sign is 1 for trans and -1 for cis
        let dihedral = compute_dihedral_angle(&p[0], &p[1], &p[2], &p[3]);
        if (dihedral - std::f64::consts::FRAC_PI_2) * f64::from(*sign) < 0.0 {
            failures.record(EmbedFailureCause::BadDoubleBondStereo);
            return false;
        }
    }
    true
}

/// The "distance matrix" chirality test.
fn dist_matrix_chirality_test(positions: &[f64], eargs: &EmbedArgs, failures: &Failures) -> bool {
    let mut atoms = std::collections::BTreeSet::new();
    for set in &eargs.chiral_centers {
        if set.idx0 != set.idx4 {
            atoms.extend([set.idx0, set.idx1, set.idx2, set.idx3, set.idx4]);
        }
    }
    let atoms: Vec<usize> = atoms.into_iter().collect();
    if !atoms.is_empty() && !bounds_fulfilled(&atoms, &eargs.mmat, positions, eargs.dim()) {
        failures.record(EmbedFailureCause::FinalChiralBounds);
        return false;
    }
    true
}

fn check_final_center_in_volume(positions: &[f64], eargs: &EmbedArgs, failures: &Failures) -> bool {
    for set in &eargs.chiral_centers {
        // the center may be outside the volume of the other four points
        if !center_in_volume(set, positions, eargs.dim(), 0.1) {
            failures.record(EmbedFailureCause::FinalCenterInVolume);
            return false;
        }
    }
    true
}

fn final_chiral_checks(positions: &[f64], eargs: &EmbedArgs, failures: &Failures) -> bool {
    // confirm the chiral volumes
    if !check_chiral_centers(&eargs.chiral_centers, positions, eargs.dim()) {
        failures.record(EmbedFailureCause::CheckChiralCenters2);
        return false;
    }
    if !dist_matrix_chirality_test(positions, eargs, failures) {
        return false;
    }
    // the "center in volume" chirality test
    check_final_center_in_volume(positions, eargs, failures)
}

/// A generator seeded with `seed`, or from the global generator for -1.
pub fn seed_rng(seed: i32) -> Result<MinStdRand> {
    if seed < -1 {
        return Err(RDKitError::invariant(
            "random seed must either be positive, zero, or negative one",
        ));
    }
    if seed > -1 {
        let mut rng = MinStdRand::new(42);
        rng.seed(seed as u32);
        Ok(rng)
    } else {
        // no seed: draw one from the global generator
        Ok(MinStdRand::new(with_global_generator(0, |g| g.next())))
    }
}

struct EmbedContext<'a> {
    eargs: &'a EmbedArgs,
    params: &'a EmbedParameters,
    basin_thresh: f64,
    end_time: Option<Instant>,
    failures: &'a Failures,
}

/// The sequential (legacy) embedding of one conformer.
fn embed_points(positions: &mut [f64], ctx: &EmbedContext<'_>, seed: i32) -> Result<bool> {
    let EmbedContext {
        eargs,
        params,
        basin_thresh,
        end_time,
        failures,
    } = *ctx;
    let dim = eargs.dim();
    let n = positions.len() / dim;
    let mut dist_mat = SymmMatrix::new(n);
    let mut rng = seed_rng(seed)?;

    let mut got_coords = false;
    let mut iter = 0;
    while !got_coords && iter < params.max_iterations {
        if timed_out(end_time) {
            break;
        }
        iter += 1;
        if let Some(callback) = params.callback {
            callback(iter);
        }
        got_coords = generate_initial_coords(positions, eargs, params, &mut dist_mat, &mut rng)?;
        if !got_coords {
            failures.record(EmbedFailureCause::InitialCoords);
            continue;
        }
        // we only return the initial embedding if successful
        if params.only_initial_embedding {
            return Ok(got_coords);
        }
        got_coords = first_minimization(positions, eargs, params, basin_thresh)?;
        if !got_coords {
            failures.record(EmbedFailureCause::FirstMinimization);
        } else {
            got_coords = check_tetrahedral_centers(
                &eargs.tetrahedral_carbons,
                positions,
                dim,
                TETRAHEDRAL_CENTERINVOLUME_TOL,
            );
            if !got_coords {
                failures.record(EmbedFailureCause::CheckTetrahedralCenters);
            }
        }

        // are any chiral centers badly out of whack?
        if got_coords && params.enforce_chirality && !eargs.chiral_centers.is_empty() {
            got_coords = check_chiral_centers(&eargs.chiral_centers, positions, dim);
            if !got_coords {
                failures.record(EmbedFailureCause::CheckChiralCenters);
            }
        }
        // redo the minimization if we have a chiral center or started from
        // random coordinates
        if got_coords && (!eargs.chiral_centers.is_empty() || params.use_random_coords) {
            got_coords =
                minimize_fourth_dimension(positions, eargs, params, basin_thresh, end_time)?;
            if !got_coords {
                if timed_out(end_time) {
                    failures.record(EmbedFailureCause::ExceededTimeout);
                }
                failures.record(EmbedFailureCause::MinimizeFourthDimension);
            }
        }

        // (ET)(K)DG
        if got_coords && (params.use_exp_torsion_angle_prefs || params.use_basic_knowledge) {
            got_coords = minimize_with_exp_torsions(positions, eargs, params)?;
            if !got_coords {
                failures.record(EmbedFailureCause::EtkMinimization);
            }
        }
        if got_coords {
            got_coords = double_bond_geometry_checks(
                positions,
                dim,
                &eargs.double_bond_ends,
                &[],
                1e-3,
                false,
            );
            if !got_coords {
                failures.record(EmbedFailureCause::LinearDoubleBond);
            }
        }
        // is the stereo correct?
        if params.enforce_chirality && got_coords {
            if !eargs.chiral_centers.is_empty() {
                // any failures are tracked in the checks
                got_coords = final_chiral_checks(positions, eargs, failures);
            }
            if got_coords && !eargs.stereo_double_bonds.is_empty() {
                got_coords = double_bond_stereo_checks(positions, eargs, failures);
            }
        }
    }
    Ok(got_coords)
}

/// The all-in-one embedding of one conformer.
fn embed_points_aio(positions: &mut [f64], ctx: &EmbedContext<'_>, seed: i32) -> Result<bool> {
    let EmbedContext {
        eargs,
        params,
        end_time,
        failures,
        ..
    } = *ctx;
    let dim = eargs.dim();
    let n = positions.len() / dim;
    let mut dist_mat = SymmMatrix::new(n);
    let mut rng = seed_rng(seed)?;
    let improper: Vec<[u32; 4]> = eargs
        .details
        .improper_atoms
        .iter()
        .map(|i| i.atoms)
        .collect();

    let mut iter = 0;
    while iter < params.max_iterations {
        iter += 1;
        if timed_out(end_time) {
            break;
        }
        if let Some(callback) = params.callback {
            callback(iter);
        }
        // initial positions
        if !generate_initial_coords(positions, eargs, params, &mut dist_mat, &mut rng)? {
            failures.record(EmbedFailureCause::InitialCoords);
            continue;
        }
        if params.only_initial_embedding {
            return Ok(true);
        }
        // minimization
        if !minimize_all_in_one(positions, eargs, params, end_time)? {
            failures.record(if timed_out(end_time) {
                EmbedFailureCause::ExceededTimeout
            } else {
                EmbedFailureCause::Minimization
            });
            continue;
        }
        // the K terms
        if params.use_basic_knowledge && !check_k_terms(positions, eargs)? {
            failures.record(EmbedFailureCause::KtermViolation);
            continue;
        }
        // the tetrahedral centers
        if !check_tetrahedral_centers(
            &eargs.tetrahedral_carbons,
            positions,
            dim,
            TETRAHEDRAL_CENTERINVOLUME_TOL,
        ) {
            failures.record(EmbedFailureCause::CheckTetrahedralCenters);
            continue;
        }
        // the chiral centers
        if params.enforce_chirality {
            if !eargs.chiral_centers.is_empty() {
                if !check_chiral_centers(&eargs.chiral_centers, positions, dim) {
                    failures.record(EmbedFailureCause::CheckChiralCenters);
                    continue;
                }
                if !dist_matrix_chirality_test(positions, eargs, failures) {
                    continue;
                }
                if !check_final_center_in_volume(positions, eargs, failures) {
                    continue;
                }
            }
            if !eargs.stereo_double_bonds.is_empty()
                && !double_bond_stereo_checks(positions, eargs, failures)
            {
                continue;
            }
        }
        // clashing atoms
        if !clash_check(positions, dim, &eargs.mmat, 0.8) {
            failures.record(EmbedFailureCause::Clash);
            continue;
        }
        // the double bond geometry
        if !double_bond_geometry_checks(
            positions,
            dim,
            &eargs.double_bond_ends,
            &improper,
            3e-3,
            true,
        ) {
            failures.record(EmbedFailureCause::LinearDoubleBond);
            continue;
        }
        return Ok(true);
    }
    Ok(false)
}

/// The seed of conformer `ci` (0-based) for the parameters' seed.
fn conformer_seed(random_seed: i32, ci: usize, sequential: bool) -> Result<i32> {
    if random_seed < -1 {
        return Err(RDKitError::invariant(
            "random seed must either be positive, zero, or negative one",
        ));
    }
    if random_seed <= -1 {
        return Ok(random_seed);
    }
    if sequential {
        return Ok(random_seed.wrapping_add(ci as i32 + 1));
    }
    if let Some(seed) = i32::try_from(ci + 1)
        .ok()
        .and_then(|c| c.checked_mul(random_seed))
    {
        // the old method of computing a new seed
        return Ok(seed);
    }
    // the simple product overflows: hash the conformer index and the seed
    // together
    let big_seed = random_seed as u64;
    let ci1 = (ci + 1) as u64;
    let max_val = ci1.max(big_seed);
    let big_num = big_seed.wrapping_add(max_val.wrapping_mul(ci1));
    // the first 31 bits xor'd with the next 31 bits, to be positive
    const POSITIVE_INT_MASK: u64 = 0x7fff_ffff;
    let folded = (big_num & POSITIVE_INT_MASK) ^ (big_num >> 31);
    Ok((folded & POSITIVE_INT_MASK) as i32)
}

/// Embeds the conformers `ci` with `ci % num_threads == thread_id`; returns
/// the positions of each, or `None` for a failure.
fn embed_helper(
    thread_id: usize,
    num_threads: usize,
    eargs: &EmbedArgs,
    params: &EmbedParameters,
    basin_thresh: f64,
    confs_ok: &[bool],
    end_time: Option<Instant>,
    failures: &Failures,
) -> Result<Vec<(usize, Option<Vec<f64>>)>> {
    let dim = eargs.dim();
    let n = eargs.mmat.num_rows();
    let mut positions = vec![0.0; n * dim];
    let ctx = EmbedContext {
        eargs,
        params,
        basin_thresh,
        end_time,
        failures,
    };
    let mut results = Vec::new();
    for ci in 0..confs_ok.len() {
        if timed_out(end_time) {
            break;
        }
        if ci % num_threads != thread_id {
            continue;
        }
        if !confs_ok[ci] {
            // one of the fragments has already failed: no sense embedding
            // this one
            continue;
        }
        let seed = conformer_seed(
            params.random_seed,
            ci,
            params.enable_sequential_random_seeds,
        )?;
        let got_coords = if params.use_legacy_implementation {
            embed_points(&mut positions, &ctx, seed)?
        } else {
            embed_points_aio(&mut positions, &ctx, seed)?
        };
        results.push((ci, got_coords.then(|| positions.clone())));
    }
    Ok(results)
}

/// The atom matches used to align conformers when pruning: with symmetry all
/// heavy-atom automorphisms, otherwise the (heavy) atoms in order.
fn get_mol_self_matches(mol: &ROMol, params: &EmbedParameters) -> Result<Vec<Vec<u32>>> {
    if params.prune_rms_thresh != 0.0 && params.use_symmetry_for_pruning {
        let mut tmol = mol.clone();
        remove_hs(&mut tmol, &RemoveHsParameters::default(), false)?;
        let mut symmetrized;
        let probe: &ROMol = if params.symmetrize_conjugated_terminal_groups_for_pruning {
            symmetrized = tmol.clone();
            symmetrize_terminal_atoms(&mut symmetrized)?;
            &symmetrized
        } else {
            &tmol
        };
        let stripped = substruct_match(
            mol,
            probe,
            &SubstructMatchParameters {
                max_matches: 1,
                ..Default::default()
            },
        )?;
        if stripped.len() != 1 {
            return Err(RDKitError::invariant("expected match not found"));
        }
        let heavy_atom_matches = substruct_match(
            &tmol,
            probe,
            &SubstructMatchParameters {
                max_matches: 1000,
                uniquify: false,
                ..Default::default()
            },
        )?;
        Ok(heavy_atom_matches
            .iter()
            .map(|m| {
                m.iter()
                    .map(|&(_, mol_idx)| stripped[0][mol_idx as usize].1 as u32)
                    .collect()
            })
            .collect())
    } else if params.only_heavy_atoms_for_rms {
        Ok(vec![
            (0..mol.num_atoms())
                .filter(|&i| mol.atom(i).atomic_num() != 1)
                .collect(),
        ])
    } else {
        Ok(vec![(0..mol.num_atoms()).collect()])
    }
}

/// Whether the conformer is at least `threshold` away from all the
/// conformers of the molecule for every match.
fn is_conf_far_from_rest(
    mol: &ROMol,
    conf: &Conformer,
    threshold: f64,
    self_matches: &[Vec<u32>],
) -> Result<bool> {
    let pick = |c: &Conformer, m: &[u32]| -> Vec<Point3D> {
        m.iter().map(|&i| c.positions()[i as usize]).collect()
    };
    let ref_points = pick(conf, &self_matches[0]);
    let ssr_thres = self_matches[0].len() as f64 * threshold * threshold;
    for m in self_matches {
        for other in mol.conformers() {
            let probe = pick(other, m);
            let (ssr, _) = align_points(&ref_points, &probe, None, false, 50)?;
            if ssr < ssr_thres {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// What [`embed_multiple_confs`] does for the fragments, apart from setting
/// up the result and the failure counters.
fn embed_fragments(
    mol: &mut ROMol,
    num_confs: u32,
    params: &mut EmbedParameters,
    failures: &Failures,
    res: &mut Vec<i32>,
    end_time: Option<Instant>,
) -> Result<()> {
    if mol.num_atoms() == 0 {
        return Err(RDKitError::Value("molecule has no atoms".into()));
    }
    if !matches!(params.et_version, 1 | 2 | 4) {
        return Err(RDKitError::Value(
            "Only version 1, 2 and 4 of the experimental torsion-angle preferences (ETversion) supported"
                .into(),
        ));
    }
    if needs_hs(mol)? {
        rd_log_warning!("Molecule does not have explicit Hs. Consider calling AddHs()");
    }
    // DEPRECATED
    if params.use_random_coords {
        params.initial_embedding_mode = InitialEmbeddingMode::RandomCoordinateEmbedding;
    }

    // initialize the conformers to create
    if params.clear_confs {
        res.clear();
        mol.clear_conformers();
    }
    let n_total = mol.num_atoms() as usize;
    let mut confs: Vec<Conformer> = (0..num_confs).map(|_| Conformer::new(n_total)).collect();
    let mut confs_ok = vec![true; num_confs as usize];

    let (frag_pieces, frag_mapping): (Vec<ROMol>, Vec<i32>) = if params.embed_fragments_separately {
        let (pieces, _) = get_mol_frags_as_mols(mol, true, false)?;
        (pieces, get_mol_frags_mapping(mol).1)
    } else {
        (vec![mol.clone()], vec![0; n_total])
    };
    let mut coord_map: Option<CoordMap> = params.coord_map.clone();
    if frag_pieces.len() > 1 && coord_map.is_some() {
        rd_log_warning!(
            "Constrained conformer generation (via the coordMap argument) does not work with molecules that have multiple fragments."
        );
        coord_map = None;
    }
    let mut constrained_atoms = vec![false; n_total];
    if let Some(map) = &coord_map {
        for (&idx, point) in map {
            if idx < 0 || idx as usize >= n_total {
                return Err(RDKitError::Value(
                    "coordMap atom index is out of range".into(),
                ));
            }
            if !(point.x.is_finite() && point.y.is_finite() && point.z.is_finite()) {
                return Err(RDKitError::Value(
                    "coordMap contains non-finite coordinates".into(),
                ));
            }
            constrained_atoms[idx as usize] = true;
        }
    }
    if frag_pieces.len() > 1 && params.bounds_mat.is_some() {
        rd_log_warning!(
            "Conformer generation using a user-provided boundsMat does not work with molecules that have multiple fragments. The boundsMat will be ignored."
        );
        coord_map = None; // FIXME (upstream): this should reset the bounds matrix
    }

    // embed the fragments one by one
    for (frag_idx, mut piece) in frag_pieces.into_iter().enumerate() {
        let n_atoms = piece.num_atoms() as usize;
        let ic_mode =
            params.initial_embedding_mode == InitialEmbeddingMode::InternalCoordinateEmbedding;
        let mut details = CrystalFFDetails {
            constrained_atoms: constrained_atoms.clone(),
            dist_mat: get_distance_mat(&piece),
            internal_coords: ic_mode
                .then(|| Box::new(InternalCoordinates::new(piece.num_bonds() as usize))),
            ..Default::default()
        };

        let mmat = match &params.bounds_mat {
            Some(user) if frag_pieces_len_is_one(&frag_mapping) => {
                // use what the user gave us, if it has the right size
                if user.num_rows() != n_atoms {
                    return Err(RDKitError::Value(
                        "size of boundsMat provided does not match the number of atoms in the molecule."
                            .into(),
                    ));
                }
                if ic_mode {
                    // collect the internal coordinates; for now the ones of
                    // the bounds matrix are used
                    let mut scratch = BoundsMatrix::new(n_atoms);
                    set_topol_bounds(
                        &piece,
                        &mut scratch,
                        params,
                        &TopolBoundsOptions {
                            embed_force_field: params.embed_force_field,
                            ..Default::default()
                        },
                        TopolBoundsOutput {
                            paths14: None,
                            internal_coords: details.internal_coords.as_deref_mut(),
                        },
                    )?;
                }
                let (bonds, angles) = collect_bonds_and_angles(&piece);
                details.bonds = bonds;
                details.angles = angles;
                (**user).clone()
            }
            _ => {
                // create and initialize the distance bounds matrix
                let mut mmat = BoundsMatrix::new(n_atoms);
                init_bounds_mat(&mut mmat, 0.0, 1000.0)?;
                if !setup_initial_bounds_matrix(
                    &piece,
                    &mut mmat,
                    coord_map.as_ref(),
                    params,
                    &mut details,
                )? {
                    // e.g. a triangle smoothing failure
                    return Ok(());
                }
                mmat
            }
        };
        init_etkdg(&piece, params, &mut details)?;

        // all the chiral centers
        assign_stereochemistry(&mut piece, false, false, false)?;
        let (chiral_centers, tetrahedral_carbons) = find_chiral_sets(&piece, coord_map.as_ref())?;

        let mut zmat = ZMatrix::new(n_atoms);
        if let Some(coords) = details.internal_coords.as_deref() {
            set_molecule_dfs(&piece, &mut zmat, coords)?;
            correct_chiral_centers(&piece, &mut zmat)?;
        }
        let (double_bond_ends, stereo_double_bonds) =
            find_double_bonds(&piece, coord_map.as_ref())?;

        // with chiral centers or random coordinates the molecule is first
        // embedded in four dimensions, otherwise in 3D
        let four_d = params.use_random_coords
            || !chiral_centers.is_empty()
            || (ic_mode && (!chiral_centers.is_empty() || !tetrahedral_carbons.is_empty()));

        let hac = piece.num_heavy_atoms() as usize;
        let eargs = EmbedArgs {
            four_d,
            mmat,
            zmat,
            chiral_centers,
            tetrahedral_carbons,
            double_bond_ends,
            stereo_double_bonds,
            details,
            hac,
        };

        // the maximum number of iterations and the basin threshold are fixed
        // by the first embedding (like upstream, they stay for later
        // fragments)
        if params.max_iterations == 0 {
            params.max_iterations = 10 * n_atoms as u32;
        }
        // the basin threshold ignores 1-4 (and higher) interactions with
        // random coordinates, folding up large flexible molecules
        if params.use_random_coords {
            params.basin_thresh = 1e8;
        }

        let num_threads = num_threads_to_use(params.num_threads);
        let thread_ids: Vec<usize> = (0..num_threads).collect();
        let snapshot = confs_ok.clone();
        let params_ref: &EmbedParameters = params;
        let results = par_map(&thread_ids, num_threads as i32, |&tid| {
            embed_helper(
                tid,
                num_threads,
                &eargs,
                params_ref,
                params_ref.basin_thresh,
                &snapshot,
                end_time,
                failures,
            )
        });
        let dim = eargs.dim();
        for result in results {
            for (ci, positions) in result? {
                match positions {
                    Some(positions) => {
                        // copy the coordinates into the conformer
                        let mut frag_atom_idx = 0;
                        for i in 0..n_total {
                            if frag_mapping[i] == frag_idx as i32 {
                                let p = &positions[dim * frag_atom_idx..dim * frag_atom_idx + 3];
                                confs[ci].set_atom_pos(i as u32, Point3D::new(p[0], p[1], p[2]))?;
                                frag_atom_idx += 1;
                            }
                        }
                    }
                    None => confs_ok[ci] = false,
                }
            }
        }
        if timed_out(end_time) {
            failures.record(EmbedFailureCause::ExceededTimeout);
            res.push(-1);
            return Ok(());
        }
    }

    let self_matches = if params.prune_rms_thresh > 0.0 {
        get_mol_self_matches(mol, params)?
    } else {
        Vec::new()
    };
    for (ci, conf) in confs.into_iter().enumerate() {
        if !confs_ok[ci] {
            continue;
        }
        // is a close-by conformation already chosen?
        if params.prune_rms_thresh <= 0.0
            || is_conf_far_from_rest(mol, &conf, params.prune_rms_thresh, &self_matches)?
        {
            let conf_id = mol.add_conformer(conf, true)?;
            res.push(conf_id as i32);
        }
    }
    Ok(())
}

/// A molecule with several fragments embeds them separately, ignoring the
/// bounds matrix of the user.
fn frag_pieces_len_is_one(frag_mapping: &[i32]) -> bool {
    frag_mapping.iter().all(|&f| f == 0)
}

/// Embeds `num_confs` conformers of the molecule; returns their ids (empty if
/// the embedding fails).
///
/// With several fragments (and `embed_fragments_separately`) they are
/// embedded separately and likely occupy the same region of space.
pub fn embed_multiple_confs(
    mol: &mut ROMol,
    num_confs: u32,
    params: &mut EmbedParameters,
) -> Result<Vec<i32>> {
    let end_time = (params.timeout > 0)
        .then(|| Instant::now() + Duration::from_secs(u64::from(params.timeout)));
    let failures = Failures::new(params.track_failures);
    let mut res = Vec::new();
    let outcome = embed_fragments(mol, num_confs, params, &failures, &mut res, end_time);
    if let Some(counts) = failures.counts {
        params.failures = counts.into_inner().unwrap_or_else(|e| e.into_inner());
    }
    outcome?;
    Ok(res)
}

/// Embeds one conformer; returns its id, or -1 if the embedding failed.
pub fn embed_molecule(mol: &mut ROMol, params: &mut EmbedParameters) -> Result<i32> {
    Ok(embed_multiple_confs(mol, 1, params)?
        .first()
        .copied()
        .unwrap_or(-1))
}
