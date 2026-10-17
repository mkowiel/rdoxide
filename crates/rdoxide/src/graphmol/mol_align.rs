//! Ported from RDKit `Code/GraphMol/MolAlign/AlignMolecules.{h,cpp}` and the
//! `randomTransform` helper of `O3AAlignMolecules.cpp`.
//!
//! Molecule alignment on a given or substructure-derived atom map, with
//! symmetry-aware best-RMS searches. The per-match and per-conformer-pair
//! loops run on the shared [`par_map`] helper, so the thread count only
//! affects speed, never the result.

use crate::geometry::{Point3D, Transform3D};
use crate::graphmol::bond::{Bond, BondType};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::mol_transforms::{compute_centroid, transform_conformer};
use crate::graphmol::queries::make_single_or_double_bond_query;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substruct::{MatchVect, SubstructMatchParameters, substruct_match};
use crate::numerics::alignment::align_points;
use crate::rdgeneral::random::{get_random_val, with_global_generator};
use crate::rdgeneral::threads::{num_threads_to_use, par_map};
use crate::{RDKitError, Result, rd_log_warning};

/// `O3_RANDOM_TRANS_COEFF`
const RANDOM_TRANS_COEFF: f64 = 5.0;
const DEFAULT_MAX_ITERS: u32 = 50;

/// `MolAlign::BestAlignmentParams`.
#[derive(Debug, Clone)]
pub struct BestAlignmentParams {
    /// Maximum number of substructure matches to consider.
    pub max_matches: u32,
    /// Treat terminal O/N of conjugated groups (e.g. carboxylates) as
    /// equivalent.
    pub symmetrize_conjugated_terminal_groups: bool,
    /// Leave hydrogens out of the atom maps.
    pub ignore_hs: bool,
    /// Worker threads for the match / conformer loops (`<= 0`: all but
    /// `|n|` hardware threads).
    pub num_threads: i32,
    /// Explicit `(probe, reference)` atom maps; empty = find them by
    /// substructure matching.
    pub map: Vec<MatchVect>,
    pub weights: Option<Vec<f64>>,
}

impl Default for BestAlignmentParams {
    fn default() -> Self {
        BestAlignmentParams {
            max_matches: 1_000_000,
            symmetrize_conjugated_terminal_groups: true,
            ignore_hs: true,
            num_threads: 0,
            map: Vec::new(),
            weights: None,
        }
    }
}

fn align_exception(message: &str) -> RDKitError {
    RDKitError::Runtime(message.to_string())
}

/// `details::symmetrizeTerminalAtoms`: makes the bond to a terminal O/N of a
/// conjugated group "single or double" and neutralizes that atom, so both
/// terminal atoms of e.g. a carboxylate match either way round.
pub fn symmetrize_terminal_atoms(mol: &mut RWMol) -> Result<()> {
    const SMARTS: &str = "[{atomPattern};$([{atomPattern}]-[*]=[{atomPattern}]),$([{atomPattern}]=[*]-[{atomPattern}])]~[*]";
    let mut params = SmartsParserParams::default();
    params
        .replacements
        .insert("{atomPattern}".into(), "O,N;D1".into());
    let query = mol_from_smarts(SMARTS, &params)?
        .ok_or_else(|| RDKitError::invariant("bad query pattern"))?;
    let matches = substruct_match(mol, &query, &SubstructMatchParameters::default())?;
    if matches.is_empty() {
        return Ok(());
    }
    let mut replacement = Bond::new(BondType::Unspecified);
    replacement.set_query(Some(make_single_or_double_bond_query()));
    for m in &matches {
        let (terminal, other) = (m[0].1 as u32, m[1].1 as u32);
        mol.atom_mut(terminal).set_formal_charge(0);
        let bond = mol
            .bond_idx_between(terminal, other)
            .ok_or_else(|| RDKitError::invariant("could not find expected bond"))?;
        mol.replace_bond(bond, &replacement, false, true)?;
    }
    Ok(())
}

/// `getAllMatchesPrbRef`: every `(probe, reference)` atom mapping of the
/// probe onto the reference.
fn all_matches_prb_ref(
    prb_mol: &ROMol,
    ref_mol: &ROMol,
    max_matches: u32,
    symmetrize: bool,
    ignore_hs: bool,
) -> Result<Vec<MatchVect>> {
    let symmetrized = if symmetrize {
        let mut m = prb_mol.clone();
        symmetrize_terminal_atoms(&mut m)?;
        Some(m)
    } else {
        None
    };
    let query = symmetrized.as_ref().unwrap_or(prb_mol);
    let params = SubstructMatchParameters {
        uniquify: false,
        recursion_possible: true,
        use_chirality: false,
        use_query_query_matches: false,
        max_matches,
        ..Default::default()
    };
    let mut matches = substruct_match(ref_mol, query, &params)?;
    if matches.is_empty() {
        return Err(align_exception(
            "No sub-structure match found between the reference and probe mol",
        ));
    }
    if matches.len() > 1_000_000 {
        let name = prb_mol.name();
        rd_log_warning!(
            "Warning in getAllMatchesPrbRef: {} matches detected for molecule {}, this may lead to a performance slowdown.",
            matches.len(),
            name
        );
    }
    if ignore_hs {
        for m in &mut matches {
            m.retain(|&(prb, _)| query.atom(prb as u32).atomic_num() != 1);
        }
    }
    Ok(matches)
}

fn atom_map_points(
    prb: &Conformer,
    reference: &Conformer,
    atom_map: &MatchVect,
) -> Result<(Vec<Point3D>, Vec<Point3D>)> {
    let mut prb_pts = Vec::with_capacity(atom_map.len());
    let mut ref_pts = Vec::with_capacity(atom_map.len());
    for &(p, r) in atom_map {
        prb_pts.push(prb.atom_pos(p as u32)?);
        ref_pts.push(reference.atom_pos(r as u32)?);
    }
    Ok((prb_pts, ref_pts))
}

/// `alignConfsOnAtomMap`: the mean squared deviation after the optimal
/// alignment, and the transformation.
fn align_confs_on_atom_map(
    prb: &Conformer,
    reference: &Conformer,
    atom_map: &MatchVect,
    weights: Option<&[f64]>,
    reflect: bool,
    max_iters: u32,
) -> Result<(f64, Transform3D)> {
    let (prb_pts, ref_pts) = atom_map_points(prb, reference, atom_map)?;
    let (ssr, trans) = align_points(&ref_pts, &prb_pts, weights, reflect, max_iters)?;
    Ok((ssr / prb_pts.len() as f64, trans))
}

/// `calcMSDInternal`: the mean squared deviation without any alignment.
fn calc_msd(
    prb: &Conformer,
    reference: &Conformer,
    atom_map: &MatchVect,
    weights: Option<&[f64]>,
) -> Result<f64> {
    let npt = atom_map.len();
    if let Some(w) = weights
        && w.len() != npt
    {
        return Err(RDKitError::precondition("Mismatch in number of weights"));
    }
    let (prb_pts, ref_pts) = atom_map_points(prb, reference, atom_map)?;
    let ssr: f64 = (0..npt)
        .map(|i| weights.map_or(1.0, |w| w[i]) * (prb_pts[i] - ref_pts[i]).length_sq())
        .sum();
    Ok(ssr / npt as f64)
}

/// The best match found by [`best_rms_internal`].
struct BestFit {
    rms: f64,
    /// Present when the matches were aligned (not just scored).
    trans: Option<Transform3D>,
    match_index: usize,
}

/// `getBestRMSInternal`: the smallest RMSD over `matches`. With `align` the
/// probe is optimally superposed for each match, otherwise it is scored in
/// place. Ties go to the earliest match whatever the thread count.
#[allow(clippy::too_many_arguments)]
fn best_rms_internal(
    prb_mol: &ROMol,
    ref_mol: &ROMol,
    prb_cid: i32,
    ref_cid: i32,
    matches: &[MatchVect],
    align: bool,
    weights: Option<&[f64]>,
    reflect: bool,
    max_iters: u32,
    num_threads: usize,
) -> Result<BestFit> {
    if matches.is_empty() {
        return Err(RDKitError::precondition("matches must not be empty"));
    }
    let prb = prb_mol.conformer(prb_cid)?;
    let reference = ref_mol.conformer(ref_cid)?;
    let target = i32::try_from(num_threads).unwrap_or(i32::MAX);
    let results = par_map(matches, target, |m| -> Result<(f64, Option<Transform3D>)> {
        if align {
            let (msd, t) = align_confs_on_atom_map(prb, reference, m, weights, reflect, max_iters)?;
            Ok((msd, Some(t)))
        } else {
            Ok((calc_msd(prb, reference, m, weights)?, None))
        }
    });
    let mut best: Option<(f64, usize, Option<Transform3D>)> = None;
    for (i, r) in results.into_iter().enumerate() {
        let (msd, trans) = r?;
        if best.as_ref().is_none_or(|b| msd < b.0) {
            best = Some((msd, i, trans));
        }
    }
    let (msd, match_index, trans) = best.unwrap_or((f64::MAX, 0, None));
    Ok(BestFit {
        rms: msd.sqrt(),
        trans,
        match_index,
    })
}

/// The matches to use: explicit ones or those found by substructure search.
fn resolve_matches(
    prb_mol: &ROMol,
    ref_mol: &ROMol,
    params: &BestAlignmentParams,
) -> Result<Vec<MatchVect>> {
    if params.map.is_empty() {
        all_matches_prb_ref(
            prb_mol,
            ref_mol,
            params.max_matches,
            params.symmetrize_conjugated_terminal_groups,
            params.ignore_hs,
        )
    } else {
        Ok(params.map.clone())
    }
}

/// `getAlignmentTransform`: the RMSD and the transformation taking the probe
/// conformer onto the reference one. Without `atom_map` the probe is matched
/// as a substructure of the reference.
#[allow(clippy::too_many_arguments)]
pub fn get_alignment_transform(
    prb_mol: &ROMol,
    ref_mol: &ROMol,
    prb_cid: i32,
    ref_cid: i32,
    atom_map: Option<&MatchVect>,
    weights: Option<&[f64]>,
    reflect: bool,
    max_iters: u32,
) -> Result<(f64, Transform3D)> {
    let prb = prb_mol.conformer(prb_cid)?;
    let reference = ref_mol.conformer(ref_cid)?;
    let found;
    let atom_map = match atom_map {
        Some(m) => m,
        None => {
            let params = SubstructMatchParameters {
                recursion_possible: true,
                use_chirality: false,
                use_query_query_matches: true,
                max_matches: 1,
                ..Default::default()
            };
            let mut ms = substruct_match(ref_mol, prb_mol, &params)?;
            if ms.is_empty() {
                return Err(align_exception(
                    "No sub-structure match found between the probe and query mol",
                ));
            }
            found = ms.swap_remove(0);
            &found
        }
    };
    let (msd, trans) =
        align_confs_on_atom_map(prb, reference, atom_map, weights, reflect, max_iters)?;
    Ok((msd.sqrt(), trans))
}

/// `alignMol`: aligns the probe conformer onto the reference one and returns
/// the RMSD.
#[allow(clippy::too_many_arguments)]
pub fn align_mol(
    prb_mol: &mut ROMol,
    ref_mol: &ROMol,
    prb_cid: i32,
    ref_cid: i32,
    atom_map: Option<&MatchVect>,
    weights: Option<&[f64]>,
    reflect: bool,
    max_iters: u32,
) -> Result<f64> {
    let (rms, trans) = get_alignment_transform(
        prb_mol, ref_mol, prb_cid, ref_cid, atom_map, weights, reflect, max_iters,
    )?;
    transform_conformer(prb_mol.conformer_mut(prb_cid)?, &trans);
    Ok(rms)
}

/// `getBestAlignmentTransform`: the lowest RMSD over all (symmetry-related)
/// atom maps, with its transformation and map.
pub fn get_best_alignment_transform(
    prb_mol: &ROMol,
    ref_mol: &ROMol,
    params: &BestAlignmentParams,
    prb_cid: i32,
    ref_cid: i32,
    reflect: bool,
    max_iters: u32,
) -> Result<(f64, Transform3D, MatchVect)> {
    let matches = resolve_matches(prb_mol, ref_mol, params)?;
    let best = best_rms_internal(
        prb_mol,
        ref_mol,
        prb_cid,
        ref_cid,
        &matches,
        true,
        params.weights.as_deref(),
        reflect,
        max_iters,
        num_threads_to_use(params.num_threads),
    )?;
    let trans = best.trans.unwrap_or_else(Transform3D::identity);
    Ok((best.rms, trans, matches[best.match_index].clone()))
}

/// `getBestRMS`: aligns the probe conformer with the lowest-RMSD atom map and
/// returns that RMSD.
pub fn get_best_rms(
    prb_mol: &mut ROMol,
    ref_mol: &ROMol,
    params: &BestAlignmentParams,
    prb_cid: i32,
    ref_cid: i32,
) -> Result<f64> {
    let matches = resolve_matches(prb_mol, ref_mol, params)?;
    if let (false, Some(w)) = (params.map.is_empty(), &params.weights)
        && w.len() != params.map[0].len()
    {
        return Err(RDKitError::Value(format!(
            "Mismatch between number of weights ({}) and number of atoms in the alignment map ({}).",
            w.len(),
            params.map[0].len()
        )));
    }
    let best = best_rms_internal(
        prb_mol,
        ref_mol,
        prb_cid,
        ref_cid,
        &matches,
        true,
        params.weights.as_deref(),
        false,
        DEFAULT_MAX_ITERS,
        num_threads_to_use(params.num_threads),
    )?;
    if let Some(trans) = best.trans {
        transform_conformer(prb_mol.conformer_mut(prb_cid)?, &trans);
    }
    Ok(best.rms)
}

fn conformer_ids(mol: &ROMol) -> Vec<i32> {
    mol.conformers().iter().map(|c| c.id() as i32).collect()
}

/// `getAllConformerBestRMS`: the best RMSD between every pair of conformers
/// of `mol`, in the order `(1,0), (2,0), (2,1), (3,0), …`.
pub fn get_all_conformer_best_rms(mol: &ROMol, params: &BestAlignmentParams) -> Result<Vec<f64>> {
    let matches = resolve_matches(mol, mol, params)?;
    let cids = conformer_ids(mol);
    let pairs: Vec<(i32, i32)> = (0..cids.len())
        .flat_map(|ci| (0..ci).map(move |cj| (ci, cj)))
        .map(|(ci, cj)| (cids[ci], cids[cj]))
        .collect();
    // pairs are the parallel unit; each one scans its matches serially
    par_map(&pairs, params.num_threads, |&(prb, reference)| {
        best_rms_internal(
            mol,
            mol,
            prb,
            reference,
            &matches,
            true,
            params.weights.as_deref(),
            false,
            DEFAULT_MAX_ITERS,
            1,
        )
        .map(|b| b.rms)
    })
    .into_iter()
    .collect()
}

/// `getAllConformerBestRMSToRef`: the best RMSD between every reference
/// conformer (outer loop) and every probe conformer (inner loop).
pub fn get_all_conformer_best_rms_to_ref(
    prb_mol: &ROMol,
    ref_mol: &ROMol,
    params: &BestAlignmentParams,
) -> Result<Vec<f64>> {
    let matches = resolve_matches(prb_mol, ref_mol, params)?;
    let prb_cids = conformer_ids(prb_mol);
    let pairs: Vec<(i32, i32)> = conformer_ids(ref_mol)
        .into_iter()
        .flat_map(|r| prb_cids.iter().map(move |&p| (p, r)))
        .collect();
    par_map(&pairs, params.num_threads, |&(prb, reference)| {
        best_rms_internal(
            prb_mol,
            ref_mol,
            prb,
            reference,
            &matches,
            true,
            params.weights.as_deref(),
            false,
            DEFAULT_MAX_ITERS,
            1,
        )
        .map(|b| b.rms)
    })
    .into_iter()
    .collect()
}

/// `CalcRMS`: the lowest RMSD over the atom maps without moving the probe.
#[allow(clippy::too_many_arguments)]
pub fn calc_rms(
    prb_mol: &ROMol,
    ref_mol: &ROMol,
    prb_cid: i32,
    ref_cid: i32,
    map: &[MatchVect],
    max_matches: u32,
    symmetrize_conjugated_terminal_groups: bool,
    weights: Option<&[f64]>,
) -> Result<f64> {
    let found;
    let matches = if map.is_empty() {
        found = all_matches_prb_ref(
            prb_mol,
            ref_mol,
            max_matches,
            symmetrize_conjugated_terminal_groups,
            false,
        )?;
        &found[..]
    } else {
        map
    };
    Ok(best_rms_internal(
        prb_mol,
        ref_mol,
        prb_cid,
        ref_cid,
        matches,
        false,
        weights,
        false,
        DEFAULT_MAX_ITERS,
        1,
    )?
    .rms)
}

fn selected_positions(conf: &Conformer, atom_ids: Option<&[u32]>) -> Result<Vec<Point3D>> {
    match atom_ids {
        None => Ok(conf.positions().to_vec()),
        Some(ids) => ids.iter().map(|&i| conf.atom_pos(i)).collect(),
    }
}

/// `alignMolConformers`: aligns the conformers of `mol` onto the first one
/// (of `conf_ids` when given) using `atom_ids` (default: all atoms). Returns
/// the per-conformer RMSDs, in alignment order.
pub fn align_mol_conformers(
    mol: &mut ROMol,
    atom_ids: Option<&[u32]>,
    conf_ids: Option<&[u32]>,
    weights: Option<&[f64]>,
    reflect: bool,
    max_iters: u32,
) -> Result<Vec<f64>> {
    let mut rms_list = Vec::new();
    if mol.num_conformers() == 0 {
        return Ok(rms_list);
    }
    let ref_id = conf_ids.and_then(|c| c.first()).map_or(-1, |&c| c as i32);
    let ref_pts = selected_positions(mol.conformer(ref_id)?, atom_ids)?;
    let targets: Vec<i32> = match conf_ids {
        None => conformer_ids(mol).into_iter().skip(1).collect(),
        Some(ids) => ids.iter().skip(1).map(|&c| c as i32).collect(),
    };
    for cid in targets {
        let conf = mol.conformer_mut(cid)?;
        let prb_pts = selected_positions(conf, atom_ids)?;
        let (ssd, trans) = align_points(&ref_pts, &prb_pts, weights, reflect, max_iters)?;
        rms_list.push((ssd / prb_pts.len() as f64).sqrt());
        transform_conformer(conf, &trans);
    }
    Ok(rms_list)
}

/// `randomTransform`: applies a random rotation about the centroid plus a
/// random translation (up to ±2.5 extents) to one conformer. A positive
/// `seed` reseeds the global generator.
pub fn random_transform(mol: &mut ROMol, cid: i32, seed: i32) -> Result<()> {
    // only reseeds; no value is drawn
    with_global_generator(seed, |_| ());
    let centroid = compute_centroid(mol, mol.conformer(cid)?, false, None)?;
    let conf = mol.conformer_mut(cid)?;
    let mut rnd_trans = [0.0; 3];
    let mut rnd_rot = [0.0; 3];
    for j in 0..3 {
        let coord = |p: &Point3D| [p.x, p.y, p.z][j];
        let max = conf.positions().iter().map(coord).fold(f64::MIN, f64::max);
        let min = conf.positions().iter().map(coord).fold(f64::MAX, f64::min);
        rnd_trans[j] = (max - min) * RANDOM_TRANS_COEFF * (get_random_val(-1) - 0.5);
        rnd_rot[j] = get_random_val(-1) * 2.0 * std::f64::consts::PI;
    }
    let axis = |x, y, z| Point3D::new(x, y, z);
    let mut to_origin = Transform3D::identity();
    to_origin.set_translation(&Point3D::new(-centroid.x, -centroid.y, -centroid.z));
    let mut x_rot = Transform3D::identity();
    x_rot.set_rotation(rnd_rot[0], &axis(1.0, 0.0, 0.0));
    let mut y_rot = Transform3D::identity();
    y_rot.set_rotation(rnd_rot[1], &axis(0.0, 1.0, 0.0));
    let mut z_rot_trans = Transform3D::identity();
    z_rot_trans.set_rotation(rnd_rot[2], &axis(0.0, 0.0, 1.0));
    z_rot_trans.set_translation(&Point3D::new(
        centroid.x + rnd_trans[0],
        centroid.y + rnd_trans[1],
        centroid.z + rnd_trans[2],
    ));
    let total = z_rot_trans
        .compose(&y_rot)
        .compose(&x_rot)
        .compose(&to_origin);
    transform_conformer(conf, &total);
    Ok(())
}
