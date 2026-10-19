//! Ported from RDKit `Code/GraphMol/MolAlign/O3AAlignMolecules.{h,cpp}`:
//! Open3DALIGN, the alignment of two molecules by atom-property based
//! pairing (MMFF94 types and charges, or Crippen logP) followed by a
//! weighted superposition.

mod lap;
mod sdm;

use std::sync::OnceLock;

pub use lap::{Lap, MolHistogram};
pub use sdm::{
    CrippenScheme, MmffScheme, O3_DEFAULT_CONSTRAINT_WEIGHT, O3_MAX_WEIGHT_COEFF,
    O3AConstraintVect, Scheme, Sdm,
};

use self::lap::{O3_DUMMY_COST, O3_LARGE_NEGATIVE_WEIGHT, O3_MAX_H_BINS};
use self::sdm::SdmContext;
use crate::geometry::{Point3D, Transform3D};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::forcefield_helpers::mmff::MmffMolProperties;
use crate::graphmol::matrices::get_3d_distance_mat;
use crate::graphmol::mol_align::{align_mol, get_alignment_transform};
use crate::graphmol::mol_transforms::transform_conformer;
use crate::graphmol::romol::ROMol;
use crate::graphmol::substruct::MatchVect;
use crate::rdgeneral::threads::par_map;
use crate::{RDKitError, Result};

const O3_MAX_SDM_ITERATIONS: u32 = 100;
const O3_MAX_SDM_THRESHOLD_ITER: u32 = 3;
const O3_SDM_THRESHOLD_START: f64 = 0.7;
const O3_SDM_THRESHOLD_STEP: f64 = 0.3;
const O3_RMSD_THRESHOLD: f64 = 1.0e-04;
const O3_SCORE_THRESHOLD: f64 = 0.01;

/// `O3_USE_MMFF_WEIGHTS`, `O3_ACCURACY_MASK` and `O3_LOCAL_ONLY`.
pub const O3_ACCURACY_MASK: u32 = 0b11;
pub const O3_LOCAL_ONLY: u32 = 1 << 2;

/// The 99x99 MMFF atom type similarity matrix.
fn mmff_similarity() -> &'static [[u8; 99]; 99] {
    static SIM: OnceLock<Box<[[u8; 99]; 99]>> = OnceLock::new();
    SIM.get_or_init(|| {
        let text = include_str!("../../../data/o3a_mmff_similarity.tsv");
        let mut sim = Box::new([[0u8; 99]; 99]);
        let rows = text.lines().filter(|l| !l.starts_with('#'));
        for (i, line) in rows.enumerate() {
            for (j, v) in line.split('\t').enumerate() {
                sim[i][j] = v.parse().expect("similarity entry");
            }
        }
        sim
    })
}

fn align_error(msg: &str) -> RDKitError {
    RDKitError::Runtime(msg.to_string())
}

/// The atom property data the two molecules are paired with.
#[allow(missing_debug_implementations)]
pub enum O3AProps<'a> {
    Mmff {
        prb: &'a MmffMolProperties,
        reference: &'a MmffMolProperties,
    },
    Crippen {
        prb: &'a [f64],
        reference: &'a [f64],
    },
}

/// Result of an O3A pairing (`O3A`).
#[derive(Debug, Clone)]
pub struct O3A {
    prb_cid: i32,
    ref_cid: i32,
    reflect: bool,
    max_iters: u32,
    /// `(probe, reference)` atom pairs
    pub matches: MatchVect,
    pub weights: Vec<f64>,
    pub score: f64,
}

/// Options shared by the constructors and the per-conformer helper.
#[derive(Debug, Clone)]
pub struct O3AParams<'a> {
    pub prb_cid: i32,
    pub ref_cid: i32,
    pub reflect: bool,
    pub max_iters: u32,
    pub options: u32,
    pub constraint_map: Option<&'a MatchVect>,
    pub constraint_weights: Option<&'a [f64]>,
}

fn atomic_nums(mol: &ROMol) -> Vec<u32> {
    (0..mol.num_atoms())
        .map(|i| mol.atom(i).atomic_num())
        .collect()
}

fn heavy_flags(nums: &[u32]) -> Vec<bool> {
    nums.iter().map(|&z| z != 1).collect()
}

fn positions(mol: &ROMol, cid: i32) -> Result<Vec<Point3D>> {
    Ok(mol.conformer(cid)?.positions().to_vec())
}

/// `LAP::computeCostMatrix`
#[allow(clippy::too_many_arguments)]
fn compute_cost_matrix(
    lap: &mut Lap,
    prb_nums: &[u32],
    prb_hist: &MolHistogram,
    ref_nums: &[u32],
    ref_hist: &MolHistogram,
    constraints: Option<&O3AConstraintVect>,
    scheme: &dyn Scheme,
) {
    let largest = ref_nums
        .iter()
        .filter(|&&z| z != 1)
        .count()
        .max(prb_nums.iter().filter(|&&z| z != 1).count());
    for i in 0..largest {
        for j in 0..largest {
            lap.set_cost(i, j, O3_DUMMY_COST);
        }
    }
    let mut map_idx = 0usize;
    let mut y = 0usize;
    for (i, &rz) in ref_nums.iter().enumerate() {
        if rz == 1 {
            continue;
        }
        let mut x = 0usize;
        for (j, &pz) in prb_nums.iter().enumerate() {
            if pz == 1 {
                continue;
            }
            let constrained = constraints.is_some_and(|c| {
                map_idx < c.len() && j == c.get(map_idx).prb_idx && i == c.get(map_idx).ref_idx
            });
            if constrained {
                lap.set_cost(y, x, O3_LARGE_NEGATIVE_WEIGHT);
                map_idx += 1;
            } else {
                let mut h_sum = 0.0;
                for k in 0..O3_MAX_H_BINS {
                    let rhyk = ref_hist.get(y, k);
                    let phxk = prb_hist.get(x, k);
                    if rhyk == 0 && phxk == 0 {
                        continue;
                    }
                    let d = f64::from(rhyk - phxk);
                    h_sum += d * d / f64::from(rhyk + phxk);
                }
                lap.set_cost(y, x, scheme.cost(j, i, h_sum));
            }
            x += 1;
        }
        y += 1;
    }
}

/// Everything one pairing run reads that does not change between runs.
struct Setup<'a> {
    ref_mol: &'a ROMol,
    prb_nums: Vec<u32>,
    ref_nums: Vec<u32>,
    prb_hvy: Vec<bool>,
    ref_hvy: Vec<bool>,
    constraints: Option<O3AConstraintVect>,
    prb_hist: Option<MolHistogram>,
    ref_hist: Option<MolHistogram>,
}

/// The core O3A constructor: one pairing run for a fixed scheme.
/// `work` is the shared scratch copy of the probe molecule, whose conformer
/// is moved around (and keeps moving between runs, as upstream).
fn run_pairing(
    scheme: &dyn Scheme,
    setup: &Setup<'_>,
    p: &O3AParams<'_>,
    work: &mut ROMol,
) -> Result<(MatchVect, Vec<f64>, f64)> {
    let accuracy = p.options & O3_ACCURACY_MASK;
    let local = p.options & O3_LOCAL_ONLY != 0;
    let constraints = setup.constraints.as_ref();
    let ref_pos = positions(setup.ref_mol, p.ref_cid)?;
    let n_prb_heavy = setup.prb_hvy.iter().filter(|&&h| h).count();
    let n_ref_heavy = setup.ref_hvy.iter().filter(|&&h| h).count();
    let largest = n_prb_heavy.max(n_ref_heavy);

    let mut start = Sdm::default();
    {
        let prb_pos = positions(work, p.prb_cid)?;
        let ctx = SdmContext {
            ref_pos: &ref_pos,
            prb_pos: &prb_pos,
            ref_atomic_nums: &setup.ref_nums,
            prb_atomic_nums: &setup.prb_nums,
            constraints,
        };
        if local {
            start.fill_from_dist(O3_SDM_THRESHOLD_START, &setup.ref_hvy, &setup.prb_hvy, &ctx);
        } else {
            let mut lap = Lap::new(largest);
            compute_cost_matrix(
                &mut lap,
                &setup.prb_nums,
                setup.prb_hist.as_ref().expect("histogram"),
                &setup.ref_nums,
                setup.ref_hist.as_ref().expect("histogram"),
                constraints,
                scheme,
            );
            lap.compute_min_cost_path(largest);
            start.fill_from_lap(&lap, &ctx);
        }
    }

    let mut pairs = [0usize; 4];
    let mut score = [0.0f64; 3];
    let mut pairs_rmsd = [0.0f64; 2];
    let mut best: [Option<Sdm>; 3] = [None, None, None];
    for i in 3..start.len() {
        let (start_matches, start_weights) = start.match_weights(i, Some(scheme), constraints);
        align_mol(
            work,
            setup.ref_mol,
            p.prb_cid,
            p.ref_cid,
            Some(&start_matches),
            Some(&start_weights),
            p.reflect,
            p.max_iters,
        )?;
        pairs[1] = 0;
        score[1] = 0.0;
        let first_it = if accuracy < 1 {
            0
        } else {
            O3_MAX_SDM_THRESHOLD_ITER - 1
        };
        for sdm_threshold_it in first_it..O3_MAX_SDM_THRESHOLD_ITER {
            pairs[2] = 0;
            let mut flag = true;
            let mut iter = 0u32;
            let threshold =
                O3_SDM_THRESHOLD_START + f64::from(sdm_threshold_it) * O3_SDM_THRESHOLD_STEP;
            while flag && iter < O3_MAX_SDM_ITERATIONS {
                let mut progress = Sdm::default();
                {
                    let prb_pos = positions(work, p.prb_cid)?;
                    let ctx = SdmContext {
                        ref_pos: &ref_pos,
                        prb_pos: &prb_pos,
                        ref_atomic_nums: &setup.ref_nums,
                        prb_atomic_nums: &setup.prb_nums,
                        constraints,
                    };
                    progress.fill_from_dist(threshold, &setup.ref_hvy, &setup.prb_hvy, &ctx);
                }
                pairs[3] = progress.len();
                if pairs[3] < 3 {
                    break;
                }
                let (matches, weights) =
                    progress.match_weights(pairs[3], Some(scheme), constraints);
                let (rmsd, trans) = get_alignment_transform(
                    work,
                    setup.ref_mol,
                    p.prb_cid,
                    p.ref_cid,
                    Some(&matches),
                    Some(&weights),
                    p.reflect,
                    p.max_iters,
                )?;
                pairs_rmsd[1] = rmsd;
                flag = pairs[3] > pairs[2]
                    || (pairs[3] == pairs[2]
                        && (iter == 0 || (pairs_rmsd[0] - pairs_rmsd[1]) > O3_RMSD_THRESHOLD));
                if flag {
                    pairs[2] = pairs[3];
                    pairs_rmsd[0] = pairs_rmsd[1];
                    best[2] = Some(progress);
                    transform_conformer(work.conformer_mut(p.prb_cid)?, &trans);
                }
                iter += 1;
            }
            score[2] = if pairs[2] >= 3 {
                let prb_pos = positions(work, p.prb_cid)?;
                best[2]
                    .as_ref()
                    .expect("best SDM")
                    .score_alignment(scheme, &ref_pos, &prb_pos)
            } else {
                0.0
            };
            if score[2] - score[1] > O3_SCORE_THRESHOLD {
                pairs[1] = pairs[2];
                score[1] = score[2];
                best[1] = best[2].clone();
            }
        }
        if score[1] - score[0] > O3_SCORE_THRESHOLD {
            pairs[0] = pairs[1];
            score[0] = score[1];
            best[0] = best[1].clone();
        }
    }
    if pairs[0] < 3 && start.len() >= 3 {
        pairs[0] = start.len();
        let prb_pos = positions(work, p.prb_cid)?;
        score[0] = start.score_alignment(scheme, &ref_pos, &prb_pos);
        best[0] = Some(start);
    }
    if pairs[0] >= 3 {
        let (matches, weights) =
            best[0]
                .as_ref()
                .expect("best SDM")
                .match_weights(pairs[0], Some(scheme), constraints);
        Ok((matches, weights, score[0]))
    } else {
        Ok((vec![(0, 0); pairs[0]], vec![0.0; pairs[0]], 0.0))
    }
}

impl O3A {
    /// `O3A::O3A(prbMol, refMol, ...)`: the best pairing over the
    /// weight/coefficient grid selected by `options`.
    pub fn new(
        prb_mol: &ROMol,
        ref_mol: &ROMol,
        props: &O3AProps<'_>,
        p: &O3AParams<'_>,
    ) -> Result<O3A> {
        let prb_nums = atomic_nums(prb_mol);
        let ref_nums = atomic_nums(ref_mol);
        let prb_hvy = heavy_flags(&prb_nums);
        let ref_hvy = heavy_flags(&ref_nums);
        let accuracy = p.options & O3_ACCURACY_MASK;
        let local = p.options & O3_LOCAL_ONLY != 0;
        let mut constraints = O3AConstraintVect::default();
        if let Some(map) = p.constraint_map {
            if let Some(w) = p.constraint_weights
                && w.len() != map.len()
            {
                return Err(align_error(
                    "The number of weights should match the number of constraints",
                ));
            }
            for (i, &(prb_idx, ref_idx)) in map.iter().enumerate() {
                if prb_idx < 0
                    || prb_idx as usize >= prb_nums.len()
                    || ref_idx < 0
                    || ref_idx as usize >= ref_nums.len()
                {
                    return Err(align_error("Constrained atom idx out of range"));
                }
                if !prb_hvy[prb_idx as usize] || !ref_hvy[ref_idx as usize] {
                    return Err(align_error("Constrained atoms must be heavy atoms"));
                }
                constraints.append(
                    prb_idx as usize,
                    ref_idx as usize,
                    p.constraint_weights
                        .map_or(O3_DEFAULT_CONSTRAINT_WEIGHT, |w| w[i]),
                );
            }
        }
        let (prb_hist, ref_hist) = if local {
            (None, None)
        } else {
            (
                Some(MolHistogram::new(
                    &prb_nums,
                    &get_3d_distance_mat(prb_mol, p.prb_cid)?,
                )),
                Some(MolHistogram::new(
                    &ref_nums,
                    &get_3d_distance_mat(ref_mol, p.ref_cid)?,
                )),
            )
        };
        let setup = Setup {
            ref_mol,
            prb_nums,
            ref_nums,
            prb_hvy,
            ref_hvy,
            constraints: Some(constraints),
            prb_hist,
            ref_hist,
        };
        let mut work = prb_mol.clone();
        let is_mmff = matches!(props, O3AProps::Mmff { .. });
        let mut best: Option<(MatchVect, Vec<f64>, f64)> = None;
        let mut score0 = 0.0f64;
        let l_max = if accuracy < 2 && is_mmff { 1 } else { 0 };
        for l in 0..=l_max {
            let step = if accuracy < 3 {
                1
            } else {
                O3_MAX_WEIGHT_COEFF + 1
            };
            let mut c = 0;
            while c <= O3_MAX_WEIGHT_COEFF {
                let result = match props {
                    O3AProps::Mmff { prb, reference } => {
                        let prb_types: Vec<u8> =
                            (0..prb_mol.num_atoms()).map(|i| prb.atom_type(i)).collect();
                        let ref_types: Vec<u8> = (0..ref_mol.num_atoms())
                            .map(|i| reference.atom_type(i))
                            .collect();
                        let prb_charges: Vec<f64> = (0..prb_mol.num_atoms())
                            .map(|i| prb.partial_charge(i))
                            .collect();
                        let ref_charges: Vec<f64> = (0..ref_mol.num_atoms())
                            .map(|i| reference.partial_charge(i))
                            .collect();
                        let scheme = MmffScheme {
                            ref_types: &ref_types,
                            prb_types: &prb_types,
                            ref_charges: &ref_charges,
                            prb_charges: &prb_charges,
                            sim: mmff_similarity(),
                            coeff: c,
                            weight: if l == 1 { c } else { 0 },
                            use_mmff_sim: l == 1,
                        };
                        run_pairing(&scheme, &setup, p, &mut work)?
                    }
                    O3AProps::Crippen { prb, reference } => {
                        let scheme = CrippenScheme {
                            ref_logp: reference,
                            prb_logp: prb,
                            coeff: c,
                        };
                        run_pairing(&scheme, &setup, p, &mut work)?
                    }
                };
                let score1 = result.2;
                if score1 - score0 > O3_SCORE_THRESHOLD || score0.abs() < 1.0e-10 {
                    best = Some(result);
                    score0 = score1;
                }
                c += step;
            }
        }
        let (matches, weights, score) = best.expect("at least one pairing");
        let score = if matches.len() >= 3 { score } else { 0.0 };
        Ok(O3A {
            prb_cid: p.prb_cid,
            ref_cid: p.ref_cid,
            reflect: p.reflect,
            max_iters: p.max_iters,
            matches,
            weights,
            score,
        })
    }

    fn rmsd(
        prb_mol: &ROMol,
        ref_mol: &ROMol,
        prb_pos: &[Point3D],
        ref_pos: &[Point3D],
        matches: Option<&MatchVect>,
    ) -> f64 {
        match matches {
            Some(m) => {
                let sum: f64 = m
                    .iter()
                    .map(|&(p, r)| (prb_pos[p as usize] - ref_pos[r as usize]).length_sq())
                    .sum();
                (sum / m.len() as f64).sqrt()
            }
            None => {
                let centroid = |mol: &ROMol, pos: &[Point3D]| {
                    let (mut c, mut n) = (Point3D::new(0.0, 0.0, 0.0), 0usize);
                    for (i, p) in pos.iter().enumerate() {
                        if mol.atom(i as u32).atomic_num() == 1 {
                            continue;
                        }
                        c.x += p.x;
                        c.y += p.y;
                        c.z += p.z;
                        n += 1;
                    }
                    c.x /= n as f64;
                    c.y /= n as f64;
                    c.z /= n as f64;
                    c
                };
                (centroid(ref_mol, ref_pos) - centroid(prb_mol, prb_pos)).length()
            }
        }
    }

    /// `O3A::align`: moves the probe conformer and returns the RMSD of the
    /// paired atoms.
    pub fn align(&self, prb_mol: &mut ROMol, ref_mol: &ROMol) -> Result<f64> {
        let paired = self.matches.len() >= 3;
        if paired {
            align_mol(
                prb_mol,
                ref_mol,
                self.prb_cid,
                self.ref_cid,
                Some(&self.matches),
                Some(&self.weights),
                self.reflect,
                self.max_iters,
            )?;
        }
        Ok(Self::rmsd(
            prb_mol,
            ref_mol,
            prb_mol.conformer(self.prb_cid)?.positions(),
            ref_mol.conformer(self.ref_cid)?.positions(),
            paired.then_some(&self.matches),
        ))
    }

    /// `O3A::trans`: the RMSD and the transformation that [`Self::align`]
    /// would apply, without moving anything.
    pub fn trans(&self, prb_mol: &ROMol, ref_mol: &ROMol) -> Result<(f64, Transform3D)> {
        let mut conf: Conformer = prb_mol.conformer(self.prb_cid)?.clone();
        let mut trans = Transform3D::identity();
        if self.matches.len() >= 3 {
            trans = get_alignment_transform(
                prb_mol,
                ref_mol,
                self.prb_cid,
                self.ref_cid,
                Some(&self.matches),
                Some(&self.weights),
                self.reflect,
                self.max_iters,
            )?
            .1;
            transform_conformer(&mut conf, &trans);
        }
        let rmsd = Self::rmsd(
            prb_mol,
            ref_mol,
            conf.positions(),
            ref_mol.conformer(self.ref_cid)?.positions(),
            Some(&self.matches),
        );
        Ok((rmsd, trans))
    }
}

/// `getO3AForProbeConfs`: one [`O3A`] per probe conformer, in conformer order.
pub fn get_o3a_for_probe_confs(
    prb_mol: &ROMol,
    ref_mol: &ROMol,
    props: &O3AProps<'_>,
    num_threads: i32,
    p: &O3AParams<'_>,
) -> Result<Vec<O3A>> {
    let ids: Vec<i32> = prb_mol.conformers().iter().map(|c| c.id() as i32).collect();
    par_map(&ids, num_threads, |&cid| {
        let params = O3AParams {
            prb_cid: cid,
            ..p.clone()
        };
        O3A::new(prb_mol, ref_mol, props, &params)
    })
    .into_iter()
    .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn similarity_matrix_is_symmetric_enough_to_load() {
        let sim = mmff_similarity();
        assert_eq!(sim[0][0], 1);
        assert_eq!(sim[98][98], 1);
    }
}
