//! Constraints, the atom-type dependent scoring schemes and the SDM matrix of
//! upstream `O3AAlignMolecules.{h,cpp}`.

use super::lap::{Lap, O3_DUMMY_COST};
use crate::geometry::Point3D;

pub const O3_DEFAULT_CONSTRAINT_WEIGHT: f64 = 100.0;
pub const O3_THRESHOLD_DIFF_DISTANCE: f64 = 0.1;
pub const O3_CHARGE_WEIGHT: f64 = 10.0;
pub const O3_SCORING_FUNCTION_ALPHA: f64 = 5.0;
pub const O3_SCORING_FUNCTION_BETA: f64 = 0.5;
pub const O3_CHARGE_COEFF: f64 = 5.0;
pub const O3_CRIPPEN_COEFF: f64 = 1.0;
pub const O3_MAX_WEIGHT_COEFF: i32 = 5;

/// A forced `(probe, reference)` atom pair with its weight.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct O3AConstraint {
    idx: usize,
    pub prb_idx: usize,
    pub ref_idx: usize,
    pub weight: f64,
}

/// `O3AConstraintVect`: constraints sorted by probe index, reference index,
/// decreasing weight and insertion order.
#[derive(Debug, Clone, Default)]
pub struct O3AConstraintVect {
    count: usize,
    items: Vec<O3AConstraint>,
}

impl O3AConstraintVect {
    pub fn append(&mut self, prb_idx: usize, ref_idx: usize, weight: f64) {
        self.items.push(O3AConstraint {
            idx: self.count,
            prb_idx,
            ref_idx,
            weight,
        });
        self.items.sort_by(|a, b| {
            a.prb_idx
                .cmp(&b.prb_idx)
                .then(a.ref_idx.cmp(&b.ref_idx))
                .then_with(|| {
                    if a.weight != b.weight {
                        b.weight.partial_cmp(&a.weight).unwrap()
                    } else {
                        a.idx.cmp(&b.idx)
                    }
                })
        });
        self.count += 1;
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    pub fn get(&self, i: usize) -> &O3AConstraint {
        &self.items[i]
    }
}

/// The atom property dependent cost, weight and scoring functions (the
/// MMFF94 and Crippen variants of upstream's `o3a*Func` callbacks).
pub trait Scheme: Sync {
    /// Cost of pairing probe atom `prb` with reference atom `rf`.
    fn cost(&self, prb: usize, rf: usize, h_sum: f64) -> i32;
    /// Weight of the pair in the final weighted alignment.
    fn weight(&self, prb: usize, rf: usize) -> f64;
    /// Contribution of the pair to the alignment score.
    fn score(&self, prb: usize, rf: usize, ref_pos: &[Point3D], prb_pos: &[Point3D]) -> f64;
}

/// MMFF94 atom types and partial charges of both molecules.
#[derive(Debug)]
pub struct MmffScheme<'a> {
    pub ref_types: &'a [u8],
    pub prb_types: &'a [u8],
    pub ref_charges: &'a [f64],
    pub prb_charges: &'a [f64],
    pub sim: &'a [[u8; 99]; 99],
    pub coeff: i32,
    pub weight: i32,
    pub use_mmff_sim: bool,
}

impl MmffScheme<'_> {
    fn similarity(&self, prb: usize, rf: usize) -> u8 {
        self.sim[self.ref_types[rf] as usize - 1][self.prb_types[prb] as usize - 1]
    }
}

impl Scheme for MmffScheme<'_> {
    fn cost(&self, prb: usize, rf: usize, h_sum: f64) -> i32 {
        let sim = self.similarity(prb, rf);
        let charge = f64::from(self.coeff)
            * O3_CHARGE_WEIGHT
            * (self.ref_charges[rf] - self.prb_charges[prb]).abs();
        let mmff = if self.use_mmff_sim {
            f64::from(O3_MAX_WEIGHT_COEFF - self.coeff) * f64::from(sim)
        } else {
            0.0
        };
        ((charge + mmff + h_sum) * 1.0e03).round() as i32
    }

    fn weight(&self, prb: usize, rf: usize) -> f64 {
        let sim = self.similarity(prb, rf);
        let (rc, pc) = (self.ref_charges[rf], self.prb_charges[prb]);
        f64::from(O3_MAX_WEIGHT_COEFF - self.weight)
            * ((1.0 + O3_CHARGE_COEFF * (rc + pc).abs()) / (1.0 + (rc - pc).abs()))
            + f64::from(self.weight) / f64::from(sim)
    }

    fn score(&self, prb: usize, rf: usize, ref_pos: &[Point3D], prb_pos: &[Point3D]) -> f64 {
        let sq_dist = (ref_pos[rf] - prb_pos[prb]).length_sq();
        let (rc, pc) = (self.ref_charges[rf], self.prb_charges[prb]);
        let charge_diff = (rc - pc).abs();
        let charge_sum = (rc + pc).abs();
        (O3_SCORING_FUNCTION_ALPHA + (1.0 + O3_CHARGE_COEFF * charge_sum) / (1.0 + charge_diff))
            * (-O3_SCORING_FUNCTION_BETA * sq_dist).exp()
    }
}

/// Crippen logP contributions of both molecules.
#[derive(Debug)]
pub struct CrippenScheme<'a> {
    pub ref_logp: &'a [f64],
    pub prb_logp: &'a [f64],
    pub coeff: i32,
}

impl Scheme for CrippenScheme<'_> {
    fn cost(&self, prb: usize, rf: usize, h_sum: f64) -> i32 {
        ((f64::from(self.coeff)
            * O3_CHARGE_WEIGHT
            * (self.ref_logp[rf] - self.prb_logp[prb]).abs()
            + h_sum)
            * 1.0e03)
            .round() as i32
    }

    fn weight(&self, prb: usize, rf: usize) -> f64 {
        let (r, p) = (self.ref_logp[rf], self.prb_logp[prb]);
        0.01 + (r + p).abs() / (1.0 + (r - p).abs())
    }

    fn score(&self, prb: usize, rf: usize, ref_pos: &[Point3D], prb_pos: &[Point3D]) -> f64 {
        let sq_dist = (ref_pos[rf] - prb_pos[prb]).length_sq();
        let (r, p) = (self.ref_logp[rf], self.prb_logp[prb]);
        (O3_SCORING_FUNCTION_ALPHA
            + (1.0 + O3_CRIPPEN_COEFF * (r + p).abs()) / (1.0 + (r - p).abs()))
            * (-O3_SCORING_FUNCTION_BETA * sq_dist).exp()
    }
}

#[derive(Debug, Clone, Copy)]
struct SdmElement {
    /// `[reference, probe]` atom indices
    idx: [usize; 2],
    score: i32,
    cost: i32,
    sq_dist: f64,
    /// index into the constraint vector
    constraint: Option<usize>,
}

/// The sorted matrix of atom equivalences (`SDM`).
#[derive(Debug, Clone, Default)]
pub struct Sdm {
    elements: Vec<SdmElement>,
}

/// Everything `SDM` reads from the two molecules.
#[derive(Debug)]
pub struct SdmContext<'a> {
    pub ref_pos: &'a [Point3D],
    pub prb_pos: &'a [Point3D],
    pub ref_atomic_nums: &'a [u32],
    pub prb_atomic_nums: &'a [u32],
    pub constraints: Option<&'a O3AConstraintVect>,
}

impl Sdm {
    pub fn len(&self) -> usize {
        self.elements.len()
    }

    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// `SDM::fillFromLAP`
    pub fn fill_from_lap(&mut self, lap: &Lap, ctx: &SdmContext<'_>) {
        let n_ref_heavy = ctx.ref_atomic_nums.iter().filter(|&&z| z != 1).count();
        self.elements.clear();
        for i in 0..n_ref_heavy {
            let sol = lap.row_sol(i);
            if lap.cost(i, sol) >= O3_DUMMY_COST {
                continue;
            }
            self.elements.push(SdmElement {
                idx: [i, sol],
                score: 0,
                cost: lap.cost(i, sol),
                sq_dist: 0.0,
                constraint: None,
            });
        }
        let n_equiv = self.elements.len();
        // map heavy atom ranks back to indices in the original molecules
        for (n, nums) in [ctx.ref_atomic_nums, ctx.prb_atomic_nums]
            .into_iter()
            .enumerate()
        {
            for e in &mut self.elements {
                let k = e.idx[n];
                let (mut j, mut n_atoms) = (0usize, 0usize);
                while n_atoms < nums.len() && j <= k {
                    if nums[n_atoms] > 1 {
                        j += 1;
                        e.idx[n] = n_atoms;
                    }
                    n_atoms += 1;
                }
            }
        }
        let mut diff = vec![0.0f64; n_equiv * n_equiv];
        for i in 0..n_equiv {
            self.elements[i].constraint = None;
            if self.elements[i].cost < 0 {
                if let Some(cv) = ctx.constraints {
                    for j in 0..cv.len() {
                        if cv.get(j).prb_idx == self.elements[i].idx[1]
                            && cv.get(j).ref_idx == self.elements[i].idx[0]
                        {
                            self.elements[i].constraint = Some(j);
                            break;
                        }
                    }
                }
            }
            for j in 0..n_equiv {
                if i == j {
                    continue;
                }
                let (a, b) = (&self.elements[i], &self.elements[j]);
                diff[i * n_equiv + j] = ((ctx.ref_pos[a.idx[0]] - ctx.ref_pos[b.idx[0]]).length()
                    - (ctx.prb_pos[a.idx[1]] - ctx.prb_pos[b.idx[1]]).length())
                .abs();
            }
        }
        for i in 0..n_equiv {
            let mut score = 0;
            for j in 0..n_equiv {
                if diff[i * n_equiv + j] > O3_THRESHOLD_DIFF_DISTANCE && self.elements[i].cost >= 0
                {
                    score += 1;
                }
            }
            self.elements[i].score = score;
        }
        self.elements.sort_by(|a, b| {
            a.score
                .cmp(&b.score)
                .then(a.cost.cmp(&b.cost))
                .then(a.idx[0].cmp(&b.idx[0]))
                .then(a.idx[1].cmp(&b.idx[1]))
        });
    }

    /// `SDM::fillFromDist`
    pub fn fill_from_dist(
        &mut self,
        threshold: f64,
        ref_hvy: &[bool],
        prb_hvy: &[bool],
        ctx: &SdmContext<'_>,
    ) {
        self.elements.clear();
        let sq_threshold = threshold * threshold;
        let (ref_n, prb_n) = (ctx.ref_pos.len(), ctx.prb_pos.len());
        let largest = ref_n.max(prb_n);
        let mut ref_used = vec![false; largest];
        let mut prb_used = vec![false; largest];
        let cv = ctx.constraints;
        let mut map_idx = 0usize;
        let is_constraint = |map_idx: usize, i: usize, j: usize| {
            cv.is_some_and(|c| {
                map_idx < c.len() && j == c.get(map_idx).prb_idx && i == c.get(map_idx).ref_idx
            })
        };
        for i in 0..ref_n {
            if !ref_hvy[i] {
                continue;
            }
            for j in 0..prb_n {
                if !prb_hvy[j] {
                    continue;
                }
                let sq_dist = (ctx.ref_pos[i] - ctx.prb_pos[j]).length_sq();
                let constrained = is_constraint(map_idx, i, j);
                if sq_dist < sq_threshold || constrained {
                    let mut constraint = None;
                    if constrained {
                        constraint = Some(map_idx);
                        // duplicates of the same pair: keep the heaviest
                        while is_constraint(map_idx, i, j) {
                            map_idx += 1;
                        }
                    }
                    self.elements.push(SdmElement {
                        idx: [i, j],
                        score: 0,
                        cost: 0,
                        sq_dist,
                        constraint,
                    });
                }
            }
        }
        let weight_of = |e: &SdmElement| {
            e.constraint
                .map(|c| cv.expect("constraints").get(c).weight)
                .unwrap_or(0.0)
        };
        self.elements.sort_by(|a, b| {
            let (aw, bw) = (weight_of(a), weight_of(b));
            if aw != bw {
                bw.partial_cmp(&aw).unwrap()
            } else if a.sq_dist != b.sq_dist {
                a.sq_dist.partial_cmp(&b.sq_dist).unwrap()
            } else {
                a.idx[0].cmp(&b.idx[0]).then(a.idx[1].cmp(&b.idx[1]))
            }
        });
        let n = self.elements.len();
        let mut pairs = 0;
        let mut i = 0;
        while i < n && !ref_used[self.elements[i].idx[0]] && !prb_used[self.elements[i].idx[1]] {
            pairs += 1;
            ref_used[self.elements[i].idx[0]] = true;
            prb_used[self.elements[i].idx[1]] = true;
            i += 1;
        }
        self.elements.truncate(pairs);
    }

    /// `SDM::scoreAlignment`
    pub fn score_alignment(
        &self,
        scheme: &dyn Scheme,
        ref_pos: &[Point3D],
        prb_pos: &[Point3D],
    ) -> f64 {
        self.elements
            .iter()
            .map(|e| scheme.score(e.idx[1], e.idx[0], ref_pos, prb_pos))
            .sum()
    }

    /// `SDM::prepareMatchWeightsVect` for the first `n` elements: the
    /// `(probe, reference)` pairs and their weights normalised by the
    /// smallest one.
    pub fn match_weights(
        &self,
        n: usize,
        scheme: Option<&dyn Scheme>,
        constraints: Option<&O3AConstraintVect>,
    ) -> (Vec<(i32, i32)>, Vec<f64>) {
        let mut matches = Vec::with_capacity(n);
        let mut weights = Vec::with_capacity(n);
        let mut min = f64::MAX;
        for (i, e) in self.elements.iter().take(n).enumerate() {
            matches.push((e.idx[1] as i32, e.idx[0] as i32));
            let mut w = scheme.map_or(1.0, |s| s.weight(e.idx[1], e.idx[0]));
            if let (Some(c), Some(cv)) = (e.constraint, constraints) {
                w += cv.get(c).weight;
            }
            if i == 0 || w < min {
                min = w;
            }
            weights.push(w);
        }
        for w in &mut weights {
            *w /= min;
        }
        (matches, weights)
    }
}
