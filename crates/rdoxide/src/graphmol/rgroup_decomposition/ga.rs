//! Genetic algorithm search for the best assignment of the R groups (upstream
//! `RGroupDecomposition/RGroupGa.{h,cpp}`).
//!
//! The runs are executed one after the other, also for `gaParallelRuns`:
//! upstream's parallel runs share one random number generator, so their
//! results depend on thread scheduling.

use std::collections::BTreeSet;
use std::time::Instant;

use super::data::RGroupMatch;
use super::fingerprint_score::{
    FingerprintVarianceScoreData, VarianceDataForLabel, fingerprint_variance_score,
};
use super::params::{RGroupDecompositionParameters, rgroup_score};
use super::score::{RGroupScorer, check_for_timeout};
use crate::Result;
use crate::external::ga::{
    DEFAULT_SEED, GaOperation, GaProblem, IntegerStringPolicy, LinkedPopLinearSel, RandomUtil,
    initialize_string, mutate_string, two_point_crossover,
};
use crate::rd_log_info;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OperationName {
    RgroupMutate,
    Crossover,
    Create,
}

/// `RGroupDecompositionChromosome`
pub struct Chromosome {
    string: Vec<i32>,
    fitness: f64,
    fingerprint_variance_score_data: FingerprintVarianceScoreData,
    operation_name: OperationName,
    permutation: Vec<usize>,
}

/// `GaResult`: the best score and the permutations that reach it.
pub struct GaResult {
    pub scorer: RGroupScorer,
}

/// `RGroupGa`
pub struct RGroupGa<'a> {
    matches: &'a [Vec<RGroupMatch>],
    labels: &'a BTreeSet<i32>,
    params: &'a RGroupDecompositionParameters,
    /// `score` mutates the scorer's working state, which upstream shares with
    /// the decomposition; here the GA gets its own
    scorer: RGroupScorer,
    policy: IntegerStringPolicy,
    popsize: usize,
    number_operations: usize,
    number_operations_without_improvement: usize,
    chrom_length: usize,
    number_permutations: u32,
    t0: Option<Instant>,
}

fn copy_variance_data(from: &FingerprintVarianceScoreData, to: &mut FingerprintVarianceScoreData) {
    for (label, data) in &from.labels_to_variance_data {
        match to.labels_to_variance_data.get_mut(label) {
            None => {
                to.labels_to_variance_data.insert(*label, data.clone());
            }
            Some(existing) => {
                existing.number_fingerprints = data.number_fingerprints;
                existing.bit_counts.clone_from(&data.bit_counts);
            }
        }
    }
    to.number_of_molecules = from.number_of_molecules;
    to.number_of_missing_user_rgroups = from.number_of_missing_user_rgroups;
}

fn clear_variance_data(data: &mut FingerprintVarianceScoreData) {
    for d in data.labels_to_variance_data.values_mut() {
        let VarianceDataForLabel {
            number_fingerprints,
            bit_counts,
        } = d;
        *number_fingerprints = 0;
        bit_counts.iter_mut().for_each(|c| *c = Default::default());
    }
    data.number_of_missing_user_rgroups = 0;
    data.number_of_molecules = 0;
}

impl Chromosome {
    /// `decode`: the match picked for every molecule
    fn decode(&mut self, matches: &[Vec<RGroupMatch>]) {
        self.permutation.clear();
        let mut pos = 0;
        for m in matches {
            if m.len() == 1 {
                self.permutation.push(0);
            } else {
                self.permutation.push(self.string[pos] as usize);
                pos += 1;
            }
        }
    }

    /// `copyGene`
    fn copy_gene(&mut self, other: &Chromosome) {
        self.string.copy_from_slice(&other.string);
        copy_variance_data(
            &other.fingerprint_variance_score_data,
            &mut self.fingerprint_variance_score_data,
        );
    }
}

impl<'a> RGroupGa<'a> {
    pub fn new(
        matches: &'a [Vec<RGroupMatch>],
        labels: &'a BTreeSet<i32>,
        params: &'a RGroupDecompositionParameters,
        t0: Option<Instant>,
    ) -> Self {
        let mut policy = IntegerStringPolicy::new(matches.len());
        let mut number_permutations = 1u32;
        let mut pos = 0;
        for m in matches {
            if m.len() == 1 {
                continue;
            }
            policy.set_max(pos, m.len() as i32);
            let count = u64::from(number_permutations) * m.len() as u64;
            number_permutations = count.min(u64::from(u32::MAX)) as u32;
            pos += 1;
        }
        let chrom_length = pos;
        let mut popsize = 100 + chrom_length / 10;
        if popsize > 200 {
            popsize = 200;
        }
        if params.ga_population_size > 0 {
            popsize = params.ga_population_size as usize;
        }
        // run for a long time and stop early when the score does not improve
        let number_operations = if params.ga_maximum_operations > 0 {
            params.ga_maximum_operations as usize
        } else {
            1_000_000
        };
        let number_operations_without_improvement =
            if params.ga_number_operations_without_improvement > 0 {
                params.ga_number_operations_without_improvement as usize
            } else {
                7500
            };
        let seed = {
            let mut rng = RandomUtil::instance();
            if params.ga_random_seed >= 0 {
                rng.seed(params.ga_random_seed as u32);
                params.ga_random_seed as u32
            } else if params.ga_random_seed == -2 {
                let seed = std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .map_or(0, |d| d.as_nanos() as u32);
                rng.seed(seed);
                seed
            } else {
                DEFAULT_SEED
            }
        };
        rd_log_info!("GA RNG seed {seed}\n");
        RGroupGa {
            matches,
            labels,
            params,
            scorer: RGroupScorer::default(),
            policy,
            popsize,
            number_operations,
            number_operations_without_improvement,
            chrom_length,
            number_permutations,
            t0,
        }
    }

    pub fn number_permutations(&self) -> u32 {
        self.number_permutations
    }

    pub fn population_size(&self) -> usize {
        self.popsize
    }

    /// `run`: one GA run
    pub fn run(&mut self, run_number: usize) -> Result<GaResult> {
        let start = Instant::now();
        let time_info = |start: Instant| format!("Time {:7.2}", start.elapsed().as_secs_f64());
        let mut rng = RandomUtil::instance();
        let mut population = LinkedPopLinearSel::new(self)?;
        rd_log_info!(
            "Running GA run {:2} number operations {:5} population size {:5} number operations without improvement {:5} chromosome length {:5} {}\n",
            run_number,
            self.number_operations,
            self.popsize,
            self.number_operations_without_improvement,
            self.chrom_length,
            time_info(start)
        );
        population.create(self, &mut rng)?;
        let mut best_score = population.best_score();
        rd_log_info!("{}\n", population.info());
        let mut n_ops = 0;
        let mut last_improvement_op = 0;
        while n_ops < self.number_operations {
            population.iterate(self, &mut rng)?;
            n_ops += 1;
            if n_ops % 1000 == 0 {
                rd_log_info!(
                    "Run {} {} {}\n",
                    run_number,
                    population.info(),
                    time_info(start)
                );
            }
            if population.best_score() > best_score {
                best_score = population.best_score();
                last_improvement_op = n_ops;
                rd_log_info!(
                    "Run {:2} OP {:5} Fit {:7.3} {}\n",
                    run_number,
                    n_ops,
                    best_score,
                    time_info(start)
                );
            }
            if n_ops - last_improvement_op > self.number_operations_without_improvement {
                rd_log_info!(
                    "Run {run_number} Op {n_ops} No improvement since {last_improvement_op} finishing..\n"
                );
                break;
            }
            if let Some(t0) = self.t0 {
                check_for_timeout(t0, self.params.timeout)?;
            }
        }
        let permutations: Vec<Vec<usize>> = population
            .tied_best(1e-6)
            .iter()
            .map(|c| c.permutation.clone())
            .collect();
        rd_log_info!("Run {run_number} Execution {}\n", time_info(start));
        Ok(GaResult {
            scorer: RGroupScorer::from_permutations(&permutations, population.best().fitness),
        })
    }

    /// `runBatch`: `gaNumberRuns` runs
    pub fn run_batch(&mut self) -> Result<Vec<GaResult>> {
        (0..self.params.ga_number_runs.max(0) as usize)
            .map(|n| self.run(n + 1))
            .collect()
    }

    /// `RGroupDecompData::score` for a chromosome
    fn data_score(
        &mut self,
        permutation: &[usize],
        data: &mut FingerprintVarianceScoreData,
    ) -> Result<f64> {
        match self.params.score_method {
            rgroup_score::MATCH => self
                .scorer
                .match_score(permutation, self.matches, self.labels),
            rgroup_score::FINGERPRINT_VARIANCE => {
                fingerprint_variance_score(permutation, self.matches, self.labels, Some(data))
            }
            _ => Ok(f64::NAN),
        }
    }
}

/// `rGroupMutateOperation`
fn mutate_operation(
    ga: &mut RGroupGa,
    rng: &mut RandomUtil,
    parents: &[&Chromosome],
    children: &mut [Chromosome],
) {
    let parent = parents[0];
    let child = &mut children[0];
    child.copy_gene(parent);
    mutate_string(&ga.policy, rng, &mut child.string, -1.0);
    child.operation_name = OperationName::RgroupMutate;
    child.decode(ga.matches);
    if child
        .fingerprint_variance_score_data
        .labels_to_variance_data
        .is_empty()
    {
        return;
    }
    // update the variance data incrementally for the changed molecules
    for pos in 0..parent.permutation.len() {
        let parent_value = parent.permutation[pos];
        let child_value = child.permutation[pos];
        if parent_value != child_value {
            let data = &mut child.fingerprint_variance_score_data;
            // the fingerprints of matched R groups were computed when the
            // parent was scored, so these cannot fail
            let _ = data.remove_variance_data(pos, parent_value, ga.matches, ga.labels);
            let _ = data.add_variance_data(pos, child_value, ga.matches, ga.labels);
        }
    }
}

/// `rGroupCrossoverOperation`
fn crossover_operation(
    _ga: &mut RGroupGa,
    rng: &mut RandomUtil,
    parents: &[&Chromosome],
    children: &mut [Chromosome],
) {
    let (child1, rest) = children.split_first_mut().expect("two children");
    let child2 = &mut rest[0];
    child1.operation_name = OperationName::Crossover;
    child2.operation_name = OperationName::Crossover;
    clear_variance_data(&mut child1.fingerprint_variance_score_data);
    clear_variance_data(&mut child2.fingerprint_variance_score_data);
    two_point_crossover(
        rng,
        &parents[0].string,
        &parents[1].string,
        &mut child1.string,
        &mut child2.string,
    );
}

impl GaProblem for RGroupGa<'_> {
    type Chromosome = Chromosome;

    fn population_size(&self) -> usize {
        self.popsize
    }

    fn selection_pressure(&self) -> f64 {
        1.0001
    }

    fn operations(&self) -> Vec<GaOperation<Self>> {
        // biased to mutation as that operator is so efficient
        vec![
            GaOperation {
                n_parents: 1,
                n_children: 1,
                weight: 75.0,
                function: mutate_operation,
            },
            GaOperation {
                n_parents: 2,
                n_children: 2,
                weight: 25.0,
                function: crossover_operation,
            },
        ]
    }

    fn create_chromosome(&self) -> Chromosome {
        Chromosome {
            string: vec![0; self.chrom_length],
            fitness: 0.0,
            fingerprint_variance_score_data: FingerprintVarianceScoreData::default(),
            operation_name: OperationName::Create,
            permutation: Vec::with_capacity(self.matches.len()),
        }
    }

    fn initialize(&mut self, rng: &mut RandomUtil, c: &mut Chromosome) {
        initialize_string(&self.policy, rng, &mut c.string);
    }

    /// `RGroupDecompositionChromosome::score`
    fn score(&mut self, c: &mut Chromosome) -> Result<f64> {
        if c.operation_name != OperationName::RgroupMutate {
            c.decode(self.matches);
        }
        if self.params.score_method == rgroup_score::FINGERPRINT_VARIANCE
            && !c
                .fingerprint_variance_score_data
                .labels_to_variance_data
                .is_empty()
            && c.operation_name == OperationName::RgroupMutate
        {
            c.fitness = c
                .fingerprint_variance_score_data
                .fingerprint_variance_group_score()?;
        } else {
            let permutation = c.permutation.clone();
            c.fitness = self.data_score(&permutation, &mut c.fingerprint_variance_score_data)?;
        }
        Ok(c.fitness)
    }

    fn fitness(&self, c: &Chromosome) -> f64 {
        c.fitness
    }

    fn is_ok(&self, _c: &Chromosome) -> bool {
        true
    }

    fn equals(&self, a: &Chromosome, b: &Chromosome) -> bool {
        a.string == b.string
    }
}
