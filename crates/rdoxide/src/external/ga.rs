//! Genetic algorithm toolkit (upstream `External/GA`, Copyright (C) 2020
//! Gareth Jones, Glysade LLC, BSD license of the RDKit): the shared random
//! number generator, integer string chromosomes and the steady state,
//! no-duplicates GA with linear-normalised parent selection used by
//! `RGroupDecomposition`.
//!
//! Only the parts RDKit instantiates are ported: `IntegerStringChromosome`
//! (its policy never allows nulls or switching) with `mutate` and
//! `twoPointCrossover`, and `LinkedPopLinearSel`.

use std::sync::{LazyLock, Mutex, MutexGuard};

use crate::Result;
use crate::rdgeneral::random::Mt19937;

/// `std::mt19937`'s default seed
pub const DEFAULT_SEED: u32 = 5489;

/// `GarethUtil::RandomUtil`
#[derive(Debug, Clone)]
pub struct RandomUtil {
    rng: Mt19937,
}

static INSTANCE: LazyLock<Mutex<RandomUtil>> = LazyLock::new(|| {
    Mutex::new(RandomUtil {
        rng: Mt19937::new(DEFAULT_SEED),
    })
});

impl RandomUtil {
    /// The process-wide generator (`RandomUtil::getInstance`); like upstream,
    /// its state carries over from one GA to the next unless it is re-seeded.
    pub fn instance() -> MutexGuard<'static, RandomUtil> {
        INSTANCE.lock().unwrap_or_else(|e| e.into_inner())
    }

    pub fn seed(&mut self, seed: u32) {
        self.rng.seed(seed);
    }

    /// `normalRand`: `std::uniform_real_distribution<double>(0, 1)` over
    /// `std::generate_canonical<double, 53>` (two 32-bit draws, low word first).
    pub fn normal_rand(&mut self) -> f64 {
        let low = f64::from(self.rng.next_u32());
        let high = f64::from(self.rng.next_u32());
        let r = (low + high * 4_294_967_296.0) / 18_446_744_073_709_551_616.0;
        if r >= 1.0 {
            1.0 - f64::EPSILON / 2.0
        } else {
            r
        }
    }

    /// `randomInt(bottom, top)`: `bottom + rng() % (top - bottom)`.
    pub fn random_int(&mut self, bottom: i32, top: i32) -> i32 {
        let step = self.rng.next_u32() % (top - bottom) as u32;
        bottom + step as i32
    }

    pub fn random_boolean(&mut self) -> bool {
        self.normal_rand() < 0.5
    }
}

/// `IntegerStringChromosomePolicy`: per position the exclusive upper bound of
/// the values (`setMax`).
#[derive(Debug, Clone)]
pub struct IntegerStringPolicy {
    maxs: Vec<i32>,
}

impl IntegerStringPolicy {
    pub fn new(size: usize) -> Self {
        IntegerStringPolicy {
            maxs: vec![10; size],
        }
    }

    pub fn set_max(&mut self, pos: usize, max: i32) {
        self.maxs[pos] = max;
    }

    /// A different value than `current` in `[0, max)`
    pub fn mutate(&self, rng: &mut RandomUtil, pos: usize, current: i32) -> i32 {
        let val = rng.random_int(0, self.maxs[pos] - 1);
        if val >= current { val + 1 } else { val }
    }

    pub fn initialize(&self, rng: &mut RandomUtil, pos: usize) -> i32 {
        rng.random_int(0, self.maxs[pos])
    }
}

/// `StringChromosomeBase::initialize`
pub fn initialize_string(policy: &IntegerStringPolicy, rng: &mut RandomUtil, string: &mut [i32]) {
    for (i, v) in string.iter_mut().enumerate() {
        *v = policy.initialize(rng, i);
    }
}

/// `StringChromosomeBase::mutate`: mutates at least one position; a negative
/// `p_mutate` means `1 / length`.
pub fn mutate_string(
    policy: &IntegerStringPolicy,
    rng: &mut RandomUtil,
    string: &mut [i32],
    p_mutate: f64,
) {
    let mut p_mutate = if p_mutate < 0.0 {
        1.0 / string.len() as f64
    } else {
        p_mutate
    };
    loop {
        let mut mutated = false;
        for (i, v) in string.iter_mut().enumerate() {
            if rng.normal_rand() < p_mutate {
                *v = policy.mutate(rng, i, *v);
                mutated = true;
            }
        }
        if mutated {
            return;
        }
        // upstream calls `mutate()` again with the default probability
        p_mutate = 1.0 / string.len() as f64;
    }
}

/// `StringChromosomeBase::twoPointCrossover` (the policy never switches the
/// children)
pub fn two_point_crossover(
    rng: &mut RandomUtil,
    parent1: &[i32],
    parent2: &[i32],
    child1: &mut [i32],
    child2: &mut [i32],
) {
    let length = parent1.len() as i32;
    let mut site1 = rng.random_int(0, length);
    let mut site2 = rng.random_int(0, length - 1);
    if site2 >= site1 {
        site2 += 1;
    } else {
        std::mem::swap(&mut site1, &mut site2);
    }
    let (site1, site2) = (site1 as usize, site2 as usize);
    for pos in 0..parent1.len() {
        let crossed = pos >= site1 && pos < site2;
        child1[pos] = if crossed { parent2[pos] } else { parent1[pos] };
        child2[pos] = if crossed { parent1[pos] } else { parent2[pos] };
    }
}

/// The function of a [`GaOperation`]: fills the children from the parents.
pub type OperationFn<P> = fn(
    &mut P,
    &mut RandomUtil,
    &[&<P as GaProblem>::Chromosome],
    &mut [<P as GaProblem>::Chromosome],
);

/// `GaOperation`: how many parents it takes, children it makes, its weight
/// and the function that performs it.
pub struct GaOperation<P: GaProblem> {
    pub n_parents: usize,
    pub n_children: usize,
    pub weight: f64,
    pub function: OperationFn<P>,
}

impl<P: GaProblem> std::fmt::Debug for GaOperation<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GaOperation")
            .field("n_parents", &self.n_parents)
            .field("n_children", &self.n_children)
            .field("weight", &self.weight)
            .finish_non_exhaustive()
    }
}

/// What `LinkedPopLinearSel` needs from the chromosome type and the
/// population policy (upstream template parameters).
pub trait GaProblem: Sized {
    type Chromosome;

    fn population_size(&self) -> usize;
    fn selection_pressure(&self) -> f64;
    fn operations(&self) -> Vec<GaOperation<Self>>;
    fn create_chromosome(&self) -> Self::Chromosome;
    /// `initialize` followed by `score`
    fn initialize(&mut self, rng: &mut RandomUtil, c: &mut Self::Chromosome);
    fn score(&mut self, c: &mut Self::Chromosome) -> Result<f64>;
    fn fitness(&self, c: &Self::Chromosome) -> f64;
    fn is_ok(&self, c: &Self::Chromosome) -> bool;
    /// equality of the genes
    fn equals(&self, a: &Self::Chromosome, b: &Self::Chromosome) -> bool;
}

const SELECT_START: f64 = 5000.0;

/// `LinkedPopLinearSel`
#[allow(missing_debug_implementations)]
pub struct LinkedPopLinearSel<P: GaProblem> {
    popsize: usize,
    operations: Vec<GaOperation<P>>,
    /// ascending by fitness, equal fitnesses in insertion order (a multimap)
    population: Vec<(f64, P::Chromosome)>,
    free_chromosomes: Vec<P::Chromosome>,
    n_operations: usize,
    n_duplicates: usize,
    n_added: usize,
    n_fail: usize,
    total_scaled_fitness: f64,
    scaled_fitness_step: f64,
    total_operator_weights: f64,
    best_score: f64,
}

impl<P: GaProblem> LinkedPopLinearSel<P> {
    pub fn new(problem: &P) -> Result<Self> {
        let popsize = problem.population_size();
        let operations = problem.operations();
        let total_operator_weights = operations.iter().map(|o| o.weight).sum();
        let selection_pressure = problem.selection_pressure();
        // the ratio between the fitness of the worst and the best individual
        if selection_pressure < 1.0 {
            return Err(crate::RDKitError::Value(format!(
                "LinkedPopLinearSel: Selection pressure should be greater than 1.0 [current value is {selection_pressure:.2}]"
            )));
        }
        let end_fitness = selection_pressure * SELECT_START;
        let total_scaled_fitness = popsize as f64 * ((SELECT_START + end_fitness) / 2.0);
        let scaled_fitness_step = (((2.0 * total_scaled_fitness) / popsize as f64)
            - (2.0 * SELECT_START))
            / (popsize as f64 - 1.0);
        Ok(LinkedPopLinearSel {
            popsize,
            operations,
            population: Vec::new(),
            free_chromosomes: Vec::with_capacity(10),
            n_operations: 0,
            n_duplicates: 0,
            n_added: 0,
            n_fail: 0,
            total_scaled_fitness,
            scaled_fitness_step,
            total_operator_weights,
            best_score: -f64::MAX,
        })
    }

    pub fn best_score(&self) -> f64 {
        self.best_score
    }

    /// `create`: the initial random population
    pub fn create(&mut self, problem: &mut P, rng: &mut RandomUtil) -> Result<()> {
        while self.population.len() < self.popsize {
            let mut chromosome = problem.create_chromosome();
            problem.initialize(rng, &mut chromosome);
            problem.score(&mut chromosome)?;
            self.add_to_population(problem, chromosome);
        }
        Ok(())
    }

    /// roulette wheel selection of the index of a parent
    fn select_parent(&self, rng: &mut RandomUtil) -> usize {
        let val = rng.normal_rand() * self.total_scaled_fitness;
        let mut sum = SELECT_START;
        let mut current_fitness = SELECT_START;
        for i in 0..self.popsize {
            if val <= sum {
                return i;
            }
            current_fitness += self.scaled_fitness_step;
            sum += current_fitness;
        }
        // numeric errors
        self.population.len() - 1
    }

    /// `iterate`: applies a single operation
    pub fn iterate(&mut self, problem: &mut P, rng: &mut RandomUtil) -> Result<()> {
        let mut total = 0.0;
        let val = rng.normal_rand() * self.total_operator_weights;
        let mut selected = self.operations.len() - 1;
        for (i, operation) in self.operations.iter().enumerate() {
            total += operation.weight;
            if val <= total {
                selected = i;
                break;
            }
        }
        let (n_children, n_parents, function) = {
            let op = &self.operations[selected];
            (op.n_children, op.n_parents, op.function)
        };
        let mut children = Vec::with_capacity(n_children);
        for _ in 0..n_children {
            children.push(
                self.free_chromosomes
                    .pop()
                    .unwrap_or_else(|| problem.create_chromosome()),
            );
        }
        let parent_indices: Vec<usize> = (0..n_parents).map(|_| self.select_parent(rng)).collect();
        {
            let parents: Vec<&P::Chromosome> = parent_indices
                .iter()
                .map(|&i| &self.population[i].1)
                .collect();
            function(problem, rng, &parents, &mut children);
        }
        for mut child in children {
            problem.score(&mut child)?;
            self.add_to_population(problem, child);
        }
        self.n_operations += 1;
        Ok(())
    }

    fn add_to_population(&mut self, problem: &P, chromosome: P::Chromosome) -> bool {
        if !problem.is_ok(&chromosome) {
            self.n_fail += 1;
            self.free_chromosomes.push(chromosome);
            return false;
        }
        let fitness = problem.fitness(&chromosome);
        // any chromosome with the same gene has the same fitness
        let duplicate = self
            .population
            .iter()
            .filter(|(f, _)| *f == fitness)
            .any(|(_, other)| problem.equals(&chromosome, other));
        if duplicate {
            self.n_duplicates += 1;
            self.free_chromosomes.push(chromosome);
            return false;
        }
        if self.population.len() >= self.popsize {
            // the worst individual is at the head
            let (_, worst) = self.population.remove(0);
            self.free_chromosomes.push(worst);
        }
        let at = self.population.partition_point(|(f, _)| *f <= fitness);
        self.population.insert(at, (fitness, chromosome));
        self.n_added += 1;
        if fitness > self.best_score {
            self.best_score = fitness;
        }
        true
    }

    /// The chromosomes whose fitness is within `tolerance` of the best one,
    /// best first (`getTiedBest`).
    pub fn tied_best(&self, tolerance: f64) -> Vec<&P::Chromosome> {
        let mut iter = self.population.iter().rev();
        let Some((best, first)) = iter.next() else {
            return Vec::new();
        };
        let mut ties = vec![first];
        for (f, c) in iter {
            if (f - best).abs() < tolerance {
                ties.push(c);
            } else {
                break;
            }
        }
        ties
    }

    pub fn best(&self) -> &P::Chromosome {
        &self.population.last().expect("empty population").1
    }

    /// `info`
    pub fn info(&self) -> String {
        format!(
            "Op {:5} Fit {:7.3} Added {:5} Dups {:5} Niche {:5} Fail {:4}",
            self.n_operations, self.best_score, self.n_added, self.n_duplicates, 0, self.n_fail
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normal_rand_is_in_unit_interval_and_reproducible() {
        let mut a = RandomUtil {
            rng: Mt19937::new(7),
        };
        let mut b = a.clone();
        for _ in 0..1000 {
            let x = a.normal_rand();
            assert!((0.0..1.0).contains(&x));
            assert_eq!(x, b.normal_rand());
        }
    }

    #[test]
    fn mutation_changes_a_position() {
        let mut rng = RandomUtil {
            rng: Mt19937::new(1),
        };
        let mut policy = IntegerStringPolicy::new(5);
        for i in 0..5 {
            policy.set_max(i, 3);
        }
        let mut s = vec![0; 5];
        initialize_string(&policy, &mut rng, &mut s);
        let before = s.clone();
        mutate_string(&policy, &mut rng, &mut s, -1.0);
        assert_ne!(before, s);
        assert!(s.iter().all(|&v| (0..3).contains(&v)));
    }

    #[test]
    fn crossover_swaps_a_segment() {
        let mut rng = RandomUtil {
            rng: Mt19937::new(3),
        };
        let p1 = [0, 0, 0, 0, 0, 0];
        let p2 = [1, 1, 1, 1, 1, 1];
        let (mut c1, mut c2) = ([9; 6], [9; 6]);
        two_point_crossover(&mut rng, &p1, &p2, &mut c1, &mut c2);
        for i in 0..6 {
            assert_eq!(c1[i] + c2[i], 1);
            assert_eq!(c1[i], 1 - c2[i]);
        }
        assert!(c1.contains(&1) && c1.contains(&0));
    }
}
