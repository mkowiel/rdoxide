//! Random number generators reproducing the ones used by RDKit, so that
//! seeded algorithms produce the same results as upstream.
//!
//! Ported from RDKit `Code/RDGeneral/utils.{h,cpp}` (BSD-3-Clause).
//!
//! * [`MinStdRand`]: `boost::minstd_rand` (the global `RDKit::getRandomGenerator`).
//! * [`CRand`]: glibc's `rand()`/`srand()` (TYPE_3 additive feedback), the
//!   generator behind `rdBase.SeedRandomNumberGenerator` on Linux builds.

use std::sync::Mutex;

/// `boost::minstd_rand`: Lehmer LCG, a = 48271, m = 2^31 - 1.
#[derive(Debug, Clone)]
pub struct MinStdRand {
    state: u32,
}

impl MinStdRand {
    pub const MIN: u32 = 1;
    pub const MAX: u32 = 2_147_483_646;
    const M: u64 = 2_147_483_647;

    pub fn new(seed: u32) -> Self {
        let mut r = MinStdRand { state: 1 };
        r.seed(seed);
        r
    }

    pub fn seed(&mut self, seed: u32) {
        let s = u64::from(seed) % Self::M;
        self.state = if s == 0 { 1 } else { s as u32 };
    }

    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> u32 {
        self.state = ((u64::from(self.state) * 48271) % Self::M) as u32;
        self.state
    }

    /// `boost::uniform_real<double>(0, 1)` driven by this engine.
    pub fn next_f64(&mut self) -> f64 {
        let range = f64::from(Self::MAX - Self::MIN) + 1.0;
        f64::from(self.next() - Self::MIN) / range
    }
}

/// glibc `random()`/`rand()` (TYPE_3, degree 31, separation 3).
#[derive(Debug, Clone)]
pub struct CRand {
    tab: Vec<i32>,
}

impl CRand {
    pub fn new(seed: u32) -> Self {
        let mut c = CRand { tab: Vec::new() };
        c.seed(seed);
        c
    }

    pub fn seed(&mut self, seed: u32) {
        let seed = if seed == 0 { 1 } else { seed };
        let mut r = [0i32; 34];
        r[0] = seed as i32;
        for i in 1..31 {
            let hi = i64::from(r[i - 1]) / 127_773;
            let lo = i64::from(r[i - 1]) % 127_773;
            let mut word = 16_807 * lo - 2_836 * hi;
            if word < 0 {
                word += 2_147_483_647;
            }
            r[i] = word as i32;
        }
        let mut tab: Vec<i32> = r[..31].to_vec();
        for i in 31..34 {
            tab.push(tab[i - 31]);
        }
        // discard the first 310 outputs
        let mut c = CRand { tab };
        for _ in 0..310 {
            c.raw();
        }
        *self = c;
    }

    fn raw(&mut self) -> u32 {
        let n = self.tab.len();
        let v = self.tab[n - 31].wrapping_add(self.tab[n - 3]);
        self.tab.push(v);
        if self.tab.len() > 4096 {
            self.tab.drain(..self.tab.len() - 34);
        }
        (v as u32) >> 1
    }

    /// `rand()`: value in `[0, RAND_MAX]` with `RAND_MAX = 2^31 - 1`.
    pub fn rand(&mut self) -> i32 {
        self.raw() as i32
    }
}

static GLOBAL: Mutex<Option<MinStdRand>> = Mutex::new(None);
static C_RAND: Mutex<Option<CRand>> = Mutex::new(None);

/// `RDKit::getRandomVal(seed)`: uniform double in [0, 1) from the global
/// generator (seeded with 42 at start-up); reseeds when `seed > 0`.
pub fn get_random_val(seed: i32) -> f64 {
    let mut g = GLOBAL.lock().unwrap_or_else(|e| e.into_inner());
    let rng = g.get_or_insert_with(|| MinStdRand::new(42));
    if seed > 0 {
        rng.seed(seed as u32);
    }
    rng.next_f64()
}

/// Run `f` with the global `boost::minstd_rand` engine, optionally reseeding.
pub fn with_global_generator<R>(seed: i32, f: impl FnOnce(&mut MinStdRand) -> R) -> R {
    let mut g = GLOBAL.lock().unwrap_or_else(|e| e.into_inner());
    let rng = g.get_or_insert_with(|| MinStdRand::new(42));
    if seed > 0 {
        rng.seed(seed as u32);
    }
    f(rng)
}

/// `std::srand(seed)` (Python: `rdBase.SeedRandomNumberGenerator`).
pub fn seed(seed: u32) {
    *C_RAND.lock().unwrap_or_else(|e| e.into_inner()) = Some(CRand::new(seed));
}

/// `std::rand()`
pub fn c_rand() -> i32 {
    let mut g = C_RAND.lock().unwrap_or_else(|e| e.into_inner());
    g.get_or_insert_with(|| CRand::new(1)).rand()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minstd_matches_reference() {
        // C++11 std::minstd_rand (same as boost) with default seed 1:
        // the 10000th value is 399268537
        let mut r = MinStdRand::new(1);
        let mut v = 0;
        for _ in 0..10_000 {
            v = r.next();
        }
        assert_eq!(v, 399_268_537);
    }

    #[test]
    fn glibc_rand_matches_reference() {
        // glibc: srand(1); rand() -> 1804289383, 846930886, 1681692777
        let mut r = CRand::new(1);
        assert_eq!(r.rand(), 1_804_289_383);
        assert_eq!(r.rand(), 846_930_886);
        assert_eq!(r.rand(), 1_681_692_777);
    }
}
