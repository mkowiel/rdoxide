//! Random number generators reproducing the ones used by RDKit, so that
//! seeded algorithms produce the same results as upstream.
//!
//! Ported from RDKit `Code/RDGeneral/utils.{h,cpp}` (BSD-3-Clause).
//!
//! * [`MinStdRand`]: `boost::minstd_rand` (the global `RDKit::getRandomGenerator`).
//! * [`MersenneTwister`]: `boost::random::mersenne_twister` (`boost::mt19937` is
//!   [`Mt19937`]) with `boost::uniform_int`.
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

/// `boost::random::mersenne_twister` with 32-bit words, mt19937's tempering
/// and `N`/`M` state words (`boost::mt19937` is `N = 624, M = 397`;
/// fingerprint bit generation uses `N = 4, M = 2`).
#[derive(Debug, Clone)]
pub struct MersenneTwister<const N: usize, const M: usize> {
    state: [u32; N],
    idx: usize,
}

/// `boost::mt19937`
pub type Mt19937 = MersenneTwister<624, 397>;

impl<const N: usize, const M: usize> MersenneTwister<N, M> {
    pub fn new(seed: u32) -> Self {
        let mut state = [0u32; N];
        state[0] = seed;
        for i in 1..N {
            state[i] = 1_812_433_253u32
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        MersenneTwister { state, idx: N }
    }

    /// Re-seeds the generator (`generator.seed(value)`).
    pub fn seed(&mut self, seed: u32) {
        *self = Self::new(seed);
    }

    fn twist(&mut self) {
        for i in 0..N {
            let y = (self.state[i] & 0x8000_0000) | (self.state[(i + 1) % N] & 0x7fff_ffff);
            let mut v = self.state[(i + M) % N] ^ (y >> 1);
            if y & 1 != 0 {
                v ^= 0x9908_b0df;
            }
            self.state[i] = v;
        }
        self.idx = 0;
    }

    /// The next raw 32-bit output.
    pub fn next_u32(&mut self) -> u32 {
        if self.idx >= N {
            self.twist();
        }
        let mut y = self.state[self.idx];
        self.idx += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    /// Bucket-rejection draw of `boost::random::detail::generate_uniform_int`
    /// for a 32-bit engine: a value in `[0, range]`.
    fn draw_up_to(&mut self, range: u32) -> u32 {
        if range == 0 {
            return 0;
        }
        if range == u32::MAX {
            return self.next_u32();
        }
        // engine range [0, 2^32 - 1] is wider than the requested one
        let mut bucket_size = u32::MAX / (range + 1);
        if u32::MAX % (range + 1) == range {
            bucket_size += 1;
        }
        loop {
            let r = self.next_u32() / bucket_size;
            if r <= range {
                return r;
            }
        }
    }

    /// `boost::uniform_int<int>(lo, hi)(engine)`.
    pub fn uniform_int(&mut self, lo: i32, hi: i32) -> i32 {
        assert!(lo <= hi, "uniform_int: empty range");
        let range = (hi as i64 - lo as i64) as u32;
        (lo as i64 + i64::from(self.draw_up_to(range))) as i32
    }

    /// `boost::random::uniform_int_distribution<uint32_t>(lo, hi)(engine)`.
    pub fn uniform_u32(&mut self, lo: u32, hi: u32) -> u32 {
        assert!(lo <= hi, "uniform_u32: empty range");
        lo + self.draw_up_to(hi - lo)
    }
}

#[cfg(test)]
mod mersenne_tests {
    use super::*;

    #[test]
    fn mt19937_reference_values() {
        // the first outputs of std::mt19937 seeded with 5489
        let mut mt = Mt19937::new(5489);
        assert_eq!(mt.next_u32(), 3_499_211_612);
        assert_eq!(mt.next_u32(), 581_869_302);
        // and the 10000th output of the default-seeded engine required by the C++ standard
        let mut mt = Mt19937::new(5489);
        let tenthousandth = (0..10_000).map(|_| mt.next_u32()).last().unwrap();
        assert_eq!(tenthousandth, 4_123_659_995);
    }

    #[test]
    fn uniform_int_stays_in_range_and_handles_degenerate_ranges() {
        let mut mt = Mt19937::new(42);
        assert_eq!(mt.uniform_int(7, 7), 7);
        for _ in 0..1000 {
            assert!((0..10).contains(&mt.uniform_int(0, 9)));
        }
        assert_eq!(
            Mt19937::new(1).uniform_int(i32::MIN, i32::MAX),
            Mt19937::new(1).next_u32() as i32 ^ i32::MIN
        );
    }
}
