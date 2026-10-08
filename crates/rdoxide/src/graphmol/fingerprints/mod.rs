//! Fingerprint generators (upstream `GraphMol/Fingerprints/FingerprintGenerator.{h,cpp}`,
//! `FingerprintUtil.{h,cpp}`).

pub mod atom_pair;
pub mod legacy;
pub mod morgan;
pub mod rdkit_fp;
pub mod torsion;

use std::collections::BTreeMap;

use crate::datastructs::sparse_int_vect::{IndexKind, SparseIntVect};
use crate::datastructs::{BitVect, ExplicitBitVect, SparseBitVect};
use crate::error::{RDKitError, Result};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;

/// `gboost::hash_combine` with RDKit's 32-bit `hash_result_t`.
#[inline]
pub fn hash_combine(seed: &mut u32, v: u32) {
    *seed ^= v
        .wrapping_add(0x9e37_79b9)
        .wrapping_add(*seed << 6)
        .wrapping_add(*seed >> 2);
}

/// `gboost::hash_value(std::pair<A, B>)` for 32-bit integral members.
#[inline]
pub fn hash_pair(a: u32, b: u32) -> u32 {
    let mut seed = 0;
    hash_combine(&mut seed, a);
    hash_combine(&mut seed, b);
    seed
}

/// `gboost::hash<std::vector<uint32_t>>`
pub fn hash_vect(v: &[u32]) -> u32 {
    let mut seed = 0;
    for &x in v {
        hash_combine(&mut seed, x);
    }
    seed
}

/// The `boost::random::mersenne_twister<uint32_t, 32, 4, 2, 31, 0x9908b0df,
/// 11, 7, 0x9d2c5680, 15, 0xefc60000, 18, 3346425566U>` used by
/// `getFingerprintHelper`: mt19937's tempering with a 4-word state (n=4, m=2).
#[derive(Debug)]
pub struct Mt19937 {
    state: [u32; Self::N],
    idx: usize,
}

impl Mt19937 {
    const N: usize = 4;
    const M: usize = 2;

    pub fn new(seed: u32) -> Self {
        let mut state = [0u32; Self::N];
        state[0] = seed;
        for i in 1..Self::N {
            state[i] = 1_812_433_253u32
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Mt19937 {
            state,
            idx: Self::N,
        }
    }
    fn twist(&mut self) {
        for i in 0..Self::N {
            let y = (self.state[i] & 0x8000_0000) | (self.state[(i + 1) % Self::N] & 0x7fff_ffff);
            let mut v = self.state[(i + Self::M) % Self::N] ^ (y >> 1);
            if y & 1 != 0 {
                v ^= 0x9908_b0df;
            }
            self.state[i] = v;
        }
        self.idx = 0;
    }
    pub fn next_u32(&mut self) -> u32 {
        if self.idx >= Self::N {
            self.twist();
        }
        let mut y = self.state[self.idx];
        self.idx += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }
    /// `boost::uniform_int<>(0, INT_MAX)` driven by this engine: the engine
    /// range is twice the target range, so values fall into buckets of two.
    pub fn uniform_0_int_max(&mut self) -> u32 {
        self.next_u32() / 2
    }
}

/// `FingerprintArguments` (common part).
#[derive(Clone, Debug)]
pub struct FingerprintArguments {
    pub count_simulation: bool,
    pub include_chirality: bool,
    pub count_bounds: Vec<u32>,
    pub fp_size: u32,
    pub num_bits_per_feature: u32,
}

impl FingerprintArguments {
    pub fn common_arguments_string(&self) -> String {
        format!(
            "Common arguments : countSimulation={} fpSize={} bitsPerFeature={} includeChirality={}",
            u8::from(self.count_simulation),
            self.fp_size,
            self.num_bits_per_feature,
            u8::from(self.include_chirality)
        )
    }
}

/// `AdditionalOutput`; each member is collected only when allocated.
#[derive(Clone, Debug, Default)]
pub struct AdditionalOutput {
    pub atom_to_bits: Option<Vec<Vec<u64>>>,
    pub bit_info_map: Option<BTreeMap<u64, Vec<(u32, u32)>>>,
    pub bit_paths: Option<BTreeMap<u64, Vec<Vec<u32>>>>,
    pub atom_counts: Option<Vec<u32>>,
    pub atoms_per_bit: Option<BTreeMap<u64, Vec<Vec<i32>>>>,
}

impl AdditionalOutput {
    /// `reinitAdditionalOutput` (note: `atomsPerBit` is not cleared upstream)
    fn reinit(&mut self, n_atoms: usize) {
        if let Some(v) = &mut self.atom_counts {
            v.clear();
            v.resize(n_atoms, 0);
        }
        if let Some(v) = &mut self.atom_to_bits {
            v.clear();
            v.resize(n_atoms, Vec::new());
        }
        if let Some(m) = &mut self.bit_info_map {
            m.clear();
        }
        if let Some(m) = &mut self.bit_paths {
            m.clear();
        }
    }
    /// `setupTempAdditionalOutput`: an empty output with the same members allocated
    fn allocated_like(&self) -> Self {
        AdditionalOutput {
            atom_to_bits: self.atom_to_bits.as_ref().map(|_| Vec::new()),
            bit_info_map: self.bit_info_map.as_ref().map(|_| BTreeMap::new()),
            bit_paths: self.bit_paths.as_ref().map(|_| BTreeMap::new()),
            atom_counts: self.atom_counts.as_ref().map(|_| Vec::new()),
            atoms_per_bit: self.atoms_per_bit.as_ref().map(|_| BTreeMap::new()),
        }
    }
    /// `duplicateAdditionalOutputBit(old=self, new=new_ao)`
    fn duplicate_bit(&self, new_ao: &mut AdditionalOutput, orig: u64, new_bit: u64) {
        if let (Some(old), Some(nw)) = (&self.atom_to_bits, &mut new_ao.atom_to_bits) {
            if nw.is_empty() {
                nw.resize(old.len(), Vec::new());
            }
            for (i, nv) in old.iter().enumerate() {
                if nv.contains(&orig) {
                    nw[i].push(new_bit);
                }
            }
        }
        if let (Some(old), Some(nw)) = (&self.bit_info_map, &mut new_ao.bit_info_map)
            && let Some(v) = old.get(&orig)
        {
            nw.insert(new_bit, v.clone());
        }
        if let (Some(old), Some(nw)) = (&self.bit_paths, &mut new_ao.bit_paths)
            && let Some(v) = old.get(&orig)
        {
            nw.insert(new_bit, v.clone());
        }
        if let (Some(old), Some(nw)) = (&self.atoms_per_bit, &mut new_ao.atoms_per_bit)
            && let Some(v) = old.get(&orig)
        {
            nw.insert(new_bit, v.clone());
        }
    }
}

/// `FingerprintFuncArguments`
#[derive(Debug, Default)]
pub struct FingerprintFuncArguments<'a> {
    pub from_atoms: Option<&'a [u32]>,
    pub ignore_atoms: Option<&'a [u32]>,
    pub conf_id: i32,
    pub additional_output: Option<&'a mut AdditionalOutput>,
    pub custom_atom_invariants: Option<&'a [u32]>,
    pub custom_bond_invariants: Option<&'a [u32]>,
}

/// One generated environment (`AtomEnvironment`).
#[derive(Debug)]
pub enum AtomEnvironment {
    Morgan(morgan::MorganAtomEnv),
    AtomPair(atom_pair::AtomPairAtomEnv),
    Torsion(torsion::TopologicalTorsionAtomEnv),
    RDKitFP(rdkit_fp::RDKitFPAtomEnv),
}

impl AtomEnvironment {
    /// `AtomEnvironment::getBitId`
    fn bit_id(
        &self,
        args: &FingerprintArguments,
        atom_invs: Option<&[u32]>,
        hash_results: bool,
    ) -> Result<u64> {
        match self {
            AtomEnvironment::Morgan(e) => Ok(e.code),
            AtomEnvironment::AtomPair(e) => e.bit_id(args, atom_invs, hash_results),
            AtomEnvironment::Torsion(e) => Ok(e.bit_id),
            AtomEnvironment::RDKitFP(e) => Ok(e.bit_id),
        }
    }
    fn update_additional_output(
        &self,
        mol: &ROMol,
        ao: &mut AdditionalOutput,
        bit_id: u64,
    ) -> Result<()> {
        match self {
            AtomEnvironment::Morgan(e) => {
                e.update_additional_output(mol, ao, bit_id);
                Ok(())
            }
            AtomEnvironment::AtomPair(e) => e.update_additional_output(ao, bit_id),
            AtomEnvironment::Torsion(e) => e.update_additional_output(ao, bit_id),
            AtomEnvironment::RDKitFP(e) => e.update_additional_output(ao, bit_id),
        }
    }
}

/// `AtomInvariantsGenerator`
#[derive(Clone, Debug)]
pub enum AtomInvariantsGenerator {
    Morgan {
        include_ring_membership: bool,
    },
    AtomPair {
        include_chirality: bool,
        topological_torsion_correction: bool,
    },
    RDKitFP,
    /// `MorganFeatureAtomInvGenerator` (`None`: the default feature SMARTS)
    MorganFeature {
        patterns: Option<Vec<ROMol>>,
    },
}

impl AtomInvariantsGenerator {
    pub fn atom_invariants(&self, mol: &ROMol) -> Result<Vec<u32>> {
        match self {
            AtomInvariantsGenerator::Morgan {
                include_ring_membership,
            } => morgan::get_connectivity_invariants(mol, *include_ring_membership),
            AtomInvariantsGenerator::AtomPair {
                include_chirality,
                topological_torsion_correction,
            } => atom_pair::atom_pair_invariants(
                mol,
                *include_chirality,
                *topological_torsion_correction,
            ),
            AtomInvariantsGenerator::RDKitFP => Ok(rdkit_fp::rdkit_fp_atom_invariants(mol)),
            AtomInvariantsGenerator::MorganFeature { patterns } => {
                morgan::get_feature_invariants(mol, patterns.as_deref())
            }
        }
    }
    pub fn info_string(&self) -> String {
        match self {
            AtomInvariantsGenerator::Morgan {
                include_ring_membership,
            } => format!(
                "MorganInvariantGenerator includeRingMembership={}",
                u8::from(*include_ring_membership)
            ),
            AtomInvariantsGenerator::AtomPair {
                topological_torsion_correction,
                ..
            } => format!(
                "AtomPairInvariantGenerator topologicalTorsionCorrection={}",
                u8::from(*topological_torsion_correction)
            ),
            AtomInvariantsGenerator::RDKitFP => "RDKitFPAtomInvGenerator".to_string(),
            AtomInvariantsGenerator::MorganFeature { .. } => {
                "MorganFeatureInvariantGenerator".to_string()
            }
        }
    }
}

/// `BondInvariantsGenerator`
#[derive(Clone, Debug)]
pub enum BondInvariantsGenerator {
    Morgan {
        use_bond_types: bool,
        use_chirality: bool,
    },
}

impl BondInvariantsGenerator {
    pub fn bond_invariants(&self, mol: &ROMol) -> Result<Vec<u32>> {
        match self {
            BondInvariantsGenerator::Morgan {
                use_bond_types,
                use_chirality,
            } => Ok(morgan::morgan_bond_invariants(
                mol,
                *use_bond_types,
                *use_chirality,
            )),
        }
    }
    pub fn info_string(&self) -> String {
        match self {
            BondInvariantsGenerator::Morgan {
                use_bond_types,
                use_chirality,
            } => format!(
                "MorganInvariantGenerator useBondTypes={} useChirality={}",
                u8::from(*use_bond_types),
                u8::from(*use_chirality)
            ),
        }
    }
}

/// The environment generator together with its type-specific arguments.
#[derive(Clone, Debug)]
pub enum EnvGenerator {
    Morgan(morgan::MorganArguments),
    AtomPair(atom_pair::AtomPairArguments),
    TopologicalTorsion(torsion::TopologicalTorsionArguments),
    RDKitFP(rdkit_fp::RDKitFPArguments),
}

/// `FingerprintGenerator<OutputType>`; `output_64` selects the OutputType.
#[derive(Clone, Debug)]
pub struct FingerprintGenerator {
    pub env: EnvGenerator,
    pub args: FingerprintArguments,
    pub atom_inv_gen: Option<AtomInvariantsGenerator>,
    pub bond_inv_gen: Option<BondInvariantsGenerator>,
    pub output_64: bool,
}

/// Result of `getFingerprintHelper`: bit id -> count
type Counts = SparseIntVect;

impl FingerprintGenerator {
    fn output_kind(&self) -> IndexKind {
        if self.output_64 {
            IndexKind::U64
        } else {
            IndexKind::U32
        }
    }
    /// `AtomEnvironmentGenerator::getResultSize`
    fn result_size(&self) -> u64 {
        if let EnvGenerator::TopologicalTorsion(t) = &self.env {
            return t.result_size(self.args.include_chirality);
        }
        if let EnvGenerator::AtomPair(_) = self.env {
            let chiral = if self.args.include_chirality {
                atom_pair::NUM_CHIRAL_BITS
            } else {
                0
            };
            return 1u64 << (atom_pair::NUM_ATOM_PAIR_FINGERPRINT_BITS + 2 * chiral);
        }
        if self.output_64 {
            u64::MAX
        } else {
            u64::from(u32::MAX)
        }
    }

    pub fn info_string(&self) -> String {
        let sep = " --- ";
        let (args, env) = match &self.env {
            EnvGenerator::Morgan(m) => (m.info_string(), "MorganEnvironmentGenerator"),
            EnvGenerator::AtomPair(a) => (a.info_string(), "AtomPairEnvironmentGenerator"),
            EnvGenerator::TopologicalTorsion(t) => {
                (t.info_string(), "TopologicalTorsionEnvGenerator")
            }
            EnvGenerator::RDKitFP(r) => (r.info_string(), "RDKitFPEnvGenerator"),
        };
        format!(
            "{}{sep}{}{sep}{}{sep}{}{sep}{}",
            self.args.common_arguments_string(),
            args,
            env,
            self.atom_inv_gen
                .as_ref()
                .map_or("No atom invariants generator".to_string(), |g| g
                    .info_string()),
            self.bond_inv_gen
                .as_ref()
                .map_or("No bond invariants generator".to_string(), |g| g
                    .info_string()),
        )
    }

    fn environments(
        &self,
        mol: &ROMol,
        fargs: &FingerprintFuncArguments<'_>,
        atom_invs: Option<&[u32]>,
        bond_invs: Option<&[u32]>,
        hash_results: bool,
    ) -> Result<Vec<AtomEnvironment>> {
        let from_atoms = fargs.from_atoms;
        match &self.env {
            EnvGenerator::RDKitFP(r) => Ok(rdkit_fp::get_environments(
                mol,
                r,
                from_atoms,
                fargs.ignore_atoms,
                atom_invs.ok_or_else(|| RDKitError::precondition("bad atomInvariants size"))?,
            )?
            .into_iter()
            .map(AtomEnvironment::RDKitFP)
            .collect()),
            EnvGenerator::TopologicalTorsion(t) => Ok(torsion::get_environments(
                mol,
                t,
                self.args.include_chirality,
                from_atoms,
                fargs.ignore_atoms,
                atom_invs.ok_or_else(|| RDKitError::precondition("bad atom invariants size"))?,
                hash_results,
            )?
            .into_iter()
            .map(AtomEnvironment::Torsion)
            .collect()),
            EnvGenerator::AtomPair(a) => Ok(atom_pair::get_environments(
                mol,
                a,
                from_atoms,
                fargs.ignore_atoms,
                fargs.conf_id,
                fargs.additional_output.as_deref(),
            )?
            .into_iter()
            .map(AtomEnvironment::AtomPair)
            .collect()),
            EnvGenerator::Morgan(m) => Ok(morgan::get_environments(
                mol,
                m,
                &self.args,
                from_atoms,
                atom_invs.ok_or_else(|| RDKitError::precondition("bad atom invariants size"))?,
                bond_invs.ok_or_else(|| RDKitError::precondition("bad bond invariants size"))?,
            )?
            .into_iter()
            .map(AtomEnvironment::Morgan)
            .collect()),
        }
    }

    /// `getFingerprintHelper`
    fn helper(
        &self,
        mol: &ROMol,
        fargs: &mut FingerprintFuncArguments<'_>,
        fp_size: u64,
    ) -> Result<Counts> {
        let stereo_mol;
        let lmol = if self.args.include_chirality && !mol.has_prop(cp::_StereochemDone) {
            let mut m = mol.clone();
            crate::graphmol::chirality::legacy::assign_stereochemistry(
                &mut m, false, false, false,
            )?;
            stereo_mol = m;
            &stereo_mol
        } else {
            mol
        };
        if let Some(ao) = fargs.additional_output.as_deref_mut() {
            ao.reinit(mol.num_atoms() as usize);
        }
        let atom_invs = match fargs.custom_atom_invariants {
            Some(v) => Some(v.to_vec()),
            None => self
                .atom_inv_gen
                .as_ref()
                .map(|g| g.atom_invariants(mol))
                .transpose()?,
        };
        let bond_invs = match fargs.custom_bond_invariants {
            Some(v) => Some(v.to_vec()),
            None => self
                .bond_inv_gen
                .as_ref()
                .map(|g| g.bond_invariants(mol))
                .transpose()?,
        };
        let envs = self.environments(
            lmol,
            fargs,
            atom_invs.as_deref(),
            bond_invs.as_deref(),
            fp_size != 0,
        )?;

        let kind = self.output_kind();
        let len = if fp_size != 0 {
            fp_size
        } else {
            self.result_size()
        };
        let mut res = SparseIntVect::new(kind, i128::from(len));
        let mask = if self.output_64 {
            u64::MAX
        } else {
            u64::from(u32::MAX)
        };
        let bump = |res: &mut SparseIntVect, bit: u64| -> Result<()> {
            let i = i128::from(bit);
            let v = res.get_val(i)?;
            res.set_val(i, v.wrapping_add(1))
        };
        let mut rng = None;
        for env in &envs {
            let seed = env.bit_id(&self.args, atom_invs.as_deref(), fp_size != 0)? & mask;
            let mut bit_id = seed;
            if fp_size != 0 {
                bit_id %= fp_size;
            }
            bump(&mut res, bit_id)?;
            if let Some(ao) = fargs.additional_output.as_deref_mut() {
                env.update_additional_output(lmol, ao, bit_id)?;
            }
            if self.args.num_bits_per_feature > 1 {
                let rg = rng.insert(Mt19937::new(seed as u32));
                for _ in 1..self.args.num_bits_per_feature {
                    let mut b = u64::from(rg.uniform_0_int_max()) & mask;
                    if fp_size != 0 {
                        b %= fp_size;
                    }
                    bump(&mut res, b)?;
                    if let Some(ao) = fargs.additional_output.as_deref_mut() {
                        env.update_additional_output(lmol, ao, b)?;
                    }
                }
            }
        }
        Ok(res)
    }

    /// Runs `helper` with a temporary AdditionalOutput when count simulation
    /// is active, then sets the (simulated) bits via `set_bit`.
    fn simulated<B>(
        &self,
        mol: &ROMol,
        fargs: &mut FingerprintFuncArguments<'_>,
        effective_size: u64,
        mut result: B,
        set_bit: impl Fn(&mut B, u64) -> Result<()>,
    ) -> Result<B> {
        let sim = self.args.count_simulation;
        let mut temp_ao = None;
        let orig = fargs.additional_output.take();
        let orig = match orig {
            Some(ao) if sim => {
                temp_ao = Some(ao.allocated_like());
                ao.reinit(mol.num_atoms() as usize);
                Some(ao)
            }
            other => {
                fargs.additional_output = other;
                None
            }
        };
        let counts = {
            let mut targs = FingerprintFuncArguments {
                from_atoms: fargs.from_atoms,
                ignore_atoms: fargs.ignore_atoms,
                conf_id: fargs.conf_id,
                additional_output: if sim {
                    temp_ao.as_mut()
                } else {
                    fargs.additional_output.as_deref_mut()
                },
                custom_atom_invariants: fargs.custom_atom_invariants,
                custom_bond_invariants: fargs.custom_bond_invariants,
            };
            self.helper(mol, &mut targs, effective_size)?
        };
        let mut orig = orig;
        let bounds = &self.args.count_bounds;
        for (&bit, &count) in counts.nonzero_elements() {
            let bit = bit as u64;
            if sim {
                for (i, &b) in bounds.iter().enumerate() {
                    if count >= b as i32 {
                        let nbit = bit.wrapping_mul(bounds.len() as u64).wrapping_add(i as u64);
                        set_bit(&mut result, nbit)?;
                        if let (Some(t), Some(o)) = (&temp_ao, orig.as_deref_mut()) {
                            t.duplicate_bit(o, bit, nbit);
                        }
                    }
                }
            } else {
                set_bit(&mut result, bit)?;
            }
        }
        if let (Some(o), Some(t)) = (orig.as_deref_mut(), temp_ao)
            && o.atom_counts.is_some()
        {
            o.atom_counts = t.atom_counts;
        }
        if let Some(o) = orig {
            fargs.additional_output = Some(o);
        }
        Ok(result)
    }

    /// `getFingerprint`
    pub fn fingerprint(
        &self,
        mol: &ROMol,
        fargs: &mut FingerprintFuncArguments<'_>,
    ) -> Result<ExplicitBitVect> {
        let mut effective = self.args.fp_size;
        if self.args.count_simulation {
            if self.args.count_bounds.is_empty() {
                return Err(RDKitError::Value("Count bounds are empty".into()));
            }
            if self.args.count_bounds.len() >= effective as usize {
                return Err(RDKitError::Value(
                    "Count bounds size is >= fingerprint size".into(),
                ));
            }
            effective /= self.args.count_bounds.len() as u32;
        }
        let result = ExplicitBitVect::new(self.args.fp_size, false);
        self.simulated(mol, fargs, u64::from(effective), result, |r, b| {
            r.set_bit(b as u32).map(|_| ())
        })
    }

    /// `getSparseFingerprint`
    pub fn sparse_fingerprint(
        &self,
        mol: &ROMol,
        fargs: &mut FingerprintFuncArguments<'_>,
    ) -> Result<SparseBitVect> {
        let result_size = self.result_size().min(u64::from(u32::MAX)) as u32;
        let mut effective = result_size;
        if self.args.count_simulation {
            effective /= self.args.count_bounds.len() as u32;
        }
        let result = SparseBitVect::new(result_size);
        self.simulated(mol, fargs, u64::from(effective), result, |r, b| {
            r.set_bit(b as u32).map(|_| ())
        })
    }

    /// `getCountFingerprint`
    pub fn count_fingerprint(
        &self,
        mol: &ROMol,
        fargs: &mut FingerprintFuncArguments<'_>,
    ) -> Result<SparseIntVect> {
        let counts = self.helper(mol, fargs, u64::from(self.args.fp_size))?;
        let mut res = SparseIntVect::new(IndexKind::U32, i128::from(self.args.fp_size));
        for (&k, &v) in counts.nonzero_elements() {
            res.set_val(i128::from(k as u32), v)?;
        }
        Ok(res)
    }

    /// `getSparseCountFingerprint`
    pub fn sparse_count_fingerprint(
        &self,
        mol: &ROMol,
        fargs: &mut FingerprintFuncArguments<'_>,
    ) -> Result<SparseIntVect> {
        self.helper(mol, fargs, 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mt19937_reference() {
        // values observed from RDKit 2026.09.1 (RDKit FP, numBitsPerFeature=2)
        assert_eq!(Mt19937::new(4_274_652_475).uniform_0_int_max(), 54_413_874);
        assert_eq!(
            Mt19937::new(1_524_090_560).uniform_0_int_max(),
            1_135_572_127
        );
    }
}
