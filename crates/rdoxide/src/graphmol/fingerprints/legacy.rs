//! Legacy fingerprint entry points built on the generators (upstream
//! `GraphMol/Fingerprints/Fingerprints.cpp`: `RDKFingerprintMol`,
//! `getUnfoldedRDKFingerprintMol`).

use std::collections::BTreeMap;

use super::rdkit_fp::RDKitFPArguments;
use super::{
    AdditionalOutput, AtomInvariantsGenerator, EnvGenerator, FingerprintArguments,
    FingerprintFuncArguments, FingerprintGenerator,
};
use crate::datastructs::bit_ops::fold_fingerprint;
use crate::datastructs::sparse_int_vect::SparseIntVect;
use crate::datastructs::{BitVect, ExplicitBitVect};
use crate::error::{RDKitError, Result};
use crate::graphmol::romol::ROMol;

/// `RDKitFP::getRDKitFPGenerator(minPath, maxPath, useHs, branchedPaths, useBondOrder)`
fn rdkit_fp_generator(
    min_path: u32,
    max_path: u32,
    use_hs: bool,
    branched_paths: bool,
    use_bond_order: bool,
    output_64: bool,
) -> FingerprintGenerator {
    FingerprintGenerator {
        env: EnvGenerator::RDKitFP(RDKitFPArguments {
            min_path,
            max_path,
            use_hs,
            branched_paths,
            use_bond_order,
        }),
        args: FingerprintArguments {
            count_simulation: false,
            include_chirality: false,
            count_bounds: vec![1, 2, 4, 8],
            fp_size: 2048,
            num_bits_per_feature: 2,
        },
        atom_inv_gen: Some(AtomInvariantsGenerator::RDKitFP),
        bond_inv_gen: None,
        output_64,
    }
}

/// Requested extra outputs and their results.
#[derive(Debug, Default)]
pub struct PathFpOutput {
    /// bits set by each atom (`atomBits`)
    pub atom_bits: Option<Vec<Vec<u64>>>,
    /// bit -> bond paths (`bitInfo`)
    pub bit_info: Option<BTreeMap<u64, Vec<Vec<u32>>>>,
}

#[allow(clippy::too_many_arguments)]
fn check_args(
    mol: &ROMol,
    min_path: u32,
    max_path: u32,
    atom_invariants: Option<&[u32]>,
) -> Result<()> {
    if min_path == 0 {
        return Err(RDKitError::precondition("minPath==0"));
    }
    if max_path < min_path {
        return Err(RDKitError::precondition("maxPath<minPath"));
    }
    if atom_invariants.is_some_and(|v| v.len() < mol.num_atoms() as usize) {
        return Err(RDKitError::precondition("bad atomInvariants size"));
    }
    Ok(())
}

fn run<R>(
    gen_: &FingerprintGenerator,
    mol: &ROMol,
    atom_invariants: Option<&[u32]>,
    from_atoms: Option<&[u32]>,
    out: &mut PathFpOutput,
    f: impl FnOnce(&FingerprintGenerator, &ROMol, &mut FingerprintFuncArguments<'_>) -> Result<R>,
) -> Result<R> {
    let want = out.atom_bits.is_some() || out.bit_info.is_some();
    let mut ao = AdditionalOutput {
        atom_to_bits: out.atom_bits.as_ref().map(|_| Vec::new()),
        bit_paths: out.bit_info.as_ref().map(|_| BTreeMap::new()),
        ..Default::default()
    };
    let res = {
        let mut fa = FingerprintFuncArguments {
            from_atoms,
            custom_atom_invariants: atom_invariants,
            additional_output: want.then_some(&mut ao),
            ..Default::default()
        };
        f(gen_, mol, &mut fa)?
    };
    if out.atom_bits.is_some() {
        out.atom_bits = ao.atom_to_bits;
    }
    if out.bit_info.is_some() {
        out.bit_info = ao.bit_paths;
    }
    Ok(res)
}

/// `RDKFingerprintMol`
#[allow(clippy::too_many_arguments)]
pub fn rdk_fingerprint_mol(
    mol: &ROMol,
    min_path: u32,
    max_path: u32,
    fp_size: u32,
    n_bits_per_hash: u32,
    use_hs: bool,
    tgt_density: f64,
    min_size: u32,
    branched_paths: bool,
    use_bond_order: bool,
    atom_invariants: Option<&[u32]>,
    from_atoms: Option<&[u32]>,
    out: &mut PathFpOutput,
) -> Result<ExplicitBitVect> {
    check_args(mol, min_path, max_path, atom_invariants)?;
    if fp_size == 0 {
        return Err(RDKitError::precondition("fpSize==0"));
    }
    if n_bits_per_hash == 0 {
        return Err(RDKitError::precondition("nBitsPerHash==0"));
    }
    let mut gen_ = rdkit_fp_generator(
        min_path,
        max_path,
        use_hs,
        branched_paths,
        use_bond_order,
        false,
    );
    gen_.args.fp_size = fp_size;
    gen_.args.num_bits_per_feature = n_bits_per_hash;
    let mut res = run(&gen_, mol, atom_invariants, from_atoms, out, |g, m, fa| {
        g.fingerprint(m, fa)
    })?;
    if let Some(bi) = &mut out.bit_info {
        // the legacy map is keyed by uint32
        *bi = std::mem::take(bi)
            .into_iter()
            .map(|(k, v)| (k & 0xffff_ffff, v))
            .collect();
    }
    if tgt_density > 0.0 {
        while f64::from(res.num_on_bits()) / f64::from(res.num_bits()) < tgt_density
            && res.num_bits() >= 2 * min_size
        {
            res = fold_fingerprint(&res, 2)?;
        }
    }
    Ok(res)
}

/// `getUnfoldedRDKFingerprintMol`
#[allow(clippy::too_many_arguments)]
pub fn unfolded_rdk_fingerprint_mol(
    mol: &ROMol,
    min_path: u32,
    max_path: u32,
    use_hs: bool,
    branched_paths: bool,
    use_bond_order: bool,
    atom_invariants: Option<&[u32]>,
    from_atoms: Option<&[u32]>,
    out: &mut PathFpOutput,
) -> Result<SparseIntVect> {
    check_args(mol, min_path, max_path, atom_invariants)?;
    let mut gen_ = rdkit_fp_generator(
        min_path,
        max_path,
        use_hs,
        branched_paths,
        use_bond_order,
        true,
    );
    gen_.args.num_bits_per_feature = 1;
    run(&gen_, mol, atom_invariants, from_atoms, out, |g, m, fa| {
        g.sparse_count_fingerprint(m, fa)
    })
}

// ---------------------------------------------------------------------------
// MorganFingerprints.cpp

use super::BondInvariantsGenerator;
use super::atom_pair::{
    AtomPairArguments, CODE_SIZE, MAX_PATH_LEN, NUM_CHIRAL_BITS, get_atom_code,
};
use super::morgan::MorganArguments;
use super::torsion::{TopologicalTorsionArguments, get_topological_torsion_code};
use crate::datastructs::sparse_int_vect::IndexKind;
use crate::graphmol::subgraphs::find_all_paths_of_length_n;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;

/// `BitInfoMap`: bit -> (atom, radius) pairs
pub type BitInfoMap = BTreeMap<u64, Vec<(u32, u32)>>;

/// `MorganFingerprint::getMorganGenerator<uint32_t>(radius, false, useChirality,
/// useBondTypes, onlyNonzeroInvariants, includeRedundantEnvironments, nullptr,
/// nullptr, fpSize)`
fn morgan_generator32(
    radius: u32,
    use_chirality: bool,
    use_bond_types: bool,
    include_redundant: bool,
    fp_size: u32,
) -> FingerprintGenerator {
    FingerprintGenerator {
        env: EnvGenerator::Morgan(MorganArguments {
            radius,
            only_nonzero_invariants: false,
            include_redundant_environments: include_redundant,
        }),
        args: FingerprintArguments {
            count_simulation: false,
            include_chirality: use_chirality,
            count_bounds: vec![1, 2, 4, 8],
            fp_size,
            num_bits_per_feature: 1,
        },
        atom_inv_gen: Some(AtomInvariantsGenerator::Morgan {
            include_ring_membership: true,
        }),
        bond_inv_gen: Some(BondInvariantsGenerator::Morgan {
            use_bond_types,
            use_chirality,
        }),
        output_64: false,
    }
}

/// Which legacy Morgan function to run.
#[derive(Clone, Copy, Debug)]
pub enum MorganKind {
    /// `getFingerprint` (unfolded; counts or presence)
    Unfolded { use_counts: bool },
    /// `getHashedFingerprint`
    Hashed(u32),
    /// `getFingerprintAsBitVect`
    BitVect(u32),
}

/// Result of a legacy Morgan call.
#[derive(Debug)]
pub enum MorganResult {
    Counts(SparseIntVect),
    Bits(ExplicitBitVect),
}

/// `MorganFingerprints::getFingerprint` / `getHashedFingerprint` / `getFingerprintAsBitVect`
#[allow(clippy::too_many_arguments)]
pub fn morgan_fingerprint(
    mol: &ROMol,
    radius: u32,
    kind: MorganKind,
    invariants: Option<&[u32]>,
    from_atoms: Option<&[u32]>,
    use_chirality: bool,
    use_bond_types: bool,
    bit_info: Option<&mut BitInfoMap>,
    include_redundant: bool,
) -> Result<MorganResult> {
    let fp_size = match kind {
        MorganKind::Unfolded { .. } => 2048,
        MorganKind::Hashed(n) | MorganKind::BitVect(n) => {
            if n == 0 {
                return Err(RDKitError::Value("nBits can not be zero".into()));
            }
            n
        }
    };
    let gen_ = morgan_generator32(
        radius,
        use_chirality,
        use_bond_types,
        include_redundant,
        fp_size,
    );
    let mut ao = AdditionalOutput {
        bit_info_map: bit_info.is_some().then(BTreeMap::new),
        ..Default::default()
    };
    let want = bit_info.is_some();
    let res = {
        let mut fa = FingerprintFuncArguments {
            from_atoms,
            custom_atom_invariants: invariants,
            additional_output: want.then_some(&mut ao),
            ..Default::default()
        };
        match kind {
            MorganKind::Unfolded { use_counts: false } => {
                let sbv = gen_.sparse_fingerprint(mol, &mut fa)?;
                let mut r = SparseIntVect::new(IndexKind::U32, i128::from(u32::MAX));
                for b in sbv.on_bits() {
                    r.set_val(i128::from(b as u32), 1)?;
                }
                MorganResult::Counts(r)
            }
            MorganKind::Unfolded { use_counts: true } => {
                MorganResult::Counts(gen_.sparse_count_fingerprint(mol, &mut fa)?)
            }
            MorganKind::Hashed(_) => MorganResult::Counts(gen_.count_fingerprint(mol, &mut fa)?),
            MorganKind::BitVect(_) => MorganResult::Bits(gen_.fingerprint(mol, &mut fa)?),
        }
    };
    if let Some(bi) = bit_info {
        *bi = ao.bit_info_map.unwrap_or_default();
    }
    Ok(res)
}

// ---------------------------------------------------------------------------
// AtomPairs.cpp

/// Copy of `mol` with stereo perceived when chirality is requested.
fn stereo_mol(mol: &ROMol, include_chirality: bool) -> Result<Option<ROMol>> {
    if include_chirality && !mol.has_prop(cp::_StereochemDone) {
        let mut m = mol.clone();
        crate::graphmol::chirality::legacy::assign_stereochemistry(&mut m, false, false, false)?;
        return Ok(Some(m));
    }
    Ok(None)
}

fn check_invariants(mol: &ROMol, inv: Option<&[u32]>) -> Result<()> {
    if inv.is_some_and(|v| v.len() < mol.num_atoms() as usize) {
        return Err(RDKitError::precondition("bad atomInvariants size"));
    }
    Ok(())
}

/// `getAtomPairFingerprintInternal` (32-bit generator); `sparse` selects the
/// unfolded sparse count fingerprint.
#[allow(clippy::too_many_arguments)]
fn atom_pair_internal(
    mol: &ROMol,
    n_bits: u32,
    min_length: u32,
    max_length: u32,
    from_atoms: Option<&[u32]>,
    ignore_atoms: Option<&[u32]>,
    atom_invariants: Option<&[u32]>,
    include_chirality: bool,
    use_2d: bool,
    conf_id: i32,
    sparse: bool,
) -> Result<SparseIntVect> {
    check_invariants(mol, atom_invariants)?;
    let sm = stereo_mol(mol, include_chirality)?;
    let lmol = sm.as_ref().unwrap_or(mol);
    if min_length > max_length {
        return Err(RDKitError::precondition("bad distances provided"));
    }
    let gen_ = FingerprintGenerator {
        env: EnvGenerator::AtomPair(AtomPairArguments {
            use_2d,
            min_distance: min_length,
            max_distance: max_length,
        }),
        args: FingerprintArguments {
            count_simulation: true,
            include_chirality,
            count_bounds: vec![1, 2, 4, 8],
            fp_size: n_bits,
            num_bits_per_feature: 1,
        },
        atom_inv_gen: Some(AtomInvariantsGenerator::AtomPair {
            include_chirality,
            topological_torsion_correction: false,
        }),
        bond_inv_gen: None,
        output_64: false,
    };
    let mut fa = FingerprintFuncArguments {
        from_atoms,
        ignore_atoms,
        custom_atom_invariants: atom_invariants,
        conf_id,
        ..Default::default()
    };
    if sparse {
        gen_.sparse_count_fingerprint(lmol, &mut fa)
    } else {
        gen_.count_fingerprint(lmol, &mut fa)
    }
}

/// Re-types a uint32 SparseIntVect as `SparseIntVect<int32_t>` (upstream's
/// `reinterpret_cast`).
fn as_i32(v: &SparseIntVect) -> Result<SparseIntVect> {
    let mut r = SparseIntVect::new(IndexKind::I32, i128::from(v.length() as u32 as i32));
    for (&k, &val) in v.nonzero_elements() {
        r.set_val(i128::from(k as u32 as i32), val)?;
    }
    Ok(r)
}

/// `getAtomPairFingerprint` (returns `SparseIntVect<int32_t>`)
#[allow(clippy::too_many_arguments)]
pub fn get_atom_pair_fingerprint(
    mol: &ROMol,
    min_length: u32,
    max_length: u32,
    from_atoms: Option<&[u32]>,
    ignore_atoms: Option<&[u32]>,
    atom_invariants: Option<&[u32]>,
    include_chirality: bool,
    use_2d: bool,
    conf_id: i32,
) -> Result<SparseIntVect> {
    let v = atom_pair_internal(
        mol,
        0,
        min_length,
        max_length,
        from_atoms,
        ignore_atoms,
        atom_invariants,
        include_chirality,
        use_2d,
        conf_id,
        true,
    )?;
    as_i32(&v)
}

/// `getHashedAtomPairFingerprint` (returns `SparseIntVect<int32_t>`)
#[allow(clippy::too_many_arguments)]
pub fn get_hashed_atom_pair_fingerprint(
    mol: &ROMol,
    n_bits: u32,
    min_length: u32,
    max_length: u32,
    from_atoms: Option<&[u32]>,
    ignore_atoms: Option<&[u32]>,
    atom_invariants: Option<&[u32]>,
    include_chirality: bool,
    use_2d: bool,
    conf_id: i32,
) -> Result<SparseIntVect> {
    let v = atom_pair_internal(
        mol,
        n_bits,
        min_length,
        max_length,
        from_atoms,
        ignore_atoms,
        atom_invariants,
        include_chirality,
        use_2d,
        conf_id,
        false,
    )?;
    let mut r = SparseIntVect::new(IndexKind::I32, i128::from(n_bits as i32));
    for (&k, &val) in v.nonzero_elements() {
        r.set_val(i128::from(k as u32 as i32), val)?;
    }
    Ok(r)
}

/// The count -> bits expansion of the `*AsBitVect` functions.
fn counts_to_bits(
    counts: &SparseIntVect,
    n_bits: u32,
    n_bits_per_entry: u32,
) -> Result<ExplicitBitVect> {
    const BOUNDS: [i32; 4] = [1, 2, 4, 8];
    let mut res = ExplicitBitVect::new(n_bits, false);
    for (&k, &val) in counts.nonzero_elements() {
        for i in 0..n_bits_per_entry {
            let set = if n_bits_per_entry != 4 {
                val > i as i32
            } else {
                val >= BOUNDS[i as usize]
            };
            if set {
                res.set_bit((k as u32).wrapping_mul(n_bits_per_entry).wrapping_add(i))?;
            }
        }
    }
    Ok(res)
}

/// `getHashedAtomPairFingerprintAsBitVect`
#[allow(clippy::too_many_arguments)]
pub fn get_hashed_atom_pair_fingerprint_as_bit_vect(
    mol: &ROMol,
    n_bits: u32,
    min_length: u32,
    max_length: u32,
    from_atoms: Option<&[u32]>,
    ignore_atoms: Option<&[u32]>,
    atom_invariants: Option<&[u32]>,
    n_bits_per_entry: u32,
    include_chirality: bool,
    use_2d: bool,
    conf_id: i32,
) -> Result<ExplicitBitVect> {
    if min_length > max_length {
        return Err(RDKitError::precondition("bad lengths provided"));
    }
    check_invariants(mol, atom_invariants)?;
    let block = n_bits
        .checked_div(n_bits_per_entry)
        .ok_or_else(|| RDKitError::Value("nBitsPerEntry can not be zero".into()))?;
    let sres = get_hashed_atom_pair_fingerprint(
        mol,
        block,
        min_length,
        max_length,
        from_atoms,
        ignore_atoms,
        atom_invariants,
        include_chirality,
        use_2d,
        conf_id,
    )?;
    counts_to_bits(&sres, n_bits, n_bits_per_entry)
}

/// `getTopologicalTorsionFingerprint` (returns `SparseIntVect<int64_t>`)
pub fn get_topological_torsion_fingerprint(
    mol: &ROMol,
    target_size: u32,
    from_atoms: Option<&[u32]>,
    ignore_atoms: Option<&[u32]>,
    atom_invariants: Option<&[u32]>,
    include_chirality: bool,
) -> Result<SparseIntVect> {
    check_invariants(mol, atom_invariants)?;
    let sm = stereo_mol(mol, include_chirality)?;
    let lmol = sm.as_ref().unwrap_or(mol);
    let bits = CODE_SIZE
        + if include_chirality {
            NUM_CHIRAL_BITS
        } else {
            0
        };
    let sz = 1u64
        .wrapping_shl(target_size.wrapping_mul(bits))
        .wrapping_sub(1);
    let mut res = SparseIntVect::new(IndexKind::I64, i128::from(sz as i64));
    let n = lmol.num_atoms() as usize;
    let atom_codes: Vec<u32> = match atom_invariants {
        None => (0..lmol.num_atoms())
            .map(|i| get_atom_code(lmol, i, 0, include_chirality))
            .collect::<Result<_>>()?,
        Some(inv) => inv[..n]
            .iter()
            .map(|&v| v % ((1 << CODE_SIZE) - 1) + 2)
            .collect(),
    };
    let to_set = |v: &[u32]| {
        let mut s = vec![false; n];
        for &a in v {
            if let Some(x) = s.get_mut(a as usize) {
                *x = true;
            }
        }
        s
    };
    let from = from_atoms.map(to_set);
    let ignore = ignore_atoms.map(to_set);
    let mut p_atoms = vec![false; n];
    for path in find_all_paths_of_length_n(lmol, target_size, false, false, -1, false, None)? {
        let mut keep = from.is_none();
        if let Some(f) = &from
            && (f[path[0] as usize] || f[*path.last().expect("non-empty") as usize])
        {
            keep = true;
        }
        if keep
            && let Some(ig) = &ignore
            && path.iter().any(|&a| ig[a as usize])
        {
            keep = false;
        }
        if !keep {
            continue;
        }
        p_atoms.iter_mut().for_each(|b| *b = false);
        let mut codes = Vec::new();
        for (k, &a) in path.iter().enumerate() {
            if k != 0 && a != path[0] && p_atoms[a as usize] {
                codes.clear();
                break;
            }
            p_atoms[a as usize] = true;
            let mut code = atom_codes[a as usize].wrapping_sub(1);
            if k != 0 && k + 1 != path.len() {
                code = code.wrapping_sub(1);
            }
            codes.push(code);
        }
        if !codes.is_empty() {
            let code = get_topological_torsion_code(&codes, include_chirality) as i64;
            let v = res.get_val(i128::from(code))?;
            res.set_val(i128::from(code), v.wrapping_add(1))?;
        }
    }
    Ok(res)
}

/// `TorsionFpCalc` with the 64-bit torsion generator; returns `SparseIntVect<int64_t>(n_bits)`.
fn torsion_fp_calc(
    mol: &ROMol,
    n_bits: u32,
    target_size: u32,
    from_atoms: Option<&[u32]>,
    ignore_atoms: Option<&[u32]>,
    atom_invariants: Option<&[u32]>,
    include_chirality: bool,
) -> Result<SparseIntVect> {
    check_invariants(mol, atom_invariants)?;
    let sm = stereo_mol(mol, include_chirality)?;
    let lmol = sm.as_ref().unwrap_or(mol);
    let gen_ = FingerprintGenerator {
        env: EnvGenerator::TopologicalTorsion(TopologicalTorsionArguments {
            torsion_atom_count: target_size,
            only_shortest_paths: false,
        }),
        args: FingerprintArguments {
            count_simulation: true,
            include_chirality,
            count_bounds: vec![1, 2, 4, 8],
            fp_size: n_bits,
            num_bits_per_feature: 1,
        },
        atom_inv_gen: Some(AtomInvariantsGenerator::AtomPair {
            include_chirality,
            topological_torsion_correction: true,
        }),
        bond_inv_gen: None,
        output_64: true,
    };
    let mut fa = FingerprintFuncArguments {
        from_atoms,
        ignore_atoms,
        custom_atom_invariants: atom_invariants,
        ..Default::default()
    };
    let siv = gen_.count_fingerprint(lmol, &mut fa)?;
    let mut res = SparseIntVect::new(IndexKind::I64, i128::from(n_bits));
    for (&k, &v) in siv.nonzero_elements() {
        res.set_val(k, v)?;
    }
    Ok(res)
}

/// `getHashedTopologicalTorsionFingerprint` (returns `SparseIntVect<int64_t>`)
pub fn get_hashed_topological_torsion_fingerprint(
    mol: &ROMol,
    n_bits: u32,
    target_size: u32,
    from_atoms: Option<&[u32]>,
    ignore_atoms: Option<&[u32]>,
    atom_invariants: Option<&[u32]>,
    include_chirality: bool,
) -> Result<SparseIntVect> {
    torsion_fp_calc(
        mol,
        n_bits,
        target_size,
        from_atoms,
        ignore_atoms,
        atom_invariants,
        include_chirality,
    )
}

/// `getHashedTopologicalTorsionFingerprintAsBitVect`
#[allow(clippy::too_many_arguments)]
pub fn get_hashed_topological_torsion_fingerprint_as_bit_vect(
    mol: &ROMol,
    n_bits: u32,
    target_size: u32,
    from_atoms: Option<&[u32]>,
    ignore_atoms: Option<&[u32]>,
    atom_invariants: Option<&[u32]>,
    n_bits_per_entry: u32,
    include_chirality: bool,
) -> Result<ExplicitBitVect> {
    let block = n_bits
        .checked_div(n_bits_per_entry)
        .ok_or_else(|| RDKitError::Value("nBitsPerEntry can not be zero".into()))?;
    let sres = torsion_fp_calc(
        mol,
        block,
        target_size,
        from_atoms,
        ignore_atoms,
        atom_invariants,
        include_chirality,
    )?;
    counts_to_bits(&sres, n_bits, n_bits_per_entry)
}

/// `maxPathLen - 1`, the default `maxLength` of the atom-pair functions
pub const DEFAULT_MAX_LENGTH: u32 = MAX_PATH_LEN - 1;
