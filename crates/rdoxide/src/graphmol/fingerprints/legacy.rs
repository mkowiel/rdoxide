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
