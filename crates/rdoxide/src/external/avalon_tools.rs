//! Avalon fingerprints for RDKit molecules (upstream `External/AvalonTools/AvalonTools.cpp`).
//! The molecule is written as a MOL block and read by the Avalon reader, as upstream does.

use crate::datastructs::bit_vect::{BitVect, ExplicitBitVect};
use crate::external::avalon::fingerprint::{
    ACCUMULATE_BITS, USE_ALL_FEATURES, USE_DY_AROMATICITY, set_fingerprint_counts_with_focus,
};
use crate::external::avalon::reaccs::Molecule;
use crate::external::avalon::reaccsio::mol_str_to_mol;
use crate::graphmol::file_parsers::mol_file_writer::{MolWriterParams, mol_to_mol_block};
use crate::graphmol::romol::ROMol;
use crate::{RDKitError, Result};

/// the default `bitFlags` of `GetAvalonFP` (all SSS features plus the non-SSS ones)
pub const DEFAULT_BIT_FLAGS: i32 = 15761407;

fn molblock_to_reaccs(molblock: &str) -> Result<Molecule> {
    mol_str_to_mol(molblock).ok_or_else(|| {
        RDKitError::Value("could not convert the molecule to the Avalon format".into())
    })
}

fn mol_to_reaccs(mol: &ROMol) -> Result<Molecule> {
    molblock_to_reaccs(&mol_to_mol_block(mol, &MolWriterParams::default(), -1)?)
}

/// `getFp`: the bytes of the fingerprint. Unless it is a query, a second pass with Daylight style
/// aromaticity is OR-ed in (the molecule keeps the changes of the first pass, as upstream).
fn fingerprint_bytes(
    mp: &mut Molecule,
    n_bytes: usize,
    is_query: bool,
    bit_flags: i32,
) -> Result<Vec<u8>> {
    let n_bytes = n_bytes.next_multiple_of(4);
    let mut bytes = vec![0u8; n_bytes];
    let passes: &[(bool, i32)] = if is_query {
        &[(true, 0)]
    } else {
        &[(false, 0), (false, ACCUMULATE_BITS | USE_DY_AROMATICITY)]
    };
    for &(as_query, fpflags) in passes {
        let mut counts = vec![0i32; n_bytes * 8];
        set_fingerprint_counts_with_focus(mp, &mut counts, bit_flags, as_query, fpflags, 0)
            .map_err(|e| RDKitError::Value(e.into()))?;
        for (i, &c) in counts.iter().enumerate() {
            if c > 0 {
                bytes[(i / 8) % n_bytes] |= 1 << (i % 8);
            }
        }
    }
    Ok(bytes)
}

/// `getAvalonFP`: sets the bits of `res` (which is cleared first if `reset_vect`)
pub fn get_avalon_fp_into(
    mol: &ROMol,
    res: &mut ExplicitBitVect,
    n_bits: u32,
    is_query: bool,
    reset_vect: bool,
    bit_flags: i32,
) -> Result<()> {
    let mut mp = mol_to_reaccs(mol)?;
    let bytes = fingerprint_bytes(&mut mp, (n_bits / 8) as usize, is_query, bit_flags)?;
    if reset_vect {
        *res = ExplicitBitVect::new(n_bits, false);
    }
    for (i, b) in bytes.iter().enumerate() {
        for j in 0..8 {
            if b & (1 << j) != 0 {
                res.set_bit((i * 8 + j) as u32)?;
            }
        }
    }
    Ok(())
}

/// `GetAvalonFP`
pub fn get_avalon_fp(
    mol: &ROMol,
    n_bits: u32,
    is_query: bool,
    bit_flags: i32,
) -> Result<ExplicitBitVect> {
    let mut res = ExplicitBitVect::new(n_bits, false);
    get_avalon_fp_into(mol, &mut res, n_bits, is_query, false, bit_flags)?;
    Ok(res)
}

/// `GetAvalonFP` of a MOL block
pub fn get_avalon_fp_from_molblock(
    molblock: &str,
    n_bits: u32,
    is_query: bool,
    bit_flags: i32,
) -> Result<ExplicitBitVect> {
    let mut mp = molblock_to_reaccs(molblock)?;
    let bytes = fingerprint_bytes(&mut mp, (n_bits / 8) as usize, is_query, bit_flags)?;
    let mut res = ExplicitBitVect::new(n_bits, false);
    for (i, b) in bytes.iter().enumerate() {
        for j in 0..8 {
            if b & (1 << j) != 0 {
                res.set_bit((i * 8 + j) as u32)?;
            }
        }
    }
    Ok(res)
}

/// `GetAvalonCountFP`: the feature multiplicities per bit position (`getCountFp`).
///
/// Upstream does not forward `bitFlags` to the counting code, which always uses `USE_ALL_FEATURES`
/// (32767); the argument is accepted and ignored here as well.
pub fn get_avalon_count_fp(
    mol: &ROMol,
    n_bits: u32,
    is_query: bool,
    bit_flags: i32,
) -> Result<Vec<u32>> {
    let _ = bit_flags;
    let bit_flags = USE_ALL_FEATURES;
    let mut mp = mol_to_reaccs(mol)?;
    let mut counts = vec![0i32; n_bits as usize];
    set_fingerprint_counts_with_focus(&mut mp, &mut counts, bit_flags, is_query, 0, 0)
        .map_err(|e| RDKitError::Value(e.into()))?;
    if !is_query {
        set_fingerprint_counts_with_focus(
            &mut mp,
            &mut counts,
            bit_flags,
            true,
            ACCUMULATE_BITS | USE_DY_AROMATICITY,
            0,
        )
        .map_err(|e| RDKitError::Value(e.into()))?;
    }
    Ok(counts.iter().map(|&c| c.max(0) as u32).collect())
}

/// `GetAvalonFPAsWords`: the fingerprint as 32 bit words
pub fn get_avalon_fp_as_words(
    mol: &ROMol,
    n_bits: u32,
    is_query: bool,
    bit_flags: i32,
) -> Result<Vec<u32>> {
    let mut mp = mol_to_reaccs(mol)?;
    let bytes = fingerprint_bytes(&mut mp, (n_bits / 8) as usize, is_query, bit_flags)?;
    Ok(bytes
        .chunks(4)
        .map(|c| {
            // upstream packs signed `char`s, so bytes >= 0x80 sign-extend into the higher bytes
            let sx = |i: usize| i32::from(c.get(i).copied().unwrap_or(0) as i8);
            (sx(0) | (sx(1) << 8) | (sx(2) << 16) | (sx(3) << 24)) as u32
        })
        .collect())
}
