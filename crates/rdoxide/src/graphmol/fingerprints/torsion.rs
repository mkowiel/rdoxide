//! Topological-torsion fingerprint generator (upstream
//! `GraphMol/Fingerprints/TopologicalTorsionGenerator.{h,cpp}` and the torsion
//! code helpers of `FingerprintUtil.cpp`).

use super::atom_pair::{CODE_SIZE, NUM_CHIRAL_BITS};
use super::{AdditionalOutput, hash_combine};
use crate::error::{RDKitError, Result};
use crate::graphmol::romol::ROMol;
use crate::graphmol::subgraphs::find_all_paths_of_length_n;

/// `TopologicalTorsionArguments` (the torsion-specific part).
#[derive(Clone, Debug)]
pub struct TopologicalTorsionArguments {
    pub torsion_atom_count: u32,
    pub only_shortest_paths: bool,
}

impl TopologicalTorsionArguments {
    pub fn info_string(&self) -> String {
        format!(
            "TopologicalTorsionArguments torsionAtomCount={} onlyShortestPaths={}",
            self.torsion_atom_count,
            u8::from(self.only_shortest_paths)
        )
    }
    /// `TopologicalTorsionEnvGenerator::getResultSize` (64-bit shift, wrapping
    /// the shift count like the hardware does)
    pub fn result_size(&self, include_chirality: bool) -> u64 {
        let bits = CODE_SIZE
            + if include_chirality {
                NUM_CHIRAL_BITS
            } else {
                0
            };
        1u64.wrapping_shl(self.torsion_atom_count.wrapping_mul(bits))
    }
}

/// `TopologicalTorsionAtomEnv`
#[derive(Clone, Debug)]
pub struct TopologicalTorsionAtomEnv {
    pub bit_id: u64,
    pub atom_path: Vec<i32>,
}

impl TopologicalTorsionAtomEnv {
    pub(super) fn update_additional_output(
        &self,
        ao: &mut AdditionalOutput,
        bit_id: u64,
    ) -> Result<()> {
        let oob = || RDKitError::Index("vector::_M_range_check".into());
        for &aid in &self.atom_path {
            if let Some(a) = &mut ao.atom_to_bits {
                a.get_mut(aid as usize).ok_or_else(oob)?.push(bit_id);
            }
            if let Some(c) = &mut ao.atom_counts {
                *c.get_mut(aid as usize).ok_or_else(oob)? += 1;
            }
        }
        if let Some(bp) = &mut ao.bit_paths {
            bp.entry(bit_id)
                .or_default()
                .push(self.atom_path.iter().map(|&a| a as u32).collect());
        }
        if let Some(apb) = &mut ao.atoms_per_bit {
            apb.entry(bit_id).or_default().push(self.atom_path.clone());
        }
        Ok(())
    }
}

/// Whether the path codes should be read back to front (smaller end first).
fn reverse_it(codes: &[u32]) -> bool {
    let (mut i, mut j) = (0usize, codes.len() - 1);
    while i < j {
        if codes[i] > codes[j] {
            return true;
        } else if codes[i] < codes[j] {
            return false;
        }
        i += 1;
        j -= 1;
    }
    false
}

/// `getTopologicalTorsionCode`
pub fn get_topological_torsion_code(codes: &[u32], include_chirality: bool) -> u64 {
    let shift = CODE_SIZE
        + if include_chirality {
            NUM_CHIRAL_BITS
        } else {
            0
        };
    let ordered: Vec<u32> = if reverse_it(codes) {
        codes.iter().rev().copied().collect()
    } else {
        codes.to_vec()
    };
    let mut res = 0u64;
    for (i, &c) in ordered.iter().enumerate() {
        res |= u64::from(c).wrapping_shl(shift * i as u32);
    }
    res
}

/// `getTopologicalTorsionHash`
pub fn get_topological_torsion_hash(codes: &[u32]) -> u32 {
    let mut res = 0u32;
    if reverse_it(codes) {
        for &c in codes.iter().rev() {
            hash_combine(&mut res, c);
        }
    } else {
        for &c in codes {
            hash_combine(&mut res, c);
        }
    }
    res
}

/// `TopologicalTorsionEnvGenerator::getEnvironments`
pub fn get_environments(
    mol: &ROMol,
    targs: &TopologicalTorsionArguments,
    include_chirality: bool,
    from_atoms: Option<&[u32]>,
    ignore_atoms: Option<&[u32]>,
    atom_invariants: &[u32],
    hash_results: bool,
) -> Result<Vec<TopologicalTorsionAtomEnv>> {
    let n = mol.num_atoms() as usize;
    let to_set = |v: &[u32]| -> Result<Vec<bool>> {
        let mut s = vec![false; n];
        for &a in v {
            *s.get_mut(a as usize)
                .ok_or_else(|| RDKitError::precondition("bad atom index"))? = true;
        }
        Ok(s)
    };
    let from = from_atoms.map(to_set).transpose()?;
    let ignore = ignore_atoms.map(to_set).transpose()?;
    let paths = find_all_paths_of_length_n(
        mol,
        targs.torsion_atom_count,
        false,
        false,
        -1,
        targs.only_shortest_paths,
        None,
    )?;
    let mut result = Vec::new();
    let mut p_atoms = vec![false; n];
    for path in paths {
        let mut keep = from.is_none();
        if let Some(f) = &from
            && (f[path[0] as usize] || f[*path.last().expect("non-empty path") as usize])
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
        let mut codes: Vec<u32> = Vec::new();
        for (k, &a) in path.iter().enumerate() {
            if k != 0 && a != path[0] && p_atoms[a as usize] {
                codes.clear();
                break;
            }
            p_atoms[a as usize] = true;
            let inv = *atom_invariants
                .get(a as usize)
                .ok_or_else(|| RDKitError::precondition("bad atom invariants size"))?;
            let mut code = inv % ((1u32 << CODE_SIZE) - 1) + 1;
            if k != 0 && k + 1 != path.len() {
                code -= 1;
            }
            codes.push(code);
        }
        if !codes.is_empty() {
            let bit_id = if hash_results {
                u64::from(get_topological_torsion_hash(&codes))
            } else {
                get_topological_torsion_code(&codes, include_chirality)
            };
            result.push(TopologicalTorsionAtomEnv {
                bit_id,
                atom_path: path,
            });
        }
    }
    Ok(result)
}
