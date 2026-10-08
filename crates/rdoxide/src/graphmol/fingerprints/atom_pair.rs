//! Atom-pair fingerprint generator (upstream `GraphMol/Fingerprints/AtomPairGenerator.{h,cpp}`
//! and the atom-code helpers of `FingerprintUtil.{h,cpp}`).

use super::{AdditionalOutput, FingerprintArguments, hash_combine};
use crate::error::{RDKitError, Result};
use crate::graphmol::atom::num_pi_electrons;
use crate::graphmol::matrices::{get_3d_distance_mat, get_distance_mat};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;

pub const NUM_TYPE_BITS: u32 = 4;
pub const ATOM_NUMBER_TYPES: [u32; 16] = [5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 53, 0];
pub const NUM_PI_BITS: u32 = 2;
pub const MAX_NUM_PI: u32 = (1 << NUM_PI_BITS) - 1;
pub const NUM_BRANCH_BITS: u32 = 3;
pub const MAX_NUM_BRANCHES: u32 = (1 << NUM_BRANCH_BITS) - 1;
pub const NUM_CHIRAL_BITS: u32 = 2;
pub const CODE_SIZE: u32 = NUM_TYPE_BITS + NUM_PI_BITS + NUM_BRANCH_BITS;
pub const NUM_PATH_BITS: u32 = 5;
pub const MAX_PATH_LEN: u32 = (1 << NUM_PATH_BITS) - 1;
pub const NUM_ATOM_PAIR_FINGERPRINT_BITS: u32 = NUM_PATH_BITS + 2 * CODE_SIZE;

/// `getAtomCode` (legacy stereo perception: uses existing `_CIPCode`s)
pub fn get_atom_code(
    mol: &ROMol,
    idx: u32,
    branch_subtract: u32,
    include_chirality: bool,
) -> Result<u32> {
    let atom = mol.atom(idx);
    let degree = mol.atom_neighbors(idx).count() as u32;
    let num_branches = degree.saturating_sub(branch_subtract);
    let mut code = num_branches % MAX_NUM_BRANCHES;
    let n_pi = num_pi_electrons(mol, idx)? % MAX_NUM_PI;
    code |= n_pi << NUM_BRANCH_BITS;

    let n_types = 1u32 << NUM_TYPE_BITS;
    let anum = atom.atomic_num();
    let mut type_idx = 0;
    while type_idx < n_types {
        // the C++ table has 15 initializers; the 16th element is zero
        let t = ATOM_NUMBER_TYPES[type_idx as usize];
        if t == anum {
            break;
        } else if t > anum {
            type_idx = n_types;
            break;
        }
        type_idx += 1;
    }
    if type_idx == n_types {
        type_idx -= 1;
    }
    code |= type_idx << (NUM_BRANCH_BITS + NUM_PI_BITS);
    if include_chirality && let Some(cip) = atom.get_prop(cp::_CIPCode) {
        let offset = NUM_BRANCH_BITS + NUM_PI_BITS + NUM_TYPE_BITS;
        match cip.to_rd_string().as_str() {
            "R" => code |= 1 << offset,
            "S" => code |= 2 << offset,
            _ => {}
        }
    }
    Ok(code)
}

/// `getAtomPairCode`
pub fn get_atom_pair_code(
    code_i: u32,
    code_j: u32,
    dist: u32,
    include_chirality: bool,
) -> Result<u32> {
    if dist >= MAX_PATH_LEN {
        return Err(RDKitError::precondition("dist too long"));
    }
    let chiral = if include_chirality {
        NUM_CHIRAL_BITS
    } else {
        0
    };
    let mut res = dist;
    res |= code_i.min(code_j) << NUM_PATH_BITS;
    res |= code_i.max(code_j) << (NUM_PATH_BITS + CODE_SIZE + chiral);
    Ok(res)
}

/// `AtomPairAtomInvGenerator::getAtomInvariants`
pub fn atom_pair_invariants(
    mol: &ROMol,
    include_chirality: bool,
    tt_correction: bool,
) -> Result<Vec<u32>> {
    (0..mol.num_atoms())
        .map(|i| {
            Ok(
                get_atom_code(mol, i, 0, include_chirality)?.wrapping_sub(if tt_correction {
                    2
                } else {
                    0
                }),
            )
        })
        .collect()
}

/// `AtomPairArguments` (the atom-pair-specific part).
#[derive(Clone, Debug)]
pub struct AtomPairArguments {
    pub use_2d: bool,
    pub min_distance: u32,
    pub max_distance: u32,
}

impl AtomPairArguments {
    pub fn info_string(&self) -> String {
        format!(
            "AtomPairArguments use2D={} minDistance={} maxDistance={}",
            u8::from(self.use_2d),
            self.min_distance,
            self.max_distance
        )
    }
}

/// `AtomPairAtomEnv`
#[derive(Clone, Debug)]
pub struct AtomPairAtomEnv {
    pub first: u32,
    pub second: u32,
    pub distance: u32,
}

impl AtomPairAtomEnv {
    pub(super) fn bit_id(
        &self,
        args: &FingerprintArguments,
        atom_invs: Option<&[u32]>,
        hash_results: bool,
    ) -> Result<u64> {
        let invs = atom_invs.ok_or_else(|| RDKitError::precondition("bad atom invariants size"))?;
        let get = |i: u32| {
            invs.get(i as usize)
                .copied()
                .ok_or_else(|| RDKitError::precondition("bad atom invariants size"))
        };
        let chiral = if args.include_chirality {
            NUM_CHIRAL_BITS
        } else {
            0
        };
        let limit = (1u32 << (CODE_SIZE + chiral)) - 1;
        let a = get(self.first)? % limit;
        let b = get(self.second)? % limit;
        let bit = if hash_results {
            let mut h = 0u32;
            hash_combine(&mut h, a.min(b));
            hash_combine(&mut h, self.distance);
            hash_combine(&mut h, a.max(b));
            h
        } else {
            get_atom_pair_code(a, b, self.distance, args.include_chirality)?
        };
        Ok(u64::from(bit))
    }

    pub(super) fn update_additional_output(
        &self,
        ao: &mut AdditionalOutput,
        bit_id: u64,
    ) -> Result<()> {
        let oob = || RDKitError::Index("vector::_M_range_check".into());
        if let Some(m) = &mut ao.bit_info_map {
            m.entry(bit_id).or_default().push((self.first, self.second));
        }
        if let Some(a) = &mut ao.atom_to_bits {
            a.get_mut(self.first as usize).ok_or_else(oob)?.push(bit_id);
            a.get_mut(self.second as usize)
                .ok_or_else(oob)?
                .push(bit_id);
        }
        if let Some(c) = &mut ao.atom_counts {
            *c.get_mut(self.first as usize).ok_or_else(oob)? += 1;
            *c.get_mut(self.second as usize).ok_or_else(oob)? += 1;
        }
        if let Some(apb) = &mut ao.atoms_per_bit {
            apb.entry(bit_id)
                .or_default()
                .push(vec![self.first as i32, self.second as i32]);
        }
        Ok(())
    }
}

/// `AtomPairEnvGenerator::getEnvironments`
pub fn get_environments(
    mol: &ROMol,
    aargs: &AtomPairArguments,
    from_atoms: Option<&[u32]>,
    ignore_atoms: Option<&[u32]>,
    conf_id: i32,
    additional_output: Option<&AdditionalOutput>,
) -> Result<Vec<AtomPairAtomEnv>> {
    let n = mol.num_atoms();
    if let Some(ao) = additional_output
        && let Some(a) = &ao.atom_to_bits
        && a.len() != n as usize
    {
        return Err(RDKitError::precondition(
            "bad atomToBits size in AdditionalOutput",
        ));
    }
    let dm = if aargs.use_2d {
        get_distance_mat(mol)
    } else {
        get_3d_distance_mat(mol, conf_id)?
    };
    let ignored = |i: u32| ignore_atoms.is_some_and(|ia| ia.contains(&i));
    let mut result = Vec::new();
    for i in 0..n {
        if ignored(i) {
            continue;
        }
        for j in i + 1..n {
            if ignored(j) {
                continue;
            }
            if let Some(fa) = from_atoms
                && !fa.contains(&i)
                && !fa.contains(&j)
            {
                continue;
            }
            let distance = dm[(i * n + j) as usize].floor() as u32;
            if distance >= aargs.min_distance && distance <= aargs.max_distance {
                result.push(AtomPairAtomEnv {
                    first: i,
                    second: j,
                    distance,
                });
            }
        }
    }
    Ok(result)
}
