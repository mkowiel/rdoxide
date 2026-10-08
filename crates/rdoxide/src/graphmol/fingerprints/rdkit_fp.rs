//! RDKit (path-based) fingerprint generator (upstream
//! `GraphMol/Fingerprints/RDKitFPGenerator.{h,cpp}` and `RDKitFPUtils` of
//! `FingerprintUtil.cpp`).

use super::{AdditionalOutput, hash_vect};
use crate::error::{RDKitError, Result};
use crate::graphmol::bond::BondType;
use crate::graphmol::queries::{is_complex_atom_query, is_complex_bond_query};
use crate::graphmol::romol::ROMol;
use crate::graphmol::subgraphs::{
    PathMap, find_all_paths_of_lengths_m_to_n, find_all_subgraphs_of_lengths_m_to_n,
};

use super::hash_combine;

/// `RDKitFPArguments` (the RDKit-FP-specific part).
#[derive(Clone, Debug)]
pub struct RDKitFPArguments {
    pub min_path: u32,
    pub max_path: u32,
    pub use_hs: bool,
    pub branched_paths: bool,
    pub use_bond_order: bool,
}

impl RDKitFPArguments {
    pub fn info_string(&self) -> String {
        format!(
            "RDKitFPArguments minPath={} maxPath={} useHs={} branchedPaths={} useBondOrder={}",
            self.min_path,
            self.max_path,
            u8::from(self.use_hs),
            u8::from(self.branched_paths),
            u8::from(self.use_bond_order)
        )
    }
}

/// `RDKitFPAtomInvGenerator::getAtomInvariants`
pub fn rdkit_fp_atom_invariants(mol: &ROMol) -> Vec<u32> {
    (0..mol.num_atoms())
        .map(|i| {
            let a = mol.atom(i);
            ((a.atomic_num() % 128) << 1) | u32::from(a.is_aromatic())
        })
        .collect()
}

/// `RDKitFPAtomEnv`
#[derive(Clone, Debug)]
pub struct RDKitFPAtomEnv {
    pub bit_id: u64,
    pub atoms_in_path: Vec<bool>,
    pub bond_path: Vec<i32>,
}

impl RDKitFPAtomEnv {
    pub(super) fn update_additional_output(
        &self,
        ao: &mut AdditionalOutput,
        bit_id: u64,
    ) -> Result<()> {
        let oob = || RDKitError::Index("vector::_M_range_check".into());
        if let Some(bp) = &mut ao.bit_paths {
            bp.entry(bit_id)
                .or_default()
                .push(self.bond_path.iter().map(|&b| b as u32).collect());
        }
        if ao.atom_to_bits.is_some() || ao.atom_counts.is_some() || ao.atoms_per_bit.is_some() {
            if let Some(apb) = &mut ao.atoms_per_bit {
                apb.entry(bit_id).or_default().push(Vec::new());
            }
            for (i, &inp) in self.atoms_in_path.iter().enumerate() {
                if !inp {
                    continue;
                }
                if let Some(apb) = &mut ao.atoms_per_bit {
                    apb.get_mut(&bit_id)
                        .and_then(|v| v.last_mut())
                        .expect("pushed above")
                        .push(i as i32);
                }
                if let Some(a) = &mut ao.atom_to_bits {
                    let l = a.get_mut(i).ok_or_else(oob)?;
                    if !l.contains(&bit_id) {
                        l.push(bit_id);
                    }
                }
                if let Some(c) = &mut ao.atom_counts {
                    *c.get_mut(i).ok_or_else(oob)? += 1;
                }
            }
        }
        Ok(())
    }
}

/// `RDKitFPUtils::enumerateAllPaths`
pub fn enumerate_all_paths(
    mol: &ROMol,
    from_atoms: Option<&[u32]>,
    branched_paths: bool,
    use_hs: bool,
    min_path: u32,
    max_path: u32,
    ignore_atoms: Option<&[bool]>,
) -> Result<PathMap> {
    let get = |root: i32| -> Result<PathMap> {
        if branched_paths {
            find_all_subgraphs_of_lengths_m_to_n(
                mol,
                min_path,
                max_path,
                use_hs,
                root,
                ignore_atoms,
            )
        } else {
            find_all_paths_of_lengths_m_to_n(
                mol,
                min_path,
                max_path,
                true,
                use_hs,
                root,
                false,
                ignore_atoms,
            )
        }
    };
    match from_atoms {
        None => get(-1),
        Some(fa) => {
            let mut all = PathMap::new();
            for &a in fa {
                for (len, paths) in get(a as i32)? {
                    // upstream inserts at the front of the existing list
                    let e = all.entry(len).or_default();
                    let mut v = paths;
                    v.append(e);
                    *e = v;
                }
            }
            Ok(all)
        }
    }
}

/// `RDKitFPUtils::identifyQueryBonds`
pub fn identify_query_bonds(mol: &ROMol) -> Vec<u8> {
    (0..mol.num_bonds())
        .map(|i| {
            let b = mol.bond(i);
            let mut q = 0u8;
            if is_complex_bond_query(b.query()) {
                q = 0x1;
            }
            if is_complex_atom_query(mol.atom(b.begin_atom_idx()).query()) {
                q |= 0x2;
            }
            if is_complex_atom_query(mol.atom(b.end_atom_idx()).query()) {
                q |= 0x4;
            }
            q
        })
        .collect()
}

/// `RDKitFPUtils::generateBondHashes`; fills `atoms_in_path`.
pub fn generate_bond_hashes(
    mol: &ROMol,
    atoms_in_path: &mut [bool],
    is_query_bond: &[u8],
    path: &[i32],
    use_bond_order: bool,
    atom_invariants: &[u32],
) -> Result<Vec<u32>> {
    if atom_invariants.len() < mol.num_atoms() as usize {
        return Err(RDKitError::precondition("bad atomInvariants size"));
    }
    atoms_in_path.iter_mut().for_each(|b| *b = false);
    let mut degrees = vec![0u32; mol.num_atoms() as usize];
    let mut query_in_path = false;
    for &b in path {
        if query_in_path {
            break;
        }
        let bi = mol.bond(b as u32);
        let (s, e) = (bi.begin_atom_idx() as usize, bi.end_atom_idx() as usize);
        degrees[s] += 1;
        degrees[e] += 1;
        atoms_in_path[s] = true;
        atoms_in_path[e] = true;
        if is_query_bond[b as usize] != 0 {
            query_in_path = true;
        }
    }
    if query_in_path {
        return Ok(Vec::new());
    }
    let mut bond_nbrs = vec![0u32; path.len()];
    let mut hashes = Vec::with_capacity(path.len() + 1);
    for i in 0..path.len() {
        let bi = mol.bond(path[i] as u32);
        let (bis, bie) = (bi.begin_atom_idx(), bi.end_atom_idx());
        for j in i + 1..path.len() {
            let bj = mol.bond(path[j] as u32);
            let (bjs, bje) = (bj.begin_atom_idx(), bj.end_atom_idx());
            if bis == bjs || bis == bje || bie == bjs || bie == bje {
                bond_nbrs[i] += 1;
                bond_nbrs[j] += 1;
            }
        }
        let mut a1 = atom_invariants[bis as usize];
        let mut a2 = atom_invariants[bie as usize];
        let mut d1 = degrees[bis as usize];
        let mut d2 = degrees[bie as usize];
        if a1 < a2 {
            std::mem::swap(&mut a1, &mut a2);
            std::mem::swap(&mut d1, &mut d2);
        } else if a1 == a2 && d1 < d2 {
            std::mem::swap(&mut d1, &mut d2);
        }
        let bond_hash = if use_bond_order {
            if bi.is_aromatic() || bi.bond_type() == BondType::Aromatic {
                BondType::Aromatic as u32
            } else {
                bi.bond_type() as u32
            }
        } else {
            1
        };
        let mut h = bond_nbrs[i];
        hash_combine(&mut h, bond_hash);
        hash_combine(&mut h, a1);
        hash_combine(&mut h, d1);
        hash_combine(&mut h, a2);
        hash_combine(&mut h, d2);
        hashes.push(h);
    }
    Ok(hashes)
}

/// `RDKitFPEnvGenerator::getEnvironments`
pub fn get_environments(
    mol: &ROMol,
    rargs: &RDKitFPArguments,
    from_atoms: Option<&[u32]>,
    ignore_atoms: Option<&[u32]>,
    atom_invariants: &[u32],
) -> Result<Vec<RDKitFPAtomEnv>> {
    let n = mol.num_atoms() as usize;
    let ignore = match ignore_atoms {
        Some(ia) => {
            let mut s = vec![false; n];
            for &a in ia {
                *s.get_mut(a as usize)
                    .ok_or_else(|| RDKitError::precondition("bad atom index"))? = true;
            }
            Some(s)
        }
        None => None,
    };
    let all = enumerate_all_paths(
        mol,
        from_atoms,
        rargs.branched_paths,
        rargs.use_hs,
        rargs.min_path,
        rargs.max_path,
        ignore.as_deref(),
    )?;
    let is_query = identify_query_bonds(mol);
    let mut atoms_in_path = vec![false; n];
    let mut result = Vec::new();
    for paths in all.values() {
        for path in paths {
            let mut hashes = generate_bond_hashes(
                mol,
                &mut atoms_in_path,
                &is_query,
                path,
                rargs.use_bond_order,
                atom_invariants,
            )?;
            if hashes.is_empty() {
                continue;
            }
            let seed = if path.len() > 1 {
                hashes.sort_unstable();
                hashes.push(atoms_in_path.iter().filter(|&&b| b).count() as u32);
                hash_vect(&hashes)
            } else {
                hashes[0]
            };
            result.push(RDKitFPAtomEnv {
                bit_id: u64::from(seed),
                atoms_in_path: atoms_in_path.clone(),
                bond_path: path.clone(),
            });
        }
    }
    Ok(result)
}
