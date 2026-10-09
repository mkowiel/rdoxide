//! Layered fingerprint (upstream `LayeredFingerprintMol` in
//! `GraphMol/Fingerprints/Fingerprints.cpp`).
//!
//! Every path/subgraph is hashed once per requested "layer"; each layer adds
//! one kind of information (topology, bond orders, atom types, ring
//! membership, ring sizes, aromaticity) so that screening can trade
//! selectivity for query tolerance.

use super::hash_vect;
use super::rdkit_fp::{enumerate_all_paths, identify_query_bonds};
use crate::datastructs::{BitVect, ExplicitBitVect};
use crate::error::{RDKitError, Result};
use crate::graphmol::bond::{Bond, BondType};
use crate::graphmol::find_rings::find_sssr;
use crate::graphmol::queries::is_atom_aromatic;
use crate::graphmol::romol::ROMol;

/// `LayeredFingerprintMolVersion`
pub const LAYERED_FINGERPRINT_MOL_VERSION: &str = "0.7.0";
/// `substructLayers`: the layers that are safe for substructure screening.
pub const SUBSTRUCT_LAYERS: u32 = 0x07;
/// `maxFingerprintLayers`
pub const MAX_FINGERPRINT_LAYERS: usize = 10;

/// Layers that are actually computed (the remaining ones are reserved).
const NUM_LAYERS: usize = 6;

/// Bits of `identify_query_bonds`: the bond itself / either atom is a query.
const QUERY_BOND: u8 = 0x1;
const QUERY_ATOMS: u8 = 0x6;

/// `LayeredFingerprintMol` arguments that do not change between paths.
#[derive(Debug, Clone)]
pub struct LayeredFingerprintArgs<'a> {
    pub layer_flags: u32,
    pub min_path: u32,
    pub max_path: u32,
    pub fp_size: u32,
    pub branched_paths: bool,
    pub from_atoms: Option<&'a [u32]>,
}

impl Default for LayeredFingerprintArgs<'_> {
    fn default() -> Self {
        LayeredFingerprintArgs {
            layer_flags: 0xFFFF_FFFF,
            min_path: 1,
            max_path: 7,
            fp_size: 2048,
            branched_paths: true,
            from_atoms: None,
        }
    }
}

fn shares_atom(a: &Bond, b: &Bond) -> bool {
    a.begin_atom_idx() == b.begin_atom_idx()
        || a.begin_atom_idx() == b.end_atom_idx()
        || a.end_atom_idx() == b.begin_atom_idx()
        || a.end_atom_idx() == b.end_atom_idx()
}

/// The molecule-wide data the per-bond layer hashes are computed from.
struct MolContext<'a> {
    mol: &'a ROMol,
    query_bond: Vec<u8>,
    aromatic_atoms: Vec<bool>,
}

/// The state of the path whose bonds are being hashed.
struct PathContext {
    atom_degrees: Vec<u32>,
    bond_nbrs: Vec<u32>,
    path_queries: u8,
}

impl MolContext<'_> {
    /// `(larger, smaller)` degree of the bond's atoms.
    fn sorted_degrees(&self, bond: &Bond, path: &PathContext) -> (u32, u32) {
        let d1 = path.atom_degrees[bond.begin_atom_idx() as usize];
        let d2 = path.atom_degrees[bond.end_atom_idx() as usize];
        (d1.max(d2), d1.min(d2))
    }

    /// Adds bond `bond` (number `i` of its path) to every requested layer.
    fn hash_bond(
        &self,
        layers: &mut [Vec<u32>; NUM_LAYERS],
        flags: u32,
        bond: &Bond,
        i: usize,
        path: &PathContext,
    ) -> Result<()> {
        let nbrs = path.bond_nbrs[i];
        let (a1_deg, a2_deg) = self.sorted_degrees(bond, path);
        let no_atom_queries = path.path_queries & QUERY_ATOMS == 0;
        if flags & 0x1 != 0 {
            // layer 1: straight topology
            layers[0].push((nbrs % 8) | (a1_deg % 8) << 3 | (a2_deg % 8) << 6);
        }
        if flags & 0x2 != 0 && path.path_queries & QUERY_BOND == 0 {
            // layer 2: bond orders; aromatic and single bonds hash alike
            let bond_hash = if !bond.is_aromatic()
                && bond.bond_type() != BondType::Single
                && bond.bond_type() != BondType::Aromatic
            {
                bond.bond_type() as u32
            } else {
                BondType::Single as u32
            };
            layers[1]
                .push((bond_hash % 8) | (nbrs % 8) << 3 | (a1_deg % 8) << 6 | (a2_deg % 8) << 9);
        }
        if flags & 0x4 != 0 && no_atom_queries {
            // layer 3: atom types
            let mol = self.mol;
            let anum = |idx: u32| mol.atom(idx).atomic_num() % 128;
            let mut a1 = (
                anum(bond.begin_atom_idx()),
                path.atom_degrees[bond.begin_atom_idx() as usize],
            );
            let mut a2 = (
                anum(bond.end_atom_idx()),
                path.atom_degrees[bond.end_atom_idx() as usize],
            );
            // order by atomic number, then degree
            if (a1.0, a1.1) < (a2.0, a2.1) {
                std::mem::swap(&mut a1, &mut a2);
            }
            layers[2]
                .push(a1.0 | a2.0 << 7 | (a1.1 % 8) << 14 | (a2.1 % 8) << 17 | (nbrs % 8) << 20);
        }
        let ring_info = self.mol.ring_info();
        if flags & 0x8 != 0 && no_atom_queries {
            // layer 4: ring membership
            if ring_info.num_bond_rings(bond.idx())? != 0 {
                layers[3].push(1);
            }
        }
        if flags & 0x10 != 0 && no_atom_queries {
            // layer 5: smallest ring size
            layers[4].push(ring_info.min_bond_ring_size(bond.idx())? % 8);
        }
        if flags & 0x20 != 0 && no_atom_queries {
            // layer 6: aromaticity of the end atoms (aromatic one first)
            let a1 = self.aromatic_atoms[bond.begin_atom_idx() as usize];
            let a2 = self.aromatic_atoms[bond.end_atom_idx() as usize];
            let (a1, a2) = if !a1 && a2 { (a2, a1) } else { (a1, a2) };
            layers[5].push(u32::from(a1) | u32::from(a2) << 1 | (nbrs % 8) << 5);
        }
        Ok(())
    }
}

fn layered_fingerprint(
    mol: &ROMol,
    args: &LayeredFingerprintArgs<'_>,
    mut atom_counts: Option<&mut [u32]>,
    set_only_bits: Option<&ExplicitBitVect>,
) -> Result<ExplicitBitVect> {
    let ctx = MolContext {
        mol,
        query_bond: identify_query_bonds(mol),
        aromatic_atoms: (0..mol.num_atoms())
            .map(|i| is_atom_aromatic(mol, i))
            .collect(),
    };
    let all_paths = enumerate_all_paths(
        mol,
        args.from_atoms,
        args.branched_paths,
        false,
        args.min_path,
        args.max_path,
        None,
    )?;

    let mut res = ExplicitBitVect::new(args.fp_size, false);
    let n_atoms = mol.num_atoms() as usize;
    let mut atom_degrees = vec![0u32; n_atoms];
    let mut atoms_in_path = vec![false; n_atoms];
    for path in all_paths.values().flatten() {
        let bonds: Vec<&Bond> = path.iter().map(|&b| mol.bond(b as u32)).collect();
        let mut pctx = PathContext {
            atom_degrees: std::mem::take(&mut atom_degrees),
            bond_nbrs: vec![0; path.len()],
            path_queries: path.iter().fold(0, |q, &b| q | ctx.query_bond[b as usize]),
        };
        pctx.atom_degrees.fill(0);
        atoms_in_path.fill(false);
        for bond in &bonds {
            for idx in [bond.begin_atom_idx(), bond.end_atom_idx()] {
                pctx.atom_degrees[idx as usize] += 1;
                atoms_in_path[idx as usize] = true;
            }
        }

        let mut layers: [Vec<u32>; NUM_LAYERS] = Default::default();
        for i in 0..bonds.len() {
            for j in i + 1..bonds.len() {
                if shares_atom(bonds[i], bonds[j]) {
                    pctx.bond_nbrs[i] += 1;
                    pctx.bond_nbrs[j] += 1;
                }
            }
            ctx.hash_bond(&mut layers, args.layer_flags, bonds[i], i, &pctx)?;
        }

        let n_path_atoms = atoms_in_path.iter().filter(|&&a| a).count() as u32;
        let mut flagged_path = false;
        for (l, layer) in layers.iter_mut().enumerate() {
            if layer.is_empty() {
                continue;
            }
            layer.sort_unstable();
            // the atom count distinguishes C1CC1 from CC(C)C
            layer.push(n_path_atoms);
            layer.push(l as u32 + 1);
            let bit_id = hash_vect(layer) % args.fp_size;
            if set_only_bits.is_none_or(|s| s.get_bit(bit_id).unwrap_or(false)) {
                res.set_bit(bit_id)?;
                if let Some(counts) = atom_counts.as_deref_mut().filter(|_| !flagged_path) {
                    for (count, _) in counts.iter_mut().zip(&atoms_in_path).filter(|(_, a)| **a) {
                        *count += 1;
                    }
                    flagged_path = true;
                }
            }
        }
        atom_degrees = pctx.atom_degrees;
    }
    Ok(res)
}

/// `LayeredFingerprintMol`
///
/// Divergence: for unbranched paths without `from_atoms` upstream calls
/// `findAllPathsOfLengthsMtoN(mol, minPath, maxPath, false)`, i.e. atom paths
/// whose atom indices are then used as bond indices (out-of-bounds reads that
/// crash the reference build). Here the bond paths that were evidently
/// intended are used instead.
///
/// `atom_counts`, when given, must have at least one entry per atom; it is
/// incremented for every atom of each path that contributes a new bit.
pub fn layered_fingerprint_mol(
    mol: &ROMol,
    args: &LayeredFingerprintArgs<'_>,
    atom_counts: Option<&mut [u32]>,
    set_only_bits: Option<&ExplicitBitVect>,
) -> Result<ExplicitBitVect> {
    if args.min_path == 0 {
        return Err(RDKitError::precondition("minPath==0"));
    }
    if args.max_path < args.min_path {
        return Err(RDKitError::precondition("maxPath<minPath"));
    }
    if args.fp_size == 0 {
        return Err(RDKitError::precondition("fpSize==0"));
    }
    if atom_counts
        .as_ref()
        .is_some_and(|c| c.len() < mol.num_atoms() as usize)
    {
        return Err(RDKitError::precondition("bad atomCounts size"));
    }
    if set_only_bits.is_some_and(|s| s.num_bits() != args.fp_size) {
        return Err(RDKitError::precondition("bad setOnlyBits size"));
    }
    if mol.ring_info().is_initialized() {
        layered_fingerprint(mol, args, atom_counts, set_only_bits)
    } else {
        // upstream perceives rings on the (mutable) ring info of a const mol
        let mut local = mol.clone();
        find_sssr(&mut local, false, false)?;
        layered_fingerprint(&local, args, atom_counts, set_only_bits)
    }
}
