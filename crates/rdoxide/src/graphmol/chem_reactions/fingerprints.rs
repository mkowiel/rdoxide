//! Ported from RDKit `Code/GraphMol/ChemReactions/ReactionFingerprints.cpp`:
//! structural (bit vector) and difference (count vector) fingerprints of
//! reactions.
//!
//! The template fingerprints of one role are independent of each other, so
//! they are computed in parallel (the result does not depend on the thread
//! count).

use super::{ChemicalReaction, ReactionMoleculeType, templates};
use crate::datastructs::sparse_int_vect::{IndexKind, SparseIntVect};
use crate::datastructs::{BitVect, ExplicitBitVect};
use crate::graphmol::find_rings::find_sssr;
use crate::graphmol::fingerprints::pattern::pattern_fingerprint_mol;
use crate::graphmol::fingerprints::{
    AtomInvariantsGenerator, BondInvariantsGenerator, EnvGenerator, FingerprintArguments,
    FingerprintFuncArguments, FingerprintGenerator, atom_pair::AtomPairArguments,
    morgan::MorganArguments, rdkit_fp::RDKitFPArguments, torsion::TopologicalTorsionArguments,
};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::threads::par_map;
use crate::{RDKitError, Result};

/// `FingerprintType`; the discriminants match upstream's.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum FingerprintType {
    AtomPair = 1,
    TopologicalTorsion = 2,
    Morgan = 3,
    RDKit = 4,
    Pattern = 5,
}

impl FingerprintType {
    pub fn from_i32(v: i32) -> Option<Self> {
        Some(match v {
            1 => Self::AtomPair,
            2 => Self::TopologicalTorsion,
            3 => Self::Morgan,
            4 => Self::RDKit,
            5 => Self::Pattern,
            _ => return None,
        })
    }
}

/// `ReactionFingerprintParams`
#[derive(Debug, Clone, PartialEq)]
pub struct ReactionFingerprintParams {
    /// include the agents of a reaction
    pub include_agents: bool,
    /// structural fingerprints: share of the bits given to the agents
    pub bit_ratio_agents: f64,
    /// difference fingerprints: weight of reactants and products
    pub non_agent_weight: u32,
    /// difference fingerprints: weight of the agents
    pub agent_weight: i32,
    pub fp_size: u32,
    pub fp_type: FingerprintType,
}

impl Default for ReactionFingerprintParams {
    fn default() -> Self {
        ReactionFingerprintParams {
            include_agents: false,
            bit_ratio_agents: 0.2,
            non_agent_weight: 10,
            agent_weight: 1,
            fp_size: 2048,
            fp_type: FingerprintType::AtomPair,
        }
    }
}

impl ReactionFingerprintParams {
    /// `DefaultStructuralFPParams`
    pub fn default_structural() -> Self {
        ReactionFingerprintParams {
            include_agents: true,
            bit_ratio_agents: 0.2,
            non_agent_weight: 1,
            agent_weight: 1,
            fp_size: 4096,
            fp_type: FingerprintType::Pattern,
        }
    }

    /// `DefaultDifferenceFPParams`
    pub fn default_difference() -> Self {
        ReactionFingerprintParams {
            include_agents: true,
            bit_ratio_agents: 0.0,
            non_agent_weight: 10,
            agent_weight: 1,
            fp_size: 2048,
            fp_type: FingerprintType::AtomPair,
        }
    }
}

fn unsupported() -> RDKitError {
    RDKitError::ChemicalReaction(">> unsupported fingerprint type".into())
}

/// The generator `generateFingerprint*` builds for `fp_type`: the default
/// arguments of the type, with the given size.
fn generator(fp_type: FingerprintType, fp_size: u32) -> Result<FingerprintGenerator> {
    let args = |count_simulation: bool, num_bits_per_feature: u32| FingerprintArguments {
        count_simulation,
        include_chirality: false,
        count_bounds: vec![1, 2, 4, 8],
        fp_size,
        num_bits_per_feature,
    };
    Ok(match fp_type {
        FingerprintType::AtomPair => FingerprintGenerator {
            env: EnvGenerator::AtomPair(AtomPairArguments {
                use_2d: true,
                min_distance: 1,
                max_distance: 30,
            }),
            args: args(true, 1),
            atom_inv_gen: Some(AtomInvariantsGenerator::AtomPair {
                include_chirality: false,
                topological_torsion_correction: false,
            }),
            bond_inv_gen: None,
            output_64: false,
        },
        FingerprintType::TopologicalTorsion => FingerprintGenerator {
            env: EnvGenerator::TopologicalTorsion(TopologicalTorsionArguments {
                torsion_atom_count: 4,
                only_shortest_paths: false,
            }),
            args: args(true, 1),
            atom_inv_gen: Some(AtomInvariantsGenerator::AtomPair {
                include_chirality: false,
                topological_torsion_correction: true,
            }),
            bond_inv_gen: None,
            output_64: true,
        },
        FingerprintType::Morgan => FingerprintGenerator {
            env: EnvGenerator::Morgan(MorganArguments {
                radius: 2,
                only_nonzero_invariants: false,
                include_redundant_environments: false,
            }),
            args: args(false, 1),
            atom_inv_gen: Some(AtomInvariantsGenerator::Morgan {
                include_ring_membership: true,
            }),
            bond_inv_gen: Some(BondInvariantsGenerator::Morgan {
                use_bond_types: true,
                use_chirality: false,
            }),
            output_64: false,
        },
        FingerprintType::RDKit => FingerprintGenerator {
            env: EnvGenerator::RDKitFP(RDKitFPArguments {
                min_path: 1,
                max_path: 7,
                use_hs: true,
                branched_paths: true,
                use_bond_order: true,
            }),
            args: args(false, 2),
            atom_inv_gen: Some(AtomInvariantsGenerator::RDKitFP),
            bond_inv_gen: None,
            output_64: false,
        },
        FingerprintType::Pattern => return Err(unsupported()),
    })
}

/// A copy of `mol` with the property cache and (for Morgan) the rings that
/// the generators need.
fn prepared(mol: &ROMol, fp_type: FingerprintType) -> Result<ROMol> {
    let mut mol = mol.clone();
    mol.update_property_cache(false)?;
    if fp_type == FingerprintType::Morgan && !mol.ring_info().is_initialized() {
        mol.update_property_cache(true)?;
        find_sssr(&mut mol, false, false)?;
    }
    Ok(mol)
}

fn count_fingerprint(mol: &ROMol, fp_size: u32, fp_type: FingerprintType) -> Result<SparseIntVect> {
    if !matches!(
        fp_type,
        FingerprintType::AtomPair | FingerprintType::TopologicalTorsion | FingerprintType::Morgan
    ) {
        return Err(unsupported());
    }
    let mol = prepared(mol, fp_type)?;
    generator(fp_type, fp_size)?.count_fingerprint(&mol, &mut FingerprintFuncArguments::default())
}

fn bit_fingerprint(mol: &ROMol, fp_size: u32, fp_type: FingerprintType) -> Result<ExplicitBitVect> {
    let mol = prepared(mol, fp_type)?;
    if fp_type == FingerprintType::Pattern {
        return pattern_fingerprint_mol(&mol, fp_size, false);
    }
    generator(fp_type, fp_size)?.fingerprint(&mol, &mut FingerprintFuncArguments::default())
}

fn check_size(fp_size: u32) -> Result<()> {
    if fp_size == 0 {
        return Err(RDKitError::precondition("fpSize==0"));
    }
    Ok(())
}

/// `generateFingerprintChemReactionAsCountVect`: the sum of the count
/// fingerprints of the templates of one role.
pub fn count_fingerprint_of_templates(
    rxn: &ChemicalReaction,
    fp_size: u32,
    fp_type: FingerprintType,
    role: ReactionMoleculeType,
) -> Result<SparseIntVect> {
    check_size(fp_size)?;
    let mut result = SparseIntVect::new(IndexKind::U32, i128::from(fp_size));
    for fp in par_map(templates(rxn, role), 0, |m| {
        count_fingerprint(m, fp_size, fp_type)
    }) {
        result.add_assign(&fp?)?;
    }
    Ok(result)
}

/// `generateFingerprintChemReactionAsBitVect`: the union of the bit
/// fingerprints of the templates of one role.
pub fn bit_fingerprint_of_templates(
    rxn: &ChemicalReaction,
    fp_size: u32,
    fp_type: FingerprintType,
    role: ReactionMoleculeType,
) -> Result<ExplicitBitVect> {
    check_size(fp_size)?;
    let mut result = ExplicitBitVect::new(fp_size, false);
    for fp in par_map(templates(rxn, role), 0, |m| {
        bit_fingerprint(m, fp_size, fp_type)
    }) {
        result = result.or(&fp?);
    }
    Ok(result)
}

/// `StructuralFingerprintChemReaction`: reactant bits, product bits and
/// (optionally) agent bits, concatenated.
pub fn structural_fingerprint(
    rxn: &ChemicalReaction,
    params: &ReactionFingerprintParams,
) -> Result<ExplicitBitVect> {
    check_size(params.fp_size)?;
    let mut size_final = params.fp_size / 2;
    if params.include_agents {
        let mut agent_size = params.fp_size / 3;
        if params.bit_ratio_agents < 1.0 {
            agent_size = (f64::from(params.fp_size) * params.bit_ratio_agents).ceil() as u32;
        }
        let scaling = if agent_size % 2 == 0 {
            agent_size
        } else {
            agent_size - 1
        };
        size_final = (params.fp_size - scaling) / 2;
    }
    let size_agent = params.fp_size - 2 * size_final;
    let role_fp = |role, size| bit_fingerprint_of_templates(rxn, size, params.fp_type, role);
    let mut res = role_fp(ReactionMoleculeType::Reactant, size_final)?;
    res.append(&role_fp(ReactionMoleculeType::Product, size_final)?);
    if size_agent > 0 {
        res.append(&role_fp(ReactionMoleculeType::Agent, size_agent)?);
    }
    Ok(res)
}

/// `DifferenceFingerprintChemReaction`: products minus reactants (plus the
/// agents), with the weights of `params`.
pub fn difference_fingerprint(
    rxn: &ChemicalReaction,
    params: &ReactionFingerprintParams,
) -> Result<SparseIntVect> {
    check_size(params.fp_size)?;
    if !(1..4).contains(&(params.fp_type as i32)) {
        return Err(RDKitError::precondition("Fingerprinttype not supported"));
    }
    let role_fp = |role| count_fingerprint_of_templates(rxn, params.fp_size, params.fp_type, role);
    let mut res = role_fp(ReactionMoleculeType::Product)?;
    let mut reactants = role_fp(ReactionMoleculeType::Reactant)?;
    if params.include_agents {
        let mut agents = role_fp(ReactionMoleculeType::Agent)?;
        res.scalar_op('*', params.non_agent_weight as i32)?;
        reactants.scalar_op('*', params.non_agent_weight as i32)?;
        agents.scalar_op('*', params.agent_weight)?;
        res.sub_assign(&reactants)?;
        res.add_assign(&agents)?;
    } else {
        res.sub_assign(&reactants)?;
    }
    Ok(res)
}
