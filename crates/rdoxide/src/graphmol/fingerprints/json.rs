//! JSON (de)serialization of fingerprint generators (upstream `toJSON` /
//! `fromJSON` of `FingerprintGenerator`, its arguments, environment and
//! invariants generators, and `generatorToJSON` / `generatorFromJSON`).
//!
//! Unlike upstream, which dereferences a null pointer when the arguments or
//! the environment generator are missing from the JSON, those cases are
//! reported as `ValueError`s. The arguments and the environment generator
//! must also describe the same fingerprint type.

use super::atom_pair::AtomPairArguments;
use super::morgan::MorganArguments;
use super::rdkit_fp::RDKitFPArguments;
use super::torsion::TopologicalTorsionArguments;
use super::{
    AtomInvariantsGenerator, BondInvariantsGenerator, EnvGenerator, FingerprintArguments,
    FingerprintGenerator,
};
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::smiles_parse::smarts_write::mol_to_smarts;
use crate::graphmol::smiles_parse::write::SmilesWriteParams;
use crate::rdgeneral::ptree::PTree;
use crate::{RDKitError, Result};

fn value_error<T>(message: impl Into<String>) -> Result<T> {
    Err(RDKitError::Value(message.into()))
}

// ------------------------------------------------------------------ writing

fn common_arguments_to_json(args: &FingerprintArguments, pt: &mut PTree) {
    pt.put("countSimulation", args.count_simulation);
    pt.put("fpSize", args.fp_size);
    pt.put("numBitsPerFeature", args.num_bits_per_feature);
    pt.put("includeChirality", args.include_chirality);
    let mut bounds = PTree::new();
    for b in &args.count_bounds {
        bounds.push(PTree::leaf(b));
    }
    pt.add_child("countBounds", bounds);
}

fn env_to_json(env: &EnvGenerator, args: &FingerprintArguments) -> (PTree, PTree) {
    let (mut a, mut e) = (PTree::new(), PTree::new());
    match env {
        EnvGenerator::Morgan(m) => {
            a.put("type", "MorganArguments");
            a.put("onlyNonzeroInvariants", m.only_nonzero_invariants);
            a.put("radius", m.radius);
            e.put("type", "MorganEnvGenerator");
        }
        EnvGenerator::AtomPair(p) => {
            a.put("type", "AtomPairArguments");
            a.put("use2D", p.use_2d);
            a.put("minDistance", p.min_distance);
            a.put("maxDistance", p.max_distance);
            e.put("type", "AtomPairEnvGenerator");
        }
        EnvGenerator::TopologicalTorsion(t) => {
            a.put("type", "TopologicalTorsionArguments");
            a.put("torsionAtomCount", t.torsion_atom_count);
            a.put("onlyShortestPaths", t.only_shortest_paths);
            e.put("type", "TopologicalTorsionEnvGenerator");
        }
        EnvGenerator::RDKitFP(r) => {
            a.put("type", "RDKitFPArguments");
            a.put("minPath", r.min_path);
            a.put("maxPath", r.max_path);
            a.put("useHs", r.use_hs);
            a.put("branchedPaths", r.branched_paths);
            a.put("useBondOrder", r.use_bond_order);
            e.put("type", "RDKitFPEnvGenerator");
        }
    }
    common_arguments_to_json(args, &mut a);
    (a, e)
}

fn atom_inv_to_json(gen_: &AtomInvariantsGenerator) -> Result<PTree> {
    let mut pt = PTree::new();
    match gen_ {
        AtomInvariantsGenerator::Morgan {
            include_ring_membership,
        } => {
            pt.put("type", "MorganAtomInvGenerator");
            pt.put("includeRingMembership", include_ring_membership);
        }
        AtomInvariantsGenerator::AtomPair {
            include_chirality,
            topological_torsion_correction,
        } => {
            pt.put("type", "AtomPairAtomInvGenerator");
            pt.put("includeChirality", include_chirality);
            pt.put(
                "topologicalTorsionCorrection",
                topological_torsion_correction,
            );
        }
        AtomInvariantsGenerator::RDKitFP => pt.put("type", "RDKitFPAtomInvGenerator"),
        AtomInvariantsGenerator::MorganFeature { patterns } => {
            pt.put("type", "MorganFeatureAtomInvGenerator");
            if let Some(patterns) = patterns {
                let mut smarts = PTree::new();
                for p in patterns {
                    smarts.push(PTree::leaf(mol_to_smarts(
                        &mut p.clone(),
                        &SmilesWriteParams::default(),
                    )?));
                }
                pt.add_child("patternSMARTS", smarts);
            }
        }
    }
    Ok(pt)
}

fn bond_inv_to_json(gen_: &BondInvariantsGenerator) -> PTree {
    let BondInvariantsGenerator::Morgan {
        use_bond_types,
        use_chirality,
    } = gen_;
    let mut pt = PTree::new();
    pt.put("type", "MorganBondInvGenerator");
    pt.put("useBondTypes", use_bond_types);
    pt.put("useChirality", use_chirality);
    pt
}

/// `generatorToJSON`
pub fn generator_to_json(generator: &FingerprintGenerator) -> Result<String> {
    let mut pt = PTree::new();
    pt.put("name", "FingerprintGenerator");
    let (args, env) = env_to_json(&generator.env, &generator.args);
    pt.add_child("fingerprintArguments", args);
    pt.add_child("atomEnvironmentGenerator", env);
    if let Some(g) = &generator.atom_inv_gen {
        pt.add_child("atomInvariantsGenerator", atom_inv_to_json(g)?);
    }
    if let Some(g) = &generator.bond_inv_gen {
        pt.add_child("bondInvariantsGenerator", bond_inv_to_json(g));
    }
    Ok(pt.to_json())
}

// ------------------------------------------------------------------ reading

/// The constructor defaults of the argument classes a JSON document starts
/// from before its values are applied.
fn read_common(pt: &PTree, defaults: FingerprintArguments) -> FingerprintArguments {
    // 'countBounds' is replaced by whatever the document lists (possibly nothing)
    let count_bounds = pt
        .child("countBounds")
        .map(|c| {
            c.children()
                .map(|b| b.value().trim().parse().unwrap_or(0))
                .collect()
        })
        .unwrap_or_default();
    FingerprintArguments {
        count_simulation: pt.get_bool("countSimulation", defaults.count_simulation),
        include_chirality: pt.get_bool("includeChirality", defaults.include_chirality),
        count_bounds,
        fp_size: pt.get_u32("fpSize", defaults.fp_size),
        num_bits_per_feature: pt.get_u32("numBitsPerFeature", defaults.num_bits_per_feature),
    }
}

fn typed<'a>(node: &'a PTree, what: &str) -> Result<&'a str> {
    node.get_str("type").map_or_else(
        || value_error(format!("{what} type not specified in JSON")),
        Ok,
    )
}

fn default_common(count_simulation: bool, num_bits_per_feature: u32) -> FingerprintArguments {
    FingerprintArguments {
        count_simulation,
        include_chirality: false,
        count_bounds: vec![1, 2, 4, 8],
        fp_size: 2048,
        num_bits_per_feature,
    }
}

fn read_arguments(pt: &PTree) -> Result<(String, EnvGenerator, FingerprintArguments)> {
    let kind = typed(pt, "FingerprintArguments")?.to_string();
    let (env, args) = match kind.as_str() {
        "MorganArguments" => (
            EnvGenerator::Morgan(MorganArguments {
                radius: pt.get_u32("radius", 3),
                only_nonzero_invariants: pt.get_bool("onlyNonzeroInvariants", false),
                include_redundant_environments: false,
            }),
            read_common(pt, default_common(false, 1)),
        ),
        "RDKitFPArguments" => (
            EnvGenerator::RDKitFP(RDKitFPArguments {
                min_path: pt.get_u32("minPath", 1),
                max_path: pt.get_u32("maxPath", 7),
                use_hs: pt.get_bool("useHs", true),
                branched_paths: pt.get_bool("branchedPaths", true),
                use_bond_order: pt.get_bool("useBondOrder", true),
            }),
            read_common(pt, default_common(false, 2)),
        ),
        "AtomPairArguments" => (
            EnvGenerator::AtomPair(AtomPairArguments {
                use_2d: pt.get_bool("use2D", true),
                min_distance: pt.get_u32("minDistance", 1),
                max_distance: pt.get_u32("maxDistance", 30),
            }),
            read_common(pt, default_common(true, 1)),
        ),
        "TopologicalTorsionArguments" => (
            EnvGenerator::TopologicalTorsion(TopologicalTorsionArguments {
                torsion_atom_count: pt.get_u32("torsionAtomCount", 4),
                only_shortest_paths: pt.get_bool("onlyShortestPaths", false),
            }),
            read_common(pt, default_common(true, 1)),
        ),
        other => return value_error(format!("Unknown FingerprintArguments type: {other}")),
    };
    Ok((kind, env, args))
}

/// The `type` an environment generator document must carry for `args_kind`.
fn check_env(pt: &PTree, args_kind: &str) -> Result<()> {
    let typ = typed(pt, "AtomEnvironmentGenerator")?;
    let expected = match typ {
        "MorganEnvGenerator" => "MorganArguments",
        "RDKitFPEnvGenerator" => "RDKitFPArguments",
        "AtomPairEnvGenerator" => "AtomPairArguments",
        "TopologicalTorsionEnvGenerator" => "TopologicalTorsionArguments",
        other => return value_error(format!("Unknown AtomEnvGenerator type: {other}")),
    };
    if expected != args_kind {
        return value_error(format!(
            "{typ} does not match the fingerprint arguments of type {args_kind}"
        ));
    }
    Ok(())
}

fn read_atom_inv(pt: &PTree) -> Result<AtomInvariantsGenerator> {
    Ok(match typed(pt, "AtomInvariantsGenerator")? {
        "MorganAtomInvGenerator" => AtomInvariantsGenerator::Morgan {
            include_ring_membership: pt.get_bool("includeRingMembership", true),
        },
        "MorganFeatureAtomInvGenerator" => AtomInvariantsGenerator::MorganFeature {
            patterns: pt.child("patternSMARTS").map(|nodes| {
                // patterns that fail to parse are skipped, as upstream does
                nodes
                    .children()
                    .filter_map(|n| {
                        mol_from_smarts(n.value(), &SmartsParserParams::default())
                            .ok()
                            .flatten()
                    })
                    .collect::<Vec<ROMol>>()
            }),
        },
        "RDKitFPAtomInvGenerator" => AtomInvariantsGenerator::RDKitFP,
        "AtomPairAtomInvGenerator" => AtomInvariantsGenerator::AtomPair {
            include_chirality: pt.get_bool("includeChirality", false),
            topological_torsion_correction: pt.get_bool("topologicalTorsionCorrection", false),
        },
        other => return value_error(format!("Unknown AtomInvariantsGenerator type: {other}")),
    })
}

fn read_bond_inv(pt: &PTree) -> Result<BondInvariantsGenerator> {
    match typed(pt, "BondInvariantsGenerator")? {
        "MorganBondInvGenerator" => Ok(BondInvariantsGenerator::Morgan {
            use_bond_types: pt.get_bool("useBondTypes", true),
            use_chirality: pt.get_bool("useChirality", false),
        }),
        other => value_error(format!("Unknown BondInvariantsGenerator type: {other}")),
    }
}

/// `generatorFromJSON`: the generator produces 64-bit results.
pub fn generator_from_json(json: &str) -> Result<FingerprintGenerator> {
    let pt = PTree::from_json(json)?;
    let Some(args_node) = pt.child("fingerprintArguments") else {
        return value_error("fingerprintArguments are missing from the JSON");
    };
    let (kind, env, args) = read_arguments(args_node)?;
    let Some(env_node) = pt.child("atomEnvironmentGenerator") else {
        return value_error("atomEnvironmentGenerator is missing from the JSON");
    };
    check_env(env_node, &kind)?;
    let atom_inv_gen = pt
        .child("atomInvariantsGenerator")
        .map(read_atom_inv)
        .transpose()?;
    let bond_inv_gen = pt
        .child("bondInvariantsGenerator")
        .map(read_bond_inv)
        .transpose()?;
    Ok(FingerprintGenerator {
        env,
        args,
        atom_inv_gen,
        bond_inv_gen,
        output_64: true,
    })
}
