//! Port of RDKit `Code/GraphMol/MolEnumerator`: enumeration of the molecules
//! described by position variation bonds, link nodes and repeat units
//! (SRUs). Random enumeration is not supported (it is not upstream either).

pub mod link_node;
mod link_node_op;
mod position_variation;
mod repeat_unit;

use std::sync::Arc;

use crate::graphmol::mol_bundle::MolBundle;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::{RDProps, common_properties as cp};
use crate::{RDKitError, Result};

pub use link_node_op::LinkNodeOp;
pub use position_variation::PositionVariationOp;
pub use repeat_unit::RepeatUnitOp;

/// `MolEnumerator::detail::idxPropName`
pub const IDX_PROP_NAME: &str = "_enumeratorOrigIdx";

/// `detail::preserveOrigIndices`
pub(crate) fn preserve_orig_indices(mol: &mut ROMol) -> Result<()> {
    mol.set_prop(IDX_PROP_NAME, 1i32, false)?;
    for idx in 0..mol.num_atoms() {
        mol.atom_mut(idx).set_prop(IDX_PROP_NAME, idx, false)?;
    }
    Ok(())
}

/// `detail::removeOrigIndices`
pub(crate) fn remove_orig_indices(mol: &mut ROMol) {
    for idx in 0..mol.num_atoms() {
        mol.atom_mut(idx).clear_prop(IDX_PROP_NAME);
    }
    mol.clear_prop(IDX_PROP_NAME);
}

/// `MolEnumeratorOp`: one kind of enumeration.
pub trait MolEnumeratorOp: Send + Sync {
    /// The number of variations at each variation point
    fn variation_counts(&self) -> Vec<usize>;
    /// The molecule with the given variation at each point
    fn apply(&self, which: &[usize]) -> Result<ROMol>;
    fn init_from_mol(&mut self, mol: &ROMol) -> Result<()>;
    /// `copy`: a fresh operation of the same kind and settings
    fn copy(&self) -> Box<dyn MolEnumeratorOp>;
}

/// `MolEnumeratorParams`
#[derive(Clone)]
pub struct MolEnumeratorParams {
    pub sanitize: bool,
    pub max_to_enumerate: usize,
    /// not implemented
    pub do_random: bool,
    /// not implemented
    pub random_seed: i32,
    pub operation: Option<Arc<dyn MolEnumeratorOp>>,
}

impl std::fmt::Debug for MolEnumeratorParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MolEnumeratorParams")
            .field("sanitize", &self.sanitize)
            .field("max_to_enumerate", &self.max_to_enumerate)
            .field("do_random", &self.do_random)
            .field("random_seed", &self.random_seed)
            .field("operation", &self.operation.is_some())
            .finish()
    }
}

impl Default for MolEnumeratorParams {
    fn default() -> Self {
        MolEnumeratorParams {
            sanitize: false,
            max_to_enumerate: 1000,
            do_random: false,
            random_seed: -1,
            operation: None,
        }
    }
}

impl MolEnumeratorParams {
    pub fn with_operation(operation: Arc<dyn MolEnumeratorOp>) -> Self {
        MolEnumeratorParams {
            operation: Some(operation),
            ..Default::default()
        }
    }
}

fn get_variations(
    level: usize,
    mut base: Vec<usize>,
    variations: &mut Vec<Vec<usize>>,
    variation_counts: &[usize],
    max_to_enumerate: usize,
) {
    for i in 0..variation_counts[level] {
        base[level] = i;
        if level + 1 >= variation_counts.len() {
            variations.push(base.clone());
        } else {
            get_variations(
                level + 1,
                base.clone(),
                variations,
                variation_counts,
                max_to_enumerate,
            );
        }
        if variations.len() >= max_to_enumerate {
            return;
        }
    }
}

fn enumerate_variations(
    variation_counts: &[usize],
    params: &MolEnumeratorParams,
) -> Result<Vec<Vec<usize>>> {
    if params.do_random {
        return Err(RDKitError::Runtime(
            "random enumeration not yet supported".into(),
        ));
    }
    let mut variations = Vec::new();
    let base = vec![0; variation_counts.len()];
    get_variations(
        0,
        base,
        &mut variations,
        variation_counts,
        params.max_to_enumerate,
    );
    Ok(variations)
}

fn clear_reaction_props(mol: &mut ROMol) {
    mol.clear_computed_props(false);
    for idx in 0..mol.num_atoms() {
        let atom = mol.atom_mut(idx);
        atom.clear_prop(cp::reactantAtomIdx);
        atom.clear_prop(cp::reactantIdx);
        atom.clear_prop("was_dummy");
        atom.clear_prop(cp::reactionMapNum);
    }
}

/// `enumerate(mol, paramsList)`: the operations are applied one after the
/// other to all the molecules the previous ones made. An empty bundle is
/// returned when no operation found variations.
pub fn enumerate(mol: &ROMol, params_list: &[MolEnumeratorParams]) -> Result<MolBundle> {
    let mut mol_copy = mol.clone();
    preserve_orig_indices(&mut mol_copy)?;
    let mut accum: Vec<ROMol> = vec![mol_copy];
    let mut variations_found = false;
    for params in params_list {
        let Some(operation) = &params.operation else {
            return Err(RDKitError::Value(
                "MolEnumeratorParams has no operation set".into(),
            ));
        };
        let mut this_round = Vec::new();
        for tmol in &accum {
            let mut op = operation.copy();
            op.init_from_mol(tmol)?;
            let variation_counts = op.variation_counts();
            if variation_counts.is_empty() {
                let mut mcp = tmol.clone();
                clear_reaction_props(&mut mcp);
                this_round.push(mcp);
                continue;
            }
            for variation in enumerate_variations(&variation_counts, params)? {
                variations_found = true;
                let mut new_mol = op.apply(&variation)?;
                new_mol.update_property_cache(false)?;
                clear_reaction_props(&mut new_mol);
                this_round.push(new_mol);
            }
        }
        accum = this_round;
    }
    if !variations_found {
        return Ok(MolBundle::new());
    }
    let mut res = MolBundle::new();
    for mut m in accum {
        remove_orig_indices(&mut m);
        res.add_mol(Arc::new(m))?;
    }
    Ok(res)
}

/// `enumerate(mol, maxPerOperation)`: position variation, repeat units and
/// link nodes, in that order (`max_per_operation == 0`: the default limits).
pub fn enumerate_default(mol: &ROMol, max_per_operation: usize) -> Result<MolBundle> {
    let limit = |p: &mut MolEnumeratorParams| {
        if max_per_operation > 0 {
            p.max_to_enumerate = max_per_operation;
        }
    };
    let mut position = MolEnumeratorParams::with_operation(Arc::new(PositionVariationOp::new()));
    limit(&mut position);
    let mut repeat_op = RepeatUnitOp::new();
    if max_per_operation > 0 {
        repeat_op.max_num_rounds = max_per_operation;
    }
    let mut repeat = MolEnumeratorParams::with_operation(Arc::new(repeat_op));
    limit(&mut repeat);
    let mut link = MolEnumeratorParams::with_operation(Arc::new(LinkNodeOp::new()));
    limit(&mut link);
    enumerate(mol, &[position, repeat, link])
}
