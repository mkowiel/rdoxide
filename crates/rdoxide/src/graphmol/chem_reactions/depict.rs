//! Reaction depiction.
//!
//! Ported from RDKit `Code/GraphMol/ChemReactions/ReactionDepict.cpp`
//! (BSD-3-Clause, Copyright (c) 2018 Novartis Institutes for BioMedical
//! Research Inc.): the templates are depicted one after the other, side by
//! side along the x axis.

use super::ChemicalReaction;
use crate::Result;
use crate::graphmol::conjug_hybrid::{set_conjugation, set_hybridization};
use crate::graphmol::depictor::{Compute2DCoordParameters, compute_2d_coords};
use crate::graphmol::romol::ROMol;

/// Parameters of [`compute_2d_coords_for_reaction`].
#[derive(Debug, Clone)]
pub struct ReactionDepictParams {
    /// gap between neighbouring templates
    pub spacing: f64,
    pub update_props: bool,
    pub canon_orient: bool,
    pub n_flips_per_sample: u32,
    pub n_samples: u32,
    pub sample_seed: i32,
    pub permute_deg4_nodes: bool,
}

/// Depicts one template and moves it to the right of `x_offset`; returns the
/// new offset. `initial_extent` is the start value of the running x range
/// upstream uses (it differs for reactants and products).
fn place(
    mol: &mut ROMol,
    p: &ReactionDepictParams,
    x_offset: f64,
    extent: (f64, f64),
) -> Result<f64> {
    if p.update_props {
        mol.update_property_cache(false)?;
        set_conjugation(mol);
        set_hybridization(mol);
    }
    compute_2d_coords(
        mol,
        &Compute2DCoordParameters {
            canon_orient: p.canon_orient,
            clear_confs: true,
            n_flips_per_sample: p.n_flips_per_sample,
            n_samples: p.n_samples,
            sample_seed: p.sample_seed,
            permute_deg4_nodes: p.permute_deg4_nodes,
            ..Default::default()
        },
    )?;
    let (mut min_x, mut max_x) = extent;
    let conf = mol.conformer_mut(-1)?;
    for pt in conf.positions() {
        min_x = min_x.min(pt.x);
        max_x = max_x.max(pt.x);
    }
    let x_offset = x_offset + (max_x - min_x) / 2.0;
    for pt in conf.positions_mut() {
        pt.x += x_offset;
        max_x = max_x.max(pt.x);
    }
    Ok(max_x + p.spacing)
}

/// `RDDepict::compute2DCoordsForReaction`
pub fn compute_2d_coords_for_reaction(
    rxn: &mut ChemicalReaction,
    p: &ReactionDepictParams,
) -> Result<()> {
    let mut x_offset = 0.0;
    for reactant in rxn.reactants_mut() {
        x_offset = place(reactant, p, x_offset, (1e8, -1e8))?;
    }
    for product in rxn.products_mut() {
        x_offset = place(product, p, x_offset, (100.0, -100.0))?;
    }
    Ok(())
}
