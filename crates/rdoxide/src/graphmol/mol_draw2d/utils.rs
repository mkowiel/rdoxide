//! Preparing molecules for drawing and drawing-option presets
//! (`MolDraw2DUtils.cpp` of RDKit 2026.09.1, without the contouring
//! functions, which need the drawer).
//!
//! Copyright (C) 2016-2021 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::helpers::{DrawColour, MolDrawOptions};
use crate::graphmol::atom::ChiralType;
use crate::graphmol::depictor::{Compute2DCoordParameters, compute_2d_coords};
use crate::graphmol::file_parsers::mol_file_stereochem::{
    StereoBondThresholds, add_wavy_bonds_for_stereo_any,
};
use crate::graphmol::kekulize::{DEFAULT_MAX_BACKTRACKS, kekulize_if_possible};
use crate::graphmol::mol_ops::add_hs::{AddHsParameters, add_hs};
use crate::graphmol::mol_transforms::get_bond_length;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::wedge_bonds::{BondWedgingParameters, wedge_mol_bonds};
use crate::rdgeneral::log::LogStateSetter;
use crate::{RDKitError, Result};

/// A candidate for an explicit H to show its stereo: chirality specified and
/// in at least two rings.
fn is_atom_cand_for_chiral_h(mol: &RWMol, idx: u32) -> bool {
    mol.ring_info().is_initialized()
        && mol.ring_info().num_atom_rings(idx).unwrap_or(0) > 1
        && matches!(
            mol.atom(idx).chiral_tag(),
            ChiralType::TetrahedralCCW | ChiralType::TetrahedralCW
        )
}

/// `prepareMolForDrawing`: kekulizes (keeping the aromatic flags), adds
/// explicit Hs on stereo centres in ring fusions, generates coordinates when
/// needed or asked for, and sets wedging and wavy bonds.
pub fn prepare_mol_for_drawing(
    mol: &mut RWMol,
    kekulize: bool,
    add_chiral_hs: bool,
    wedge_bonds: bool,
    force_coords: bool,
    wavy_bonds: bool,
) -> Result<()> {
    if kekulize {
        let _blocker = LogStateSetter::new();
        kekulize_if_possible(mol, false, true, DEFAULT_MAX_BACKTRACKS)?;
    }
    if add_chiral_hs {
        let chiral_ats: Vec<u32> = (0..mol.num_atoms())
            .filter(|&i| is_atom_cand_for_chiral_h(mol, i))
            .collect();
        if !chiral_ats.is_empty() {
            let add_coords = !force_coords && mol.num_conformers() > 0;
            let params = AddHsParameters {
                explicit_only: false,
                add_coords,
                add_residue_info: false,
                skip_queries: false,
            };
            add_hs(mol, &params, Some(&chiral_ats))?;
        }
    }
    if force_coords || mol.num_conformers() == 0 {
        // 2D coordinates in a standard orientation
        let params = Compute2DCoordParameters {
            canon_orient: true,
            ..Default::default()
        };
        compute_2d_coords(mol, &params)?;
    }
    if wedge_bonds {
        wedge_mol_bonds(mol, None, &BondWedgingParameters::default())?;
    }
    if wavy_bonds {
        add_wavy_bonds_for_stereo_any(mol, true, StereoBondThresholds::DBL_BOND_NO_STEREO)?;
    }
    Ok(())
}

/// `meanBondLength`: the mean length of the bonds in the conformer.
pub fn mean_bond_length(mol: &ROMol, conf_id: i32) -> Result<f64> {
    let n = mol.num_bonds();
    if n == 0 {
        return Ok(0.0);
    }
    let conf = mol.conformer(conf_id)?;
    let mut total = 0.0;
    for b in mol.bonds() {
        total += get_bond_length(conf, b.begin_atom_idx(), b.end_atom_idx())?;
    }
    Ok(total / f64::from(n))
}

/// `setACS1996Options`: sets the options for the ACS 1996 drawing guidelines,
/// which fix a bond length of 14.4 pixels.
pub fn set_acs1996_options(opts: &mut MolDrawOptions, mean_bond_len: f64) -> Result<()> {
    if mean_bond_len <= 0.0 {
        return Err(RDKitError::Value(
            "ACS1996Options requires mean bond length > 0.0.".into(),
        ));
    }
    opts.bond_line_width = 0.6;
    opts.scale_bond_width = false;
    // the guideline is for a bond length of 14.4px, and things are set up in
    // pixels per Angstrom
    opts.scaling_factor = 14.4 / mean_bond_len;
    // a fixed bond length means the drawing won't be scaled up in a drawer of
    // defined size
    opts.fixed_bond_length = 14.4 / mean_bond_len;
    // the offset for multiple bonds is 18% of the bond length
    opts.multiple_bond_offset = 0.18;
    opts.highlight_bond_width_multiplier = 32;
    opts.set_monochrome_mode(
        DrawColour::new(0.0, 0.0, 0.0),
        DrawColour::new(1.0, 1.0, 1.0),
    );
    opts.fixed_font_size = 10;
    opts.additional_atom_label_padding = 0.066;
    // The guidelines say Arial, which is not free. FreeSans is close but GPL,
    // so it can't be embedded: use it if it is in $RDBASE/Data/Fonts, else
    // fall back on Roboto (Apache 2.0). A font file the user set is kept.
    if opts.font_file.is_empty() {
        let free_sans = std::env::var("RDBASE")
            .ok()
            .map(|b| format!("{b}/Data/Fonts/FreeSans.ttf"))
            .filter(|p| std::path::Path::new(p).exists());
        opts.font_file = free_sans.unwrap_or_else(|| "BuiltinRobotoRegular".to_string());
    }
    Ok(())
}
