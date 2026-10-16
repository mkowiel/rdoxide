//! Removal of the minor fragments of a molecule.
//!
//! Ported from RDKit `Code/GraphMol/StructChecker/StripSmallFragments.cpp`
//! (2026.09.1).
//!
//! Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::Result;
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::BondDir;
use crate::graphmol::chirality::clear_single_bond_dir_flags;
use crate::graphmol::chirality::from_coords::detect_bond_stereochemistry;
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::{get_mol_frags_as_mols, sanitize_mol};
use crate::graphmol::mol_props::{get_exact_mol_wt, get_mol_formula};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::{RDProps, RDValue, common_properties as cp};

/// `AddMWMF`: sets the molecular formula and mass properties
/// ("MF_PRE"/"MW_PRE" or "MF_POST"/"MW_POST").
pub fn add_mw_mf(mol: &mut ROMol, pre: bool) -> Result<()> {
    let mass = get_exact_mol_wt(mol, false)?;
    let formula = get_mol_formula(mol, false, true)?;
    if !formula.is_empty() {
        mol.set_prop(if pre { "MF_PRE" } else { "MF_POST" }, formula, false)?;
    }
    mol.set_prop(if pre { "MW_PRE" } else { "MW_POST" }, mass, false)?;
    Ok(())
}

/// `StripSmallFragments`: keeps only the biggest fragment (the last one of
/// the biggest); true if there was more than one.
pub fn strip_small_fragments(mol: &mut ROMol, verbose: bool) -> Result<bool> {
    let (mut frags, _) = get_mol_frags_as_mols(mol, false, true)?;
    if frags.len() <= 1 {
        return Ok(false);
    }
    let mut max_frag_size = 0;
    let mut max_frag_idx = 0;
    for (i, f) in frags.iter().enumerate() {
        let size = f.num_atoms();
        if size >= max_frag_size {
            max_frag_size = size;
            max_frag_idx = i;
        }
    }
    if verbose {
        let name = match mol.get_prop(cp::_Name) {
            Some(RDValue::String(s)) => s.clone(),
            Some(v) => v.to_rd_string(),
            None => "<no name>".to_string(),
        };
        for (i, f) in frags.iter().enumerate() {
            if i != max_frag_idx {
                crate::rd_log_warning!(
                    "{name} removed fragment i={i} with {} atoms\n",
                    f.num_atoms()
                );
            }
        }
    }

    // the chirality has to be saved for checking later
    let mut check_chiral = false;
    let mut kept = frags.swap_remove(max_frag_idx);
    if let Some(chiral_flag) = mol
        .get_prop(cp::_MolFileChiralFlag)
        .and_then(RDValue::as_uint)
    {
        kept.set_prop(cp::_MolFileChiralFlag, chiral_flag, false)?;
        check_chiral = chiral_flag != 0;
    }
    *mol = kept;

    // see if the chirality of the mol file possibly came from this fragment
    if check_chiral {
        let mut ischiral = false;
        let mut copy = mol.clone();
        let mut analyse = || -> Result<bool> {
            sanitize_mol(&mut copy, SANITIZE_ALL).map_err(|(_, e)| e)?;
            clear_single_bond_dir_flags(&mut copy, false);
            detect_bond_stereochemistry(&mut copy, -1)?;
            assign_stereochemistry(&mut copy, true, true, true)?;
            Ok(copy
                .atoms()
                .iter()
                .any(|a| a.has_prop(cp::_ChiralityPossible)))
        };
        if let Ok(true) = analyse() {
            ischiral = true;
            check_chiral = false;
        }
        // are chiral tags set?
        if check_chiral {
            if mol.atoms().iter().any(|a| {
                matches!(
                    a.chiral_tag(),
                    ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW
                )
            }) {
                ischiral = true;
            }
            if mol
                .bonds()
                .iter()
                .any(|b| matches!(b.bond_dir(), BondDir::BeginDash | BondDir::BeginWedge))
            {
                ischiral = true;
            }
        }
        if !ischiral {
            mol.set_prop(cp::_MolFileChiralFlag, 0u32, false)?;
        }
    }
    Ok(true)
}
