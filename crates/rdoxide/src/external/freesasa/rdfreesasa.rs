//! `External/FreeSASA/RDFreeSASA.cpp`: SASA of RDKit molecules (conformer coordinates, radii by atom classification from the
//! PDB residue info) with the results stored as atom and molecule properties.

use super::classifier::{AtomClass, Classifier};
use super::{Algorithm, Parameters, calc_coord};
use crate::graphmol::monomer_info::AtomMonomerType;
use crate::graphmol::queries::{AtomQuery, PropValue, make_prop_query};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::dict::RDProps;
use crate::{RDKitError, Result};

/// `FreeSASA::SASAOpts`
#[derive(Debug, Clone, Copy)]
pub struct SasaOpts {
    pub algorithm: Algorithm,
    pub classifier: Classifier,
    pub probe_radius: f64,
}

impl Default for SasaOpts {
    fn default() -> Self {
        SasaOpts {
            algorithm: Algorithm::LeeRichards,
            classifier: Classifier::Protor,
            probe_radius: 1.4,
        }
    }
}

/// `classifyAtoms`: the radii (0 for atoms without PDB residue info) and the atom properties `SASAClass` (int) and
/// `SASAClassName`; false when an atom with residue info could not be classified
pub fn classify_atoms(mol: &mut ROMol, opts: &SasaOpts) -> Result<(bool, Vec<f64>)> {
    let mut radii = Vec::new();
    let mut success = true;
    for idx in 0..mol.num_atoms() {
        let mut cls = AtomClass::Unknown;
        let mut classification = "Unclassified".to_string();
        let mut radius = 0.0;
        if let Some(info) = mol.atom(idx).monomer_info()
            && info.monomer_type == AtomMonomerType::PdbResidue
        {
            radius = opts.classifier.radius(&info.residue_name, &info.name);
            if radius == 0.0 {
                crate::rd_log_warning!("Atom {idx} has zero radius");
            }
            cls = opts.classifier.class(&info.residue_name, &info.name);
            if cls == AtomClass::Unknown {
                crate::rd_log_warning!("Atom {idx} could not be classified");
                success = false;
            } else {
                classification = cls.name().to_string();
            }
        }
        radii.push(radius);
        let atom = mol.atom_mut(idx);
        atom.set_prop("SASAClass", cls as i32, false)?;
        atom.set_prop("SASAClassName", classification, false)?;
    }
    Ok((success, radii))
}

/// `calcSASA`: the total area (or of the atoms matching `query`); the atom properties `SASA` and the molecule property `SASA`
/// are set as a side effect
pub fn calc_sasa(
    mol: &mut ROMol,
    radii: &[f64],
    conf_idx: i32,
    query: Option<&AtomQuery>,
    opts: &SasaOpts,
) -> Result<f64> {
    if mol.num_conformers() == 0 {
        return Err(RDKitError::precondition("No conformers in molecule"));
    }
    if mol.num_atoms() == 0 {
        return Err(RDKitError::precondition("Empty molecule"));
    }
    if radii.len() < mol.num_atoms() as usize {
        return Err(RDKitError::Value("Need one radius per atom".into()));
    }
    let mut xyz = Vec::with_capacity(3 * mol.num_atoms() as usize);
    for p in mol.conformer(conf_idx)?.positions() {
        xyz.extend_from_slice(&[p.x, p.y, p.z]);
    }
    let params = Parameters {
        alg: opts.algorithm,
        probe_radius: opts.probe_radius,
        ..Default::default()
    };
    let res = calc_coord(&xyz, &radii[..mol.num_atoms() as usize], &params);
    mol.set_prop("SASA", res.total, false)?;
    for (i, s) in res.sasa.iter().enumerate() {
        mol.atom_mut(i as u32).set_prop("SASA", *s, false)?;
    }
    let mut result = res.total;
    if let Some(q) = query {
        result = 0.0;
        for i in 0..mol.num_atoms() {
            if q.matches(mol, i)? {
                result += res.sasa[i as usize];
            }
        }
    }
    Ok(result)
}

/// `makeFreeSasaAPolarAtomQuery`
pub fn make_apolar_atom_query() -> AtomQuery {
    make_prop_query("SASAClassName", PropValue::String("Apolar".into()), 0.0)
}

/// `makeFreeSasaPolarAtomQuery`
pub fn make_polar_atom_query() -> AtomQuery {
    make_prop_query("SASAClassName", PropValue::String("Polar".into()), 0.0)
}
