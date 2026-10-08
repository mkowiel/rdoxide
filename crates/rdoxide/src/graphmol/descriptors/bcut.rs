//! BCUT2D descriptors (upstream `GraphMol/Descriptors/BCUT.cpp`), from
//! J. Chem. Inf. Comput. Sci., Vol. 39, No. 1, 1999.
//!
//! The Burden matrix holds an atom property on the diagonal and
//! `1/sqrt(bond order)` for bonded atom pairs (0.001 for all other pairs); the
//! descriptor is its highest and lowest eigenvalue.
//!
//! Copyright (C) 2020 Brian P. Kelley.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::graphmol::bond::BondType;
use crate::graphmol::descriptors::crippen::get_crippen_atom_contribs;
use crate::graphmol::mol_ops::add_hs::remove_all_hs;
use crate::graphmol::partial_charges::gasteiger::compute_gasteiger_charges;
use crate::graphmol::romol::ROMol;
use crate::numerics::symmetric_eigenvalues;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::dict::{RDProps, RDValue};
use crate::{RDKitError, Result};

pub const BCUT2D_VERSION: &str = "1.0.0";
const NUM_BCUTS: usize = 8;

/// Off-diagonal value of atoms that are not bonded.
const UNBONDED: f64 = 0.001;

/// `make_burden`: the row-major matrix with the diagonal still unset (0).
fn make_burden(mol: &ROMol) -> Result<Vec<f64>> {
    let n = mol.num_atoms() as usize;
    let mut burden = vec![UNBONDED; n * n];
    for bond in mol.bonds() {
        let score = match bond.bond_type() {
            BondType::Aromatic => 0.816_496_580_927_726_1, // 1/sqrt(1.5)
            BondType::Single => 1.0,
            BondType::Double => 0.707_106_781_186_547_5, // 1/sqrt(2)
            BondType::Triple => 0.577_350_269_189_625_8, // 1/sqrt(3)
            _ => {
                return Err(RDKitError::invariant(
                    "Bond order must be Single, Double, Triple or Aromatic",
                ));
            }
        };
        let (i, j) = (bond.begin_atom_idx() as usize, bond.end_atom_idx() as usize);
        burden[i * n + j] = score;
        burden[j * n + i] = score;
    }
    Ok(burden)
}

/// (highest, lowest) eigenvalue of the Burden matrix with `atom_props` on the
/// diagonal.
fn bcut_from_burden(burden: &mut [f64], atom_props: &[f64]) -> Result<(f64, f64)> {
    let n = atom_props.len();
    if n == 0 {
        return Ok((0.0, 0.0));
    }
    for (i, &p) in atom_props.iter().enumerate() {
        burden[i * n + i] = p;
    }
    let eigenvalues = symmetric_eigenvalues(n, burden)?;
    Ok((eigenvalues[n - 1], eigenvalues[0]))
}

/// `BCUT2D(mol, atom_props)`: (highest, lowest) eigenvalue.
pub fn bcut2d_with_props(mol: &ROMol, atom_props: &[f64]) -> Result<(f64, f64)> {
    if atom_props.len() != mol.num_atoms() as usize {
        return Err(RDKitError::precondition(
            "Number of atom props not equal to number of atoms",
        ));
    }
    if mol.num_atoms() == 0 {
        return Ok((0.0, 0.0));
    }
    let mut burden = make_burden(mol)?;
    bcut_from_burden(&mut burden, atom_props)
}

/// `BCUT2D(mol, atom_double_prop)`: uses the double property of every atom.
pub fn bcut2d_with_prop_name(mol: &ROMol, prop: &str) -> Result<(f64, f64)> {
    let props = mol
        .atoms()
        .iter()
        .map(|atom| match atom.get_prop(prop) {
            Some(RDValue::Double(d)) => Ok(*d),
            Some(_) => Err(RDKitError::Runtime("bad any cast".into())),
            None => Err(RDKitError::Key(prop.to_string())),
        })
        .collect::<Result<Vec<f64>>>()?;
    bcut2d_with_props(mol, &props)
}

/// `BCUT2D(mol)`: the eight descriptors
/// `[mass hi, mass lo, gasteiger hi, gasteiger lo, logp hi, logp lo, mr hi, mr lo]`
/// computed on the molecule without hydrogens.
pub fn bcut2d(mol: &ROMol) -> Result<Vec<f64>> {
    if mol.num_atoms() == 0 {
        return Ok(vec![0.0; NUM_BCUTS]);
    }
    let mut mol = ROMol::copy_from(mol, false, -1);
    remove_all_hs(&mut mol, true)?;
    compute_gasteiger_charges(&mut mol, 12, true)?;
    let masses: Vec<f64> = mol.atoms().iter().map(|a| a.mass()).collect();
    let charges = mol
        .atoms()
        .iter()
        .map(|a| match a.get_prop(cp::_GasteigerCharge) {
            Some(RDValue::Double(d)) => Ok(*d),
            _ => Err(RDKitError::Key(cp::_GasteigerCharge.to_string())),
        })
        .collect::<Result<Vec<f64>>>()?;
    let contribs = get_crippen_atom_contribs(&mut mol, true)?;
    let mut burden = make_burden(&mol)?;
    let mut res = Vec::with_capacity(NUM_BCUTS);
    for props in [&masses, &charges, &contribs.logp, &contribs.mr] {
        let (hi, lo) = bcut_from_burden(&mut burden, props)?;
        res.push(hi);
        res.push(lo);
    }
    Ok(res)
}
