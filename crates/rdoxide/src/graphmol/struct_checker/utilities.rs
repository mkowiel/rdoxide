//! Neighbour lists, atom coordinates and logging helpers of the structure
//! checker.
//!
//! Ported from RDKit `Code/GraphMol/StructChecker/Utilites.cpp` (2026.09.1).
//!
//! Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::Result;
use crate::geometry::Point3D;
use crate::graphmol::bond::BondType;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::{RDProps, RDValue, common_properties as cp};

/// A set of the neighbours of an atom.
#[derive(Debug, Clone, Default)]
pub struct Neighbourhood {
    /// Indices of the neighbour atoms.
    pub atoms: Vec<u32>,
    /// Indices of the bonds to them.
    pub bonds: Vec<u32>,
}

/// `SetupNeighbourhood`
pub fn setup_neighbourhood(mol: &ROMol) -> Vec<Neighbourhood> {
    let mut neighbours = vec![Neighbourhood::default(); mol.num_atoms() as usize];
    for i in 0..mol.num_bonds() {
        let bond = mol.bond(i);
        let a1 = bond.begin_atom_idx() as usize;
        let a2 = bond.end_atom_idx() as usize;
        neighbours[a1].atoms.push(a2 as u32);
        neighbours[a1].bonds.push(i);
        neighbours[a2].atoms.push(a1 as u32);
        neighbours[a2].bonds.push(i);
    }
    neighbours
}

/// `getMolAtomPoints`: the coordinates of the atoms (taken from the first 3D
/// conformer, or from the first conformer when `twod`) and whether any of
/// them has a non-zero z.
pub fn get_mol_atom_points(mol: &ROMol, twod: bool) -> Result<(Vec<Point3D>, bool)> {
    let mut non_zero_z = false;
    let mut points = vec![Point3D::default(); mol.num_atoms() as usize];
    for conf in mol.conformers() {
        if twod || conf.is_3d() {
            for (i, p) in points.iter_mut().enumerate() {
                *p = conf.atom_pos(i as u32)?;
                if p.z.abs() >= 1.0e-7 {
                    non_zero_z = true;
                }
            }
            break;
        }
    }
    Ok((points, non_zero_z))
}

/// `LogNeighbourhood`
pub fn log_neighbourhood(mol: &ROMol, idx: u32, neighbours: &[Neighbourhood]) -> String {
    let mut out = String::new();
    let name = match mol.get_prop(cp::_Name) {
        Some(RDValue::String(s)) => s.clone(),
        Some(v) => v.to_rd_string(),
        None => String::new(),
    };
    let atm = mol.atom(idx);
    out.push_str(&format!("atom '{name}' idx={idx} AA: {}", atm.symbol()));
    if atm.formal_charge() != 0 {
        if atm.formal_charge() > 0 {
            out.push('+');
        }
        out.push_str(&atm.formal_charge().to_string());
    }
    if atm.num_radical_electrons() != 0 {
        out.push_str(&atm.num_radical_electrons().to_string());
    }

    let nbr = &neighbours[idx as usize];
    let mut nbrs: Vec<(String, i32, i32, i32)> = nbr
        .atoms
        .iter()
        .zip(&nbr.bonds)
        .map(|(&a, &b)| {
            let na = mol.atom(a);
            (
                na.symbol(),
                mol.bond(b).bond_type() as i32,
                na.formal_charge(),
                na.num_radical_electrons() as i32,
            )
        })
        .collect();
    nbrs.sort();
    for (symbol, bt, charge, radicals) in nbrs {
        let bs = match BondType::from_u8(bt as u8) {
            Some(BondType::Single) => "-",
            Some(BondType::Double) => "=",
            Some(BondType::Triple) => "#",
            Some(BondType::Aromatic) => "~",
            _ => "",
        };
        if bs.is_empty() {
            out.push_str(&format!("(?{bt}?{symbol}"));
        } else {
            out.push_str(&format!("({bs}{symbol}"));
        }
        if charge != 0 {
            if charge > 0 {
                out.push('+');
            }
            out.push_str(&charge.to_string());
        }
        if radicals != 0 {
            out.push_str(&radicals.to_string());
        }
        out.push(')');
    }
    out
}
