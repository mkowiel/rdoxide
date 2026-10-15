//! Building a Z-matrix of a molecule from internal coordinates (upstream
//! `ZMatrixBuilder.cpp`).

use std::f64::consts::PI;

use super::details::{BitSet, unified_id2, unified_id3};
use super::internal_coordinates::InternalCoordinates;
use crate::distgeom::{InternalCoordinateReference, TorsionCandidates, TorsionDependence, ZMatrix};
use crate::error::{RDKitError, Result};
use crate::graphmol::atom::ChiralType;
use crate::graphmol::romol::ROMol;
use crate::rd_log_warning;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;

/// Per atom: the atom that depends on it for torsions, the bond (angle)
/// reference and the torsion reference.
type Type14Reference = (Option<usize>, usize, Option<usize>);

#[derive(Debug, Clone, Copy)]
struct StackElem {
    atom_idx: u32,
    precursor_idx: u32,
}

fn bond_between(mol: &ROMol, a: u32, b: u32) -> Result<u32> {
    mol.bond_idx_between(a, b)
        .ok_or_else(|| RDKitError::invariant("no bond between the atoms"))
}

#[allow(clippy::too_many_arguments)]
fn add_neighbors_to_stack(
    atom_idx: u32,
    precursor_idx: u32,
    mol: &ROMol,
    stack: &mut Vec<StackElem>,
    visited_atoms: &BitSet,
    zmat: &ZMatrix,
    internal_coords: &InternalCoordinates,
) -> Result<()> {
    let mut bonds: Vec<u32> = Vec::new();
    let mut ring_closures: Vec<u32> = Vec::new();

    for bid in mol.atom_bond_indices(atom_idx) {
        let nbr = mol.bond(bid).other_atom_idx(atom_idx)?;
        if nbr == precursor_idx {
            continue;
        }
        if visited_atoms.get(nbr as usize) {
            ring_closures.push(nbr);
        } else {
            bonds.push(bid);
        }
    }
    if bonds.is_empty() && ring_closures.is_empty() {
        return Ok(());
    }

    let references = zmat.references(atom_idx as usize)?;
    if let (Some(angle_ref), false) = (references.angle_ref, bonds.is_empty()) {
        // if torsions are given, the most constrained one goes last onto the
        // stack
        let bond_idx1 = bond_between(mol, angle_ref as u32, precursor_idx)?;
        let bond_idx2 = bond_between(mol, precursor_idx, atom_idx)?;
        let nb = mol.num_bonds();
        let torsion_of = |bond: u32| -> Result<&TorsionCandidates> {
            internal_coords
                .torsion_range
                .get(&unified_id3(bond_idx1, bond_idx2, bond, nb))
                .ok_or_else(|| RDKitError::Key("torsion not found".into()))
        };
        let mut first = 0;
        for i in 1..bonds.len() {
            if torsion_of(bonds[i])?.less(torsion_of(bonds[first])?) {
                first = i;
            }
        }
        let torsion = torsion_of(bonds[first])?.clone();
        // it is handled first, so it is added last to the stack
        let next = (first + 1) % bonds.len();
        bonds.rotate_left(next);

        if let TorsionCandidates::Range(range) = &torsion {
            if range.upper - range.lower >= PI * 2.0 - 1e-6 {
                // completely unconstrained: non-ring bonds are handled first
                let rinfo = mol.ring_info();
                let sizes: Vec<(u32, u32)> = bonds
                    .iter()
                    .map(|&b| Ok((b, rinfo.min_bond_ring_size(b)?)))
                    .collect::<Result<_>>()?;
                let mut sizes = sizes;
                libcxx_sort_by(&mut sizes, |a, b| a.1 > b.1);
                bonds = sizes.into_iter().map(|(b, _)| b).collect();
            }
        }
    }

    for &bid in &bonds {
        stack.push(StackElem {
            atom_idx: mol.bond(bid).other_atom_idx(atom_idx)?,
            precursor_idx: atom_idx,
        });
    }
    // ring closures are handled first
    for nbr in ring_closures {
        stack.push(StackElem {
            atom_idx: nbr,
            precursor_idx: atom_idx,
        });
    }
    Ok(())
}

fn ring_closure(atom_idx: usize, precursor_idx: usize, references: &mut [Option<Type14Reference>]) {
    // the precursor was set already (the first atom cannot close a ring)
    let (_, ref1, ref2) = references[precursor_idx].expect("the precursor is set");
    if ref2 == Some(atom_idx) {
        // small rings (3 atoms) or highly fused systems
        return;
    }
    references[precursor_idx] = Some((Some(atom_idx), ref1, ref2));
    let (_, angle_ref, torsion_ref) = references[atom_idx].expect("references are set first");
    // only matters if the ring closure is more or less chemically valid
    references[atom_idx] = Some((Some(precursor_idx), angle_ref, torsion_ref));
}

fn add_element(
    atom_idx: u32,
    precursor_idx: u32,
    mol: &ROMol,
    zmat: &mut ZMatrix,
    internal_coords: &InternalCoordinates,
    references: &mut [Option<Type14Reference>],
) -> Result<()> {
    let (t_ref, ref1, ref2) = references[precursor_idx as usize].expect("the precursor is set");
    let nb = mol.num_bonds();
    let bnd_idx1 = bond_between(mol, atom_idx, precursor_idx)?;
    let bnd_idx2 = bond_between(mol, ref1 as u32, precursor_idx)?;
    let bond_length = internal_coords.lengths[bnd_idx1 as usize];

    let angle_key = unified_id2(bnd_idx1, bnd_idx2, nb);
    let Some(&angle) = internal_coords.angles.get(&angle_key) else {
        rd_log_warning!("Angle {bnd_idx1} {bnd_idx2} not found\n");
        return Err(RDKitError::Key("angle not found".into()));
    };
    let mut angle = angle;
    if (angle - PI).abs() <= 1e-8 {
        // to avoid colinearity
        angle += 1e-6;
    }

    let mut torsion = None;
    let mut torsion_dependence = None;
    if let Some(t_ref) = t_ref {
        // this also means that a torsion reference is not needed
        let degree = mol.atom_degree(precursor_idx).min(5);
        let offset = 2.0 * PI / f64::from(degree.wrapping_sub(1));
        torsion_dependence = Some(TorsionDependence {
            offset,
            reference: t_ref,
        });
    } else if let Some(ref2) = ref2 {
        // the 3rd element of the matrix has no torsion reference
        let bnd3 = bond_between(mol, ref1 as u32, ref2 as u32)?;
        let key = unified_id3(bnd_idx1, bnd_idx2, bnd3, nb);
        let Some(range) = internal_coords.torsion_range.get(&key) else {
            rd_log_warning!("Torsion {bnd_idx1} {bnd_idx2} {bnd3} not found\n");
            return Err(RDKitError::Key("torsion not found".into()));
        };
        torsion = Some(range.clone());
    }

    zmat.add_element(
        atom_idx as usize,
        InternalCoordinateReference {
            length: Some(bond_length),
            bond_ref: Some(precursor_idx as usize),
            angle: Some(angle),
            angle_ref: Some(ref1),
            torsion,
            torsion_ref: if torsion_dependence.is_some() {
                None
            } else {
                ref2
            },
        },
        torsion_dependence,
    )?;

    references[precursor_idx as usize].as_mut().expect("set").0 = Some(atom_idx as usize);
    debug_assert!(references[atom_idx as usize].is_none());
    references[atom_idx as usize] = Some((None, precursor_idx as usize, Some(ref1)));
    Ok(())
}

/// Fills the Z-matrix by a depth-first walk through the molecule; it starts
/// with a non-ring bond at the end of a chain, if there is one.
pub fn set_molecule_dfs(
    mol: &ROMol,
    zmat: &mut ZMatrix,
    internal_coords: &InternalCoordinates,
) -> Result<()> {
    if mol.num_atoms() == 1 {
        return set_molecule_dfs_from(mol, zmat, internal_coords, 0, 0);
    }
    // prefer bonds where neither atom is in a ring, then low degrees
    let rinfo = mol.ring_info();
    let key = |bond: &crate::graphmol::bond::Bond| -> Result<(u32, u32, u32)> {
        let rings = rinfo
            .num_atom_rings(bond.begin_atom_idx())?
            .max(rinfo.num_atom_rings(bond.end_atom_idx())?);
        let degree = mol
            .atom_degree(bond.begin_atom_idx())
            .min(mol.atom_degree(bond.end_atom_idx()));
        Ok((rings, degree, bond.idx()))
    };
    let mut start: Option<((u32, u32, u32), u32)> = None;
    for bond in mol.bonds() {
        let k = key(bond)?;
        if start.is_none_or(|(best, _)| k < best) {
            start = Some((k, bond.idx()));
        }
    }
    let bond = mol.bond(
        start
            .ok_or_else(|| RDKitError::precondition("molecule has no bonds"))?
            .1,
    );
    let (mut first, mut second) = (bond.begin_atom_idx(), bond.end_atom_idx());
    if mol.atom_degree(first) > mol.atom_degree(second) {
        std::mem::swap(&mut first, &mut second);
    }
    set_molecule_dfs_from(mol, zmat, internal_coords, first, second)
}

/// Like [`set_molecule_dfs`] with the first two rows given.
pub fn set_molecule_dfs_from(
    mol: &ROMol,
    zmat: &mut ZMatrix,
    internal_coords: &InternalCoordinates,
    mut first_atom_idx: u32,
    mut second_atom_idx: u32,
) -> Result<()> {
    let n_atoms = mol.num_atoms() as usize;
    let mut references: Vec<Option<Type14Reference>> = vec![None; n_atoms];
    let mut visited_atoms = BitSet::new(n_atoms);
    let mut visited_bonds = BitSet::new(mol.num_bonds() as usize);
    let mut stack: Vec<StackElem> = Vec::with_capacity(n_atoms / 2);

    if n_atoms > 2 && mol.atom_degree(second_atom_idx) <= 1 {
        // the second atom should have another neighbor than the first one
        // (for a bond angle reference)
        std::mem::swap(&mut first_atom_idx, &mut second_atom_idx);
    }

    zmat.add_element(
        first_atom_idx as usize,
        InternalCoordinateReference::default(),
        None,
    )?;
    if n_atoms == 1 {
        return Ok(());
    }
    let first_bond = bond_between(mol, first_atom_idx, second_atom_idx)?;
    zmat.add_element(
        second_atom_idx as usize,
        InternalCoordinateReference {
            length: Some(internal_coords.lengths[first_bond as usize]),
            bond_ref: Some(first_atom_idx as usize),
            ..Default::default()
        },
        None,
    )?;
    if n_atoms == 2 {
        return Ok(());
    }

    add_neighbors_to_stack(
        first_atom_idx,
        second_atom_idx,
        mol,
        &mut stack,
        &visited_atoms,
        zmat,
        internal_coords,
    )?;
    visited_atoms.set(first_atom_idx as usize);
    add_neighbors_to_stack(
        second_atom_idx,
        first_atom_idx,
        mol,
        &mut stack,
        &visited_atoms,
        zmat,
        internal_coords,
    )?;
    visited_atoms.set(second_atom_idx as usize);
    visited_bonds.set(first_bond as usize);

    // the stack holds at least one element (more than 2 atoms) and it is not a
    // ring closure since the matrix has only two rows
    let last = stack
        .last()
        .expect("a neighbor of the first atoms")
        .atom_idx;
    references[first_atom_idx as usize] =
        Some((None, second_atom_idx as usize, Some(last as usize)));
    references[second_atom_idx as usize] = Some((None, first_atom_idx as usize, None));

    while let Some(StackElem {
        atom_idx,
        precursor_idx,
    }) = stack.pop()
    {
        let bnd_idx = bond_between(mol, atom_idx, precursor_idx)?;
        if visited_bonds.get(bnd_idx as usize) {
            continue;
        }
        if !visited_atoms.get(atom_idx as usize) {
            add_element(
                atom_idx,
                precursor_idx,
                mol,
                zmat,
                internal_coords,
                &mut references,
            )?;
            visited_atoms.set(atom_idx as usize);
            add_neighbors_to_stack(
                atom_idx,
                precursor_idx,
                mol,
                &mut stack,
                &visited_atoms,
                zmat,
                internal_coords,
            )?;
        } else {
            ring_closure(atom_idx as usize, precursor_idx as usize, &mut references);
        }
        visited_bonds.set(bnd_idx as usize);
    }
    Ok(())
}

/// Corrects the sign of the offsets of improper torsions at tetrahedral
/// stereocenters so that the chirality is reproduced.
pub fn correct_chiral_centers(mol: &ROMol, zmat: &mut ZMatrix) -> Result<()> {
    // the first four atoms can only span a torsion, which cannot depend on
    // another one
    let rows: Vec<(usize, usize, usize, usize)> = zmat
        .rows()
        .skip(3)
        .filter_map(|row| {
            let dep = row.torsion_dependence?;
            Some((
                row.atom_idx,
                row.internal.bond_ref?,
                row.internal.angle_ref?,
                dep.reference,
            ))
        })
        .collect();
    for (current, center, axis, anchor) in rows.into_iter().rev() {
        let chiral_tag = mol.atom(center as u32).chiral_tag();
        if !matches!(
            chiral_tag,
            ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW
        ) {
            continue;
        }
        let b = |other: usize| bond_between(mol, center as u32, other as u32);
        let (bnd1, bnd2, bnd3) = (b(axis)?, b(anchor)?, b(current)?);
        let mut order = vec![bnd1 as i32, bnd2 as i32, bnd3 as i32];
        // the remaining bond
        for bid in mol.atom_bond_indices(center as u32) {
            if bid != bnd1 && bid != bnd2 && bid != bnd3 {
                order.push(bid as i32);
            }
        }
        let is_ccw = mol.perturbation_order(center as u32, &order)? % 2 != 0;
        if is_ccw != (chiral_tag == ChiralType::TetrahedralCCW) {
            // wrong order: invert the offsets
            zmat.invert_torsion_dependence(current)?;
        }
    }
    Ok(())
}
