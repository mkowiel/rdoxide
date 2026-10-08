//! Atom/bond query predicates (in progress).
//!
//! Ported from RDKit `Code/GraphMol/QueryOps.{h,cpp}` (2026.09.1).
//!
//! Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use crate::graphmol::romol::ROMol;

/// `queryIsAtomBridgehead`: at least three ring bonds, and every ring
/// containing the atom shares at least two bonds with another ring that also
/// contains the atom.
pub fn query_is_atom_bridgehead(mol: &ROMol, atom: u32) -> bool {
    if mol.atom_degree(atom) < 3 {
        return false;
    }
    let ri = mol.ring_info();
    if !ri.is_initialized() {
        return false;
    }
    let nb = mol.num_bonds() as usize;
    let mut atom_ring_bonds = vec![false; nb];
    let mut count = 0;
    for b in mol.atom_bond_indices(atom) {
        if ri.num_bond_rings(b).unwrap_or(0) > 0 {
            atom_ring_bonds[b as usize] = true;
            count += 1;
        }
    }
    if count < 3 {
        return false;
    }
    let bond_rings = ri.bond_rings();
    let mut bonds_in_ring_i = vec![false; nb];
    let mut rings_overlap = vec![false; bond_rings.len()];
    for i in 0..bond_rings.len() {
        bonds_in_ring_i.iter_mut().for_each(|x| *x = false);
        let mut atom_in_ring_i = false;
        for &bidx in &bond_rings[i] {
            bonds_in_ring_i[bidx as usize] = true;
            if atom_ring_bonds[bidx as usize] {
                atom_in_ring_i = true;
            }
        }
        if !atom_in_ring_i {
            continue;
        }
        for j in i + 1..bond_rings.len() {
            let mut overlap = 0;
            let mut atom_in_ring_j = false;
            for &bidx in &bond_rings[j] {
                if atom_ring_bonds[bidx as usize] {
                    atom_in_ring_j = true;
                }
                if bonds_in_ring_i[bidx as usize] {
                    overlap += 1;
                }
                if overlap >= 2 && atom_in_ring_j {
                    rings_overlap[i] = true;
                    rings_overlap[j] = true;
                    break;
                }
            }
        }
        if !rings_overlap[i] {
            return false;
        }
    }
    true
}

// ------------------------------------------------------------------ query atoms/bonds

use crate::Result;
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::Bond;
use crate::graphmol::queries::{
    AtomQuery, BondQuery, CompositeQueryType, expand_query, make_atom_formal_charge_query,
    make_atom_isotope_query, make_atom_mass_query, make_atom_num_query,
    make_atom_num_radical_electrons_query, make_bond_order_equals_query,
};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;

/// `QueryAtom::expandQuery(what, how, maintainOrder)` on an atom that
/// already has a query.
pub fn atom_expand_query(
    atom: &mut Atom,
    what: AtomQuery,
    how: CompositeQueryType,
    maintain_order: bool,
) -> Result<()> {
    let current = atom
        .query
        .take()
        .ok_or_else(|| crate::RDKitError::precondition("atom has no query"))?;
    let q = expand_query(*current, what, how, maintain_order, "Atom")?;
    atom.set_query(Some(q));
    Ok(())
}

/// `QueryBond::expandQuery(what, how, maintainOrder)`
pub fn bond_expand_query(
    bond: &mut Bond,
    what: BondQuery,
    how: CompositeQueryType,
    maintain_order: bool,
) -> Result<()> {
    let current = bond
        .query
        .take()
        .ok_or_else(|| crate::RDKitError::precondition("bond has no query"))?;
    let q = expand_query(*current, what, how, maintain_order, "Bond")?;
    bond.set_query(Some(q));
    Ok(())
}

/// `QueryAtom(const Atom &other)`: a query atom equivalent to `other`.
pub fn query_atom_from_atom(other: &Atom) -> Result<Atom> {
    let mut a = other.copy_detached();
    if other.has_query() {
        return Ok(a);
    }
    a.set_query(Some(make_atom_num_query(other.atomic_num() as i32)));
    let and = CompositeQueryType::And;
    if other.isotope() != 0 {
        atom_expand_query(
            &mut a,
            make_atom_isotope_query(other.isotope() as i32),
            and,
            true,
        )?;
    }
    if other.formal_charge() != 0 {
        atom_expand_query(
            &mut a,
            make_atom_formal_charge_query(other.formal_charge()),
            and,
            true,
        )?;
    }
    if other.num_radical_electrons() != 0 {
        atom_expand_query(
            &mut a,
            make_atom_num_radical_electrons_query(other.num_radical_electrons() as i32),
            and,
            true,
        )?;
    }
    Ok(a)
}

/// `QueryBond(const Bond &other)`
pub fn query_bond_from_bond(other: &Bond) -> Bond {
    let mut b = other.copy_detached();
    if !other.has_query() {
        b.set_query(Some(make_bond_order_equals_query(other.bond_type())));
    }
    b
}

/// `QueryOps::replaceAtomWithQueryAtom`
pub fn replace_atom_with_query_atom(mol: &mut ROMol, idx: u32) -> Result<()> {
    if mol.get_atom(idx)?.has_query() {
        return Ok(());
    }
    let atom = mol.atom(idx);
    let mut qa = query_atom_from_atom(atom)?;
    if atom.has_prop(cp::_hasMassQuery) {
        atom_expand_query(
            &mut qa,
            make_atom_mass_query(atom.mass() as i32),
            CompositeQueryType::And,
            true,
        )?;
    }
    mol.replace_atom(idx, &qa, false)
}
