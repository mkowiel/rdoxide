//! Port of RDKit `Code/GraphMol/GenericGroups/GenericGroups.{h,cpp}`: matchers
//! for generic (Reaxys-style) atom labels such as `ALK` or `G*`, used by
//! substructure matching with `useGenericMatchers`.
//!
//! Matchers take the index of the atom the generic atom is matched to and the
//! atoms to ignore (the ones already matched by the query); their `ignore`
//! argument is by value upstream, so nothing is shared between calls.
//! Upstream perceives rings on the (mutable) ring info of a const molecule;
//! here a copy with rings is made when the molecule has none.

use std::borrow::Cow;
use std::collections::{BTreeSet, VecDeque};

use super::adjust_query::{AdjustQueryParameters, adjust_query_properties_in_place};
use super::atom::Atom;
use super::bond::{Bond, BondType};
use super::find_rings::{fast_find_rings, find_sssr};
use super::queries::make_atom_null_query;
use super::romol::{ROMol, RWMol};
use super::substance_group::{SubstanceGroup, add_substance_group};
use crate::Result;
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;

/// The signature of all the matchers.
pub type GenericMatcher = fn(&ROMol, u32, &[bool]) -> bool;

type AtomFn<'a> = &'a dyn Fn(&Atom) -> bool;
type BondFn<'a> = &'a dyn Fn(&Bond) -> bool;

/// Rings of `mol` (a copy if it has none, or no SSSR where `sssr` says so).
fn with_rings(mol: &ROMol, sssr: bool) -> Cow<'_, ROMol> {
    let ok = if sssr {
        mol.ring_info().is_initialized()
    } else {
        mol.ring_info().is_find_fast_or_better()
    };
    if ok {
        return Cow::Borrowed(mol);
    }
    let mut local = mol.clone();
    let found = if sssr {
        find_sssr(&mut local, false, false).map(|_| ())
    } else {
        fast_find_rings(&mut local)
    };
    // a molecule whose rings cannot be perceived has no usable ring info
    found.map_or(Cow::Borrowed(mol), |()| Cow::Owned(local))
}

fn num_atom_rings(mol: &ROMol, idx: u32) -> u32 {
    mol.ring_info().num_atom_rings(idx).unwrap_or(0)
}

/// `queryIsBondInRing`
fn bond_in_ring(mol: &ROMol, bond: &Bond) -> bool {
    mol.ring_info().num_bond_rings(bond.idx()).unwrap_or(0) > 0
}

fn is_hydrogen(mol: &ROMol, idx: u32) -> bool {
    mol.atom(idx).atomic_num() == 1 && mol.atom_degree(idx) == 1
}

fn is_carbon_or_h(at: &Atom) -> bool {
    at.atomic_num() == 6 || at.atomic_num() == 1
}

fn is_aromatic_bond(b: &Bond) -> bool {
    b.is_aromatic() || b.bond_type() == BondType::Aromatic
}

/// `AllAtomsMatch`: breadth-first over the side chain starting at `start`.
fn all_atoms_match(
    mol: &ROMol,
    start: u32,
    mut ignore: Vec<bool>,
    matcher: Option<AtomFn>,
    bond_matcher: Option<BondFn>,
    at_least_one_atom: Option<AtomFn>,
    at_least_one_bond: Option<BondFn>,
) -> bool {
    if matcher.is_some_and(|m| !m(mol.atom(start))) {
        return false;
    }
    let mut atom_at_least = at_least_one_atom.is_none();
    let mut bond_at_least = at_least_one_bond.is_none();
    let mut queue = VecDeque::from([start]);
    while let Some(idx) = queue.pop_front() {
        if !atom_at_least && at_least_one_atom.is_some_and(|f| f(mol.atom(idx))) {
            atom_at_least = true;
        }
        ignore[idx as usize] = true;
        for nbr in mol.atom_neighbors(idx) {
            if ignore[nbr as usize] {
                continue;
            }
            if matcher.is_some_and(|m| !m(mol.atom(nbr))) {
                return false;
            }
            if bond_matcher.is_some() || !bond_at_least {
                let Some(bi) = mol.bond_idx_between(idx, nbr) else {
                    continue;
                };
                let bond = mol.bond(bi);
                if bond_matcher.is_some_and(|m| !m(bond)) {
                    return false;
                }
                if !bond_at_least && at_least_one_bond.is_some_and(|f| f(bond)) {
                    bond_at_least = true;
                }
            }
            queue.push_back(nbr);
        }
    }
    atom_at_least && bond_at_least
}

/// Wraps a matcher so that an H atom (and nothing else) also matches.
macro_rules! with_h {
    ($(#[$m:meta])* $name:ident, $base:ident) => {
        $(#[$m])*
        pub fn $name(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
            is_hydrogen(mol, atom) || $base(mol, atom, ignore)
        }
    };
}

/// Single, non-aromatic, non-ring bond (alkyl chains).
fn chain_single_bond(mol: &ROMol) -> impl Fn(&Bond) -> bool + '_ {
    move |b| b.bond_type() == BondType::Single && !b.is_aromatic() && !bond_in_ring(mol, b)
}

fn aliphatic_c_or_h(at: &Atom) -> bool {
    !at.is_aromatic() && is_carbon_or_h(at)
}

fn is_carbon(at: &Atom) -> bool {
    at.atomic_num() == 6
}

fn is_not_h(at: &Atom) -> bool {
    at.atomic_num() != 1
}

/// Reaxys query type G: any sidechain with at least one non-hydrogen atom.
pub fn group_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    all_atoms_match(
        mol,
        atom,
        ignore.to_vec(),
        None,
        None,
        Some(&is_not_h),
        None,
    )
}
with_h!(
    /// Reaxys query type GH
    group_h_atom_matcher,
    group_atom_matcher
);

/// Reaxys query type G*: a sidechain with a ring closure.
pub fn group_star_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let mol = with_rings(mol, false);
    let in_ring = |b: &Bond| bond_in_ring(&mol, b);
    all_atoms_match(
        &mol,
        atom,
        ignore.to_vec(),
        None,
        None,
        Some(&is_not_h),
        Some(&in_ring),
    )
}
with_h!(
    /// Reaxys query type GH*
    group_star_h_atom_matcher,
    group_star_atom_matcher
);

pub fn alkyl_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let mol = with_rings(mol, false);
    let bonds = chain_single_bond(&mol);
    all_atoms_match(
        &mol,
        atom,
        ignore.to_vec(),
        Some(&aliphatic_c_or_h),
        Some(&bonds),
        Some(&is_carbon),
        None,
    )
}

pub fn alkyl_h_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let mol = with_rings(mol, false);
    let bonds = chain_single_bond(&mol);
    all_atoms_match(
        &mol,
        atom,
        ignore.to_vec(),
        Some(&aliphatic_c_or_h),
        Some(&bonds),
        None,
        None,
    )
}

pub fn acyclic_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let mol = with_rings(mol, false);
    let acyclic = |at: &Atom| num_atom_rings(&mol, at.idx()) == 0;
    all_atoms_match(
        &mol,
        atom,
        ignore.to_vec(),
        Some(&acyclic),
        None,
        Some(&is_not_h),
        None,
    )
}

pub fn acyclic_h_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let mol = with_rings(mol, false);
    let acyclic = |at: &Atom| num_atom_rings(&mol, at.idx()) == 0;
    all_atoms_match(
        &mol,
        atom,
        ignore.to_vec(),
        Some(&acyclic),
        None,
        None,
        None,
    )
}

pub fn carboacyclic_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let mol = with_rings(mol, false);
    let m = |at: &Atom| is_carbon_or_h(at) && num_atom_rings(&mol, at.idx()) == 0;
    all_atoms_match(
        &mol,
        atom,
        ignore.to_vec(),
        Some(&m),
        None,
        Some(&is_carbon),
        None,
    )
}

pub fn carboacyclic_h_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let mol = with_rings(mol, false);
    let m = |at: &Atom| is_carbon_or_h(at) && num_atom_rings(&mol, at.idx()) == 0;
    all_atoms_match(&mol, atom, ignore.to_vec(), Some(&m), None, None, None)
}

pub fn heteroacyclic_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let mol = with_rings(mol, false);
    let acyclic = |at: &Atom| num_atom_rings(&mol, at.idx()) == 0;
    let hetero = |at: &Atom| !is_carbon_or_h(at);
    all_atoms_match(
        &mol,
        atom,
        ignore.to_vec(),
        Some(&acyclic),
        None,
        Some(&hetero),
        None,
    )
}
with_h!(heteroacyclic_h_atom_matcher, heteroacyclic_atom_matcher);

pub fn alkoxyacyclic_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let mol = with_rings(mol, false);
    if mol.atom_degree(atom) != 2 || mol.atom(atom).atomic_num() != 8 {
        return false;
    }
    let Some(nnbr) = mol.atom_neighbors(atom).find(|&n| !ignore[n as usize]) else {
        return false;
    };
    let bonds = chain_single_bond(&mol);
    all_atoms_match(
        &mol,
        nnbr,
        ignore.to_vec(),
        Some(&aliphatic_c_or_h),
        Some(&bonds),
        Some(&is_carbon),
        None,
    )
}
with_h!(alkoxyacyclic_h_atom_matcher, alkoxyacyclic_atom_matcher);

/// Alkenyl/alkynyl: carbon and hydrogen only, with at least one bond of
/// `extra` type and no ring bonds.
fn unsaturated_alk_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool], extra: BondType) -> bool {
    let mol = with_rings(mol, false);
    let bonds = |b: &Bond| {
        (b.bond_type() == BondType::Single || b.bond_type() == extra)
            && !b.is_aromatic()
            && !bond_in_ring(&mol, b)
    };
    let at_least = |b: &Bond| b.bond_type() == extra;
    all_atoms_match(
        &mol,
        atom,
        ignore.to_vec(),
        Some(&aliphatic_c_or_h),
        Some(&bonds),
        None,
        Some(&at_least),
    )
}

pub fn alkenyl_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    unsaturated_alk_atom_matcher(mol, atom, ignore, BondType::Double)
}
with_h!(alkenyl_h_atom_matcher, alkenyl_atom_matcher);

pub fn alkynyl_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    unsaturated_alk_atom_matcher(mol, atom, ignore, BondType::Triple)
}
with_h!(alkynyl_h_atom_matcher, alkynyl_atom_matcher);

fn check_atom_ring(
    mol: &ROMol,
    atom: u32,
    ignore: &[bool],
    ring: &[i32],
    matcher: Option<AtomFn>,
    at_least_one: Option<AtomFn>,
) -> bool {
    let mut at_least = at_least_one.is_none();
    for &aidx in ring {
        let a = mol.atom(aidx as u32);
        if aidx as u32 != atom && (ignore[aidx as usize] || matcher.is_some_and(|m| !m(a))) {
            return false;
        }
        if !at_least && at_least_one.is_some_and(|f| f(a)) {
            at_least = true;
        }
    }
    at_least
}

fn check_bond_ring(
    mol: &ROMol,
    ring: &[i32],
    matcher: Option<BondFn>,
    at_least_one: Option<BondFn>,
) -> bool {
    let mut at_least = at_least_one.is_none();
    for &bidx in ring {
        let b = mol.bond(bidx as u32);
        if matcher.is_some_and(|m| !m(b)) {
            return false;
        }
        if !at_least && at_least_one.is_some_and(|f| f(b)) {
            at_least = true;
        }
    }
    at_least
}

/// The checks `FusedRingMatch` makes on the fused ring system of an atom.
#[derive(Default)]
struct FusedChecks<'a> {
    atom: Option<AtomFn<'a>>,
    bond: Option<BondFn<'a>>,
    one_atom_per_ring: Option<AtomFn<'a>>,
    one_bond_per_ring: Option<BondFn<'a>>,
    one_atom: Option<AtomFn<'a>>,
}

/// `FusedRingMatch`: all the rings fused (sharing at least two atoms) to the
/// ring of `atom` must pass the checks.
fn fused_ring_match(mol: &ROMol, atom: u32, ignore: &[bool], c: &FusedChecks) -> bool {
    if c.atom.is_some_and(|m| !m(mol.atom(atom))) {
        return false;
    }
    let mol = with_rings(mol, true);
    if num_atom_rings(&mol, atom) == 0 {
        return false;
    }
    let ri = mol.ring_info();
    let (atom_rings, bond_rings) = (ri.atom_rings(), ri.bond_rings());
    let mut ring_atoms: BTreeSet<i32> = BTreeSet::new();
    if let Some(i) = atom_rings.iter().position(|r| r.contains(&(atom as i32))) {
        let ring = &atom_rings[i];
        if !check_atom_ring(&mol, atom, ignore, ring, c.atom, c.one_atom_per_ring)
            || !check_bond_ring(&mol, &bond_rings[i], c.bond, c.one_bond_per_ring)
        {
            return false;
        }
        ring_atoms.extend(ring.iter().copied());
    }
    for (i, ring) in atom_rings.iter().enumerate() {
        let sring: BTreeSet<i32> = ring.iter().copied().collect();
        let diff: Vec<i32> = sring.difference(&ring_atoms).copied().collect();
        if diff.is_empty() || sring.len() - diff.len() < 2 {
            continue;
        }
        if !check_atom_ring(&mol, atom, ignore, ring, c.atom, c.one_atom_per_ring)
            || !check_bond_ring(&mol, &bond_rings[i], c.bond, c.one_bond_per_ring)
        {
            return false;
        }
        ring_atoms.extend(diff);
    }
    match c.one_atom {
        Some(f) => ring_atoms.iter().any(|&i| f(mol.atom(i as u32))),
        None => true,
    }
}

pub fn carbocycloalkyl_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let atom_ok = |at: &Atom| !at.is_aromatic() && at.atomic_num() == 6;
    let bond_ok = |b: &Bond| !b.is_aromatic() && b.bond_type() == BondType::Single;
    let checks = FusedChecks {
        atom: Some(&atom_ok),
        bond: Some(&bond_ok),
        ..Default::default()
    };
    fused_ring_match(mol, atom, ignore, &checks)
}
with_h!(carbocycloalkyl_h_atom_matcher, carbocycloalkyl_atom_matcher);

pub fn carbocycloalkenyl_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let at_least = |b: &Bond| is_aromatic_bond(b) || b.bond_type() == BondType::Double;
    let checks = FusedChecks {
        atom: Some(&is_carbon),
        one_bond_per_ring: Some(&at_least),
        ..Default::default()
    };
    fused_ring_match(mol, atom, ignore, &checks)
}
with_h!(
    carbocycloalkenyl_h_atom_matcher,
    carbocycloalkenyl_atom_matcher
);

pub fn carboaryl_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let atom_ok = |at: &Atom| at.is_aromatic() && at.atomic_num() == 6;
    let checks = FusedChecks {
        atom: Some(&atom_ok),
        bond: Some(&is_aromatic_bond),
        ..Default::default()
    };
    fused_ring_match(mol, atom, ignore, &checks)
}
with_h!(carboaryl_h_atom_matcher, carboaryl_atom_matcher);

pub fn carbocyclic_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let checks = FusedChecks {
        atom: Some(&is_carbon),
        ..Default::default()
    };
    fused_ring_match(mol, atom, ignore, &checks)
}
with_h!(carbocyclic_h_atom_matcher, carbocyclic_atom_matcher);

pub fn no_carbon_ring_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let not_c = |at: &Atom| at.atomic_num() != 6;
    let checks = FusedChecks {
        atom: Some(&not_c),
        ..Default::default()
    };
    fused_ring_match(mol, atom, ignore, &checks)
}
with_h!(no_carbon_ring_h_atom_matcher, no_carbon_ring_atom_matcher);

pub fn heterocyclic_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let hetero = |at: &Atom| !is_carbon_or_h(at);
    let checks = FusedChecks {
        one_atom: Some(&hetero),
        ..Default::default()
    };
    fused_ring_match(mol, atom, ignore, &checks)
}
with_h!(heterocyclic_h_atom_matcher, heterocyclic_atom_matcher);

pub fn heteroaryl_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let aromatic = |at: &Atom| at.is_aromatic();
    let hetero = |at: &Atom| !is_carbon_or_h(at);
    let checks = FusedChecks {
        atom: Some(&aromatic),
        bond: Some(&is_aromatic_bond),
        one_atom: Some(&hetero),
        ..Default::default()
    };
    fused_ring_match(mol, atom, ignore, &checks)
}
with_h!(heteroaryl_h_atom_matcher, heteroaryl_atom_matcher);

pub fn cyclic_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    let mol = with_rings(mol, false);
    let in_ring = |at: &Atom| num_atom_rings(&mol, at.idx()) > 0;
    let checks = FusedChecks {
        atom: Some(&in_ring),
        ..Default::default()
    };
    fused_ring_match(&mol, atom, ignore, &checks)
}
with_h!(cyclic_h_atom_matcher, cyclic_atom_matcher);

pub fn d_atom_matcher(mol: &ROMol, atom: u32, _ignore: &[bool]) -> bool {
    let a = mol.atom(atom);
    a.atomic_num() == 1 && a.isotope() == 2
}

pub fn t_atom_matcher(mol: &ROMol, atom: u32, _ignore: &[bool]) -> bool {
    let a = mol.atom(atom);
    a.atomic_num() == 1 && a.isotope() == 3
}

pub fn hplus_atom_matcher(mol: &ROMol, atom: u32, _ignore: &[bool]) -> bool {
    let a = mol.atom(atom);
    a.atomic_num() == 1 && a.formal_charge() == 1
}

pub fn pol_atom_matcher(mol: &ROMol, atom: u32, _ignore: &[bool]) -> bool {
    string_prop(mol.atom(atom), cp::atomLabel).is_some_and(|l| l == "Pol")
}

pub fn r_atom_matcher(mol: &ROMol, atom: u32, ignore: &[bool]) -> bool {
    group_h_atom_matcher(mol, atom, ignore)
}

fn string_prop(p: &impl RDProps, key: &str) -> Option<String> {
    p.get_prop(key).map(|v| v.to_rd_string())
}

/// `genericMatchers`: the matcher of a generic label.
pub fn generic_matcher(label: &str) -> Option<GenericMatcher> {
    Some(match label {
        "Group" | "G" => group_atom_matcher,
        "GroupH" | "GH" => group_h_atom_matcher,
        "Group*" | "G*" => group_star_atom_matcher,
        "GroupH*" | "GH*" => group_star_h_atom_matcher,
        "Alkyl" | "ALK" => alkyl_atom_matcher,
        "AlkylH" | "ALH" => alkyl_h_atom_matcher,
        "Alkenyl" | "AEL" => alkenyl_atom_matcher,
        "AlkenylH" | "AEH" => alkenyl_h_atom_matcher,
        "Alkynyl" | "AYL" => alkynyl_atom_matcher,
        "AlkynylH" | "AYH" => alkynyl_h_atom_matcher,
        "Carbocyclic" | "CBC" => carbocyclic_atom_matcher,
        "CarbocyclicH" | "CBH" => carbocyclic_h_atom_matcher,
        "Carbocycloalkyl" | "CAL" => carbocycloalkyl_atom_matcher,
        "CarbocycloalkylH" | "CAH" => carbocycloalkyl_h_atom_matcher,
        "Carbocycloalkenyl" | "CEL" => carbocycloalkenyl_atom_matcher,
        "CarbocycloalkenylH" | "CEH" => carbocycloalkenyl_h_atom_matcher,
        "Carboaryl" | "ARY" => carboaryl_atom_matcher,
        "CarboarylH" | "ARH" => carboaryl_h_atom_matcher,
        "Cyclic" | "CYC" => cyclic_atom_matcher,
        "CyclicH" | "CYH" => cyclic_h_atom_matcher,
        "Acyclic" | "ACY" => acyclic_atom_matcher,
        "AcyclicH" | "ACH" => acyclic_h_atom_matcher,
        "Carboacyclic" | "ABC" => carboacyclic_atom_matcher,
        "CarboacyclicH" | "ABH" => carboacyclic_h_atom_matcher,
        "Heteroacyclic" | "AHC" => heteroacyclic_atom_matcher,
        "HeteroacyclicH" | "AHH" => heteroacyclic_h_atom_matcher,
        "Alkoxy" | "AOX" => alkoxyacyclic_atom_matcher,
        "AlkoxyH" | "AOH" => alkoxyacyclic_h_atom_matcher,
        "Heterocyclic" | "CHC" => heterocyclic_atom_matcher,
        "HeterocyclicH" | "CHH" => heterocyclic_h_atom_matcher,
        "Heteroaryl" | "HAR" => heteroaryl_atom_matcher,
        "HeteroarylH" | "HAH" => heteroaryl_h_atom_matcher,
        "NoCarbonRing" | "CXX" => no_carbon_ring_atom_matcher,
        "NoCarbonRingH" | "CXH" => no_carbon_ring_h_atom_matcher,
        _ => return None,
    })
}

/// `genericAtomMatcher`: false if any of the generic atoms of `query` is not
/// satisfied by the atoms of `mol` it is matched to.
pub fn generic_atom_matcher(mol: &ROMol, query: &ROMol, matched: &[u32]) -> bool {
    let mut ignore = vec![false; mol.num_atoms() as usize];
    for &idx in matched {
        ignore[idx as usize] = true;
    }
    // perceive rings once for all the matchers
    let mol = with_rings(mol, true);
    for atom in query.atoms() {
        if query.atom_degree(atom.idx()) > 1 {
            continue;
        }
        let Some(label) = string_prop(atom, cp::_QueryAtomGenericLabel) else {
            continue;
        };
        if let Some(m) = generic_matcher(&label) {
            if !m(&mol, matched[atom.idx() as usize], &ignore) {
                return false;
            }
        }
    }
    true
}

/// `adjustQueryPropertiesWithGenericGroups`
pub fn adjust_query_properties_with_generic_groups(
    mol: &ROMol,
    params: &AdjustQueryParameters,
) -> Result<RWMol> {
    let mut res = mol.clone();
    adjust_query_properties_in_place(&mut res, params)?;
    set_generic_queries_from_properties(&mut res, true, true);
    Ok(res)
}

/// `convertGenericQueriesToSubstanceGroups`
pub fn convert_generic_queries_to_substance_groups(mol: &mut ROMol) -> Result<()> {
    for idx in 0..mol.num_atoms() {
        let Some(label) = string_prop(mol.atom(idx), cp::_QueryAtomGenericLabel) else {
            continue;
        };
        let mut sg = SubstanceGroup::new("SUP");
        sg.set_prop("LABEL", label, false)?;
        sg.add_atom_with_idx(mol, idx)?;
        add_substance_group(mol, sg);
        mol.atom_mut(idx).clear_prop(cp::_QueryAtomGenericLabel);
    }
    Ok(())
}

/// `setGenericQueriesFromProperties`: tags atoms with generic labels taken
/// from their `atomLabel` properties and/or from SUP substance groups.
pub fn set_generic_queries_from_properties(
    mol: &mut ROMol,
    use_atom_labels: bool,
    use_sgroups: bool,
) {
    if use_atom_labels {
        for idx in 0..mol.num_atoms() {
            let Some(mut label) = string_prop(mol.atom(idx), cp::atomLabel) else {
                continue;
            };
            let atom = mol.atom_mut(idx);
            // special case for generic groups from mol blocks
            if atom.atomic_num() == 0
                && atom
                    .query()
                    .is_some_and(|q| q.description() == "AtomAtomicNum" && !q.negate)
            {
                atom.set_query(Some(make_atom_null_query()));
            }
            // pseudoatom labels from CXSMILES end with "_p"
            if label.len() > 4 && label.ends_with("_p") {
                label.truncate(label.len() - 2);
            }
            if generic_matcher(&label).is_some() {
                // `label` is a valid property value
                let _ = atom.set_prop(cp::_QueryAtomGenericLabel, label, false);
                atom.clear_prop(cp::atomLabel);
            }
        }
    }
    if use_sgroups {
        let mut tagged: Vec<(u32, String)> = Vec::new();
        mol.substance_groups_mut().retain(|sg| {
            if string_prop(sg, "TYPE").as_deref() != Some("SUP") {
                return true;
            }
            match string_prop(sg, "LABEL") {
                Some(label) if generic_matcher(&label).is_some() => {
                    tagged.extend(sg.atoms().iter().map(|&a| (a, label.clone())));
                    false
                }
                _ => true,
            }
        });
        for (idx, label) in tagged {
            let _ = mol
                .atom_mut(idx)
                .set_prop(cp::_QueryAtomGenericLabel, label, false);
        }
    }
}
