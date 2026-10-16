//! SLN attributes of atoms, bonds and molecules.
//!
//! Ported from RDKit `Code/GraphMol/SLNParse/SLNAttribs.{h,cpp}` (2026.09.1).
//!
//! Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::VecDeque;

use super::SlnError;
use crate::graphmol::atom::{Atom, ChiralType};
use crate::graphmol::bond::{Bond, BondType};
use crate::graphmol::queries::{
    AtomFn, AtomQuery, CompositeQueryType, Query, QueryKind, make_atom_in_ring_query,
    make_bond_is_in_ring_query, make_bond_order_equals_query,
};
use crate::graphmol::query_ops::{atom_expand_query, bond_expand_query};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::{RDProps, RDValue, common_properties as cp};

/// `AttribCombineOp`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AttribCombineOp {
    LowPriAnd,
    Or,
    And,
    Not,
}

/// `AttribType`
#[derive(Debug, Clone, Default)]
pub struct AttribType {
    pub first: String,
    pub second: String,
    pub op: String,
    pub negated: bool,
    /// The (OR of the) recursive queries of an `is=`/`not=` attribute.
    pub struct_query: Option<AtomQuery>,
}

pub type AttribListType = Vec<(AttribCombineOp, AttribType)>;

fn sln_err(msg: String) -> SlnError {
    SlnError::Sln(msg)
}

/// `parseIntAttribVal`
fn parse_int_attrib_val(name: &str, val: &str) -> Result<i32, SlnError> {
    // (`boost::lexical_cast<int>`)
    val.to_lowercase().parse::<i32>().map_err(|_| {
        sln_err(format!(
            "SLN Parser error: bad integer value ({}) provided for property: {name}",
            val.to_lowercase()
        ))
    })
}

/// `makeQueryFromOp`
pub fn make_query_from_op(
    op: &str,
    val: i32,
    func: AtomFn,
    description: &str,
) -> Result<AtomQuery, SlnError> {
    // (the RDKit `Greater` query is true if the value to compare to is
    // greater than the target, so the sense is reversed here)
    let (kind, negate) = match op {
        "=" => (QueryKind::Equality, false),
        "!=" => (QueryKind::Equality, true),
        ">" => (QueryKind::Less, false),
        ">=" => (QueryKind::LessEqual, false),
        "<" => (QueryKind::Greater, false),
        "<=" => (QueryKind::GreaterEqual, false),
        _ => {
            return Err(sln_err(format!(
                "SLN Parser error: bad attribute operator ({op}) provided."
            )));
        }
    };
    let mut q = Query::with(kind, val, func, description);
    q.negate = negate;
    Ok(q)
}

fn composite(kind: QueryKind, description: &str, children: Vec<AtomQuery>) -> AtomQuery {
    let mut q = Query::new(kind);
    q.description = description.into();
    q.children = children;
    q
}

/// `parseAtomAttribs`: sets the properties/queries given by the attributes.
/// Some query values can't be set until the molecule is complete; those are
/// handled by [`parse_final_atom_attribs`].
pub fn parse_atom_attribs(
    atom: &mut Atom,
    attribs: AttribListType,
    doing_query: bool,
) -> Result<(), SlnError> {
    let mut atom_query: Option<AtomQuery> = None;
    let mut last_was_low_pri_and = false;
    for (how, mut attrib) in attribs {
        let mut query: Option<AtomQuery> = None;
        let attrib_name = attrib.first.to_lowercase();
        let attrib_val = attrib.second.clone();
        match attrib_name.as_str() {
            "charge" => {
                let chg = match attrib_val.as_str() {
                    "-" => -1,
                    "+" => 1,
                    _ => parse_int_attrib_val(&attrib_name, &attrib_val)?,
                };
                if !doing_query {
                    atom.set_formal_charge(chg);
                } else {
                    query = Some(make_query_from_op(
                        &attrib.op,
                        chg,
                        AtomFn::FormalCharge,
                        "AtomFormalCharge",
                    )?);
                }
            }
            "i" => {
                let val = parse_int_attrib_val(&attrib_name, &attrib_val)?;
                if !doing_query {
                    atom.set_isotope(val as u32);
                } else {
                    query = Some(make_query_from_op(
                        &attrib.op,
                        val,
                        AtomFn::Isotope,
                        "AtomIsotope",
                    )?);
                }
            }
            "r" => {
                if !attrib_val.is_empty() {
                    crate::rd_log_warning!("Query value '{attrib_val}' ignored for r query\n");
                }
                if !doing_query {
                    crate::rd_log_warning!(
                        "Query property '{attrib_name}' ignored on non-query atom\n"
                    );
                } else {
                    query = Some(make_atom_in_ring_query());
                }
            }
            "is" => {
                // recursive queries
                match attrib.struct_query.take() {
                    Some(q) => query = Some(q),
                    None => return Err(sln_err("failed recursive query".into())),
                }
            }
            "s" => {
                if attrib.op != "=" {
                    return Err(sln_err(format!(
                        "SLN Parser error: comparison operator '{}' not supported for chiral attributes.\n",
                        attrib.op
                    )));
                }
                let val = attrib_val.to_lowercase();
                let b = val.as_bytes();
                if b.first().is_some_and(|&c| c == b'i' || c == b'n') {
                    if b.get(1) == Some(&b'*') {
                        crate::rd_log_warning!(
                            "Chiral modifier * ignored, chiral spec {} will be used\n",
                            b[0] as char
                        );
                    }
                    if b.get(1) == Some(&b'm') {
                        crate::rd_log_warning!(
                            "Chiral modifier m ignored, chiral spec {} will be used\n",
                            b[0] as char
                        );
                    }
                } else {
                    crate::rd_log_warning!(
                        "Unsupported stereochemistry specifier '{val}' ignored.\n"
                    );
                }
            }
            _ => {
                // a block of properties that can have "f" values, and so need
                // special handling
                let mut f_tag = "";
                let mut val = 0;
                if attrib_val == "f" || attrib_name == "f" {
                    f_tag = "_SLN_";
                    atom.set_prop(cp::_Unfinished_SLN_, 1i32, false)
                        .map_err(SlnError::Other)?;
                    val = -666;
                }
                let spec: Option<(AtomFn, &str, bool)> = match attrib_name.as_str() {
                    "rbc" => Some((AtomFn::RingBondCount, "AtomRingBondCount", false)),
                    "tbo" => Some((AtomFn::TotalValence, "AtomTotalValence", false)),
                    "tac" => Some((AtomFn::TotalDegree, "AtomTotalDegree", false)),
                    "hc" => Some((AtomFn::HCount, "AtomHCount", false)),
                    "hac" => Some((AtomFn::NonHydrogenDegree, "AtomHeavyAtomDegree", false)),
                    "f" => Some((AtomFn::AllBondProduct, "AtomBondEnvironment", true)),
                    _ => None,
                };
                if let Some((func, descr, always_equal)) = spec {
                    if f_tag.is_empty() {
                        val = parse_int_attrib_val(&attrib_name, &attrib_val)?;
                    }
                    let op = if always_equal {
                        "="
                    } else {
                        attrib.op.as_str()
                    };
                    query = Some(make_query_from_op(
                        op,
                        val,
                        func,
                        &format!("{f_tag}{descr}"),
                    )?);
                } else {
                    // anything we don't know how to deal with is stored raw
                    atom.set_prop(&attrib_name, attrib_val.clone(), false)
                        .map_err(SlnError::Other)?;
                }
            }
        }

        // if we've constructed a query from all that, then we need to add it
        // to the atom query
        if let Some(mut query) = query {
            if !doing_query {
                crate::rd_log_warning!(
                    "Query property '{attrib_name}' ignored on non-query atom\n"
                );
            } else {
                if attrib.negated {
                    query.negate = !query.negate;
                }
                atom_query = Some(match atom_query.take() {
                    None => query,
                    Some(prev) => match how {
                        AttribCombineOp::And => {
                            // high-priority and
                            last_was_low_pri_and = false;
                            composite(QueryKind::And, "AtomAnd", vec![prev, query])
                        }
                        AttribCombineOp::LowPriAnd => {
                            last_was_low_pri_and = true;
                            composite(QueryKind::And, "AtomAnd", vec![prev, query])
                        }
                        AttribCombineOp::Or => {
                            let res = if last_was_low_pri_and {
                                // if the last query was a low-priority AND the
                                // tree is restructured: (a AND b) OR c becomes
                                // a AND (b OR c). (Upstream never adds c to
                                // the new OR.)
                                let mut children = prev.children.into_iter();
                                let a = children.next();
                                let b = children.next();
                                let or =
                                    composite(QueryKind::Or, "AtomOr", b.into_iter().collect());
                                composite(
                                    QueryKind::And,
                                    "AtomAnd",
                                    a.into_iter().chain(Some(or)).collect(),
                                )
                            } else {
                                composite(QueryKind::Or, "AtomOr", vec![prev, query])
                            };
                            last_was_low_pri_and = false;
                            res
                        }
                        AttribCombineOp::Not => {
                            return Err(sln_err("unrecognized query composition operator".into()));
                        }
                    },
                });
            }
        }
    }
    if let Some(q) = atom_query {
        atom_expand_query(atom, q, CompositeQueryType::And, true).map_err(SlnError::Other)?;
    }
    Ok(())
}

/// `parseFinalAtomAttribs`: finalizes the attributes that had "f" in the
/// SLN (their query descriptions start with "_SLN_") once the molecule is
/// complete.
pub fn parse_final_atom_attribs(
    mol: &mut ROMol,
    idx: u32,
    doing_query: bool,
) -> Result<(), SlnError> {
    let atom = mol.atom(idx);
    if !doing_query || !atom.has_query() || !atom.has_prop(cp::_Unfinished_SLN_) {
        return Ok(());
    }
    mol.atom_mut(idx).clear_prop(cp::_Unfinished_SLN_);
    let mut query = mol.atom_mut(idx).take_query().expect("checked above");
    let res = (|| -> Result<(), SlnError> {
        let mut todo: VecDeque<&mut AtomQuery> = VecDeque::new();
        todo.push_back(&mut query);
        while let Some(q) = todo.pop_front() {
            if q.description.len() > 5 && q.description.starts_with("_SLN_") {
                q.description = q.description[5..].to_string();
                if let Some(f) = &q.data_func {
                    use crate::graphmol::queries::QueryTarget;
                    q.val = f.eval(mol, idx).map_err(SlnError::Other)?;
                }
            }
            todo.extend(q.children.iter_mut());
        }
        Ok(())
    })();
    mol.atom_mut(idx).set_query(Some(query));
    res
}

/// `parseBondAttribs`
pub fn parse_bond_attribs(
    bond: &mut Bond,
    attribs: AttribListType,
    doing_query: bool,
) -> Result<(), SlnError> {
    // FIX (upstream): the query tree reordering of atoms is missing here
    let mut seen_type_query = false;
    for (op, attrib) in attribs {
        let how = match op {
            AttribCombineOp::And | AttribCombineOp::LowPriAnd => CompositeQueryType::And,
            AttribCombineOp::Or => CompositeQueryType::Or,
            AttribCombineOp::Not => {
                return Err(sln_err("unrecognized query composition operator".into()));
            }
        };
        let attrib_name = attrib.first.to_lowercase();
        let attrib_val = attrib.second.clone();
        if attrib_name == "type" {
            let val = attrib_val.to_lowercase();
            let bond_type = match val.as_str() {
                "-" | "1" => BondType::Single,
                "=" | "2" => BondType::Double,
                "#" | "3" => BondType::Triple,
                ":" | "aromatic" => BondType::Aromatic,
                _ => {
                    bond.set_prop("SLN_Type", val.clone(), false)
                        .map_err(SlnError::Other)?;
                    BondType::Other
                }
            };
            if !doing_query {
                bond.set_bond_type(bond_type);
            } else {
                let mut query = make_bond_order_equals_query(bond_type);
                if attrib.negated {
                    query.negate = !query.negate;
                }
                if seen_type_query {
                    bond_expand_query(bond, query, how, true).map_err(SlnError::Other)?;
                } else {
                    // if this is the first type query, any existing bond
                    // order queries are replaced
                    // FIX (upstream): this replaces too much, ring queries
                    // also get blown out
                    bond.set_query(Some(query));
                }
                seen_type_query = true;
            }
        } else if attrib_name == "r" {
            if !attrib_val.is_empty() {
                crate::rd_log_warning!("Query value '{attrib_val}' ignored for r query\n");
            }
            if !doing_query {
                crate::rd_log_warning!(
                    "Query property '{attrib_name}' ignored on non-query bond\n"
                );
            } else {
                let mut query = make_bond_is_in_ring_query();
                if attrib.negated {
                    query.negate = true;
                }
                bond_expand_query(bond, query, how, true).map_err(SlnError::Other)?;
            }
        } else {
            // anything we don't know how to deal with is stored raw
            bond.set_prop(&attrib_name, attrib_val, false)
                .map_err(SlnError::Other)?;
        }
    }
    Ok(())
}

/// `parseMolAttribs`
pub fn parse_mol_attribs(mol: &mut ROMol, attribs: AttribListType) -> Result<(), SlnError> {
    for (op, attrib) in attribs {
        assert_eq!(op, AttribCombineOp::And, "bad attrib type");
        let name = attrib.first.to_lowercase();
        let mut val = attrib.second;
        if val.starts_with('"') && val.ends_with('"') {
            // (upstream erases the first and the last character, also of a
            // single `"`)
            val.remove(0);
            val.pop();
        }
        let key = if name == "name" {
            cp::_Name
        } else {
            name.as_str()
        };
        mol.set_prop(key, val, false).map_err(SlnError::Other)?;
    }
    Ok(())
}

/// `adjustAtomChiralities`: translates the SLN chirality marks to RDKit
/// chiral tags.
pub fn adjust_atom_chiralities(mol: &mut ROMol) -> Result<(), SlnError> {
    for idx in 0..mol.num_atoms() {
        let Some(RDValue::String(attrib_val)) = mol.atom(idx).get_prop(cp::_SLN_s).cloned() else {
            continue;
        };
        // the SLN chirality maps straight to the RDKit one; a few pairs:
        //   C[s=n]H(Cl)(F)Br  <->  [C@@H](Cl)(F)Br (CHI_TETRAHEDRAL_CW)
        //   ClC[s=n]H(F)Br  <->  Cl[C@H](F)Br (CHI_TETRAHEDRAL_CCW)
        //   FC[1:s=n](Cl)OCH2@1   <->  F[C@@]1(Cl)OC1 (CHI_TETRAHEDRAL_CW)
        match attrib_val.as_bytes().first() {
            Some(b'n') => mol.atom_mut(idx).set_chiral_tag(ChiralType::TetrahedralCW),
            Some(b'i') => mol.atom_mut(idx).set_chiral_tag(ChiralType::TetrahedralCCW),
            _ => {}
        }
        let mut neighbors: Vec<(i32, i32)> = mol
            .atom_bond_indices(idx)
            .map(|b| (mol.bond(b).other_atom_idx(idx).map(|a| a as i32), b as i32))
            .map(|(a, b)| a.map(|a| (a, b)))
            .collect::<Result<_, _>>()
            .map_err(SlnError::Other)?;
        // sort by neighbour index
        neighbors.sort();
        let bond_ordering: Vec<i32> = neighbors.iter().map(|n| n.1).collect();
        // the permutation order relative to the atom numbering (sln chirality)
        let n_swaps = mol
            .perturbation_order(idx, &bond_ordering)
            .map_err(SlnError::Other)?;
        if n_swaps % 2 != 0 {
            let flipped = if mol.atom(idx).chiral_tag() == ChiralType::TetrahedralCW {
                ChiralType::TetrahedralCCW
            } else {
                ChiralType::TetrahedralCW
            };
            mol.atom_mut(idx).set_chiral_tag(flipped);
        }
    }
    Ok(())
}
