//! MolJSON parser (upstream `Parser.cpp`).

use std::collections::BTreeSet;
use std::sync::Arc;

use super::{
    CURRENT_CHARGE_REPRESENTATION_VERSION, CURRENT_MOL_JSON_VERSION,
    CURRENT_QUERY_REPRESENTATION_VERSION, CURRENT_RDKIT_JSON_VERSION,
    CURRENT_RDKIT_REPRESENTATION_VERSION, JsonParseParameters, bond_stereo_from_name,
    bond_type_from_order, chiral_tag_from_name, stereo_group_type_from_name,
};
use crate::geometry::Point3D;
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::{Bond, BondType};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::pickle::{finalize_atom_query, finalize_bond_query, tags};
use crate::graphmol::queries::{AtomQuery, BondQuery, Query, QueryKind, QueryMol, QueryTarget};
use crate::graphmol::ring_info::FindRingType;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::stereo_group::StereoGroup;
use crate::graphmol::substance_group::{
    SubstanceGroup, add_substance_group, checks as SubstanceGroupChecks,
};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::json::Json;
use crate::rdgeneral::threads::par_map;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

fn bad_format(msg: impl AsRef<str>) -> RDKitError {
    RDKitError::FileParse(format!("Bad Format: {}", msg.as_ref()))
}

/// The `defaults` objects of the document (empty when absent).
struct Defaults<'a>(Option<&'a Json>);

/// `getIntDefaultValue`: the member of `from`, else of the defaults, else 0.
fn int_default(key: &str, from: &Json, defaults: &Defaults<'_>) -> Result<i64> {
    let not_int = || RDKitError::FileParse(format!("Bad format: value of {key} is not an int"));
    for source in [Some(from), defaults.0] {
        if let Some(val) = source
            .and_then(|s| s.as_object().ok())
            .and_then(|o| o.get(key))
        {
            return match val {
                Json::Int(_) | Json::UInt(_) => val.as_i64(),
                _ => Err(not_int()),
            };
        }
    }
    Ok(0)
}

/// `getStringDefaultValue`
fn string_default(key: &str, from: &Json, defaults: &Defaults<'_>) -> Result<String> {
    for source in [Some(from), defaults.0] {
        if let Some(val) = source
            .and_then(|s| s.as_object().ok())
            .and_then(|o| o.get(key))
        {
            return match val {
                Json::String(s) => Ok(s.clone()),
                _ => Err(RDKitError::FileParse(format!(
                    "Bad format: value of {key} is not a string"
                ))),
            };
        }
    }
    Ok(String::new())
}

/// Whether the object `val` has the member `key`.
fn has(val: &Json, key: &str) -> Result<bool> {
    Ok(val.as_object()?.contains(key))
}

fn as_u32(val: &Json) -> Result<u32> {
    Ok(val.as_i64()? as u32)
}

fn indices(val: &Json) -> Result<Vec<u32>> {
    val.as_array()?.iter().map(as_u32).collect()
}

fn read_atom(
    mol: &mut RWMol,
    val: &Json,
    defaults: &Defaults<'_>,
    params: &JsonParseParameters,
) -> Result<()> {
    let stereo = string_default("stereo", val, defaults)?;
    let chiral_tag =
        chiral_tag_from_name(&stereo).ok_or_else(|| bad_format("bad stereo value for atom"))?;
    let z = int_default("z", val, defaults)?;
    if !(0..=i64::from(PeriodicTable::get().get_max_atomic_number())).contains(&z) {
        return Err(RDKitError::precondition("Atomic number not found"));
    }
    let mut atom = Atom::new(z as u32);
    if params.use_h_counts {
        atom.set_no_implicit(true);
        atom.set_num_explicit_hs(int_default("impHs", val, defaults)? as u32);
    }
    atom.set_formal_charge(int_default("chg", val, defaults)? as i32);
    atom.set_num_radical_electrons(int_default("nRad", val, defaults)? as u32);
    atom.set_isotope(int_default("isotope", val, defaults)? as u32);
    atom.set_chiral_tag(chiral_tag);
    mol.add_atom(atom, false);
    Ok(())
}

/// Returns whether the bond carries stereo information to read later.
fn read_bond(mol: &mut RWMol, val: &Json, defaults: &Defaults<'_>) -> Result<bool> {
    let order = int_default("bo", val, defaults)?;
    let bond_type =
        bond_type_from_order(order).ok_or_else(|| bad_format("bad bond order for bond"))?;
    let ids = val.at("atoms")?;
    let begin = ids.at_idx(0)?.as_i64()? as i32;
    let end = ids.at_idx(1)?.as_i64()? as i32;
    let mut bond = Bond::new(bond_type);
    for (name, idx, line) in [
        ("bond_pin->getBeginAtomIdx()", begin, 408),
        ("bond_pin->getEndAtomIdx()", end, 409),
    ] {
        if idx < 0 || idx as u32 >= mol.num_atoms() {
            return Err(RDKitError::range_at(
                name,
                idx as u32,
                mol.num_atoms(),
                "Code/GraphMol/ROMol.cpp",
                line,
            ));
        }
    }
    bond.set_begin_atom_idx_unchecked(begin as u32);
    bond.set_end_atom_idx_unchecked(end as u32);
    mol.add_bond_object(bond)?;
    Ok(string_default("stereo", val, defaults)? != "unspecified")
}

fn read_bond_stereo(mol: &mut RWMol, idx: u32, val: &Json, defaults: &Defaults<'_>) -> Result<()> {
    let stereo = string_default("stereo", val, defaults)?;
    if stereo == "unspecified" {
        return Ok(());
    }
    let bond_stereo =
        bond_stereo_from_name(&stereo).ok_or_else(|| bad_format("bond stereo value for bond"))?;
    if has(val, "stereoAtoms")? {
        let ids = val.at("stereoAtoms")?.as_array()?;
        let stereo_atom = |i: usize| -> Result<u32> {
            // upstream indexes the array unchecked; a missing entry is
            // not an int
            let a = ids
                .get(i)
                .ok_or_else(|| RDKitError::Runtime("value is not a std::int64_t number".into()))?
                .as_i64()? as u32;
            if a >= mol.num_atoms() {
                return Err(RDKitError::range_at(
                    "idx2",
                    a,
                    mol.num_atoms(),
                    "Code/GraphMol/ROMol.cpp",
                    345,
                ));
            }
            Ok(a)
        };
        let (a, b) = (stereo_atom(0)?, stereo_atom(1)?);
        mol.set_bond_stereo_atoms(idx, a, b)?;
    } else if stereo == "cis" || stereo == "trans" {
        return Err(bad_format("bond stereo provided without stereoAtoms"));
    }
    mol.bond_mut(idx).set_stereo(bond_stereo)
}

/// `parseProperties`: ints, doubles and strings become properties, other
/// JSON values are ignored.
fn parse_properties<T: RDProps>(obj: &mut T, props: &Json) -> Result<()> {
    for (key, val) in props.as_object()?.iter() {
        match val {
            Json::Int(_) | Json::UInt(_) => obj.set_prop(key, val.as_i64()? as i32, false)?,
            Json::Double(d) => obj.set_prop(key, *d, false)?,
            Json::String(s) => obj.set_prop(key, s.as_str(), false)?,
            _ => {}
        }
    }
    Ok(())
}

fn read_stereo_groups(mol: &mut RWMol, groups: &Json) -> Result<()> {
    let mut mol_groups = mol.stereo_groups().to_vec();
    for val in groups.as_array()? {
        if !has(val, "type")? {
            return Err(bad_format("stereogroup does not have a type"));
        }
        if !has(val, "atoms")? && !has(val, "bonds")? {
            return Err(bad_format(
                "stereogroup does not have either atoms or bonds",
            ));
        }
        let group_type = stereo_group_type_from_name(val.at("type")?.as_str()?)
            .ok_or_else(|| bad_format("bad stereogroup type"))?;
        let mut group_id = 0;
        if group_type != crate::graphmol::stereo_group::StereoGroupType::Absolute && has(val, "id")?
        {
            group_id = as_u32(val.at("id")?)?;
        }
        let atoms = if has(val, "atoms")? {
            indices(val.at("atoms")?)?
        } else {
            Vec::new()
        };
        let bonds = if has(val, "bonds")? {
            indices(val.at("bonds")?)?
        } else {
            Vec::new()
        };
        for &a in &atoms {
            mol.get_atom(a)?;
        }
        for &b in &bonds {
            mol.get_bond(b)?;
        }
        if !atoms.is_empty() || !bonds.is_empty() {
            mol_groups.push(StereoGroup::new(group_type, atoms, bonds, group_id)?);
        }
    }
    mol.set_stereo_groups(mol_groups)
}

fn read_point(val: &Json, what: &str) -> Result<Point3D> {
    let p = val.as_array()?;
    if p.len() != 3 {
        return Err(bad_format(format!("{what} doesn't have three coordinates")));
    }
    Ok(Point3D::new(p[0].as_f64()?, p[1].as_f64()?, p[2].as_f64()?))
}

fn read_substance_groups(mol: &mut RWMol, groups: &Json) -> Result<()> {
    for val in groups.as_array()? {
        let sg_type = match val.as_object()?.get("properties") {
            Some(props) if has(props, "TYPE")? => props.at("TYPE")?.as_str()?.to_string(),
            _ => return Err(bad_format("substance group does not have TYPE property")),
        };
        if !SubstanceGroupChecks::is_valid_type(&sg_type) {
            return Err(bad_format(format!(
                "substance group TYPE '{sg_type}' not recognized"
            )));
        }
        let mut sg = SubstanceGroup::new(&sg_type);
        parse_properties(&mut sg, val.at("properties")?)?;
        if let Some(RDValue::String(sub)) = sg.get_prop("SUBTYPE") {
            if !SubstanceGroupChecks::is_valid_sub_type(sub) {
                return Err(bad_format(format!(
                    "substance group SUBTYPE '{sub}' not recognized"
                )));
            }
        }
        if let Some(RDValue::String(connect)) = sg.get_prop("CONNECT") {
            if !SubstanceGroupChecks::is_valid_connect_type(connect) {
                return Err(bad_format(format!(
                    "substance group CONNECT type '{connect}' not recognized"
                )));
            }
        }
        if has(val, "atoms")? {
            sg.set_atoms(mol, indices(val.at("atoms")?)?)?;
        }
        if has(val, "bonds")? {
            sg.set_bonds(mol, indices(val.at("bonds")?)?)?;
        }
        if has(val, "parentAtoms")? {
            sg.set_parent_atoms(mol, indices(val.at("parentAtoms")?)?)?;
        }
        if has(val, "brackets")? {
            for bracket in val.at("brackets")?.as_array()? {
                let mut points = [Point3D::default(); 3];
                for (slot, pt) in bracket.as_array()?.iter().enumerate() {
                    *points
                        .get_mut(slot)
                        .ok_or_else(|| RDKitError::Runtime("out of range".into()))? =
                        read_point(pt, "bracket point")?;
                }
                sg.add_bracket(points);
            }
        }
        if has(val, "cstates")? {
            for cstate in val.at("cstates")?.as_array()? {
                let bond_idx = as_u32(cstate.at("bond")?)?;
                let vector = if has(cstate, "vector")? {
                    read_point(cstate.at("vector")?, "cstate vector")?
                } else {
                    Point3D::default()
                };
                sg.cstates_mut()
                    .push(crate::graphmol::substance_group::CState { bond_idx, vector });
            }
        }
        if has(val, "attachPoints")? {
            for ap in val.at("attachPoints")?.as_array()? {
                let a_idx = as_u32(ap.at("aIdx")?)?;
                let lv_idx = if has(ap, "lvIdx")? {
                    as_u32(ap.at("lvIdx")?)? as i32
                } else {
                    -1
                };
                let id = if has(ap, "id")? {
                    ap.at("id")?.as_str()?.to_string()
                } else {
                    String::new()
                };
                sg.attach_points_mut()
                    .push(crate::graphmol::substance_group::AttachPoint { a_idx, lv_idx, id });
            }
        }
        add_substance_group(mol, sg);
    }
    Ok(())
}

fn read_conformer(num_atoms: usize, val: &Json) -> Result<Conformer> {
    if !has(val, "dim")? {
        return Err(bad_format("no conformer dimension"));
    }
    let dim = val.at("dim")?.as_i64()? as usize;
    let mut conf = Conformer::new(num_atoms);
    match dim {
        2 => conf.set_3d(false),
        3 => conf.set_3d(true),
        _ => return Err(bad_format("conformer dimension != 2 or 3")),
    }
    if !has(val, "coords")? {
        return Err(bad_format("no conformer coords"));
    }
    let mut count = 0;
    for pt in val.at("coords")?.as_array()? {
        let coords = pt.as_array()?;
        if coords.len() != dim {
            return Err(RDKitError::FileParse(
                "coordinate contains wrong number of values".into(),
            ));
        }
        let z = if dim == 3 { coords[2].as_f64()? } else { 0.0 };
        let p = Point3D::new(coords[0].as_f64()?, coords[1].as_f64()?, z);
        conf.set_atom_pos(count, p)?;
        count += 1;
    }
    // `setAtomPos` grows the conformer, so too many positions only fail once
    // the conformer is added to the molecule
    if count != conf.num_atoms() {
        return Err(bad_format(
            "conformer doesn't contain coordinates for all atoms",
        ));
    }
    Ok(conf)
}

fn read_partial_charges(mol: &mut RWMol, rep: &Json) -> Result<()> {
    if !has(rep, "formatVersion")? {
        return Err(bad_format("missing version"));
    }
    let version = rep.at("formatVersion")?.as_i64()?;
    if version > CURRENT_CHARGE_REPRESENTATION_VERSION {
        crate::rd_log_warning!("partialCharges version {version} too recent. Ignoring it.");
        return Ok(());
    }
    if has(rep, "values")? {
        let values = rep.at("values")?.as_array()?;
        if values.len() != mol.num_atoms() as usize {
            return Err(bad_format("size of values array != num atoms"));
        }
        for (idx, val) in values.iter().enumerate() {
            let Json::Double(charge) = val else {
                return Err(bad_format("partial charge not double"));
            };
            mol.atom_mut(idx as u32)
                .set_prop(cp::_GasteigerCharge, *charge, false)?;
        }
    }
    Ok(())
}

// ------------------------------------------------------------------ queries

/// `readBaseQuery`: the query node for a tag, before it is finished.
fn read_base_query<F: QueryTarget>(val: &Json) -> Result<Query<F>> {
    let tag = val.at("tag")?.as_i64()?;
    if !has(val, "descr")? {
        return Err(bad_format("missing query description"));
    }
    let optional_tol = |q: &mut Query<F>| -> Result<()> {
        if has(val, "tol")? {
            q.tol = val.at("tol")?.as_i64()? as i32;
        }
        Ok(())
    };
    let val_query = |kind: QueryKind| -> Result<Query<F>> {
        let mut q = Query::new(kind);
        q.val = val.at("val")?.as_i64()? as i32;
        optional_tol(&mut q)?;
        Ok(q)
    };
    let tag = u8::try_from(tag).unwrap_or(u8::MAX);
    match tag {
        tags::QUERY_AND => Ok(Query::new(QueryKind::And)),
        tags::QUERY_OR => Ok(Query::new(QueryKind::Or)),
        tags::QUERY_XOR => Ok(Query::new(QueryKind::XOr)),
        tags::QUERY_EQUALS => val_query(QueryKind::Equality),
        tags::QUERY_GREATER => val_query(QueryKind::Greater),
        tags::QUERY_GREATEREQUAL => val_query(QueryKind::GreaterEqual),
        tags::QUERY_LESS => val_query(QueryKind::Less),
        tags::QUERY_LESSEQUAL => val_query(QueryKind::LessEqual),
        tags::QUERY_NULL => Ok(Query::new(QueryKind::Base { always_true: false })),
        tags::QUERY_RANGE => {
            let lower = val.at("lower")?.as_i64()? as i32;
            let upper = val.at("upper")?.as_i64()? as i32;
            let (mut lower_open, mut upper_open) = (false, false);
            if has(val, "ends")? {
                let ends = val.at("ends")?.as_i64()?;
                lower_open = ends & 0x2 != 0;
                upper_open = ends & 0x1 != 0;
            }
            let mut q = Query::new(QueryKind::Range {
                lower,
                upper,
                lower_open,
                upper_open,
            });
            optional_tol(&mut q)?;
            Ok(q)
        }
        tags::QUERY_SET => {
            let mut members = BTreeSet::new();
            if has(val, "set")? {
                for member in val.at("set")?.as_array()? {
                    members.insert(member.as_i64()? as i32);
                }
            }
            Ok(Query::new(QueryKind::Set(members)))
        }
        _ => Err(bad_format(format!(
            "unknown query tag {}",
            val.at("tag")?.as_i64()?
        ))),
    }
}

/// The parts of `finishQuery` that do not depend on the owner's type.
fn finish_common<F: QueryTarget>(q: &mut Query<F>, val: &Json) -> Result<()> {
    q.description = val.at("descr")?.as_str()?.to_string();
    if has(val, "type")? {
        q.type_label = val.at("type")?.as_str()?.to_string();
    }
    q.negate = has(val, "negated")? && val.at("negated")?.as_bool()?;
    Ok(())
}

fn read_children<F: QueryTarget>(
    q: &mut Query<F>,
    val: &Json,
    read: impl Fn(&Json) -> Result<Query<F>>,
) -> Result<()> {
    if has(val, "children")? {
        for child in val.at("children")?.as_array()? {
            q.children.push(read(child)?);
        }
    }
    Ok(())
}

fn read_atom_query(
    val: &Json,
    atom_defaults: &Defaults<'_>,
    bond_defaults: &Defaults<'_>,
    params: &JsonParseParameters,
) -> Result<AtomQuery> {
    if !has(val, "tag")? {
        return Err(bad_format("missing atom query tag"));
    }
    let tag = val.at("tag")?.as_i64()?;
    let mut q: AtomQuery = if tag == i64::from(tags::QUERY_RECURSIVE) {
        if !has(val, "subquery")? {
            return Err(bad_format("missing subquery"));
        }
        let mut sub = RWMol::new();
        process_mol(
            &mut sub,
            val.at("subquery")?,
            atom_defaults,
            bond_defaults,
            params,
        )?;
        Query::new(QueryKind::Recursive {
            mol: QueryMol(Arc::new(sub)),
            serial_number: 0,
            set: BTreeSet::new(),
            initialized: false,
        })
    } else if tag == i64::from(tags::QUERY_ATOMRING) {
        let mut q = Query::new(QueryKind::AtomRing);
        q.val = val.at("val")?.as_i64()? as i32;
        if has(val, "tol")? {
            q.tol = val.at("tol")?.as_i64()? as i32;
        }
        q.data_func = Some(crate::graphmol::queries::AtomFn::RingMembership);
        q
    } else {
        read_base_query(val)?
    };
    finish_common(&mut q, val)?;
    finalize_atom_query(&mut q)?;
    read_children(&mut q, val, |c| {
        read_atom_query(c, atom_defaults, bond_defaults, params)
    })?;
    Ok(q)
}

fn read_bond_query(val: &Json) -> Result<BondQuery> {
    if !has(val, "tag")? {
        return Err(bad_format("missing bond query tag"));
    }
    let mut q: BondQuery = read_base_query(val)?;
    finish_common(&mut q, val)?;
    finalize_bond_query(&mut q)?;
    read_children(&mut q, val, read_bond_query)?;
    Ok(q)
}

/// The entries of an `atomQueries`/`bondQueries` array that hold a query.
///
/// Entries without a `tag` are skipped *without* advancing the target index,
/// exactly as upstream does, so a molecule mixing query and plain atoms does
/// not round-trip.
fn query_entries<'a>(
    rep: &'a Json,
    key: &str,
    what: &str,
    count: u32,
) -> Result<Vec<(u32, &'a Json)>> {
    let Some(entries) = rep.as_object()?.get(key) else {
        return Ok(Vec::new());
    };
    let mut out = Vec::new();
    let mut idx = 0;
    for entry in entries.as_array()? {
        if !entry.is_object() {
            return Err(bad_format(format!("{what}Query not object")));
        }
        if !has(entry, "tag")? {
            continue;
        }
        if idx >= count {
            return Err(RDKitError::FileParse(format!("too much {what} data found")));
        }
        out.push((idx, entry));
        idx += 1;
    }
    Ok(out)
}

fn read_queries(
    mol: &mut RWMol,
    rep: &Json,
    atom_defaults: &Defaults<'_>,
    bond_defaults: &Defaults<'_>,
    params: &JsonParseParameters,
) -> Result<()> {
    if !has(rep, "formatVersion")? {
        return Err(bad_format("missing format_version"));
    }
    let version = rep.at("formatVersion")?.as_i64()?;
    if version > CURRENT_QUERY_REPRESENTATION_VERSION {
        crate::rd_log_warning!(
            "RDKit query representation format version {version} too recent. Ignoring it."
        );
        return Ok(());
    }
    for (idx, entry) in query_entries(rep, "atomQueries", "atom", mol.num_atoms())? {
        let q = read_atom_query(entry, atom_defaults, bond_defaults, params)?;
        mol.atom_mut(idx).set_query(Some(q));
    }
    for (idx, entry) in query_entries(rep, "bondQueries", "bond", mol.num_bonds())? {
        let q = read_bond_query(entry)?;
        mol.bond_mut(idx).set_query(Some(q));
    }
    Ok(())
}

// ------------------------------------------------------------------ molecules

fn read_rdkit_representation(
    mol: &mut RWMol,
    rep: &Json,
    params: &JsonParseParameters,
) -> Result<()> {
    if !has(rep, "formatVersion")? {
        return Err(bad_format("missing format_version"));
    }
    let version = rep.at("formatVersion")?.as_i64()?;
    if version > CURRENT_RDKIT_REPRESENTATION_VERSION {
        crate::rd_log_warning!(
            "RDKit representation format version {version} too recent. Ignoring it."
        );
        return Ok(());
    }
    let member = |key: &str| rep.as_object().map(|o| o.get(key));
    let int_entries = |key: &str, what: &str| -> Result<Vec<u32>> {
        let Some(vals) = member(key)? else {
            return Ok(Vec::new());
        };
        vals.as_array()?
            .iter()
            .map(|v| match v {
                Json::Int(_) | Json::UInt(_) => as_u32(v),
                _ => Err(bad_format(format!("{what} not int"))),
            })
            .collect()
    };
    for idx in int_entries("aromaticAtoms", "aromaticAtom")? {
        mol.get_atom_mut(idx)?.set_is_aromatic(true);
    }
    for idx in int_entries("aromaticBonds", "aromaticBond")? {
        let bond = mol.get_bond_mut(idx)?;
        bond.set_is_aromatic(true);
        if params.set_aromatic_bonds {
            bond.set_bond_type(BondType::Aromatic);
        }
    }
    for (idx, rank) in int_entries("cipRanks", "ciprank")?.into_iter().enumerate() {
        mol.get_atom_mut(idx as u32)?
            .set_prop(cp::_CIPRank, rank, false)?;
    }
    if let Some(codes) = member("cipCodes")? {
        for code in codes.as_array()? {
            if !code.is_array() {
                return Err(bad_format("CIPCode not string"));
            }
            let idx = code.at_idx(0)?.as_i64()? as u32;
            let label = code.at_idx(1)?.as_str()?;
            mol.get_atom_mut(idx)?
                .set_prop(cp::_CIPCode, label, false)?;
        }
    }
    if let Some(rings) = member("atomRings")? {
        if mol.ring_info().is_initialized() {
            return Err(RDKitError::invariant("rings already initialized"));
        }
        mol.ring_info_mut().initialize(FindRingType::OtherOrUnknown);
        for ring in rings.as_array()? {
            if !ring.is_array() {
                return Err(bad_format("atomRing not array"));
            }
            let atoms = ring
                .as_array()?
                .iter()
                .map(|a| a.as_i64().map(|a| a as i32))
                .collect::<Result<Vec<_>>>()?;
            let mut bonds = Vec::with_capacity(atoms.len());
            for (i, &a) in atoms.iter().enumerate() {
                let b = atoms[(i + 1) % atoms.len()];
                let Some(bond) = mol.bond_idx_between(a as u32, b as u32) else {
                    return Err(RDKitError::invariant("no bond found for ring"));
                };
                bonds.push(bond as i32);
            }
            mol.ring_info_mut().add_ring(atoms, bonds)?;
        }
    }
    Ok(())
}

/// `processMol`: fills `mol` from the JSON description `val`.
fn process_mol(
    mol: &mut RWMol,
    val: &Json,
    atom_defaults: &Defaults<'_>,
    bond_defaults: &Defaults<'_>,
    params: &JsonParseParameters,
) -> Result<()> {
    if has(val, "name")? {
        mol.set_prop(cp::_Name, val.at("name")?.as_str()?, false)?;
    }
    if !has(val, "atoms")? {
        return Err(bad_format("missing atoms in JSON"));
    }
    if !has(val, "bonds")? {
        return Err(bad_format("missing bonds in JSON"));
    }
    for atom in val.at("atoms")?.as_array()? {
        read_atom(mol, atom, atom_defaults, params)?;
    }
    let bonds = val.at("bonds")?.as_array()?;
    let mut need_stereo_loop = false;
    for bond in bonds {
        need_stereo_loop |= read_bond(mol, bond, bond_defaults)?;
    }
    if need_stereo_loop {
        // the stereo atoms need all the bonds to be there
        for (idx, bond) in bonds.iter().enumerate() {
            read_bond_stereo(mol, idx as u32, bond, bond_defaults)?;
        }
    }
    if has(val, "stereoGroups")? {
        read_stereo_groups(mol, val.at("stereoGroups")?)?;
    }
    if has(val, "substanceGroups")? {
        read_substance_groups(mol, val.at("substanceGroups")?)?;
    }
    if params.parse_conformers && has(val, "conformers")? {
        for conf in val.at("conformers")?.as_array()? {
            let conf = read_conformer(mol.num_atoms() as usize, conf)?;
            mol.add_conformer(conf, true)?;
        }
    }
    if params.parse_properties && has(val, "properties")? {
        parse_properties(mol, val.at("properties")?)?;
    }
    if has(val, "extensions")? {
        for ext in val.at("extensions")?.as_array()? {
            if !has(ext, "name")? {
                return Err(bad_format("representation has no name member"));
            }
            match ext.at("name")?.as_str()? {
                "rdkitRepresentation" => read_rdkit_representation(mol, ext, params)?,
                "partialCharges" => read_partial_charges(mol, ext)?,
                "rdkitQueries" => read_queries(mol, ext, atom_defaults, bond_defaults, params)?,
                _ => {}
            }
        }
    }
    mol.set_prop(cp::_StereochemDone, 1, false)
}

/// The checks of the document header.
fn check_header(doc: &Json) -> Result<()> {
    let obj = doc.as_object()?;
    if let Some(header) = obj.get("commonchem") {
        let Some(version) = header.as_object().ok().and_then(|h| h.get("version")) else {
            return Err(bad_format("missing version in JSON"));
        };
        if version.as_i64()? != CURRENT_MOL_JSON_VERSION {
            return Err(bad_format("bad version in JSON"));
        }
    } else if let Some(header) = obj.get("rdkitjson") {
        let Some(version) = header.as_object().ok().and_then(|h| h.get("version")) else {
            return Err(bad_format("missing version in JSON"));
        };
        // version 10 files can be read by 11 and 12, but not vice versa
        let version = version.as_i64()?;
        if !(10..=CURRENT_RDKIT_JSON_VERSION).contains(&version) {
            return Err(bad_format("bad version in JSON"));
        }
    } else {
        return Err(bad_format("missing header in JSON"));
    }
    Ok(())
}

/// The `defaults.<key>` object, which must be an object when present.
fn read_defaults<'a>(doc: &'a Json, key: &str) -> Result<Defaults<'a>> {
    if let Some(defaults) = doc.as_object()?.get("defaults") {
        if has(defaults, key)? {
            let d = defaults.at(key)?;
            if !d.is_object() {
                return Err(bad_format(format!("{key}Defaults is not an object")));
            }
            return Ok(Defaults(Some(d)));
        }
    }
    Ok(Defaults(None))
}

/// `JSONDataToMols`: the molecules are built on all cores.
pub fn json_to_mols(text: &str, params: &JsonParseParameters) -> Result<Vec<ROMol>> {
    let doc = Json::parse(text)?;
    if !doc.is_object() {
        return Err(bad_format("JSON should be an object"));
    }
    check_header(&doc)?;
    let atom_defaults = read_defaults(&doc, "atom")?;
    let bond_defaults = read_defaults(&doc, "bond")?;
    let Some(molecules) = doc.as_object()?.get("molecules") else {
        return Ok(Vec::new());
    };
    if !molecules.is_array() {
        return Err(bad_format("molecules is not an array"));
    }
    let built = par_map(molecules.as_array()?, 0, |val| {
        let mut mol = RWMol::new();
        process_mol(&mut mol, val, &atom_defaults, &bond_defaults, params)?;
        mol.update_property_cache(params.strict_valence_check)?;
        mol.set_prop(cp::_StereochemDone, 1, false)?;
        Ok(mol)
    });
    built.into_iter().collect()
}
