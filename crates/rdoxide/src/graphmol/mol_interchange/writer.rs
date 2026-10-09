//! MolJSON writer (upstream `Writer.cpp`).

use super::{
    CURRENT_CHARGE_REPRESENTATION_VERSION, CURRENT_MOL_JSON_VERSION,
    CURRENT_QUERY_REPRESENTATION_VERSION, CURRENT_RDKIT_JSON_VERSION,
    CURRENT_RDKIT_REPRESENTATION_VERSION, JsonWriteParameters, bond_order, bond_stereo_name,
    chiral_tag_name, stereo_group_type_name,
};
use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, symmetrize_sssr};
use crate::graphmol::kekulize::kekulize;
use crate::graphmol::pickle::{QueryDetails, query_details, tags};
use crate::graphmol::queries::{Query, QueryKind, QueryTarget};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::stereo_group::StereoGroup;
use crate::graphmol::substance_group::SubstanceGroup;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::json::{Json, JsonObject};
use crate::rdgeneral::threads::par_map;
use crate::rdgeneral::{Dict, RDKIT_VERSION, RDProps, RDValue};
use crate::{RDKitError, Result, SanitizeError};

/// the default of `MolOps::Kekulize`
const KEKULIZE_MAX_BACK_TRACKS: u32 = 100;

fn init_atom_defaults() -> JsonObject {
    let mut d = JsonObject::new();
    d.insert("z", 6);
    d.insert("impHs", 0);
    d.insert("chg", 0);
    d.insert("nRad", 0);
    d.insert("isotope", 0);
    d.insert("stereo", "unspecified");
    d
}

fn init_bond_defaults() -> JsonObject {
    let mut d = JsonObject::new();
    d.insert("bo", 1);
    d.insert("stereo", "unspecified");
    d
}

/// Writes `val` unless it equals the default.
fn add_int_val(dest: &mut JsonObject, defaults: &JsonObject, tag: &str, val: i64) -> Result<()> {
    match defaults.get(tag) {
        Some(default) if default.as_i64()? == val => {}
        _ => dest.insert(tag, val),
    }
    Ok(())
}

/// Writes `val` unless it is empty or equals the default.
fn add_string_val(
    dest: &mut JsonObject,
    defaults: &JsonObject,
    tag: &str,
    val: &str,
) -> Result<()> {
    match defaults.get(tag) {
        Some(default) => {
            if !val.is_empty() && default.as_str()? != val {
                dest.insert(tag, val);
            }
        }
        None => dest.insert(tag, val),
    }
    Ok(())
}

fn add_atom(mol: &ROMol, idx: u32, defaults: &JsonObject) -> Result<JsonObject> {
    let atom = mol.atom(idx);
    let mut out = JsonObject::new();
    add_int_val(&mut out, defaults, "z", i64::from(atom.atomic_num()))?;
    if !atom.has_query() {
        add_int_val(
            &mut out,
            defaults,
            "impHs",
            i64::from(mol.total_num_hs(idx, false)?),
        )?;
    }
    add_int_val(&mut out, defaults, "chg", i64::from(atom.formal_charge()))?;
    add_int_val(&mut out, defaults, "isotope", i64::from(atom.isotope()))?;
    add_int_val(
        &mut out,
        defaults,
        "nRad",
        i64::from(atom.num_radical_electrons()),
    )?;
    let stereo = chiral_tag_name(atom.chiral_tag()).unwrap_or_else(|| {
        crate::rd_log_warning!(" unrecognized atom chirality set to default while writing");
        ""
    });
    add_string_val(&mut out, defaults, "stereo", stereo)?;
    Ok(out)
}

fn add_bond(
    mol: &ROMol,
    idx: u32,
    defaults: &JsonObject,
    params: &JsonWriteParameters,
) -> Result<JsonObject> {
    let bond = mol.bond(idx);
    let mut out = JsonObject::new();
    let mut order = bond_order(bond.bond_type());
    // commonchem only supports a few bond orders
    if !params.use_rdkit_extensions && order.is_some_and(|o| o > 3) {
        order = None;
    }
    let order = order.unwrap_or_else(|| {
        crate::rd_log_warning!(
            " unrecognized bond type {} set to zero while writing",
            bond.bond_type() as u8
        );
        0
    });
    add_int_val(&mut out, defaults, "bo", order)?;
    out.insert(
        "atoms",
        vec![
            Json::from(bond.begin_atom_idx()),
            Json::from(bond.end_atom_idx()),
        ],
    );
    let stereo = bond_stereo_name(bond.stereo());
    add_string_val(&mut out, defaults, "stereo", stereo)?;
    if stereo != "unspecified" && bond.stereo_atoms().len() == 2 {
        out.insert(
            "stereoAtoms",
            bond.stereo_atoms()
                .iter()
                .map(|&a| Json::from(a))
                .collect::<Vec<_>>(),
        );
    }
    Ok(out)
}

/// `addProperties`: the named public properties, as ints, doubles and
/// strings; anything else is skipped with a warning.
fn add_properties(dict: &Dict, prop_names: &[String]) -> JsonObject {
    let mut props = JsonObject::new();
    for (key, val) in dict.iter() {
        if !prop_names.iter().any(|n| n == key) {
            continue;
        }
        match val {
            RDValue::Int(i) => props.insert(key, *i),
            RDValue::UnsignedInt(u) => props.insert(key, *u as i32),
            RDValue::Double(d) => props.insert(key, *d),
            RDValue::Float(f) => props.insert(key, f64::from(*f)),
            RDValue::String(s) => props.insert(key, s.as_str()),
            _ => crate::rd_log_warning!(
                "Warning: Could not convert property {key} to a recognized type. Skipping it."
            ),
        }
    }
    props
}

fn index_array<I: IntoIterator<Item = u32>>(indices: I) -> Json {
    Json::Array(indices.into_iter().map(Json::from).collect())
}

fn point_array(x: f64, y: f64, z: f64) -> Json {
    Json::Array(vec![Json::from(x), Json::from(y), Json::from(z)])
}

fn add_substance_group(sg: &SubstanceGroup) -> JsonObject {
    let mut out = JsonObject::new();
    let prop_names = sg.get_prop_list(false, false);
    if !prop_names.is_empty() {
        out.insert("properties", add_properties(sg.dict(), &prop_names));
    }
    if !sg.atoms().is_empty() {
        out.insert("atoms", index_array(sg.atoms().iter().copied()));
    }
    if !sg.bonds().is_empty() {
        out.insert("bonds", index_array(sg.bonds().iter().copied()));
    }
    if !sg.parent_atoms().is_empty() {
        out.insert(
            "parentAtoms",
            index_array(sg.parent_atoms().iter().copied()),
        );
    }
    if !sg.brackets().is_empty() {
        let brackets = sg
            .brackets()
            .iter()
            .map(|bracket| {
                Json::Array(bracket.iter().map(|p| point_array(p.x, p.y, p.z)).collect())
            })
            .collect::<Vec<_>>();
        out.insert("brackets", brackets);
    }
    if !sg.cstates().is_empty() {
        let is_sup = matches!(sg.get_prop("TYPE"), Some(RDValue::String(t)) if t == "SUP");
        let cstates = sg
            .cstates()
            .iter()
            .map(|cs| {
                let mut o = JsonObject::new();
                o.insert("bond", cs.bond_idx);
                if is_sup {
                    o.insert("vector", point_array(cs.vector.x, cs.vector.y, cs.vector.z));
                }
                Json::Object(o)
            })
            .collect::<Vec<_>>();
        out.insert("cstates", cstates);
    }
    if !sg.attach_points().is_empty() {
        let points = sg
            .attach_points()
            .iter()
            .map(|ap| {
                let mut o = JsonObject::new();
                o.insert("aIdx", ap.a_idx);
                if ap.lv_idx != -1 {
                    o.insert("lvIdx", ap.lv_idx);
                }
                if !ap.id.is_empty() {
                    o.insert("id", ap.id.as_str());
                }
                Json::Object(o)
            })
            .collect::<Vec<_>>();
        out.insert("attachPoints", points);
    }
    out
}

fn add_stereo_group(sg: &StereoGroup) -> JsonObject {
    let mut out = JsonObject::new();
    out.insert("type", stereo_group_type_name(sg.group_type()));
    if sg.group_type() != crate::graphmol::stereo_group::StereoGroupType::Absolute {
        out.insert("id", sg.write_id());
    }
    if !sg.atoms().is_empty() {
        out.insert("atoms", index_array(sg.atoms().iter().copied()));
    }
    if !sg.bonds().is_empty() {
        out.insert("bonds", index_array(sg.bonds().iter().copied()));
    }
    out
}

fn add_conformer(conf: &crate::graphmol::conformer::Conformer) -> JsonObject {
    let dim = if conf.is_3d() { 3 } else { 2 };
    let coords = conf
        .positions()
        .iter()
        .map(|p| {
            let mut pos = vec![Json::from(p.x), Json::from(p.y)];
            if dim == 3 {
                pos.push(Json::from(p.z));
            }
            Json::Array(pos)
        })
        .collect::<Vec<_>>();
    let mut out = JsonObject::new();
    out.insert("dim", dim);
    out.insert("coords", coords);
    out
}

fn add_query<F: QueryTarget>(q: &Query<F>, params: &JsonWriteParameters) -> Result<JsonObject> {
    let mut out = JsonObject::new();
    out.insert("descr", q.description.as_str());
    if !q.type_label.is_empty() {
        out.insert("type", q.type_label.as_str());
    }
    if q.negate {
        out.insert("negated", true);
    }
    if let QueryKind::Recursive { mol, .. } = &q.kind {
        let mut sub = JsonObject::new();
        add_mol(
            &mol.0,
            &mut sub,
            &init_atom_defaults(),
            &init_bond_defaults(),
            params,
        )?;
        out.insert("tag", i32::from(tags::QUERY_RECURSIVE));
        out.insert("subquery", sub);
    } else {
        match query_details(q) {
            Some(QueryDetails::Tag(tag)) => out.insert("tag", i32::from(tag)),
            Some(QueryDetails::ValTol { tag, val, tol }) => {
                out.insert("tag", i32::from(tag));
                out.insert("val", val);
                if tol != 0 {
                    out.insert("tol", tol);
                }
            }
            Some(QueryDetails::Range {
                lower,
                upper,
                tol,
                ends,
            }) => {
                out.insert("tag", i32::from(tags::QUERY_RANGE));
                out.insert("lower", lower);
                out.insert("upper", upper);
                if tol != 0 {
                    out.insert("tol", tol);
                }
                out.insert("ends", i32::from(ends));
            }
            Some(QueryDetails::Set(set)) => {
                out.insert("tag", i32::from(tags::QUERY_SET));
                out.insert("set", set.into_iter().map(Json::from).collect::<Vec<_>>());
            }
            // the property queries have no JSON form
            _ => {
                return Err(RDKitError::Runtime(
                    "do not know how to pickle part of the query.".into(),
                ));
            }
        }
    }
    if !q.children.is_empty() {
        let children = q
            .children
            .iter()
            .map(|c| add_query(c, params).map(Json::Object))
            .collect::<Result<Vec<_>>>()?;
        out.insert("children", children);
    }
    Ok(out)
}

/// The molecule with symmetrized SSSR and a property cache, kekulized when
/// possible (`RWMol mol(imol)` and the following steps of `addMol`).
///
/// When kekulization fails the input molecule is restored, so an unsanitized
/// molecule ends up without ring info and writing it fails, as upstream.
fn prepare(imol: &ROMol) -> Result<RWMol> {
    let mut mol = imol.clone();
    if !mol.ring_info().is_symm_sssr() {
        symmetrize_sssr(
            &mut mol,
            SymmetrizeSSSRAlgorithm::Default,
            true,
            false,
            false,
        )?;
    }
    if mol.needs_update_property_cache() {
        mol.update_property_cache(false)?;
    }
    match kekulize(&mut mol, false, false, KEKULIZE_MAX_BACK_TRACKS) {
        Ok(()) => {}
        Err(RDKitError::Sanitize(
            SanitizeError::Kekulize { .. } | SanitizeError::AtomKekulize { .. },
        )) => {
            mol = imol.clone();
            if mol.needs_update_property_cache() {
                mol.update_property_cache(false)?;
            }
        }
        Err(e) => return Err(e),
    }
    Ok(mol)
}

fn rdkit_representation(mol: &ROMol) -> Result<JsonObject> {
    let mut rep = JsonObject::new();
    rep.insert("name", "rdkitRepresentation");
    rep.insert("formatVersion", CURRENT_RDKIT_REPRESENTATION_VERSION);
    rep.insert("toolkitVersion", RDKIT_VERSION);
    if mol.atoms().iter().any(|a| a.is_aromatic()) {
        let atoms = mol
            .atoms()
            .iter()
            .filter(|a| a.is_aromatic())
            .map(|a| a.idx());
        rep.insert("aromaticAtoms", index_array(atoms));
        let bonds = mol
            .bonds()
            .iter()
            .filter(|b| b.is_aromatic())
            .map(|b| b.idx());
        rep.insert("aromaticBonds", index_array(bonds));
    }
    if mol.get_atom(0)?.has_prop(cp::_CIPRank) {
        let ranks = mol
            .atoms()
            .iter()
            .map(|a| {
                a.get_prop(cp::_CIPRank)
                    .and_then(RDValue::as_uint)
                    .map(Json::from)
                    .ok_or_else(|| RDKitError::Key(cp::_CIPRank.to_string()))
            })
            .collect::<Result<Vec<_>>>()?;
        rep.insert("cipRanks", ranks);
    }
    let codes = mol
        .atoms()
        .iter()
        .filter_map(|a| {
            a.get_prop(cp::_CIPCode)
                .map(|code| Json::Array(vec![Json::from(a.idx()), Json::from(code.to_rd_string())]))
        })
        .collect::<Vec<_>>();
    if !codes.is_empty() {
        rep.insert("cipCodes", codes);
    }
    if mol.ring_info().num_rings()? > 0 {
        let rings = mol
            .ring_info()
            .atom_rings()
            .iter()
            .map(|ring| Json::Array(ring.iter().map(|&a| Json::from(a)).collect()))
            .collect::<Vec<_>>();
        rep.insert("atomRings", rings);
    }
    Ok(rep)
}

fn partial_charges(mol: &ROMol) -> JsonObject {
    let mut rep = JsonObject::new();
    rep.insert("name", "partialCharges");
    rep.insert("generator", "RDKit");
    rep.insert("formatVersion", CURRENT_CHARGE_REPRESENTATION_VERSION);
    rep.insert("generatorVersion", RDKIT_VERSION);
    let values = mol
        .atoms()
        .iter()
        .map(|a| {
            Json::from(
                a.get_prop(cp::_GasteigerCharge)
                    .and_then(RDValue::as_double)
                    .unwrap_or(0.0),
            )
        })
        .collect::<Vec<_>>();
    rep.insert("values", values);
    rep
}

/// One entry per atom/bond: the query, or an empty object.
fn query_array<'a, F: QueryTarget + 'a>(
    queries: impl Iterator<Item = Option<&'a Query<F>>>,
    params: &JsonWriteParameters,
) -> Result<Json> {
    queries
        .map(|q| match q {
            Some(q) => add_query(q, params).map(Json::Object),
            None => Ok(Json::Object(JsonObject::new())),
        })
        .collect::<Result<Vec<_>>>()
        .map(Json::Array)
}

fn query_representation(
    mol: &ROMol,
    has_query_atoms: bool,
    has_query_bonds: bool,
    params: &JsonWriteParameters,
) -> Result<JsonObject> {
    let mut rep = JsonObject::new();
    rep.insert("name", "rdkitQueries");
    rep.insert("formatVersion", CURRENT_QUERY_REPRESENTATION_VERSION);
    rep.insert("toolkitVersion", RDKIT_VERSION);
    if has_query_atoms {
        rep.insert(
            "atomQueries",
            query_array(mol.atoms().iter().map(|a| a.query()), params)?,
        );
    }
    if has_query_bonds {
        rep.insert(
            "bondQueries",
            query_array(mol.bonds().iter().map(|b| b.query()), params)?,
        );
    }
    Ok(rep)
}

/// `addMol`: fills `out` with the description of `imol`.
fn add_mol(
    imol: &ROMol,
    out: &mut JsonObject,
    atom_defaults: &JsonObject,
    bond_defaults: &JsonObject,
    params: &JsonWriteParameters,
) -> Result<()> {
    let mol = prepare(imol)?;
    if let Some(name) = mol.get_prop(cp::_Name) {
        out.insert("name", name.to_rd_string());
    }
    let atoms = (0..mol.num_atoms())
        .map(|i| add_atom(&mol, i, atom_defaults).map(Json::Object))
        .collect::<Result<Vec<_>>>()?;
    out.insert("atoms", atoms);
    let bonds = (0..mol.num_bonds())
        .map(|i| add_bond(&mol, i, bond_defaults, params).map(Json::Object))
        .collect::<Result<Vec<_>>>()?;
    out.insert("bonds", bonds);

    if params.use_rdkit_extensions && !mol.stereo_groups().is_empty() {
        let groups = mol
            .stereo_groups()
            .iter()
            .map(|sg| Json::Object(add_stereo_group(sg)))
            .collect::<Vec<_>>();
        out.insert("stereoGroups", groups);
    }
    if params.use_rdkit_extensions && !mol.substance_groups().is_empty() {
        let groups = mol
            .substance_groups()
            .iter()
            .map(|sg| Json::Object(add_substance_group(sg)))
            .collect::<Vec<_>>();
        out.insert("substanceGroups", groups);
    }
    if mol.num_conformers() > 0 {
        let confs = mol
            .conformers()
            .iter()
            .map(|c| Json::Object(add_conformer(c)))
            .collect::<Vec<_>>();
        out.insert("conformers", confs);
    }
    let prop_names = mol.get_prop_list(false, false);
    if !prop_names.is_empty() {
        out.insert("properties", add_properties(mol.dict(), &prop_names));
    }

    let mut extensions = vec![Json::Object(rdkit_representation(&mol)?)];
    if mol.get_atom(0)?.has_prop(cp::_GasteigerCharge) {
        extensions.push(Json::Object(partial_charges(&mol)));
    }
    let has_query_atoms = mol.atoms().iter().any(|a| a.has_query());
    let has_query_bonds = mol.bonds().iter().any(|b| b.has_query());
    if has_query_atoms || has_query_bonds {
        extensions.push(Json::Object(query_representation(
            &mol,
            has_query_atoms,
            has_query_bonds,
            params,
        )?));
    }
    out.insert("extensions", extensions);
    Ok(())
}

/// `MolsToJSONData`: the molecules are described on all cores.
pub fn mols_to_json(mols: &[Option<&ROMol>], params: &JsonWriteParameters) -> Result<String> {
    let mut doc = JsonObject::new();
    let mut header = JsonObject::new();
    let (key, version) = if params.use_rdkit_extensions {
        ("rdkitjson", CURRENT_RDKIT_JSON_VERSION)
    } else {
        ("commonchem", CURRENT_MOL_JSON_VERSION)
    };
    header.insert("version", version);
    doc.insert(key, header);

    let atom_defaults = init_atom_defaults();
    let bond_defaults = init_bond_defaults();
    let mut defaults = JsonObject::new();
    defaults.insert("atom", atom_defaults.clone());
    defaults.insert("bond", bond_defaults.clone());
    doc.insert("defaults", defaults);

    if mols.iter().any(Option::is_none) {
        return Err(RDKitError::Value(
            "null molecule passed to MolsToJSONData".into(),
        ));
    }
    let described = par_map(mols, 0, |mol| {
        let mut out = JsonObject::new();
        let mol = mol.expect("null molecules were rejected above");
        add_mol(mol, &mut out, &atom_defaults, &bond_defaults, params).map(|()| Json::Object(out))
    });
    doc.insert(
        "molecules",
        described.into_iter().collect::<Result<Vec<_>>>()?,
    );
    Ok(Json::Object(doc).serialize())
}

/// `MolToJSONData`
pub fn mol_to_json(mol: &ROMol, params: &JsonWriteParameters) -> Result<String> {
    mols_to_json(&[Some(mol)], params)
}
