//! Port of RDKit `MarvinParser.cpp`'s reader (`MolFromMrvBlock`,
//! `MolFromMrvFile`): Marvin CML to an `RWMol`.

use super::defs::{MarvinAtom, MarvinMol, NO_COORD, NodeId, ROOT, SgroupKind, radical_electrons};
use super::xml::{XmlTree, read_xml};
use crate::geometry::Point3D;
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::{Bond, BondDir, BondStereo, BondType};
use crate::graphmol::chirality::clear_single_bond_dir_flags;
use crate::graphmol::chirality::from_coords::{
    assign_chiral_types_from_3d, assign_chiral_types_from_bond_dirs, detect_bond_stereochemistry,
};
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, remove_hs};
use crate::graphmol::mol_ops::sanitize_flags::{SANITIZE_ALL, SANITIZE_CLEANUP};
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::queries::{
    complete_mol_queries, make_atom_null_query, make_bond_null_query,
    make_double_or_aromatic_bond_query, make_single_or_aromatic_bond_query,
    make_single_or_double_bond_query,
};
use crate::graphmol::romol::RWMol;
use crate::graphmol::stereo_group::{StereoGroup, StereoGroupType};
use crate::graphmol::substance_group::{SubstanceGroup, add_substance_group};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::{RDKitError, Result};

/// `MrvParserParams`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MrvParserParams {
    pub sanitize: bool,
    pub remove_hs: bool,
}

impl Default for MrvParserParams {
    fn default() -> Self {
        MrvParserParams {
            sanitize: true,
            remove_hs: true,
        }
    }
}

fn file_parse(msg: impl Into<String>) -> RDKitError {
    RDKitError::FileParse(msg.into())
}

/// `MarvinCMLReader::molAtomFromMarvinAtom`
fn mol_atom_from_marvin_atom(marvin: &MarvinMol, node: NodeId, atom: &MarvinAtom) -> Result<Atom> {
    let symb = atom.element_type.trim().to_string();
    let mut res;
    if symb == "*" {
        res = Atom::new(0);
        res.set_query(Some(make_atom_null_query()));
        res.props.set_val(cp::dummyLabel, "*".to_string());
        // queries have no implicit Hs
        res.set_no_implicit(true);
    } else if symb == "R" {
        if atom.rgroup_ref < 0 {
            return Err(file_parse(
                "Expected an R group atom to have an rgroupRef designation",
            ));
        }
        res = Atom::new(0);
        res.props
            .set_val(cp::dummyLabel, format!("R{}", atom.rgroup_ref));
        if !atom.mrv_alias.is_empty() {
            res.props.set_val(cp::molFileAlias, atom.mrv_alias.clone());
        }
        res.set_query(Some(make_atom_null_query()));
        res.set_isotope(atom.rgroup_ref as u32);
        res.props
            .set_val(cp::_MolFileRLabel, atom.rgroup_ref as u32);
    } else if symb.chars().count() <= 2 {
        let mut s: Vec<char> = symb.chars().collect();
        if s.len() == 2 && s[1].is_ascii_uppercase() {
            s[1] = s[1].to_ascii_lowercase();
        }
        let s: String = s.into_iter().collect();
        res = Atom::new(0);
        let num = PeriodicTable::get()
            .get_atomic_number(&s)
            .map_err(|e| file_parse(e.to_string()))?;
        res.set_atomic_num(num as i32);
    } else {
        return Err(file_parse("Unrecognized atom type in MRV file"));
    }

    if atom.formal_charge != 0 {
        res.set_formal_charge(atom.formal_charge);
    }
    if atom.isotope != 0 {
        res.set_isotope(atom.isotope as u32);
    }
    // mrvValence and hydrogenCount both set the number of explicit Hs; if
    // both are present they must agree
    if atom.hydrogen_count >= 0 {
        res.set_num_explicit_hs(atom.hydrogen_count as u32);
    }
    if atom.mrv_valence >= 0 {
        let explicit = marvin.explicit_valence(node, atom)?;
        let h_count = atom.mrv_valence - explicit;
        if h_count >= 0 {
            if atom.hydrogen_count >= 0 {
                if atom.hydrogen_count != h_count {
                    return Err(file_parse(
                        "mrvValence and hydrogenCount both specified for an atom, and they do not agree",
                    ));
                }
            } else {
                res.set_no_implicit(true);
                res.set_num_explicit_hs(h_count as u32);
            }
        }
    }
    if !atom.radical.is_empty() {
        let n = radical_electrons(&atom.radical)
            .ok_or_else(|| RDKitError::Runtime("invalid map<K, T> key".into()))?;
        res.set_num_radical_electrons(n);
    }
    if atom.mrv_map != 0 {
        res.props.set_val(cp::molAtomMapNumber, atom.mrv_map);
    }
    Ok(res)
}

/// `MarvinCMLReader::molBondFromMarvinBond`
fn mol_bond_from_marvin_bond(
    marvin: &MarvinMol,
    bond: &super::defs::MarvinBond,
    mol: &mut RWMol,
    chirality_possible: &mut bool,
) -> Result<()> {
    let idx1 = marvin.atom_index(ROOT, &bond.atom_refs2[0]);
    let idx2 = marvin.atom_index(ROOT, &bond.atom_refs2[1]);

    let mut b_type: u32;
    let mut btype;
    let mut res;
    match bond.bond_type()?.as_str() {
        "SD" => {
            b_type = 5;
            btype = BondType::Unspecified;
            res = Bond::new(btype);
            res.set_query(Some(make_single_or_double_bond_query()));
        }
        "SA" => {
            b_type = 6;
            btype = BondType::Unspecified;
            res = Bond::new(btype);
            res.set_query(Some(make_single_or_aromatic_bond_query()));
        }
        "DA" => {
            b_type = 7;
            btype = BondType::Unspecified;
            res = Bond::new(btype);
            res.set_query(Some(make_double_or_aromatic_bond_query()));
        }
        "ANY" => {
            b_type = 8;
            btype = BondType::Unspecified;
            res = Bond::new(btype);
            res.set_query(Some(make_bond_null_query()));
        }
        "DATIVE" => {
            b_type = 1;
            btype = BondType::Dative;
            res = Bond::new(btype);
        }
        "HYDROGEN" => {
            b_type = 1;
            btype = BondType::Hydrogen;
            res = Bond::new(btype);
        }
        "1" => {
            b_type = 1;
            btype = BondType::Single;
            res = Bond::new(btype);
        }
        "2" => {
            b_type = 2;
            btype = BondType::Double;
            res = Bond::new(btype);
        }
        "3" => {
            b_type = 3;
            btype = BondType::Triple;
            res = Bond::new(btype);
        }
        "A" => {
            b_type = 4;
            btype = BondType::Aromatic;
            res = Bond::new(btype);
        }
        _ => return Err(RDKitError::invariant("unexpected marvin bond type")),
    }
    res.set_begin_atom_idx_unchecked(idx1 as u32);
    res.set_end_atom_idx_unchecked(idx2 as u32);

    let mut bond_dir = BondDir::None;
    let mut bond_stereo = BondStereo::StereoNone;
    let mut cfg = 0u32;
    let stereo = &bond.bond_stereo;
    let value = stereo.value.to_lowercase();
    if !value.is_empty() {
        match value.as_str() {
            "w" => bond_dir = BondDir::BeginWedge,
            "h" => bond_dir = BondDir::BeginDash,
            _ => {
                return Err(file_parse(format!(
                    "unrecognized bond stereo value{}",
                    stereo.value
                )));
            }
        }
    } else if !stereo.convention.is_empty() {
        if stereo.convention != "MDL" {
            return Err(file_parse(format!(
                "unrecognized bond stereo conventrion{}",
                stereo.convention
            )));
        }
        let mdl: i32 = super::defs::clean_number(&stereo.convention_value)
            .ok_or_else(|| file_parse("MDL Convention Value must be one of: 1, 3, 4, 6"))?;
        bond_dir = match mdl {
            1 => BondDir::BeginWedge,
            6 => BondDir::BeginDash,
            3 => BondDir::EitherDouble,
            4 => BondDir::Unknown,
            _ => {
                return Err(file_parse(
                    "MDL Convention Value must be one of: 1, 3, 4, 6",
                ));
            }
        };
    } else if !stereo.dict_ref.is_empty() {
        bond_dir = match stereo.dict_ref.as_str() {
            "cml:W" => BondDir::BeginWedge,
            "cml:H" => BondDir::BeginDash,
            _ => return Err(file_parse("dictRef must be one of: cml:W or cml:H")),
        };
    }

    match bond_dir {
        BondDir::BeginWedge => {
            b_type = 1;
            btype = BondType::Single;
            *chirality_possible = true;
            cfg = 1;
        }
        BondDir::BeginDash => {
            b_type = 1;
            btype = BondType::Single;
            cfg = 3;
            *chirality_possible = true;
        }
        BondDir::EitherDouble => {
            bond_stereo = BondStereo::StereoAny;
            b_type = 2;
            btype = BondType::Double;
            cfg = 2;
        }
        BondDir::Unknown => {
            b_type = 1;
            btype = BondType::Single;
            cfg = 2;
        }
        _ => {}
    }
    res.set_bond_type(btype);
    res.props.set_val(cp::_MolFileBondType, b_type);
    res.set_bond_dir(bond_dir);
    if cfg != 0 {
        res.props.set_val(cp::_MolFileBondCfg, cfg);
    }
    if bond_stereo != BondStereo::StereoNone {
        res.set_stereo(bond_stereo)?;
    }
    if res.bond_type() == BondType::Aromatic {
        res.set_is_aromatic(true);
    }
    mol.add_bond_object(res)?;
    Ok(())
}

fn conformer_position(a: &MarvinAtom, three_d: bool) -> Point3D {
    if three_d {
        if a.x3 != NO_COORD && a.y3 != NO_COORD && a.z3 != NO_COORD {
            Point3D::new(a.x3, a.y3, a.z3)
        } else {
            Point3D::new(0.0, 0.0, 0.0)
        }
    } else if a.x2 != NO_COORD && a.y2 != NO_COORD {
        Point3D::new(a.x2, a.y2, 0.0)
    } else {
        Point3D::new(0.0, 0.0, 0.0)
    }
}

/// The substance group a child of the root becomes
/// (`parseMoleculeSpecific`).
fn substance_group_of(
    marvin: &MarvinMol,
    mol: &RWMol,
    node: NodeId,
    sequence_id: u32,
) -> Result<SubstanceGroup> {
    let n = marvin.node(node);
    let index_of_atom = |a: usize| marvin.atom_index(ROOT, &marvin.atoms[a].id) as u32;
    let index_of_bond = |b: usize| marvin.bond_index(ROOT, &marvin.bonds[b].id) as u32;
    let add_atoms = |sg: &mut SubstanceGroup| -> Result<()> {
        for &a in &n.atoms {
            sg.add_atom_with_idx(mol, index_of_atom(a))?;
        }
        Ok(())
    };
    let add_bonds = |sg: &mut SubstanceGroup| -> Result<()> {
        for &b in &n.bonds {
            sg.add_bond_with_idx(mol, index_of_bond(b))?;
        }
        Ok(())
    };
    let mut sg;
    match &n.kind {
        SgroupKind::SuperatomExpanded { title } => {
            sg = SubstanceGroup::new("SUP");
            sg.set_prop("index", sequence_id, false)?;
            add_atoms(&mut sg)?;
            add_bonds(&mut sg)?;
            sg.set_prop("LABEL", title.clone(), false)?;
        }
        SgroupKind::Multiple {
            title,
            parent_atoms,
            bonds_to_atoms_not_in_expanded_group,
            ..
        } => {
            sg = SubstanceGroup::new("MUL");
            sg.set_prop("index", sequence_id, false)?;
            add_atoms(&mut sg)?;
            for &a in parent_atoms {
                sg.add_parent_atom_with_idx(mol, index_of_atom(a))?;
            }
            for &b in bonds_to_atoms_not_in_expanded_group {
                sg.add_bond_with_idx(mol, index_of_bond(b))?;
            }
            sg.set_prop("MULT", title.clone(), false)?;
        }
        SgroupKind::SruCoMod {
            role_name,
            title,
            connect,
            ..
        } => {
            let typ = match role_name.as_str() {
                "SruSgroup" => "SRU",
                "CopolymerSgroup" => "COP",
                "ModificationSgroup" => "MOD",
                _ => {
                    return Err(file_parse(
                        "Internal error: unrecognized role in a MarvinSruCoPolSgroup",
                    ));
                }
            };
            sg = SubstanceGroup::new(typ);
            sg.set_prop("index", sequence_id, false)?;
            sg.set_prop("CONNECT", connect.clone(), false)?;
            add_atoms(&mut sg)?;
            add_bonds(&mut sg)?;
            sg.set_prop("LABEL", title.clone(), false)?;
        }
        SgroupKind::Data {
            context,
            field_name,
            placement,
            units_displayed,
            query_type,
            query_op,
            field_data,
            units,
            x,
            y,
        } => {
            sg = SubstanceGroup::new("DAT");
            sg.set_prop("index", sequence_id, false)?;
            add_atoms(&mut sg)?;
            sg.set_prop("FIELDNAME", field_name.clone(), false)?;
            if !query_type.is_empty() {
                sg.set_prop("QUERYTYPE", query_type.clone(), false)?;
            }
            if !query_op.is_empty() {
                sg.set_prop("QUERYOP", query_op.clone(), false)?;
            }
            // (not used by RDKit)
            sg.set_prop(
                "FIELDDISP",
                format!("{x:10.4}{y:10.4}    DRU   ALL  0       0"),
                false,
            )?;
            sg.set_prop("DATAFIELDS", vec![field_data.clone()], false)?;
            // not part of the MOL format, kept to round-trip MRV
            sg.set_prop("UNITS", units.clone(), false)?;
            sg.set_prop("UNITSDISPLAYED", units_displayed.clone(), false)?;
            sg.set_prop("CONTEXT", context.clone(), false)?;
            sg.set_prop("PLACEMENT", placement.clone(), false)?;
            sg.set_prop("X", *x, false)?;
            sg.set_prop("Y", *y, false)?;
        }
        SgroupKind::Multicenter { .. } => {
            return Err(RDKitError::precondition(
                "Internal error:   a MarvinMulticenterSgroup has not been removed",
            ));
        }
        SgroupKind::Generic { .. } => {
            sg = SubstanceGroup::new("GEN");
            sg.set_prop("index", sequence_id, false)?;
            add_atoms(&mut sg)?;
        }
        SgroupKind::Monomer { title, .. } => {
            sg = SubstanceGroup::new("MON");
            sg.set_prop("index", sequence_id, false)?;
            // (sic) the index within the group of every atom of the parent
            let parent = n.parent.expect("a group has a parent");
            for &a in &marvin.node(parent).atoms {
                let idx = marvin.atom_index(node, &marvin.atoms[a].id);
                sg.add_atom_with_idx(mol, idx as u32)?;
            }
            sg.set_prop("LABEL", title.clone(), false)?;
        }
        SgroupKind::Mol | SgroupKind::Superatom { .. } => {
            return Err(RDKitError::invariant("unexpected kind of sgroup"));
        }
    }
    Ok(sg)
}

/// `MarvinCMLReader::parseMolecule(MarvinMol *, ...)`
pub fn mol_from_marvin_mol(marvin: &MarvinMol, params: &MrvParserParams) -> Result<RWMol> {
    let mut mol = RWMol::new();
    mol.props
        .set_val("_MolFileComments", "Generated by RDKit".to_string());
    let atoms = &marvin.node(ROOT).atoms;

    let mut conf2d = marvin
        .has_any_2d_coords(ROOT)
        .then(|| Conformer::new(atoms.len()));
    if let Some(c) = &mut conf2d {
        c.set_3d(false);
    }
    let mut conf3d = marvin
        .has_any_3d_coords(ROOT)
        .then(|| Conformer::new(atoms.len()));
    if let Some(c) = &mut conf3d {
        c.set_3d(true);
    }

    // (type, number, atoms)
    let mut stereo_groups: Vec<(StereoGroupType, i32, Vec<u32>)> = Vec::new();
    for &a in atoms {
        let atom = &marvin.atoms[a];
        let aid = mol.add_atom(mol_atom_from_marvin_atom(marvin, ROOT, atom)?, false);
        if let Some(c) = &mut conf2d {
            c.set_atom_pos(aid, conformer_position(atom, false))?;
        }
        if let Some(c) = &mut conf3d {
            c.set_atom_pos(aid, conformer_position(atom, true))?;
        }
        if !atom.mrv_stereo_group.is_empty() {
            let temp = atom.mrv_stereo_group.to_lowercase();
            let (group_type, number) = if temp == "abs" {
                (StereoGroupType::Absolute, -1)
            } else if let Some(rest) = temp.strip_prefix("and") {
                let n = super::defs::clean_number(rest).ok_or_else(|| {
                    file_parse(
                        "Group Number must be an integer in a stereo group AND# in a MRV file",
                    )
                })?;
                (StereoGroupType::And, n)
            } else if let Some(rest) = temp.strip_prefix("or") {
                let n = super::defs::clean_number(rest).ok_or_else(|| {
                    file_parse(
                        "Group Number must be an integer in a stereo group OR# in a MRV file",
                    )
                })?;
                (StereoGroupType::Or, n)
            } else {
                return Err(file_parse("Unrecognized group definition"));
            };
            let index = marvin.atom_index(ROOT, &atom.id) as u32;
            match stereo_groups
                .iter_mut()
                .find(|(t, n, _)| *n == number && *t == group_type)
            {
                Some((_, _, list)) => list.push(index),
                None => stereo_groups.push((group_type, number, vec![index])),
            }
        }
    }
    let conf2d_id = match conf2d {
        Some(c) => Some(mol.add_conformer(c, true)?),
        None => None,
    };
    let conf3d_id = match conf3d {
        Some(c) => Some(mol.add_conformer(c, true)?),
        None => None,
    };

    let mut chirality_possible = false;
    for &b in &marvin.node(ROOT).bonds {
        mol_bond_from_marvin_bond(marvin, &marvin.bonds[b], &mut mol, &mut chirality_possible)?;
    }

    let mut groups = Vec::new();
    for (group_type, number, list) in stereo_groups {
        groups.push(StereoGroup::new(
            group_type,
            list,
            Vec::new(),
            number as u32,
        )?);
    }
    if !groups.is_empty() {
        mol.set_stereo_groups(groups)?;
    }

    // only the groups directly below the molecule are written
    for (sequence_id, &child) in marvin.node(ROOT).children.iter().enumerate() {
        let sg = substance_group_of(marvin, &mol, child, sequence_id as u32)?;
        if sg.is_valid() {
            add_substance_group(&mut mol, sg);
        }
    }

    mol.clear_all_atom_bookmarks();
    mol.clear_all_bond_bookmarks();
    for i in 0..mol.num_atoms() {
        mol.calc_explicit_valence(i, false)?;
        mol.calc_implicit_valence(i, false)?;
    }

    // NOTE: the chirality is detected before sanitizing/removing Hs because
    // removing them may remove the wedged bond
    if let (true, Some(id)) = (chirality_possible, conf2d_id) {
        assign_chiral_types_from_bond_dirs(&mut mol, id as i32, true)?;
    } else if let Some(id) = conf3d_id {
        mol.update_property_cache(false)?;
        assign_chiral_types_from_3d(&mut mol, id as i32, true)?;
    }
    clear_single_bond_dir_flags(&mut mol, false);

    if params.sanitize {
        if params.remove_hs {
            sanitize_mol(&mut mol, SANITIZE_CLEANUP).map_err(|(_, e)| e)?;
            detect_bond_stereochemistry(&mut mol, -1)?;
            remove_hs(&mut mol, &RemoveHsParameters::default(), true)?;
        } else {
            sanitize_mol(&mut mol, SANITIZE_ALL).map_err(|(_, e)| e)?;
            detect_bond_stereochemistry(&mut mol, -1)?;
        }
        assign_stereochemistry(&mut mol, true, true, true)?;
    } else {
        detect_bond_stereochemistry(&mut mol, -1)?;
    }

    if mol.has_prop(cp::_NeedsQueryScan) {
        mol.clear_prop(cp::_NeedsQueryScan);
        complete_mol_queries(&mut mol, 0xDEADBEEF)?;
    }
    Ok(mol)
}

/// `parseMolecule(ptree, ...)`: finds the molecule in the document.
fn parse_document(tree: &XmlTree, params: &MrvParserParams) -> Result<RWMol> {
    let Some(section) = tree.get_child("cml.MDocument.MChemicalStruct.molecule") else {
        if tree.get_child("cml.MDocument").is_some() || tree.get_child("cml").is_some() {
            return Ok(RWMol::new());
        }
        return Err(file_parse("Expected \"molecule\" in MRV file"));
    };
    let mut marvin = MarvinMol::parse(section)?;
    marvin.prep_sgroups_for_rdkit(ROOT)?;
    mol_from_marvin_mol(&marvin, params)
}

/// `MolFromMrvBlock`
pub fn mol_from_mrv_block(text: &str, params: &MrvParserParams) -> Result<RWMol> {
    parse_document(&read_xml(text)?, params)
}

/// `MolFromMrvFile`: `Ok(None)` for an empty file.
pub fn mol_from_mrv_file(path: &str, params: &MrvParserParams) -> Result<Option<RWMol>> {
    let bytes =
        std::fs::read(path).map_err(|_| RDKitError::BadFile(format!("Bad input file {path}")))?;
    if bytes.is_empty() {
        return Ok(None);
    }
    mol_from_mrv_block(&String::from_utf8_lossy(&bytes), params).map(Some)
}
