//! Port of RDKit `MarvinWriter.cpp`'s molecule writer (`MolToMrvBlock`,
//! `MolToMrvFile`): an `ROMol` to Marvin CML.

use super::defs::{
    MarvinAtom, MarvinBond, MarvinMol, MarvinNode, NO_COORD, NodeId, ROOT, SgroupKind,
};
use super::transform::NumberingCounts;
use super::xml::write_xml;
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::{Bond, BondDir, BondType};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::file_parsers::mol_file_writer::{
    all_atom_ids_for_stereo_group, get_mol_file_bond_stereo_info,
};
use crate::graphmol::kekulize::kekulize;
use crate::graphmol::queries::{AtomQuery, BondQuery};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::write::in_organic_subset;
use crate::graphmol::stereo_group::StereoGroupType;
use crate::graphmol::substance_group::SubstanceGroup;
use crate::graphmol::wedge_bonds::pick_bonds_to_wedge;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// `MrvWriterParams`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MrvWriterParams {
    pub include_stereo: bool,
    pub kekulize: bool,
    pub pretty_print: bool,
    pub precision: u32,
}

impl Default for MrvWriterParams {
    fn default() -> Self {
        MrvWriterParams {
            include_stereo: true,
            kekulize: true,
            pretty_print: false,
            precision: 6,
        }
    }
}

fn writer_error(msg: &str) -> RDKitError {
    RDKitError::Runtime(msg.to_string())
}

fn no_queries() -> RDKitError {
    writer_error("Query atoms are not supported for MarvinWriter")
}

fn string_prop(props: &impl RDProps, key: &str) -> Option<String> {
    match props.get_prop(key) {
        Some(RDValue::String(s)) => Some(s.clone()),
        _ => None,
    }
}

/// `hasComplexQuery`
fn has_complex_query(atom: &Atom) -> bool {
    let Some(q) = atom.query() else {
        return false;
    };
    match q.description.as_str() {
        "AtomAtomicNum" => false,
        "AtomAnd" => q
            .children
            .first()
            .is_none_or(|c| c.description != "AtomAtomicNum"),
        _ => true,
    }
}

fn is_x_x_query(q: &AtomQuery, label: &str) -> bool {
    q.type_label == label
}

/// `GetMarvinAtomInfo`
fn marvin_atom_info(atom: &Atom, out: &mut MarvinAtom) -> Result<()> {
    if let Some(r) = atom
        .props
        .get(cp::_MolFileRLabel)
        .and_then(RDValue::as_uint)
    {
        out.element_type = "R".into();
        out.rgroup_ref = r as i32;
        if let Some(alias) = string_prop(&atom.props, cp::molFileAlias) {
            out.mrv_alias = alias;
        }
    } else if atom.atomic_num() != 0 {
        out.element_type = atom.symbol();
    } else if let Some(symb) = string_prop(&atom.props, cp::dummyLabel) {
        let (group, alias) = match symb.as_str() {
            "*" => {
                out.element_type = "*".into();
                return Ok(());
            }
            "X" => (0, "R"),
            "Xa" => (1, "R1"),
            "Xb" => (2, "R2"),
            "Xc" => (3, "R3"),
            "Xd" => (4, "R4"),
            "Xf" => (5, "R5"),
            "Xg" => (6, "R6"),
            "Xh" => (7, "R7"),
            "Xi" => (8, "R8"),
            "Xj" => (9, "R9"),
            _ => return Err(no_queries()),
        };
        out.element_type = "R".into();
        out.rgroup_ref = group;
        out.mrv_alias = alias.into();
    } else if let Some(q) = atom.query() {
        let a_like = q.type_label == "A"
            || q.type_label.is_empty()
            || (q.negate && q.description == "AtomAtomicNum" && q.val == 1);
        if a_like {
            out.element_type = "*".into();
        } else {
            let q_like = is_x_x_query(q, "Q")
                || (q.negate
                    && q.description == "AtomOr"
                    && q.children.len() == 2
                    && q.children[0].description == "AtomAtomicNum"
                    && q.children[0].val == 6
                    && q.children[1].description == "AtomAtomicNum"
                    && q.children[1].val == 1);
            if q_like
                || ["X", "M", "AH", "QH", "XH", "MH"]
                    .iter()
                    .any(|l| is_x_x_query(q, l))
                || has_complex_query(atom)
            {
                return Err(no_queries());
            }
            out.element_type = "*".into();
        }
    } else {
        out.element_type = "*".into();
    }
    Ok(())
}

fn bond_in_ring(q: &BondQuery) -> bool {
    let mut q = q;
    if q.description == "BondAnd" && !q.negate && q.children.len() == 2 {
        let (c1, c2) = (&q.children[0], &q.children[1]);
        if (c1.description == "BondInRing") != (c2.description == "BondInRing") {
            q = if c1.description != "BondInRing" {
                c2
            } else {
                c1
            };
        }
    }
    q.description == "BondInRing"
}

/// `getMarvinQueryBondSymbol`
fn query_bond_symbol(q: &BondQuery) -> Result<&'static str> {
    if q.description == "BondOrder" || bond_in_ring(q) {
        return Ok("");
    }
    let mut q = q;
    if q.description == "BondAnd" && !q.negate && q.children.len() == 2 {
        let (c1, c2) = (&q.children[0], &q.children[1]);
        if c2.description == "BondInRing" {
            q = c1;
        } else if c1.description == "BondInRing" {
            q = c2;
        }
    }
    if q.description == "BondOr" && !q.negate {
        if q.children.len() == 2 {
            let (c1, c2) = (&q.children[0], &q.children[1]);
            if c1.description == "BondOrder"
                && !c1.negate
                && c2.description == "BondOrder"
                && !c2.negate
            {
                let (mut t1, mut t2) = (c1.val, c2.val);
                if t1 > t2 {
                    std::mem::swap(&mut t1, &mut t2);
                }
                let (single, double, aromatic) = (
                    BondType::Single as i32,
                    BondType::Double as i32,
                    BondType::Aromatic as i32,
                );
                if t1 == single && t2 == double {
                    return Ok("SD");
                } else if t1 == single && t2 == aromatic {
                    return Ok("SA");
                } else if t1 == double && t2 == aromatic {
                    return Ok("DA");
                }
            }
        }
    } else if !q.negate {
        match q.description.as_str() {
            "SingleOrAromaticBond" => return Ok("SA"),
            "SingleOrDoubleBond" => return Ok("SD"),
            "DoubleOrAromaticBond" => return Ok("DA"),
            "BondNull" => return Ok("Any"),
            _ => {}
        }
    }
    Err(writer_error(
        "Only SA, DA, SD, and Any query bond are supported for MarvinWriter",
    ))
}

const BOND_SYMBOL_ERROR: &str =
    "Only 1,2,3,Aromatic, and query bonds SA, DA, and SD are supported for MarvinWriter";

/// `GetMarvinBondSymbol`: `(order, queryType, convention)`
fn marvin_bond_symbol(bond: &Bond) -> Result<(String, String, String)> {
    if let Some(q) = bond.query() {
        let query_type = query_bond_symbol(q)?;
        if query_type.is_empty() {
            return Err(writer_error(BOND_SYMBOL_ERROR));
        }
        return Ok(("1".into(), query_type.into(), String::new()));
    }
    let (order, convention) = match bond.bond_type() {
        BondType::Single => (if bond.is_aromatic() { "A" } else { "1" }, ""),
        BondType::Double => (if bond.is_aromatic() { "A" } else { "2" }, ""),
        BondType::Triple => ("3", ""),
        BondType::Aromatic => ("A", ""),
        BondType::Dative => ("", "cxn:coord"),
        BondType::Hydrogen => ("", "cxn:hydrogen"),
        _ => return Err(writer_error(BOND_SYMBOL_ERROR)),
    };
    Ok((order.into(), String::new(), convention.into()))
}

/// `hasNonDefaultValence`
fn has_non_default_valence(atom: &Atom) -> bool {
    if atom.num_radical_electrons() != 0 {
        return true;
    }
    if atom.has_query() {
        return false;
    }
    !(atom.atomic_num() == 1 || in_organic_subset(atom.atomic_num() as i32))
}

fn radical_name(n: u32) -> Option<&'static str> {
    Some(match n {
        1 => "monovalent",
        2 => "divalent",
        3 => "trivalent4",
        4 => "4",
        _ => return None,
    })
}

struct Builder {
    marvin: MarvinMol,
    mol_count: u32,
    atom_count: u32,
    bond_count: u32,
    sg_count: u32,
}

impl Builder {
    fn add_group(&mut self, kind: SgroupKind) -> NodeId {
        self.sg_count += 1;
        self.mol_count += 1;
        let mut node = MarvinNode::new(kind, Some(ROOT));
        node.id = format!("sg{}", self.sg_count);
        node.mol_id = format!("m{}", self.mol_count);
        self.marvin.nodes.push(node);
        let id = self.marvin.nodes.len() - 1;
        self.marvin.nodes[ROOT].children.push(id);
        id
    }

    /// Adds the atoms of the substance group to the node; for the kinds
    /// that refer to their group the atoms get the reference.
    fn add_atoms(&mut self, node: NodeId, atoms: &[u32], set_ref: bool) {
        for &i in atoms {
            let a = self.marvin.nodes[ROOT].atoms[i as usize];
            self.marvin.nodes[node].atoms.push(a);
            if set_ref && !self.marvin.atoms[a].sgroup_ref_is_superatom {
                self.marvin.atoms[a].sgroup_ref = self.marvin.nodes[node].id.clone();
            }
        }
    }
}

fn sgroup_to_marvin(b: &mut Builder, sg: &SubstanceGroup) -> Result<()> {
    let typ =
        string_prop(sg, "TYPE").ok_or_else(|| writer_error("TYPE not found for an Sgroup"))?;
    match typ.as_str() {
        "SRU" | "MOD" | "COP" => {
            let role = match typ.as_str() {
                "SRU" => "SruSgroup",
                "MOD" => "ModificationSgroup",
                _ => "CopolymerSgroup",
            };
            let title = string_prop(sg, "LABEL").ok_or_else(|| {
                writer_error("Expected a LABEL attribute for an SRU, MOD, or COP group")
            })?;
            let connect = string_prop(sg, "CONNECT").ok_or_else(|| {
                writer_error("Expected a CONNECT attribute for an SRU, MOD, or COP group")
            })?;
            let node = b.add_group(SgroupKind::SruCoMod {
                role_name: role.into(),
                title,
                connect,
                correspondence: String::new(),
            });
            b.add_atoms(node, sg.atoms(), true);
            for &i in sg.bonds() {
                let bond = b.marvin.nodes[ROOT].bonds[i as usize];
                b.marvin.nodes[node].bonds.push(bond);
            }
        }
        "DAT" => {
            let fields = match sg.get_prop("DATAFIELDS") {
                Some(RDValue::VecString(v)) => v.join("\n"),
                _ => String::new(),
            };
            let double = |key: &str| match sg.get_prop(key) {
                Some(v) => v.as_double().unwrap_or(0.0),
                None => 0.0,
            };
            let node = b.add_group(SgroupKind::Data {
                context: string_prop(sg, "CONTEXT").unwrap_or_default(),
                field_name: string_prop(sg, "FIELDNAME").unwrap_or_default(),
                placement: string_prop(sg, "PLACEMENT").unwrap_or_default(),
                units_displayed: string_prop(sg, "UNITSDISPLAYED").unwrap_or_default(),
                query_type: string_prop(sg, "QUERYTYPE").unwrap_or_default(),
                query_op: string_prop(sg, "QUERYOP").unwrap_or_default(),
                field_data: fields,
                units: string_prop(sg, "UNITS").unwrap_or_default(),
                x: double("X"),
                y: double("Y"),
            });
            b.add_atoms(node, sg.atoms(), true);
        }
        "SUP" => {
            let title = string_prop(sg, "LABEL")
                .ok_or_else(|| writer_error("LABEL not found for a SuperatomSgroup"))?;
            let node = b.add_group(SgroupKind::SuperatomExpanded { title });
            b.add_atoms(node, sg.atoms(), false);
        }
        "MUL" => {
            let title = string_prop(sg, "MULT")
                .or_else(|| string_prop(sg, "LABEL"))
                .ok_or_else(|| writer_error("Title not found for a MultipleSgroup"))?;
            let node = b.add_group(SgroupKind::Multiple {
                title,
                is_expanded: true,
                parent_atoms: Vec::new(),
                bonds_to_atoms_not_in_expanded_group: Vec::new(),
            });
            b.add_atoms(node, sg.atoms(), true);
            let parents: Vec<usize> = sg
                .parent_atoms()
                .iter()
                .map(|&i| b.marvin.nodes[ROOT].atoms[i as usize])
                .collect();
            if let SgroupKind::Multiple { parent_atoms, .. } = &mut b.marvin.nodes[node].kind {
                *parent_atoms = parents;
            }
        }
        "GEN" => {
            // RDKit has no place for the charge location: assume onAtoms
            let node = b.add_group(SgroupKind::Generic {
                charge: "onAtoms".into(),
            });
            b.add_atoms(node, sg.atoms(), true);
        }
        "MON" => {
            let title = string_prop(sg, "MULT")
                .or_else(|| string_prop(sg, "LABEL"))
                .ok_or_else(|| writer_error("Title not found for a MultipleSgroup"))?;
            // (sic) upstream does not add the atoms
            b.add_group(SgroupKind::Monomer {
                title,
                charge: "onAtoms".into(),
            });
        }
        _ => {}
    }
    Ok(())
}

/// `MolToMarvinMol`: the molecule must have its property cache; the
/// numbering counters continue from `counts` (reactions chain several
/// molecules).
pub fn mol_to_marvin_mol(
    mol: &mut RWMol,
    counts: &mut NumberingCounts,
    conf_id: i32,
) -> Result<MarvinMol> {
    let mut b = Builder {
        marvin: MarvinMol::default(),
        mol_count: 1,
        atom_count: 0,
        bond_count: 0,
        sg_count: 0,
    };
    b.marvin.nodes[ROOT].mol_id = "m1".into();

    // pick the conformers
    let n_confs = mol.num_conformers() as i32;
    let mut conf2d: Option<Conformer> = None;
    let mut conf3d: Option<Conformer> = None;
    if n_confs > 0 {
        if conf_id >= 0 && conf_id < n_confs {
            let c = mol.conformer(conf_id)?.clone();
            if c.is_3d() {
                conf3d = Some(c);
            } else {
                conf2d = Some(c);
            }
        }
        if conf2d.is_none() && conf3d.is_none() {
            for c in mol.conformers() {
                if c.is_3d() {
                    conf3d.get_or_insert_with(|| c.clone());
                } else {
                    conf2d.get_or_insert_with(|| c.clone());
                }
                if conf2d.is_some() && conf3d.is_some() {
                    break;
                }
            }
        }
    }

    if mol.needs_update_property_cache() {
        mol.update_property_cache(false)?;
    }
    for idx in 0..mol.num_atoms() {
        let atom = mol.atom(idx);
        let mut ma = MarvinAtom::default();
        b.atom_count += 1;
        ma.id = format!("a{}", b.atom_count);
        marvin_atom_info(atom, &mut ma)?;
        ma.formal_charge = atom.formal_charge();
        let n_rad = atom.num_radical_electrons();
        if n_rad != 0
            && let Some(name) = radical_name(n_rad)
        {
            ma.radical = name.into();
        }
        if ma.is_element() {
            ma.isotope = atom.isotope() as i32;
            if ma.radical.is_empty() && has_non_default_valence(atom) {
                if mol.total_degree(idx)? == 0 {
                    // zero valence for elements/metals without neighbors
                    ma.mrv_valence = -1;
                } else if atom.no_implicit() && atom.num_explicit_hs() > 0 {
                    ma.hydrogen_count = atom.num_explicit_hs() as i32;
                } else {
                    let total = atom.total_valence()?;
                    if total != 15 {
                        ma.mrv_valence = (total % 15) as i32;
                    }
                }
            }
        }
        ma.mrv_map = atom
            .props
            .get(cp::molAtomMapNumber)
            .and_then(RDValue::as_int)
            .unwrap_or(0);
        if let Some(c) = &conf2d {
            let p = c.atom_pos(idx)?;
            ma.x2 = p.x;
            ma.y2 = p.y;
        }
        if let Some(c) = &conf3d {
            let p = c.atom_pos(idx)?;
            ma.x3 = p.x;
            ma.y3 = p.y;
            ma.z3 = p.z;
        }
        b.marvin.atoms.push(ma);
        let id = b.marvin.atoms.len() - 1;
        b.marvin.nodes[ROOT].atoms.push(id);
    }

    let conf_to_use = conf2d.as_ref().or(conf3d.as_ref());
    let wedge = pick_bonds_to_wedge(mol, conf_to_use)?;
    for bidx in 0..mol.num_bonds() {
        let bond = mol.bond(bidx);
        let mut mb = MarvinBond::default();
        b.bond_count += 1;
        mb.id = format!("b{}", b.bond_count);
        (mb.order, mb.query_type, mb.convention) = marvin_bond_symbol(bond)?;
        let (code, reverse) = match conf_to_use {
            Some(c) => get_mol_file_bond_stereo_info(mol, bidx, &wedge, Some(c))?,
            None => (0, false),
        };
        let atoms = &b.marvin.nodes[ROOT].atoms;
        let (first, second) = if reverse {
            (bond.end_atom_idx(), bond.begin_atom_idx())
        } else {
            (bond.begin_atom_idx(), bond.end_atom_idx())
        };
        mb.atom_refs2 = [
            b.marvin.atoms[atoms[first as usize]].id.clone(),
            b.marvin.atoms[atoms[second as usize]].id.clone(),
        ];
        let dir = match code {
            1 => BondDir::BeginWedge,
            6 => BondDir::BeginDash,
            4 => BondDir::Unknown,
            3 => BondDir::EitherDouble,
            _ => BondDir::None,
        };
        match dir {
            BondDir::BeginWedge => mb.bond_stereo.value = "W".into(),
            BondDir::BeginDash => mb.bond_stereo.value = "H".into(),
            BondDir::Unknown => {
                mb.bond_stereo.convention = "MDL".into();
                mb.bond_stereo.convention_value = "4".into();
            }
            BondDir::EitherDouble => {
                mb.bond_stereo.convention = "MDL".into();
                mb.bond_stereo.convention_value = "3".into();
            }
            _ => {}
        }
        b.marvin.bonds.push(mb);
        let id = b.marvin.bonds.len() - 1;
        b.marvin.nodes[ROOT].bonds.push(id);
    }

    // stereo groups
    let (mut or_count, mut and_count) = (0, 0);
    for group in mol.stereo_groups() {
        let name = match group.group_type() {
            StereoGroupType::Absolute => "abs".to_string(),
            StereoGroupType::Or => {
                or_count += 1;
                format!("or{or_count}")
            }
            StereoGroupType::And => {
                and_count += 1;
                format!("and{and_count}")
            }
        };
        for id in all_atom_ids_for_stereo_group(mol, group, &wedge) {
            let a = b.marvin.nodes[ROOT].atoms[id as usize];
            b.marvin.atoms[a].mrv_stereo_group = name.clone();
        }
    }

    for sg in mol.substance_groups() {
        sgroup_to_marvin(&mut b, sg)?;
    }

    b.marvin.process_sgroups_from_rdkit(ROOT)?;
    b.marvin.clean_up_numbering(counts);
    Ok(b.marvin)
}

/// `MolToMrvBlock`
pub fn mol_to_mrv_block(mol: &ROMol, params: &MrvWriterParams, conf_id: i32) -> Result<String> {
    let mut m = mol.clone();
    if m.needs_update_property_cache() {
        m.update_property_cache(false)?;
    }
    if params.kekulize {
        kekulize(&mut m, true, true, 100)?;
    }
    if params.include_stereo && m.num_conformers() == 0 {
        // coordinates, so that the stereo we write makes sense
        crate::graphmol::depictor::compute_2d_coords(&mut m, &Default::default())?;
    }
    crate::graphmol::generic_groups::convert_generic_queries_to_substance_groups(&mut m)?;
    let mut counts = NumberingCounts::default();
    let mut marvin = mol_to_marvin_mol(&mut m, &mut counts, conf_id)?;
    marvin.nodes[ROOT].coordinate_precision = params.precision;
    Ok(write_xml(&marvin.document_tree(), params.pretty_print))
}

/// `MolToMrvFile`
pub fn mol_to_mrv_file(
    mol: &ROMol,
    path: &str,
    params: &MrvWriterParams,
    conf_id: i32,
) -> Result<()> {
    let text = mol_to_mrv_block(mol, params, conf_id)?;
    std::fs::write(path, text).map_err(|_| RDKitError::BadFile(format!("Bad output file {path}")))
}

#[allow(dead_code)]
fn unused(_: f64) -> f64 {
    NO_COORD
}
