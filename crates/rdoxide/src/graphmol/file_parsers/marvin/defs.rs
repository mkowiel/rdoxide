//! Port of RDKit `Code/GraphMol/MarvinParse/MarvinDefs.{h,cpp}`: the Marvin
//! (CML) molecule model with its superatom, multiple, SRU, data ... groups
//! and the transformations that make them palatable to RDKit.
//!
//! Upstream links atoms, bonds and groups with raw pointers; here a
//! [`MarvinMol`] owns all of them in arenas and the groups form a tree of
//! [`MarvinNode`]s that refer to atoms and bonds by index.  Where upstream
//! identifies groups by their `id` string (and not by identity) this port
//! does too, so that its quirks are reproduced.

use super::xml::XmlTree;
use crate::{RDKitError, Result};

pub type AtomId = usize;
pub type BondId = usize;
pub type NodeId = usize;

/// The marker for "no such coordinate" (`DBL_MAX`).
pub const NO_COORD: f64 = f64::MAX;

pub const SRU_SGROUP_CONNECT_CHOICES: [&str; 3] = ["hh", "ht", "eu"];
pub const MARVIN_BOND_ORDERS: [&str; 4] = ["1", "2", "3", "A"];
pub const MARVIN_QUERY_BOND_TYPES: [&str; 4] = ["SD", "SA", "DA", "Any"];
pub const MARVIN_CONVENTION_TYPES: [&str; 2] = ["cxn:coord", "cxn:hydrogen"];
pub const MARVIN_STEREO_DICT_REF_TYPES: [&str; 2] = ["cml:W", "cml:H"];
pub const MARVIN_STEREO_CONVENTION_TYPES: [&str; 4] = ["1", "3", "4", "6"];
pub const MARVIN_RADICAL_VALS: [&str; 8] = [
    "monovalent",
    "divalent",
    "divalent1",
    "divalent3",
    "trivalent",
    "trivalent2",
    "trivalent4",
    "4",
];

/// `marvinRadicalToRadicalElectrons`
pub fn radical_electrons(radical: &str) -> Option<u32> {
    Some(match radical {
        "monovalent" => 1,
        "divalent" | "divalent1" | "divalent3" => 2,
        "trivalent" | "trivalent2" | "trivalent4" => 3,
        "4" => 4,
        _ => return None,
    })
}

fn file_parse(msg: impl Into<String>) -> RDKitError {
    RDKitError::FileParse(msg.into())
}

/// `getCleanNumber`: a number with no surrounding white space.
pub fn clean_number<T: std::str::FromStr>(s: &str) -> Option<T> {
    if s.trim() != s {
        return None;
    }
    s.parse().ok()
}

/// `boost::algorithm::split(.., is_space())`: splits at every white space
/// character, keeping empty fields (an empty string gives one empty field).
fn split_space(s: &str) -> Vec<String> {
    s.split(|c: char| c.is_ascii_whitespace())
        .map(str::to_string)
        .collect()
}

/// A Marvin atom.
#[derive(Debug, Clone)]
pub struct MarvinAtom {
    pub id: String,
    pub element_type: String,
    pub x2: f64,
    pub y2: f64,
    pub x3: f64,
    pub y3: f64,
    pub z3: f64,
    pub formal_charge: i32,
    pub radical: String,
    pub isotope: i32,
    pub mrv_valence: i32,
    pub hydrogen_count: i32,
    pub mrv_alias: String,
    pub mrv_stereo_group: String,
    pub mrv_map: i32,
    pub sgroup_ref: String,
    /// when set, the `sgroup_ref` is that of a superatom and is not changed
    pub sgroup_ref_is_superatom: bool,
    pub sgroup_attachment_point: String,
    /// -1: not specified
    pub rgroup_ref: i32,
}

impl Default for MarvinAtom {
    fn default() -> Self {
        MarvinAtom {
            id: String::new(),
            element_type: String::new(),
            x2: NO_COORD,
            y2: NO_COORD,
            x3: NO_COORD,
            y3: NO_COORD,
            z3: NO_COORD,
            formal_charge: 0,
            radical: String::new(),
            isotope: 0,
            mrv_valence: -1,
            hydrogen_count: -1,
            mrv_alias: String::new(),
            mrv_stereo_group: String::new(),
            mrv_map: 0,
            sgroup_ref: String::new(),
            sgroup_ref_is_superatom: false,
            sgroup_attachment_point: String::new(),
            rgroup_ref: -1,
        }
    }
}

impl MarvinAtom {
    /// `MarvinAtom(atomToCopy, newId)`
    pub fn copy_with_id(&self, new_id: String) -> MarvinAtom {
        MarvinAtom {
            id: new_id,
            ..self.clone()
        }
    }

    pub fn is_element(&self) -> bool {
        self.element_type != "R" && self.element_type != "*"
    }
}

#[derive(Debug, Clone, Default)]
pub struct MarvinBondStereo {
    pub value: String,
    pub convention: String,
    pub convention_value: String,
    pub dict_ref: String,
}

/// A Marvin bond.
#[derive(Debug, Clone, Default)]
pub struct MarvinBond {
    pub id: String,
    pub atom_refs2: [String; 2],
    pub order: String,
    pub bond_stereo: MarvinBondStereo,
    pub query_type: String,
    pub convention: String,
}

impl MarvinBond {
    /// `MarvinBond(bondToCopy, newId, atomRef1, atomRef2)`
    pub fn copy_with(&self, new_id: String, ref1: String, ref2: String) -> MarvinBond {
        MarvinBond {
            id: new_id,
            atom_refs2: [ref1, ref2],
            order: self.order.clone(),
            bond_stereo: self.bond_stereo.clone(),
            query_type: self.query_type.clone(),
            convention: self.convention.clone(),
        }
    }

    /// `getBondType`: SD, SA, DA, ANY, DATIVE, HYDROGEN, 1, 2, 3 or A.
    pub fn bond_type(&self) -> Result<String> {
        let query_type = self.query_type.to_uppercase();
        let order = self.order.to_uppercase();
        let convention = self.convention.to_uppercase();
        if !query_type.is_empty() {
            if ["SD", "SA", "DA", "ANY"].contains(&query_type.as_str()) {
                Ok(query_type)
            } else {
                Err(file_parse(format!(
                    "unrecognized query bond type {} in MRV File ",
                    self.query_type
                )))
            }
        } else if !convention.is_empty() {
            match convention.as_str() {
                "CXN:COORD" => Ok("DATIVE".into()),
                "CXN:HYDROGEN" => Ok("HYDROGEN".into()),
                _ => Err(file_parse(format!(
                    "unrecognized convention {} in MRV File ",
                    self.convention
                ))),
            }
        } else if !order.is_empty() {
            if ["1", "2", "3", "A"].contains(&order.as_str()) {
                Ok(order)
            } else {
                Err(file_parse(format!(
                    "unrecognized bond type {} in MRV File ",
                    self.order
                )))
            }
        } else {
            Err(file_parse(
                "bond must have one of:  order, queryType, or convention in MRV File ",
            ))
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct MarvinAttachmentPoint {
    pub atom: String,
    pub bond: String,
    pub order: String,
}

/// What kind of group a node is, with the data specific to it.
#[derive(Debug, Clone)]
pub enum SgroupKind {
    /// the molecule itself (the root of the tree)
    Mol,
    /// `MarvinSruCoModSgroup`: SruSgroup, CopolymerSgroup or ModificationSgroup
    SruCoMod {
        role_name: String,
        title: String,
        connect: String,
        correspondence: String,
    },
    Data {
        context: String,
        field_name: String,
        placement: String,
        units_displayed: String,
        query_type: String,
        query_op: String,
        field_data: String,
        units: String,
        x: f64,
        y: f64,
    },
    SuperatomExpanded {
        title: String,
    },
    Multiple {
        title: String,
        is_expanded: bool,
        parent_atoms: Vec<AtomId>,
        /// only when expanded
        bonds_to_atoms_not_in_expanded_group: Vec<BondId>,
    },
    Multicenter {
        center: AtomId,
    },
    Generic {
        charge: String,
    },
    Monomer {
        title: String,
        charge: String,
    },
    Superatom {
        title: String,
        attachment_points: Vec<MarvinAttachmentPoint>,
    },
}

impl SgroupKind {
    pub fn role(&self) -> &str {
        match self {
            SgroupKind::Mol => "MarvinMol",
            SgroupKind::SruCoMod { role_name, .. } => role_name,
            SgroupKind::Data { .. } => "DataSgroup",
            SgroupKind::SuperatomExpanded { .. } => "SuperatomSgroupExpanded",
            SgroupKind::Multiple { .. } => "MultipleSgroup",
            SgroupKind::Multicenter { .. } => "MulticenterSgroup",
            SgroupKind::Generic { .. } => "GenericSgroup",
            SgroupKind::Monomer { .. } => "MonomerSgroup",
            SgroupKind::Superatom { .. } => "SuperatomSgroup",
        }
    }
}

/// How the atoms of a group relate to a set of atoms
/// (`IsSgroupInAtomSetResult`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InAtomSet {
    In,
    NotIn,
    BothInAndNotIn,
}

/// A molecule or a group of one (`MarvinMolBase`).
#[derive(Debug, Clone)]
pub struct MarvinNode {
    pub mol_id: String,
    pub id: String,
    pub coordinate_precision: u32,
    pub atoms: Vec<AtomId>,
    pub bonds: Vec<BondId>,
    /// the groups (and for the root, nothing else) below this one
    pub children: Vec<NodeId>,
    pub parent: Option<NodeId>,
    pub kind: SgroupKind,
}

impl MarvinNode {
    pub fn new(kind: SgroupKind, parent: Option<NodeId>) -> Self {
        MarvinNode {
            mol_id: String::new(),
            id: String::new(),
            coordinate_precision: 6,
            atoms: Vec::new(),
            bonds: Vec::new(),
            children: Vec::new(),
            parent,
            kind,
        }
    }

    pub fn role(&self) -> &str {
        self.kind.role()
    }

    pub fn has_atom_bond_blocks(&self) -> bool {
        matches!(self.kind, SgroupKind::Mol | SgroupKind::Superatom { .. })
    }

    /// `isPassiveRoleForExpansion`
    pub fn is_passive_role_for_expansion(&self) -> bool {
        !matches!(
            self.kind,
            SgroupKind::Superatom { .. } | SgroupKind::Multiple { .. }
        )
    }

    /// `isPassiveRoleForContraction`
    pub fn is_passive_role_for_contraction(&self) -> bool {
        !matches!(
            self.kind,
            SgroupKind::Mol | SgroupKind::SuperatomExpanded { .. } | SgroupKind::Multiple { .. }
        )
    }
}

/// A Marvin molecule: the arenas and the tree of groups; node 0 is the
/// molecule itself.
#[derive(Debug, Clone)]
pub struct MarvinMol {
    pub atoms: Vec<MarvinAtom>,
    pub bonds: Vec<MarvinBond>,
    pub nodes: Vec<MarvinNode>,
}

pub const ROOT: NodeId = 0;

impl Default for MarvinMol {
    fn default() -> Self {
        MarvinMol {
            atoms: Vec::new(),
            bonds: Vec::new(),
            nodes: vec![MarvinNode::new(SgroupKind::Mol, None)],
        }
    }
}

fn required(xml: &XmlTree, name: &str, msg: &str) -> Result<String> {
    let v = xml.attr(name, "");
    if v.is_empty() {
        Err(file_parse(msg))
    } else {
        Ok(v)
    }
}

impl MarvinMol {
    /// `MarvinMol(molTree)` and the sub-molecules below it.
    pub fn parse(mol_tree: &XmlTree) -> Result<MarvinMol> {
        let mut mol = MarvinMol::default();
        mol.parse_atoms_and_bonds(ROOT, mol_tree)?;
        mol.parse_children(ROOT, mol_tree)?;
        Ok(mol)
    }

    /// The sub-molecules of `tree` become the groups of `node`.
    fn parse_children(&mut self, node: NodeId, tree: &XmlTree) -> Result<()> {
        for sub in tree.children_named("molecule") {
            let child = self.parse_sub_molecule(node, sub)?;
            self.nodes[node].children.push(child);
            self.parse_children(child, sub)?;
        }
        Ok(())
    }

    // ------------------------------------------------------------------
    // lookups

    pub fn node(&self, id: NodeId) -> &MarvinNode {
        &self.nodes[id]
    }

    fn add_node(&mut self, node: MarvinNode) -> NodeId {
        self.nodes.push(node);
        self.nodes.len() - 1
    }

    /// `getAtomIndex`: the position of the atom in the atoms of `node`.
    pub fn atom_index(&self, node: NodeId, id: &str) -> i32 {
        self.nodes[node]
            .atoms
            .iter()
            .position(|&a| self.atoms[a].id == id)
            .map_or(-1, |p| p as i32)
    }

    pub fn bond_index(&self, node: NodeId, id: &str) -> i32 {
        self.nodes[node]
            .bonds
            .iter()
            .position(|&b| self.bonds[b].id == id)
            .map_or(-1, |p| p as i32)
    }

    /// `findAtomByRef`: in the atoms of `node` or, failing that, of its
    /// ancestors.
    pub fn find_atom_by_ref(&self, node: NodeId, id: &str) -> Option<AtomId> {
        let mut cur = Some(node);
        while let Some(n) = cur {
            if let Some(&a) = self.nodes[n]
                .atoms
                .iter()
                .find(|&&a| self.atoms[a].id == id)
            {
                return Some(a);
            }
            cur = self.nodes[n].parent;
        }
        None
    }

    pub fn find_bond_by_ref(&self, node: NodeId, id: &str) -> Option<BondId> {
        let mut cur = Some(node);
        while let Some(n) = cur {
            if let Some(&b) = self.nodes[n]
                .bonds
                .iter()
                .find(|&&b| self.bonds[b].id == id)
            {
                return Some(b);
            }
            cur = self.nodes[n].parent;
        }
        None
    }

    pub fn atom_list(&self, node: NodeId) -> Vec<String> {
        self.nodes[node]
            .atoms
            .iter()
            .map(|&a| self.atoms[a].id.clone())
            .collect()
    }

    pub fn bond_list(&self, node: NodeId) -> Vec<String> {
        self.nodes[node]
            .bonds
            .iter()
            .map(|&b| self.bonds[b].id.clone())
            .collect()
    }

    pub fn has_2d_coords(&self, node: NodeId) -> bool {
        self.nodes[node]
            .atoms
            .iter()
            .all(|&a| self.atoms[a].x2 != NO_COORD && self.atoms[a].y2 != NO_COORD)
    }

    pub fn has_any_2d_coords(&self, node: NodeId) -> bool {
        self.nodes[node]
            .atoms
            .iter()
            .any(|&a| self.atoms[a].x2 != NO_COORD && self.atoms[a].y2 != NO_COORD)
    }

    pub fn has_any_3d_coords(&self, node: NodeId) -> bool {
        self.nodes[node].atoms.iter().any(|&a| {
            let at = &self.atoms[a];
            at.x3 != NO_COORD && at.y3 != NO_COORD && at.z3 != NO_COORD
        })
    }

    /// `getExplicitValence`
    pub fn explicit_valence(&self, node: NodeId, atom: &MarvinAtom) -> Result<i32> {
        let mut times10 = 0u32;
        for &b in &self.nodes[node].bonds {
            let bond = &self.bonds[b];
            if bond.atom_refs2[0] != atom.id && bond.atom_refs2[1] != atom.id {
                continue;
            }
            times10 += match bond.bond_type()?.as_str() {
                "SD" | "SA" | "DA" | "A" => 15,
                // no good answer for Any bonds: a single bond
                "ANY" | "1" => 10,
                "2" => 20,
                "3" => 30,
                _ => 0, // dative and hydrogen bonds
            };
        }
        Ok((times10 / 10) as i32)
    }

    /// `getActualParent`: the first ancestor that is not a MultipleSgroup.
    pub fn actual_parent(&self, child: NodeId) -> NodeId {
        let mut cur = self.nodes[child].parent.expect("a group has a parent");
        while self.nodes[cur].role() == "MultipleSgroup" {
            cur = self.nodes[cur].parent.expect("a group has a parent");
        }
        cur
    }

    // ------------------------------------------------------------------
    // parsing

    fn atoms_from_refs(
        &self,
        node: NodeId,
        parent: NodeId,
        refs: &str,
        what: &str,
    ) -> Result<Vec<AtomId>> {
        let _ = node;
        split_space(refs)
            .iter()
            .map(|r| {
                self.find_atom_by_ref(parent, r).ok_or_else(|| {
                    file_parse(format!(
                        "AtomRef specification for {what} definition was not found in any parent"
                    ))
                })
            })
            .collect()
    }

    /// `parseAtomsAndBonds`
    fn parse_atoms_and_bonds(&mut self, node: NodeId, mol_tree: &XmlTree) -> Result<()> {
        let atom_array = mol_tree
            .child("atomArray")
            .ok_or_else(|| RDKitError::FileParse("No such node (atomArray)".into()))?;
        let atom_id = atom_array.attr("atomID", "UseLongForm");
        if atom_id == "UseLongForm" {
            for (name, v) in &atom_array.children {
                if name != "atom" {
                    continue;
                }
                let atom = self.parse_long_form_atom(v)?;
                self.atoms.push(atom);
                self.nodes[node].atoms.push(self.atoms.len() - 1);
            }
        } else {
            self.parse_short_form_atoms(node, atom_array, &atom_id)?;
        }
        if let Some(bond_array) = mol_tree.child("bondArray") {
            for (name, v) in &bond_array.children {
                if name != "bond" {
                    continue;
                }
                let bond = self.parse_bond(node, v)?;
                self.bonds.push(bond);
                self.nodes[node].bonds.push(self.bonds.len() - 1);
            }
        }
        Ok(())
    }

    fn parse_long_form_atom(&self, v: &XmlTree) -> Result<MarvinAtom> {
        let mut atom = MarvinAtom {
            id: v.attr("id", ""),
            element_type: v.attr("elementType", ""),
            ..MarvinAtom::default()
        };
        if atom.id.is_empty() || atom.element_type.is_empty() {
            return Err(file_parse(
                "Expected id, elementType for an atom definition in MRV file",
            ));
        }
        let (x2, y2) = (v.attr("x2", ""), v.attr("y2", ""));
        if !x2.is_empty() && !y2.is_empty() {
            match (clean_number(&x2), clean_number(&y2)) {
                (Some(x), Some(y)) => (atom.x2, atom.y2) = (x, y),
                _ => {
                    return Err(file_parse(
                        "The values x2 and y2 must be large floats in MRV file",
                    ));
                }
            }
        }
        let (x3, y3, z3) = (v.attr("x3", ""), v.attr("y3", ""), v.attr("z3", ""));
        if !x3.is_empty() && !y3.is_empty() && !z3.is_empty() {
            match (clean_number(&x3), clean_number(&y3), clean_number(&z3)) {
                (Some(x), Some(y), Some(z)) => (atom.x3, atom.y3, atom.z3) = (x, y, z),
                _ => {
                    return Err(file_parse(
                        "The values x3, y3,  and z2 must be large floats in MRV file",
                    ));
                }
            }
        }
        let charge = v.attr("formalCharge", "");
        if !charge.is_empty() {
            atom.formal_charge = clean_number(&charge).ok_or_else(|| {
                file_parse("The value for formalCharge must be an integer in MRV file")
            })?;
        }
        atom.radical = v.attr("radical", "");
        if !atom.radical.is_empty() && !MARVIN_RADICAL_VALS.contains(&atom.radical.as_str()) {
            return Err(file_parse(format!(
                "The value for radical must be one of {} in MRV file",
                MARVIN_RADICAL_VALS.join(", ")
            )));
        }
        let isotope = v.attr("isotope", "");
        if !isotope.is_empty() {
            atom.isotope = clean_number(&isotope).filter(|&i| i > 0).ok_or_else(|| {
                file_parse("The value for isotope must be a positive number in MRV file")
            })?;
        }
        let valence = v.attr("mrvValence", "");
        if !valence.is_empty() {
            atom.mrv_valence = clean_number(&valence).filter(|&i| i >= 0).ok_or_else(|| {
                file_parse("The value for mrvValence must be a positive number in MRV file")
            })?;
        }
        let hcount = v.attr("hydrogenCount", "");
        if !hcount.is_empty() {
            atom.hydrogen_count = clean_number(&hcount).filter(|&i| i >= 0).ok_or_else(|| {
                file_parse("The value for hydrogenCount must be a non-negative number in MRV file")
            })?;
        }
        atom.mrv_alias = v.attr("mrvAlias", "");
        // (boost's get<int> throws a ptree_bad_data for a non-integer)
        let rgroup = v.attr("rgroupRef", "");
        atom.rgroup_ref = if v
            .child("<xmlattr>")
            .and_then(|a| a.child("rgroupRef"))
            .is_some()
        {
            clean_number(&rgroup).ok_or_else(|| {
                RDKitError::FileParse("conversion of data to type \"int\" failed".into())
            })?
        } else {
            -1
        };
        atom.mrv_stereo_group = v.attr("mrvStereoGroup", "");
        if atom.mrv_stereo_group == "0" {
            atom.mrv_stereo_group.clear();
        }
        let map = v.attr("mrvMap", "");
        if !map.is_empty() {
            atom.mrv_map = clean_number(&map).filter(|&i| i > 0).ok_or_else(|| {
                file_parse("The value for mrvMap must be an non-=negative integer in MRV file")
            })?;
        }
        atom.sgroup_ref = v.attr("sgroupRef", "");
        atom.sgroup_attachment_point = v.attr("sgroupAttachmentPoint", "");
        Ok(atom)
    }

    fn parse_short_form_atoms(
        &mut self,
        node: NodeId,
        atom_array: &XmlTree,
        atom_id: &str,
    ) -> Result<()> {
        let atom_ids: Vec<String> = if atom_id.is_empty() {
            Vec::new()
        } else {
            split_space(atom_id)
        };
        let atom_count = atom_ids.len();
        let column = |name: &str| -> (String, Vec<String>) {
            let text = atom_array.attr(name, "");
            let parts = split_space(&text);
            (text, parts)
        };
        let (element_type, element_types) = column("elementType");
        let (x2, x2s) = column("x2");
        let (y2, y2s) = column("y2");
        let (formal_charge, formal_charges) = column("formalCharge");
        let (isotope, isotopes) = column("isotope");
        let (radical, radicals) = column("radical");
        let (hydrogen_count, hydrogen_counts) = column("hydrogenCount");
        let (mrv_valence, mrv_valences) = column("mrvValence");
        let (mrv_alias, mrv_aliases) = column("mrvAlias");
        let (rgroup_ref, rgroup_refs) = column("rgroupRef");
        let (mrv_stereo_group, mrv_stereo_groups) = column("mrvStereoGroup");
        let (mrv_map, mrv_maps) = column("mrvMap");
        let (sgroup_ref, sgroup_refs) = column("sgroupRef");
        let (sgroup_attachment_point, sgroup_attachment_points) = column("sgroupAttachmentPoint");

        if !atom_id.is_empty() {
            if element_type.is_empty() {
                return Err(file_parse(
                    "Expected an elementType array for an atomArray definition in MRV file",
                ));
            }
            if element_types.len() < atom_count {
                return Err(file_parse("There must be an element type for each atom id"));
            }
        }
        for i in 0..atom_count {
            let mut atom = MarvinAtom {
                id: atom_ids[i].clone(),
                element_type: element_types[i].clone(),
                ..MarvinAtom::default()
            };
            if !x2.is_empty() && !y2.is_empty() && x2s.len() > i && y2s.len() > i {
                match (clean_number(&x2s[i]), clean_number(&y2s[i])) {
                    (Some(x), Some(y)) => (atom.x2, atom.y2) = (x, y),
                    _ => {
                        return Err(file_parse(
                            "The values x2 and y2 must be large floats in MRV file",
                        ));
                    }
                }
            }
            if !formal_charge.is_empty() && formal_charges.len() > i {
                atom.formal_charge = clean_number(&formal_charges[i]).ok_or_else(|| {
                    file_parse("The value for formalCharge must be an integer in MRV file")
                })?;
            }
            if !isotope.is_empty() && isotopes.len() > i {
                atom.isotope = clean_number(&isotopes[i]).ok_or_else(|| {
                    file_parse("The value for formalCharge must be an integer in MRV file")
                })?;
            }
            if !mrv_valence.is_empty() && mrv_valences.len() > i {
                if mrv_valences[i] != "-" {
                    atom.mrv_valence = clean_number(&mrv_valences[i]).ok_or_else(|| {
                        file_parse("The value for mrvValences must be an integer in MRV file")
                    })?;
                }
            }
            if !hydrogen_count.is_empty() && hydrogen_counts.len() > i && hydrogen_counts[i] != "-"
            {
                atom.hydrogen_count = clean_number(&hydrogen_counts[i]).ok_or_else(|| {
                    file_parse("The value for hydrogenCount must be an integer in MRV file")
                })?;
            }
            if !radical.is_empty() && radicals.len() > i && radicals[i] != "0" {
                atom.radical = radicals[i].clone();
                if !MARVIN_RADICAL_VALS.contains(&atom.radical.as_str()) {
                    return Err(file_parse(format!(
                        "The value for radical must be one of {} in MRV file",
                        MARVIN_RADICAL_VALS.join(", ")
                    )));
                }
            }
            if !mrv_alias.is_empty() && mrv_aliases.len() > i {
                atom.mrv_alias = mrv_aliases[i].clone();
            }
            if !rgroup_ref.is_empty() && rgroup_refs.len() > i {
                atom.rgroup_ref = clean_number(&rgroup_refs[i])
                    .ok_or_else(|| file_parse("rgroupRef value must be an integer in MRV file"))?;
            }
            if !mrv_stereo_group.is_empty()
                && mrv_stereo_groups.len() > i
                && mrv_stereo_groups[i] != "0"
            {
                atom.mrv_stereo_group = mrv_stereo_groups[i].clone();
            }
            if !mrv_map.is_empty() && mrv_maps.len() > i {
                atom.mrv_map = clean_number(&mrv_maps[i])
                    .filter(|&m| m >= 0)
                    .ok_or_else(|| {
                        file_parse(
                            "The value for mrvMap must be an non-negative integer in MRV file",
                        )
                    })?;
            }
            if !sgroup_ref.is_empty() && sgroup_refs.len() > i && sgroup_refs[i] != "0" {
                atom.sgroup_ref = sgroup_refs[i].clone();
            }
            if !sgroup_attachment_point.is_empty()
                && sgroup_attachment_points.len() > i
                && sgroup_attachment_points[i] != "0"
            {
                atom.sgroup_attachment_point = sgroup_attachment_points[i].clone();
            }
            self.atoms.push(atom);
            self.nodes[node].atoms.push(self.atoms.len() - 1);
        }
        Ok(())
    }

    fn parse_bond(&self, node: NodeId, v: &XmlTree) -> Result<MarvinBond> {
        let mut bond = MarvinBond {
            id: v.attr("id", ""),
            ..MarvinBond::default()
        };
        if bond.id.is_empty() {
            return Err(file_parse("Expected id for an bond definition in MRV file"));
        }
        let refs = split_space(&v.attr("atomRefs2", ""));
        let known = |r: &String| {
            self.nodes[node]
                .atoms
                .iter()
                .any(|&a| self.atoms[a].id == *r)
        };
        if refs.len() != 2 || !known(&refs[0]) || !known(&refs[1]) {
            return Err(file_parse(
                "atomRefs2 must contain two atom refs that must appear in the atoms array in MRV file",
            ));
        }
        bond.atom_refs2 = [refs[0].clone(), refs[1].clone()];
        bond.order = v.attr("order", "");
        if !bond.order.is_empty() && !MARVIN_BOND_ORDERS.contains(&bond.order.as_str()) {
            return Err(file_parse(format!(
                "Expected one of  {} for order for an bond definition in MRV file",
                MARVIN_BOND_ORDERS.join(", ")
            )));
        }
        bond.query_type = v.attr("queryType", "");
        if !bond.query_type.is_empty()
            && !MARVIN_QUERY_BOND_TYPES.contains(&bond.query_type.as_str())
        {
            return Err(file_parse(format!(
                "Expected one of  {} for queryType for an bond definition in MRV file",
                MARVIN_QUERY_BOND_TYPES.join(", ")
            )));
        }
        bond.convention = v.attr("convention", "");
        if !bond.convention.is_empty()
            && !MARVIN_CONVENTION_TYPES.contains(&bond.convention.as_str())
        {
            return Err(file_parse(format!(
                "Expected one of  {} for convention for an bond definition in MRV file",
                MARVIN_CONVENTION_TYPES.join(", ")
            )));
        }
        let mut stereo_decl_count = 0;
        bond.bond_stereo.value = v.get_str("bondStereo", "");
        if !bond.bond_stereo.value.is_empty() {
            stereo_decl_count += 1;
            match bond.bond_stereo.value.to_lowercase().as_str() {
                "w" | "h" => {}
                // cis and trans are ignored
                "c" | "t" => bond.bond_stereo.value.clear(),
                _ => {
                    return Err(file_parse(
                        "The value for bondStereo must be \"H\", \"W\", \"C\" or \"T\" in MRV file (\"C\" and \"T\" are ignored)",
                    ));
                }
            }
        }
        if v.child("bondStereo").is_some() {
            // upstream loops over all the children of the bond; the last one
            // that is not an attribute node decides (the others just reset
            // the values to empty)
            for (_, ww) in &v.children {
                bond.bond_stereo.convention = ww.attr("convention", "");
                if !bond.bond_stereo.convention.is_empty() {
                    stereo_decl_count += 1;
                    if bond.bond_stereo.convention != "MDL" {
                        return Err(file_parse(
                            "Expected MDL as value for the bond convention attribute",
                        ));
                    }
                    bond.bond_stereo.convention_value = ww.attr("conventionValue", "");
                    if !MARVIN_STEREO_CONVENTION_TYPES
                        .contains(&bond.bond_stereo.convention_value.as_str())
                    {
                        return Err(file_parse(format!(
                            "Expected one of  {} for a bond convention for an bond stereo def",
                            MARVIN_STEREO_CONVENTION_TYPES.join(", ")
                        )));
                    }
                }
                bond.bond_stereo.dict_ref = ww.attr("dictRef", "");
                if !bond.bond_stereo.dict_ref.is_empty() {
                    stereo_decl_count += 1;
                    if !MARVIN_STEREO_DICT_REF_TYPES.contains(&bond.bond_stereo.dict_ref.as_str()) {
                        return Err(file_parse(format!(
                            "Expected one of  {} for a discRef value for an bond stereo def",
                            MARVIN_STEREO_DICT_REF_TYPES.join(", ")
                        )));
                    }
                }
            }
        }
        if stereo_decl_count > 1 {
            return Err(file_parse(
                "bondStereo can either have only one of: a value, dictRef, or convention value",
            ));
        }
        Ok(bond)
    }

    /// `parseMarvinMolecule` for a sub-molecule (a group).
    fn parse_sub_molecule(&mut self, parent: NodeId, tree: &XmlTree) -> Result<NodeId> {
        let role = tree.attr("role", "");
        if role.is_empty() {
            return Err(file_parse("Expected a role for a sub-molecule in MRV file"));
        }
        let ids = |tree: &XmlTree| -> Result<(String, String)> {
            let id = required(tree, "id", "Expected a iD in MRV file")?;
            let mol_id = required(tree, "molID", "Expected a molID in MRV file")?;
            Ok((id, mol_id))
        };
        let atom_refs_required = |tree: &XmlTree, what: &str| -> Result<String> {
            required(
                tree,
                "atomRefs",
                &format!("Expected  atomRefs for a {what} definition in MRV file"),
            )
        };
        let (id, mol_id) = ids(tree)?;
        let mut node;
        match role.as_str() {
            "SuperatomSgroup" => {
                let atom_refs = tree.attr("atomRefs", "");
                if atom_refs.is_empty() && tree.child("atomArray").is_some() {
                    // a regular (contracted) superatom
                    let title = required(
                        tree,
                        "title",
                        "Expected  title for a SuperatomSgroup definition in MRV file",
                    )?;
                    node = MarvinNode::new(
                        SgroupKind::Superatom {
                            title,
                            attachment_points: Vec::new(),
                        },
                        Some(parent),
                    );
                    node.id = id;
                    node.mol_id = mol_id;
                    let nid = self.add_node(node);
                    self.parse_atoms_and_bonds(nid, tree)?;
                    self.parse_attachment_points(nid, parent, tree)?;
                    return Ok(nid);
                }
                node = MarvinNode::new(
                    SgroupKind::SuperatomExpanded {
                        title: String::new(),
                    },
                    Some(parent),
                );
                if !atom_refs.is_empty() {
                    node.atoms = self.atoms_from_refs(
                        ROOT,
                        parent,
                        &atom_refs,
                        "an SuperatomSgroupExpanded group",
                    )?;
                } else {
                    // no atoms nor atomRefs: the atoms that refer to this group
                    node.atoms = self.nodes[parent]
                        .atoms
                        .iter()
                        .copied()
                        .filter(|&a| self.atoms[a].sgroup_ref == id)
                        .collect();
                }
                let title = tree.attr("title", "");
                if title.is_empty() {
                    return Err(file_parse(
                        "Expected  title for a SuperatomSgroupExpanded definition in MRV file",
                    ));
                }
                node.kind = SgroupKind::SuperatomExpanded { title };
            }
            "SruSgroup" | "CopolymerSgroup" | "ModificationSgroup" => {
                let refs = atom_refs_required(tree, "SruSgroup")?;
                let atoms = self.atoms_from_refs(
                    ROOT,
                    parent,
                    &refs,
                    "an SRU, Copolymer, or Modification group",
                )?;
                let title = tree.attr("title", "");
                if title.is_empty() {
                    return Err(file_parse(
                        "Expected title for a SruSgroup definition in MRV file",
                    ));
                }
                let mut connect = tree.attr("connect", "");
                if connect.is_empty() {
                    connect = "ht".into();
                }
                if !SRU_SGROUP_CONNECT_CHOICES.contains(&connect.as_str()) {
                    return Err(file_parse(format!(
                        "Expected a connect  string of \"{}\" for a SruSgroup definition in MRV file",
                        SRU_SGROUP_CONNECT_CHOICES.join(", ")
                    )));
                }
                let mut bonds = Vec::new();
                let bond_list = tree.attr("bondList", "");
                if !bond_list.is_empty() {
                    for b in split_space(&bond_list) {
                        bonds.push(self.find_bond_by_ref(parent, &b).ok_or_else(|| {
                            file_parse(
                                "BondList specification for an SRU, Copolymer, or Modification group definition was not found in any parent",
                            )
                        })?);
                    }
                }
                node = MarvinNode::new(
                    SgroupKind::SruCoMod {
                        role_name: role.clone(),
                        title,
                        connect,
                        correspondence: tree.attr("correspondence", ""),
                    },
                    Some(parent),
                );
                node.atoms = atoms;
                node.bonds = bonds;
            }
            "MultipleSgroup" => {
                let refs = atom_refs_required(tree, "MultipleSgroup")?;
                let atoms = self.atoms_from_refs(ROOT, parent, &refs, "a MultipleSgroup")?;
                let title = tree.attr("title", "");
                if clean_number::<i32>(&title).is_none_or(|n| n <= 0) {
                    return Err(file_parse(
                        "Expected a positive integer title for a MultipleSgroup definition in MRV file",
                    ));
                }
                node = MarvinNode::new(
                    SgroupKind::Multiple {
                        title,
                        is_expanded: false,
                        parent_atoms: atoms.clone(),
                        bonds_to_atoms_not_in_expanded_group: Vec::new(),
                    },
                    Some(parent),
                );
                node.atoms = atoms;
            }
            "MulticenterSgroup" => {
                let refs = atom_refs_required(tree, "MulticenterSgroup")?;
                let atoms = self.atoms_from_refs(ROOT, parent, &refs, "a MulticenterSgroup")?;
                let center_id = tree.attr("center", "");
                let center = self.find_atom_by_ref(parent, &center_id).ok_or_else(|| {
                    file_parse(
                        "Center specification for a MulticenterSgroup definition was not found in any parent",
                    )
                })?;
                node = MarvinNode::new(SgroupKind::Multicenter { center }, Some(parent));
                node.atoms = atoms;
            }
            "GenericSgroup" => {
                let refs = atom_refs_required(tree, "GenericSgroup")?;
                let atoms = self.atoms_from_refs(ROOT, parent, &refs, "a GenericSgroup")?;
                let charge = tree.attr("charge", "");
                if charge != "onAtoms" && charge != "onBracket" {
                    return Err(file_parse(
                        "Expected  omAtoms or onBracket for a charge attr of a GenericSgroup definition in MRV file",
                    ));
                }
                node = MarvinNode::new(SgroupKind::Generic { charge }, Some(parent));
                node.atoms = atoms;
            }
            "MonomerSgroup" => {
                let refs = atom_refs_required(tree, "MonomerSgroup")?;
                let atoms = self.atoms_from_refs(ROOT, parent, &refs, "a MonomerSgroup")?;
                let title = tree.attr("title", "");
                if title.is_empty() {
                    return Err(file_parse(
                        "Expected  title for a MonomerSgroup definition in MRV file",
                    ));
                }
                let charge = tree.attr("charge", "");
                if charge.is_empty() {
                    return Err(file_parse(
                        "Expected  omAtoms or onBracket for a charge attr of a MonomerSgroup definition in MRV file",
                    ));
                }
                node = MarvinNode::new(SgroupKind::Monomer { title, charge }, Some(parent));
                node.atoms = atoms;
            }
            "DataSgroup" => {
                let refs = atom_refs_required(tree, "DataSgroup")?;
                let atoms = self.atoms_from_refs(ROOT, parent, &refs, "a DataSgroup")?;
                let mut units_displayed = tree.attr("unitsDisplayed", "Unit not displayed");
                if units_displayed.is_empty() {
                    units_displayed = "Unit not displayed".into();
                }
                let lower = units_displayed.to_lowercase();
                if lower != "unit displayed" && lower != "unit not displayed" {
                    return Err(file_parse(
                        "Expected unitsDisplayed to be either \"Unit displayed\" or \"Unit not displayed\" for a DataSgroup definition in MRV file",
                    ));
                }
                let x_text = tree.attr("x", "0.0");
                let y_text = tree.attr("y", "0.0");
                if x_text.is_empty() {
                    return Err(file_parse(
                        "Expected x for a DataSgroup definition in MRV file",
                    ));
                }
                let x = clean_number(&x_text).ok_or_else(|| {
                    file_parse("The value for x must be a floating point value in MRV file")
                })?;
                if y_text.is_empty() {
                    return Err(file_parse(
                        "Expected y for a DataSgroup definition in MRV file",
                    ));
                }
                let y = clean_number(&y_text).ok_or_else(|| {
                    file_parse("The value for y must be a floating point value in MRV file")
                })?;
                node = MarvinNode::new(
                    SgroupKind::Data {
                        context: tree.attr("context", ""),
                        field_name: tree.attr("fieldName", ""),
                        placement: tree.attr("placement", ""),
                        units_displayed,
                        query_type: tree.attr("queryType", ""),
                        query_op: tree.attr("queryOp", ""),
                        field_data: tree.attr("fieldData", ""),
                        units: tree.attr("units", ""),
                        x,
                        y,
                    },
                    Some(parent),
                );
                node.atoms = atoms;
            }
            other => {
                return Err(file_parse(format!("Unexpected role {other} in MRV file")));
            }
        }
        node.id = id;
        node.mol_id = mol_id;
        Ok(self.add_node(node))
    }

    /// The attachment points of a contracted superatom.
    fn parse_attachment_points(
        &mut self,
        node: NodeId,
        parent: NodeId,
        tree: &XmlTree,
    ) -> Result<()> {
        let Some(array) = tree.child("AttachmentPointArray") else {
            return Ok(());
        };
        let mut points = Vec::new();
        for (_, v) in &array.children {
            let bond = v.attr("bond", "");
            if bond.is_empty() {
                // can happen if the attachment point is not connected
                continue;
            }
            let point = MarvinAttachmentPoint {
                atom: v.attr("atom", ""),
                order: v.attr("order", ""),
                bond,
            };
            if point.atom.is_empty() || point.order.is_empty() {
                return Err(file_parse(
                    "Expected atom, order and bond, for an AttachmentPoint definition in MRV file",
                ));
            }
            if !self.nodes[node]
                .atoms
                .iter()
                .any(|&a| self.atoms[a].id == point.atom)
            {
                return Err(file_parse(
                    "Atom specification for an AttachmentPoint definition must be in the parent's atom array in MRV file",
                ));
            }
            if self.find_bond_by_ref(parent, &point.bond).is_none() {
                return Err(file_parse(
                    "Bond specification for an AttachmentPoint definition must be in the bond array in MRV file",
                ));
            }
            if clean_number::<i32>(&point.order).is_none() {
                return Err(file_parse(
                    "Order for an AttachmentPoint definition must be an integer in MRV file",
                ));
            }
            points.push(point);
        }
        if let SgroupKind::Superatom {
            attachment_points, ..
        } = &mut self.nodes[node].kind
        {
            *attachment_points = points;
        }
        Ok(())
    }
}
