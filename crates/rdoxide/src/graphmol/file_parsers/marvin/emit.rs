//! The XML trees of a Marvin molecule (`toPtree` in upstream
//! `MarvinDefs.cpp`).

use super::defs::{MarvinAtom, MarvinBond, MarvinMol, NO_COORD, NodeId, ROOT, SgroupKind};
use super::xml::XmlTree;
use crate::rdgeneral::strings::format_g;

/// `getDoubleAsText`
pub(super) fn double_as_text(val: f64, precision: u32) -> String {
    // (kept for backwards compatibility upstream)
    if precision == 6 && val.abs() < 0.00001 {
        return "0.00000".to_string();
    }
    format_g(val, precision as usize)
}

fn atom_tree(a: &MarvinAtom, precision: u32) -> XmlTree {
    let mut out = XmlTree::default();
    out.put_attr("id", a.id.clone());
    out.put_attr("elementType", a.element_type.clone());
    if a.x2 != NO_COORD && a.y2 != NO_COORD {
        out.put_attr("x2", double_as_text(a.x2, precision));
        out.put_attr("y2", double_as_text(a.y2, precision));
    }
    if a.x3 != NO_COORD && a.y3 != NO_COORD && a.z3 != NO_COORD {
        out.put_attr("x3", double_as_text(a.x3, precision));
        out.put_attr("y3", double_as_text(a.y3, precision));
        out.put_attr("z3", double_as_text(a.z3, precision));
    }
    if a.formal_charge != 0 {
        out.put_attr("formalCharge", a.formal_charge.to_string());
    }
    if !a.radical.is_empty() {
        out.put_attr("radical", a.radical.clone());
    }
    if a.is_element() && a.isotope != 0 {
        out.put_attr("isotope", a.isotope.to_string());
    }
    if !a.mrv_alias.is_empty() {
        out.put_attr("mrvAlias", a.mrv_alias.clone());
    }
    if a.mrv_valence >= 0 {
        out.put_attr("mrvValence", a.mrv_valence.to_string());
    }
    if a.hydrogen_count > 0 {
        out.put_attr("hydrogenCount", a.hydrogen_count.to_string());
    }
    if !a.mrv_stereo_group.is_empty() {
        out.put_attr("mrvStereoGroup", a.mrv_stereo_group.clone());
    }
    if a.mrv_map != 0 {
        out.put_attr("mrvMap", a.mrv_map.to_string());
    }
    if !a.sgroup_ref.is_empty() {
        out.put_attr("sgroupRef", a.sgroup_ref.clone());
    }
    if a.rgroup_ref >= 0 {
        out.put_attr("rgroupRef", a.rgroup_ref.to_string());
    }
    if !a.sgroup_attachment_point.is_empty() {
        out.put_attr("sgroupAttachmentPoint", a.sgroup_attachment_point.clone());
    }
    out
}

fn bond_tree(b: &MarvinBond) -> XmlTree {
    let mut out = XmlTree::default();
    out.put_attr("id", b.id.clone());
    out.put_attr(
        "atomRefs2",
        format!("{} {}", b.atom_refs2[0], b.atom_refs2[1]),
    );
    if !b.order.is_empty() {
        out.put_attr("order", b.order.clone());
    }
    if !b.query_type.is_empty() {
        out.put_attr("queryType", b.query_type.clone());
    }
    if !b.convention.is_empty() {
        out.put_attr("convention", b.convention.clone());
    }
    let s = &b.bond_stereo;
    if !s.value.is_empty()
        || !s.dict_ref.is_empty()
        || (!s.convention.is_empty() && !s.convention_value.is_empty())
    {
        let mut st = XmlTree::default();
        if !s.value.is_empty() {
            st.data = s.value.clone();
        } else if !s.dict_ref.is_empty() {
            st.put_attr("dictRef", s.dict_ref.clone());
        } else {
            st.put_attr("convention", s.convention.clone());
            st.put_attr("conventionValue", s.convention_value.clone());
        }
        out.put_child("bondStereo", st);
    }
    out
}

impl MarvinMol {
    /// `MarvinMolBase::toPtree` and the overrides
    pub fn node_tree(&self, node: NodeId) -> XmlTree {
        let n = self.node(node);
        let mut out = XmlTree::default();
        out.put_attr("molID", n.mol_id.clone());
        if n.parent.is_some() {
            out.put_attr("id", n.id.clone());
            out.put_attr("role", n.role().to_string());
        }
        if n.has_atom_bond_blocks() {
            let mut atoms = XmlTree::default();
            for &a in &n.atoms {
                atoms.add_child("atom", atom_tree(&self.atoms[a], n.coordinate_precision));
            }
            if !n.atoms.is_empty() {
                out.add_child("atomArray", atoms);
            }
            let mut bonds = XmlTree::default();
            for &b in &n.bonds {
                bonds.add_child("bond", bond_tree(&self.bonds[b]));
            }
            if !n.bonds.is_empty() {
                out.add_child("bondArray", bonds);
            }
        }
        let atom_refs = || self.atom_list(node).join(" ");
        let mut with_children = false;
        match &n.kind {
            SgroupKind::Mol => with_children = true,
            SgroupKind::SruCoMod {
                title,
                connect,
                correspondence,
                ..
            } => {
                out.put_attr("atomRefs", atom_refs());
                out.put_attr("title", title.clone());
                out.put_attr("connect", connect.clone());
                out.put_attr("correspondence", correspondence.clone());
                out.put_attr("bondList", self.bond_list(node).join(" "));
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
                out.put_attr("atomRefs", atom_refs());
                out.put_attr("context", context.clone());
                out.put_attr("fieldName", field_name.clone());
                out.put_attr("placement", placement.clone());
                out.put_attr("unitsDisplayed", units_displayed.clone());
                out.put_attr("fieldData", field_data.clone());
                if !units.is_empty() {
                    out.put_attr("units", units.clone());
                }
                if !query_type.is_empty() && !query_op.is_empty() {
                    out.put_attr("queryType", query_type.clone());
                    out.put_attr("queryOp", query_op.clone());
                }
                out.put_attr("x", format_g(*x, 6));
                out.put_attr("y", format_g(*y, 6));
            }
            SgroupKind::SuperatomExpanded { title } => {
                out.put_attr("atomRefs", atom_refs());
                out.put_attr("title", title.clone());
            }
            SgroupKind::Multiple { title, .. } => {
                out.put_attr("atomRefs", atom_refs());
                out.put_attr("title", title.clone());
                with_children = true;
            }
            SgroupKind::Multicenter { center } => {
                out.put_attr("atomRefs", atom_refs());
                out.put_attr("center", self.atoms[*center].id.clone());
            }
            SgroupKind::Generic { .. } => {
                out.put_attr("atomRefs", atom_refs());
            }
            SgroupKind::Monomer { title, charge } => {
                out.put_attr("atomRefs", atom_refs());
                out.put_attr("title", title.clone());
                out.put_attr("charge", charge.clone());
            }
            SgroupKind::Superatom {
                title,
                attachment_points,
            } => {
                out.put_attr("title", title.clone());
                if !attachment_points.is_empty() {
                    let mut arr = XmlTree::default();
                    for p in attachment_points {
                        let mut t = XmlTree::default();
                        t.put_attr("atom", p.atom.clone());
                        t.put_attr("order", p.order.clone());
                        t.put_attr("bond", p.bond.clone());
                        arr.add_child("attachmentPoint", t);
                    }
                    out.put_child("AttachmentPointArray", arr);
                }
                with_children = true;
            }
        }
        if with_children {
            for &c in &n.children {
                out.add_child("molecule", self.node_tree(c));
            }
        }
        out
    }

    /// `toMolPtree`: the whole document.
    pub fn document_tree(&self) -> XmlTree {
        let mut cml = XmlTree::default();
        cml.put_attr("xmlns", "http://www.chemaxon.com");
        cml.put_attr("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
        cml.put_attr(
            "xsi:schemaLocation",
            "http://www.chemaxon.com http://www.chemaxon.com/marvin/schema/mrvSchema_20_20_0.xsd",
        );
        let mut chem = XmlTree::default();
        chem.add_child("molecule", self.node_tree(ROOT));
        let mut doc = XmlTree::default();
        doc.add_child("MChemicalStruct", chem);
        cml.add_child("MDocument", doc);
        let mut root = XmlTree::default();
        root.add_child("cml", cml);
        root
    }
}
