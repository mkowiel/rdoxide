//! The SVG specifics of the drawer (`MolDraw2DSVG.cpp` of RDKit 2026.09.1):
//! the header and background, finishing off, the molecule metadata and the
//! invisible shapes that tag atoms and bonds for interactive use.
//!
//! Copyright (C) 2015-2021 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;
use std::fmt::Write;

use super::details::format_double;
use super::drawer::MolDraw2DSVG;
use super::helpers::colour_to_svg;
use super::svg_canvas::stream_double;
use crate::Result;
use crate::geometry::Point2D;
use crate::graphmol::ROMol;
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, get_atom_smiles, get_bond_smiles};
use crate::rdgeneral::{RDProps, RDValue};

const RDKIT_SVG_VERSION: &str = "0.9";

/// The classes from the `_tagClass` property of an atom or bond.
fn output_tag_classes<T: RDProps>(obj: &T, os: &mut String, active_class: &str) {
    if !active_class.is_empty() {
        os.push(' ');
        os.push_str(active_class);
    }
    if let Some(RDValue::String(value)) = obj.get_prop("_tagClass") {
        os.push(' ');
        os.push_str(&value.replace(['"', '\'', '.'], "_"));
    }
}

/// The ` name="value"` attributes for the `_metaData-` properties.
fn output_meta_data<T: RDProps>(obj: &T, os: &mut String) {
    for prop in obj.get_prop_list(true, true) {
        if prop.len() < 11 || !prop.starts_with("_metaData-") {
            continue;
        }
        if let Some(value) = obj.get_prop(&prop) {
            let value = value.to_rd_string().replace('"', "&quot;");
            let _ = write!(os, " {}=\"{}\"", &prop[10..], value);
        }
    }
}

impl MolDraw2DSVG {
    /// Writes the start of the SVG.
    pub(super) fn init_drawing(&mut self) {
        let out = &mut self.canvas.out;
        out.push_str("<?xml version='1.0' encoding='iso-8859-1'?>\n");
        out.push_str(
            "<svg version='1.1' baseProfile='full'\n              \
             xmlns='http://www.w3.org/2000/svg'\n                      \
             xmlns:rdkit='http://www.rdkit.org/xml'\n                      \
             xmlns:xlink='http://www.w3.org/1999/xlink'\n                  \
             xml:space='preserve'\n",
        );
        let _ = writeln!(
            out,
            "width='{w}px' height='{h}px' viewBox='0 0 {w} {h}'>",
            w = self.width,
            h = self.height
        );
        out.push_str("<!-- END OF HEADER -->\n");
    }

    /// Paints the background.
    pub fn clear_drawing(&mut self) {
        if self.needs_init {
            self.init_drawing();
            self.needs_init = false;
        }
        let col = colour_to_svg(&self.canvas.options.background_colour)
            .unwrap_or_else(|_| "#FFFFFF".to_string());
        let offset = self.offset();
        let _ = writeln!(
            self.canvas.out,
            "<rect style='opacity:1.0;fill:{col};stroke:none' width='{}' height='{}' x='{}' y='{}'> </rect>",
            format_double(f64::from(self.width)),
            format_double(f64::from(self.height)),
            format_double(offset.x),
            format_double(offset.y)
        );
    }

    /// Ends the SVG.
    pub fn finish_drawing(&mut self) {
        self.canvas.out.push_str("</svg>\n");
    }

    /// The SVG written so far.
    pub fn get_drawing_text(&self) -> String {
        self.canvas.out.clone()
    }

    /// Adds the `<metadata>` for a molecule, with its atoms and bonds as
    /// SMILES and their drawing and molecule coordinates.
    pub fn add_molecule_metadata(&mut self, mol: &ROMol, conf_id: i32) -> Result<()> {
        let mut os = String::new();
        os.push_str("<metadata>\n");
        let _ = writeln!(
            os,
            "<rdkit:mol xmlns:rdkit = \"http://www.rdkit.org/xml\" version=\"{RDKIT_SVG_VERSION}\">"
        );
        let atom_params = SmilesWriteParams {
            do_kekule: false,
            all_hs_explicit: true,
            do_isomeric_smiles: true,
            ..Default::default()
        };
        let conf = mol.conformer(conf_id)?;
        for idx in 0..mol.num_atoms() {
            let _ = write!(os, "<rdkit:atom idx=\"{}\"", idx + 1);
            let _ = write!(
                os,
                " atom-smiles=\"{}\"",
                get_atom_smiles(mol, idx, &atom_params)?
            );
            let pos = conf.atom_pos(idx)?;
            let dpos = self.get_draw_coords(Point2D::new(pos.x, pos.y))?;
            let _ = write!(
                os,
                " drawing-x=\"{}\" drawing-y=\"{}\"",
                format_double(dpos.x),
                format_double(dpos.y)
            );
            let _ = write!(
                os,
                " x=\"{}\" y=\"{}\" z=\"{}\"",
                stream_double(pos.x),
                stream_double(pos.y),
                stream_double(pos.z)
            );
            output_meta_data(mol.atom(idx), &mut os);
            os.push_str(" />\n");
        }
        let bond_params = SmilesWriteParams {
            do_kekule: false,
            all_bonds_explicit: true,
            ..Default::default()
        };
        // writing a bond's SMILES touches the bond, so it works on a copy
        let mut smiles_mol = mol.clone();
        for idx in 0..mol.num_bonds() {
            let bond = mol.bond(idx);
            let _ = write!(os, "<rdkit:bond idx=\"{}\"", idx + 1);
            let _ = write!(os, " begin-atom-idx=\"{}\"", bond.begin_atom_idx() + 1);
            let _ = write!(os, " end-atom-idx=\"{}\"", bond.end_atom_idx() + 1);
            let _ = write!(
                os,
                " bond-smiles=\"{}\"",
                get_bond_smiles(&mut smiles_mol, idx, &bond_params, -1)?
            );
            output_meta_data(bond, &mut os);
            os.push_str(" />\n");
        }
        os.push_str("</rdkit:mol></metadata>\n");
        self.canvas.out.push_str(&os);
        Ok(())
    }

    /// Adds the metadata of several molecules, one after the other.
    pub fn add_molecules_metadata(&mut self, mols: &[&ROMol], conf_ids: &[i32]) -> Result<()> {
        for (i, mol) in mols.iter().enumerate() {
            let conf_id = if conf_ids.len() == mols.len() {
                conf_ids[i]
            } else {
                -1
            };
            self.set_active_mol_idx(i as i32)?;
            self.add_molecule_metadata(mol, conf_id)?;
        }
        Ok(())
    }

    /// Adds invisible shapes over the bonds and atoms of the active molecule,
    /// with classes by which they can be found, and the given event handlers
    /// (`name` -> function called with the atom index) on the atoms.
    pub fn tag_atoms(
        &mut self,
        mol: &ROMol,
        radius: f64,
        events: &BTreeMap<String, String>,
    ) -> Result<()> {
        let active_class = self.canvas.active_class.clone();
        let split_bonds = self.canvas.options.split_bonds;
        let width = 2.0 + self.canvas.options.bond_line_width;
        let mut os = String::new();
        // first bonds so that they are under the atoms
        for idx in 0..mol.num_bonds() {
            let bond = mol.bond(idx);
            let (a_idx1, a_idx2) = (bond.begin_atom_idx(), bond.end_atom_idx());
            let a1pos = self.get_draw_coords_atom(a_idx1 as usize)?;
            let a2pos = self.get_draw_coords_atom(a_idx2 as usize)?;
            let style = format!(
                " style='fill:#fff;stroke:#fff;stroke-width:{}px;fill-opacity:0;stroke-opacity:0' ",
                format_double(width)
            );
            let segment = |from: Point2D, to: Point2D, atoms: &str, os: &mut String| {
                let _ = write!(
                    os,
                    "<path  d='M {},{} L {},{}' class='bond-selector bond-{idx} {atoms}",
                    format_double(from.x),
                    format_double(from.y),
                    format_double(to.x),
                    format_double(to.y),
                );
                output_tag_classes(bond, os, &active_class);
                os.push('\'');
                os.push_str(&style);
                os.push_str("/>\n");
            };
            if split_bonds {
                let midp = (a1pos + a2pos) * 0.5;
                segment(a1pos, midp, &format!("atom-{a_idx1}"), &mut os);
                segment(midp, a2pos, &format!("atom-{a_idx2}"), &mut os);
            } else {
                segment(
                    a1pos,
                    a2pos,
                    &format!("atom-{a_idx1} atom-{a_idx2}"),
                    &mut os,
                );
            }
        }
        let scale = self.active_scale_for_tags()?;
        for idx in 0..mol.num_atoms() {
            let pos = self.get_draw_coords_atom(idx as usize)?;
            let _ = write!(
                os,
                "<circle  cx='{}' cy='{}' r='{}' class='atom-selector atom-{idx}",
                format_double(pos.x),
                format_double(pos.y),
                format_double(scale * radius)
            );
            output_tag_classes(mol.atom(idx), &mut os, &active_class);
            os.push('\'');
            os.push_str(
                " style='fill:#fff;stroke:#fff;stroke-width:1px;fill-opacity:0;stroke-opacity:0' ",
            );
            for (name, func) in events {
                let _ = write!(os, " {name}='{func}({idx});'");
            }
            os.push_str("/>\n");
        }
        self.canvas.out.push_str(&os);
        Ok(())
    }
}
