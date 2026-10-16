//! Drawing chemical reactions (the reaction functions of `MolDraw2D.cpp` of
//! RDKit 2026.09.1): reactants, agents and products side by side with plus
//! signs and an arrow, all to one scale.
//!
//! Copyright (C) 2015-2022 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeMap;
use std::f64::consts::PI;

use super::draw_mol::{DrawMol, DrawMolHighlights};
use super::draw_mol_utils::center_mol_for_drawing;
use super::drawer::MolDraw2DSVG;
use super::helpers::DrawColour;
use super::utils::prepare_mol_for_drawing;
use crate::Result;
use crate::geometry::Point2D;
use crate::graphmol::chem_reactions::ChemicalReaction;
use crate::graphmol::conjug_hybrid::set_hybridization;
use crate::graphmol::depictor::{Compute2DCoordParameters, compute_2d_coords};
use crate::graphmol::kekulize::{DEFAULT_MAX_BACKTRACKS, kekulize_if_possible};
use crate::graphmol::{ROMol, RWMol};
use crate::rdgeneral::log::LogStateSetter;

/// The three kinds of component of a reaction, as ranges of the drawer's
/// molecules (which are in this order: reactants, products, agents).
struct Components {
    reagents: std::ops::Range<usize>,
    products: std::ops::Range<usize>,
    agents: std::ops::Range<usize>,
}

fn reaction_width(
    dms: &[DrawMol],
    c: &Components,
    padding: f64,
    arrow_mult: i32,
    gap_width: i32,
) -> i32 {
    let mut tot = 0i32;
    if !c.reagents.is_empty() {
        tot += c.reagents.clone().map(|i| dms[i].width).sum::<i32>();
        tot += gap_width * (c.reagents.len() as i32 - 1);
        tot += gap_width / 2;
    }
    if c.agents.is_empty() {
        tot += arrow_mult * gap_width;
    } else {
        // the agent doesn't start at front of arrow
        tot += gap_width / 2;
        for i in c.agents.clone() {
            tot += dms[i].width + gap_width / 2;
        }
    }
    if !c.products.is_empty() {
        tot += gap_width / 2;
        tot += c.products.clone().map(|i| dms[i].width).sum::<i32>();
        // we don't want a plus after the last product
        tot += gap_width * (c.products.len() as i32 - 1);
    }
    (f64::from(tot) + 2.0 * padding * f64::from(tot)) as i32
}

impl MolDraw2DSVG {
    /// Draws a reaction: reactants, agents above the arrow and products. With
    /// `highlight_by_reactant` atoms and bonds are coloured by the reactant
    /// they come from, and atom maps are not shown.
    pub fn draw_reaction(
        &mut self,
        rxn: &ChemicalReaction,
        highlight_by_reactant: bool,
        highlight_colours: Option<&[DrawColour]>,
        conf_ids: Option<&[i32]>,
    ) -> Result<()> {
        // Copy the reaction because processing for drawing alters it.
        let mut nrxn = rxn.clone();
        let (comps, mut plus_width) = self.get_reaction_draw_mols(
            &mut nrxn,
            highlight_by_reactant,
            highlight_colours,
            conf_ids,
        )?;

        if self.height == -1 {
            for dm in &self.draw_mols {
                self.height = self.height.max(dm.height);
            }
            self.height += (f64::from(self.height) * 2.0 * self.canvas.options.padding) as i32;
            self.panel_height = self.height;
        }

        let (offsets, arrow_beg, arrow_end) =
            self.calc_reaction_offsets(&comps, &mut plus_width)?;
        self.active_mol_idx = -1;
        self.start_drawing();
        let mut x_offset = 0;
        x_offset = self.draw_reaction_part(&comps.reagents, plus_width, x_offset, &offsets)?;
        // make sure the arrowhead is big enough
        let mut frac = 0.05;
        let delta = (arrow_end - arrow_beg).length() * frac;
        // an arbitrary 5 pixel minimum
        if delta < 5.0 {
            frac *= 5.0 / delta;
        }
        x_offset = self.draw_reaction_part(&comps.agents, 0, x_offset, &offsets)?;
        self.draw_reaction_part(&comps.products, plus_width, x_offset, &offsets)?;
        let osbw = self.canvas.options.scale_bond_width;
        if comps.reagents.is_empty() && comps.products.is_empty() && comps.agents.is_empty() {
            // an empty reaction has no DrawMol with a scale, so can't have a
            // scaled bond width
            self.canvas.options.scale_bond_width = false;
        }
        let symbol_colour = self.canvas.options.symbol_colour;
        let res = self.draw_arrow(
            arrow_beg,
            arrow_end,
            false,
            frac,
            PI / 6.0,
            symbol_colour,
            true,
        );
        self.canvas.options.scale_bond_width = osbw;
        res
    }

    fn get_reaction_draw_mols(
        &mut self,
        rxn: &mut ChemicalReaction,
        highlight_by_reactant: bool,
        highlight_colours: Option<&[DrawColour]>,
        conf_ids: Option<&[i32]>,
    ) -> Result<(Components, i32)> {
        let agent_frac = 0.4;
        let mut min_scale = f64::MAX;
        let mut min_font_scale = 0.0;

        let atom_colours =
            self.find_reaction_highlights(rxn, highlight_by_reactant, highlight_colours);
        // reactants and products; at the end min_scale and min_font_scale are
        // the smallest used in any of the drawings
        let padding = self.canvas.options.padding;
        let use_height = if self.panel_height > 0 {
            (f64::from(self.panel_height) * (1.0 - 2.0 * padding)) as i32
        } else {
            self.panel_height
        };
        let first = self.draw_mols.len();
        let n_react = rxn.reactants().len();
        let n_prod = rxn.products().len();
        let n_agent = rxn.agents().len();
        let mut reactants: Vec<ROMol> = rxn.reactants().to_vec();
        self.make_reaction_components(
            &mut reactants,
            conf_ids,
            use_height,
            &atom_colours,
            &mut min_scale,
            &mut min_font_scale,
        )?;
        let mut products: Vec<ROMol> = rxn.products().to_vec();
        self.make_reaction_components(
            &mut products,
            conf_ids,
            use_height,
            &atom_colours,
            &mut min_scale,
            &mut min_font_scale,
        )?;
        let comps = Components {
            reagents: first..first + n_react,
            products: first + n_react..first + n_react + n_prod,
            agents: first + n_react + n_prod..first + n_react + n_prod + n_agent,
        };
        // Now scale everything to min_scale, so all the elements are drawn to
        // the same scale.
        for i in comps.reagents.clone().chain(comps.products.clone()) {
            self.draw_mols[i].set_scale(min_scale, min_font_scale, false, &mut self.text)?;
            self.draw_mols[i].shrink_to_fit(false, &mut self.text)?;
        }
        // set the spacing for the plus signs to be 0.5 Angstrom
        let plus_width = min_scale as i32;

        // agents
        let agent_height = (agent_frac * f64::from(use_height)) as i32;
        min_scale = f64::MAX;
        let mut agents: Vec<ROMol> = rxn.agents().to_vec();
        self.make_reaction_components(
            &mut agents,
            conf_ids,
            agent_height,
            &atom_colours,
            &mut min_scale,
            &mut min_font_scale,
        )?;
        for i in comps.agents.clone() {
            self.draw_mols[i].set_scale(min_scale, min_font_scale, false, &mut self.text)?;
            self.draw_mols[i].shrink_to_fit(false, &mut self.text)?;
        }

        // set the active atom/bond indices so they run in series for all
        // pieces: reagents, products then agents
        let (mut atom_idx_offset, mut bond_idx_offset) = (0, 0);
        for dm in &mut self.draw_mols {
            dm.active_atm_idx_offset = atom_idx_offset;
            dm.active_bnd_idx_offset = bond_idx_offset;
            atom_idx_offset += dm.draw_mol.num_atoms() as i32;
            bond_idx_offset += dm.draw_mol.num_bonds() as i32;
        }
        Ok((comps, plus_width))
    }

    fn make_reaction_components(
        &mut self,
        bits: &mut [ROMol],
        conf_ids: Option<&[i32]>,
        height_to_use: i32,
        atom_colours: &BTreeMap<i32, DrawColour>,
        min_scale: &mut f64,
        min_font_scale: &mut f64,
    ) -> Result<()> {
        for (midx, frag) in bits.iter_mut().enumerate() {
            let cid = conf_ids.map_or(-1, |c| c[midx]);
            let mut hl = DrawMolHighlights::default();
            for ai in 0..frag.num_atoms() {
                let colour = atom_colours.get(&frag.atom(ai).atom_map_num()).copied();
                for bi in 0..frag.num_bonds() {
                    let bond = frag.bond(bi);
                    let beg = frag.atom(bond.begin_atom_idx()).atom_map_num();
                    if let Some(begc) = atom_colours.get(&beg) {
                        let end = frag.atom(bond.end_atom_idx()).atom_map_num();
                        if atom_colours.get(&end).is_some_and(|endc| endc == begc) {
                            hl.bonds.push(bi as i32);
                            hl.bond_map.entry(bi as i32).or_insert(*begc);
                        }
                    }
                }
                if let Some(c) = colour {
                    hl.atoms.push(ai as i32);
                    hl.atom_map.entry(ai as i32).or_insert(c);
                    frag.atom_mut(ai).set_atom_map_num(0, false)?;
                }
            }
            self.make_reaction_draw_mol(frag, cid, height_to_use, hl)?;
            let dm = self.draw_mols.last().expect("just added");
            if dm.scale < *min_scale {
                *min_scale = dm.scale;
                *min_font_scale = dm.font_scale;
            }
        }
        Ok(())
    }

    fn make_reaction_draw_mol(
        &mut self,
        mol: &mut RWMol,
        conf_id: i32,
        mol_height: i32,
        highlights: DrawMolHighlights,
    ) -> Result<()> {
        mol.update_property_cache(false)?;
        if self.canvas.options.prepare_mols_before_drawing {
            let _blocker = LogStateSetter::new();
            kekulize_if_possible(mol, false, false, DEFAULT_MAX_BACKTRACKS)?;
            set_hybridization(mol);
        }
        if mol.num_conformers() == 0 {
            let params = Compute2DCoordParameters {
                canon_orient: true,
                ..Default::default()
            };
            compute_2d_coords(mol, &params)?;
        } else {
            // we need to center the molecule
            center_mol_for_drawing(mol, conf_id)?;
        }
        // when preparing a reaction component to be drawn we should neither
        // kekulize (we did that above if required) nor add chiral Hs
        prepare_mol_for_drawing(mol, false, false, true, false, false)?;
        // The height is fixed, but the width is allowed to be as large as the
        // height and molecule dimensions dictate. Padding is turned off to
        // create the molecule - it is handled for the whole picture.
        let padding = self.canvas.options.padding;
        let mut opts = self.canvas.options.clone();
        opts.padding = 0.0;
        let mut dm = DrawMol::new(
            mol,
            "",
            -1,
            mol_height,
            &opts,
            &mut self.text,
            highlights,
            true,
            conf_id,
            true,
        )?;
        dm.opts.padding = padding;
        dm.create_draw_objects(&mut self.text)?;
        self.draw_mols.push(dm);
        self.active_mol_idx += 1;
        Ok(())
    }

    fn calc_reaction_offsets(
        &mut self,
        c: &Components,
        plus_width: &mut i32,
    ) -> Result<(Vec<Point2D>, Point2D, Point2D)> {
        const ARROW_MULT: i32 = 2; // number of plus widths for an empty arrow
        let padding = self.canvas.options.padding;
        let mut width_to_use = if self.panel_width == -1 {
            self.width
        } else {
            self.panel_width
        };
        if width_to_use == -1 {
            width_to_use = reaction_width(&self.draw_mols, c, padding, ARROW_MULT, *plus_width);
            if self.width == -1 {
                self.width = width_to_use;
            }
        }
        // The DrawMols are sized according to the panel height, which is all
        // there is to go on initially, so they may be wider than the width to
        // use in total. If so, shrink them to fit. Because the shrinking
        // imposes min and max font sizes, we may not get the smaller size we
        // want first go, so iterate until we do or we give up.
        let mut tot_width = reaction_width(&self.draw_mols, c, padding, ARROW_MULT, *plus_width);
        for _ in 0..5 {
            let Some(max_width) = self.draw_mols.iter().map(|dm| dm.width).reduce(i32::max) else {
                break;
            };
            *plus_width = max_width / 4;
            *plus_width = (*plus_width).min(width_to_use / 20);
            let old_tot_width = tot_width;
            let stretch = f64::from(width_to_use) * (1.0 - 2.0 * padding) / f64::from(tot_width);
            // If stretch < 1, the DrawMols need to shrink to fit; otherwise
            // they only get stretched along the panel and already fit for
            // height.
            if stretch >= 1.0 {
                break;
            }
            for i in c
                .reagents
                .clone()
                .chain(c.agents.clone())
                .chain(c.products.clone())
            {
                let (scale, font_scale) = (self.draw_mols[i].scale, self.draw_mols[i].font_scale);
                self.draw_mols[i].set_scale(
                    stretch * scale,
                    stretch * font_scale,
                    false,
                    &mut self.text,
                )?;
                self.draw_mols[i].shrink_to_fit(false, &mut self.text)?;
            }
            tot_width = reaction_width(&self.draw_mols, c, padding, ARROW_MULT, *plus_width);
            if f64::from((tot_width - old_tot_width).abs()) < 0.01 * f64::from(width_to_use) {
                break;
            }
        }

        // make sure plus_width remains 1A, in pixels
        if !c.reagents.is_empty() {
            *plus_width = self.draw_mols[c.reagents.start].scale as i32;
            tot_width = reaction_width(&self.draw_mols, c, padding, ARROW_MULT, *plus_width);
        } else if !c.products.is_empty() {
            *plus_width = self.draw_mols[c.products.start].scale as i32;
            tot_width = reaction_width(&self.draw_mols, c, padding, ARROW_MULT, *plus_width);
        }
        // if there's space, we can afford make the extras a bit bigger
        let mut num_gaps = if c.reagents.is_empty() {
            0
        } else {
            c.reagents.len() as i32 - 1
        };
        num_gaps += if c.products.is_empty() {
            0
        } else {
            c.products.len() as i32 - 1
        };
        num_gaps += 1; // for either side of the arrow
        num_gaps += if c.agents.is_empty() {
            2
        } else {
            c.agents.len() as i32 - 1
        };
        if width_to_use - tot_width > num_gaps * 5 {
            *plus_width += 5;
        }
        tot_width = reaction_width(&self.draw_mols, c, padding, ARROW_MULT, *plus_width);
        // And finally work out where to put all the pieces, centring them.
        // The padding is already taken care of, so take it out.
        let mut x_offset = self.x_offset
            + ((f64::from(width_to_use) - f64::from(tot_width) / (1.0 + 2.0 * padding)) / 2.0)
                .round() as i32;
        let y_offset = f64::from(self.y_offset);
        let panel_height = f64::from(self.panel_height);
        let mut offsets = Vec::new();
        for (n, i) in c.reagents.clone().enumerate() {
            let h_offset = y_offset + (panel_height - f64::from(self.draw_mols[i].height)) / 2.0;
            offsets.push(Point2D::new(f64::from(x_offset), h_offset));
            x_offset += self.draw_mols[i].width;
            if n + 1 < c.reagents.len() {
                x_offset += *plus_width;
            }
        }
        // only half a plus width to the arrow
        if !c.reagents.is_empty() {
            x_offset += *plus_width / 2;
        }
        let arrow_beg = Point2D::new(f64::from(x_offset), y_offset + panel_height / 2.0);
        let arrow_end = if c.agents.is_empty() {
            Point2D::new(
                arrow_beg.x + f64::from(ARROW_MULT * *plus_width),
                arrow_beg.y,
            )
        } else {
            x_offset += *plus_width / 2;
            for i in c.agents.clone() {
                offsets.push(Point2D::new(
                    f64::from(x_offset),
                    y_offset + 0.45 * panel_height - f64::from(self.draw_mols[i].height),
                ));
                x_offset += self.draw_mols[i].width + *plus_width / 2;
            }
            // the overlap at the end of the arrow has already been added
            Point2D::new(f64::from(x_offset), arrow_beg.y)
        };
        // arrow_end.x is an int upstream
        x_offset = arrow_end.x as i32 + *plus_width / 2;
        for i in c.products.clone() {
            let h_offset = y_offset + (panel_height - f64::from(self.draw_mols[i].height)) / 2.0;
            offsets.push(Point2D::new(f64::from(x_offset), h_offset));
            x_offset += self.draw_mols[i].width + *plus_width;
        }
        Ok((offsets, arrow_beg, arrow_end))
    }

    fn draw_reaction_part(
        &mut self,
        part: &std::ops::Range<usize>,
        plus_width: i32,
        init_offset: usize,
        offsets: &[Point2D],
    ) -> Result<usize> {
        if part.is_empty() {
            return Ok(init_offset);
        }
        let mut init_offset = init_offset;
        let mut plus_pos = Point2D::new(
            0.0,
            f64::from(self.y_offset) + f64::from(self.panel_height) / 2.0,
        );
        let component_padding = self.canvas.options.component_padding;
        for (n, i) in part.clone().enumerate() {
            // Adjust the scale to take account of any reagent padding. The
            // molecules will be re-centred in their panel, so don't need to
            // fiddle with offsets.
            let (w, h, scale) = {
                let dm = &self.draw_mols[i];
                (f64::from(dm.width), f64::from(dm.height), dm.scale)
            };
            let new_x_scale = (w - 2.0 * (w * component_padding)) / w;
            let new_y_scale = (h - 2.0 * (h * component_padding)) / h;
            let new_scale = new_x_scale.min(new_y_scale) * scale;
            self.draw_mols[i].set_scale(new_scale, new_scale, true, &mut self.text)?;
            self.active_mol_idx += 1;
            let off = offsets[init_offset];
            self.draw_mols[i].set_offsets(off.x, off.y, &mut self.text)?;
            self.canvas.scale = self.draw_mols[i].scale;
            self.draw_mols[i].draw(&mut self.canvas, &mut self.text)?;
            if plus_width != 0 && n + 1 < part.len() {
                plus_pos.x =
                    (off.x + f64::from(self.draw_mols[i].width) + offsets[init_offset + 1].x) / 2.0;
                let plus_stroke = if plus_width > 30 { 10 } else { plus_width / 3 };
                let col = self.canvas.options.symbol_colour;
                self.draw_plus(plus_pos, plus_stroke, col, true)?;
            }
            init_offset += 1;
        }
        Ok(init_offset)
    }

    /// The colours of the atoms by atom map number, when highlighting by
    /// reactant: atoms on the reactant side with an atom map number.
    fn find_reaction_highlights(
        &self,
        rxn: &ChemicalReaction,
        highlight_by_reactant: bool,
        highlight_colours: Option<&[DrawColour]>,
    ) -> BTreeMap<i32, DrawColour> {
        let mut atom_colours = BTreeMap::new();
        if !highlight_by_reactant {
            return atom_colours;
        }
        let colours: &[DrawColour] =
            highlight_colours.unwrap_or(&self.canvas.options.highlight_colour_palette);
        if colours.is_empty() {
            return atom_colours;
        }
        for (midx, frag) in rxn.reactants().iter().enumerate() {
            // upstream marks every atom of the reactants with the reaction role
            // 1 first (`ChemicalReactionToRxnMol`), so only the atom map matters
            for ai in 0..frag.num_atoms() {
                let map_num = frag.atom(ai).atom_map_num();
                if map_num != 0 {
                    atom_colours
                        .entry(map_num)
                        .or_insert(colours[midx % colours.len()]);
                }
            }
        }
        atom_colours
    }
}
