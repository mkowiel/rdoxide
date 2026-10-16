//! `SynthonShapeInput.cpp`: the shapes of a synthon (or query fragment) with
//! the volumes of its dummy atoms.

use crate::graphmol::gaussian_shape::{
    ShapeInput, ShapeInputOptions, ShapeOverlayOptions, Tokens, calc_shape_vol_and_grads,
    fmt_double,
};
use crate::graphmol::romol::ROMol;
use crate::{RDKitError, Result};

#[derive(Debug, Clone)]
pub struct SynthonShapeInput {
    shapes: ShapeInput,
    /// The volumes of the dummy atoms in each of the shapes
    dummy_volumes: Vec<f64>,
    /// The dummy atoms and their neighbours. On rare occasions a dummy may
    /// have 2 neighbours, so the first numbers may not be unique.
    dummy_atoms_and_nbrs: Vec<(u32, u32)>,
}

impl SynthonShapeInput {
    pub fn new(
        mol: &ROMol,
        conf_id: i32,
        opts: &ShapeInputOptions,
        overlay_opts: &ShapeOverlayOptions,
    ) -> Result<Self> {
        let shapes = ShapeInput::new(mol, conf_id, opts, overlay_opts)?;
        let mut s = SynthonShapeInput {
            shapes,
            dummy_volumes: Vec::new(),
            dummy_atoms_and_nbrs: Vec::new(),
        };
        s.build_dummy_atoms_and_nbrs()?;
        s.calculate_dummy_vols(overlay_opts)?;
        Ok(s)
    }

    /// Merges the other, assuming it has the correct number of atoms etc.
    pub fn merge(&mut self, mut other: SynthonShapeInput) {
        self.shapes.merge(&mut other.shapes);
        self.dummy_volumes.append(&mut other.dummy_volumes);
        // The dummy atoms and neighbours are properties of the synthon, not
        // of the individual shapes.
    }

    pub fn shapes(&self) -> &ShapeInput {
        &self.shapes
    }
    pub fn shapes_mut(&mut self) -> &mut ShapeInput {
        &mut self.shapes
    }

    pub fn dummy_volume(&self, shape_num: usize) -> Result<f64> {
        if shape_num >= self.shapes.num_shapes() {
            return Err(RDKitError::precondition(format!(
                "Invalid shape number ({shape_num} vs {}).",
                self.shapes.num_shapes()
            )));
        }
        Ok(self.dummy_volumes[shape_num])
    }

    pub fn num_dummy_atoms(&self) -> usize {
        let mut d: Vec<u32> = self.dummy_atoms_and_nbrs.iter().map(|p| p.0).collect();
        d.sort_unstable();
        d.dedup();
        d.len()
    }

    pub fn num_dummy_atom_nbors(&self) -> usize {
        self.dummy_atoms_and_nbrs.len()
    }

    pub fn dummy_atoms_and_nbrs(&self) -> &[(u32, u32)] {
        &self.dummy_atoms_and_nbrs
    }

    fn build_dummy_atoms_and_nbrs(&mut self) -> Result<()> {
        let tmp = self.shapes.shape_to_mol(false, true)?;
        for atom in tmp.atoms() {
            if atom.atomic_num() == 0 {
                for nbr in tmp.atom_neighbors(atom.idx()) {
                    self.dummy_atoms_and_nbrs.push((atom.idx(), nbr));
                }
            }
        }
        Ok(())
    }

    fn calculate_dummy_vols(&mut self, opts: &ShapeOverlayOptions) -> Result<()> {
        // Set the alphas for the dummy atoms to -1.0 so they are skipped in
        // the volume calculation.
        let mut alphas = self.shapes.alphas().to_vec();
        for p in &self.dummy_atoms_and_nbrs {
            alphas[p.0 as usize] = -1.0;
        }
        let na = self.shapes.num_atoms();
        let mut conv = vec![0.0; 12 * (na + self.shapes.num_features())];
        let cutoff2 = opts.dist_cutoff * opts.dist_cutoff;
        self.dummy_volumes.clear();
        for i in 0..self.shapes.num_shapes() {
            self.shapes.set_active_shape(i)?;
            let full = self.shapes.shape_volume();
            let carbon = self.shapes.carbon_radii();
            let coords = self.shapes.coords();
            let part = calc_shape_vol_and_grads(
                coords,
                &alphas,
                na,
                carbon,
                coords,
                &alphas,
                na,
                carbon,
                &mut conv,
                opts.use_dist_cutoff,
                cutoff2,
                None,
            );
            self.dummy_volumes.push(full - part);
        }
        Ok(())
    }

    /// The Boost text archive form, as stored in the database files.
    pub fn to_text(&self) -> String {
        let inner = self.shapes.to_text();
        let mut out = String::from("22 serialization::archive 19 0 0");
        out.push_str(&format!(" {} {inner}", inner.len()));
        out.push_str(&format!(" {} 0", self.dummy_volumes.len()));
        for v in &self.dummy_volumes {
            out.push(' ');
            out.push_str(&fmt_double(*v));
        }
        out.push_str(&format!(" 0 0 {} 0", self.dummy_atoms_and_nbrs.len()));
        if !self.dummy_atoms_and_nbrs.is_empty() {
            out.push_str(" 0 0");
            for (a, b) in &self.dummy_atoms_and_nbrs {
                out.push_str(&format!(" {a} {b}"));
            }
        }
        out
    }

    pub fn from_text(text: &str) -> Result<Self> {
        let mut t = Tokens::new(text);
        t.expect_header()?;
        t.class_info()?;
        let inner = t.string()?;
        let shapes = ShapeInput::from_text(&inner)?;
        let dummy_volumes = t.doubles()?;
        t.class_info()?;
        let n = t.usize()?;
        t.token()?;
        let mut pairs = Vec::with_capacity(n);
        if n > 0 {
            t.class_info()?;
        }
        for _ in 0..n {
            pairs.push((t.usize()? as u32, t.usize()? as u32));
        }
        Ok(SynthonShapeInput {
            shapes,
            dummy_volumes,
            dummy_atoms_and_nbrs: pairs,
        })
    }
}
