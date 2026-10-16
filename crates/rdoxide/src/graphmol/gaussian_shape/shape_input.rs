//! `ShapeInput`: the Gaussian representation of a molecule (one shape per
//! conformer), with the colour features and the canonical orientation.

use std::sync::OnceLock;

use super::alignment::{ShapeView, calc_color_vol_and_grads, calc_shape_vol_and_grads};
use super::options::ShapeOverlayOptions;
use super::overlay::align_shape;
use crate::geometry::{Point3D, Transform3D};
use crate::graphmol::atom::Atom;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::find_rings::find_sssr;
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, remove_hs};
use crate::graphmol::mol_ops::structure::renumber_atoms;
use crate::graphmol::mol_transforms::{
    compute_canonical_transform, compute_principal_axes_and_moments,
};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::RDValue;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::sim_div_pickers::leader;
use crate::{RDKitError, Result, rd_log_warning};

pub const CARBON_RAD: f64 = 1.70;
/// same as Xe
pub const DUMMY_RAD: f64 = 2.16;
pub const KAPPA: f64 = 2.41798793102;
/// same radius for all feature/colour "atoms", as used by the PubChem code
const RADIUS_COLOR: f64 = 1.08265;

/// The vdW radius of an element, for those that have one.
fn vdw_radius(z: u32) -> Option<f64> {
    Some(match z {
        0 => DUMMY_RAD, // Dummy, same as Xe.
        1 => 1.10,
        2 => 1.40,
        3 => 1.81,
        4 => 1.53,
        5 => 1.92,
        6 => CARBON_RAD,
        7 => 1.55,
        8 => 1.52,
        9 => 1.47,
        10 => 1.54,
        11 => 2.27,
        12 => 1.73,
        13 => 1.84,
        14 => 2.10,
        15 => 1.80,
        16 => 1.80,
        17 => 1.75,
        18 => 1.88,
        19 => 2.75,
        20 => 2.31,
        31 => 1.87,
        32 => 2.11,
        33 => 1.85,
        34 => 1.90,
        35 => 1.83,
        36 => 2.02,
        37 => 3.03,
        38 => 2.49,
        49 => 1.93,
        50 => 2.17,
        51 => 2.06,
        52 => 2.06,
        53 => 1.98,
        54 => 2.16,
        55 => 3.43,
        56 => 2.68,
        81 => 1.96,
        82 => 2.02,
        83 => 2.07,
        84 => 1.97,
        85 => 2.02,
        86 => 2.20,
        87 => 3.48,
        88 => 2.83,
        _ => return None,
    })
}

fn standard_atom_radius(z: u32) -> Result<f64> {
    // Mostly they will be carbons, so don't look up.
    if z == 6 {
        return Ok(CARBON_RAD);
    }
    vdw_radius(z).ok_or_else(|| RDKitError::Value(format!("No VdW radius for atom with Z={z}")))
}

/// A colour feature supplied by the caller.
#[derive(Debug, Clone, PartialEq)]
pub struct CustomFeature {
    pub kind: u32,
    pub pos: Point3D,
    pub rad: f64,
    /// The atoms that the feature was derived from. May be empty.
    pub atoms: Vec<u32>,
}

/// `ShapeInputOptions`
#[derive(Debug, Clone, PartialEq)]
pub struct ShapeInputOptions {
    /// Whether to build the colour features (by default with the RDKit
    /// pharmacophore definitions).
    pub use_colors: bool,
    /// Custom colour features used verbatim, one inner vector for each
    /// conformation in the molecule.
    pub custom_features: Vec<Vec<CustomFeature>>,
    /// If not empty, use just these atoms of the molecule.
    pub atom_subset: Vec<u32>,
    /// Non-standard radii for atoms (by atom index in the molecule); atoms not
    /// listed keep the standard radius.
    pub atom_radii: Vec<(u32, f64)>,
    /// Whether to use carbon radii for all atoms (quicker but less accurate)
    /// or the vdW radii of the elements.
    pub all_carbon_radii: bool,
    /// If there is more than 1 conformer, prune the shapes so that none of
    /// them is more similar to another than this. Negative: no pruning.
    pub shape_prune_threshold: f64,
    /// Sort the shapes in descending order of total volume.
    pub sort_shapes: bool,
    /// Whether to include dummy atoms in the shape.
    pub include_dummies: bool,
}

impl Default for ShapeInputOptions {
    fn default() -> Self {
        ShapeInputOptions {
            use_colors: true,
            custom_features: Vec::new(),
            atom_subset: Vec::new(),
            atom_radii: Vec::new(),
            all_carbon_radii: true,
            shape_prune_threshold: -1.0,
            sort_shapes: true,
            include_dummies: true,
        }
    }
}

/// The shapes responsible for the best similarity of `best_similarity` and
/// the transformation that gave it.
#[derive(Debug, Clone)]
pub struct ShapePair {
    pub this_shape: usize,
    pub fit_shape: usize,
    pub xform: Transform3D,
}

/// The Gaussian shapes of the conformers of a molecule.
#[derive(Debug, Clone, Default)]
pub struct ShapeInput {
    pub(super) active_shape: usize,
    /// The coordinates of the atoms and features, 3 values for each, for
    /// each shape.
    pub(super) coords: Vec<Vec<f64>>,
    /// `KAPPA / (r * r)` for the atoms and features. Negative values are
    /// switched off.
    pub(super) alphas: Vec<f64>,
    /// The feature types, 0 for atoms.
    pub(super) types: Vec<i32>,
    pub(super) num_atoms: usize,
    pub(super) num_feats: usize,
    pub(super) self_overlap_shape_vols: Vec<f64>,
    pub(super) self_overlap_color_vols: Vec<f64>,
    /// The points at the extremes of x, y and z: min x, y, z then max.
    pub(super) extreme_points: Vec<[usize; 6]>,
    /// Flags those atoms with a carbon radius, for faster calculation later.
    pub(super) carbon_radii: Option<Vec<bool>>,
    pub(super) smiles: String,
    /// Whether the shape has been put in its canonical orientation.
    pub(super) normalizeds: Vec<bool>,
    /// If the shape is moved, the normalization is no longer valid.
    pub(super) normalization_oks: Vec<bool>,
    pub(super) canon_rots: Vec<[f64; 9]>,
    pub(super) canon_transs: Vec<[f64; 3]>,
    /// The sorted eigenvalues of the principal axes.
    pub(super) eigen_values: Vec<[f64; 3]>,
}

/// The atoms of `mol` to keep, with aromatic flags fixed on those that are no
/// longer in a ring.
fn extract_subset(mol: &RWMol, atoms_to_keep: &[u32]) -> Result<RWMol> {
    let mut ret = mol.clone();
    let mut keep = vec![false; mol.num_atoms() as usize];
    for &a in atoms_to_keep {
        keep[a as usize] = true;
    }
    ret.begin_batch_edit()?;
    for idx in 0..ret.num_atoms() {
        if !keep[idx as usize] {
            ret.remove_atom(idx)?;
        }
    }
    ret.commit_batch_edit()?;

    // Fix the aromatic flags on anything that is no longer in a ring.
    find_sssr(&mut ret, false, false)?;
    for idx in 0..ret.num_atoms() {
        if ret.atom(idx).is_aromatic() && ret.ring_info().num_atom_rings(idx)? == 0 {
            ret.atom_mut(idx).set_is_aromatic(false);
        }
    }
    Ok(ret)
}

impl ShapeInput {
    /// Creates the shapes of `mol`: for the conformer `conf_id`, or all of
    /// them if it is negative. The distance cutoff of `overlay_opts` is used in
    /// the self-overlap calculations.
    pub fn new(
        mol: &ROMol,
        conf_id: i32,
        opts: &ShapeInputOptions,
        overlay_opts: &ShapeOverlayOptions,
    ) -> Result<Self> {
        if mol.num_conformers() == 0 {
            let mut name = mol.name();
            let mut tmp = mol.clone();
            let smi = mol_to_smiles(&mut tmp, &SmilesWriteParams::default())?;
            name.push_str("  ");
            name.push_str(&smi);
            return Err(RDKitError::precondition(format!(
                "ShapeInput object needs the molecule to have conformers.  {name}"
            )));
        }
        // Subsetting the molecule makes any bespoke atom radii, identified by
        // atom index, incorrect so stash them as atom properties.
        let mut work = mol.clone();
        for &(idx, radius) in &opts.atom_radii {
            work.atom_mut(idx)
                .set_prop("BespokeRadius", radius, false)?;
        }
        let mut tmp_mol = if opts.atom_subset.is_empty() {
            work
        } else {
            extract_subset(&work, &opts.atom_subset)?
        };
        let mut res = ShapeInput {
            smiles: mol_to_smiles(&mut tmp_mol, &SmilesWriteParams::default())?,
            ..Default::default()
        };
        let at_order: Vec<u32> = match tmp_mol.get_prop(cp::_smilesAtomOutputOrder) {
            Some(RDValue::VecUnsignedInt(v)) => v.clone(),
            Some(RDValue::VecInt(v)) => v.iter().map(|&i| i as u32).collect(),
            _ => return Err(RDKitError::Key(cp::_smilesAtomOutputOrder.into())),
        };
        let mut tmp_mol = renumber_atoms(&tmp_mol, &at_order)?;

        if !tmp_mol.ring_info().is_initialized() {
            // Query molecules don't seem to have the ring info generated on
            // creation.
            find_sssr(&mut tmp_mol, false, false)?;
        }
        let n_confs = tmp_mol.num_conformers() as usize;
        res.normalizeds = vec![false; n_confs];
        res.normalization_oks = vec![false; n_confs];

        // Molecules don't necessarily have consecutive conformer IDs from 0
        // upwards. conf_num is just the sequential number of the iterated
        // conformers.
        let process_conf = |res: &mut ShapeInput, conf: &Conformer, conf_num: usize, fa: bool| {
            res.extract_atoms(&tmp_mol, conf, opts, fa)?;
            if opts.use_colors {
                res.extract_features(&tmp_mol, conf, conf_num, opts, fa)?;
            }
            Ok::<(), RDKitError>(())
        };
        if conf_id >= 0 {
            process_conf(&mut res, tmp_mol.conformer(conf_id)?, 0, true)?;
        } else {
            for (i, conf) in tmp_mol.conformers().iter().enumerate() {
                process_conf(&mut res, conf, i, i == 0)?;
            }
        }
        res.active_shape = 0;
        res.calc_normalization()?;
        res.calc_extremes();
        res.calculate_self_overlaps(overlay_opts);
        if opts.shape_prune_threshold > 0.0 && tmp_mol.num_conformers() > 1 {
            res.prune_shapes(opts.shape_prune_threshold)?;
        }
        if opts.sort_shapes {
            res.sort_shapes_by_volumes();
        }
        res.active_shape = 0;
        Ok(res)
    }

    /// A new `ShapeInput` with the single shape number `shape_num` of `other`.
    pub fn from_shape(other: &ShapeInput, shape_num: usize) -> Result<Self> {
        if shape_num >= other.num_shapes() {
            return Err(RDKitError::precondition(format!(
                "Invalid shape number in makeSingleShape : {} vs {}.",
                shape_num,
                other.num_shapes()
            )));
        }
        Ok(ShapeInput {
            active_shape: 0,
            coords: vec![other.coords[shape_num].clone()],
            alphas: other.alphas.clone(),
            types: other.types.clone(),
            num_atoms: other.num_atoms,
            num_feats: other.num_feats,
            self_overlap_shape_vols: vec![other.self_overlap_shape_vols[shape_num]],
            self_overlap_color_vols: vec![other.self_overlap_color_vols[shape_num]],
            extreme_points: vec![other.extreme_points[shape_num]],
            carbon_radii: other.carbon_radii.clone(),
            smiles: other.smiles.clone(),
            normalizeds: vec![other.normalizeds[shape_num]],
            normalization_oks: vec![other.normalization_oks[shape_num]],
            canon_rots: vec![other.canon_rots[shape_num]],
            canon_transs: vec![other.canon_transs[shape_num]],
            eigen_values: vec![other.eigen_values[shape_num]],
        })
    }

    /// Merges the shapes of `other` into this one, assuming it has the
    /// correct number of atoms etc. Empties `other`, unless they can't be
    /// merged in which case it is left as it was. The merge can only be done
    /// if `other` has the same number of coordinates per conformer, and the
    /// feature types of the two match.
    pub fn merge(&mut self, other: &mut ShapeInput) {
        if let (Some(a), Some(b)) = (self.coords.first(), other.coords.first())
            && a.len() != b.len()
        {
            rd_log_warning!("Can't merge shapes as different sizes.");
            return;
        }
        if self.types != other.types {
            rd_log_warning!("Can't merge shapes as different types.");
            return;
        }
        if other.coords.is_empty() {
            return;
        }
        self.coords.append(&mut other.coords);
        self.self_overlap_shape_vols
            .append(&mut other.self_overlap_shape_vols);
        self.self_overlap_color_vols
            .append(&mut other.self_overlap_color_vols);
        self.extreme_points.append(&mut other.extreme_points);
        self.canon_rots.append(&mut other.canon_rots);
        self.canon_transs.append(&mut other.canon_transs);
        self.eigen_values.append(&mut other.eigen_values);
        self.normalizeds.append(&mut other.normalizeds);
        self.normalization_oks.append(&mut other.normalization_oks);
    }

    pub fn smiles(&self) -> &str {
        &self.smiles
    }

    pub fn active_shape(&self) -> usize {
        self.active_shape
    }

    /// Sets the conformation used for future calculations.
    pub fn set_active_shape(&mut self, new_shape: usize) -> Result<()> {
        if new_shape >= self.coords.len() {
            return Err(RDKitError::precondition(format!(
                "Invalid shape number ({} vs {}).",
                new_shape,
                self.coords.len()
            )));
        }
        self.active_shape = new_shape;
        Ok(())
    }

    /// The coordinates of the active shape, 3 values for each atom and
    /// feature.
    pub fn coords(&self) -> &[f64] {
        &self.coords[self.active_shape]
    }

    /// The alpha values for the atoms and colour features.
    pub fn alphas(&self) -> &[f64] {
        &self.alphas
    }

    /// Multiplies the alpha value for the given atom/feature by -1.0, which
    /// toggles whether it takes part in the volume calculation.
    pub fn negate_alpha(&mut self, alpha_num: usize) -> Result<()> {
        if alpha_num >= self.alphas.len() {
            return Err(RDKitError::precondition("Bad number for alpha negation"));
        }
        self.alphas[alpha_num] *= -1.0;
        Ok(())
    }

    /// The coordinates of the atoms and optionally the features.
    pub fn atom_points(&self, include_colors: bool) -> Vec<Point3D> {
        let mut n = self.num_atoms;
        if include_colors {
            n += self.num_feats;
        }
        self.coords[self.active_shape]
            .chunks_exact(3)
            .take(n)
            .map(|c| Point3D::new(c[0], c[1], c[2]))
            .collect()
    }

    /// Whether the active shape is in its canonical orientation.
    pub fn is_normalized(&self) -> bool {
        self.normalizeds[self.active_shape]
    }

    /// The feature types of all atoms and features, 0 for atoms.
    pub fn feature_types(&self) -> &[i32] {
        &self.types
    }

    pub fn num_atoms(&self) -> usize {
        self.num_atoms
    }

    pub fn num_features(&self) -> usize {
        self.num_feats
    }

    /// The number of shapes, which may be smaller than the number of
    /// conformations of the input molecule if pruning was done.
    pub fn num_shapes(&self) -> usize {
        self.coords.len()
    }

    /// The volume of the atoms of the active shape.
    pub fn shape_volume(&self) -> f64 {
        self.self_overlap_shape_vols[self.active_shape]
    }

    pub fn shape_volume_of(&self, shape_num: usize) -> Result<f64> {
        self.check_shape_num(shape_num)?;
        Ok(self.self_overlap_shape_vols[shape_num])
    }

    /// The volume of the colour features of the active shape.
    pub fn color_volume(&self) -> f64 {
        self.self_overlap_color_vols[self.active_shape]
    }

    pub fn color_volume_of(&self, shape_num: usize) -> Result<f64> {
        self.check_shape_num(shape_num)?;
        Ok(self.self_overlap_color_vols[shape_num])
    }

    fn check_shape_num(&self, shape_num: usize) -> Result<()> {
        if shape_num >= self.coords.len() {
            return Err(RDKitError::precondition(format!(
                "Invalid shape number ({} vs {}).",
                shape_num,
                self.coords.len()
            )));
        }
        Ok(())
    }

    /// Flags for the atoms that have a carbon radius, if they are tracked.
    pub fn carbon_radii(&self) -> Option<&[bool]> {
        self.carbon_radii.as_deref()
    }

    /// What the alignment needs of the active shape.
    pub fn view(&self) -> ShapeView<'_> {
        ShapeView {
            coords: self.coords(),
            alphas: &self.alphas,
            types: &self.types,
            carbon_radii: self.carbon_radii(),
            num_shape: self.num_atoms,
            num_color: self.num_feats,
            shape_vol: self.shape_volume(),
            color_vol: self.color_volume(),
        }
    }

    /// The canonical rotation of the active shape.
    pub fn calc_canonical_rotation(&mut self) -> Result<[f64; 9]> {
        self.ensure_normalization()?;
        Ok(self.canon_rots[self.active_shape])
    }

    /// The canonical translation of the active shape.
    pub fn calc_canonical_translation(&mut self) -> Result<[f64; 3]> {
        self.ensure_normalization()?;
        Ok(self.canon_transs[self.active_shape])
    }

    /// The eigenvalues of the coordinates matrix.
    pub fn calc_eigen_values(&mut self) -> Result<[f64; 3]> {
        self.ensure_normalization()?;
        Ok(self.eigen_values[self.active_shape])
    }

    fn ensure_normalization(&mut self) -> Result<()> {
        if !self.normalization_oks[self.active_shape] {
            self.calc_normalization()?;
        }
        Ok(())
    }

    /// The numbers of the points at the extremes of x, y and z of the active
    /// shape: minimum x, y, z then the maxima.
    pub fn calc_extremes(&mut self) -> [usize; 6] {
        if !self.normalization_oks[self.active_shape] {
            self.calculate_extremes();
        }
        self.extreme_points[self.active_shape]
    }

    /// The principal moments of inertia of the active shape, with the colour
    /// features (as xenon atoms) if requested.
    pub fn calc_moments_of_inertia(&self, include_colors: bool) -> Result<[f64; 3]> {
        let mut tmp = self.shape_to_mol(include_colors, true)?;
        let axes = compute_principal_axes_and_moments(&mut tmp, -1, true, false, None)?
            .ok_or_else(|| RDKitError::Runtime("eigenvalue calculation did not converge".into()))?;
        Ok(axes.moments)
    }

    /// Aligns the principal axes to the cartesian axes and centres the active
    /// shape on the origin.
    pub fn normalize_coords(&mut self) -> Result<()> {
        let a = self.active_shape;
        if self.normalizeds[a] {
            return Ok(());
        }
        self.ensure_normalization()?;
        let mut canon_rot = Transform3D::identity();
        for i in 0..3 {
            for j in 0..3 {
                canon_rot.set_val(i, j, self.canon_rots[a][3 * i + j]);
            }
        }
        let t = self.canon_transs[a];
        let trans = canon_rot.transform_point(&Point3D::new(t[0], t[1], t[2]));
        canon_rot.set_translation(&trans);
        self.transform_coords(&canon_rot);
        self.normalizeds[a] = true;
        // Recalculate the extremes now we've changed the coordinates.
        self.calc_extremes();
        Ok(())
    }

    /// Applies the transformation to the active shape.
    pub fn transform_coords(&mut self, xform: &Transform3D) {
        let a = self.active_shape;
        apply_transform_to_shape(&mut self.coords[a], xform);
        self.normalizeds[a] = false;
        self.normalization_oks[a] = false;
    }

    /// A molecule made from the active shape. If requested the features are
    /// added as xenon atoms. Without bonds it is just carbon atoms, otherwise
    /// a full molecule built from the SMILES string.
    pub fn shape_to_mol(&self, include_colors: bool, with_bonds: bool) -> Result<RWMol> {
        let mut mol = if with_bonds {
            // The SMILES string and the atom coordinates are in the same
            // order.
            let params = SmilesParserParams {
                sanitize: false,
                ..Default::default()
            };
            let mut mol = mol_from_smiles(&self.smiles, &params)?
                .ok_or_else(|| RDKitError::Value("bad SMILES in the shape".into()))?;
            // Construction of the shape removes all Hs except isotopes, so do
            // that from the molecule too.
            let rhps = RemoveHsParameters {
                remove_degree_zero: true,
                remove_higher_degrees: true,
                remove_only_h_neighbors: true,
                remove_isotopes: false,
                remove_dummy_neighbors: true,
                remove_defining_bond_stereo: true,
                remove_with_wedged_bond: true,
                remove_with_query: true,
                remove_in_sgroups: true,
                show_warnings: false,
                remove_nonimplicit: true,
                remove_hydrides: true,
                remove_nontetrahedral_neighbors: true,
                ..Default::default()
            };
            remove_hs(&mut mol, &rhps, false)?;
            mol
        } else {
            let mut mol = RWMol::new();
            for _ in 0..self.num_atoms {
                mol.add_atom(Atom::new(6), true);
            }
            mol
        };
        if include_colors {
            for _ in 0..self.num_feats {
                mol.add_atom(Atom::new(54), true);
            }
        }
        let mut num = self.num_atoms;
        if include_colors {
            num += self.num_feats;
        }
        let mut conf = Conformer::new(num);
        let cds = self.coords();
        for (i, pos) in conf.positions_mut().iter_mut().enumerate() {
            *pos = Point3D::new(cds[3 * i], cds[3 * i + 1], cds[3 * i + 2]);
        }
        mol.add_conformer(conf, true)?;
        Ok(mol)
    }

    /// The best similarity score between all shapes of this one and the
    /// other, stopping as soon as one above `threshold` is found. The scores
    /// are `[-1, -1, -1]` if nothing is above the threshold. Also the numbers
    /// of the two shapes that gave the best score and the transformation.
    pub fn best_similarity(
        &mut self,
        fit_shape: &ShapeInput,
        threshold: f64,
        overlay_opts: &ShapeOverlayOptions,
    ) -> Result<([f64; 3], Option<ShapePair>)> {
        let mut best_sim = [-1.0; 3];
        if self.max_possible_similarity(fit_shape, &ShapeOverlayOptions::default())? < threshold {
            return Ok((best_sim, None));
        }
        let mut best = None;
        let curr_active = self.active_shape;
        'outer: for i in 0..self.num_shapes() {
            self.set_active_shape(i)?;
            for j in 0..fit_shape.num_shapes() {
                let max_sim = max_score(
                    self.shape_volume(),
                    fit_shape.shape_volume_of(j)?,
                    self.color_volume(),
                    fit_shape.color_volume_of(j)?,
                    overlay_opts,
                );
                if max_sim > threshold {
                    let mut single_fit = ShapeInput::from_shape(fit_shape, j)?;
                    let (scores, xform) = align_shape(self, &mut single_fit, overlay_opts)?;
                    if scores[0] > best_sim[0] {
                        best_sim = scores;
                        best = Some(ShapePair {
                            this_shape: i,
                            fit_shape: j,
                            xform,
                        });
                    }
                }
                // Floating point cruft means we sometimes get a similarity
                // slightly above 1.0. 1.0 is the maximum possible, so stop if
                // we hit it.
                if best_sim[0] > 1.0 || (best_sim[0] - 1.0).abs() < 1.0e-6 {
                    break 'outer;
                }
            }
        }
        self.set_active_shape(curr_active)?;
        Ok((best_sim, best))
    }

    /// The maximum similarity achievable between the 2 shapes, when one shape
    /// is entirely inside the other: the upper bound of what is achievable.
    pub fn max_possible_similarity(
        &self,
        fit_shape: &ShapeInput,
        overlay_opts: &ShapeOverlayOptions,
    ) -> Result<f64> {
        let mut max_sim = 0.0;
        for i in 0..self.num_shapes() {
            for j in 0..fit_shape.num_shapes() {
                let sim = max_score(
                    fit_shape.self_overlap_shape_vols[j],
                    self.self_overlap_shape_vols[i],
                    fit_shape.self_overlap_color_vols[j],
                    self.self_overlap_color_vols[i],
                    overlay_opts,
                );
                if sim > max_sim {
                    max_sim = sim;
                }
            }
        }
        Ok(max_sim)
    }

    /// Prunes the shapes so none is more similar to another than the
    /// threshold.
    pub fn prune_shapes(&mut self, sim_threshold: f64) -> Result<()> {
        if self.coords.len() < 2 || sim_threshold < 0.0 {
            return Ok(());
        }
        let shapes = self.clone();
        // The picker works on distances.
        let dist = |i: u32, j: u32| -> Result<f64> {
            let mut a = shapes.clone();
            let mut b = shapes.clone();
            a.set_active_shape(i as usize)?;
            b.set_active_shape(j as usize)?;
            let (scores, _) = align_shape(&a, &mut b, &ShapeOverlayOptions::default())?;
            Ok(1.0 - scores[0])
        };
        let mut picks = leader::lazy_pick(
            dist,
            self.coords.len() as u32,
            0,
            &[],
            1.0 - sim_threshold,
            1,
        )?;
        // Allow for mysterious LeaderPicker behaviour where it returns a
        // vector of 0s when it should only pick 1 shape.
        picks.sort_unstable();
        picks.dedup();
        self.select_conformations(&picks);
        self.active_shape = 0;
        Ok(())
    }

    fn extract_atoms(
        &mut self,
        mol: &ROMol,
        conf: &Conformer,
        opts: &ShapeInputOptions,
        fill_alphas: bool,
    ) -> Result<()> {
        if !opts.all_carbon_radii {
            let n = if opts.atom_subset.is_empty() {
                conf.num_atoms() as usize
            } else {
                opts.atom_subset.len()
            };
            self.carbon_radii = Some(vec![false; n]);
        }
        let mut these = Vec::with_capacity(conf.num_atoms() as usize * 3);
        for atom in mol.atoms() {
            // Ignore H atoms except deuterium and tritium which are treated
            // elsewhere as explicit atoms, and do use dummies if requested.
            // The latter two are set to be ignored in the volume calculation,
            // however. This is just to keep the atom numberings correct.
            // Dummy atoms can also have isotope number 1.
            let z = atom.atomic_num();
            if z == 1 && atom.isotope() <= 1 {
                continue;
            }
            if !opts.include_dummies && z == 0 {
                continue;
            }
            let at_idx = atom.idx();
            let pos = conf.atom_pos(at_idx)?;
            these.extend([pos.x, pos.y, pos.z]);
            if !fill_alphas {
                continue;
            }
            let alpha = if opts.all_carbon_radii && self.carbon_radii.is_none() {
                KAPPA / (CARBON_RAD * CARBON_RAD)
            } else {
                let bespoke = atom.get_prop("BespokeRadius").and_then(RDValue::as_double);
                let rad = match bespoke {
                    Some(r) if !opts.atom_radii.is_empty() => r,
                    _ => {
                        if z == 6
                            && let Some(flags) = self.carbon_radii.as_mut()
                        {
                            flags[at_idx as usize] = true;
                        }
                        standard_atom_radius(z)?
                    }
                };
                // Check there hasn't been a logic error leaving a radius of
                // 0.0.
                if rad == 0.0 {
                    return Err(RDKitError::Value("Atom has radius 0.0.".into()));
                }
                KAPPA / (rad * rad)
            };
            // H atoms (isotopes) are not included in the volume calculations.
            self.alphas.push(if z == 1 { -1.0 } else { alpha });
        }
        self.active_shape = self.coords.len();
        self.coords.push(these);
        if fill_alphas {
            self.num_atoms = self.coords[0].len() / 3;
            self.types.resize(self.num_atoms, 0);
            self.num_feats = 0;
        }
        Ok(())
    }

    fn extract_features(
        &mut self,
        mol: &ROMol,
        conf: &Conformer,
        conf_num: usize,
        opts: &ShapeInputOptions,
        fill_alphas: bool,
    ) -> Result<()> {
        let feats = if opts.custom_features.is_empty() {
            let mut feats = Vec::new();
            find_features(mol, conf, &mut feats, None)?;
            feats
        } else {
            opts.custom_features
                .get(conf_num)
                .ok_or_else(|| {
                    RDKitError::precondition(format!(
                        "Conformer number {} too large for custom features size {}",
                        conf_num,
                        opts.custom_features.len()
                    ))
                })?
                .clone()
        };
        // Copy them directly except that the last is a radius, so convert to
        // alpha.
        let a = self.active_shape;
        for f in &feats {
            self.coords[a].extend([f.pos.x, f.pos.y, f.pos.z]);
            if fill_alphas {
                self.alphas.push(KAPPA / (f.rad * f.rad));
                self.types.push(f.kind as i32);
                self.num_feats += 1;
            }
        }
        Ok(())
    }

    /// The rotation and translation that align the principal axes to the
    /// cartesian axes and centre on the origin, for the active shape.
    fn calc_normalization(&mut self) -> Result<()> {
        let n = self.coords.len();
        if self.eigen_values.len() < n {
            self.eigen_values.resize(n, [0.0; 3]);
            self.canon_rots.resize(n, [0.0; 9]);
            self.canon_transs.resize(n, [0.0; 3]);
            self.extreme_points.resize(n, [0; 6]);
        }
        let a = self.active_shape;
        // Build a "molecule" just of the shape points, not the colour
        // features, with which to calculate the canonical transformation.
        // Doesn't ever use the input molecule in case the shape was built from
        // a subset of its atoms.
        let tmp = self.shape_to_mol(false, false)?;
        let (canon_xform, eigen_values) =
            compute_canonical_transform(&tmp, tmp.conformer(-1)?, None, false, true)?;
        self.eigen_values[a] = eigen_values;
        for i in 0..3 {
            for j in 0..3 {
                self.canon_rots[a][3 * i + j] = canon_xform.get_val(i, j);
            }
        }
        let mut trans = [0.0; 3];
        for c in self.coords[a].chunks_exact(3).take(self.num_atoms) {
            trans[0] -= c[0];
            trans[1] -= c[1];
            trans[2] -= c[2];
        }
        let n_atoms = self.num_atoms as f64;
        self.canon_transs[a] = [trans[0] / n_atoms, trans[1] / n_atoms, trans[2] / n_atoms];
        self.normalization_oks[a] = true;
        Ok(())
    }

    fn calculate_extremes(&mut self) {
        let n = self.coords.len();
        if self.extreme_points.len() < n {
            self.extreme_points.resize(n, [0; 6]);
        }
        let a = self.active_shape;
        let cds = &self.coords[a];
        let mut ext = [0usize; 6];
        for (j, c) in cds.chunks_exact(3).enumerate() {
            for k in 0..3 {
                if c[k] < cds[3 * ext[k] + k] {
                    ext[k] = j;
                }
                if c[k] > cds[3 * ext[k + 3] + k] {
                    ext[k + 3] = j;
                }
            }
        }
        self.extreme_points[a] = ext;
    }

    fn select_conformations(&mut self, picks: &[i32]) {
        fn pick<T: Clone>(v: &[T], picks: &[i32]) -> Vec<T> {
            picks.iter().map(|&p| v[p as usize].clone()).collect()
        }
        self.coords = pick(&self.coords, picks);
        self.self_overlap_shape_vols = pick(&self.self_overlap_shape_vols, picks);
        self.self_overlap_color_vols = pick(&self.self_overlap_color_vols, picks);
        self.extreme_points = pick(&self.extreme_points, picks);
        self.normalizeds = pick(&self.normalizeds, picks);
        self.normalization_oks = pick(&self.normalization_oks, picks);
        self.canon_rots = pick(&self.canon_rots, picks);
        self.canon_transs = pick(&self.canon_transs, picks);
        self.eigen_values = pick(&self.eigen_values, picks);
    }

    fn calculate_self_overlaps(&mut self, opts: &ShapeOverlayOptions) {
        let cutoff2 = opts.dist_cutoff * opts.dist_cutoff;
        let mut conv = vec![0.0; 12 * (self.num_atoms + self.num_feats)];
        let (na, nf) = (self.num_atoms, self.num_feats);
        let carbon = self.carbon_radii.as_deref();
        for cds in &self.coords {
            self.self_overlap_shape_vols.push(calc_shape_vol_and_grads(
                cds,
                &self.alphas,
                na,
                carbon,
                cds,
                &self.alphas,
                na,
                carbon,
                &mut conv,
                opts.use_dist_cutoff,
                cutoff2,
                None,
            ));
            self.self_overlap_color_vols.push(calc_color_vol_and_grads(
                &cds[3 * na..],
                &self.alphas[na..],
                nf,
                &self.types[na..],
                &cds[3 * na..],
                &self.alphas[na..],
                nf,
                &self.types[na..],
                na,
                &mut conv,
                opts.use_dist_cutoff,
                cutoff2,
                None,
            ));
        }
    }

    /// Sorts the shapes in descending order of the sum of the shape and colour
    /// volumes.
    fn sort_shapes_by_volumes(&mut self) {
        let mut vals: Vec<(f64, usize)> = (0..self.coords.len())
            .map(|i| {
                (
                    self.self_overlap_shape_vols[i] + self.self_overlap_color_vols[i],
                    i,
                )
            })
            .collect();
        libcxx_sort_by(&mut vals, |a, b| a.0 > b.0);
        let picks: Vec<i32> = vals.iter().map(|v| v.1 as i32).collect();
        self.select_conformations(&picks);
    }
}

/// The SMARTS of the pharmacophore features, one group per feature type; the
/// type is the (1-based) position of the group.
const SMARTS_PATTERNS: &[&[&str]] = &[
    // Donor
    &["[$([N;!H0;v3,v4&+1]),\
$([O,S;H1;+0]),\
n&H1&+0]"],
    // Acceptor
    &["[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),\
$([O,S;H0;v2]),\
$([O,S;-]),\
$([N;v3;!$(N-*=[O,N,P,S])]),\
n&H0&+0,\
$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]"],
    // rings
    &[
        "[r]1[r][r]1",
        "[r]1[r][r][r]1",
        "[r]1[r][r][r][r]1",
        "[r]1[r][r][r][r][r]1",
        "[r]1[r][r][r][r][r][r]1",
    ],
    // Basic
    &["[#7;+,\
$([N;H2&+0][$([C,a]);!$([C,a](=O))]),\
$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),\
$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]"],
    // Acidic
    &["[$([C,S](=[O,S,P])-[O;H1,-1])]"],
];

fn ph4_patterns() -> &'static Vec<Vec<ROMol>> {
    static PATTERNS: OnceLock<Vec<Vec<ROMol>>> = OnceLock::new();
    PATTERNS.get_or_init(|| {
        SMARTS_PATTERNS
            .iter()
            .map(|group| {
                group
                    .iter()
                    .map(|s| {
                        mol_from_smarts(s, &SmartsParserParams::default())
                            .ok()
                            .flatten()
                            .expect("the pharmacophore SMARTS are valid")
                    })
                    .collect()
            })
            .collect()
    })
}

/// The pharmacophore colour features of the molecule, appended to
/// `features`. If given and non-empty, only features all of whose atoms are in
/// `atom_subset` are kept.
pub fn find_features(
    mol: &ROMol,
    conf: &Conformer,
    features: &mut Vec<CustomFeature>,
    atom_subset: Option<&[u32]>,
) -> Result<()> {
    let params = SubstructMatchParameters::default();
    for (patt_idx, patts) in ph4_patterns().iter().enumerate() {
        for patt in patts {
            for m in substruct_match(mol, patt, &params)? {
                let ats: Vec<u32> = m.iter().map(|&(_, a)| a as u32).collect();
                if let Some(subset) = atom_subset
                    && !subset.is_empty()
                    && !ats.iter().all(|a| subset.contains(a))
                {
                    continue;
                }
                features.push(CustomFeature {
                    kind: patt_idx as u32 + 1,
                    pos: compute_feature_pos(conf, &ats)?,
                    rad: RADIUS_COLOR,
                    atoms: ats,
                });
            }
        }
    }
    Ok(())
}

/// The centroid of the atoms.
pub fn compute_feature_pos(conf: &Conformer, ats: &[u32]) -> Result<Point3D> {
    let mut pos = Point3D::new(0.0, 0.0, 0.0);
    for &a in ats {
        pos += conf.atom_pos(a)?;
    }
    pos /= ats.len() as f64;
    Ok(pos)
}

/// Applies the transformation to each point of the shape in place.
pub fn apply_transform_to_shape(shape: &mut [f64], xform: &Transform3D) {
    for c in shape.chunks_exact_mut(3) {
        let p = xform.transform_point(&Point3D::new(c[0], c[1], c[2]));
        c.copy_from_slice(&[p.x, p.y, p.z]);
    }
}

/// Applies the transformation to the first `num_points` points of `input`,
/// putting them into `output`.
pub(crate) fn apply_transform_to_shape_into(
    input: &[f64],
    output: &mut [f64],
    num_points: usize,
    xform: &Transform3D,
) {
    for (i, o) in input
        .chunks_exact(3)
        .zip(output.chunks_exact_mut(3))
        .take(num_points)
    {
        let p = xform.transform_point(&Point3D::new(i[0], i[1], i[2]));
        o.copy_from_slice(&[p.x, p.y, p.z]);
    }
}

/// Translates each point of the shape in place.
pub fn translate_shape(shape: &mut [f64], translation: Point3D) {
    for c in shape.chunks_exact_mut(3) {
        c[0] += translation.x;
        c[1] += translation.y;
        c[2] += translation.z;
    }
}

pub(crate) fn translate_shape_into(
    input: &[f64],
    output: &mut [f64],
    num_points: usize,
    translation: Point3D,
) {
    for (i, o) in input
        .chunks_exact(3)
        .zip(output.chunks_exact_mut(3))
        .take(num_points)
    {
        o[0] = i[0] + translation.x;
        o[1] = i[1] + translation.y;
        o[2] = i[2] + translation.z;
    }
}

/// The best Tversky score achievable (with the smaller volume entirely inside
/// the larger one), mixing shape and colour with the options' weight.
pub fn max_score(v1: f64, v2: f64, c1: f64, c2: f64, opts: &ShapeOverlayOptions) -> f64 {
    // Tversky s = O / (A * (R - O) + B * (F - O) + O). There are 2 cases to
    // handle, where p1 < p2 in which case the max overlap is p1, and the
    // opposite.
    let max_part = |p1: f64, p2: f64| -> f64 {
        // If either is 0, the maximum overlap must be 0.
        if p1 == 0.0 || p2 == 0.0 {
            return 0.0;
        }
        if p1 < p2 {
            return p1 / (opts.sim_beta * (p2 - p1) + p1);
        }
        p2 / (opts.sim_alpha * (p1 - p2) + p2)
    };
    let max_st = max_part(v1, v2);
    let max_ct = max_part(c1, c2);
    max_st * (1.0 - opts.opt_param) + max_ct * opts.opt_param
}
