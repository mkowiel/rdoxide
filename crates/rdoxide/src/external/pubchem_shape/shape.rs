//! `External/pubchem_shape/PubChemShape.cpp`: shapes of RDKit molecules (atoms with Bondi radii plus pharmacophore feature
//! points) for the pubchem-align3d optimizer, and the alignment and scoring of molecules and shapes.

use std::collections::BTreeMap;
use std::sync::OnceLock;

use super::functions::*;
use super::neighbor::{ShapeData, neighbor_conformers};
use crate::geometry::{Point3D, Transform3D};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::mol_transforms::{compute_canonical_transform, transform_conformer};
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::dict::RDProps;
use crate::{RDKitError, Result};

const PUBCHEM_FEATURE_NAME: &str = "PUBCHEM_PHARMACOPHORE_FEATURES";
/// the same radius for all feature "atoms"
const RADIUS_COLOR: f64 = 1.08265;

/// Bondi radii (more in Table 12 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3658832/)
fn vdw_radius(z: u32) -> Option<f64> {
    Some(match z {
        1 => 1.10,
        2 => 1.40,
        3 => 1.81,
        4 => 1.53,
        5 => 1.92,
        6 => 1.70,
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

/// Definitions for feature points adapted from Gobbi and Poppinger, Biotech. Bioeng. 61, 47-54 (1998)
const SMARTS_PATTERNS: [&[&str]; 5] = [
    &["[$([N;!H0;v3,v4&+1]),\
$([O,S;H1;+0]),\
n&H1&+0]"], // Donor
    &["[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),\
$([O,S;H0;v2]),\
$([O,S;-]),\
$([N;v3;!$(N-*=[O,N,P,S])]),\
n&H0&+0,\
$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]"], // Acceptor
    &[
        "[r]1[r][r]1",
        "[r]1[r][r][r]1",
        "[r]1[r][r][r][r]1",
        "[r]1[r][r][r][r][r]1",
        "[r]1[r][r][r][r][r][r]1",
    ], // rings
    &["[#7;+,\
$([N;H2&+0][$([C,a]);!$([C,a](=O))]),\
$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),\
$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]"], // Basic
    &["[$([C,S](=[O,S,P])-[O;H1,-1])]"], // Acidic
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
                            .expect("bad smarts")
                    })
                    .collect()
            })
            .collect()
    })
}

/// `ShapeInput`
#[derive(Debug, Clone)]
pub struct ShapeInput {
    pub coord: Vec<f32>,
    pub alpha_vector: Vec<f64>,
    pub atom_type_vector: Vec<u32>,
    pub volume_atom_index_vector: Vec<u32>,
    pub color_atom_type_to_index_vector_map: ColorMap,
    pub shift: Vec<f64>,
    pub inertial_rot: Vec<f64>,
    pub sov: f64,
    pub sof: f64,
}

impl Default for ShapeInput {
    fn default() -> Self {
        ShapeInput {
            coord: Vec::new(),
            alpha_vector: Vec::new(),
            atom_type_vector: Vec::new(),
            volume_atom_index_vector: Vec::new(),
            color_atom_type_to_index_vector_map: ColorMap::new(),
            shift: Vec::new(),
            inertial_rot: vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            sov: 0.0,
            sof: 0.0,
        }
    }
}

/// `ShapeInputOptions`
#[derive(Debug, Clone)]
pub struct ShapeInputOptions {
    /// use the pharmacophore features
    pub use_colors: bool,
    /// hydrogens and dummy atoms are normally skipped
    pub include_dummies: bool,
    /// the radius of dummy atoms (that of Xe)
    pub dummy_radius: f64,
    /// if not empty, use just these atoms
    pub atom_subset: Vec<u32>,
    /// atoms that must not take part in features
    pub not_color_atoms: Vec<u32>,
    /// non-standard radii
    pub atom_radii: Vec<(u32, f64)>,
    /// feature definitions replacing the default ones: (type, position, radius)
    pub custom_features: Vec<(u32, Point3D, f64)>,
    /// put the shape at the origin along its principal axes
    pub normalize: bool,
}

impl Default for ShapeInputOptions {
    fn default() -> Self {
        ShapeInputOptions {
            use_colors: true,
            include_dummies: false,
            dummy_radius: 2.16,
            atom_subset: Vec::new(),
            not_color_atoms: Vec::new(),
            atom_radii: Vec::new(),
            custom_features: Vec::new(),
            normalize: true,
        }
    }
}

type Features = Vec<(Vec<u32>, u32)>;

fn value_error<T>(msg: &str) -> Result<T> {
    Err(RDKitError::Value(msg.to_string()))
}

/// `extractFeatures`: the PubChem feature property of an SDF, else the pharmacophore SMARTS
fn extract_features(mol: &ROMol, opts: &ShapeInputOptions) -> Result<Features> {
    let mut feature_idx_type: Features = Vec::new();
    if !opts.use_colors {
        return Ok(feature_idx_type);
    }
    if let Some(v) = mol.get_prop(PUBCHEM_FEATURE_NAME) {
        let features = v.to_rd_string();
        let atom_type = |t: &str| -> Option<u32> {
            Some(match t {
                "acceptor" => 1,
                "anion" => 2,
                "cation" => 3,
                "donor" => 4,
                "hydrophobe" => 5,
                "rings" => 6,
                _ => return None,
            })
        };
        let mut n = 0usize;
        for line in features.lines() {
            if n == 0 {
                let count: usize = line
                    .trim()
                    .parse()
                    .map_err(|_| RDKitError::Value("bad feature count".into()))?;
                feature_idx_type = vec![(Vec::new(), 0); count];
            } else {
                let f = n - 1;
                if f >= feature_idx_type.len() {
                    return value_error("Too many features");
                }
                let mut t = 0usize;
                for token in line.split(' ') {
                    if t == 0 {
                        let k: usize = token
                            .parse()
                            .map_err(|_| RDKitError::Value("bad feature size".into()))?;
                        feature_idx_type[f].0 = vec![0; k];
                    } else if t <= feature_idx_type[f].0.len() {
                        let idx: u32 = token
                            .parse()
                            .map_err(|_| RDKitError::Value("bad feature atom".into()))?;
                        feature_idx_type[f].0[t - 1] = idx.wrapping_sub(1);
                    } else {
                        match atom_type(token) {
                            Some(ty) => feature_idx_type[f].1 = ty,
                            None => {
                                return value_error(&format!("Invalid feature type: {token}"));
                            }
                        }
                    }
                    t += 1;
                }
                if t != feature_idx_type[f].0.len() + 2 {
                    return value_error("Wrong number of tokens in feature");
                }
            }
            n += 1;
        }
        if n != feature_idx_type.len() + 1 {
            return value_error("Wrong number of features");
        }
    } else {
        let ps = SubstructMatchParameters::default();
        for (patt_idx, patts) in ph4_patterns().iter().enumerate() {
            for patt in patts {
                for m in substruct_match(mol, patt, &ps)? {
                    let ats: Vec<u32> = m.iter().map(|p| p.1 as u32).collect();
                    feature_idx_type.push((ats, patt_idx as u32 + 1));
                }
            }
        }
    }
    Ok(feature_idx_type)
}

fn atom_in_subset(idx: u32, opts: &ShapeInputOptions) -> bool {
    opts.atom_subset.is_empty() || opts.atom_subset.contains(&idx)
}

fn atom_allowed_in_color(idx: u32, opts: &ShapeInputOptions) -> bool {
    !opts.not_color_atoms.contains(&idx)
}

fn atom_radius(idx: u32, opts: &ShapeInputOptions) -> f64 {
    opts.atom_radii
        .iter()
        .find(|p| p.0 == idx)
        .map_or(-1.0, |p| p.1)
}

/// `PrepareConformer`: the shape of a conformer, translated to the origin and aligned along its principal axes
pub fn prepare_conformer(
    mol: &ROMol,
    conf_id: i32,
    opts: &ShapeInputOptions,
) -> Result<ShapeInput> {
    let mut res = ShapeInput::default();
    let feature_idx_type = extract_features(mol, opts)?;
    let conformer = mol.conformer(conf_id)?;
    if !conformer.is_3d() {
        return value_error("Conformer must be 3D");
    }
    let n_atoms = mol.num_atoms();
    let mut n_alignment_atoms = n_atoms as usize;
    n_alignment_atoms += if opts.custom_features.is_empty() {
        feature_idx_type.len()
    } else {
        opts.custom_features.len()
    };
    let mut rad_vector = vec![0.0f64; n_alignment_atoms];
    res.atom_type_vector = vec![0; n_alignment_atoms];

    let mut conf_cp = conformer.clone();
    let mut trans: Option<Transform3D> = None;
    let mut initial_centroid = Point3D::default();
    if opts.normalize {
        // this ignores hydrogens by default, which is what is needed here
        let (t, _) = compute_canonical_transform(mol, &conf_cp, None, false, true)?;
        transform_conformer(&mut conf_cp, &t);
        let mut k = 0;
        for i in 0..3 {
            for j in 0..3 {
                res.inertial_rot[k] = t.get_val(i, j);
                k += 1;
            }
        }
        // the initial centroid is the translation of the canonical transform rotated back by its inverse rotation
        initial_centroid = Point3D::new(t.get_val(0, 3), t.get_val(1, 3), t.get_val(2, 3));
        let mut inv = Transform3D::identity();
        for i in 0..3 {
            for j in 0..3 {
                inv.set_val(i, j, t.get_val(j, i));
            }
        }
        initial_centroid = inv.transform_point(&initial_centroid);
        trans = Some(t);
    }

    // radii, and the average coordinates of the selected atoms
    let mut ave = Point3D::default();
    let mut n_selected = 0usize;
    for i in 0..n_atoms {
        if !atom_in_subset(i, opts) {
            continue;
        }
        let rad = atom_radius(i, opts);
        if rad > 0.0 {
            rad_vector[n_selected] = rad;
            n_selected += 1;
            ave += conf_cp.atom_pos(i)?;
        } else {
            let z = mol.atom(i).atomic_num();
            if z > 1 {
                ave += conf_cp.atom_pos(i)?;
                match vdw_radius(z) {
                    Some(r) => {
                        rad_vector[n_selected] = r;
                        n_selected += 1;
                    }
                    None => return value_error(&format!("No VdW radius for atom with Z={z}")),
                }
            } else if opts.include_dummies && z == 0 {
                ave += conf_cp.atom_pos(i)?;
                rad_vector[n_selected] = opts.dummy_radius;
                n_selected += 1;
            }
        }
    }
    ave /= n_selected as f64;
    res.shift = vec![-ave.x, -ave.y, -ave.z];
    res.coord = vec![0.0; n_alignment_atoms * 3];

    // coordinates of the atoms
    let mut j = 0;
    for i in 0..n_atoms {
        if !atom_in_subset(i, opts) {
            continue;
        }
        let z = mol.atom(i).atomic_num();
        if z > 1 || (opts.include_dummies && z == 0) {
            let mut pos = conf_cp.atom_pos(i)?;
            if opts.normalize {
                pos -= ave;
            }
            res.coord[j * 3] = pos.x as f32;
            res.coord[j * 3 + 1] = pos.y as f32;
            res.coord[j * 3 + 2] = pos.z as f32;
            j += 1;
        }
    }

    let mut num_features = 0usize;
    if opts.custom_features.is_empty() {
        // the feature position is the average of the coordinates of its atoms
        for (atoms, ftype) in &feature_idx_type {
            let mut floc = Point3D::default();
            let mut n_sel = 0usize;
            for &idx in atoms {
                if !atom_in_subset(idx, opts) || !atom_allowed_in_color(idx, opts) {
                    continue;
                }
                if idx >= n_atoms || mol.atom(idx).atomic_num() <= 1 {
                    return value_error("Invalid feature atom index");
                }
                floc += conf_cp.atom_pos(idx)?;
                n_sel += 1;
            }
            if n_sel == atoms.len() {
                floc /= n_sel as f64;
                if opts.normalize {
                    floc -= ave;
                }
                let k = n_selected + num_features;
                res.coord[k * 3] = floc.x as f32;
                res.coord[k * 3 + 1] = floc.y as f32;
                res.coord[k * 3 + 2] = floc.z as f32;
                rad_vector[k] = RADIUS_COLOR;
                res.atom_type_vector[k] = *ftype;
                num_features += 1;
            }
        }
    } else if opts.use_colors {
        for (ftype, pos, radius) in &opts.custom_features {
            let mut floc = *pos;
            if let Some(t) = &trans {
                floc = t.transform_point(&floc);
            }
            floc -= ave;
            let k = n_selected + num_features;
            res.coord[k * 3] = floc.x as f32;
            res.coord[k * 3 + 1] = floc.y as f32;
            res.coord[k * 3 + 2] = floc.z as f32;
            rad_vector[k] = *radius;
            res.atom_type_vector[k] = *ftype;
            num_features += 1;
        }
    }

    // cut the vectors down to the atoms and features used
    n_alignment_atoms = n_selected + num_features;
    res.coord.truncate(n_alignment_atoms * 3);
    rad_vector.truncate(n_alignment_atoms);
    res.atom_type_vector.truncate(n_alignment_atoms);
    res.alpha_vector = set_alpha(&rad_vector);

    res.volume_atom_index_vector = get_volume_atom_index_vector(&res.atom_type_vector);
    res.sov = shape_overlap(
        &res.coord,
        &res.alpha_vector,
        &res.volume_atom_index_vector,
        &res.coord,
        &res.alpha_vector,
        &res.volume_atom_index_vector,
    );
    if !feature_idx_type.is_empty() {
        res.color_atom_type_to_index_vector_map =
            get_color_atom_type_to_index_vector_map(&res.atom_type_vector);
        res.sof = feature_overlap(
            &res.coord,
            &res.alpha_vector,
            &res.color_atom_type_to_index_vector_map,
            &res.coord,
            &res.alpha_vector,
            &res.color_atom_type_to_index_vector_map,
            None,
        );
    }
    res.shift = vec![initial_centroid.x, initial_centroid.y, initial_centroid.z];
    Ok(res)
}

/// `AlignShape`: (shape Tanimoto, colour Tanimoto); `fit_shape.coord` is transformed and `matrix` receives the overlay
pub fn align_shape(
    ref_shape: &ShapeInput,
    fit_shape: &mut ShapeInput,
    matrix: &mut [f32; 12],
    opt_param: f64,
    max_preiters: u32,
    max_postiters: u32,
) -> (f64, f64) {
    let joint = get_joint_color_type_set(&ref_shape.atom_type_vector, &fit_shape.atom_type_vector);
    let mut ref_map = ref_shape.color_atom_type_to_index_vector_map.clone();
    restrict_color_map(&mut ref_map, &joint);
    let mut fit_map = fit_shape.color_atom_type_to_index_vector_map.clone();
    restrict_color_map(&mut fit_map, &joint);
    let r = ShapeData {
        coord: &ref_shape.coord,
        alpha: &ref_shape.alpha_vector,
        volume_atoms: &ref_shape.volume_atom_index_vector,
        color_map: &ref_map,
        sov: ref_shape.sov,
        sof: ref_shape.sof,
    };
    let f = ShapeData {
        coord: &fit_shape.coord,
        alpha: &fit_shape.alpha_vector,
        volume_atoms: &fit_shape.volume_atom_index_vector,
        color_map: &fit_map,
        sov: fit_shape.sov,
        sof: fit_shape.sof,
    };
    let res = neighbor_conformers(
        &r,
        &f,
        !joint.is_empty(),
        max_preiters,
        max_postiters,
        opt_param,
    );
    *matrix = res.matrix;
    fit_shape.coord = apply_rot_trans_matrix(&fit_shape.coord, &res.matrix);
    (res.st, res.ct)
}

/// `TransformConformer`
pub fn transform_conformer_by_overlay(
    final_trans: &[f64],
    final_rot: &[f64],
    matrix: &[f32],
    fit_shape: &ShapeInput,
    fit_conf: &mut Conformer,
) {
    let mut transform0 = Transform3D::identity();
    transform0.set_translation(&Point3D::new(
        fit_shape.shift[0],
        fit_shape.shift[1],
        fit_shape.shift[2],
    ));
    let mut transform1 = Transform3D::identity();
    for i in 0..3 {
        for j in 0..3 {
            transform1.set_val(i, j, fit_shape.inertial_rot[i * 3 + j]);
        }
    }
    let m = |k: usize| f64::from(matrix[k]);
    let mut transform2 = Transform3D::identity();
    for (i, row) in [[0, 1, 2, 9], [3, 4, 5, 10], [6, 7, 8, 11]]
        .iter()
        .enumerate()
    {
        for (j, &k) in row.iter().enumerate() {
            transform2.set_val(i, j, m(k));
        }
    }
    let mut transform3 = Transform3D::identity();
    for i in 0..3 {
        for j in 0..3 {
            transform3.set_val(i, j, final_rot[i * 3 + j]);
        }
        transform3.set_val(i, 3, final_trans[i]);
    }
    let final_transform = transform3
        .compose(&transform2)
        .compose(&transform1)
        .compose(&transform0);
    transform_conformer(fit_conf, &final_transform);
}

/// `AlignMolecule(refShape, fit, ...)`: aligns the conformer of `fit`; returns the Tanimotos and sets the properties
/// `shape_align_shape_tanimoto` and `shape_align_color_tanimoto` on `fit`
#[allow(clippy::too_many_arguments)]
pub fn align_molecule_to_shape(
    ref_shape: &ShapeInput,
    fit: &mut ROMol,
    matrix: &mut [f32; 12],
    shape_opts: &ShapeInputOptions,
    fit_conf_id: i32,
    opt_param: f64,
    max_preiters: u32,
    max_postiters: u32,
    apply_ref_shift: bool,
) -> Result<(f64, f64)> {
    let mut fit_shape = prepare_conformer(fit, fit_conf_id, shape_opts)?;
    let tanis = align_shape(
        ref_shape,
        &mut fit_shape,
        matrix,
        opt_param,
        max_preiters,
        max_postiters,
    );
    let mut final_trans = vec![0.0, 0.0, 0.0];
    let mut final_rot = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    if apply_ref_shift {
        final_trans = vec![
            -ref_shape.shift[0],
            -ref_shape.shift[1],
            -ref_shape.shift[2],
        ];
        let r = &ref_shape.inertial_rot;
        final_rot = vec![r[0], r[3], r[6], r[1], r[4], r[7], r[2], r[5], r[8]];
    }
    let conf = fit.conformer_mut(fit_conf_id)?;
    transform_conformer_by_overlay(&final_trans, &final_rot, matrix, &fit_shape, conf);
    fit.set_prop("shape_align_shape_tanimoto", tanis.0, false)?;
    fit.set_prop("shape_align_color_tanimoto", tanis.1, false)?;
    Ok(tanis)
}

/// `AlignMolecule(ref, fit, ...)`
#[allow(clippy::too_many_arguments)]
pub fn align_molecule(
    reference: &ROMol,
    fit: &mut ROMol,
    matrix: &mut [f32; 12],
    ref_opts: &ShapeInputOptions,
    fit_opts: &ShapeInputOptions,
    ref_conf_id: i32,
    fit_conf_id: i32,
    opt_param: f64,
    max_preiters: u32,
    max_postiters: u32,
) -> Result<(f64, f64)> {
    if ref_opts.use_colors != fit_opts.use_colors {
        crate::rd_log_warning!(
            "useColor values inconsistent between the reference and fit molecules. Colors will not be used in the alignment."
        );
    }
    let ref_shape = prepare_conformer(reference, ref_conf_id, ref_opts)?;
    align_molecule_to_shape(
        &ref_shape,
        fit,
        matrix,
        fit_opts,
        fit_conf_id,
        opt_param,
        max_preiters,
        max_postiters,
        true,
    )
}

/// `ScoreShape`: (shape Tanimoto, colour Tanimoto) of two shapes as they are placed
pub fn score_shape(shape1: &ShapeInput, shape2: &ShapeInput, use_colors: bool) -> (f64, f64) {
    let shape_ov = shape_overlap(
        &shape1.coord,
        &shape1.alpha_vector,
        &shape1.volume_atom_index_vector,
        &shape2.coord,
        &shape2.alpha_vector,
        &shape2.volume_atom_index_vector,
    );
    let shape = shape_ov / (shape1.sov + shape2.sov - shape_ov);
    let joint = get_joint_color_type_set(&shape1.atom_type_vector, &shape2.atom_type_vector);
    let mut color = 0.0;
    if use_colors {
        let mut m1 = shape1.color_atom_type_to_index_vector_map.clone();
        restrict_color_map(&mut m1, &joint);
        let mut m2 = shape2.color_atom_type_to_index_vector_map.clone();
        restrict_color_map(&mut m2, &joint);
        let feature = feature_overlap(
            &shape1.coord,
            &shape1.alpha_vector,
            &m1,
            &shape2.coord,
            &shape2.alpha_vector,
            &m2,
            None,
        );
        color = feature / (shape1.sof + shape2.sof - feature);
    }
    (shape, color)
}

/// `ScoreMolecule(shape, mol, ...)`
pub fn score_molecule_to_shape(
    shape: &ShapeInput,
    mol: &ROMol,
    mol_opts: &ShapeInputOptions,
    mol_conf_id: i32,
) -> Result<(f64, f64)> {
    let mut opts = mol_opts.clone();
    opts.normalize = false;
    let shape2 = prepare_conformer(mol, mol_conf_id, &opts)?;
    Ok(score_shape(shape, &shape2, opts.use_colors))
}

/// `ScoreMolecule(mol1, mol2, ...)`
pub fn score_molecule(
    mol1: &ROMol,
    mol2: &ROMol,
    opts1: &ShapeInputOptions,
    opts2: &ShapeInputOptions,
    conf1: i32,
    conf2: i32,
) -> Result<(f64, f64)> {
    let mut o1 = opts1.clone();
    o1.normalize = false;
    let shape1 = prepare_conformer(mol1, conf1, &o1)?;
    let mut o2 = opts2.clone();
    o2.normalize = false;
    let shape2 = prepare_conformer(mol2, conf2, &o2)?;
    let use_colors = opts1.use_colors && opts2.use_colors;
    Ok(score_shape(&shape1, &shape2, use_colors))
}

/// a map of the colour feature types to their atoms, for tests
pub fn color_types(shape: &ShapeInput) -> BTreeMap<u32, usize> {
    shape
        .color_atom_type_to_index_vector_map
        .iter()
        .map(|(k, v)| (*k, v.len()))
        .collect()
}
