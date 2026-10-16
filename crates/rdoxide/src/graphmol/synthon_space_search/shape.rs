//! `SynthonSpaceShapeSearcher.cpp`: searching a SynthonSpace by Gaussian
//! shape similarity.

use std::collections::{BTreeMap, HashMap};
use std::f64::consts::PI;
use std::sync::Arc;
use std::time::Instant;

use super::details;
use super::hit_set::SynthonSpaceHitSet;
use super::params::{MAX_CONNECTOR_NUM, SynthonSpaceSearchParams, ThreadMode};
use super::search_results::SearchResults;
use super::searcher::{
    SearchResultCallback, SearcherBase, SynthonSearcher, ToTry, check_possible_hits,
    make_hits_from_to_try, search, search_with_callback,
};
use super::shape_build::{generate_isomer_conformers, split_dummy_dummy_bonds};
use super::shape_input::SynthonShapeInput;
use super::space::SynthonSpace;
use super::synthon::Synthon;
use super::synthon_set::SynthonSet;
use crate::geometry::{Point3D, Transform3D};
use crate::graphmol::bond::BondType;
use crate::graphmol::chirality::from_coords::assign_stereochemistry_from_3d;
use crate::graphmol::dist_geom_helpers::EmbedParameters;
use crate::graphmol::find_rings::get_mol_frags;
use crate::graphmol::gaussian_shape::{
    CARBON_RAD, CustomFeature, ShapeInput, ShapeInputOptions, ShapeOverlayOptions, StartMode,
    align_molecule_to_shape, align_shape, find_features, max_score, score_shape,
};
use crate::graphmol::mol_transforms::{get_bond_length, transform_conformer};
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::rdgeneral::threads::par_map;
use crate::{RDKitError, Result};

use super::fingerprint::{check_lines, check_query};

/// The best overlay of a synthon onto a fragment: the combination score, the
/// number of the synthon shape and the transformation to apply to the shape
/// to get the overlay that gave the score.
#[derive(Clone, Debug)]
pub struct SynthonOverlay {
    pub score: f64,
    pub shape_num: usize,
    pub xform: Option<Arc<Transform3D>>,
}

impl SynthonOverlay {
    fn none() -> Self {
        SynthonOverlay {
            score: -1.0,
            shape_num: 0,
            xform: None,
        }
    }
}

/// (the number of the fragment shape, the address of the synthon)
type SimKey = (usize, usize);

fn synthon_key(synthon: &Synthon) -> usize {
    synthon as *const Synthon as usize
}

fn frag_key(frag: &Arc<ROMol>) -> usize {
    Arc::as_ptr(frag) as usize
}

pub struct ShapeSearcher<'a> {
    base: SearcherBase<'a>,
    /// Shape for the query - only ever 1 conformer.
    query_shape: Option<SynthonShapeInput>,
    /// The fragment shapes for this search, derived from the query.
    frag_shapes_pool: Vec<Option<SynthonShapeInput>>,
    /// Number of the fragment's shape in `frag_shapes_pool`, keyed by the
    /// address of the fragment.
    frag_shapes: HashMap<usize, usize>,
    /// Precomputed similarities between fragments and synthons. This speeds
    /// things up because synthons are re-used and for the shape search we
    /// have to calculate the similarity between every query fragment and all
    /// the synthons in each SynthonSet. There's an entry only if the
    /// fragment->synthon similarity exceeded the threshold.
    frag_synthon_sims: HashMap<SimKey, SynthonOverlay>,
}

impl<'a> ShapeSearcher<'a> {
    pub fn new(
        query: &ROMol,
        params: &SynthonSpaceSearchParams,
        space: &'a SynthonSpace,
        has_space: bool,
    ) -> Result<Self> {
        if has_space && space.num_conformers() == 0 {
            return Err(RDKitError::Runtime(
                "No conformers found in SynthonSpaceSearch".into(),
            ));
        }
        let mut base = SearcherBase::new(query, params, space)?;
        base.frag_set_uniquify_mode = super::params::FragSetUniquifyMode::ByAtoms;
        let mut s = ShapeSearcher {
            base,
            query_shape: None,
            frag_shapes_pool: Vec::new(),
            frag_shapes: HashMap::new(),
            frag_synthon_sims: HashMap::new(),
        };
        if !has_space {
            let q = s.base.query.clone();
            s.build_query_shape(&q)?;
        }
        Ok(s)
    }

    fn query_shape(&self) -> &SynthonShapeInput {
        self.query_shape.as_ref().expect("query shape built")
    }

    pub fn has_precomputed_sims(&self) -> bool {
        !self.frag_synthon_sims.is_empty()
    }

    /// Use the precomputed sims to decide if the fragment matched the
    /// synthon.
    pub fn frag_matched_synthon(
        &self,
        frag_shape: usize,
        synthon: &Synthon,
    ) -> Option<&SynthonOverlay> {
        self.frag_synthon_sims
            .get(&(frag_shape, synthon_key(synthon)))
    }

    fn build_query_shape(&mut self, mol: &ROMol) -> Result<Vec<CustomFeature>> {
        crate::rd_log_info!("Generating query shapes for {}", details::smiles_of(mol));
        let mut all_features = Vec::new();
        find_features(mol, mol.conformer(-1)?, &mut all_features, None)?;
        let opts = ShapeInputOptions {
            custom_features: vec![all_features.clone()],
            ..Default::default()
        };
        self.query_shape = Some(SynthonShapeInput::new(
            mol,
            -1,
            &opts,
            &self.base.params.shape_overlay_options,
        )?);
        Ok(all_features)
    }

    fn frag_shape_of(&self, frag: &Arc<ROMol>) -> Option<usize> {
        let idx = *self.frag_shapes.get(&frag_key(frag))?;
        self.frag_shapes_pool[idx].as_ref().map(|_| idx)
    }

    /// Fills in the fragment->synthon similarities.
    fn compute_frag_synthon_sims(
        &mut self,
        end_time: Option<Instant>,
        min_frag_set_size: &[usize],
    ) -> Result<()> {
        let threshold =
            self.base.params.similarity_cutoff - self.base.params.frag_similarity_adjuster;
        let space = self.base.space;
        let overlay_opts = self.base.params.shape_overlay_options.clone();
        let num_threads = self.base.params.num_threads;
        let mut to_do: Vec<(usize, usize)> = Vec::new();
        // The pairs are done in chunks, to limit the memory used.
        const CHUNK: usize = 2_500_000;
        let pool = space.pool();
        let mut sims = HashMap::new();
        let mut process = |to_do: &mut Vec<(usize, usize)>| -> Result<()> {
            let frag_shapes = &self.frag_shapes_pool;
            let counter = std::sync::atomic::AtomicUsize::new(0);
            let timed_out = std::sync::atomic::AtomicBool::new(false);
            let res = par_map(
                to_do,
                num_threads,
                |&(f, s)| -> Result<Option<(SimKey, SynthonOverlay)>> {
                    use std::sync::atomic::Ordering::Relaxed;
                    if timed_out.load(Relaxed) {
                        return Ok(None);
                    }
                    if counter.fetch_add(1, Relaxed) % 100 == 0 && details::check_time_out(end_time)
                    {
                        timed_out.store(true, Relaxed);
                        return Ok(None);
                    }
                    let Some(frag_shape) = &frag_shapes[f] else {
                        return Ok(None);
                    };
                    let synthon = &pool[s];
                    if !synthon.has_shapes() {
                        return Ok(None);
                    }
                    let sim = best_sim_synthon_onto_fragment(
                        frag_shape,
                        synthon,
                        threshold,
                        &overlay_opts,
                    )?;
                    Ok((sim.score >= threshold).then(|| ((f, synthon_key(synthon)), sim)))
                },
            );
            for r in res {
                if let Some((k, v)) = r? {
                    sims.insert(k, v);
                }
            }
            to_do.clear();
            Ok(())
        };
        for f in 0..self.frag_shapes_pool.len() {
            for (s, synthon) in pool.iter().enumerate() {
                // If the smallest fragment set that this fragment is in is
                // bigger than the largest SynthonSet the synthon is in, we
                // won't ever need to know the similarity between them, so
                // skip.
                if min_frag_set_size[f] <= synthon.max_synthon_set_size() as usize {
                    to_do.push((f, s));
                }
                if to_do.len() == CHUNK {
                    process(&mut to_do)?;
                }
            }
        }
        process(&mut to_do)?;
        self.frag_synthon_sims = sims;
        Ok(())
    }
}

/// The atoms in `frag` that are in the query have the `ORIG_IDX` property.
fn orig_idx(props: &impl RDProps) -> Option<u32> {
    props.get_prop("ORIG_IDX").and_then(|v| v.as_uint())
}

fn select_frag_features(
    all_features: &[CustomFeature],
    query: &ROMol,
    frag_atoms: &[u32],
) -> Vec<CustomFeature> {
    let mut in_frag = vec![false; query.num_atoms() as usize];
    for &f in frag_atoms {
        in_frag[f as usize] = true;
    }
    all_features
        .iter()
        .filter(|feat| {
            // Single atom features that are dummy atoms in the query
            // shouldn't be used, but larger features, which probably just
            // means rings, should be kept.
            if feat.atoms.len() == 1
                && in_frag[feat.atoms[0] as usize]
                && query.atom(feat.atoms[0]).atomic_num() == 0
            {
                return false;
            }
            feat.atoms.iter().all(|&a| in_frag[a as usize])
        })
        .cloned()
        .collect()
}

fn generate_shapes(
    query: &ROMol,
    frag: &ROMol,
    all_features: &[CustomFeature],
) -> Result<SynthonShapeInput> {
    // The fragSets molecules will have their atoms labelled with ORIG_IDX
    // apart from the dummy atoms, but we need coords for them, too. They are
    // normally copied from the atom at the other end of the broken bond so
    // find that atom too. Work on a copy the query.
    let mut query_cp = query.clone();
    let mut frag_atoms: Vec<u32> = Vec::with_capacity(frag.num_atoms() as usize);
    let mut in_frag = vec![false; query_cp.num_atoms() as usize];
    for atom in frag.atoms() {
        if let Some(orig) = orig_idx(atom) {
            frag_atoms.push(orig);
            in_frag[orig as usize] = true;
        }
    }
    frag_atoms.sort_unstable();
    frag_atoms.dedup();
    let mut dummy_radii: Vec<(u32, f64)> = Vec::new();
    for atom in frag.atoms() {
        if atom.atomic_num() == 0
            && atom.isotope() >= 1
            && atom.isotope() as usize <= MAX_CONNECTOR_NUM
        {
            let Some(nbr) = frag.atom_neighbors(atom.idx()).next() else {
                continue;
            };
            let orig_nbr =
                orig_idx(frag.atom(nbr)).ok_or_else(|| RDKitError::Key("ORIG_IDX".into()))?;
            let nbr_nbrs: Vec<u32> = query_cp.atom_neighbors(orig_nbr).collect();
            for nbr_nbr in nbr_nbrs {
                if !in_frag[nbr_nbr as usize] {
                    dummy_radii.push((nbr_nbr, 2.16));
                    // Make it a dummy atom in the copy so it's correct in the
                    // Shape.
                    query_cp.atom_mut(nbr_nbr).set_atomic_num(0);
                }
            }
        }
    }
    // Put the dummy atoms into the fragAtoms as we want them in the shape,
    // too. If the fragment was formed by taking a single atom out of a ring,
    // the dummy atom will appear twice in the list, so fix that. It will mean
    // that in the shape, there will be a dummy in a ring. A synthon would
    // have 2, but in the same place.
    dummy_radii.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    dummy_radii.dedup();
    frag_atoms.extend(dummy_radii.iter().map(|p| p.0));
    // Break any dummy-dummy bonds
    split_dummy_dummy_bonds(&mut query_cp)?;
    let frag_features = select_frag_features(all_features, &query_cp, &frag_atoms);
    let opts = ShapeInputOptions {
        atom_subset: frag_atoms,
        atom_radii: dummy_radii,
        custom_features: vec![frag_features],
        ..Default::default()
    };
    SynthonShapeInput::new(&query_cp, -1, &opts, &ShapeOverlayOptions::default())
}

/// This is very similar to `details::map_frags_by_smiles`: the key is an
/// underscore-separated list of the indices of the query atoms in the
/// fragment, and the value the (fragment set, fragment) positions.
fn map_frags_by_atoms(frag_sets: &[Vec<Arc<ROMol>>]) -> BTreeMap<String, Vec<(usize, usize)>> {
    let mut atoms_to_frags: BTreeMap<String, Vec<(usize, usize)>> = BTreeMap::new();
    const DUMMY_IDX: u32 = 10000;
    for (i, frag_set) in frag_sets.iter().enumerate() {
        for (j, frag) in frag_set.iter().enumerate() {
            let mut at_idxs: Vec<u32> = Vec::new();
            for a in frag.atoms() {
                if let Some(orig) = orig_idx(a) {
                    at_idxs.push(orig);
                } else if a.atomic_num() == 0 && a.isotope() > 0 {
                    at_idxs.push(DUMMY_IDX + a.isotope());
                }
            }
            at_idxs.sort_unstable();
            let key = at_idxs
                .iter()
                .map(u32::to_string)
                .collect::<Vec<_>>()
                .join("_");
            atoms_to_frags.entry(key).or_default().push((i, j));
        }
    }
    atoms_to_frags
}

/// Aligns the synth dummy and its neighbour with the frag equivalents and
/// returns the transformation that did it.
fn align_dummies(
    frag_dummy: Point3D,
    frag_dummy_nbr: Point3D,
    synth_dummy: Point3D,
    synth_dummy_nbr: Point3D,
) -> Result<Transform3D> {
    let frag_vec = frag_dummy.direction_vector(&frag_dummy_nbr)?;
    let synth_vec = synth_dummy.direction_vector(&synth_dummy_nbr)?;
    let mut axis = synth_vec.cross_product(&frag_vec);
    axis.normalize()?;
    let cos_t = frag_vec.dot_product(&synth_vec);
    let sin_t = (1.0 - cos_t * cos_t).sqrt();
    let mut rot = Transform3D::identity();
    rot.set_rotation_cs(cos_t, sin_t, &axis);
    let mut to_origin = Transform3D::identity();
    to_origin.set_translation(&Point3D::new(
        -synth_dummy.x,
        -synth_dummy.y,
        -synth_dummy.z,
    ));
    let mut from_origin = Transform3D::identity();
    from_origin.set_translation(&frag_dummy);
    Ok(from_origin.compose(&rot).compose(&to_origin))
}

fn shape_point(shape: &ShapeInput, idx: u32) -> Point3D {
    let c = &shape.coords()[3 * idx as usize..3 * idx as usize + 3];
    Point3D::new(c[0], c[1], c[2])
}

/// Makes sure the synthon dummy atom hasn't moved too far from the fragment
/// one or the directions haven't diverged too much.
fn check_dummies(
    frag_dummy: Point3D,
    frag_dummy_nbr: Point3D,
    synth_shape: &ShapeInput,
    synth_dummy_idx: u32,
    synth_dummy_nbr_idx: u32,
    dist_threshold_sq: f64,
    angle_threshold: f64,
) -> Result<bool> {
    let sd = shape_point(synth_shape, synth_dummy_idx);
    if (frag_dummy - sd).length_sq() > dist_threshold_sq {
        return Ok(false);
    }
    let sdn = shape_point(synth_shape, synth_dummy_nbr_idx);
    let v1 = frag_dummy.direction_vector(&frag_dummy_nbr)?;
    let v2 = sd.direction_vector(&sdn)?;
    Ok(v1.angle_to(&v2) <= angle_threshold)
}

/// Finds the best similarity of a fragment shape with a synthon shape, within
/// the threshold, making sure that they start with dummy atoms aligned and
/// the synthon's dummy atom doesn't move too far from the fragment's in the
/// final overlay. Assumes the fragment shape hasn't been normalised, so is in
/// its original position as it was in the input molecule conformation. Keeps
/// track of the transformation of the best synthon shape onto the fragment
/// for later use.
fn best_sim_synthon_onto_fragment(
    frag_shape: &SynthonShapeInput,
    synthon: &Synthon,
    threshold: f64,
    shape_overlay_opts: &ShapeOverlayOptions,
) -> Result<SynthonOverlay> {
    // This should only be called if the synthon has shapes.
    let synth_shapes = synthon.shapes().expect("synthon has shapes");
    if frag_shape
        .shapes()
        .max_possible_similarity(synth_shapes.shapes(), &ShapeOverlayOptions::default())?
        < threshold
    {
        return Ok(SynthonOverlay::none());
    }
    // Work on a copy so that alpha values can be changed if need be.
    let mut frag_shape_cp = frag_shape.clone();
    let mut opts = shape_overlay_opts.clone();
    opts.start_mode = StartMode::Rotate0;
    opts.normalize = false;
    // Rotations of 0, 90, 180, 270 degrees.
    let sin_t = [0.0, (PI / 2.0).sin(), PI.sin(), (PI * 1.5).sin()];
    let cos_t = [1.0, (PI / 2.0).cos(), PI.cos(), (PI * 1.5).cos()];
    // Tolerances for how far the synthon shape moves when overlaid onto the
    // fragment.
    const DIST_TOL_SQ: f64 = 4.0; // Slightly over a bond's length
    // There's a big lever, potentially, so make the tolerance quite shallow.
    const ANGLE_TOL: f64 = PI / 4.0;
    let mut best_synth_shape = 0;
    let mut best_xform = Transform3D::identity();

    // The frag_shape will only have one shape, since the query is only
    // allowed to have 1 conformation. If it doesn't have any dummies, just
    // do a normal alignment.
    if frag_shape_cp.num_dummy_atoms() == 0 {
        let synth_shape_cp = synth_shapes.shapes().clone();
        let (sim, pair) = frag_shape_cp.shapes_mut().best_similarity(
            &synth_shape_cp,
            threshold,
            &ShapeOverlayOptions::default(),
        )?;
        if let Some(pair) = pair {
            best_synth_shape = pair.fit_shape;
            best_xform = pair.xform;
        }
        return Ok(SynthonOverlay {
            score: sim[0],
            shape_num: best_synth_shape,
            xform: Some(Arc::new(best_xform)),
        });
    }
    // For every dummy atom in the frag_shape
    //    For every shape in the synthon:
    //       For every dummy atom in the synthon shape
    //          Align the dummy atoms and their neighbours
    //          Rotate the synthon shape by 0, 90, 180 and 270 degrees
    //          Do the alignment from that start point
    //          Keep the best score and the transformation that caused it.
    let mut best_score = 0.0;
    let frag_dummies: Vec<(u32, u32)> = frag_shape_cp.dummy_atoms_and_nbrs().to_vec();
    for &(frag_dummy_idx, frag_dummy_nbr_idx) in &frag_dummies {
        let frag_dummy = shape_point(frag_shape_cp.shapes(), frag_dummy_idx);
        let frag_dummy_nbr = shape_point(frag_shape_cp.shapes(), frag_dummy_nbr_idx);
        // Set the alpha values of the other dummies -ve so they are ignored.
        for &(other, _) in &frag_dummies {
            if other != frag_dummy_idx && frag_shape_cp.shapes().alphas()[other as usize] > 0.0 {
                frag_shape_cp.shapes_mut().negate_alpha(other as usize)?;
            }
        }
        for ssn in 0..synth_shapes.shapes().num_shapes() {
            for &(synth_dummy_idx, synth_dummy_nbr_idx) in synth_shapes.dummy_atoms_and_nbrs() {
                // Fresh working copy every time.
                let mut single_shape = ShapeInput::from_shape(synth_shapes.shapes(), ssn)?;
                let xform = align_dummies(
                    frag_dummy,
                    frag_dummy_nbr,
                    shape_point(&single_shape, synth_dummy_idx),
                    shape_point(&single_shape, synth_dummy_nbr_idx),
                )?;
                single_shape.transform_coords(&xform);
                let single_dummy = shape_point(&single_shape, synth_dummy_idx);
                let rot_axis = single_dummy
                    .direction_vector(&shape_point(&single_shape, synth_dummy_nbr_idx))?;
                // Get the new dummy->neighbour axis to rotate about.
                let mut to_origin = Transform3D::identity();
                to_origin.set_translation(&Point3D::new(
                    -single_dummy.x,
                    -single_dummy.y,
                    -single_dummy.z,
                ));
                let mut from_origin = Transform3D::identity();
                from_origin.set_translation(&single_dummy);
                // Rotate 4 times around the axis
                for k in 0..4 {
                    // Work on a copy so the transformation will be relevant at
                    // the end.
                    let mut single_shape_cp = single_shape.clone();
                    let mut rot = Transform3D::identity();
                    // upstream passes these as (sin, cos) to a function that
                    // takes (cos, sin)
                    rot.set_rotation_cs(sin_t[k], cos_t[k], &rot_axis);
                    let full_xform = from_origin.compose(&rot).compose(&to_origin);
                    single_shape_cp.transform_coords(&full_xform);
                    let (scores, ovly_xform) =
                        align_shape(frag_shape_cp.shapes(), &mut single_shape_cp, &opts)?;
                    if scores[0] > best_score
                        // Reject anything where the synthon dummies have gone
                        // too far.
                        && check_dummies(
                            frag_dummy,
                            frag_dummy_nbr,
                            &single_shape_cp,
                            synth_dummy_idx,
                            synth_dummy_nbr_idx,
                            DIST_TOL_SQ,
                            ANGLE_TOL,
                        )?
                    {
                        best_score = scores[0];
                        // The total transform that the synthon undergoes to
                        // get this score.
                        best_xform = ovly_xform.compose(&full_xform).compose(&xform);
                        best_synth_shape = ssn;
                    }
                }
            }
        }
        // Set the alpha values of the other dummies back again.
        for &(other, _) in &frag_dummies {
            if other != frag_dummy_idx && frag_shape_cp.shapes().alphas()[other as usize] < 0.0 {
                frag_shape_cp.shapes_mut().negate_alpha(other as usize)?;
            }
        }
    }
    Ok(SynthonOverlay {
        score: best_score,
        shape_num: best_synth_shape,
        xform: Some(Arc::new(best_xform)),
    })
}

/// Takes the fragged mol ShapeSets and flags all those synthons that have a
/// fragment as a similarity match.
fn get_hit_synthons(
    searcher: &ShapeSearcher<'_>,
    frag_shapes: &[usize],
    similarity_cutoff: f64,
    reaction: &SynthonSet,
    synthon_set_order: &[usize],
) -> Vec<Vec<usize>> {
    let pool = searcher.base.space.pool();
    let n_sets = reaction.synthons().len();
    let mut synthons_to_use: Vec<Vec<bool>> = reaction
        .synthons()
        .iter()
        .map(|s| vec![false; s.len()])
        .collect();
    let mut frag_sims: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_sets];

    // It makes sense to match fragments against synthon sets in order of
    // smallest synthon set first because if a fragment doesn't have a match
    // in a synthon set, the whole thing's a bust. So if frag_shapes[0] is
    // matched against 1000 synthons and then frag_shapes[1] is matched
    // against 10 synthons and doesn't match any of them, the first set of
    // matches was wasted time. This doesn't change the order of either, just
    // the order of comparison.
    let mut frag_orders: Vec<(usize, usize)> = synthon_set_order
        .iter()
        .enumerate()
        .map(|(i, &so)| (i, reaction.synthons()[so].len()))
        .collect();
    libcxx_sort_by(&mut frag_orders, |a, b| a.1 < b.1);
    for &(frag_num, _) in &frag_orders {
        let set_num = synthon_set_order[frag_num];
        let synthons = &reaction.synthons()[set_num];
        let frag_shape = searcher.frag_shapes_pool[frag_shapes[frag_num]]
            .as_ref()
            .expect("fragment shape");
        let mut frag_matched = false;
        // Because the combination score is the sum of 2 Tversky scores, it's
        // not possible to use the threshold to set upper and lower bounds on
        // the search space as is done with fingerprints and Rascal
        // similarity. We just need to plough through them in order.
        for (j, (_, synthon_idx)) in synthons.iter().enumerate() {
            let synthon = &pool[*synthon_idx];
            let Some(synth_shapes) = synthon.shapes() else {
                continue;
            };
            if searcher.has_precomputed_sims() {
                if let Some(sim) = searcher.frag_matched_synthon(frag_shapes[frag_num], synthon) {
                    synthons_to_use[set_num][j] = true;
                    frag_sims[set_num].push((j, sim.score));
                    frag_matched = true;
                }
            } else {
                let mut frag_cp = frag_shape.shapes().clone();
                let best = frag_cp.best_similarity(
                    synth_shapes.shapes(),
                    similarity_cutoff,
                    &ShapeOverlayOptions::default(),
                );
                if let Ok((csim, _)) = best
                    && csim[0] >= similarity_cutoff
                {
                    synthons_to_use[set_num][j] = true;
                    frag_sims[set_num].push((j, csim[0]));
                    frag_matched = true;
                }
            }
        }
        if !frag_matched {
            // No synthons matched this fragment, so the whole fragment set
            // is a bust.
            return Vec::new();
        }
    }

    // Fill in any synthons where they all didn't match because there were
    // fewer fragments than synthons.
    details::expand_bit_set(&mut synthons_to_use);
    let mut ret = details::bit_sets_to_vectors(&synthons_to_use);

    // Now order the synthons in descending order of their similarity to the
    // corresponding fragment.
    for i in 0..frag_shapes.len() {
        if frag_sims[i].is_empty() {
            // This one will have been filled in by expandBitSet so we need to
            // use all the synthons and a dummy similarity.
            frag_sims[i] = (0..synthons_to_use[i].len()).map(|j| (j, 0.0)).collect();
        } else {
            libcxx_sort_by(&mut frag_sims[i], |a, b| a.1 > b.1);
        }
        ret[i] = frag_sims[i].iter().map(|fs| fs.0).collect();
    }
    ret
}

/// Aligns `mol` conformer-for-conformer onto the query via the shape: the
/// volume of the hit that is inside the excluded volume.
fn calc_excluded_volume(
    mol: &ROMol,
    exc_vol: &ShapeInput,
    overlay_opts: &ShapeOverlayOptions,
) -> Result<f64> {
    let hit_shape = ShapeInput::new(mol, -1, &ShapeInputOptions::default(), overlay_opts)?;
    let (_, (a, b)) = score_shape(&hit_shape, exc_vol, overlay_opts);
    Ok(a + b)
}

fn calc_num_clashes(mol: &ROMol, exc_vol: &ShapeInput) -> Result<usize> {
    let double_c_sq = 4.0 * CARBON_RAD * CARBON_RAD;
    let shp_coords = exc_vol.coords();
    let conf = mol.conformer(-1)?;
    let mut clash_atoms = vec![false; mol.num_atoms() as usize];
    for atom in mol.atoms() {
        let a_pos = conf.atom_pos(atom.idx())?;
        for c in shp_coords.chunks_exact(3) {
            let s_pos = Point3D::new(c[0], c[1], c[2]);
            if (a_pos - s_pos).length_sq() < double_c_sq {
                clash_atoms[atom.idx() as usize] = true;
                break;
            }
        }
    }
    Ok(clash_atoms.iter().filter(|&&c| c).count())
}

/// Checks the hit against the excluded volume, if there is one. Returns
/// false if it fails one of the checks.
fn check_excluded_vols(
    mol: &ROMol,
    params: &SynthonSpaceSearchParams,
    excluded_vol: &mut f64,
    mean_excluded_vol: &mut f64,
) -> Result<bool> {
    let Some(exc) = &params.excluded_volume else {
        return Ok(true);
    };
    *excluded_vol = calc_excluded_volume(mol, exc, &params.shape_overlay_options)?;
    if params.max_excluded_volume > -0.5 && *excluded_vol > params.max_excluded_volume {
        return Ok(false);
    }
    let num_clashes = calc_num_clashes(mol, exc)?;
    if num_clashes > 0 {
        *mean_excluded_vol = *excluded_vol / num_clashes as f64;
        if params.max_mean_excluded_volume > -0.5
            && *mean_excluded_vol > params.max_mean_excluded_volume
        {
            return Ok(false);
        }
    } else {
        *mean_excluded_vol = 0.0;
    }
    Ok(true)
}

/// Sets the scores and name of the hit and the excluded volumes if there is
/// one.
fn set_hit_props(
    hit: &mut ROMol,
    scores: &[f64; 3],
    rxn_id: &str,
    synth_names: &[&str],
    params: &SynthonSpaceSearchParams,
    excluded_vol: f64,
    mean_excluded_vol: f64,
) -> Result<()> {
    hit.set_prop("Similarity", scores[0], false)?;
    hit.set_prop("ShapeScore", scores[1], false)?;
    hit.set_prop("ColorScore", scores[2], false)?;
    hit.set_prop(
        cp::_Name,
        details::build_product_name(rxn_id, synth_names),
        false,
    )?;
    if params.excluded_volume.is_some() {
        hit.set_prop("ExcludedVolume", excluded_vol, false)?;
        hit.set_prop("MeanExcludedVolume", mean_excluded_vol, false)?;
    }
    Ok(())
}

/// Assuming the hit shape is aligned onto the query, extracts that
/// conformation into `hit` and adds the scores to it as properties. Returns
/// false if it fails one of the excluded volume checks, true otherwise.
#[allow(clippy::too_many_arguments)]
fn finalise_hit(
    hit_shapes: &mut ShapeInput,
    hit_shape_num: usize,
    scores: &[f64; 3],
    rxn_id: &str,
    synth_names: &[&str],
    hit: &mut ROMol,
    params: &SynthonSpaceSearchParams,
    excluded_vol: &mut f64,
    mean_excluded_vol: &mut f64,
) -> Result<bool> {
    // Copy the conformer into the hit.
    hit_shapes.set_active_shape(hit_shape_num)?;
    let mut poss_hit = hit_shapes.shape_to_mol(false, true)?;
    if !check_excluded_vols(&poss_hit, params, excluded_vol, mean_excluded_vol)? {
        return Ok(false);
    }
    set_hit_props(
        &mut poss_hit,
        scores,
        rxn_id,
        synth_names,
        params,
        *excluded_vol,
        *mean_excluded_vol,
    )?;
    assign_stereochemistry_from_3d(&mut poss_hit, -1, true)?;
    *hit = poss_hit;
    Ok(true)
}

/// Makes a molecule from the synthon's shape, and transforms it to the
/// overlaid coords.
fn get_shape_mol(
    synthon: &Synthon,
    shape_num: usize,
    xform: Option<&Transform3D>,
) -> Result<RWMol> {
    let synth_shape = ShapeInput::from_shape(
        synthon.shapes().expect("synthon has shapes").shapes(),
        shape_num,
    )?;
    let mut shape_mol = synth_shape.shape_to_mol(false, true)?;
    if let Some(x) = xform {
        transform_conformer(shape_mol.conformer_mut(-1)?, x);
    }
    Ok(shape_mol)
}

fn check_bond_lengths(mol: &ROMol) -> Result<bool> {
    // DetermineBonds::connectivityVdw uses a covalent factor of 1.3.
    const RAD_FACTOR: f64 = 1.3;
    let conf = mol.conformer(-1)?;
    let pt = PeriodicTable::get();
    for bond in mol.bonds() {
        let (b, e) = (bond.begin_atom_idx(), bond.end_atom_idx());
        let (z1, z2) = (mol.atom(b).atomic_num(), mol.atom(e).atomic_num());
        if z1 == 0 || z2 == 0 {
            continue;
        }
        let bond_len = get_bond_length(conf, b, e)?;
        if bond_len > RAD_FACTOR * (pt.get_rcovalent(z1)? + pt.get_rcovalent(z2)?) {
            return Ok(false);
        }
    }
    Ok(true)
}

impl<'a> ShapeSearcher<'a> {
    fn shape_hit_data<'h>(hitset: &'h SynthonSpaceHitSet<'a>) -> &'h super::hit_set::ShapeHitData {
        hitset.shape.as_ref().expect("shape hit set")
    }

    fn verify_hit_impl(&self, hit: &mut ROMol, rxn_id: &str, synth_names: &[&str]) -> Result<bool> {
        let params = &self.base.params;
        let query_shape = self.query_shape();
        let mut init_scores = [-1.0; 3];
        let mut init_hit: Option<ROMol> = None;
        if hit.num_conformers() > 0 {
            init_scores = align_molecule_to_shape(
                query_shape.shapes(),
                hit,
                &ShapeInputOptions::default(),
                &params.shape_overlay_options,
                -1,
            )?
            .0;
            let (mut excluded_vol, mut mean_excluded_vol) = (-1.0, -1.0);
            if check_excluded_vols(hit, params, &mut excluded_vol, &mut mean_excluded_vol)? {
                // If the excluded vol is also ok, take this hit. Otherwise,
                // pass it through to conformation expansion to see if it wins
                // there.
                set_hit_props(
                    hit,
                    &init_scores,
                    rxn_id,
                    synth_names,
                    params,
                    excluded_vol,
                    mean_excluded_vol,
                )?;
                if !params.best_hit
                    && check_bond_lengths(hit)?
                    && init_scores[0] >= params.similarity_cutoff
                {
                    self.base.update_best_hit_so_far(hit, init_scores[0]);
                    return Ok(true);
                }
                if init_scores[0] >= params.similarity_cutoff {
                    init_hit = Some(hit.clone());
                }
            } else {
                init_scores[0] = -1.0; // To flag it as a bad hit.
            }
        }
        // If the run is multi-threaded, this will already be running on the
        // maximum number of threads, so do the embedding on a single thread.
        let mut dg_params = EmbedParameters::etkdg_v3();
        dg_params.prune_rms_thresh = params.conf_rms_threshold;
        dg_params.random_seed = params.random_seed;
        dg_params.timeout = params.time_out as u32;
        let hit_confs = generate_isomer_conformers(
            hit,
            params.num_conformers,
            true,
            &params.stereo_enum_opts,
            &dg_params,
            params.user_conformer_generator.as_ref(),
            u32::MAX,
        )?;
        let mut found_hit = false;
        // Prefer a hit from a generated conf rather than the initial hit
        // which might have iffy bond lengths because it's the output from
        // molzip.
        let mut best_sim = params.similarity_cutoff;
        let overlay_options = &params.shape_overlay_options;
        // Don't prune the hit shapes, it just wastes time.
        let opts = ShapeInputOptions {
            shape_prune_threshold: -1.0,
            ..Default::default()
        };
        for isomer in &hit_confs {
            let mut isomer_shapes = ShapeInput::new(isomer, -1, &opts, overlay_options)?;
            // Do all the shapes one by one rather than via bestSimilarity
            // because excluded volumes might be too high for the best scoring
            // overlay, but there might be an acceptable one with a lower
            // score.
            for shape_num in 0..isomer_shapes.num_shapes() {
                let (mut excluded_vol, mut mean_excluded_vol) = (-1.0, -1.0);
                let mut finalised_hit = false;
                isomer_shapes.set_active_shape(shape_num)?;
                let (this_score, _) =
                    align_shape(query_shape.shapes(), &mut isomer_shapes, overlay_options)?;
                if this_score[0] > self.base.best_similarity_so_far() {
                    if !finalise_hit(
                        &mut isomer_shapes,
                        shape_num,
                        &this_score,
                        rxn_id,
                        synth_names,
                        hit,
                        params,
                        &mut excluded_vol,
                        &mut mean_excluded_vol,
                    )? {
                        // It failed excluded vols so don't use it.
                        continue;
                    }
                    finalised_hit = true;
                    self.base.update_best_hit_so_far(hit, this_score[0]);
                }
                if this_score[0] >= best_sim && this_score[0] > params.similarity_cutoff {
                    if !finalised_hit
                        && !finalise_hit(
                            &mut isomer_shapes,
                            shape_num,
                            &this_score,
                            rxn_id,
                            synth_names,
                            hit,
                            params,
                            &mut excluded_vol,
                            &mut mean_excluded_vol,
                        )?
                    {
                        continue;
                    }
                    // If we're only interested in whether there's a shape
                    // match, and not in finding the best shape, we're done.
                    if !params.best_hit {
                        return Ok(true);
                    }
                    found_hit = true;
                    best_sim = this_score[0];
                }
            }
        }
        if !found_hit
            && let Some(init) = init_hit
            && init_scores[0] >= params.similarity_cutoff
        {
            // Stick with what we found at the start because the conformational
            // expansion didn't do any better.
            *hit = init;
            self.base.update_best_hit_so_far(hit, init_scores[0]);
            found_hit = true;
        }
        Ok(found_hit)
    }
}

impl<'a> SynthonSearcher<'a> for ShapeSearcher<'a> {
    fn base(&self) -> &SearcherBase<'a> {
        &self.base
    }

    fn search_frag_set(
        &self,
        frag_set: &[Arc<ROMol>],
        reaction: &'a SynthonSet,
    ) -> Result<Vec<SynthonSpaceHitSet<'a>>> {
        let mut results = Vec::new();
        if frag_set.len() > reaction.synthons().len() {
            return Ok(results);
        }
        // Collect the ShapeSets for the fragSet
        let Some(frag_shapes) = frag_set
            .iter()
            .map(|f| self.frag_shape_of(f))
            .collect::<Option<Vec<usize>>>()
        else {
            return Ok(results);
        };
        let params = &self.base.params;
        let pool = self.base.space.pool();
        // Need to try all combinations of synthon orders.
        for synthon_order in details::perm_m_from_n(frag_set.len(), reaction.synthons().len()) {
            let these = get_hit_synthons(
                self,
                &frag_shapes,
                params.similarity_cutoff - params.frag_similarity_adjuster,
                reaction,
                &synthon_order,
            );
            if !these.is_empty() {
                let hs = SynthonSpaceHitSet::new_shape(
                    reaction,
                    pool,
                    &these,
                    frag_set,
                    frag_shapes.clone(),
                    synthon_order,
                );
                if hs.num_hits > 0 {
                    results.push(hs);
                }
            }
        }
        Ok(results)
    }

    fn approx_similarity(&self, hitset: &SynthonSpaceHitSet<'a>, synth_nums: &[usize]) -> f64 {
        let params = &self.base.params;
        let mut approx_sim = 0.0;
        if self.has_precomputed_sims() {
            let hs = Self::shape_hit_data(hitset);
            let sso = &hs.synthon_set_order;
            // Take the approximate similarity as the mean of the individual
            // frag->synthon similarities weighted by the synthon volumes.
            let mut tot_vol = 0.0;
            let mut synth_vols = vec![0.0; synth_nums.len()];
            let mut scores = vec![0.0; synth_nums.len()];
            let m = params.shape_overlay_options.opt_param;
            // Not all fragShapes will have a matching synthon. For example,
            // if a 2 fragment split matched 2 synthons from a 3 synthon
            // SynthonSet. All synthons of the 3rd set will have been included
            // as potential hits. Do the overlap vol calcs with the fragments
            // to matching synthon.
            let mut synths_matched = vec![false; synth_nums.len()];
            for (i, &frag_shape) in hs.frag_shapes.iter().enumerate() {
                let synthon = hitset.synthons_to_use[sso[i]][synth_nums[sso[i]]].1;
                // If one of the synthons doesn't have a shape, the hit is
                // clearly a bust.
                let Some(shapes) = synthon.shapes() else {
                    return 0.0;
                };
                // Look up the similarity for this shape/synthon combination.
                let Some(sim) = self.frag_matched_synthon(frag_shape, synthon) else {
                    return 0.0;
                };
                // Get the volume of the synthon for the shape conformer that
                // gave the similarity values.
                let (Ok(shape_vol), Ok(dummy_vol), Ok(colour_vol)) = (
                    shapes.shapes().shape_volume_of(sim.shape_num),
                    shapes.dummy_volume(sim.shape_num),
                    shapes.shapes().color_volume_of(sim.shape_num),
                ) else {
                    return 0.0;
                };
                synth_vols[sso[i]] = (1.0 - m) * (shape_vol - dummy_vol) + m * colour_vol;
                scores[sso[i]] = sim.score;
                tot_vol += synth_vols[sso[i]];
                synths_matched[sso[i]] = true;
            }
            for i in 0..synth_nums.len() {
                if synths_matched[i] {
                    continue;
                }
                let synthon = hitset.synthons_to_use[i][synth_nums[i]].1;
                // Use the first shape, which should have the largest volume.
                let Some(shapes) = synthon.shapes() else {
                    return 0.0;
                };
                let (Ok(shape_vol), Ok(dummy_vol), Ok(colour_vol)) = (
                    shapes.shapes().shape_volume_of(0),
                    shapes.dummy_volume(0),
                    shapes.shapes().color_volume_of(0),
                ) else {
                    return 0.0;
                };
                synth_vols[i] = (1.0 - m) * (shape_vol - dummy_vol) + m * colour_vol;
                tot_vol += synth_vols[i];
            }
            // scores[i] will be zero for un-matched synthons.
            for i in 0..synth_nums.len() {
                approx_sim += scores[i] * synth_vols[i] / tot_vol;
            }
        } else {
            // This is the best we can do without pre-computed similarities
            let mut max_vol = 0.0;
            let mut feature_vol = 0.0;
            // The synthon shapes are sorted in descending order of sov + sof.
            // Assume therefore that the maximum volume of the hit is the sum
            // of the sov's of the first shape in each synthon, minus the
            // volume of their dummy atoms.
            for (i, &n) in synth_nums.iter().enumerate() {
                let synth = hitset.synthons_to_use[i][n].1;
                let Some(shapes) = synth.shapes() else {
                    return 0.0;
                };
                let (Ok(shape_vol), Ok(dummy_vol), Ok(colour_vol)) = (
                    shapes.shapes().shape_volume_of(0),
                    shapes.dummy_volume(0),
                    shapes.shapes().color_volume_of(0),
                ) else {
                    return 0.0;
                };
                max_vol += shape_vol - dummy_vol;
                feature_vol += colour_vol;
            }
            let query = self.query_shape().shapes();
            let (Ok(q_shape_vol), Ok(q_colour_vol)) =
                (query.shape_volume_of(0), query.color_volume_of(0))
            else {
                return 0.0;
            };
            approx_sim = max_score(
                max_vol,
                q_shape_vol,
                feature_vol,
                q_colour_vol,
                &params.shape_overlay_options,
            );
        }
        approx_sim
    }

    /// Builds the hit, doing its best to line the synthons up correctly. It
    /// may not do a great job, however for at least 2 reasons. The synthon
    /// may not match a fragment, in which case there is no good way of lining
    /// it up on the query. Also, molzip lines the fragments up using the bond
    /// vectors to the dummy atoms, but this fails if it's building a ring. The
    /// product geometry should be checked for bad bond lengths.
    fn build_hit(
        &self,
        hitset: &SynthonSpaceHitSet<'a>,
        synth_nums: &[usize],
        synth_names: &mut Vec<&'a str>,
    ) -> Option<ROMol> {
        let hs = Self::shape_hit_data(hitset);
        let sso = &hs.synthon_set_order;
        // shape_synths are the molecules produced by the shapes, with coords.
        let mut shape_synths: Vec<Option<ROMol>> = vec![None; synth_nums.len()];
        let mut synths_matched = vec![false; synth_nums.len()];
        // We need to build all the synthons into a molecule, with them put
        // into the position the fragment->shape alignment gave. Not all
        // fragments will match a synthon, such as when the reaction had 3
        // synthon sets and there were only 2 fragments, and both those
        // fragments had a synthon that matched. The third synthon will have
        // been filled in from all possibilities.
        for (i, &frag_shape) in hs.frag_shapes.iter().enumerate() {
            let (name, synthon) = hitset.synthons_to_use[sso[i]][synth_nums[sso[i]]];
            // If one of the synthons doesn't have a shape, the hit is clearly
            // a bust.
            synthon.shapes()?;
            // Look up the similarity for this shape/synthon combination.
            if let Some(sim) = self.frag_matched_synthon(frag_shape, synthon) {
                // We need to build a copy of the synthon and give it the
                // coords of the associated shape that is similar to the
                // fragment, and transform it to the overlaid coords.
                let shape_mol = get_shape_mol(synthon, sim.shape_num, sim.xform.as_deref()).ok()?;
                shape_synths[sso[i]] = Some(shape_mol);
                synth_names[sso[i]] = name;
                synths_matched[sso[i]] = true;
            }
        }
        for i in 0..synth_nums.len() {
            if !synths_matched[i] {
                let (name, synthon) = hitset.synthons_to_use[i][synth_nums[i]];
                shape_synths[i] = Some(get_shape_mol(synthon, 0, None).ok()?);
                synth_names[i] = name;
            }
        }
        let plain_synths: Vec<&ROMol> = (0..synth_nums.len())
            .map(|i| hitset.synthons_to_use[i][synth_nums[i]].1.orig_mol())
            .collect();
        let shape_refs: Vec<&ROMol> = shape_synths.iter().flatten().collect();
        let mut prod = details::build_product(&shape_refs);
        if prod.is_none() {
            // Try it with the plain synthons, which will obviously need
            // coordinates creating later on.
            prod = details::build_product(&plain_synths);
            if prod.is_none() {
                println!("Failed building ");
                for s in synth_names.iter() {
                    println!("{s} ");
                }
                println!(" for reaction {}", hitset.reaction.id());
            }
        }
        prod
    }

    fn extra_search_setup(
        &mut self,
        frag_sets: &mut Vec<Vec<Arc<ROMol>>>,
        end_time: Option<Instant>,
    ) -> Result<bool> {
        // Use the given conformer unless it looks like a 2D molecule.
        let query = self.base.query.clone();
        if query.num_conformers() == 0 || !query.conformer(-1)?.is_3d() {
            crate::rd_log_error!("The query molecule needs a 3D conformer.");
            return Ok(false);
        }
        if query.num_conformers() > 1 {
            crate::rd_log_warning!(
                "The query molecule has multiple conformers.  Just using the default."
            );
        }
        let all_features = self.build_query_shape(&query)?;

        // Organise the fragments into a map where each entry is a set of
        // fragments from the same atoms in the query. SMILES strings aren't
        // good enough for this, as they are for the other synthon searches,
        // because 2 fragments can have the same SMILES but be of different
        // atoms and therefore potentially different shapes. The key is an
        // underscore-separated list of the atom indices.
        let frag_atoms_to_frag = map_frags_by_atoms(frag_sets);

        // Compute ShapeSets for the fragments
        let frags_for_shape: Vec<&ROMol> = frag_atoms_to_frag
            .values()
            .map(|members| &*frag_sets[members[0].0][members[0].1])
            .collect();
        let counter = std::sync::atomic::AtomicUsize::new(0);
        let timed_out = std::sync::atomic::AtomicBool::new(false);
        let shapes = par_map(&frags_for_shape, self.base.params.num_threads, |frag| {
            use std::sync::atomic::Ordering::Relaxed;
            if timed_out.load(Relaxed) {
                return Ok(None);
            }
            if counter.fetch_add(1, Relaxed) % 100 == 0 && details::check_time_out(end_time) {
                timed_out.store(true, Relaxed);
                return Ok(None);
            }
            generate_shapes(&query, frag, &all_features).map(Some)
        });
        self.frag_shapes_pool = shapes.into_iter().collect::<Result<Vec<_>>>()?;

        // Keep track of the minimum fragSet size that each frag is in, which
        // will be used in computeFragSynthonSims.
        let mut min_frag_set_size: HashMap<usize, usize> = HashMap::new();
        for frag_set in frag_sets.iter() {
            for frag in frag_set {
                let e = min_frag_set_size
                    .entry(frag_key(frag))
                    .or_insert(frag_set.len());
                if *e < frag_set.len() {
                    *e = frag_set.len();
                }
            }
        }
        // Use the pooled ShapeSets to populate the map for each fragment.
        // Use the smallest min_frag_set_size as the minShapeSetSize as we are
        // only doing 1 shape for all the identical frags. frag_shapes will
        // hold the number of the shape in the pool, keyed on the
        // corresponding fragment.
        self.frag_shapes.clear();
        let mut min_shape_set_size = vec![10usize; self.frag_shapes_pool.len()];
        for (frag_num, members) in frag_atoms_to_frag.values().enumerate() {
            for &(i, j) in members {
                let key = frag_key(&frag_sets[i][j]);
                min_shape_set_size[frag_num] =
                    min_shape_set_size[frag_num].min(min_frag_set_size[&key]);
                self.frag_shapes.insert(key, frag_num);
            }
        }

        // It's likely that a fragment occurs more than once in the fragment
        // sets, so compute all the fragment->synthon similarities up front.
        // This makes maximum use of the parallel environment, doesn't do
        // anything that wouldn't be done at some point and may prevent
        // duplicated comparisons.
        self.compute_frag_synthon_sims(end_time, &min_shape_set_size)?;
        Ok(true)
    }

    fn verify_hit(&self, hit: &mut ROMol, rxn_id: &str, synth_names: &[&str]) -> bool {
        match self.verify_hit_impl(hit, rxn_id, synth_names) {
            Ok(b) => b,
            Err(e) => {
                crate::rd_log_error!("{e}");
                false
            }
        }
    }

    fn process_to_try_set(
        &self,
        to_try: &mut ToTry<'_, 'a>,
        end_time: Option<Instant>,
        results: &mut Vec<ROMol>,
        num_hits_found: &std::sync::atomic::AtomicI64,
        num_poss_hits_written: &mut u64,
    ) -> Result<()> {
        // This needs to be done a bit differently because the approximate
        // similarities depend on the fragments as well as the synthons
        // selected. For the other searches, such as fingerprints, the
        // synthons are compared with the full query so it doesn't matter
        // which entry in the hit set is used. If we uniquify on
        // reaction/synthons straight away, we may end up with a
        // reaction/synthon pair that doesn't hit the threshold even though
        // one of the other hitsets might have the same reaction/synthon pair
        // coupled with a different fragSet that does go over the threshold.
        // First calculate the approximate similarities and discard any that
        // don't meet the threshold.
        let params = &self.base.params;
        let mut approx_sims: Vec<(usize, f64)> = Vec::with_capacity(to_try.len());
        for (i, (hitset, nums)) in to_try.iter().enumerate() {
            let sim = self.approx_similarity(hitset, nums);
            if sim >= params.similarity_cutoff - params.approx_similarity_adjuster {
                approx_sims.push((i, sim));
            }
        }
        if approx_sims.is_empty() {
            return Ok(());
        }
        libcxx_sort_by(&mut approx_sims, |a, b| a.1 > b.1);
        let new_to_try: ToTry<'_, 'a> = approx_sims.iter().map(|a| to_try[a.0].clone()).collect();
        *to_try = new_to_try;
        details::sort_and_uniquify_to_try(to_try);
        make_hits_from_to_try(
            self,
            to_try,
            end_time,
            results,
            num_hits_found,
            num_poss_hits_written,
        )
    }
}

/// Join the fragments by a zero-order bond between the 2 nearest atoms. This
/// is intended only for the query.
fn add_zero_order_bonds_between_frags(query: &ROMol, frags: &[Vec<i32>]) -> Result<ROMol> {
    if frags.len() != 2 {
        return Err(RDKitError::precondition(
            "The query can have at most 2 fragments.",
        ));
    }
    let mut ret = query.clone();
    let conf = query.conformer(-1)?;
    let mut nearest = (0u32, 0u32);
    let mut nearest_dist = 1.0e6;
    for &a0 in &frags[0] {
        let pos0 = conf.atom_pos(a0 as u32)?;
        for &a1 in &frags[1] {
            let d = (pos0 - conf.atom_pos(a1 as u32)?).length_sq();
            if d < nearest_dist {
                nearest_dist = d;
                nearest = (a0 as u32, a1 as u32);
            }
        }
    }
    ret.add_bond(nearest.0, nearest.1, BondType::Zero)?;
    crate::graphmol::mol_ops::sanitize_mol(
        &mut ret,
        crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL,
    )
    .map_err(|(_, e)| e)?;
    Ok(ret)
}

/// If there's 2 fragments in the query, join them together with a zero-order
/// bond between the closest atoms in each. It makes the fragmentation work,
/// but otherwise doesn't affect the final outcome. More than 2 fragments
/// aren't allowed for now.
fn join_query_fragments(query: &ROMol) -> Result<ROMol> {
    let frags = get_mol_frags(query);
    if frags.len() > 1 {
        add_zero_order_bonds_between_frags(query, &frags)
    } else {
        Ok(query.clone())
    }
}

impl SynthonSpace {
    /// `shapeSearch`
    pub fn shape_search(
        &self,
        query: &ROMol,
        params: &SynthonSpaceSearchParams,
    ) -> Result<SearchResults> {
        check_query(query)?;
        let query_cp = join_query_fragments(query)?;
        let mut s = ShapeSearcher::new(&query_cp, params, self, true)?;
        search(&mut s, ThreadMode::ThreadFragments)
    }

    /// `shapeSearch` returning the hits to a callback in chunks.
    pub fn shape_search_incremental(
        &self,
        query: &ROMol,
        cb: SearchResultCallback<'_>,
        params: &SynthonSpaceSearchParams,
    ) -> Result<()> {
        check_query(query)?;
        let query_cp = join_query_fragments(query)?;
        let mut s = ShapeSearcher::new(&query_cp, params, self, true)?;
        search_with_callback(&mut s, cb, ThreadMode::ThreadFragments)
    }

    /// `shapeSearch` on the possible hits between two lines of the possible
    /// hits file.
    pub fn shape_search_lines(
        &self,
        query: &ROMol,
        params: &SynthonSpaceSearchParams,
        start_line: u64,
        finish_line: u64,
    ) -> Result<SearchResults> {
        check_query(query)?;
        check_lines(start_line, finish_line)?;
        let empty = SynthonSpace::new();
        let s = ShapeSearcher::new(query, params, &empty, false)?;
        check_possible_hits(&s, start_line, finish_line)
    }
}
