//! Building the shapes of the synthons: `trimSampleMol`,
//! `generateIsomerConformers` and `makeShapesFromMol(s)` of
//! `SynthonSpaceSearch_details.cpp`, and `SynthonSpace::buildSynthonShapes`
//! of `SynthonSpace.cpp`.

use super::params::{ShapeBuildParams, UserConfGenerator};
use super::shape_input::SynthonShapeInput;
use super::space::SynthonSpace;
use super::synthon::Synthon;
use super::synthon_set::{SampleMolRec, SynthonSet};
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::BondType;
use crate::graphmol::chem_transforms::mol_fragmenter::fragment_on_bonds;
use crate::graphmol::dist_geom_helpers::{EmbedParameters, embed_multiple_confs};
use crate::graphmol::enumerate_stereoisomers::{StereoEnumerationOptions, StereoisomerEnumerator};
use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, get_mol_frags, symmetrize_sssr};
use crate::graphmol::gaussian_shape::{
    DUMMY_RAD, ShapeInputOptions, ShapeOverlayOptions, find_features,
};
use crate::graphmol::mol_ops::add_hs::{AddHsParameters, add_hs};
use crate::graphmol::mol_ops::{RemoveHsParameters, remove_hs};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::cx_write::{RestoreBondDirOption, cx_fields, mol_to_cx_smiles};
use crate::graphmol::smiles_parse::write::SmilesWriteParams;
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::rdgeneral::dict::Dict;
use crate::rdgeneral::threads::par_map;
use crate::{RDKitError, Result};

/// The most stereo centres that are enumerated when building the shapes.
const MAX_STEREO_CENTRES_TO_DO: u32 = 3;

fn cx_smiles(mol: &ROMol) -> String {
    mol_to_cx_smiles(
        mol,
        &SmilesWriteParams::default(),
        cx_fields::ALL,
        RestoreBondDirOption::Clear,
    )
    .unwrap_or_default()
}

fn mol_num_of(props: &impl RDProps) -> Option<u32> {
    props.get_prop("molNum").and_then(|v| v.as_uint())
}

/// Carries the sample molecule properties from an atom or bond to the
/// equivalent in the conformer generator's output.
fn carry_props(old: &Dict, new: &mut Dict) {
    for key in ["molNum", "idx", "HoldsJoin"] {
        if let Some(v) = old.get(key) {
            new.set_val(key, v.clone());
        }
    }
}

/// `useUserConfGen`: the molecule is replaced by the generator's output with
/// the sample molecule properties copied over. `None` if that failed.
fn use_user_conf_gen(
    mol: &ROMol,
    num_conformers: u32,
    generator: &UserConfGenerator,
) -> Result<Option<ROMol>> {
    let tmp = mol.clone();
    let smiles = super::details::smiles_of(&tmp);
    let Some(mut new_mol) = generator(&smiles, num_conformers) else {
        return Ok(None);
    };
    // The mol has information that won't have survived, and it's not safe to
    // assume the atom orders have survived.
    let ps = SubstructMatchParameters {
        uniquify: false,
        max_matches: 1,
        ..Default::default()
    };
    let matches = substruct_match(&new_mol, mol, &ps)?;
    let Some(res) = matches.first() else {
        // Bit weird!
        crate::rd_log_error!(
            "Couldn't generate conformers for isomer {} of molecule.",
            cx_smiles(mol)
        );
        return Ok(None);
    };
    // query (mol) atom -> new_mol atom
    let query2mol: std::collections::BTreeMap<u32, u32> =
        res.iter().map(|&(q, m)| (q as u32, m as u32)).collect();
    for (&curr, &new) in &query2mol {
        carry_props(mol.atom(curr).dict(), new_mol.atom_mut(new).dict_mut());
        for nbr in mol.atom_neighbors(curr) {
            let Some(&new_nbr) = query2mol.get(&nbr) else {
                continue;
            };
            let (Some(old_bond), Some(new_bond)) = (
                mol.get_bond_between_atoms(curr, nbr)?,
                new_mol
                    .get_bond_between_atoms(new, new_nbr)?
                    .map(|b| b.idx()),
            ) else {
                continue;
            };
            carry_props(old_bond.dict(), new_mol.bond_mut(new_bond).dict_mut());
        }
    }
    if let Some(name) = mol.get_prop(cp::_Name) {
        new_mol.set_prop(cp::_Name, name.clone(), false)?;
    }
    Ok(Some(new_mol))
}

/// `useDefaultConfGen`: ETKDG on the molecule with Hs, which are removed
/// again. `None` if no conformer was made.
fn use_default_conf_gen(
    mol: &ROMol,
    num_conformers: u32,
    dg_params: &EmbedParameters,
) -> Result<Option<ROMol>> {
    let mut mol = mol.clone();
    add_hs(&mut mol, &AddHsParameters::default(), None)?;
    let mut params = dg_params.clone();
    let cids = embed_multiple_confs(&mut mol, num_conformers, &mut params)?;
    if cids.is_empty() || mol.num_conformers() == 0 {
        crate::rd_log_error!(
            "Couldn't generate conformers for isomer {} of molecule.",
            cx_smiles(&mol)
        );
        return Ok(None);
    }
    remove_hs(&mut mol, &RemoveHsParameters::default(), true)?;
    Ok(Some(mol))
}

/// `generateIsomerConformers`: conformers for the molecule, including
/// enumerating stereoisomers if requested (at most `max_stereo_centres` of
/// them). The molecules that no conformers could be made for are dropped.
pub fn generate_isomer_conformers(
    mol: &ROMol,
    num_conformers: u32,
    enumerate_stereo: bool,
    enum_opts: &StereoEnumerationOptions,
    dg_params: &EmbedParameters,
    user_conf_generator: Option<&UserConfGenerator>,
    max_stereo_centres: u32,
) -> Result<Vec<ROMol>> {
    let mut isomers = Vec::new();
    if enumerate_stereo {
        let mut enu = StereoisomerEnumerator::new(mol, enum_opts.clone())?;
        let mut i = 0;
        while let Some(isomer) = enu.next()? {
            isomers.push(isomer);
            i += 1;
            if i == max_stereo_centres {
                break;
            }
        }
    } else {
        isomers.push(mol.clone());
    }
    let mut conf_mols = Vec::with_capacity(isomers.len());
    for isomer in isomers {
        let made = match user_conf_generator {
            Some(g) => use_user_conf_gen(&isomer, num_conformers, g)?,
            None => use_default_conf_gen(&isomer, num_conformers, dg_params)?,
        };
        if let Some(m) = made
            && m.num_conformers() > 0
        {
            conf_mols.push(m);
        }
    }
    Ok(conf_mols)
}

/// `trimSampleMol`: trims the molecule down to the minimum needed for a
/// sensible conformation expansion of the atoms with property `molNum`.
pub fn trim_sample_mol(mol: &mut ROMol, mol_num: usize) -> Result<RWMol> {
    let mut mol_num_atoms = vec![false; mol.num_atoms() as usize];
    for atom in mol.atoms() {
        if mol_num_of(atom).is_some_and(|n| n as usize == mol_num) {
            mol_num_atoms[atom.idx() as usize] = true;
        }
    }
    if !mol.ring_info().is_initialized() {
        symmetrize_sssr(mol, SymmetrizeSSSRAlgorithm::Default, true, false, false)?;
    }
    for ring in mol.ring_info().atom_rings() {
        if ring.iter().any(|&a| mol_num_atoms[a as usize]) {
            for &ap in ring {
                mol_num_atoms[ap as usize] = true;
            }
        }
    }
    let mut bonds_to_go = Vec::new();
    for bond in mol.bonds() {
        // Remove single non-cyclic bonds that aren't in mol_num_atoms at
        // either end.
        if bond.bond_type() != BondType::Single || mol.ring_info().num_bond_rings(bond.idx())? != 0
        {
            continue;
        }
        let (beg, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
        let mut bond_to_go = !(mol_num_atoms[beg as usize] || mol_num_atoms[end as usize]);
        // Don't break a S(=O)-[C,N] bond as that gives an odd-looking [SH](=O)
        // product.
        let mut s_atom = None;
        if bond_to_go && mol.atom(beg).atomic_num() == 16 {
            s_atom = Some((beg, end));
        }
        if bond_to_go && mol.atom(end).atomic_num() == 16 {
            s_atom = Some((end, beg));
        }
        if let Some((s, other)) = s_atom
            && matches!(mol.atom(other).atomic_num(), 6 | 7)
        {
            let num_oxy = mol
                .atom_neighbors(s)
                .filter(|&nbr| {
                    mol.atom(nbr).atomic_num() == 8
                        && mol
                            .get_bond_between_atoms(s, nbr)
                            .is_ok_and(|b| b.is_some())
                })
                .count();
            if num_oxy > 0 {
                bond_to_go = false;
            }
        }
        if bond_to_go {
            bonds_to_go.push(bond.idx());
        }
    }
    let mut new_mol = fragment_on_bonds(mol, &bonds_to_go, false, None, None, None)?;
    let mol_frags = get_mol_frags(&new_mol);
    new_mol.begin_batch_edit()?;
    for frag in &mol_frags {
        let hose_frag = !frag
            .iter()
            .any(|&fa| mol_num_of(new_mol.atom(fa as u32)).is_some_and(|n| n as usize == mol_num));
        if hose_frag {
            for &fa in frag {
                new_mol.remove_atom(fa as u32)?;
            }
        }
    }
    new_mol.commit_batch_edit()?;
    Ok(new_mol)
}

/// For deciphering the contents of the HoldsJoin property on 2 atoms. Used
/// to decide the isotope number a dummy atom must be set to.
fn join_in_common(joins1: &str, joins2: &str) -> u32 {
    // Each of the join strings will be a underscore-separated list of ints,
    // with no more than 3 ints per string. The ints will be between 1 and 4.
    // We want the single number in common to both lists.
    for c1 in joins1.chars().filter(|&c| c != '_') {
        if joins2.contains(c1) {
            return c1.to_digit(10).unwrap_or(0);
        }
    }
    0
}

fn holds_join(mol: &ROMol, atom: u32) -> String {
    match mol.atom(atom).get_prop("HoldsJoin") {
        Some(crate::rdgeneral::RDValue::String(s)) => s.clone(),
        _ => String::new(),
    }
}

/// Sometimes the atom subsetting leaves bonds between atoms not in molNum.
/// Once the dummy isotopes have been set, they will be between dummy atoms
/// and so can safely be removed.
pub(super) fn split_dummy_dummy_bonds(mol: &mut RWMol) -> Result<()> {
    let to_go: Vec<(u32, u32)> = mol
        .bonds()
        .iter()
        .filter(|b| {
            mol.atom(b.begin_atom_idx()).atomic_num() == 0
                && mol.atom(b.end_atom_idx()).atomic_num() == 0
        })
        .map(|b| (b.begin_atom_idx(), b.end_atom_idx()))
        .collect();
    mol.begin_batch_edit()?;
    for (b, e) in to_go {
        mol.remove_bond(b, e)?;
    }
    mol.commit_batch_edit()
}

/// (atom made into a dummy, the atom it is attached to, the isotope number)
type AtomIsotope = (u32, u32, u32);

fn set_join_isotope(
    mol: &mut RWMol,
    atom_isotopes: &[AtomIsotope],
    dummy_radii: &mut Vec<(u32, f64)>,
) {
    for &(idx, _, iso) in atom_isotopes {
        let atom = mol.atom_mut(idx);
        atom.set_atomic_num(0);
        atom.set_isotope(iso);
        atom.set_num_explicit_hs(0);
        dummy_radii.push((idx, DUMMY_RAD));
    }
}

/// Deals with something of an edge case. The synthon was most likely a
/// single atom that finished a ring, so had 2 connection points e.g.
/// `R1-[14C]([1*])[2*]` and `R2C([1*])CCC[2*]` to give `R1[14C]1C(R2)CCCC2`.
/// The shape needs to have 2 dummy atoms, so not `R1[*]1C(R2)CCCC2` but
/// `R1[*]CCC[*]R2` with the 2 dummy atoms at the same coordinates. This needs
/// to be done last thing before the shapes are created as it might break an
/// aromatic system that will cause problems.
fn duplicate_join_isotope(
    mol: &mut RWMol,
    atom_isotopes: &[AtomIsotope],
    dummy_radii: &mut Vec<(u32, f64)>,
    frag_atoms: &mut Vec<u32>,
) -> Result<()> {
    for &(idx, other, iso) in atom_isotopes {
        if mol.atom(idx).atomic_num() != 0 || mol.atom(idx).isotope() == iso {
            continue;
        }
        let mut dummy = Atom::new(0);
        dummy.set_isotope(iso);
        let new_idx = mol.add_atom(dummy, true);
        // Remove the bond between this atom and the other atom and make a
        // bond between the new atom and the other atom.
        let bt = mol
            .get_bond_between_atoms(idx, other)?
            .ok_or_else(|| RDKitError::Value("missing join bond".into()))?
            .bond_type();
        mol.remove_bond(idx, other)?;
        mol.add_bond(other, new_idx, bt)?;
        for conf in mol.conformers_mut() {
            let pos = conf.atom_pos(idx)?;
            conf.set_atom_pos(new_idx, pos)?;
        }
        frag_atoms.push(new_idx);
        dummy_radii.push((new_idx, DUMMY_RAD));
    }
    Ok(())
}

/// The shapes for an isomer of the sample molecule, for the atoms of the
/// synthon with `synthon_set_num`.
fn make_isomer_shapes(
    isomer: &mut RWMol,
    synthon_set_num: u32,
    shape_build_params: &ShapeBuildParams,
) -> Result<Option<SynthonShapeInput>> {
    let mut frag_atoms: Vec<u32> = Vec::new();
    let mut dummy_radii: Vec<(u32, f64)> = Vec::new();
    let mut atom_isotopes: Vec<AtomIsotope> = Vec::new();
    for bond in isomer.bonds() {
        let (beg, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
        match mol_num_of(bond) {
            None => {
                let (Some(beg_num), Some(end_num)) =
                    (mol_num_of(isomer.atom(beg)), mol_num_of(isomer.atom(end)))
                else {
                    // Probably an explicit H added by the stereoisomer
                    // enumerator.
                    continue;
                };
                if (beg_num == synthon_set_num) == (end_num == synthon_set_num) {
                    continue;
                }
                frag_atoms.push(beg);
                frag_atoms.push(end);
                // Set the atom as a dummy, so it shows up in the shape.
                let join_isotope =
                    join_in_common(&holds_join(isomer, end), &holds_join(isomer, beg));
                if beg_num == synthon_set_num {
                    atom_isotopes.push((end, beg, join_isotope));
                } else {
                    atom_isotopes.push((beg, end, join_isotope));
                }
            }
            Some(n) if n == synthon_set_num => {
                frag_atoms.push(beg);
                frag_atoms.push(end);
            }
            Some(_) => {}
        }
    }
    set_join_isotope(isomer, &atom_isotopes, &mut dummy_radii);
    frag_atoms.sort_unstable();
    frag_atoms.dedup();
    let mut shape_opts = ShapeInputOptions {
        shape_prune_threshold: shape_build_params.shape_sim_threshold,
        ..Default::default()
    };
    // Make custom features for the atom subset. Doing it from the parent
    // molecule means that the atom types will be correct which may not be the
    // case after the fragment is extracted. There might be a split aromatic
    // ring, for example.
    for conf in isomer.conformers() {
        let mut feats = Vec::new();
        find_features(isomer, conf, &mut feats, Some(&frag_atoms))?;
        shape_opts.custom_features.push(feats);
    }
    duplicate_join_isotope(isomer, &atom_isotopes, &mut dummy_radii, &mut frag_atoms)?;
    shape_opts.atom_radii = dummy_radii;
    shape_opts.atom_subset = frag_atoms;
    split_dummy_dummy_bonds(isomer)?;
    match SynthonShapeInput::new(isomer, -1, &shape_opts, &ShapeOverlayOptions::default()) {
        Ok(s) => Ok(Some(s)),
        Err(RDKitError::Value(m)) => {
            // It fails if it doesn't have a radius for an atom in the
            // molecule.
            crate::rd_log_warning!("{m}");
            Ok(None)
        }
        Err(e) => Err(e),
    }
}

/// `makeShapesFromMol` for one sample molecule: the shapes of the synthon it
/// was made for, `None` if they couldn't be made.
pub fn make_shapes_from_sample(
    pool: &[Synthon],
    reaction: &SynthonSet,
    sample: &SampleMolRec,
    dg_params: &EmbedParameters,
    shape_build_params: &ShapeBuildParams,
) -> Result<Option<SynthonShapeInput>> {
    let Ok(mut mol) = reaction.build_molecule(pool, &sample.synthon_nums) else {
        return Ok(None);
    };
    let mut mol = trim_sample_mol(&mut mol, sample.synthon_set_num)?;
    let user_gen = shape_build_params.user_conformer_generator.as_ref();
    let mut isomer_confs = generate_isomer_conformers(
        &mol,
        shape_build_params.num_confs,
        true,
        &shape_build_params.stereo_enum_opts,
        dg_params,
        user_gen,
        MAX_STEREO_CENTRES_TO_DO,
    )?;
    if isomer_confs.is_empty() {
        // Sometimes we can get something from the un-enumerated isomers
        // which will be better than nothing.
        isomer_confs = generate_isomer_conformers(
            &mol,
            shape_build_params.num_confs,
            false,
            &shape_build_params.stereo_enum_opts,
            dg_params,
            user_gen,
            MAX_STEREO_CENTRES_TO_DO,
        )?;
        if isomer_confs.is_empty() {
            crate::rd_log_warning!(
                "No conformers generated for sample molecule {} : {} when generating conformers for synthon {}",
                mol.name(),
                super::details::smiles_of(&mol),
                pool[sample.synthon_idx].smiles()
            );
            return Ok(None);
        }
    }
    let num_isomers = isomer_confs.len();
    let mut all_shapes: Option<SynthonShapeInput> = None;
    for isomer in &mut isomer_confs {
        if isomer.num_conformers() == 0 {
            continue;
        }
        let shapes = make_isomer_shapes(isomer, sample.synthon_set_num as u32, shape_build_params)
            .map_err(|e| {
                RDKitError::Runtime(format!(
                    "{e}\nNo shapes generated for sample molecule {} : {} when generating conformers for synthon {} with isomer : {}",
                    mol.name(),
                    super::details::smiles_of(&mol),
                    pool[sample.synthon_idx].smiles(),
                    super::details::smiles_of(isomer)
                ))
            })?;
        let Some(shapes) = shapes else { continue };
        match &mut all_shapes {
            // Because stereoisomers should all have the same number of atoms
            // and bonds, we can just combine the shapes into one set. We
            // don't need to keep track of which stereoisomer they came from.
            Some(all) => all.merge(shapes),
            None => all_shapes = Some(shapes),
        }
    }
    if let Some(all) = &mut all_shapes
        && num_isomers > 1
    {
        // all_shapes will be the product of at least 1 merge, so prune.
        all.shapes_mut()
            .prune_shapes(shape_build_params.shape_sim_threshold)?;
    }
    mol.clear_conformers();
    Ok(all_shapes)
}

/// A sample molecule for a synthon and the reaction it was made from.
#[derive(Clone)]
struct SampleJob {
    reaction: usize,
    rec: SampleMolRec,
}

fn write_interim_file(space: &SynthonSpace, filename: &str) -> Result<()> {
    // Write to a temporary file and then move it onto filename. That way if
    // there's a crash during the write, we reduce the chance of losing
    // everything.
    let temp_file = format!("{filename}.tmp");
    space.write_db_file(&temp_file)?;
    std::fs::rename(&temp_file, filename)
        .map_err(|e| RDKitError::Runtime(format!("Couldn't rename {temp_file}: {e}")))
}

impl SynthonSpace {
    /// `buildSynthonSampleMolecules`: for each synthon that needs shapes, the
    /// sample molecules of the reactions it is in.
    fn build_synthon_sample_molecules(&self, max_synthon_atoms: u32) -> Vec<Vec<SampleJob>> {
        let mut all = Vec::with_capacity(self.synthon_reactions.len());
        for (idx, reactions) in self.synthon_reactions.iter().enumerate() {
            let synthon = &self.pool[idx];
            if reactions.is_empty()
                || (max_synthon_atoms != 0 && synthon.num_heavy_atoms() > max_synthon_atoms)
                || synthon.has_shapes()
            {
                continue;
            }
            let mut these: Vec<SampleJob> = reactions
                .iter()
                .filter_map(|&r| {
                    // In the unlikely event that we didn't get anything, drop
                    // it.
                    self.reactions[r]
                        .make_sample_molecule(&self.pool, idx)
                        .filter(|rec| rec.num_atoms > 0)
                        .map(|rec| SampleJob { reaction: r, rec })
                })
                .collect();
            // Sort in descending order of SynthonSet size, so that if a
            // synthon is used more than once, the sample is most
            // representative of the type of chemistry it's involved with.
            libcxx_sort_by(&mut these, |a, b| {
                self.reactions[a.reaction].num_products()
                    > self.reactions[b.reaction].num_products()
            });
            all.push(these);
        }
        all
    }

    /// `buildSynthonShapes`: makes the shapes of the synthons that don't have
    /// any yet, writing the interim file as it goes if asked.
    pub fn build_synthon_shapes(&mut self, shape_build_params: &ShapeBuildParams) -> Result<()> {
        if self.num_conformers == shape_build_params.num_confs {
            if self.pool.iter().any(|s| !s.has_shapes()) {
                crate::rd_log_warning!("Resuming building SynthonSpace shapes.");
            } else {
                crate::rd_log_warning!(
                    "SynthonSpace has already been built with {} conformers.",
                    shape_build_params.num_confs
                );
                return Ok(());
            }
        }
        crate::rd_log_warning!(
            "Building the conformers of {} synthons may take some time.",
            self.pool.len()
        );
        if self.synthon_reactions.is_empty() {
            self.fill_synthon_reactions();
        }
        self.num_conformers = shape_build_params.num_confs;
        let mut all_samples =
            self.build_synthon_sample_molecules(shape_build_params.max_synthon_atoms);
        // Do the big molecules first so there is less chance of sitting at the
        // end waiting for 1 big embedding to finish.
        libcxx_sort_by(&mut all_samples, |a, b| {
            a[0].rec.num_atoms > b[0].rec.num_atoms
        });
        let interim_write =
            shape_build_params.interim_writes != 0 && !shape_build_params.interim_file.is_empty();
        loop {
            // Loop around until all synthons have some shapes or we've run
            // out of synthon/reaction combinations.
            let mut samples: Vec<SampleJob> = Vec::new();
            for all in &mut all_samples {
                // If we have a shape object in the synthon with some shapes,
                // don't do anything. If there's no shapes then probably the
                // embedding failed last time, so try again. It might be the
                // fault of this synthon, but it might conceivably have been
                // a problem with whatever it was attached to, and another
                // reaction might have bolted on something more amenable.
                if let Some(last) = all.last()
                    && !self.pool[last.rec.synthon_idx].has_shapes()
                {
                    samples.extend(all.pop());
                }
                if interim_write && shape_build_params.interim_writes == samples.len() as u64 {
                    break;
                }
            }
            if samples.is_empty() {
                break;
            }
            libcxx_sort_by(&mut samples, |a, b| a.rec.num_atoms > b.rec.num_atoms);
            let mut dg_params = EmbedParameters::etkdg_v3();
            dg_params.num_threads = 1;
            dg_params.prune_rms_thresh = shape_build_params.rms_threshold;
            dg_params.random_seed = shape_build_params.random_seed;
            dg_params.max_iterations = shape_build_params.max_embed_attempts;
            dg_params.timeout = shape_build_params.time_out;
            let with_shapes_before = self.num_synthons_with_shapes();
            let results = par_map(&samples, shape_build_params.num_threads, |job| {
                make_shapes_from_sample(
                    &self.pool,
                    &self.reactions[job.reaction],
                    &job.rec,
                    &dg_params,
                    shape_build_params,
                )
            });
            for (job, shapes) in samples.iter().zip(results) {
                if let Some(shapes) = shapes? {
                    self.pool[job.rec.synthon_idx].set_shapes(shapes);
                }
            }
            let with_shapes = self.num_synthons_with_shapes();
            if interim_write && with_shapes > with_shapes_before {
                crate::rd_log_info!(
                    "Writing interim file {} with {} shapes.",
                    shape_build_params.interim_file,
                    with_shapes
                );
                write_interim_file(self, &shape_build_params.interim_file)?;
            }
        }
        Ok(())
    }
}
