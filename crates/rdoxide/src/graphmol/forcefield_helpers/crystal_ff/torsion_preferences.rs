//! Collecting the torsion, inversion and planarity knowledge terms of a
//! molecule (upstream `TorsionPreferences.cpp`).

use std::f64::consts::PI;

use super::params::{
    ExpTorsionAngle, GaussianExpTorsionAngle, cosine_torsion_params, gaussian_torsion_params,
};
use crate::error::{RDKitError, Result};
use crate::forcefield::crystal_ff::{
    LOOKUP_GRID_SIZE, gaussian_torsion_de_dphi, gaussian_torsion_energy,
};
use crate::graphmol::atom::HybridizationType;
use crate::graphmol::dist_geom_helpers::{InternalCoordinates, Path14Configuration, TorsionType};
use crate::graphmol::romol::ROMol;
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rd_log_info;

/// the "macrocycle" patterns for ETKDGv3 use a minimum ring size of 9
const MIN_MACROCYCLE_SIZE: usize = 9;

// parameters for ETKDGv4
const AROMATIC_HEIGHTS: [f64; 1] = [2.305];
const AROMATIC_POSITIONS: [f64; 1] = [0.0];
const TRANS_POSITION: [f64; 1] = [PI];
const CIS_POSITION: [f64; 1] = [0.0];
const AROMATIC_WIDTHS: [f64; 1] = [0.251];

/// the torsion indices of the forced trans/cis amides and of aromatic rings
pub const AMIDE_TRANS_IDX: usize = 2000;
pub const AMIDE_CIS_IDX: usize = 2001;
pub const AROM_TORSION_IDX: usize = 2002;

/// Force constants of the ETKDG terms (`ETKDGForceConsts::Params`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EtkdgForceConsts {
    pub distance: f64,
    pub fourth_dim: f64,
    pub chiral: f64,
    pub k_term_angle: f64,
    pub k_term_improper: f64,
    pub k_term_torsion: f64,
    pub et_term_scaling: f64,
}

impl Default for EtkdgForceConsts {
    fn default() -> Self {
        Self::SEQ_COSINE
    }
}

impl EtkdgForceConsts {
    /// sequential embedding with cosine torsions
    pub const SEQ_COSINE: Self = Self {
        distance: 1.0,
        fourth_dim: 1.0,
        chiral: 1.0,
        k_term_angle: 1.0,
        k_term_improper: 10.0,
        k_term_torsion: 100.0,
        et_term_scaling: 1.0,
    };
    /// sequential embedding with Gaussian torsions
    pub const SEQ_GAUSSIAN: Self = Self::SEQ_COSINE;
    /// all-in-one embedding with cosine torsions
    pub const AIO_COSINE: Self = Self {
        distance: 2.15,
        fourth_dim: 2.15,
        chiral: 1.0,
        k_term_angle: 0.1,
        k_term_improper: 0.001,
        k_term_torsion: 2.15,
        et_term_scaling: 0.05,
    };
    /// all-in-one embedding with Gaussian torsions
    pub const AIO_GAUSSIAN: Self = Self {
        k_term_torsion: 0.75,
        ..Self::AIO_COSINE
    };
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TorsionParamKind {
    #[default]
    Cosine,
    Gaussian,
}

/// The parameters of one torsion term.
#[derive(Debug, Clone, PartialEq)]
pub enum TorsionParams {
    /// signs and force constants of the six cosines
    Cosine { signs: Vec<i32>, v: Vec<f64> },
    /// heights, positions and widths of the Gaussians and a scaling factor
    Gaussian {
        heights: Vec<f64>,
        positions: Vec<f64>,
        widths: Vec<f64>,
        scaling: f64,
    },
}

/// The atoms of an inversion term: the neighbours of the central sp2 atom
/// (slots 0, 2 and 3 of the term, the centre being slot 1), its atomic number
/// and whether a carbon centre is bound to an sp2 oxygen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ImproperAtoms {
    pub atoms: [u32; 4],
    pub center_atomic_num: u32,
    pub bound_to_sp2_o: bool,
}

/// The pattern that matched a torsion.
#[derive(Debug, Clone, Copy)]
pub enum TorsionAngleRef {
    Cosine(&'static ExpTorsionAngle),
    Gaussian(&'static GaussianExpTorsionAngle),
}

/// A matched torsion: (central bond index, the four atoms, the pattern).
pub type TorsionBond = (u32, [u32; 4], TorsionAngleRef);

/// Everything the ETKDG force field needs to know about a molecule.
#[derive(Debug, Default)]
pub struct CrystalFFDetails {
    pub torsion_param_kind: TorsionParamKind,
    pub exp_torsion_atoms: Vec<[u32; 4]>,
    pub exp_torsion_angles: Vec<TorsionParams>,
    pub torsion_idx: Vec<usize>,
    pub improper_atoms: Vec<ImproperAtoms>,
    pub bonds: Vec<(i32, i32)>,
    pub angles: Vec<[i32; 4]>,
    pub atom_nums: Vec<i32>,
    pub bounds_mat_force_scaling: f64,
    /// empty if no atoms are constrained
    pub constrained_atoms: Vec<bool>,
    /// the topological distance matrix (row-major)
    pub dist_mat: Vec<f64>,
    pub force_consts: EtkdgForceConsts,
    pub path14_configs: Vec<Path14Configuration>,
    pub phi_to_energy: Vec<Vec<f64>>,
    pub phi_to_grad: Vec<Vec<f64>>,
    pub internal_coords: Option<Box<InternalCoordinates>>,
}

/// Options of [`get_experimental_torsions`].
#[derive(Debug, Clone, Copy)]
pub struct TorsionOptions {
    pub use_exp_torsions: bool,
    pub use_small_ring_torsions: bool,
    pub use_macrocycle_torsions: bool,
    pub use_basic_knowledge: bool,
    pub version: u32,
    pub verbose: bool,
}

impl Default for TorsionOptions {
    fn default() -> Self {
        TorsionOptions {
            use_exp_torsions: false,
            use_small_ring_torsions: false,
            use_macrocycle_torsions: false,
            use_basic_knowledge: false,
            version: 2,
            verbose: false,
        }
    }
}

/// What differs between cosine and Gaussian torsion parameters.
trait AngleKind: Sized + 'static {
    fn collection(options: &TorsionOptions) -> Result<&'static [Self]>;
    fn torsion_idx(&self) -> usize;
    fn pattern(&self) -> &ROMol;
    fn atom_indices(&self) -> [u32; 4];
    fn as_ref(&'static self) -> TorsionAngleRef;
    /// the parameters of the experimental term, scaled by `et_term_scaling`
    fn term_params(&self, et_term_scaling: f64) -> TorsionParams;
    fn log_match(&self, atoms: [u32; 4]) -> String;
    /// the parameters of the forced cis or trans amide term
    fn amide_params(is_cis: bool, is_aio: bool) -> TorsionParams;
    /// the parameters of the term of a flat ring torsion
    fn flat_ring_params(k_term_torsion: f64) -> TorsionParams;
}

impl AngleKind for ExpTorsionAngle {
    fn collection(options: &TorsionOptions) -> Result<&'static [Self]> {
        cosine_torsion_params(
            options.version,
            options.use_small_ring_torsions,
            options.use_macrocycle_torsions,
        )
    }
    fn torsion_idx(&self) -> usize {
        self.torsion_idx
    }
    fn pattern(&self) -> &ROMol {
        &self.pattern
    }
    fn atom_indices(&self) -> [u32; 4] {
        self.idx
    }
    fn as_ref(&'static self) -> TorsionAngleRef {
        TorsionAngleRef::Cosine(self)
    }
    fn term_params(&self, et_term_scaling: f64) -> TorsionParams {
        TorsionParams::Cosine {
            signs: self.signs.clone(),
            v: self.v.iter().map(|v| v * et_term_scaling).collect(),
        }
    }
    fn log_match(&self, a: [u32; 4]) -> String {
        let terms: Vec<String> = self
            .signs
            .iter()
            .zip(&self.v)
            .map(|(s, v)| format!("({s} {v})"))
            .collect();
        format!(
            "{}: {} {} {} {}, [{}] ",
            self.smarts,
            a[0],
            a[1],
            a[2],
            a[3],
            terms.join(", ")
        )
    }
    fn amide_params(is_cis: bool, is_aio: bool) -> TorsionParams {
        let mut v = vec![0.0; 6];
        let mut signs = vec![1; 6];
        v[0] = if is_aio { 4.0 } else { 75.0 };
        if is_cis {
            signs[0] = -1;
        }
        TorsionParams::Cosine { signs, v }
    }
    fn flat_ring_params(k_term_torsion: f64) -> TorsionParams {
        let mut signs = vec![1; 6];
        signs[1] = -1; // MMFF sign for m = 2
        let mut v = vec![0.0; 6];
        v[1] = k_term_torsion; // 7.0 is the MMFF force constant for aromatic rings
        TorsionParams::Cosine { signs, v }
    }
}

impl AngleKind for GaussianExpTorsionAngle {
    fn collection(options: &TorsionOptions) -> Result<&'static [Self]> {
        gaussian_torsion_params(
            options.version,
            options.use_small_ring_torsions,
            options.use_macrocycle_torsions,
        )
    }
    fn torsion_idx(&self) -> usize {
        self.torsion_idx
    }
    fn pattern(&self) -> &ROMol {
        &self.pattern
    }
    fn atom_indices(&self) -> [u32; 4] {
        self.idx
    }
    fn as_ref(&'static self) -> TorsionAngleRef {
        TorsionAngleRef::Gaussian(self)
    }
    fn term_params(&self, et_term_scaling: f64) -> TorsionParams {
        TorsionParams::Gaussian {
            heights: self.heights.clone(),
            positions: self.positions.clone(),
            widths: self.widths.clone(),
            scaling: et_term_scaling,
        }
    }
    fn log_match(&self, a: [u32; 4]) -> String {
        format!("{}: {} {} {} {}, [", self.smarts, a[0], a[1], a[2], a[3])
    }
    fn amide_params(is_cis: bool, is_aio: bool) -> TorsionParams {
        TorsionParams::Gaussian {
            heights: AROMATIC_HEIGHTS.to_vec(),
            positions: if is_cis { CIS_POSITION } else { TRANS_POSITION }.to_vec(),
            widths: AROMATIC_WIDTHS.to_vec(),
            scaling: if is_aio { 0.05 } else { 1.0 },
        }
    }
    fn flat_ring_params(k_term_torsion: f64) -> TorsionParams {
        TorsionParams::Gaussian {
            heights: AROMATIC_HEIGHTS.to_vec(),
            positions: AROMATIC_POSITIONS.to_vec(),
            widths: AROMATIC_WIDTHS.to_vec(),
            scaling: k_term_torsion,
        }
    }
}

/// The bonds of rings that share more than one bond with another ring (bridged
/// systems) are excluded, unless both rings are macrocycles.
fn excluded_bonds(mol: &ROMol) -> Vec<bool> {
    let nb = mol.num_bonds() as usize;
    let mut excluded = vec![false; nb];
    let bond_rings = mol.ring_info().bond_rings();
    for (i, ring_i) in bond_rings.iter().enumerate() {
        let mut in_ring_i = vec![false; nb];
        for &b in ring_i {
            in_ring_i[b as usize] = true;
        }
        for ring_j in &bond_rings[i + 1..] {
            // we don't worry about the overlap if both rings are macrocycles
            if ring_i.len() >= MIN_MACROCYCLE_SIZE && ring_j.len() >= MIN_MACROCYCLE_SIZE {
                continue;
            }
            let in_common = ring_j
                .iter()
                .filter(|&&b| in_ring_i[b as usize])
                .take(2)
                .count();
            if in_common > 1 {
                // more than one bond in common: exclude the bonds of the
                // rings that are not macrocycles
                for ring in [ring_i, ring_j] {
                    if ring.len() < MIN_MACROCYCLE_SIZE {
                        for &b in ring {
                            excluded[b as usize] = true;
                        }
                    }
                }
            }
        }
    }
    excluded
}

fn all_constrained(details: &CrystalFFDetails, atoms: [u32; 4]) -> bool {
    !details.constrained_atoms.is_empty()
        && atoms.iter().all(|&a| details.constrained_atoms[a as usize])
}

fn get_experimental_torsions_impl<A: AngleKind>(
    mol: &ROMol,
    details: &mut CrystalFFDetails,
    torsion_bonds: &mut Vec<TorsionBond>,
    options: &TorsionOptions,
) -> Result<()> {
    torsion_bonds.clear();
    let na = mol.num_atoms();
    if na == 0 {
        return Err(RDKitError::Value("molecule has no atoms".into()));
    }
    details.exp_torsion_atoms.clear();
    details.exp_torsion_angles.clear();
    details.improper_atoms.clear();
    details.torsion_idx.clear();

    let excluded = excluded_bonds(mol);
    let rinfo = mol.ring_info();
    let mut done_bonds = vec![false; mol.num_bonds() as usize];

    // apply basic knowledge such as flat aromatic rings, other sp2-centers,
    // straight triple bonds, etc.
    if options.use_basic_knowledge {
        // torsions for forced trans amides / esters
        let forced: Vec<Path14Configuration> = details
            .path14_configs
            .iter()
            .filter(|c| {
                c.value.is_forced && matches!(c.value.kind, TorsionType::Trans | TorsionType::Cis)
            })
            .copied()
            .collect();
        for config in forced {
            let atoms = [config.aid1, config.aid2, config.aid3, config.aid4];
            let bnd_idx = mol
                .bond_idx_between(atoms[1], atoms[2])
                .ok_or_else(|| RDKitError::invariant("no bond between the central atoms"))?
                as usize;
            if excluded[bnd_idx] || rinfo.num_bond_rings(bnd_idx as u32)? > 3 {
                done_bonds[bnd_idx] = true;
            }
            if done_bonds[bnd_idx] || all_constrained(details, atoms) {
                continue;
            }
            let is_cis = config.value.kind == TorsionType::Cis;
            details.exp_torsion_atoms.push(atoms);
            let is_aio = (details.force_consts.et_term_scaling - 1.0).abs() > 1e-3;
            details.torsion_idx.push(if is_cis {
                AMIDE_CIS_IDX
            } else {
                AMIDE_TRANS_IDX
            });
            details
                .exp_torsion_angles
                .push(A::amide_params(is_cis, is_aio));
        }
    }

    if options.use_exp_torsions {
        // the torsion angles with experimental data
        for param in A::collection(options)? {
            let matches = substruct_match(
                mol,
                param.pattern(),
                &SubstructMatchParameters {
                    uniquify: false,
                    recursion_possible: true,
                    ..Default::default()
                },
            )?;
            let idx = param.atom_indices();
            for m in matches {
                let atoms = idx.map(|i| m[i as usize].1 as u32);
                let bid2 = mol
                    .bond_idx_between(atoms[1], atoms[2])
                    .ok_or_else(|| RDKitError::invariant("bond between central atoms not found"))?;
                // a bond may be part of at most one ring
                if excluded[bid2 as usize] || rinfo.num_bond_rings(bid2)? > 3 {
                    done_bonds[bid2 as usize] = true;
                }
                if done_bonds[bid2 as usize] {
                    continue;
                }
                // no ET terms between constrained atoms
                if all_constrained(details, atoms) {
                    continue;
                }
                torsion_bonds.push((bid2, atoms, param.as_ref()));
                done_bonds[bid2 as usize] = true;
                details.exp_torsion_atoms.push(atoms);
                details.torsion_idx.push(param.torsion_idx());
                details
                    .exp_torsion_angles
                    .push(param.term_params(details.force_consts.et_term_scaling));
                if options.verbose {
                    rd_log_info!("{}", param.log_match(atoms));
                }
            }
        }
    }

    if options.use_basic_knowledge {
        // inversion terms (improper torsions / out-of-plane bends / inversion)
        for aid2 in 0..na {
            let atom2 = mol.atom(aid2);
            let z = atom2.atomic_num();
            // a N, O or C that is sp2 hybridized and has three neighbors
            if (z == 6 || z == 7 || z == 8)
                && atom2.hybridization() == HybridizationType::SP2
                && mol.atom_degree(aid2) == 3
            {
                let mut atoms = [0u32; 4];
                atoms[1] = aid2;
                let mut bound_to_sp2_o = false;
                // the neighbours go to the slots 0, 2 and 3
                for (slot, nbr) in [0usize, 2, 3].into_iter().zip(mol.atom_neighbors(aid2)) {
                    atoms[slot] = nbr;
                    // an sp2 carbon bound to an sp2 oxygen
                    if !bound_to_sp2_o {
                        bound_to_sp2_o = z == 6
                            && mol.atom(nbr).atomic_num() == 8
                            && mol.atom(nbr).hybridization() == HybridizationType::SP2;
                    }
                }
                details.improper_atoms.push(ImproperAtoms {
                    atoms,
                    center_atomic_num: z,
                    bound_to_sp2_o,
                });
            }
        }

        // torsions for flat rings
        if !rinfo.is_initialized() {
            return Err(RDKitError::invariant("ring info not initialized"));
        }
        for ring in rinfo.atom_rings() {
            let r_size = ring.len();
            // 3-membered rings need no treatment, rings above 6 are not treated
            if !(4..=6).contains(&r_size) {
                continue;
            }
            for i in 0..r_size {
                let atoms = [
                    ring[i] as u32,
                    ring[(i + 1) % r_size] as u32,
                    ring[(i + 2) % r_size] as u32,
                    ring[(i + 3) % r_size] as u32,
                ];
                let bid2 = mol
                    .bond_idx_between(atoms[1], atoms[2])
                    .ok_or_else(|| RDKitError::invariant("no bond between ring atoms"))?
                    as usize;
                // if all 4 atoms are sp2, add the torsion
                if !done_bonds[bid2]
                    && atoms
                        .iter()
                        .all(|&a| mol.atom(a).hybridization() == HybridizationType::SP2)
                {
                    done_bonds[bid2] = true;
                    details.exp_torsion_atoms.push(atoms);
                    details.torsion_idx.push(AROM_TORSION_IDX);
                    details
                        .exp_torsion_angles
                        .push(A::flat_ring_params(details.force_consts.k_term_torsion));
                }
            }
        }
    }
    Ok(())
}

/// Collects the torsion terms of the molecule: forced amide torsions, the
/// experimental preferences, inversion centres and flat-ring torsions, and
/// the matched torsions with their patterns.
pub fn get_experimental_torsions_with_bonds(
    mol: &ROMol,
    details: &mut CrystalFFDetails,
    options: &TorsionOptions,
) -> Result<Vec<TorsionBond>> {
    let mut torsion_bonds = Vec::new();
    details.torsion_param_kind = if options.version == 4 {
        TorsionParamKind::Gaussian
    } else {
        TorsionParamKind::Cosine
    };
    match details.torsion_param_kind {
        TorsionParamKind::Gaussian => get_experimental_torsions_impl::<GaussianExpTorsionAngle>(
            mol,
            details,
            &mut torsion_bonds,
            options,
        )?,
        TorsionParamKind::Cosine => {
            get_experimental_torsions_impl::<ExpTorsionAngle>(
                mol,
                details,
                &mut torsion_bonds,
                options,
            )?;
        }
    }
    Ok(torsion_bonds)
}

/// [`get_experimental_torsions_with_bonds`] without the matched torsions.
pub fn get_experimental_torsions(
    mol: &ROMol,
    details: &mut CrystalFFDetails,
    options: &TorsionOptions,
) -> Result<()> {
    get_experimental_torsions_with_bonds(mol, details, options).map(|_| ())
}

/// Computes the energy and gradient lookup tables of the Gaussian torsion
/// terms, once for each distinct torsion index, and replaces the torsion
/// indices by the rank of the index among the distinct ones.
///
/// The order of the terms is kept, which matters for the handling of trans
/// amides.
pub fn populate_ref_table(details: &mut CrystalFFDetails) -> Result<()> {
    let mut sorted = details.torsion_idx.clone();
    sorted.sort_unstable();
    sorted.dedup();
    for idx in &mut details.torsion_idx {
        *idx = sorted.partition_point(|&s| s < *idx);
    }
    details.phi_to_energy = vec![vec![0.0; LOOKUP_GRID_SIZE]; sorted.len()];
    details.phi_to_grad = vec![vec![0.0; LOOKUP_GRID_SIZE]; sorted.len()];
    for torsion_idx in 0..sorted.len() {
        // the first occurrence of the index locates its parameters
        let term_idx = details
            .torsion_idx
            .iter()
            .position(|&t| t == torsion_idx)
            .expect("every rank occurs");
        let TorsionParams::Gaussian {
            heights,
            positions,
            widths,
            ..
        } = &details.exp_torsion_angles[term_idx]
        else {
            return Err(RDKitError::invariant(
                "lookup tables need Gaussian torsion parameters",
            ));
        };
        for grid_point in 0..LOOKUP_GRID_SIZE {
            let phi = grid_point as f64 * PI / (LOOKUP_GRID_SIZE - 1) as f64;
            details.phi_to_energy[torsion_idx][grid_point] =
                gaussian_torsion_energy(heights, positions, widths, phi);
            details.phi_to_grad[torsion_idx][grid_point] =
                gaussian_torsion_de_dphi(heights, positions, widths, phi);
        }
    }
    Ok(())
}
