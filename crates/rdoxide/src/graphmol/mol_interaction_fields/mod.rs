//! Molecular interaction field (MIF) descriptors: Coulomb, distance
//! dependent dielectric Coulomb, van der Waals (MMFF and UFF probes),
//! evaluated on a uniform grid around a molecule.
//!
//! Ported from RDKit `Code/GraphMol/MolInteractionFields` (2026.09.1).
//! Note that, as upstream says, this functionality is experimental.
//!
//! Copyright (c) 2014-2024, Novartis Institutes for BioMedical Research and
//! other RDKit contributors. Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-BSD-3-Clause.txt).

mod cube;

use crate::geometry::Point3D;
use crate::geometry::grid::UniformRealValueGrid3D;
use crate::graphmol::forcefield_helpers::mmff::MmffMolProperties;
use crate::graphmol::forcefield_helpers::uff::get_atom_types;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDValue;
use crate::rdgeneral::threads::par_map;
use crate::{RDKitError, Result};

pub use cube::{
    read_from_cube_file, read_from_cube_stream, write_to_cube_file, write_to_cube_stream,
};

pub(crate) const CUTOFF: f64 = 0.001;
pub(crate) const MIN_CUTOFF_VAL: f64 = CUTOFF * CUTOFF;

/// A function of the position of a probe (`x`, `y`, `z`) and the squared
/// distance up to which atoms contribute.
pub trait Descriptor: Sync {
    fn value(&self, x: f64, y: f64, z: f64, thres: f64) -> f64;
}

/// `constructGrid`: a grid around the conformer with a margin on each side.
pub fn construct_grid(
    mol: &ROMol,
    conf_id: i32,
    margin: f64,
    spacing: f64,
) -> Result<UniformRealValueGrid3D> {
    if mol.num_conformers() == 0 {
        return Err(RDKitError::precondition(
            "No conformers available for molecule.",
        ));
    }
    let pts = mol.conformer(conf_id)?.positions();
    let first = pts.first().copied().unwrap_or_default();
    let (mut lo, mut hi) = (first, first);
    for p in pts {
        lo = Point3D::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
        hi = Point3D::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
    }
    lo = Point3D::new(lo.x - margin, lo.y - margin, lo.z - margin);
    hi = Point3D::new(hi.x + margin, hi.y + margin, hi.z + margin);
    UniformRealValueGrid3D::new(hi.x - lo.x, hi.y - lo.y, hi.z - lo.z, spacing, Some(lo))
}

/// `calculateDescriptors`: fills the grid with the values of the descriptor;
/// a negative threshold means "no threshold". The planes of the grid are
/// computed on `num_threads` threads (the values do not depend on it).
pub fn calculate_descriptors<D: Descriptor>(
    grd: &mut UniformRealValueGrid3D,
    descriptor: &D,
    thres: f64,
    num_threads: i32,
) -> Result<()> {
    let off = grd.offset();
    let spacing = grd.spacing();
    let thres = if thres < 0.0 {
        spacing * f64::from(grd.size())
    } else {
        thres
    };
    let thres = thres * thres; // comparing squared distances
    // the coordinates are accumulated like upstream does
    let axis = |start: f64, n: u32| -> Vec<f64> {
        let mut v = Vec::with_capacity(n as usize);
        let mut x = start;
        for _ in 0..n {
            v.push(x);
            x += spacing;
        }
        v
    };
    let xs = axis(off.x, grd.num_x());
    let ys = axis(off.y, grd.num_y());
    let zs = axis(off.z, grd.num_z());
    let planes = par_map(&zs, num_threads, |&z| {
        let mut plane = Vec::with_capacity(xs.len() * ys.len());
        for &y in &ys {
            for &x in &xs {
                plane.push(descriptor.value(x, y, z, thres));
            }
        }
        plane
    });
    let mut id = 0;
    for plane in planes {
        for v in plane {
            grd.set_val(id, v)?;
            id += 1;
        }
    }
    Ok(())
}

fn positions_and_charges(mol: &ROMol, conf_id: i32, prop: &str) -> Result<(Vec<f64>, Vec<f64>)> {
    let conf = mol.conformer(conf_id)?;
    let mut charges = Vec::new();
    let mut pos = Vec::new();
    for i in 0..mol.num_atoms() {
        let c = match mol.atom(i).props.get(prop) {
            Some(RDValue::Double(d)) => *d,
            Some(v) => v
                .as_double()
                .ok_or_else(|| RDKitError::Key(prop.to_string()))?,
            None => return Err(RDKitError::Key(prop.to_string())),
        };
        charges.push(c);
        let p = conf.atom_pos(i)?;
        pos.extend([p.x, p.y, p.z]);
    }
    Ok((charges, pos))
}

fn flatten(positions: &[Point3D]) -> Vec<f64> {
    positions.iter().flat_map(|p| [p.x, p.y, p.z]).collect()
}

#[allow(clippy::approx_constant)] // (upstream's value of pi)
const PREFACTOR: f64 = 1.0 / (4.0 * 3.141592 * 8.854188) * 1.602 * 1.602 * 6.02214129 * 10000.0;

/// Coulomb interaction between the charges of a molecule and a probe in
/// vacuum.
#[derive(Debug, Clone)]
pub struct Coulomb {
    n_atoms: usize,
    softcore: bool,
    abs_val: bool,
    cutoff: f64,
    probe: f64,
    alpha: f64,
    charges: Vec<f64>,
    pos: Vec<f64>,
}

impl Coulomb {
    fn finish(mut self) -> Self {
        if self.alpha.abs() < MIN_CUTOFF_VAL {
            self.softcore = false;
            if self.cutoff < MIN_CUTOFF_VAL {
                self.cutoff = CUTOFF;
            }
        } else {
            self.softcore = true;
        }
        self
    }

    /// From charges and positions.
    pub fn new(
        charges: &[f64],
        positions: &[Point3D],
        probe_charge: f64,
        abs_val: bool,
        alpha: f64,
        cutoff: f64,
    ) -> Result<Self> {
        if charges.len() != positions.len() {
            return Err(RDKitError::precondition(
                "Lengths of positions and charges vectors do not match.",
            ));
        }
        Ok(Coulomb {
            n_atoms: charges.len(),
            softcore: false,
            abs_val,
            cutoff: cutoff * cutoff,
            probe: probe_charge,
            alpha,
            charges: charges.to_vec(),
            pos: flatten(positions),
        }
        .finish())
    }

    /// From a molecule whose atoms have the property `prop` (a partial
    /// charge).
    pub fn from_mol(
        mol: &ROMol,
        conf_id: i32,
        probe_charge: f64,
        abs_val: bool,
        prop: &str,
        alpha: f64,
        cutoff: f64,
    ) -> Result<Self> {
        let (charges, pos) = positions_and_charges(mol, conf_id, prop)?;
        Ok(Coulomb {
            n_atoms: charges.len(),
            softcore: false,
            abs_val,
            cutoff: cutoff * cutoff,
            probe: probe_charge,
            alpha,
            charges,
            pos,
        }
        .finish())
    }
}

impl Descriptor for Coulomb {
    fn value(&self, x: f64, y: f64, z: f64, thres: f64) -> f64 {
        let mut res = 0.0;
        for i in 0..self.n_atoms {
            let j = 3 * i;
            let (dx, dy, dz) = (x - self.pos[j], y - self.pos[j + 1], z - self.pos[j + 2]);
            let dist2 = dx * dx + dy * dy + dz * dz;
            if dist2 < thres {
                if self.softcore {
                    res += self.charges[i] * (1.0 / (self.alpha + dist2).sqrt());
                } else {
                    res += self.charges[i] * (1.0 / dist2.max(self.cutoff).sqrt());
                }
            }
        }
        res *= PREFACTOR * self.probe;
        if self.abs_val {
            res = -res.abs();
        }
        res
    }
}

/// The screening factor for a number of neighbours within 4 A.
fn screening(neighbors: u32) -> f64 {
    match neighbors {
        0..=6 => 0.0,
        7 => 0.4,
        8 => 0.9,
        9 => 1.4,
        10 => 1.9,
        11 => 2.6,
        _ => 4.0,
    }
}

/// Coulomb interaction with a distance dependent dielectric (the energy term
/// of GRID; J. Med. Chem. 1985, 28, 849 and J. Comp. Chem. 1983, 4, 187).
#[derive(Debug, Clone)]
pub struct CoulombDielectric {
    n_atoms: usize,
    softcore: bool,
    abs_val: bool,
    cutoff: f64,
    probe: f64,
    dielectric: f64,
    xi: f64,
    alpha: f64,
    charges: Vec<f64>,
    pos: Vec<f64>,
    sp: Vec<f64>,
}

impl CoulombDielectric {
    #[allow(clippy::too_many_arguments)]
    fn build(
        charges: Vec<f64>,
        pos: Vec<f64>,
        probe_charge: f64,
        abs_val: bool,
        alpha: f64,
        cutoff: f64,
        epsilon: f64,
        xi: f64,
    ) -> Self {
        let n = charges.len();
        let mut neighbors = vec![0u32; n];
        for i in 0..n.saturating_sub(1) {
            for j in i + 1..n {
                let (dx, dy, dz) = (
                    pos[3 * j] - pos[3 * i],
                    pos[3 * j + 1] - pos[3 * i + 1],
                    pos[3 * j + 2] - pos[3 * i + 2],
                );
                if dx * dx + dy * dy + dz * dz < 16.0 {
                    neighbors[i] += 1;
                    neighbors[j] += 1;
                }
            }
        }
        let mut res = CoulombDielectric {
            n_atoms: n,
            softcore: false,
            abs_val,
            cutoff: cutoff * cutoff,
            probe: probe_charge,
            dielectric: (xi - epsilon) / (xi + epsilon),
            xi,
            alpha,
            charges,
            pos,
            sp: neighbors.into_iter().map(screening).collect(),
        };
        if alpha.abs() < MIN_CUTOFF_VAL {
            if res.cutoff < MIN_CUTOFF_VAL * MIN_CUTOFF_VAL {
                res.cutoff = CUTOFF * CUTOFF;
            }
        } else {
            res.softcore = true;
        }
        res
    }

    /// From charges and positions.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        charges: &[f64],
        positions: &[Point3D],
        probe_charge: f64,
        abs_val: bool,
        alpha: f64,
        cutoff: f64,
        epsilon: f64,
        xi: f64,
    ) -> Result<Self> {
        if charges.len() != positions.len() {
            return Err(RDKitError::precondition(
                "Lengths of position and charge vectors do not match.",
            ));
        }
        Ok(Self::build(
            charges.to_vec(),
            flatten(positions),
            probe_charge,
            abs_val,
            alpha,
            cutoff,
            epsilon,
            xi,
        ))
    }

    /// From a molecule whose atoms have the property `prop`.
    #[allow(clippy::too_many_arguments)]
    pub fn from_mol(
        mol: &ROMol,
        conf_id: i32,
        probe_charge: f64,
        abs_val: bool,
        prop: &str,
        alpha: f64,
        cutoff: f64,
        epsilon: f64,
        xi: f64,
    ) -> Result<Self> {
        if mol.num_conformers() == 0 {
            return Err(RDKitError::precondition("No Conformers for Molecule"));
        }
        let (charges, pos) = positions_and_charges(mol, conf_id, prop)?;
        Ok(Self::build(
            charges,
            pos,
            probe_charge,
            abs_val,
            alpha,
            cutoff,
            epsilon,
            xi,
        ))
    }
}

impl Descriptor for CoulombDielectric {
    fn value(&self, x: f64, y: f64, z: f64, thres: f64) -> f64 {
        let mut dists = Vec::with_capacity(self.n_atoms);
        let mut neigh = 0;
        for i in 0..self.n_atoms {
            let j = 3 * i;
            let (dx, dy, dz) = (x - self.pos[j], y - self.pos[j + 1], z - self.pos[j + 2]);
            let dist2 = dx * dx + dy * dy + dz * dz;
            dists.push(dist2);
            if dist2 < 16.0 {
                neigh += 1;
            }
        }
        let sq = screening(neigh);
        let mut res = 0.0;
        for i in 0..self.n_atoms {
            if dists[i] < thres {
                if self.softcore {
                    res += self.charges[i]
                        * (1.0 / (self.alpha + dists[i]).sqrt()
                            + self.dielectric
                                / (self.alpha + dists[i] + 4.0 * self.sp[i] * sq).sqrt());
                } else {
                    let dist2 = dists[i].max(self.cutoff);
                    res += self.charges[i]
                        * (1.0 / dist2.sqrt()
                            + self.dielectric / (dist2 + 4.0 * self.sp[i] * sq).sqrt());
                }
            }
        }
        res *= PREFACTOR * (1.0 / self.xi) * self.probe;
        if self.abs_val {
            res = -res.abs();
        }
        res
    }
}

/// Which force field a [`VdWaals`] descriptor uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VdwKind {
    Mmff,
    Uff,
}

/// van der Waals interaction between a molecule and a probe atom, with
/// MMFF or UFF parameters.
#[derive(Debug, Clone)]
pub struct VdWaals {
    kind: VdwKind,
    n_atoms: usize,
    cutoff: f64,
    pos: Vec<f64>,
    r_star_ij: Vec<f64>,
    well_depth: Vec<f64>,
}

impl VdWaals {
    fn start(mol: &ROMol, conf_id: i32, cutoff: f64, kind: VdwKind) -> Result<Self> {
        let conf = mol.conformer(conf_id)?;
        let mut pos = Vec::with_capacity(3 * mol.num_atoms() as usize);
        for i in 0..mol.num_atoms() {
            let p = conf.atom_pos(i)?;
            pos.extend([p.x, p.y, p.z]);
        }
        Ok(VdWaals {
            kind,
            n_atoms: mol.num_atoms() as usize,
            cutoff: (cutoff * cutoff).max(MIN_CUTOFF_VAL),
            pos,
            r_star_ij: Vec::new(),
            well_depth: Vec::new(),
        })
    }

    /// `MMFFVdWaals`: `probe_atom_type` is an MMFF94 atom type (6: sp3 O);
    /// `scaling` accounts for undirected hydrogen bonds.
    pub fn mmff(
        mol: &ROMol,
        conf_id: i32,
        probe_atom_type: u32,
        scaling: bool,
        cutoff: f64,
    ) -> Result<Self> {
        use crate::forcefield::mmff::params::MmffParams;
        use crate::forcefield::mmff::terms::{
            calc_unscaled_vdw_minimum, calc_unscaled_vdw_well_depth, scale_vdw_params,
        };
        let mut res = Self::start(mol, conf_id, cutoff, VdwKind::Mmff)?;
        let mut copy = mol.clone();
        let props = MmffMolProperties::new(&mut copy, "MMFF94")?;
        if !props.is_valid() {
            return Err(RDKitError::Value(
                "No MMFF atom types available for at least one atom in molecule.".into(),
            ));
        }
        let params = MmffParams::get(false);
        let probe = params
            .vdw(probe_atom_type)
            .ok_or_else(|| RDKitError::Value("Unknown MMFF probe atom type.".into()))?;
        for i in 0..mol.num_atoms() {
            let p = params
                .vdw(u32::from(props.atom_type(i)))
                .ok_or_else(|| RDKitError::Value("Unknown MMFF atom type.".into()))?;
            let mut r = calc_unscaled_vdw_minimum(params.vdw_constants(), p, probe);
            let mut depth = calc_unscaled_vdw_well_depth(r, p, probe);
            if scaling {
                scale_vdw_params(&mut r, &mut depth, params.vdw_constants(), p, probe);
            }
            res.r_star_ij.push(r);
            res.well_depth.push(depth);
        }
        Ok(res)
    }

    /// `UFFVdWaals`: `probe_atom_type` is a UFF atom label ("O_3").
    pub fn uff(mol: &ROMol, conf_id: i32, probe_atom_type: &str, cutoff: f64) -> Result<Self> {
        use crate::forcefield::uff::ParamCollection;
        use crate::forcefield::uff::calc_nonbonded_depth;
        let mut res = Self::start(mol, conf_id, cutoff, VdwKind::Uff)?;
        let probe = ParamCollection::default_params()
            .get(probe_atom_type)
            .copied()
            .ok_or_else(|| RDKitError::Value("Unknown UFF probe atom type.".into()))?;
        let (params, found_all) = get_atom_types(mol)?;
        if !found_all {
            return Err(RDKitError::Value(
                "No UFF atom types available for at least one atom in molecule.".into(),
            ));
        }
        for p in params.into_iter().flatten() {
            res.r_star_ij.push(probe.x1 * p.x1);
            res.well_depth.push(calc_nonbonded_depth(&probe, &p));
        }
        Ok(res)
    }

    fn energy(&self, dist2: f64, r_star_ij: f64, well_depth: f64) -> f64 {
        match self.kind {
            VdwKind::Uff => {
                let mut r6 = r_star_ij / dist2;
                r6 = r6 * r6 * r6;
                well_depth * (r6 * r6 - 2.0 * r6)
            }
            VdwKind::Mmff => {
                crate::forcefield::mmff::terms::calc_vdw_energy(dist2.sqrt(), r_star_ij, well_depth)
            }
        }
    }
}

impl Descriptor for VdWaals {
    fn value(&self, x: f64, y: f64, z: f64, thres: f64) -> f64 {
        let mut res = 0.0;
        for i in 0..self.n_atoms {
            let j = 3 * i;
            let (dx, dy, dz) = (x - self.pos[j], y - self.pos[j + 1], z - self.pos[j + 2]);
            let mut dist2 = dx * dx + dy * dy + dz * dz;
            if dist2 < thres {
                dist2 = dist2.max(self.cutoff);
                res += self.energy(dist2, self.r_star_ij[i], self.well_depth[i]);
            }
        }
        res
    }
}
