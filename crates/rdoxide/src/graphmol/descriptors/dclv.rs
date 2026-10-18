//! Double cubic lattice volume: Shrake-Rupley surface area and volume of a
//! molecule or protein (upstream `Descriptors/DCLV.cpp`).
//!
//! Ported from RDKit 2026.09.1 `Code/GraphMol/Descriptors/DCLV.cpp`.
//! Copyright (C) 2024 Greg Landrum and other RDKit contributors, contributed
//! by NextMove Software, Cambridge, UK. Licensed under the BSD 3-Clause
//! license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! Frank Eisenhaber, Philip Lijnzaad, Patrick Argos, Chris Sander and Michael
//! Scharf, "The Double Cubic Lattice Method: Efficient Approaches to Numerical
//! Integration of Surface Area and Volume and to Dot Surface Contouring of
//! Molecular Assemblies", J. Comput. Chem. 16, 273-284 (1995).

use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::sync::OnceLock;

use crate::geometry::Point3D;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::ROMol;
use crate::{RDKitError, Result};

/// Voxels per axis.
const VOXORDER: usize = 16;

/// The 320 standard dots and the areas of the polyhedron and of one dot
/// (`DCLV_dots.h`).
struct Dots {
    dots: Vec<Point3D>,
    standard_area: f64,
    dot_area: f64,
}

fn standard_dots() -> &'static Dots {
    static DOTS: OnceLock<Dots> = OnceLock::new();
    DOTS.get_or_init(|| {
        let text = include_str!("../../../data/dclv_dots.tsv");
        let (mut standard_area, mut dot_area) = (0.0, 0.0);
        let mut dots = Vec::new();
        for line in text.lines() {
            if let Some(v) = line.strip_prefix("# standard_area\t") {
                standard_area = v.parse().expect("standard area");
            } else if let Some(v) = line.strip_prefix("# dot_area\t") {
                dot_area = v.parse().expect("dot area");
            } else if !line.starts_with('#') && !line.is_empty() {
                // the header stores `float` literals
                let c: Vec<f64> = line
                    .split('\t')
                    .map(|v| f64::from(v.parse::<f32>().expect("dot coordinate")))
                    .collect();
                dots.push(Point3D::new(c[0], c[1], c[2]));
            }
        }
        Dots {
            dots,
            standard_area,
            dot_area,
        }
    })
}

/// Whether an atom belongs to the solvent or ligand that is left out.
fn is_excluded_atom(atom: &crate::graphmol::atom::Atom, include_ligand: bool) -> bool {
    let Some(info) = atom.monomer_info() else {
        return false;
    };
    let res_name = info.residue_name.as_str();
    let solvent = match res_name.chars().next() {
        Some('D') => matches!(res_name, "DOD" | "D20"),
        Some('H') => matches!(res_name, "HOH" | "H20"),
        Some('S') => matches!(res_name, "SOL" | "SO4" | "SUL"),
        Some('W') => res_name == "WAT",
        Some('T') => res_name == "TIP",
        Some('P') => res_name == "P04",
        _ => false,
    };
    solvent || (!include_ligand && info.pdb.as_ref().is_some_and(|p| p.is_hetero_atom))
}

/// Peter Ertl's definition: O and N, optionally P and S, and attached Hs.
fn is_polar(mol: &ROMol, idx: u32, include_s_and_p: bool, include_hs: bool) -> bool {
    match mol.atom(idx).atomic_num() {
        7 | 8 => true,
        15 | 16 => include_s_and_p,
        1 => {
            include_hs
                && mol
                    .atom_neighbors(idx)
                    .any(|n| is_polar(mol, n, include_s_and_p, include_hs))
        }
        _ => false,
    }
}

/// Atoms sorted into a 16 x 16 x 16 grid over the molecule.
struct Voxels {
    vox: [f64; 3],
    origin: [f64; 3],
    cells: Vec<Vec<u32>>,
}

impl Voxels {
    fn new(min: &Point3D, max: &Point3D, positions: &[Point3D], radii: &[f64]) -> Self {
        let vox = [
            VOXORDER as f64 / ((max.x - min.x) + 0.1),
            VOXORDER as f64 / ((max.y - min.y) + 0.1),
            VOXORDER as f64 / ((max.z - min.z) + 0.1),
        ];
        let origin = [min.x, min.y, min.z];
        let mut cells = vec![Vec::new(); VOXORDER * VOXORDER * VOXORDER];
        for (i, pos) in positions.iter().enumerate() {
            if radii[i] != 0.0 {
                let x = (vox[0] * (pos.x - origin[0])) as usize;
                let y = (vox[1] * (pos.y - origin[1])) as usize;
                let z = (vox[2] * (pos.z - origin[2])) as usize;
                cells[(x * VOXORDER + y) * VOXORDER + z].push(i as u32);
            }
        }
        Voxels { vox, origin, cells }
    }

    /// For each atom with a radius, the atoms close enough to hide part of
    /// its probe-inflated sphere. Like upstream, atoms without radius add no
    /// entry, so the list is indexed by position among the atoms that count.
    fn neighbours(
        &self,
        positions: &[Point3D],
        radii: &[f64],
        probe_radius: f64,
        max_radius: f64,
    ) -> Vec<Vec<u32>> {
        let last = (VOXORDER - 1) as i32;
        let mut nbrs = Vec::new();
        for (atm_idx, pos) in positions.iter().enumerate() {
            if radii[atm_idx] == 0.0 {
                continue;
            }
            let range = radii[atm_idx] + probe_radius + probe_radius;
            let max_dist = range + max_radius;
            let lo =
                |p: f64, k: usize| ((self.vox[k] * (p - max_dist - self.origin[k])) as i32).max(0);
            let hi = |p: f64, k: usize| {
                ((self.vox[k] * (p + max_dist - self.origin[k])) as i32).min(last)
            };
            let (lx, ly, lz) = (lo(pos.x, 0), lo(pos.y, 1), lo(pos.z, 2));
            let (ux, uy, uz) = (hi(pos.x, 0), hi(pos.y, 1), hi(pos.z, 2));
            let mut atom_neighbours = Vec::new();
            for x in lx..=ux {
                for y in ly..=uy {
                    for z in lz..=uz {
                        let cell = &self.cells
                            [((x as usize) * VOXORDER + y as usize) * VOXORDER + z as usize];
                        for &idx in cell {
                            if idx as usize != atm_idx {
                                let dist = range + radii[idx as usize];
                                if (*pos - positions[idx as usize]).length_sq() < dist * dist {
                                    atom_neighbours.push(idx);
                                }
                            }
                        }
                    }
                }
            }
            nbrs.push(atom_neighbours);
        }
        nbrs
    }
}

/// Surface area and volume by the double cubic lattice method.
#[derive(Debug, Clone)]
pub struct DoubleCubicLatticeVolume {
    mol: ROMol,
    radii: Vec<f64>,
    probe_radius: f64,
    positions: Vec<Point3D>,
    centre_of_gravity: Point3D,
    neighbours: Vec<Vec<u32>>,
    /// atom whose sphere hid the previous probe point (`recordCache`)
    record_cache: Option<usize>,
    surface_area: f64,
    total_volume: f64,
    vdw_volume: f64,
    surface_points: BTreeMap<u32, Vec<Point3D>>,
    all_surface_points: BTreeMap<u32, Vec<Point3D>>,
}

impl DoubleCubicLatticeVolume {
    /// `DoubleCubicLatticeVolume(mol, radii, isProtein, includeLigand, probeRadius, confId)`;
    /// without `radii` the van der Waals radii of the elements are used.
    pub fn new(
        mol: &ROMol,
        radii: Vec<f64>,
        is_protein: bool,
        include_ligand: bool,
        probe_radius: f64,
        conf_id: i32,
    ) -> Result<Self> {
        let n = mol.num_atoms() as usize;
        let mut radii = radii;
        if radii.is_empty() {
            let table = PeriodicTable::get();
            radii = mol
                .atoms()
                .iter()
                .map(|a| table.get_rvdw(a.atomic_num()))
                .collect::<Result<_>>()?;
        }
        if radii.len() != n || n == 0 {
            return Err(RDKitError::Value(
                "the number of radii must match the number of atoms".into(),
            ));
        }
        let positions = mol.conformer(conf_id)?.positions().to_vec();
        let max_radius = radii.iter().copied().fold(f64::MIN, f64::max);

        let (mut num_atoms, mut centre) = (0u32, Point3D::new(0.0, 0.0, 0.0));
        let mut bounds: Option<(Point3D, Point3D)> = None;
        for (idx, atom) in mol.atoms().iter().enumerate() {
            num_atoms += 1;
            if is_protein && is_excluded_atom(atom, include_ligand) {
                radii[idx] = 0.0;
                num_atoms -= 1;
            }
            if radii[idx] != 0.0 {
                let p = positions[idx];
                centre.x += p.x;
                centre.y += p.y;
                centre.z += p.z;
                match &mut bounds {
                    None => bounds = Some((p, p)),
                    Some((min, max)) => {
                        // upstream tests each axis with `if ... else if`
                        if p.x > max.x {
                            max.x = p.x;
                        } else if p.x < min.x {
                            min.x = p.x;
                        }
                        if p.y > max.y {
                            max.y = p.y;
                        } else if p.y < min.y {
                            min.y = p.y;
                        }
                        if p.z > max.z {
                            max.z = p.z;
                        } else if p.z < min.z {
                            min.z = p.z;
                        }
                    }
                }
            }
        }
        let centre_of_gravity = Point3D::new(
            centre.x / f64::from(num_atoms),
            centre.y / f64::from(num_atoms),
            centre.z / f64::from(num_atoms),
        );
        let neighbours = match &bounds {
            Some((min, max)) => Voxels::new(min, max, &positions, &radii).neighbours(
                &positions,
                &radii,
                probe_radius,
                max_radius,
            ),
            None => Vec::new(),
        };
        Ok(DoubleCubicLatticeVolume {
            mol: mol.clone(),
            radii,
            probe_radius,
            positions,
            centre_of_gravity,
            neighbours,
            record_cache: None,
            surface_area: 0.0,
            total_volume: 0.0,
            vdw_volume: 0.0,
            surface_points: BTreeMap::new(),
            all_surface_points: BTreeMap::new(),
        })
    }

    pub fn num_atoms(&self) -> u32 {
        self.mol.num_atoms()
    }

    /// Whether the probe of radius `solvrad` at `point` touches none of `nbrs`.
    fn test_point(&mut self, point: &Point3D, solvrad: f64, nbrs: &[u32]) -> bool {
        if let Some(cached) = self.record_cache {
            let dist = self.radii[cached] + solvrad;
            if (self.positions[cached] - *point).length_sq() < dist * dist {
                return false;
            }
            self.record_cache = None;
        }
        for &i in nbrs {
            let i = i as usize;
            let dist = self.radii[i] + solvrad;
            if (self.positions[i] - *point).length_sq() < dist * dist {
                self.record_cache = Some(i);
                return false;
            }
        }
        true
    }

    /// The neighbours of an atom; upstream indexes the compacted list by atom
    /// index, so an atom after masked ones may have none.
    fn neighbours_of(&self, atom_idx: u32) -> Vec<u32> {
        self.neighbours
            .get(atom_idx as usize)
            .cloned()
            .unwrap_or_default()
    }

    fn check_atom(&self, atom_idx: u32) -> Result<()> {
        self.mol.get_atom(atom_idx).map(|_| ())
    }

    /// `getAtomSurfaceArea`
    pub fn atom_surface_area(&mut self, atom_idx: u32) -> Result<f64> {
        self.check_atom(atom_idx)?;
        if let Some(points) = self.surface_points.get_mut(&atom_idx) {
            points.clear();
        }
        let rad = self.radii[atom_idx as usize];
        if rad == 0.0 {
            return Ok(0.0);
        }
        let pos = self.positions[atom_idx as usize];
        let nbr = self.neighbours_of(atom_idx);
        let factor = rad + self.probe_radius;
        self.record_cache = None;
        let dots = standard_dots();
        let mut area = 0.0;
        for dot in &dots.dots {
            let vect = pos + *dot * factor;
            if self.test_point(&vect, self.probe_radius, &nbr) {
                area += dots.dot_area;
                self.surface_points.entry(atom_idx).or_default().push(vect);
            }
        }
        area *= ((4.0 * PI) * (factor * factor)) / dots.standard_area;
        Ok(area)
    }

    /// `getSurfaceArea`: the solvent accessible surface area.
    pub fn surface_area(&mut self) -> Result<f64> {
        if self.surface_area != 0.0 {
            return Ok(self.surface_area);
        }
        for idx in 0..self.num_atoms() {
            if self.radii[idx as usize] != 0.0 {
                self.surface_area += self.atom_surface_area(idx)?;
            }
        }
        Ok(self.surface_area)
    }

    /// `getPartialSurfaceArea`: the surface area of a set of atoms.
    pub fn partial_surface_area(&mut self, included: &[bool]) -> Result<f64> {
        let mut area = 0.0;
        for idx in 0..self.num_atoms() {
            if included[idx as usize] && self.radii[idx as usize] != 0.0 {
                area += self.atom_surface_area(idx)?;
            }
        }
        Ok(area)
    }

    fn polar_atoms(&self, include_s_and_p: bool, include_hs: bool) -> Vec<bool> {
        (0..self.num_atoms())
            .map(|i| is_polar(&self.mol, i, include_s_and_p, include_hs))
            .collect()
    }

    /// `getPolarSurfaceArea`
    pub fn polar_surface_area(&mut self, include_s_and_p: bool, include_hs: bool) -> Result<f64> {
        let polar = self.polar_atoms(include_s_and_p, include_hs);
        self.partial_surface_area(&polar)
    }

    /// `getSurfacePoints`: the dots on the surface, or with `all_points` all
    /// dots around every atom. The result of the first call is kept.
    pub fn surface_points(&mut self, all_points: bool) -> Result<&BTreeMap<u32, Vec<Point3D>>> {
        if !all_points && !self.surface_points.is_empty() {
            return Ok(&self.surface_points);
        }
        if all_points && !self.all_surface_points.is_empty() {
            return Ok(&self.all_surface_points);
        }
        let dots = standard_dots();
        for idx in 0..self.num_atoms() {
            if all_points {
                let points = self.all_surface_points.entry(idx).or_default();
                points.clear();
                let rad = self.radii[idx as usize];
                if rad != 0.0 {
                    let pos = self.positions[idx as usize];
                    let factor = rad + self.probe_radius;
                    points.extend(dots.dots.iter().map(|d| pos + *d * factor));
                }
            } else {
                self.atom_surface_area(idx)?;
            }
        }
        Ok(if all_points {
            &self.all_surface_points
        } else {
            &self.surface_points
        })
    }

    /// `getAtomVolume`: the volume of one atom for a probe of `solvent_radius`.
    pub fn atom_volume(&mut self, atom_idx: u32, solvent_radius: f64) -> Result<f64> {
        self.check_atom(atom_idx)?;
        let rad = self.radii[atom_idx as usize];
        if rad == 0.0 {
            return Ok(0.0);
        }
        let pos = self.positions[atom_idx as usize];
        let nbr = self.neighbours_of(atom_idx);
        let factor = rad + solvent_radius;
        self.record_cache = None;
        let dots = standard_dots();
        let mut p = Point3D::new(0.0, 0.0, 0.0);
        let mut count = 0u32;
        for dot in &dots.dots {
            let vect = pos + *dot * factor;
            if self.test_point(&vect, solvent_radius, &nbr) {
                p += *dot;
                count += 1;
            }
        }
        // the Gauss-Ostrogradskii theorem
        let mut volume = (pos - self.centre_of_gravity).dot_product(&p);
        volume = factor * factor * (volume + factor * f64::from(count));
        volume *= ((4.0 / 3.0) * PI) / dots.dots.len() as f64;
        Ok(volume)
    }

    /// `getVolume`: the volume enclosed by the probe sphere.
    pub fn volume(&mut self) -> Result<f64> {
        if self.total_volume != 0.0 {
            return Ok(self.total_volume);
        }
        for idx in 0..self.num_atoms() {
            if self.radii[idx as usize] != 0.0 {
                self.total_volume += self.atom_volume(idx, self.probe_radius)?;
            }
        }
        Ok(self.total_volume)
    }

    /// `getVDWVolume`: the van der Waals volume (probe radius 0).
    pub fn vdw_volume(&mut self) -> Result<f64> {
        if self.vdw_volume != 0.0 {
            return Ok(self.vdw_volume);
        }
        for idx in 0..self.num_atoms() {
            if self.radii[idx as usize] != 0.0 {
                self.vdw_volume += self.atom_volume(idx, 0.0)?;
            }
        }
        Ok(self.vdw_volume)
    }

    /// `getPartialVolume`: the volume of a set of atoms.
    pub fn partial_volume(&mut self, included: &[bool]) -> Result<f64> {
        let mut volume = 0.0;
        for idx in 0..self.num_atoms() {
            if included[idx as usize] && self.radii[idx as usize] != 0.0 {
                volume += self.atom_volume(idx, self.probe_radius)?;
            }
        }
        Ok(volume)
    }

    /// `getPolarVolume`
    pub fn polar_volume(&mut self, include_s_and_p: bool, include_hs: bool) -> Result<f64> {
        let polar = self.polar_atoms(include_s_and_p, include_hs);
        self.partial_volume(&polar)
    }

    /// `getCompactness`
    pub fn compactness(&mut self) -> Result<f64> {
        let area = self.surface_area()?;
        let volume = self.volume()?;
        Ok(area / (36.0 * PI * volume * volume).cbrt())
    }

    /// `getPackingDensity`
    pub fn packing_density(&mut self) -> Result<f64> {
        let vdw = self.vdw_volume()?;
        let volume = self.volume()?;
        Ok(vdw / volume)
    }

    /// The molecule the calculation was made for.
    pub fn mol(&self) -> &ROMol {
        &self.mol
    }
}
