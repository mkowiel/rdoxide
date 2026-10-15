//! Trajectories: series of coordinate snapshots (upstream
//! `GraphMol/Trajectory/{Snapshot.h,Trajectory.h,Trajectory.cpp}`), as
//! produced by force-field minimizations or read from AMBER and GROMOS
//! files.

use std::collections::BTreeSet;
use std::sync::Arc;

use super::conformer::Conformer;
use super::romol::ROMol;
use crate::geometry::{Point2D, Point3D};
use crate::numerics::bfgs;
use crate::{RDKitError, Result};

/// `Snapshot`: the coordinates (`numPoints * dimension` values) and energy
/// of one step. Once a snapshot is in a [`Trajectory`] it knows the
/// trajectory's shape and can return points.
#[derive(Debug, Clone, PartialEq)]
pub struct Snapshot {
    /// (dimension, number of points) of the trajectory it belongs to
    shape: Option<(u32, u32)>,
    energy: f64,
    pos: Arc<[f64]>,
}

impl Snapshot {
    pub fn new(pos: Vec<f64>, energy: f64) -> Self {
        Snapshot {
            shape: None,
            energy,
            pos: pos.into(),
        }
    }

    pub fn energy(&self) -> f64 {
        self.energy
    }

    pub fn set_energy(&mut self, energy: f64) {
        self.energy = energy;
    }

    /// The coordinates, all points one after the other.
    pub fn coordinates(&self) -> &[f64] {
        &self.pos
    }

    fn checked_shape(&self, point_num: u32, min_dimension: u32) -> Result<(u32, usize)> {
        let Some((dimension, num_points)) = self.shape else {
            return Err(RDKitError::precondition("d_trajectory must not be NULL"));
        };
        if self.pos.is_empty() {
            return Err(RDKitError::precondition("d_pos must not be NULL"));
        }
        if dimension < min_dimension {
            return Err(RDKitError::precondition(
                if min_dimension == 2 && dimension != 2 {
                    "d_dimension must be >= 2"
                } else {
                    "d_dimension must be == 2"
                },
            ));
        }
        if num_points == 0 {
            return Err(RDKitError::precondition("d_numPoints must be > 0"));
        }
        if point_num >= num_points {
            return Err(RDKitError::range("pointNum", point_num, num_points));
        }
        Ok((dimension, (point_num * dimension) as usize))
    }

    /// `getPoint2D`; the trajectory dimension must be 2.
    pub fn point_2d(&self, point_num: u32) -> Result<Point2D> {
        let (dimension, i) = self.checked_shape(point_num, 2)?;
        if dimension != 2 {
            return Err(RDKitError::precondition("d_dimension must be == 2"));
        }
        Ok(Point2D::new(self.pos[i], self.pos[i + 1]))
    }

    /// `getPoint3D`; the trajectory dimension must be at least 2.
    pub fn point_3d(&self, point_num: u32) -> Result<Point3D> {
        let (dimension, i) = self.checked_shape(point_num, 2)?;
        let z = if dimension == 3 { self.pos[i + 2] } else { 0.0 };
        Ok(Point3D::new(self.pos[i], self.pos[i + 1], z))
    }
}

impl From<bfgs::Snapshot> for Snapshot {
    fn from(s: bfgs::Snapshot) -> Self {
        Snapshot::new(s.coords, s.energy)
    }
}

/// `Trajectory`
#[derive(Debug, Clone)]
pub struct Trajectory {
    dimension: u32,
    num_points: u32,
    snapshots: Vec<Snapshot>,
}

impl Trajectory {
    pub fn new(dimension: u32, num_points: u32) -> Self {
        Trajectory {
            dimension,
            num_points,
            snapshots: Vec::new(),
        }
    }

    pub fn with_snapshots(dimension: u32, num_points: u32, snapshots: Vec<Snapshot>) -> Self {
        let mut traj = Trajectory::new(dimension, num_points);
        for s in snapshots {
            traj.add_snapshot(s);
        }
        traj
    }

    pub fn dimension(&self) -> u32 {
        self.dimension
    }

    pub fn num_points(&self) -> u32 {
        self.num_points
    }

    pub fn len(&self) -> usize {
        self.snapshots.len()
    }

    pub fn is_empty(&self) -> bool {
        self.snapshots.is_empty()
    }

    /// `addSnapshot`: returns the index of the added snapshot.
    pub fn add_snapshot(&mut self, s: Snapshot) -> u32 {
        self.insert_snapshot(self.snapshots.len() as u32, s)
            .expect("appending is always in range")
    }

    /// `getSnapshot`
    pub fn snapshot(&self, snapshot_num: u32) -> Result<&Snapshot> {
        self.snapshots.get(snapshot_num as usize).ok_or_else(|| {
            RDKitError::range("snapshotNum", snapshot_num, self.snapshots.len() as u32)
        })
    }

    /// `insertSnapshot`: inserts before `snapshot_num`, returns its index.
    pub fn insert_snapshot(&mut self, snapshot_num: u32, mut s: Snapshot) -> Result<u32> {
        if snapshot_num as usize > self.snapshots.len() {
            return Err(RDKitError::range(
                "snapshotNum",
                snapshot_num,
                self.snapshots.len() as u32 + 1,
            ));
        }
        s.shape = Some((self.dimension, self.num_points));
        self.snapshots.insert(snapshot_num as usize, s);
        Ok(snapshot_num)
    }

    /// `removeSnapshot`: returns the index of the snapshot after the removed one.
    pub fn remove_snapshot(&mut self, snapshot_num: u32) -> Result<u32> {
        if snapshot_num as usize >= self.snapshots.len() {
            return Err(RDKitError::range(
                "snapshotNum",
                snapshot_num,
                self.snapshots.len() as u32,
            ));
        }
        self.snapshots.remove(snapshot_num as usize);
        Ok(snapshot_num)
    }

    /// `clear`
    pub fn clear(&mut self) {
        self.snapshots.clear();
    }

    /// `addConformersToMol`: adds the snapshots `from..=to` (-1: first/last)
    /// as conformers; returns how many were added.
    pub fn add_conformers_to_mol(&self, mol: &mut ROMol, from: i32, to: i32) -> Result<u32> {
        if self.num_points != mol.num_atoms() {
            return Err(RDKitError::precondition(
                "Number of atom mismatch between ROMol and Trajectory",
            ));
        }
        let size = self.len() as i32;
        if from >= size {
            return Err(RDKitError::precondition("from must be < size()"));
        }
        if to >= size {
            return Err(RDKitError::precondition("to must be < size()"));
        }
        let from = from.max(0);
        let to = if to < 0 { size - 1 } else { to };
        if size > 0 && from > to {
            return Err(RDKitError::precondition("from must be <= to"));
        }
        let mut added = 0;
        if size > 0 {
            for n in from..=to {
                let snapshot = &self.snapshots[n as usize];
                let mut conf = Conformer::new(mol.num_atoms() as usize);
                for i in 0..mol.num_atoms() {
                    conf.set_atom_pos(i, snapshot.point_3d(i)?)?;
                }
                mol.add_conformer(conf, true)?;
                added += 1;
            }
        }
        Ok(added)
    }
}

/// A `std::istream` positioned on text, with its eof flag.
struct Stream<'a> {
    data: &'a str,
    pos: usize,
    eof: bool,
}

impl<'a> Stream<'a> {
    fn new(data: &'a str) -> Self {
        Stream {
            data,
            pos: 0,
            eof: false,
        }
    }

    /// `std::getline`: the next line; eof is set when it has no newline.
    fn get_line(&mut self) -> String {
        let rest = &self.data[self.pos..];
        match rest.find('\n') {
            Some(n) => {
                self.pos += n + 1;
                rest[..n].to_string()
            }
            None => {
                self.pos = self.data.len();
                self.eof = true;
                rest.to_string()
            }
        }
    }

    /// `stream >> double`: `None` on failure; reading up to the end of the
    /// data sets eof.
    fn read_f64(&mut self) -> Option<f64> {
        let rest = &self.data[self.pos..];
        let trimmed = rest.trim_start();
        self.pos += rest.len() - trimmed.len();
        if trimmed.is_empty() {
            self.eof = true;
            return None;
        }
        let end = trimmed.find(char::is_whitespace).unwrap_or(trimmed.len());
        self.pos += end;
        if end == trimmed.len() {
            self.eof = true;
        }
        trimmed[..end].parse().ok()
    }
}

fn read_trajectory_file(name: &str) -> Result<String> {
    std::fs::read_to_string(name)
        .map_err(|_| RDKitError::BadFile(format!("Bad input file: {name}")))
}

/// `readAmberTrajectory`: reads the snapshots of an AMBER trajectory file
/// into `traj`; returns how many were read.
pub fn read_amber_trajectory(file_name: &str, traj: &mut Trajectory) -> Result<u32> {
    if traj.dimension() != 3 {
        return Err(RDKitError::precondition(
            "The trajectory must have dimension == 3",
        ));
    }
    let text = read_trajectory_file(file_name)?;
    let mut stream = Stream::new(&text);
    stream.get_line(); // title
    let n_coords = (traj.num_points() * 3) as usize;
    let mut n_snapshots = 0;
    while !stream.eof {
        let mut c = vec![0.0; n_coords];
        let mut i = 0;
        while i < n_coords {
            match stream.read_f64() {
                Some(v) => c[i] = v,
                None => {
                    if !stream.eof {
                        return Err(RDKitError::Value(format!(
                            "Error while reading file: {file_name}"
                        )));
                    }
                    if i > 0 && i + 1 < n_coords {
                        return Err(RDKitError::Value(format!(
                            "Premature end of file: {file_name}"
                        )));
                    }
                    break;
                }
            }
            i += 1;
        }
        if !stream.eof {
            traj.add_snapshot(Snapshot::new(c, 0.0));
            n_snapshots += 1;
        }
    }
    Ok(n_snapshots)
}

/// `readGromosTrajectory`: reads the `POSITION`/`POSITIONRED` blocks of a
/// GROMOS trajectory file (coordinates in nm, stored in Angstrom) into
/// `traj`; returns how many snapshots were read.
pub fn read_gromos_trajectory(file_name: &str, traj: &mut Trajectory) -> Result<u32> {
    if traj.dimension() != 3 {
        return Err(RDKitError::precondition(
            "The trajectory must have dimension == 3",
        ));
    }
    let text = read_trajectory_file(file_name)?;
    let mut stream = Stream::new(&text);
    let n_coords = (traj.num_points() * 3) as usize;
    let ignored: BTreeSet<&str> = [
        "TITLE",
        "TIMESTEP",
        "VELOCITYRED",
        "VELOCITY",
        "GENBOX",
        "BOX",
    ]
    .into_iter()
    .collect();
    let mut n_snapshots = 0;
    while !stream.eof {
        let line = stream.get_line();
        if stream.eof {
            continue;
        }
        if ignored.contains(line.as_str()) {
            let mut line = line;
            while !stream.eof && line != "END" {
                line = stream.get_line();
            }
        } else if line == "POSITIONRED" || line == "POSITION" {
            let mut c = vec![0.0; n_coords];
            let mut j = 0;
            let mut i = 0;
            while i < traj.num_points() {
                let line = stream.get_line();
                if stream.eof || line == "END" {
                    return Err(RDKitError::Value("Wrong number of coordinates".into()));
                }
                if line.contains('#') {
                    continue;
                }
                let values: Vec<f64> = line
                    .split_whitespace()
                    .take(3)
                    .map(|t| t.parse::<f64>())
                    .collect::<std::result::Result<_, _>>()
                    .map_err(|_| RDKitError::Value("Error while reading file".into()))?;
                if values.len() != 3 {
                    return Err(RDKitError::Value("Error while reading file".into()));
                }
                for v in values {
                    c[j] = v * 10.0;
                    j += 1;
                }
                i += 1;
            }
            let end = stream.get_line();
            if stream.eof || end != "END" {
                return Err(RDKitError::Value("Wrong number of coordinates".into()));
            }
            traj.add_snapshot(Snapshot::new(c, 0.0));
            n_snapshots += 1;
        } else {
            let mut supported = String::from("POSITIONRED, POSITION");
            for kw in &ignored {
                supported.push_str(", ");
                supported.push_str(kw);
            }
            return Err(RDKitError::Value(format!(
                "Unsupported block: {line}. Supported blocks are {supported}"
            )));
        }
    }
    Ok(n_snapshots)
}
