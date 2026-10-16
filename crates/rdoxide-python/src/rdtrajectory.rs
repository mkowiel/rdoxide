//! `rdoxide.Chem.rdtrajectory` (upstream `GraphMol/Wrap/{rdTrajectory,Trajectory}.cpp`):
//! trajectories and snapshots.

use pyo3::prelude::*;
use pyo3::types::PyList;
use rdoxide::graphmol::trajectory::{
    self, Snapshot, Trajectory, read_amber_trajectory, read_gromos_trajectory,
};

use crate::errors::IntoPyResult;
use crate::geometry::{PyPoint2D, PyPoint3D};
use crate::mol::PyMol;

/// A class which allows storing coordinates from a trajectory
#[pyclass(
    name = "Snapshot",
    module = "rdoxide.Chem.rdtrajectory",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
pub struct PySnapshot {
    pub snapshot: Snapshot,
}

impl PySnapshot {
    pub fn wrap(snapshot: Snapshot) -> Self {
        PySnapshot { snapshot }
    }
}

#[pymethods]
impl PySnapshot {
    /// Constructor;
    /// coordList: list of floats containing the coordinates for this Snapshot;
    /// energy:    the energy for this Snapshot.
    ///
    /// Alternatively, a copy constructor taking another Snapshot.
    #[new]
    #[pyo3(signature = (coordList, energy = 0.0))]
    #[allow(non_snake_case)]
    fn py_new(coordList: &Bound<'_, PyAny>, energy: f64) -> PyResult<Self> {
        if let Ok(other) = coordList.extract::<PyRef<'_, PySnapshot>>() {
            return Ok(PySnapshot {
                snapshot: other.snapshot.clone(),
            });
        }
        let coords: Vec<f64> = coordList.extract()?;
        Ok(PySnapshot::wrap(Snapshot::new(coords, energy)))
    }

    /// return the coordinates at pointNum as a Point2D object; requires the Trajectory dimension to be == 2
    #[pyo3(name = "GetPoint2D", signature = (pointNum))]
    #[allow(non_snake_case)]
    fn point_2d(&self, pointNum: u32) -> PyResult<PyPoint2D> {
        Ok(PyPoint2D {
            p: self.snapshot.point_2d(pointNum).py()?,
        })
    }

    /// return the coordinates at pointNum as a Point3D object; requires the Trajectory dimension to be >= 2
    #[pyo3(name = "GetPoint3D", signature = (pointNum))]
    #[allow(non_snake_case)]
    fn point_3d(&self, pointNum: u32) -> PyResult<PyPoint3D> {
        Ok(PyPoint3D {
            p: self.snapshot.point_3d(pointNum).py()?,
        })
    }

    /// returns the energy for this Snapshot
    #[pyo3(name = "GetEnergy")]
    fn energy(&self) -> f64 {
        self.snapshot.energy()
    }

    /// sets the energy for this Snapshot
    #[pyo3(name = "SetEnergy", signature = (energy))]
    fn set_energy(&mut self, energy: f64) {
        self.snapshot.set_energy(energy);
    }
}

/// A class which allows storing Snapshots from a trajectory
#[pyclass(
    name = "Trajectory",
    module = "rdoxide.Chem.rdtrajectory",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
pub struct PyTrajectory {
    traj: Trajectory,
}

#[pymethods]
impl PyTrajectory {
    /// Constructor;
    /// dimension:    dimensionality of this Trajectory's coordinate tuples;
    /// numPoints:    number of coordinate tuples associated to each Snapshot;
    /// snapshotList: list of Snapshot objects used to initialize the Trajectory (optional; defaults to []).
    ///
    /// Alternatively, a copy constructor taking another Trajectory.
    #[new]
    #[pyo3(signature = (dimension, numPoints = None, snapshotList = None))]
    #[allow(non_snake_case)]
    fn py_new(
        dimension: &Bound<'_, PyAny>,
        numPoints: Option<u32>,
        snapshotList: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        if let Ok(other) = dimension.extract::<PyRef<'_, PyTrajectory>>() {
            return Ok(PyTrajectory {
                traj: other.traj.clone(),
            });
        }
        let dimension: u32 = dimension.extract()?;
        let num_points = numPoints.ok_or_else(|| {
            pyo3::exceptions::PyTypeError::new_err("Trajectory() missing argument 'numPoints'")
        })?;
        let mut traj = Trajectory::new(dimension, num_points);
        if let Some(list) = snapshotList {
            for s in list.try_iter()? {
                traj.add_snapshot(s?.extract::<PyRef<'_, PySnapshot>>()?.snapshot.clone());
            }
        }
        Ok(PyTrajectory { traj })
    }

    /// returns the dimensionality of this Trajectory's coordinate tuples
    #[pyo3(name = "Dimension")]
    fn dimension(&self) -> u32 {
        self.traj.dimension()
    }

    /// returns the number of coordinate tuples associated to each Snapshot
    #[pyo3(name = "NumPoints")]
    fn num_points(&self) -> u32 {
        self.traj.num_points()
    }

    fn __len__(&self) -> usize {
        self.traj.len()
    }

    /// appends Snapshot s to this Trajectory; returns the zero-based index position of the added snapshot
    #[pyo3(name = "AddSnapshot", signature = (s))]
    fn add_snapshot(&mut self, s: PyRef<'_, PySnapshot>) -> u32 {
        self.traj.add_snapshot(s.snapshot.clone())
    }

    /// returns the Snapshot snapshotNum, where the latter is the zero-based index of the retrieved Snapshot
    #[pyo3(name = "GetSnapshot", signature = (snapshotNum))]
    #[allow(non_snake_case)]
    fn get_snapshot(&self, snapshotNum: u32) -> PyResult<PySnapshot> {
        Ok(PySnapshot::wrap(
            self.traj.snapshot(snapshotNum).py()?.clone(),
        ))
    }

    /// inserts Snapshot s into the Trajectory at the position snapshotNum, where the latter is the zero-based index of the Trajectory's Snapshot before which the Snapshot s will be inserted; returns the zero-based index position of the inserted snapshot
    #[pyo3(name = "InsertSnapshot", signature = (snapshotNum, s))]
    #[allow(non_snake_case)]
    fn insert_snapshot(&mut self, snapshotNum: u32, s: PyRef<'_, PySnapshot>) -> PyResult<u32> {
        self.traj
            .insert_snapshot(snapshotNum, s.snapshot.clone())
            .py()
    }

    /// removes Snapshot snapshotNum from the Trajectory, where snapshotNum is the zero-based index of Snapshot to be removed
    #[pyo3(name = "RemoveSnapshot", signature = (snapshotNum))]
    #[allow(non_snake_case)]
    fn remove_snapshot(&mut self, snapshotNum: u32) -> PyResult<u32> {
        self.traj.remove_snapshot(snapshotNum).py()
    }

    /// removes all Snapshots from the Trajectory
    #[pyo3(name = "Clear")]
    fn clear(&mut self) {
        self.traj.clear();
    }

    /// adds conformations from the Trajectory to mol
    /// fromCid is the first Snapshot that will be added as a Conformer; defaults to -1 (first available)
    /// toCid is the last Snapshot that will be added as a Conformer; defaults to -1 (all)
    #[pyo3(name = "AddConformersToMol", signature = (mol, fromCid = -1, toCid = -1))]
    #[allow(non_snake_case)]
    fn add_conformers_to_mol(
        &self,
        mol: &Bound<'_, PyMol>,
        fromCid: i32,
        toCid: i32,
    ) -> PyResult<u32> {
        self.traj
            .add_conformers_to_mol(&mut mol.borrow_mut().mol, fromCid, toCid)
            .py()
    }
}

/// reads coordinates from an AMBER trajectory file into the Trajectory object; returns the number of Snapshot objects read in
#[pyfunction]
#[pyo3(name = "ReadAmberTrajectory", signature = (fName, traj))]
#[allow(non_snake_case)]
fn py_read_amber(fName: &str, traj: &Bound<'_, PyTrajectory>) -> PyResult<u32> {
    read_amber_trajectory(fName, &mut traj.borrow_mut().traj).py()
}

/// reads coordinates from a GROMOS trajectory file into the Trajectory object; returns the number of Snapshot objects read in
#[pyfunction]
#[pyo3(name = "ReadGromosTrajectory", signature = (fName, traj))]
#[allow(non_snake_case)]
fn py_read_gromos(fName: &str, traj: &Bound<'_, PyTrajectory>) -> PyResult<u32> {
    read_gromos_trajectory(fName, &mut traj.borrow_mut().traj).py()
}

/// A Python list of snapshots.
pub fn snapshot_list<'py>(
    py: Python<'py>,
    snapshots: Vec<rdoxide::numerics::bfgs::Snapshot>,
) -> PyResult<Bound<'py, PyList>> {
    let list = PyList::empty(py);
    for s in snapshots {
        list.append(Py::new(
            py,
            PySnapshot::wrap(trajectory::Snapshot::from(s)),
        )?)?;
    }
    Ok(list)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.setattr(
        "__doc__",
        "Module containing Trajectory and Snapshot objects\n",
    )?;
    m.add_class::<PyTrajectory>()?;
    m.add_class::<PySnapshot>()?;
    m.add_function(wrap_pyfunction!(py_read_amber, m)?)?;
    m.add_function(wrap_pyfunction!(py_read_gromos, m)?)?;
    Ok(())
}
