//! `rdchem.Conformer` — port of upstream `Code/GraphMol/Wrap/Conformer.cpp`.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rdoxide::geometry::Point3D;
use rdoxide::graphmol::Conformer;
use rdoxide::rdgeneral::Dict;

use crate::errors::IntoPyResult;
use crate::geometry::PyPoint3D;
use crate::mol::PyMol;
use crate::props::impl_props;

#[derive(Debug)]
pub enum ConfHandle {
    Owned(Conformer),
    InMol { mol: Py<PyMol>, id: u32 },
}

/// The class to store 2D or 3D conformation of a molecule
#[pyclass(name = "Conformer", module = "rdoxide.Chem.rdchem", subclass)]
#[derive(Debug)]
pub struct PyConformer {
    pub handle: ConfHandle,
}

fn gone() -> PyErr {
    PyValueError::new_err("Bad Conformer Id")
}

impl PyConformer {
    pub fn in_mol(py: Python<'_>, mol: Py<PyMol>, id: u32) -> PyResult<Py<PyConformer>> {
        Py::new(
            py,
            PyConformer {
                handle: ConfHandle::InMol { mol, id },
            },
        )
    }

    pub fn with<R>(&self, py: Python<'_>, f: impl FnOnce(&Conformer) -> R) -> PyResult<R> {
        match &self.handle {
            ConfHandle::Owned(c) => Ok(f(c)),
            ConfHandle::InMol { mol, id } => {
                let m = mol.bind(py).borrow();
                m.mol
                    .conformers()
                    .iter()
                    .find(|c| c.id() == *id)
                    .map(f)
                    .ok_or_else(gone)
            }
        }
    }

    pub fn with_mut<R>(
        &mut self,
        py: Python<'_>,
        f: impl FnOnce(&mut Conformer) -> R,
    ) -> PyResult<R> {
        match &mut self.handle {
            ConfHandle::Owned(c) => Ok(f(c)),
            ConfHandle::InMol { mol, id } => {
                let mut m = mol.bind(py).borrow_mut();
                let id = *id;
                m.mol
                    .conformers_mut()
                    .iter_mut()
                    .find(|c| c.id() == id)
                    .map(f)
                    .ok_or_else(gone)
            }
        }
    }

    pub fn snapshot(&self, py: Python<'_>) -> PyResult<Conformer> {
        self.with(py, Conformer::clone)
    }

    pub fn with_dict<R>(&self, py: Python<'_>, f: impl FnOnce(&Dict) -> R) -> PyResult<R> {
        self.with(py, |c| f(rdoxide::rdgeneral::RDProps::dict(c)))
    }
    pub fn with_dict_mut<R>(
        &mut self,
        py: Python<'_>,
        f: impl FnOnce(&mut Dict) -> R,
    ) -> PyResult<R> {
        self.with_mut(py, |c| f(rdoxide::rdgeneral::RDProps::dict_mut(c)))
    }
}

impl_props!(PyConformer, computed = true, dict_defaults = (false, false));

/// A position given as a Point3D or as a sequence of 2 or 3 numbers.
fn extract_point(obj: &Bound<'_, PyAny>) -> PyResult<Point3D> {
    if let Ok(p) = obj.cast::<PyPoint3D>() {
        return Ok(p.borrow().p);
    }
    let v: Vec<f64> = obj.extract()?;
    match v.len() {
        2 => Ok(Point3D::new(v[0], v[1], 0.0)),
        3 => Ok(Point3D::new(v[0], v[1], v[2])),
        _ => Err(PyValueError::new_err(
            "Position must have 2 or 3 dimensions",
        )),
    }
}

#[derive(FromPyObject)]
enum ConfInit<'py> {
    NumAtoms(u32),
    Other(PyRef<'py, PyConformer>),
}

#[pymethods]
impl PyConformer {
    #[new]
    #[pyo3(signature = (arg = None))]
    fn py_new(py: Python<'_>, arg: Option<ConfInit<'_>>) -> PyResult<Self> {
        let c = match arg {
            None => Conformer::new(0),
            Some(ConfInit::NumAtoms(n)) => Conformer::new(n as usize),
            Some(ConfInit::Other(o)) => o.snapshot(py)?.copy_detached(),
        };
        Ok(PyConformer {
            handle: ConfHandle::Owned(c),
        })
    }

    /// Get the number of atoms in the conformer
    #[pyo3(name = "GetNumAtoms")]
    fn get_num_atoms(&self, py: Python<'_>) -> PyResult<u32> {
        self.with(py, Conformer::num_atoms)
    }

    /// Returns whether or not this instance belongs to a molecule
    #[pyo3(name = "HasOwningMol")]
    fn has_owning_mol(&self) -> bool {
        matches!(self.handle, ConfHandle::InMol { .. })
    }

    /// Get the owning molecule
    #[pyo3(name = "GetOwningMol")]
    fn get_owning_mol(&self, py: Python<'_>) -> PyResult<Py<PyMol>> {
        match &self.handle {
            ConfHandle::InMol { mol, .. } => Ok(mol.clone_ref(py)),
            ConfHandle::Owned(_) => Err(crate::errors::to_py(rdoxide::RDKitError::precondition(
                "no owner",
            ))),
        }
    }

    /// Get the ID of the conformer
    #[pyo3(name = "GetId")]
    fn get_id(&self, py: Python<'_>) -> PyResult<u32> {
        self.with(py, Conformer::id)
    }

    /// Set the ID of the conformer
    #[pyo3(name = "SetId")]
    fn set_id(&mut self, py: Python<'_>, id: u32) -> PyResult<()> {
        self.with_mut(py, |c| c.set_id(id))?;
        if let ConfHandle::InMol { id: hid, .. } = &mut self.handle {
            *hid = id;
        }
        Ok(())
    }

    /// Get the posistion of an atom
    #[pyo3(name = "GetAtomPosition")]
    fn get_atom_position(&self, py: Python<'_>, aid: u32) -> PyResult<PyPoint3D> {
        Ok(PyPoint3D {
            p: self.with(py, |c| c.atom_pos(aid))?.py()?,
        })
    }

    /// Get positions of all the atoms
    #[pyo3(name = "GetPositions")]
    fn get_positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let rows: Vec<[f64; 3]> = self.with(py, |c| {
            c.positions().iter().map(|p| [p.x, p.y, p.z]).collect()
        })?;
        let np = py.import("numpy")?;
        let arr = np.call_method1("array", (rows,))?;
        if arr.getattr("size")?.extract::<usize>()? == 0 {
            return np.call_method1("zeros", ((0usize, 3usize),));
        }
        Ok(arr)
    }

    /// Set positions of all the atoms given a 2D or 3D numpy array of type double
    #[pyo3(name = "SetPositions")]
    fn set_positions(&mut self, py: Python<'_>, positions: &Bound<'_, PyAny>) -> PyResult<()> {
        let rows: Vec<Vec<f64>> = positions
            .call_method0("tolist")
            .unwrap_or_else(|_| positions.clone())
            .extract()?;
        let n = self.with(py, Conformer::num_atoms)? as usize;
        if rows.len() != n {
            return Err(PyValueError::new_err(
                "Number of positions does not match number of atoms",
            ));
        }
        let mut pts = Vec::with_capacity(n);
        for r in rows {
            pts.push(match r.len() {
                2 => Point3D::new(r[0], r[1], 0.0),
                3 => Point3D::new(r[0], r[1], r[2]),
                _ => return Err(PyValueError::new_err("Positions must have 2 or 3 columns")),
            });
        }
        self.with_mut(py, |c| *c.positions_mut() = pts)
    }

    /// Set the position of the specified atom
    #[pyo3(name = "SetAtomPosition")]
    fn set_atom_position(
        &mut self,
        py: Python<'_>,
        aid: u32,
        loc: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let p = extract_point(loc)?;
        self.with_mut(py, |c| c.set_atom_pos(aid, p))?.py()
    }

    /// Set the 3D flag of the conformer
    #[pyo3(name = "Set3D")]
    fn set_3d(&mut self, py: Python<'_>, v: bool) -> PyResult<()> {
        self.with_mut(py, |c| c.set_3d(v))
    }

    /// returns the 3D flag of the conformer
    #[pyo3(name = "Is3D")]
    fn is_3d(&self, py: Python<'_>) -> PyResult<bool> {
        self.with(py, Conformer::is_3d)
    }

    /// Removes all computed properties from the conformer.
    #[pyo3(name = "ClearComputedProps")]
    fn clear_computed_props(&mut self, py: Python<'_>) -> PyResult<()> {
        self.with_dict_mut(py, rdoxide::rdgeneral::RDProps::clear_computed_props)
    }
}
