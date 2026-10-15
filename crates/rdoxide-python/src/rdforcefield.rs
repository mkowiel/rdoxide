//! `rdkit.ForceField.rdForceField` (upstream `ForceField/Wrap/ForceField.cpp`)
//! and `rdkit.Chem.rdForceFieldHelpers` (upstream
//! `ForceFieldHelpers/Wrap/rdForceFields.cpp`), UFF part.
//!
//! A force field made from a molecule stays linked to it: upstream's points
//! are pointers into the conformer, so the force field sees later coordinate
//! changes and `Minimize` moves the atoms of the molecule.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::forcefield::ForceField;
use rdoxide::graphmol::forcefield_helpers::{
    flat_positions, optimize_molecule, uff, uff_optimize_molecule, uff_optimize_molecule_confs,
};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// A force field
#[pyclass(name = "ForceField", module = "rdkit.ForceField.rdForceField")]
pub struct PyForceField {
    field: ForceField,
    /// the molecule and conformer the points come from
    link: Option<(Py<PyMol>, i32)>,
}

impl PyForceField {
    /// takes over the coordinates of the linked conformer
    fn pull(&mut self, py: Python<'_>) -> PyResult<()> {
        if let Some((mol, conf_id)) = &self.link {
            let mol = mol.borrow(py);
            let conf = mol.mol.conformer(*conf_id).py()?;
            self.field
                .set_positions(&flat_positions(conf.positions()))
                .py()?;
        }
        Ok(())
    }

    /// writes the coordinates to the linked conformer
    fn push(&self, py: Python<'_>) -> PyResult<()> {
        if let Some((mol, conf_id)) = &self.link {
            let mut mol = mol.borrow_mut(py);
            let conf = mol.mol.conformer_mut(*conf_id).py()?;
            for (p, c) in conf
                .positions_mut()
                .iter_mut()
                .zip(self.field.positions().chunks_exact(3))
            {
                *p = rdoxide::geometry::Point3D::new(c[0], c[1], c[2]);
            }
        }
        Ok(())
    }

    fn coords(&self, pos: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Vec<f64>>> {
        let Some(pos) = pos.filter(|p| !p.is_none()) else {
            return Ok(None);
        };
        let expected = self.field.dimension() * self.field.num_points();
        if pos.len()? != expected {
            return Err(PyValueError::new_err(
                "The Python container must have length equal to Dimension() * NumPoints()",
            ));
        }
        let c = pos
            .try_iter()?
            .map(|v| v?.extract())
            .collect::<PyResult<Vec<f64>>>()?;
        Ok(Some(c))
    }
}

#[pymethods]
impl PyForceField {
    /// Returns the energy (in kcal/mol) of the current arrangement
    /// or of the supplied coordinate list (if non-empty)
    #[pyo3(name = "CalcEnergy", signature = (pos = None))]
    fn calc_energy(&mut self, py: Python<'_>, pos: Option<&Bound<'_, PyAny>>) -> PyResult<f64> {
        match self.coords(pos)? {
            Some(c) => Ok(self.field.energy_at(&c)),
            None => {
                self.pull(py)?;
                Ok(self.field.calc_energy())
            }
        }
    }

    /// Returns a tuple filled with the per-coordinate gradients
    /// of the current arrangement or of the supplied coordinate list (if non-empty)
    #[pyo3(name = "CalcGrad", signature = (pos = None))]
    fn calc_grad<'py>(
        &mut self,
        py: Python<'py>,
        pos: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let grad = match self.coords(pos)? {
            Some(c) => {
                let mut g = vec![0.0; c.len()];
                self.field.grad_at(&c, &mut g);
                g
            }
            None => {
                self.pull(py)?;
                self.field.calc_grad()
            }
        };
        PyTuple::new(py, grad)
    }

    /// Returns a tuple filled with the coordinates of the
    /// points the ForceField is handling
    #[pyo3(name = "Positions")]
    fn positions<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        self.pull(py)?;
        PyTuple::new(py, self.field.positions().to_vec())
    }

    /// Returns the dimension of the ForceField
    #[pyo3(name = "Dimension")]
    fn dimension(&self) -> usize {
        self.field.dimension()
    }

    /// Returns the number of points the ForceField is handling
    #[pyo3(name = "NumPoints")]
    fn num_points(&self) -> usize {
        self.field.num_points()
    }

    /// Runs some minimization iterations.
    ///
    ///   Returns 0 if the minimization succeeded.
    #[pyo3(name = "Minimize", signature = (maxIts = 200, forceTol = 1e-4, energyTol = 1e-6))]
    #[allow(non_snake_case)]
    fn minimize(
        &mut self,
        py: Python<'_>,
        maxIts: u32,
        forceTol: f64,
        energyTol: f64,
    ) -> PyResult<i32> {
        self.pull(py)?;
        let field = &mut self.field;
        let res = py
            .detach(|| field.minimize(maxIts, forceTol, energyTol))
            .py()?;
        self.push(py)?;
        Ok(res)
    }

    /// Initializes the force field (the points are read again from the
    /// molecule, so there is nothing else to do).
    #[pyo3(name = "Initialize")]
    fn initialize(&mut self, py: Python<'_>) -> PyResult<()> {
        self.pull(py)
    }

    /// Adds a fixed point to the force field.
    #[pyo3(name = "AddFixedPoint", signature = (idx))]
    fn add_fixed_point(&mut self, idx: usize) {
        self.field.fixed_points_mut().push(idx);
    }
}

fn float_tuple<'py>(py: Python<'py>, values: &[f64]) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(py, values.to_vec())
}

fn opt_tuple(py: Python<'_>, values: Option<&[f64]>) -> PyResult<Py<PyAny>> {
    Ok(match values {
        Some(v) => float_tuple(py, v)?.into_any().unbind(),
        None => py.None(),
    })
}

/// uses UFF to optimize a molecule's structure
///
///  ARGUMENTS:
///
///     - mol : the molecule of interest
///     - maxIters : the maximum number of iterations (defaults to 200)
///     - vdwThresh : used to exclude long-range van der Waals interactions
///                   (defaults to 10.0)
///     - confId : indicates which conformer to optimize
///     - ignoreInterfragInteractions : if true, nonbonded terms between
///                   fragments will not be added to the forcefield.
///
///  RETURNS: 0 if the optimization converged, 1 if more iterations are required.
#[pyfunction]
#[pyo3(name = "UFFOptimizeMolecule", signature = (mol, maxIters = 200, vdwThresh = 10.0, confId = -1, ignoreInterfragInteractions = true))]
#[allow(non_snake_case)]
fn py_uff_optimize_molecule(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    maxIters: u32,
    vdwThresh: f64,
    confId: i32,
    ignoreInterfragInteractions: bool,
) -> PyResult<i32> {
    let mut mol = mol.borrow_mut();
    let mol = &mut mol.mol;
    py.detach(|| {
        uff_optimize_molecule(
            mol,
            maxIters,
            vdwThresh,
            confId,
            ignoreInterfragInteractions,
        )
    })
    .py()
    .map(|(res, _)| res)
}

/// uses UFF to optimize all of a molecule's conformations
///
///  ARGUMENTS:
///
///     - mol : the molecule of interest
///     - numThreads : the number of threads to use, only has an effect if the RDKit
///                    was built with thread support (defaults to 1)
///                    If set to zero, the max supported by the system will be used.
///     - maxIters : the maximum number of iterations (defaults to 200)
///     - vdwThresh : used to exclude long-range van der Waals interactions
///                   (defaults to 10.0)
///     - ignoreInterfragInteractions : if true, nonbonded terms between
///                   fragments will not be added to the forcefield.
///
///  RETURNS: a list of (not_converged, energy) 2-tuples.
///      If not_converged is 0 the optimization converged for that conformer.
#[pyfunction]
#[pyo3(name = "UFFOptimizeMoleculeConfs", signature = (mol, numThreads = 1, maxIters = 200, vdwThresh = 10.0, ignoreInterfragInteractions = true))]
#[allow(non_snake_case)]
fn py_uff_optimize_molecule_confs(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    numThreads: i32,
    maxIters: u32,
    vdwThresh: f64,
    ignoreInterfragInteractions: bool,
) -> PyResult<Vec<(i32, f64)>> {
    let mut mol = mol.borrow_mut();
    let mol = &mut mol.mol;
    py.detach(|| {
        uff_optimize_molecule_confs(
            mol,
            numThreads,
            maxIters,
            vdwThresh,
            ignoreInterfragInteractions,
        )
    })
    .py()
}

/// returns a UFF force field for a molecule
///
///  ARGUMENTS:
///
///     - mol : the molecule of interest
///     - vdwThresh : used to exclude long-range van der Waals interactions
///                   (defaults to 10.0)
///     - confId : indicates which conformer to optimize
///     - ignoreInterfragInteractions : if true, nonbonded terms between
///                   fragments will not be added to the forcefield.
#[pyfunction]
#[pyo3(name = "UFFGetMoleculeForceField", signature = (mol, vdwThresh = 10.0, confId = -1, ignoreInterfragInteractions = true))]
#[allow(non_snake_case)]
fn py_uff_get_molecule_force_field(
    mol: &Bound<'_, PyMol>,
    vdwThresh: f64,
    confId: i32,
    ignoreInterfragInteractions: bool,
) -> PyResult<PyForceField> {
    let field = uff::construct_force_field(
        &mol.borrow().mol,
        vdwThresh,
        confId,
        ignoreInterfragInteractions,
    )
    .py()?;
    Ok(PyForceField {
        field,
        link: Some((mol.clone().unbind(), confId)),
    })
}

/// checks if UFF parameters are available for all of a molecule's atoms
///
///  ARGUMENTS:
///
///     - mol : the molecule of interest.
#[pyfunction]
#[pyo3(name = "UFFHasAllMoleculeParams", signature = (mol))]
fn py_uff_has_all_molecule_params(mol: &PyMol) -> PyResult<bool> {
    Ok(uff::get_atom_types(&mol.mol).py()?.1)
}

/// Retrieves UFF bond stretch parameters for atoms with indices idx1, idx2 as a
/// (kb, r0) tuple, or None if no parameters could be found
#[pyfunction]
#[pyo3(name = "GetUFFBondStretchParams", signature = (mol, idx1, idx2))]
fn py_uff_bond_params(py: Python<'_>, mol: &PyMol, idx1: u32, idx2: u32) -> PyResult<Py<PyAny>> {
    let p = uff::get_uff_bond_stretch_params(&mol.mol, idx1, idx2).py()?;
    opt_tuple(py, p.map(|b| [b.kb, b.r0]).as_ref().map(|a| a.as_slice()))
}

/// Retrieves UFF angle bend parameters for atoms with indices idx1, idx2, idx3 as
/// a (ka, theta0) tuple, or None if no parameters could be found
#[pyfunction]
#[pyo3(name = "GetUFFAngleBendParams", signature = (mol, idx1, idx2, idx3))]
fn py_uff_angle_params(
    py: Python<'_>,
    mol: &PyMol,
    idx1: u32,
    idx2: u32,
    idx3: u32,
) -> PyResult<Py<PyAny>> {
    let p = uff::get_uff_angle_bend_params(&mol.mol, idx1, idx2, idx3).py()?;
    opt_tuple(
        py,
        p.map(|a| [a.ka, a.theta0]).as_ref().map(|a| a.as_slice()),
    )
}

/// Retrieves UFF torsion parameters for atoms with indices idx1, idx2, idx3, idx4
/// as a V float value, or None if no parameters could be found
#[pyfunction]
#[pyo3(name = "GetUFFTorsionParams", signature = (mol, idx1, idx2, idx3, idx4))]
fn py_uff_torsion_params(
    mol: &PyMol,
    idx1: u32,
    idx2: u32,
    idx3: u32,
    idx4: u32,
) -> PyResult<Option<f64>> {
    uff::get_uff_torsion_params(&mol.mol, idx1, idx2, idx3, idx4).py()
}

/// Retrieves UFF inversion parameters for atoms with indices idx1, idx2, idx3,
/// idx4 as a K float value, or None if no parameters could be found
#[pyfunction]
#[pyo3(name = "GetUFFInversionParams", signature = (mol, idx1, idx2, idx3, idx4))]
fn py_uff_inversion_params(
    mol: &PyMol,
    idx1: u32,
    idx2: u32,
    idx3: u32,
    idx4: u32,
) -> PyResult<Option<f64>> {
    uff::get_uff_inversion_params(&mol.mol, idx1, idx2, idx3, idx4).py()
}

/// Retrieves UFF van der Waals parameters for atoms with indices idx1, idx2 as a
/// (x_ij, D_ij) tuple, or None if no parameters could be found
#[pyfunction]
#[pyo3(name = "GetUFFVdWParams", signature = (mol, idx1, idx2))]
fn py_uff_vdw_params(py: Python<'_>, mol: &PyMol, idx1: u32, idx2: u32) -> PyResult<Py<PyAny>> {
    let p = uff::get_uff_vdw_params(&mol.mol, idx1, idx2).py()?;
    opt_tuple(py, p.map(|(x, d)| [x, d]).as_ref().map(|a| a.as_slice()))
}

/// `rdkit.ForceField.rdForceField`
pub fn init_force_field(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__doc__", "Exposes the ForceField class")?;
    m.add_class::<PyForceField>()?;
    Ok(())
}

/// `rdkit.Chem.rdForceFieldHelpers`
pub fn init_helpers(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions to handle force fields",
    )?;
    m.add_function(wrap_pyfunction!(py_uff_optimize_molecule, m)?)?;
    m.add_function(wrap_pyfunction!(py_uff_optimize_molecule_confs, m)?)?;
    m.add_function(wrap_pyfunction!(py_uff_get_molecule_force_field, m)?)?;
    m.add_function(wrap_pyfunction!(py_uff_has_all_molecule_params, m)?)?;
    m.add_function(wrap_pyfunction!(py_uff_bond_params, m)?)?;
    m.add_function(wrap_pyfunction!(py_uff_angle_params, m)?)?;
    m.add_function(wrap_pyfunction!(py_uff_torsion_params, m)?)?;
    m.add_function(wrap_pyfunction!(py_uff_inversion_params, m)?)?;
    m.add_function(wrap_pyfunction!(py_uff_vdw_params, m)?)?;
    // referenced so that optimize_molecule stays linked for MMFF later
    let _ = optimize_molecule;
    Ok(())
}
