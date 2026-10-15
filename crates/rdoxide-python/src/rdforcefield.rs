//! `rdkit.ForceField.rdForceField` (upstream `ForceField/Wrap/ForceField.cpp`)
//! and `rdkit.Chem.rdForceFieldHelpers` (upstream
//! `ForceFieldHelpers/Wrap/rdForceFields.cpp`), UFF part.
//!
//! A force field made from a molecule stays linked to it: upstream's points
//! are pointers into the conformer, so the force field sees later coordinate
//! changes and `Minimize` moves the atoms of the molecule.

use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::forcefield::ForceField;
use rdoxide::forcefield::constraints::{
    AngleConstraintContrib, DistanceConstraintContrib, PositionConstraintContrib,
    TorsionConstraintContrib,
};
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
    /// the number of extra points after the points of the molecule
    extra_points: usize,
}

impl PyForceField {
    /// takes over the coordinates of the linked conformer
    fn pull(&mut self, py: Python<'_>) -> PyResult<()> {
        if let Some((mol, conf_id)) = &self.link {
            let mol = mol.borrow(py);
            let conf = mol.mol.conformer(*conf_id).py()?;
            self.field
                .set_leading_positions(&flat_positions(conf.positions()))
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

    #[allow(clippy::too_many_arguments)]
    fn distance_constraint(
        &mut self,
        py: Python<'_>,
        idx1: usize,
        idx2: usize,
        relative: bool,
        min_len: f64,
        max_len: f64,
        k: f64,
    ) -> PyResult<()> {
        self.pull(py)?;
        let c =
            DistanceConstraintContrib::new(&self.field, idx1, idx2, relative, min_len, max_len, k)
                .py()?;
        self.field.push_contrib(Box::new(c));
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn angle_constraint(
        &mut self,
        py: Python<'_>,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        relative: bool,
        min_deg: f64,
        max_deg: f64,
        k: f64,
    ) -> PyResult<()> {
        self.pull(py)?;
        let c = AngleConstraintContrib::new(
            &self.field,
            idx1,
            idx2,
            idx3,
            relative,
            min_deg,
            max_deg,
            k,
        )
        .py()?;
        self.field.push_contrib(Box::new(c));
        Ok(())
    }

    fn torsion_constraint(
        &mut self,
        py: Python<'_>,
        idx: [usize; 4],
        relative: bool,
        min_deg: f64,
        max_deg: f64,
        k: f64,
    ) -> PyResult<()> {
        self.pull(py)?;
        let c =
            TorsionConstraintContrib::new(&self.field, idx, relative, min_deg, max_deg, k).py()?;
        self.field.push_contrib(Box::new(c));
        Ok(())
    }

    fn position_constraint(
        &mut self,
        py: Python<'_>,
        idx: usize,
        max_displ: f64,
        k: f64,
    ) -> PyResult<()> {
        self.pull(py)?;
        let c = PositionConstraintContrib::new(&self.field, idx, max_displ, k).py()?;
        self.field.push_contrib(Box::new(c));
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

    /// Adds an extra point, this can be useful for adding constraints.
    #[pyo3(name = "AddExtraPoint", signature = (x, y, z, fixed = true))]
    fn add_extra_point(
        &mut self,
        py: Python<'_>,
        x: f64,
        y: f64,
        z: f64,
        fixed: bool,
    ) -> PyResult<usize> {
        self.pull(py)?;
        let idx = self.field.add_point(&[x, y, z]);
        self.extra_points += 1;
        if fixed {
            self.field.fixed_points_mut().push(idx);
        }
        // like upstream, the number of points rather than the index
        Ok(idx + 1)
    }

    /// returns the location of an extra point as a tuple
    #[pyo3(name = "GetExtraPointPos", signature = (idx))]
    fn get_extra_point_pos<'py>(
        &self,
        py: Python<'py>,
        idx: usize,
    ) -> PyResult<Bound<'py, PyTuple>> {
        if idx >= self.extra_points {
            return Err(PyIndexError::new_err(format!("{idx}")));
        }
        let first = self.field.num_points() - self.extra_points;
        let p = &self.field.positions()[3 * (first + idx)..3 * (first + idx) + 3];
        PyTuple::new(py, p.to_vec())
    }

    /// Adds a distance constraint to the UFF force field (deprecated, use
    /// UFFAddDistanceConstraint instead).
    #[pyo3(name = "AddDistanceConstraint", signature = (idx1, idx2, minLen, maxLen, forceConstant))]
    #[allow(non_snake_case)]
    fn add_distance_constraint(
        &mut self,
        py: Python<'_>,
        idx1: usize,
        idx2: usize,
        minLen: f64,
        maxLen: f64,
        forceConstant: f64,
    ) -> PyResult<()> {
        self.distance_constraint(py, idx1, idx2, false, minLen, maxLen, forceConstant)
    }

    /// Adds a distance constraint to the UFF force field; if relative == True,
    /// then minLen and maxLen are intended as relative to the current distance.
    #[pyo3(name = "UFFAddDistanceConstraint", signature = (idx1, idx2, relative, minLen, maxLen, forceConstant))]
    #[allow(non_snake_case)]
    fn uff_add_distance_constraint(
        &mut self,
        py: Python<'_>,
        idx1: usize,
        idx2: usize,
        relative: bool,
        minLen: f64,
        maxLen: f64,
        forceConstant: f64,
    ) -> PyResult<()> {
        self.distance_constraint(py, idx1, idx2, relative, minLen, maxLen, forceConstant)
    }

    /// Adds a distance constraint to the MMFF force field; if relative == True,
    /// then minLen and maxLen are intended as relative to the current distance.
    #[pyo3(name = "MMFFAddDistanceConstraint", signature = (idx1, idx2, relative, minLen, maxLen, forceConstant))]
    #[allow(non_snake_case)]
    fn mmff_add_distance_constraint(
        &mut self,
        py: Python<'_>,
        idx1: usize,
        idx2: usize,
        relative: bool,
        minLen: f64,
        maxLen: f64,
        forceConstant: f64,
    ) -> PyResult<()> {
        self.distance_constraint(py, idx1, idx2, relative, minLen, maxLen, forceConstant)
    }

    /// Adds an angle constraint to the UFF force field; if relative == True,
    /// then minAngleDeg and maxAngleDeg are intended as relative to the current
    /// angle.
    #[pyo3(name = "UFFAddAngleConstraint", signature = (idx1, idx2, idx3, relative, minAngleDeg, maxAngleDeg, forceConstant))]
    #[allow(non_snake_case)]
    fn uff_add_angle_constraint(
        &mut self,
        py: Python<'_>,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        relative: bool,
        minAngleDeg: f64,
        maxAngleDeg: f64,
        forceConstant: f64,
    ) -> PyResult<()> {
        self.angle_constraint(
            py,
            idx1,
            idx2,
            idx3,
            relative,
            minAngleDeg,
            maxAngleDeg,
            forceConstant,
        )
    }

    /// Adds an angle constraint to the MMFF force field; if relative == True,
    /// then minAngleDeg and maxAngleDeg are intended as relative to the current
    /// angle.
    #[pyo3(name = "MMFFAddAngleConstraint", signature = (idx1, idx2, idx3, relative, minAngleDeg, maxAngleDeg, forceConstant))]
    #[allow(non_snake_case)]
    fn mmff_add_angle_constraint(
        &mut self,
        py: Python<'_>,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        relative: bool,
        minAngleDeg: f64,
        maxAngleDeg: f64,
        forceConstant: f64,
    ) -> PyResult<()> {
        self.angle_constraint(
            py,
            idx1,
            idx2,
            idx3,
            relative,
            minAngleDeg,
            maxAngleDeg,
            forceConstant,
        )
    }

    /// Adds a dihedral angle constraint to the UFF force field; if relative ==
    /// True, then minDihedralDeg and maxDihedralDeg are intended as relative to
    /// the current dihedral angle.
    #[pyo3(name = "UFFAddTorsionConstraint", signature = (idx1, idx2, idx3, idx4, relative, minDihedralDeg, maxDihedralDeg, forceConstant))]
    #[allow(non_snake_case)]
    fn uff_add_torsion_constraint(
        &mut self,
        py: Python<'_>,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        relative: bool,
        minDihedralDeg: f64,
        maxDihedralDeg: f64,
        forceConstant: f64,
    ) -> PyResult<()> {
        self.torsion_constraint(
            py,
            [idx1, idx2, idx3, idx4],
            relative,
            minDihedralDeg,
            maxDihedralDeg,
            forceConstant,
        )
    }

    /// Adds a dihedral angle constraint to the MMFF force field; if relative ==
    /// True, then minDihedralDeg and maxDihedralDeg are intended as relative to
    /// the current dihedral angle.
    #[pyo3(name = "MMFFAddTorsionConstraint", signature = (idx1, idx2, idx3, idx4, relative, minDihedralDeg, maxDihedralDeg, forceConstant))]
    #[allow(non_snake_case)]
    fn mmff_add_torsion_constraint(
        &mut self,
        py: Python<'_>,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        relative: bool,
        minDihedralDeg: f64,
        maxDihedralDeg: f64,
        forceConstant: f64,
    ) -> PyResult<()> {
        self.torsion_constraint(
            py,
            [idx1, idx2, idx3, idx4],
            relative,
            minDihedralDeg,
            maxDihedralDeg,
            forceConstant,
        )
    }

    /// Adds a position constraint to the UFF force field.
    #[pyo3(name = "UFFAddPositionConstraint", signature = (idx, maxDispl, forceConstant))]
    #[allow(non_snake_case)]
    fn uff_add_position_constraint(
        &mut self,
        py: Python<'_>,
        idx: usize,
        maxDispl: f64,
        forceConstant: f64,
    ) -> PyResult<()> {
        self.position_constraint(py, idx, maxDispl, forceConstant)
    }

    /// Adds a position constraint to the MMFF force field.
    #[pyo3(name = "MMFFAddPositionConstraint", signature = (idx, maxDispl, forceConstant))]
    #[allow(non_snake_case)]
    fn mmff_add_position_constraint(
        &mut self,
        py: Python<'_>,
        idx: usize,
        maxDispl: f64,
        forceConstant: f64,
    ) -> PyResult<()> {
        self.position_constraint(py, idx, maxDispl, forceConstant)
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
        extra_points: 0,
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
