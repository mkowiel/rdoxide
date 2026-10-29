//! `rdoxide.Chem.rdEHTTools` (upstream `External/YAeHMOP/Wrap/rdEHTTools.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rdoxide::external::yaehmop::{EhtResults, run_mol};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::pyconv::float64_array;

/// Results from a YAeHMOP calculation
#[pyclass(
    name = "EHTResults",
    module = "rdoxide.Chem.rdEHTTools",
    skip_from_py_object
)]
struct PyEhtResults {
    res: EhtResults,
}

fn matrix<'py>(
    py: Python<'py>,
    values: &[f64],
    rows: usize,
    cols: usize,
) -> PyResult<Bound<'py, PyAny>> {
    float64_array(py, values)?.call_method1("reshape", (rows, cols))
}

#[pymethods]
impl PyEhtResults {
    #[getter(numOrbitals)]
    fn num_orbitals(&self) -> usize {
        self.res.num_orbitals
    }
    #[getter(numElectrons)]
    fn num_electrons(&self) -> usize {
        self.res.num_electrons
    }
    #[getter(fermiEnergy)]
    fn fermi_energy(&self) -> f64 {
        self.res.fermi_energy
    }
    #[getter(totalEnergy)]
    fn total_energy(&self) -> f64 {
        self.res.total_energy
    }

    /// returns the reduced charge matrix
    #[pyo3(name = "GetReducedChargeMatrix")]
    fn reduced_charge_matrix<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        matrix(
            py,
            &self.res.reduced_charge_matrix,
            self.res.num_atoms,
            self.res.num_orbitals,
        )
    }

    /// returns the reduced overlap population matrix
    #[pyo3(name = "GetReducedOverlapPopulationMatrix")]
    fn reduced_overlap_population_matrix<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, PyAny>> {
        float64_array(py, &self.res.reduced_overlap_population_matrix)
    }

    /// returns the calculated atomic charges
    #[pyo3(name = "GetAtomicCharges")]
    fn atomic_charges<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        float64_array(py, &self.res.atomic_charges)
    }

    /// returns the symmetric Hamiltonian matrix
    #[pyo3(name = "GetHamiltonian")]
    fn hamiltonian<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        match &self.res.hamiltonian_matrix {
            Some(h) => matrix(py, h, self.res.num_orbitals, self.res.num_orbitals),
            None => Err(PyValueError::new_err(
                "Hamiltonian not available, set keepOverlapAndHamiltonianMatrices=True to preserve it.",
            )),
        }
    }

    /// returns the symmetric overlap matrix
    #[pyo3(name = "GetOverlapMatrix")]
    fn overlap_matrix<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        match &self.res.overlap_matrix {
            Some(h) => matrix(py, h, self.res.num_orbitals, self.res.num_orbitals),
            None => Err(PyValueError::new_err(
                "Overlap matrix not available, set keepOverlapAndHamiltonianMatrices=True to preserve it.",
            )),
        }
    }

    /// returns the energies of the molecular orbitals as a vector
    #[pyo3(name = "GetOrbitalEnergies")]
    fn orbital_energies<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        float64_array(py, &self.res.orbital_energies)
    }
}

/// Runs an extended Hueckel calculation for a molecule.
///
/// ARGUMENTS:
///    - mol: molecule to use
///    - confId: (optional) conformer ID to use
///    - keepOverlapAndHamiltonianMatrices: (optional) triggers storing the overlap and hamiltonian matrices in the EHTResults object
///
/// RETURNS: a 2-tuple:
///    - a boolean indicating whether or not the calculation succeeded
///    - an EHTResults object with the results
#[pyfunction]
#[pyo3(name = "RunMol", signature = (mol, confId = -1, keepOverlapAndHamiltonianMatrices = false))]
#[allow(non_snake_case)]
fn run_mol_py(
    mol: PyRef<'_, PyMol>,
    confId: i32,
    keepOverlapAndHamiltonianMatrices: bool,
) -> PyResult<(bool, PyEhtResults)> {
    let res = run_mol(&mol.mol, confId, keepOverlapAndHamiltonianMatrices).py()?;
    Ok((true, PyEhtResults { res }))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing interface to the YAeHMOP extended Hueckel library.\n\
         **NOTE** This functionality is experimental and the API\nor results may change from one release to the next.",
    )?;
    m.add_class::<PyEhtResults>()?;
    m.add_function(wrap_pyfunction!(run_mol_py, m)?)?;
    Ok(())
}
