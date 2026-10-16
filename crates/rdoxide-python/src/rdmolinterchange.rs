//! `rdoxide.Chem.rdMolInterchange` (upstream
//! `GraphMol/MolInterchange/Wrap/rdMolInterchange.cpp`).

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::RDKitError;
use rdoxide::graphmol::mol_interchange::{
    JsonParseParameters, JsonWriteParameters, json_to_mols, mol_to_json, mols_to_json,
};
use rdoxide::graphmol::romol::ROMol;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// `FileParseException` text as the Boost.Python translator reports it.
fn map_error(err: RDKitError) -> PyErr {
    match err {
        RDKitError::FileParse(m) => PyRuntimeError::new_err(format!("File parsing error: {m}")),
        other => Err::<(), _>(other).py().unwrap_err(),
    }
}

/// Parameters controlling the JSON parser
#[pyclass(
    name = "JSONParseParameters",
    module = "rdoxide.Chem.rdMolInterchange",
    skip_from_py_object
)]
#[derive(Clone, Default)]
struct PyJsonParseParameters {
    inner: JsonParseParameters,
}

#[pymethods]
impl PyJsonParseParameters {
    #[new]
    fn new() -> Self {
        PyJsonParseParameters::default()
    }

    /// set bond types to aromatic for bonds flagged aromatic
    #[getter(setAromaticBonds)]
    fn get_set_aromatic_bonds(&self) -> bool {
        self.inner.set_aromatic_bonds
    }
    #[setter(setAromaticBonds)]
    fn set_set_aromatic_bonds(&mut self, v: bool) {
        self.inner.set_aromatic_bonds = v;
    }

    /// be strict when checking atom valences
    #[getter(strictValenceCheck)]
    fn get_strict_valence_check(&self) -> bool {
        self.inner.strict_valence_check
    }
    #[setter(strictValenceCheck)]
    fn set_strict_valence_check(&mut self, v: bool) {
        self.inner.strict_valence_check = v;
    }

    /// parse conformers in the JSON
    #[getter(parseConformers)]
    fn get_parse_conformers(&self) -> bool {
        self.inner.parse_conformers
    }
    #[setter(parseConformers)]
    fn set_parse_conformers(&mut self, v: bool) {
        self.inner.parse_conformers = v;
    }

    /// parse molecular properties in the JSON
    #[getter(parseProperties)]
    fn get_parse_properties(&self) -> bool {
        self.inner.parse_properties
    }
    #[setter(parseProperties)]
    fn set_parse_properties(&mut self, v: bool) {
        self.inner.parse_properties = v;
    }

    /// use atomic H counts from the JSON. You may want to set this to False
    /// when parsing queries.
    #[getter(useHCounts)]
    fn get_use_h_counts(&self) -> bool {
        self.inner.use_h_counts
    }
    #[setter(useHCounts)]
    fn set_use_h_counts(&mut self, v: bool) {
        self.inner.use_h_counts = v;
    }
}

/// Parameters controlling the JSON writer
#[pyclass(
    name = "JSONWriteParameters",
    module = "rdoxide.Chem.rdMolInterchange",
    skip_from_py_object
)]
#[derive(Clone, Default)]
struct PyJsonWriteParameters {
    inner: JsonWriteParameters,
}

#[pymethods]
impl PyJsonWriteParameters {
    #[new]
    fn new() -> Self {
        PyJsonWriteParameters::default()
    }

    /// use RDKit extensions to the commonchem format
    #[getter(useRDKitExtensions)]
    fn get_use_rdkit_extensions(&self) -> bool {
        self.inner.use_rdkit_extensions
    }
    #[setter(useRDKitExtensions)]
    fn set_use_rdkit_extensions(&mut self, v: bool) {
        self.inner.use_rdkit_extensions = v;
    }
}

fn write_params(params: Option<PyRef<'_, PyJsonWriteParameters>>) -> JsonWriteParameters {
    params.map(|p| p.inner).unwrap_or_default()
}

/// Convert a single molecule to JSON
///
///     ARGUMENTS:
///       - mol: the molecule to work with
///     RETURNS:
///       a string
#[pyfunction]
#[pyo3(name = "MolToJSON", signature = (mol, params = None))]
fn py_mol_to_json(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    params: Option<PyRef<'_, PyJsonWriteParameters>>,
) -> PyResult<String> {
    let params = write_params(params);
    let mol: &ROMol = &mol.mol;
    py.detach(|| mol_to_json(mol, &params)).map_err(map_error)
}

/// Convert a set of molecules to JSON
///
///     ARGUMENTS:
///       - mols: the molecules to work with
///     RETURNS:
///       a string
#[pyfunction]
#[pyo3(name = "MolsToJSON", signature = (mols, params = None))]
fn py_mols_to_json(
    py: Python<'_>,
    mols: &Bound<'_, PyAny>,
    params: Option<PyRef<'_, PyJsonWriteParameters>>,
) -> PyResult<String> {
    // an empty (falsy) sequence yields an empty string, as upstream
    if !mols.is_truthy()? {
        return Ok(String::new());
    }
    let held = mols
        .try_iter()?
        .map(|item| {
            let item = item?;
            if item.is_none() {
                Ok(None)
            } else {
                item.extract::<PyRef<'_, PyMol>>()
                    .map(Some)
                    .map_err(PyErr::from)
            }
        })
        .collect::<PyResult<Vec<_>>>()?;
    let refs: Vec<Option<&ROMol>> = held.iter().map(|m| m.as_ref().map(|m| &m.mol)).collect();
    let params = write_params(params);
    py.detach(|| mols_to_json(&refs, &params))
        .map_err(map_error)
}

/// Convert JSON to a tuple of molecules
///
///     ARGUMENTS:
///       - jsonBlock: the molecule to work with
///       - params: (optional) JSONParseParameters controlling the JSON parsing
///     RETURNS:
///       a tuple of Mols
#[pyfunction]
#[pyo3(name = "JSONToMols", signature = (jsonBlock, params = None))]
#[allow(non_snake_case)]
fn py_json_to_mols<'py>(
    py: Python<'py>,
    jsonBlock: &Bound<'py, PyAny>,
    params: Option<PyRef<'_, PyJsonParseParameters>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let params = params.map(|p| p.inner).unwrap_or_default();
    // std::string arguments accept str and bytes
    let text = match jsonBlock.extract::<String>() {
        Ok(s) => s,
        Err(_) => String::from_utf8_lossy(&jsonBlock.extract::<Vec<u8>>()?).into_owned(),
    };
    let mols = py
        .detach(|| json_to_mols(&text, &params))
        .map_err(map_error)?;
    PyTuple::new(py, mols.into_iter().map(|mol| PyMol { mol }))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions for interchange of molecules.\n\
         Note that this should be considered beta and that the format\n  \
         and API will very likely change in future releases.",
    )?;
    m.add_class::<PyJsonParseParameters>()?;
    m.add_class::<PyJsonWriteParameters>()?;
    m.add_function(wrap_pyfunction!(py_mol_to_json, m)?)?;
    m.add_function(wrap_pyfunction!(py_mols_to_json, m)?)?;
    m.add_function(wrap_pyfunction!(py_json_to_mols, m)?)?;
    Ok(())
}
