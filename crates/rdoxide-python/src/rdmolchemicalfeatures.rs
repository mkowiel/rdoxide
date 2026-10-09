//! `rdkit.Chem.rdMolChemicalFeatures` (upstream
//! `MolChemicalFeatures/Wrap/`): feature definitions, the factory that finds
//! the features of a molecule, and the features themselves.

use std::sync::{Arc, Mutex};

use pyo3::exceptions::{PyIOError, PyIndexError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::graphmol::mol_chemical_features::parser::{FeatureError, build_feature_factory};
use rdoxide::graphmol::mol_chemical_features::{MolChemicalFeature, MolChemicalFeatureFactory};

use crate::errors::{IntoPyResult, to_py};
use crate::geometry::PyPoint3D;
use crate::mol::PyMol;

fn feature_error(e: FeatureError) -> PyErr {
    match e {
        FeatureError::Parse(e) => PyValueError::new_err(e.to_string()),
        FeatureError::Rdkit(e) => to_py(e),
    }
}

/// Class to featurize a molecule
#[pyclass(
    name = "MolChemicalFeatureFactory",
    module = "rdkit.Chem.rdMolChemicalFeatures"
)]
#[derive(Debug)]
struct PyFactory {
    factory: Arc<MolChemicalFeatureFactory>,
}

/// Class to represent a chemical feature.
/// These chemical features may or may not have been derived from molecule object;
/// i.e. it is possible to have a chemical feature that was created just from its type
/// and location.
#[pyclass(
    name = "MolChemicalFeature",
    module = "rdkit.Chem.rdMolChemicalFeatures"
)]
#[derive(Debug)]
struct PyFeature {
    feat: MolChemicalFeature,
    mol: Py<PyMol>,
    factory: Py<PyFactory>,
}

/// The features computed by the last `GetMolFeature(…, recompute=True)` call:
/// upstream keeps them in a function-level static, so `recompute=False`
/// indexes whatever the previous call produced.
static LAST_FEATURES: Mutex<Vec<Py<PyFeature>>> = Mutex::new(Vec::new());

fn features_of(
    py: Python<'_>,
    factory: &Bound<'_, PyFactory>,
    mol: &Bound<'_, PyMol>,
    include_only: &str,
    conf_id: i32,
) -> PyResult<Vec<Py<PyFeature>>> {
    let core = Arc::clone(&factory.borrow().factory);
    let found = {
        let m = mol.borrow();
        let mol_core = &m.mol;
        py.detach(|| core.features_for_mol(mol_core, include_only, conf_id))
            .py()?
    };
    found
        .into_iter()
        .map(|feat| {
            Py::new(
                py,
                PyFeature {
                    feat,
                    mol: mol.clone().unbind(),
                    factory: factory.clone().unbind(),
                },
            )
        })
        .collect()
}

#[pymethods]
impl PyFactory {
    /// Get the number of feature definitions
    #[pyo3(name = "GetNumFeatureDefs")]
    fn get_num_feature_defs(&self) -> usize {
        self.factory.num_feature_defs()
    }

    /// Get a tuple of feature types
    #[pyo3(name = "GetFeatureFamilies")]
    fn get_feature_families<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let mut families: Vec<&str> = Vec::new();
        for def in self.factory.feature_defs() {
            if !families.contains(&def.family()) {
                families.push(def.family());
            }
        }
        PyTuple::new(py, families)
    }

    /// Get a dictionary with SMARTS definitions for each feature type
    #[pyo3(name = "GetFeatureDefs")]
    fn get_feature_defs<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let res = PyDict::new(py);
        for def in self.factory.feature_defs() {
            res.set_item(
                format!("{}.{}", def.family(), def.feature_type()),
                def.smarts(),
            )?;
        }
        Ok(res)
    }

    /// Get the number of features the molecule has
    #[pyo3(name = "GetNumMolFeatures", signature = (mol, includeOnly = ""))]
    #[allow(non_snake_case)]
    fn get_num_mol_features(
        slf: &Bound<'_, Self>,
        mol: &Bound<'_, PyMol>,
        includeOnly: &str,
    ) -> PyResult<usize> {
        Ok(features_of(slf.py(), slf, mol, includeOnly, -1)?.len())
    }

    /// returns a particular feature (by index)
    #[pyo3(name = "GetMolFeature", signature = (mol, idx, includeOnly = "", recompute = true, confId = -1))]
    #[allow(non_snake_case)]
    fn get_mol_feature(
        slf: &Bound<'_, Self>,
        mol: &Bound<'_, PyMol>,
        idx: i64,
        includeOnly: &str,
        recompute: bool,
        confId: i32,
    ) -> PyResult<Py<PyFeature>> {
        let py = slf.py();
        let fresh = if recompute {
            Some(features_of(py, slf, mol, includeOnly, confId)?)
        } else {
            None
        };
        let mut last = LAST_FEATURES.lock().unwrap_or_else(|e| e.into_inner());
        if let Some(fresh) = fresh {
            *last = fresh;
        }
        usize::try_from(idx)
            .ok()
            .and_then(|i| last.get(i))
            .map(|f| f.clone_ref(py))
            .ok_or_else(|| PyIndexError::new_err(idx.to_string()))
    }
}

#[pymethods]
impl PyFeature {
    /// Returns the identifier of the feature
    #[pyo3(name = "GetId")]
    fn get_id(&self) -> i32 {
        self.feat.id()
    }

    /// Get the family to which the feature belongs; donor, acceptor, etc.
    #[pyo3(name = "GetFamily")]
    fn get_family(&self) -> &str {
        self.feat.family()
    }

    /// Get the specific type for the feature
    #[pyo3(name = "GetType")]
    fn get_type(&self) -> &str {
        self.feat.feature_type()
    }

    /// Get the location of the chemical feature (in the active conformer
    /// when no conformer id is given)
    #[pyo3(name = "GetPos", signature = (confId = None))]
    #[allow(non_snake_case)]
    fn get_pos(&mut self, py: Python<'_>, confId: Option<i32>) -> PyResult<PyPoint3D> {
        let mol = self.mol.bind(py).borrow();
        let p = match confId {
            Some(c) => self.feat.pos_in(&mol.mol, c),
            None => self.feat.pos(&mol.mol),
        }
        .py()?;
        Ok(PyPoint3D { p })
    }

    /// Get the IDs of the atoms that participate in the feature
    #[pyo3(name = "GetAtomIds")]
    fn get_atom_ids<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.feat.atoms())
    }

    /// Get the molecule used to derive the features
    #[pyo3(name = "GetMol")]
    fn get_mol(&self, py: Python<'_>) -> Py<PyMol> {
        self.mol.clone_ref(py)
    }

    /// Get the factory used to generate this feature
    #[pyo3(name = "GetFactory")]
    fn get_factory(&self, py: Python<'_>) -> Py<PyFactory> {
        self.factory.clone_ref(py)
    }

    /// Clears the cache used to store position information.
    #[pyo3(name = "ClearCache")]
    fn clear_cache(&mut self) {
        self.feat.clear_cache();
    }

    /// Sets the conformer to use (must be associated with a molecule).
    #[pyo3(name = "SetActiveConformer")]
    #[allow(non_snake_case)]
    fn set_active_conformer(&mut self, confId: i32) {
        self.feat.set_active_conformer(confId);
    }

    /// Gets the conformer to use.
    #[pyo3(name = "GetActiveConformer")]
    fn get_active_conformer(&self) -> i32 {
        self.feat.active_conformer()
    }
}

/// Construct a feature factory given a feature definition in a file
#[pyfunction]
#[pyo3(name = "BuildFeatureFactory", signature = (fileName))]
#[allow(non_snake_case)]
fn build_feature_factory_from_file(py: Python<'_>, fileName: &str) -> PyResult<PyFactory> {
    let text = std::fs::read(fileName)
        .map_err(|_| PyIOError::new_err(format!("File: {fileName} could not be opened.")))?;
    factory_from_text(py, &String::from_utf8_lossy(&text))
}

/// Construct a feature factory given a feature definition block
#[pyfunction]
#[pyo3(name = "BuildFeatureFactoryFromString", signature = (fdefString))]
#[allow(non_snake_case)]
fn build_feature_factory_from_string(py: Python<'_>, fdefString: &str) -> PyResult<PyFactory> {
    factory_from_text(py, fdefString)
}

fn factory_from_text(py: Python<'_>, text: &str) -> PyResult<PyFactory> {
    let factory = py
        .detach(|| build_feature_factory(text))
        .map_err(feature_error)?;
    Ok(PyFactory {
        factory: Arc::new(factory),
    })
}

/// Returns an empty list if any of the features passed in share an atom.
///  Otherwise a list of lists of atom indices is returned.
#[pyfunction]
#[pyo3(name = "GetAtomMatch", signature = (featMatch, maxAts = 1024))]
#[allow(non_snake_case)]
fn get_atom_match<'py>(
    py: Python<'py>,
    featMatch: &Bound<'py, PyAny>,
    maxAts: usize,
) -> PyResult<Bound<'py, pyo3::types::PyList>> {
    let mut used = vec![false; maxAts];
    let res = pyo3::types::PyList::empty(py);
    for i in 0..featMatch.len()? {
        let feat = featMatch.get_item(i)?;
        let feat = feat.cast::<PyFeature>()?.borrow();
        let local = pyo3::types::PyList::empty(py);
        for &idx in feat.feat.atoms() {
            // upstream's fixed-size bitset is indexed unchecked
            let slot = used
                .get_mut(idx as usize)
                .ok_or_else(|| PyIndexError::new_err(format!("Index Error: {idx}")))?;
            if *slot {
                return Ok(pyo3::types::PyList::empty(py));
            }
            *slot = true;
            local.append(idx)?;
        }
        res.append(local)?;
    }
    Ok(res)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing from chemical feature and functions to generate the",
    )?;
    m.add_class::<PyFactory>()?;
    m.add_class::<PyFeature>()?;
    m.add_function(wrap_pyfunction!(build_feature_factory_from_file, m)?)?;
    m.add_function(wrap_pyfunction!(build_feature_factory_from_string, m)?)?;
    m.add_function(wrap_pyfunction!(get_atom_match, m)?)?;
    Ok(())
}
