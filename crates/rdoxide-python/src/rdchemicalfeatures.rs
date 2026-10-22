//! `rdoxide.Chem.rdChemicalFeatures` (upstream `ChemicalFeatures/Wrap/`):
//! chemical features that are not associated with a molecule.

use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyTuple};
use rdoxide::chemical_features::FreeChemicalFeature;

use crate::errors::IntoPyResult;
use crate::geometry::PyPoint3D;

/// Class to represent free chemical features.
/// These chemical features are not associated with a molecule, though they can be matched
/// to molecular features
#[pyclass(
    name = "FreeChemicalFeature",
    module = "rdoxide.Chem.rdChemicalFeatures",
    skip_from_py_object,
    subclass,
    dict
)]
#[derive(Debug, Clone)]
struct PyFreeChemicalFeature {
    feat: FreeChemicalFeature,
}

#[pymethods]
impl PyFreeChemicalFeature {
    /// Subclasses call `FreeChemicalFeature.__init__(self, ...)`; the feature is
    /// built by `__new__`, so there is nothing left to do.
    #[pyo3(signature = (*_args, **_kwargs))]
    fn __init__(&self, _args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) {}

    /// `FreeChemicalFeature(pickle)`, `FreeChemicalFeature()`,
    /// `FreeChemicalFeature(family, loc)` or
    /// `FreeChemicalFeature(family, type, loc, id=-1)`
    #[new]
    #[pyo3(signature = (*args, **kwargs))]
    fn new(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        // keyword form: family, type, loc, id (as in the boost signature)
        let merged;
        let args = match kwargs {
            Some(kw) if !kw.is_empty() => {
                let mut items: Vec<Bound<'_, PyAny>> = args.iter().collect();
                for name in ["family", "type", "loc", "id"].iter().skip(items.len()) {
                    match kw.get_item(name)? {
                        Some(v) => items.push(v),
                        None => break,
                    }
                }
                merged = PyTuple::new(args.py(), items)?;
                &merged
            }
            _ => args,
        };
        let point =
            |obj: Bound<'_, PyAny>| -> PyResult<_> { Ok(obj.extract::<PyRef<'_, PyPoint3D>>()?.p) };
        let feat = match args.len() {
            0 => FreeChemicalFeature::default(),
            1 => {
                let arg = args.get_item(0)?;
                let bytes: Vec<u8> = match arg.extract::<Vec<u8>>() {
                    Ok(b) if arg.is_instance_of::<PyBytes>() => b,
                    _ => arg.extract::<String>()?.into_bytes(),
                };
                FreeChemicalFeature::from_bytes(&bytes).py()?
            }
            2 => FreeChemicalFeature::new(
                &args.get_item(0)?.extract::<String>()?,
                "",
                point(args.get_item(1)?)?,
                -1,
            ),
            3 | 4 => {
                let id = if args.len() == 4 {
                    args.get_item(3)?.extract()?
                } else {
                    -1
                };
                FreeChemicalFeature::new(
                    &args.get_item(0)?.extract::<String>()?,
                    &args.get_item(1)?.extract::<String>()?,
                    point(args.get_item(2)?)?,
                    id,
                )
            }
            n => {
                return Err(pyo3::exceptions::PyTypeError::new_err(format!(
                    "FreeChemicalFeature() takes at most 4 arguments ({n} given)"
                )));
            }
        };
        Ok(PyFreeChemicalFeature { feat })
    }

    /// Set the id of the feature
    #[pyo3(name = "SetId")]
    fn set_id(&mut self, id: i32) {
        self.feat.id = id;
    }
    /// Set the family of the feature
    #[pyo3(name = "SetFamily")]
    fn set_family(&mut self, family: &str) {
        self.feat.family = family.to_string();
    }
    /// Set the sepcific type for the feature
    #[pyo3(name = "SetType")]
    fn set_type(&mut self, r#type: &str) {
        self.feat.feature_type = r#type.to_string();
    }
    /// Get the id of the feature
    #[pyo3(name = "GetId")]
    fn get_id(&self) -> i32 {
        self.feat.id
    }
    /// Get the family of the feature
    #[pyo3(name = "GetFamily")]
    fn get_family(&self) -> &str {
        &self.feat.family
    }
    /// Get the sepcific type for the feature
    #[pyo3(name = "GetType")]
    fn get_type(&self) -> &str {
        &self.feat.feature_type
    }
    /// Set the feature position
    #[pyo3(name = "SetPos")]
    fn set_pos(&mut self, loc: PyRef<'_, PyPoint3D>) {
        self.feat.position = loc.p;
    }
    /// Get the position of the feature
    #[pyo3(name = "GetPos")]
    fn get_pos(&self) -> PyPoint3D {
        PyPoint3D {
            p: self.feat.position,
        }
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let py = slf.py();
        let pickle = PyBytes::new(py, &slf.borrow().feat.to_bytes());
        PyTuple::new(
            py,
            [
                slf.get_type().into_any(),
                PyTuple::new(py, [pickle])?.into_any(),
            ],
        )
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing free chemical feature functionality\n     These are features that are not associated with molecules. They are \n     typically derived from pharmacophores and site-maps.\n",
    )?;
    m.add_class::<PyFreeChemicalFeature>()?;
    Ok(())
}
