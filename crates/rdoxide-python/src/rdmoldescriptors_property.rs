//! `rdoxide.Chem.rdMolDescriptors`: the property registry — `Properties`,
//! `PythonPropertyFunctor`, `PropertyRangeQuery` and `MakePropertyRangeQuery`
//! (upstream `Descriptors/Property.cpp` and the `nbWrap` wrapper).
//!
//! The registry is a process-wide list of functor objects. The built-in
//! descriptors are registered the first time it is read, exactly like upstream's
//! `registerDescriptors`; Python subclasses of `PythonPropertyFunctor` that
//! override `__call__` can be added with `Properties.RegisterProperty`.

use std::sync::{Mutex, MutexGuard, Once};

use pyo3::exceptions::{PyKeyError, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::graphmol::descriptors::property::{
    DataFunc, PropertyFunctor, PropertyRangeQuery, builtin_properties,
};
use rdoxide::rdgeneral::dict::RDProps;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdmoldescriptors::bool_arg;
use crate::vects::{PyVectD, PyVectStr};

/// A named, versioned molecule -> number function. Built-in descriptors are
/// instances of this class; Python subclasses override `__call__`.
#[pyclass(
    name = "PythonPropertyFunctor",
    module = "rdoxide.Chem.rdMolDescriptors",
    subclass,
    weakref
)]
#[derive(Debug)]
pub struct PyPropertyFunctor {
    functor: PropertyFunctor,
}

impl PyPropertyFunctor {
    fn named(name: &str, version: &str) -> Self {
        PyPropertyFunctor {
            functor: PropertyFunctor::new(name, version, None),
        }
    }
}

#[pymethods]
impl PyPropertyFunctor {
    // Python subclasses are created by `Sub()` and then initialised with
    // `PythonPropertyFunctor.__init__(self, name, version)`, so `__new__` takes
    // whatever the subclass constructor is given.
    #[new]
    #[pyo3(signature = (*args, **_kwargs))]
    fn new(args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> Self {
        match (args.len(), args.get_item(0), args.get_item(1)) {
            (2, Ok(name), Ok(version)) => {
                match (name.extract::<String>(), version.extract::<String>()) {
                    (Ok(name), Ok(version)) => Self::named(&name, &version),
                    _ => Self::named("", ""),
                }
            }
            _ => Self::named("", ""),
        }
    }

    #[pyo3(name = "__init__", signature = (name, version))]
    fn init(&mut self, name: &str, version: &str) {
        self.functor.name = name.to_string();
        self.functor.version = version.to_string();
    }

    /// Compute the property for the specified molecule
    #[pyo3(name = "__call__", signature = (mol))]
    fn call(&self, mol: &Bound<'_, PyMol>) -> PyResult<f64> {
        self.functor.compute(&mut mol.borrow_mut().mol).py()
    }

    /// Return the name of the property to calculate
    #[pyo3(name = "GetName")]
    fn get_name(&self) -> String {
        self.functor.name.clone()
    }

    /// Return the version of the calculated property
    #[pyo3(name = "GetVersion")]
    fn get_version(&self) -> String {
        self.functor.version.clone()
    }
}

static REGISTRY: Mutex<Vec<Py<PyPropertyFunctor>>> = Mutex::new(Vec::new());
static BUILTINS: Once = Once::new();

fn lock_registry() -> MutexGuard<'static, Vec<Py<PyPropertyFunctor>>> {
    REGISTRY
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

/// `Properties::registerProperty`: replaces the functor of the same name or
/// appends, and returns the index.
fn register(py: Python<'_>, prop: Py<PyPropertyFunctor>) -> usize {
    let name = prop.borrow(py).functor.name.clone();
    let mut registry = lock_registry();
    if let Some(i) = registry
        .iter()
        .position(|p| p.borrow(py).functor.name == name)
    {
        registry[i] = prop;
        return i;
    }
    registry.push(prop);
    registry.len() - 1
}

/// `registerDescriptors`: the built-ins are added once, on first use of the
/// registry (but not by `registerProperty` itself, as upstream).
fn register_builtins(py: Python<'_>) -> PyResult<()> {
    let mut result = Ok(());
    BUILTINS.call_once(|| {
        for functor in builtin_properties() {
            match Py::new(py, PyPropertyFunctor { functor }) {
                Ok(prop) => {
                    register(py, prop);
                }
                Err(e) => {
                    result = Err(e);
                    return;
                }
            }
        }
    });
    result
}

/// A snapshot of the registered functors.
fn registered(py: Python<'_>) -> PyResult<Vec<Py<PyPropertyFunctor>>> {
    register_builtins(py)?;
    Ok(lock_registry().iter().map(|p| p.clone_ref(py)).collect())
}

/// `Properties::getProperty`
fn get_property(py: Python<'_>, name: &str) -> PyResult<Py<PyPropertyFunctor>> {
    registered(py)?
        .into_iter()
        .find(|p| p.borrow(py).functor.name == name)
        .ok_or_else(|| PyKeyError::new_err(name.to_string()))
}

/// Property computation and registry system.  To compute all registered properties:
/// mol = Chem.MolFromSmiles('c1ccccc1')
/// properties = rdMolDescriptors.Properties()
/// for name, value in zip(properties.GetPropertyNames(), properties.ComputeProperties(mol)):
///   print(name, value)
///
/// To compute a subset
/// properties = rdMolDescriptors.Properties(['exactmw', 'lipinskiHBA'])
/// for name, value in zip(properties.GetPropertyNames(), properties.ComputeProperties(mol)):
///   print(name, value)
#[pyclass(name = "Properties", module = "rdoxide.Chem.rdMolDescriptors")]
#[derive(Debug)]
pub struct PyProperties {
    properties: Vec<Py<PyPropertyFunctor>>,
}

/// Calls `prop(mol)` through Python so that overridden `__call__`s are honoured.
fn compute(prop: &Py<PyPropertyFunctor>, mol: &Bound<'_, PyMol>) -> PyResult<f64> {
    prop.bind(mol.py())
        .call_method1("__call__", (mol,))?
        .extract()
}

#[pymethods]
impl PyProperties {
    #[new]
    #[pyo3(signature = (propNames = None))]
    #[allow(non_snake_case)]
    fn new(py: Python<'_>, propNames: Option<Vec<String>>) -> PyResult<Self> {
        let properties = match propNames {
            None => registered(py)?,
            Some(names) => {
                register_builtins(py)?;
                names
                    .iter()
                    .map(|name| get_property(py, name))
                    .collect::<PyResult<Vec<_>>>()?
            }
        };
        Ok(PyProperties { properties })
    }

    /// Return the property names computed by this instance
    #[pyo3(name = "GetPropertyNames")]
    fn get_property_names(&self, py: Python<'_>) -> PyVectStr {
        PyVectStr {
            v: self
                .properties
                .iter()
                .map(|p| p.borrow(py).functor.name.clone())
                .collect(),
        }
    }

    /// Return a list of computed properties, if annotateMol==True, annotate the molecule with the computed properties.
    #[pyo3(name = "ComputeProperties", signature = (mol, annotateMol = None))]
    #[allow(non_snake_case)]
    fn compute_properties(
        &self,
        py: Python<'_>,
        mol: &Bound<'_, PyMol>,
        annotateMol: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let annotate = bool_arg(annotateMol, false)?;
        let mut res = Vec::with_capacity(self.properties.len());
        for prop in &self.properties {
            res.push(compute(prop, mol)?);
            if annotate {
                // upstream evaluates the property a second time to annotate
                let value = compute(prop, mol)?;
                let name = prop.borrow(py).functor.name.clone();
                mol.borrow_mut().mol.set_prop(&name, value, false).py()?;
            }
        }
        PyVectD::new_py(py, res)
    }

    /// Annotate the molecule with the computed properties.  These properties will be available as SDData or from mol.GetProp(prop)
    #[pyo3(name = "AnnotateProperties", signature = (mol))]
    fn annotate_properties(&self, py: Python<'_>, mol: &Bound<'_, PyMol>) -> PyResult<()> {
        for prop in &self.properties {
            let value = compute(prop, mol)?;
            let name = prop.borrow(py).functor.name.clone();
            mol.borrow_mut().mol.set_prop(&name, value, false).py()?;
        }
        Ok(())
    }

    /// Return all available property names that can be computed
    #[staticmethod]
    #[pyo3(name = "GetAvailableProperties")]
    fn get_available_properties(py: Python<'_>) -> PyResult<PyVectStr> {
        let v = registered(py)?
            .iter()
            .map(|p| p.borrow(py).functor.name.clone())
            .collect();
        Ok(PyVectStr { v })
    }

    /// Return the named property if it exists
    #[staticmethod]
    #[pyo3(name = "GetProperty", signature = (propName))]
    #[allow(non_snake_case)]
    fn get_property(py: Python<'_>, propName: &str) -> PyResult<Py<PyPropertyFunctor>> {
        get_property(py, propName)
    }

    /// Register a new property object (not thread safe)
    #[staticmethod]
    #[pyo3(name = "RegisterProperty", signature = (propertyFunctor))]
    #[allow(non_snake_case)]
    fn register_property(py: Python<'_>, propertyFunctor: &Bound<'_, PyAny>) -> PyResult<usize> {
        let prop = propertyFunctor
            .cast::<PyPropertyFunctor>()
            .map_err(|_| PyTypeError::new_err("propertyFunctor must be a PythonPropertyFunctor"))?;
        Ok(register(py, prop.clone().unbind()))
    }
}

/// Property Range Query for a molecule.  Match(mol) -> true if in range
#[pyclass(name = "PropertyRangeQuery", module = "rdoxide.Chem.rdMolDescriptors")]
#[derive(Debug)]
pub struct PyPropertyRangeQuery {
    query: PropertyRangeQuery,
}

#[pymethods]
impl PyPropertyRangeQuery {
    #[pyo3(name = "Match", signature = (what))]
    fn matches(&self, what: &Bound<'_, PyMol>) -> PyResult<bool> {
        self.query.matches(&mut what.borrow_mut().mol).py()
    }
}

/// Generates a Range property for the specified property, between min and max
/// query = MakePropertyRangeQuery('exactmw', 0, 500)
/// query.Match( mol )
#[pyfunction]
#[pyo3(name = "MakePropertyRangeQuery", signature = (name, min, max))]
fn make_property_range_query(
    py: Python<'_>,
    name: &str,
    min: f64,
    max: f64,
) -> PyResult<PyPropertyRangeQuery> {
    let data_func: Option<DataFunc> = get_property(py, name)?.borrow(py).functor.data_func;
    Ok(PyPropertyRangeQuery {
        query: PropertyRangeQuery::new(min, max, data_func),
    })
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPropertyFunctor>()?;
    m.add_class::<PyProperties>()?;
    m.add_class::<PyPropertyRangeQuery>()?;
    m.add_function(wrap_pyfunction!(make_property_range_query, m)?)
}
