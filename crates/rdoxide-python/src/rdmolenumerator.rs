//! `rdkit.Chem.rdMolEnumerator` (upstream `GraphMol/MolEnumerator/Wrap/`).

use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rdoxide::graphmol::mol_enumerator::{
    self, LinkNodeOp, MolEnumeratorOp, MolEnumeratorParams, PositionVariationOp, RepeatUnitOp,
};

use crate::enums::add_enum;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::molbundle::PyMolBundle;

const ENUMERATOR_TYPE: &str = "EnumeratorType";

fn operation(typ: i64) -> PyResult<Arc<dyn MolEnumeratorOp>> {
    Ok(match typ {
        0 => Arc::new(LinkNodeOp::new()),
        1 => Arc::new(PositionVariationOp::new()),
        2 => Arc::new(RepeatUnitOp::new()),
        _ => return Err(PyValueError::new_err("unrecognized operator type")),
    })
}

/// Molecular enumerator parameters
#[pyclass(
    name = "MolEnumeratorParams",
    module = "rdkit.Chem.rdMolEnumerator",
    skip_from_py_object
)]
#[derive(Clone)]
#[allow(non_snake_case)]
struct PyMolEnumeratorParams {
    params: MolEnumeratorParams,
}

#[pymethods]
#[allow(non_snake_case)]
impl PyMolEnumeratorParams {
    #[new]
    #[pyo3(signature = (typ = None))]
    fn new(typ: Option<i64>) -> PyResult<Self> {
        let mut params = MolEnumeratorParams::default();
        if let Some(typ) = typ {
            params.operation = Some(operation(typ)?);
        }
        Ok(PyMolEnumeratorParams { params })
    }

    /// sanitize molecules after enumeration
    #[getter]
    fn get_sanitize(&self) -> bool {
        self.params.sanitize
    }
    #[setter]
    fn set_sanitize(&mut self, v: bool) {
        self.params.sanitize = v;
    }

    /// maximum number of molecules to enumerate
    #[getter]
    fn get_maxToEnumerate(&self) -> usize {
        self.params.max_to_enumerate
    }
    #[setter]
    fn set_maxToEnumerate(&mut self, v: usize) {
        self.params.max_to_enumerate = v;
    }

    /// do random enumeration (not yet implemented
    #[getter]
    fn get_doRandom(&self) -> bool {
        self.params.do_random
    }
    #[setter]
    fn set_doRandom(&mut self, v: bool) {
        self.params.do_random = v;
    }

    /// seed for the random enumeration (not yet implemented
    #[getter]
    fn get_randomSeed(&self) -> i32 {
        self.params.random_seed
    }
    #[setter]
    fn set_randomSeed(&mut self, v: i32) {
        self.params.random_seed = v;
    }

    /// set the operator to be used for enumeration
    fn SetEnumerationOperator(&mut self, typ: i64) -> PyResult<()> {
        self.params.operation = Some(operation(typ)?);
        Ok(())
    }

    fn __setattr__(slf: &Bound<'_, Self>, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        crate::pyconv::set_declared_attribute(slf.as_any(), name, value)
    }
}

/// do an enumeration and return a MolBundle.
///   If maxPerOperation is >0 that will be used as the maximum number of molecules which
///     can be returned by any given operation.
/// Limitations:
///   - the current implementation does not support molecules which include both
///     SRUs and LINKNODEs
///   - Overlapping SRUs, i.e. where one monomer is contained within another, are
///     not supported
///
/// or, with the parameters of an enumeration instead of maxPerOperation, do an
/// enumeration for the supplied parameter type and return a MolBundle
#[pyfunction]
#[pyo3(name = "Enumerate", signature = (mol, maxPerOperation = None, *, enumParams = None))]
#[allow(non_snake_case)]
fn py_enumerate(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    maxPerOperation: Option<&Bound<'_, PyAny>>,
    enumParams: Option<PyRef<'_, PyMolEnumeratorParams>>,
) -> PyResult<PyMolBundle> {
    let mut params = enumParams.map(|p| p.params.clone());
    if params.is_none()
        && let Some(arg) = maxPerOperation
        && let Ok(p) = arg.extract::<PyRef<'_, PyMolEnumeratorParams>>()
    {
        params = Some(p.params.clone());
    }
    let source = &mol.mol;
    let bundle = match params {
        Some(p) => py.detach(|| mol_enumerator::enumerate(source, &[p])),
        None => {
            let max = match maxPerOperation {
                Some(v) => v.extract::<usize>()?,
                None => 0,
            };
            py.detach(|| mol_enumerator::enumerate_default(source, max))
        }
    }
    .py()?;
    Ok(PyMolBundle::wrap(bundle))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing classes and functions for enumerating molecules",
    )?;
    add_enum(
        m,
        ENUMERATOR_TYPE,
        &[("LinkNode", 0), ("PositionVariation", 1), ("RepeatUnit", 2)],
        None,
    )?;
    m.add_class::<PyMolEnumeratorParams>()?;
    m.add_function(wrap_pyfunction!(py_enumerate, m)?)?;
    Ok(())
}
