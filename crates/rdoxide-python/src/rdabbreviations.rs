//! `rdkit.Chem.rdAbbreviations` (upstream
//! `GraphMol/Abbreviations/Wrap/rdAbbreviations.cpp`).

use std::sync::Arc;

use pyo3::prelude::*;
use pyo3::types::PyList;
use rdoxide::graphmol::abbreviations::{
    AbbreviationDefinition, condense_abbreviation_substance_groups, condense_mol_abbreviations,
    default_abbreviations, default_linkers, label_mol_abbreviations, parse_abbreviations,
    parse_linkers,
};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Abbreviation Definition
#[pyclass(
    name = "AbbreviationDefinition",
    module = "rdkit.Chem.rdAbbreviations",
    skip_from_py_object
)]
#[derive(Clone, Default)]
pub struct PyAbbreviationDefinition {
    defn: AbbreviationDefinition,
}

#[pymethods]
impl PyAbbreviationDefinition {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// the label
    #[getter]
    fn label(&self) -> String {
        self.defn.label.clone()
    }
    #[setter(label)]
    fn set_label(&mut self, v: String) {
        self.defn.label = v;
    }

    /// the label in a drawing when the bond comes from the right
    #[getter(displayLabel)]
    fn display_label(&self) -> String {
        self.defn.display_label.clone()
    }
    #[setter(displayLabel)]
    fn set_display_label(&mut self, v: String) {
        self.defn.display_label = v;
    }

    /// the label in a drawing when the bond comes from the west
    #[getter(displayLabelW)]
    fn display_label_w(&self) -> String {
        self.defn.display_label_w.clone()
    }
    #[setter(displayLabelW)]
    fn set_display_label_w(&mut self, v: String) {
        self.defn.display_label_w = v;
    }

    /// the query molecule (should have a dummy as the first atom if includesXBonds is true)
    #[getter]
    fn mol(&self) -> Option<PyMol> {
        self.defn.mol.as_ref().map(|m| PyMol {
            mol: m.as_ref().clone(),
        })
    }
    #[setter(mol)]
    fn set_mol(&mut self, v: Option<PyRef<'_, PyMol>>) {
        self.defn.mol = v.map(|m| Arc::new(m.mol.clone()));
    }

    /// whether or not the abbreviation definition includes bonds to non-abbreviation atoms
    #[getter(includesXBonds)]
    fn includes_x_bonds(&self) -> bool {
        self.defn.includes_x_bonds
    }
    #[setter(includesXBonds)]
    fn set_includes_x_bonds(&mut self, v: bool) {
        self.defn.includes_x_bonds = v;
    }
}

fn to_py_list(py: Python<'_>, defs: Vec<AbbreviationDefinition>) -> PyResult<Bound<'_, PyList>> {
    PyList::new(
        py,
        defs.into_iter()
            .map(|defn| Py::new(py, PyAbbreviationDefinition { defn }))
            .collect::<PyResult<Vec<_>>>()?,
    )
}

fn definitions(abbrevs: &Bound<'_, PyAny>) -> PyResult<Vec<AbbreviationDefinition>> {
    abbrevs
        .try_iter()?
        .map(|item| {
            Ok(item?
                .extract::<PyRef<'_, PyAbbreviationDefinition>>()?
                .defn
                .clone())
        })
        .collect()
}

/// returns a list of the default abbreviation definitions
#[pyfunction]
#[pyo3(name = "GetDefaultAbbreviations")]
fn py_get_default_abbreviations(py: Python<'_>) -> PyResult<Bound<'_, PyList>> {
    to_py_list(py, default_abbreviations())
}

/// returns a list of the default linker definitions
#[pyfunction]
#[pyo3(name = "GetDefaultLinkers")]
fn py_get_default_linkers(py: Python<'_>) -> PyResult<Bound<'_, PyList>> {
    to_py_list(py, default_linkers())
}

/// Returns a set of abbreviation definitions from a string.  Format of the text data:  A series of lines, each of which contains: label SMARTS displayLabel displayLabelW  Where label is the label used for the abbreviation, SMARTS is the SMARTS definition of the abbreviation, displayLabel is used in drawings to render the abbreviations and displayLabelW is the display label if a bond comes in from the right.  The 'displayLabel' and 'displayLabelW' fields are optional.  Use dummies in the SMARTS to indicate attachment points. The assumption is that the first atom is a dummy (one will be added if this is not true) and that the second atom is the surrogate for the rest of the group.
#[pyfunction]
#[pyo3(name = "ParseAbbreviations", signature = (text, removeExtraDummies = false, allowConnectionToDummies = false))]
#[allow(non_snake_case)]
fn py_parse_abbreviations<'py>(
    py: Python<'py>,
    text: &str,
    removeExtraDummies: bool,
    allowConnectionToDummies: bool,
) -> PyResult<Bound<'py, PyList>> {
    to_py_list(
        py,
        parse_abbreviations(text, removeExtraDummies, allowConnectionToDummies).py()?,
    )
}

/// Returns a set of linker definitions from a string.  Equivalent to calling ParseAbbreviations(text, True True).
#[pyfunction]
#[pyo3(name = "ParseLinkers", signature = (text))]
fn py_parse_linkers<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyList>> {
    to_py_list(py, parse_linkers(text).py()?)
}

/// Finds and replaces abbreviations in a molecule. The result is not sanitized.
#[pyfunction]
#[pyo3(name = "CondenseMolAbbreviations", signature = (mol, abbrevs, maxCoverage = 0.4, sanitize = true))]
#[allow(non_snake_case)]
fn py_condense_mol_abbreviations(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    abbrevs: &Bound<'_, PyAny>,
    maxCoverage: f64,
    sanitize: bool,
) -> PyResult<PyMol> {
    let abbrevs = definitions(abbrevs)?;
    let mut res = mol.mol.clone();
    py.detach(|| condense_mol_abbreviations(&mut res, &abbrevs, maxCoverage, sanitize))
        .py()?;
    Ok(PyMol { mol: res })
}

/// Finds abbreviations and adds to them to a molecule as "SUP" SubstanceGroups
#[pyfunction]
#[pyo3(name = "LabelMolAbbreviations", signature = (mol, abbrevs, maxCoverage = 0.4))]
#[allow(non_snake_case)]
fn py_label_mol_abbreviations(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    abbrevs: &Bound<'_, PyAny>,
    maxCoverage: f64,
) -> PyResult<PyMol> {
    let abbrevs = definitions(abbrevs)?;
    let mut res = mol.mol.clone();
    py.detach(|| label_mol_abbreviations(&mut res, &abbrevs, maxCoverage))
        .py()?;
    Ok(PyMol { mol: res })
}

/// Finds and replaces abbreviation (i.e. "SUP") substance groups in a molecule. The result is not sanitized.
#[pyfunction]
#[pyo3(name = "CondenseAbbreviationSubstanceGroups", signature = (mol))]
fn py_condense_abbreviation_substance_groups(mol: PyRef<'_, PyMol>) -> PyResult<PyMol> {
    let mut res = mol.mol.clone();
    condense_abbreviation_substance_groups(&mut res).py()?;
    Ok(PyMol { mol: res })
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions for working with molecular abbreviations",
    )?;
    m.add_class::<PyAbbreviationDefinition>()?;
    m.add_function(wrap_pyfunction!(py_get_default_abbreviations, m)?)?;
    m.add_function(wrap_pyfunction!(py_get_default_linkers, m)?)?;
    m.add_function(wrap_pyfunction!(py_parse_abbreviations, m)?)?;
    m.add_function(wrap_pyfunction!(py_parse_linkers, m)?)?;
    m.add_function(wrap_pyfunction!(py_condense_mol_abbreviations, m)?)?;
    m.add_function(wrap_pyfunction!(py_label_mol_abbreviations, m)?)?;
    m.add_function(wrap_pyfunction!(
        py_condense_abbreviation_substance_groups,
        m
    )?)?;
    Ok(())
}
