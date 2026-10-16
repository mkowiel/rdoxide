//! The pipeline part of `rdoxide.Chem.MolStandardize.rdMolStandardize`
//! (upstream `GraphMol/MolStandardize/nbWrap/Pipeline.cpp`).

use pyo3::prelude::*;
use rdoxide::graphmol::mol_standardize::pipeline::{
    Pipeline, PipelineLogEntry, PipelineOptions, PipelineResult, status,
};

use crate::enums::{add_enum, to_py};
use crate::errors::IntoPyResult;
use crate::vector_suite::{VectorSuite, no_python_init, vector_suite_methods};

const STATUS: [(&str, i64); 22] = [
    ("NO_EVENT", status::NO_EVENT as i64),
    ("INPUT_ERROR", status::INPUT_ERROR as i64),
    (
        "PREPARE_FOR_VALIDATION_ERROR",
        status::PREPARE_FOR_VALIDATION_ERROR as i64,
    ),
    (
        "FEATURES_VALIDATION_ERROR",
        status::FEATURES_VALIDATION_ERROR as i64,
    ),
    (
        "BASIC_VALIDATION_ERROR",
        status::BASIC_VALIDATION_ERROR as i64,
    ),
    (
        "IS2D_VALIDATION_ERROR",
        status::IS2D_VALIDATION_ERROR as i64,
    ),
    (
        "LAYOUT2D_VALIDATION_ERROR",
        status::LAYOUT2D_VALIDATION_ERROR as i64,
    ),
    (
        "STEREO_VALIDATION_ERROR",
        status::STEREO_VALIDATION_ERROR as i64,
    ),
    ("VALIDATION_ERROR", status::VALIDATION_ERROR as i64),
    (
        "PREPARE_FOR_STANDARDIZATION_ERROR",
        status::PREPARE_FOR_STANDARDIZATION_ERROR as i64,
    ),
    (
        "METAL_STANDARDIZATION_ERROR",
        status::METAL_STANDARDIZATION_ERROR as i64,
    ),
    (
        "NORMALIZER_STANDARDIZATION_ERROR",
        status::NORMALIZER_STANDARDIZATION_ERROR as i64,
    ),
    (
        "FRAGMENT_STANDARDIZATION_ERROR",
        status::FRAGMENT_STANDARDIZATION_ERROR as i64,
    ),
    (
        "CHARGE_STANDARDIZATION_ERROR",
        status::CHARGE_STANDARDIZATION_ERROR as i64,
    ),
    (
        "STANDARDIZATION_ERROR",
        status::STANDARDIZATION_ERROR as i64,
    ),
    ("OUTPUT_ERROR", status::OUTPUT_ERROR as i64),
    ("PIPELINE_ERROR", status::PIPELINE_ERROR as i64),
    ("METALS_DISCONNECTED", status::METALS_DISCONNECTED as i64),
    (
        "NORMALIZATION_APPLIED",
        status::NORMALIZATION_APPLIED as i64,
    ),
    ("FRAGMENTS_REMOVED", status::FRAGMENTS_REMOVED as i64),
    ("PROTONATION_CHANGED", status::PROTONATION_CHANGED as i64),
    (
        "STRUCTURE_MODIFICATION",
        status::STRUCTURE_MODIFICATION as i64,
    ),
];

const STAGE: [(&str, i64); 7] = [
    ("PARSING_INPUT", 1),
    ("PREPARE_FOR_VALIDATION", 2),
    ("VALIDATION", 3),
    ("PREPARE_FOR_STANDARDIZATION", 4),
    ("STANDARDIZATION", 5),
    ("SERIALIZING_OUTPUT", 9),
    ("COMPLETED", 10),
];

#[pyclass(
    name = "PipelineOptions",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
pub struct PyPipelineOptions {
    inner: PipelineOptions,
}

macro_rules! options {
    ($($py:literal, $get:ident, $set:ident, $field:ident: $ty:ty;)*) => {
        #[pymethods]
        impl PyPipelineOptions {
            #[new]
            fn new() -> Self {
                Self::default()
            }
            $(
                #[getter($py)]
                fn $get(&self) -> $ty {
                    self.inner.$field.clone()
                }
                #[setter($py)]
                fn $set(&mut self, value: $ty) {
                    self.inner.$field = value;
                }
            )*
        }
    };
}

options! {
    "strictParsing", get_strict_parsing, set_strict_parsing, strict_parsing: bool;
    "reportAllFailures", get_report_all_failures, set_report_all_failures, report_all_failures: bool;
    "allowEmptyMolecules", get_allow_empty_molecules, set_allow_empty_molecules, allow_empty_molecules: bool;
    "allowEnhancedStereo", get_allow_enhanced_stereo, set_allow_enhanced_stereo, allow_enhanced_stereo: bool;
    "allowAromaticBondType", get_allow_aromatic_bond_type, set_allow_aromatic_bond_type, allow_aromatic_bond_type: bool;
    "allowDativeBondType", get_allow_dative_bond_type, set_allow_dative_bond_type, allow_dative_bond_type: bool;
    "is2DZeroThreshold", get_is2d_zero_threshold, set_is2d_zero_threshold, is_2d_zero_threshold: f64;
    "atomClashLimit", get_atom_clash_limit, set_atom_clash_limit, atom_clash_limit: f64;
    "minMedianBondLength", get_min_median_bond_length, set_min_median_bond_length, min_median_bond_length: f64;
    "bondLengthLimit", get_bond_length_limit, set_bond_length_limit, bond_length_limit: f64;
    "allowLongBondsInRings", get_allow_long_bonds_in_rings, set_allow_long_bonds_in_rings, allow_long_bonds_in_rings: bool;
    "allowAtomBondClashExemption", get_allow_atom_bond_clash_exemption, set_allow_atom_bond_clash_exemption, allow_atom_bond_clash_exemption: bool;
    "metalNof", get_metal_nof, set_metal_nof, metal_nof: String;
    "metalNon", get_metal_non, set_metal_non, metal_non: String;
    "normalizerData", get_normalizer_data, set_normalizer_data, normalizer_data: String;
    "normalizerMaxRestarts", get_normalizer_max_restarts, set_normalizer_max_restarts, normalizer_max_restarts: u32;
    "scaledMedianBondLength", get_scaled_median_bond_length, set_scaled_median_bond_length, scaled_median_bond_length: f64;
    "outputV2000", get_output_v2000, set_output_v2000, output_v2000: bool;
}

#[pyclass(
    name = "PipelineLogEntry",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
pub struct PyPipelineLogEntry {
    inner: PipelineLogEntry,
}

#[pymethods]
impl PyPipelineLogEntry {
    #[getter]
    fn status(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(py, "PipelineStatus", self.inner.status as i64)
    }

    #[getter]
    fn detail(&self) -> String {
        self.inner.detail.clone()
    }

    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        other
            .cast::<PyPipelineLogEntry>()
            .is_ok_and(|o| o.borrow().inner == self.inner)
    }
}

/// The log entries of a pipeline run (upstream's `vector_indexing_suite`).
#[pyclass(
    name = "PipelineLog",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    sequence,
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
pub struct PyPipelineLog {
    entries: Vec<PipelineLogEntry>,
}

impl VectorSuite for PyPipelineLog {
    type Item = PipelineLogEntry;

    fn items(&self) -> &Vec<PipelineLogEntry> {
        &self.entries
    }

    fn items_mut(&mut self) -> &mut Vec<PipelineLogEntry> {
        &mut self.entries
    }

    fn to_py<'py>(py: Python<'py>, item: &PipelineLogEntry) -> PyResult<Bound<'py, PyAny>> {
        Ok(Bound::new(
            py,
            PyPipelineLogEntry {
                inner: item.clone(),
            },
        )?
        .into_any())
    }

    fn from_py(obj: &Bound<'_, PyAny>) -> Option<PipelineLogEntry> {
        obj.cast::<PyPipelineLogEntry>()
            .ok()
            .map(|e| e.borrow().inner.clone())
    }

    fn new_vector<'py>(
        py: Python<'py>,
        items: Vec<PipelineLogEntry>,
    ) -> PyResult<Bound<'py, PyAny>> {
        Ok(Bound::new(py, PyPipelineLog { entries: items })?.into_any())
    }
}

vector_suite_methods!(PyPipelineLog);

no_python_init!(PyPipelineLogEntry, PyPipelineLog, PyPipelineResult);

#[pyclass(
    name = "PipelineResult",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize"
)]
#[derive(Debug)]
pub struct PyPipelineResult {
    inner: PipelineResult,
    /// the log, shared with Python (upstream returns it by reference)
    log: Py<PyPipelineLog>,
}

#[pymethods]
impl PyPipelineResult {
    #[getter]
    fn status(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(py, "PipelineStatus", self.inner.status as i64)
    }

    #[getter]
    fn stage(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(py, "PipelineStage", self.inner.stage as i64)
    }

    #[getter]
    fn log(&self, py: Python<'_>) -> Py<PyPipelineLog> {
        self.log.clone_ref(py)
    }

    #[getter(inputMolData)]
    fn input_mol_data(&self) -> String {
        self.inner.input_mol_data.clone()
    }

    #[getter(outputMolData)]
    fn output_mol_data(&self) -> String {
        self.inner.output_mol_data.clone()
    }

    #[getter(parentMolData)]
    fn parent_mol_data(&self) -> String {
        self.inner.parent_mol_data.clone()
    }
}

#[pyclass(
    name = "Pipeline",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize"
)]
#[derive(Debug)]
pub struct PyPipeline {
    inner: Pipeline,
}

#[pymethods]
impl PyPipeline {
    #[new]
    #[pyo3(signature = (options = None))]
    fn new(options: Option<PyRef<'_, PyPipelineOptions>>) -> Self {
        PyPipeline {
            inner: Pipeline::new(options.map(|o| o.inner.clone()).unwrap_or_default()),
        }
    }

    #[pyo3(signature = (molData))]
    #[allow(non_snake_case)]
    fn run(&self, py: Python<'_>, molData: &str) -> PyResult<PyPipelineResult> {
        let mut inner = py.detach(|| self.inner.run(molData)).py()?;
        let entries = std::mem::take(&mut inner.log);
        let log = Py::new(py, PyPipelineLog { entries })?;
        Ok(PyPipelineResult { inner, log })
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPipelineOptions>()?;
    add_enum(m, "PipelineStatus", &STATUS, None)?;
    add_enum(m, "PipelineStage", &STAGE, None)?;
    m.add_class::<PyPipelineLogEntry>()?;
    m.add_class::<PyPipelineLog>()?;
    m.add_class::<PyPipelineResult>()?;
    m.add_class::<PyPipeline>()?;
    Ok(())
}
