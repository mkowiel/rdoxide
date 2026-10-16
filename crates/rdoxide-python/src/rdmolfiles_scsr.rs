//! SCSR mol files (upstream `Wrap/rdmolfiles.cpp`): `MolFromSCSRParams`,
//! `SCSRTemplateNames`, `SCSRBaseHbondOptions`, `MolFromSCSRBlock` and
//! `MolFromSCSRFile`.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::MolFileParserParams;
use rdoxide::graphmol::file_parsers::scsr::{
    MolFromScsrParams, ScsrBaseHbondOptions, ScsrTemplateNames, mol_from_scsr_block,
    mol_from_scsr_file,
};

use crate::enums::{add_enum, to_py};
use crate::errors::to_py as error_to_py;
use crate::mol::PyMol;

const TEMPLATE_NAMES: [(&str, ScsrTemplateNames); 3] = [
    ("UseFirstName", ScsrTemplateNames::UseFirstName),
    ("UseSecondName", ScsrTemplateNames::UseSecondName),
    ("AsEntered", ScsrTemplateNames::AsEntered),
];

const HBOND_OPTIONS: [(&str, ScsrBaseHbondOptions); 4] = [
    ("Ignore", ScsrBaseHbondOptions::Ignore),
    ("UseSapAll", ScsrBaseHbondOptions::UseSapAll),
    ("UseSapOne", ScsrBaseHbondOptions::UseSapOne),
    ("Auto", ScsrBaseHbondOptions::Auto),
];

fn template_names_from(value: i64) -> PyResult<ScsrTemplateNames> {
    [
        ScsrTemplateNames::AsEntered,
        ScsrTemplateNames::UseFirstName,
        ScsrTemplateNames::UseSecondName,
    ]
    .into_iter()
    .find(|v| *v as i64 == value)
    .ok_or_else(|| PyValueError::new_err("invalid SCSRTemplateNames value"))
}

fn hbond_options_from(value: i64) -> PyResult<ScsrBaseHbondOptions> {
    HBOND_OPTIONS
        .iter()
        .map(|(_, v)| *v)
        .find(|v| *v as i64 == value)
        .ok_or_else(|| PyValueError::new_err("invalid SCSRBaseHbondOptions value"))
}

/// Parameters controlling conversion of an SCSRMol to a Mol
#[pyclass(
    name = "MolFromSCSRParams",
    module = "rdoxide.Chem.rdmolfiles",
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
pub struct PyMolFromScsrParams {
    params: MolFromScsrParams,
}

#[pymethods]
impl PyMolFromScsrParams {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// include leaving groups atoms if not substited at that position
    #[getter(includeLeavingGroups)]
    fn get_include_leaving_groups(&self) -> bool {
        self.params.include_leaving_groups
    }

    #[setter(includeLeavingGroups)]
    fn set_include_leaving_groups(&mut self, value: bool) {
        self.params.include_leaving_groups = value;
    }

    /// If True, the first template name in the Sgroup is used as the Sgroup label
    #[getter(scsrTemplateNames)]
    fn get_scsr_template_names(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(
            py,
            "SCSRTemplateNames",
            self.params.scsr_template_names as i64,
        )
    }

    #[setter(scsrTemplateNames)]
    fn set_scsr_template_names(&mut self, value: i64) -> PyResult<()> {
        self.params.scsr_template_names = template_names_from(value)?;
        Ok(())
    }

    /// One of Ignore, UseSapAll(default) , UseSapOne, Auto
    #[getter(scsrBaseHbondOptions)]
    fn get_scsr_base_hbond_options(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(
            py,
            "SCSRBaseHbondOptions",
            self.params.scsr_base_hbond_options as i64,
        )
    }

    #[setter(scsrBaseHbondOptions)]
    fn set_scsr_base_hbond_options(&mut self, value: i64) -> PyResult<()> {
        self.params.scsr_base_hbond_options = hbond_options_from(value)?;
        Ok(())
    }
}

/// Upstream logs `FileParseException`s as warnings and swallows everything
/// else; only `BadFileException` reaches Python (as `OSError`) from the file
/// reader.
fn finish(
    res: rdoxide::Result<Option<rdoxide::graphmol::RWMol>>,
    file: bool,
) -> PyResult<Option<PyMol>> {
    match res {
        Ok(Some(mol)) => Ok(Some(PyMol { mol })),
        Ok(None) => Ok(None),
        Err(e @ RDKitError::BadFile(_)) if file => Err(error_to_py(e)),
        Err(RDKitError::FileParse(m)) => {
            rdoxide::rd_log_warning!("{m}");
            Ok(None)
        }
        Err(_) => Ok(None),
    }
}

fn file_params(sanitize: bool, remove_hs: bool) -> MolFileParserParams {
    MolFileParserParams {
        sanitize,
        remove_hs,
        strict_parsing: false,
        ..MolFileParserParams::default()
    }
}

fn scsr_params(params: Option<PyRef<'_, PyMolFromScsrParams>>) -> MolFromScsrParams {
    params.map_or_else(MolFromScsrParams::default, |p| p.params.clone())
}

/// Construct a molecule from an SCSR Mol block.
///
///       ARGUMENTS:
///
///         - molBlock: string containing the SCSR Mol block
///
///         - sanitize: (optional) toggles sanitization of the molecule.
///           Defaults to True.
///
///         - removeHs: (optional) toggles removing hydrogens from the molecule.
///           This only make sense when sanitization is done.
///           Defaults to true.
///
///         - molFromSCSRParams : MolFromSCSRParams to control conversion
///
///     RETURNS :
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromSCSRBlock", signature = (molBlock, sanitize = true, removeHs = true, molFromSCSRParams = None))]
#[allow(non_snake_case)]
fn mol_from_scsr_block_py(
    py: Python<'_>,
    molBlock: &str,
    sanitize: bool,
    removeHs: bool,
    molFromSCSRParams: Option<PyRef<'_, PyMolFromScsrParams>>,
) -> PyResult<Option<PyMol>> {
    let params = scsr_params(molFromSCSRParams);
    let p = file_params(sanitize, removeHs);
    finish(
        py.detach(|| mol_from_scsr_block(molBlock, &p, &params)),
        false,
    )
}

/// Construct a molecule from an SCSR Mol file.
///
///       ARGUMENTS:
///
///         - fileName: name of the file to read
///
///         - sanitize: (optional) toggles sanitization of the molecule.
///           Defaults to True.
///
///         - removeHs: (optional) toggles removing hydrogens from the molecule.
///           This only make sense when sanitization is done.
///           Defaults to true.
///
///         - molFromSCSRParams : MolFromSCSRParams to control conversion
///
///     RETURNS :
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromSCSRFile", signature = (molFileName, sanitize = true, removeHs = true, molFromSCSRParams = None))]
#[allow(non_snake_case)]
fn mol_from_scsr_file_py(
    py: Python<'_>,
    molFileName: String,
    sanitize: bool,
    removeHs: bool,
    molFromSCSRParams: Option<PyRef<'_, PyMolFromScsrParams>>,
) -> PyResult<Option<PyMol>> {
    let params = scsr_params(molFromSCSRParams);
    let p = file_params(sanitize, removeHs);
    finish(
        py.detach(|| mol_from_scsr_file(&molFileName, &p, &params)),
        true,
    )
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let hbond: Vec<(&str, i64)> = HBOND_OPTIONS.iter().map(|(n, v)| (*n, *v as i64)).collect();
    add_enum(m, "SCSRBaseHbondOptions", &hbond, None)?;
    let names: Vec<(&str, i64)> = TEMPLATE_NAMES
        .iter()
        .map(|(n, v)| (*n, *v as i64))
        .collect();
    add_enum(m, "SCSRTemplateNames", &names, None)?;
    m.add_class::<PyMolFromScsrParams>()?;
    m.add_function(wrap_pyfunction!(mol_from_scsr_block_py, m)?)?;
    m.add_function(wrap_pyfunction!(mol_from_scsr_file_py, m)?)?;
    Ok(())
}
