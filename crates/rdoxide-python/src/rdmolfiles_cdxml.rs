//! CDXML reading (upstream `Wrap/rdmolfiles.cpp`): `MolsFromCDXML`,
//! `MolsFromCDXMLFile`, `MolToCDXMLBlock`, `HasChemDrawCDXSupport`,
//! `CDXMLFormat` and `CDXMLParserParams`.

use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::RDKitError;
use rdoxide::graphmol::RWMol;
use rdoxide::graphmol::file_parsers::cdxml::{
    CdxmlFormat, CdxmlParserParams, has_chem_draw_cdx_support, mol_to_cdxml_block, mols_from_cdxml,
    mols_from_cdxml_file,
};

use crate::enums::add_enum;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdmolfiles_sequence::text_arg;

/// Parameters controlling conversion of a CDXML document to molecules
#[pyclass(
    name = "CDXMLParserParams",
    module = "rdoxide.Chem.rdmolfiles",
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
pub struct PyCdxmlParserParams {
    params: CdxmlParserParams,
}

fn format_from_obj(obj: &Bound<'_, PyAny>) -> PyResult<CdxmlFormat> {
    Ok(match obj.getattr("real")?.extract::<i64>()? {
        0 => CdxmlFormat::Cdxml,
        1 => CdxmlFormat::Cdx,
        _ => CdxmlFormat::Auto,
    })
}

#[pymethods]
impl PyCdxmlParserParams {
    /// Construct a default CDXMLFormat, or one from explicit values
    #[new]
    #[pyo3(signature = (sanitize = None, removeHs = None, format = None, parseQueries = false, strictQueryParsing = false))]
    #[allow(non_snake_case)]
    fn new(
        sanitize: Option<bool>,
        removeHs: Option<bool>,
        format: Option<&Bound<'_, PyAny>>,
        parseQueries: bool,
        strictQueryParsing: bool,
    ) -> PyResult<Self> {
        let mut params = CdxmlParserParams::default();
        if let (Some(s), Some(r), Some(f)) = (sanitize, removeHs, format) {
            params = CdxmlParserParams {
                sanitize: s,
                remove_hs: r,
                format: format_from_obj(f)?,
                parse_queries: parseQueries,
                strict_query_parsing: strictQueryParsing,
            };
        }
        Ok(Self { params })
    }

    /// controls whether or not the molecule is sanitized before being returned
    #[getter]
    fn sanitize(&self) -> bool {
        self.params.sanitize
    }

    #[setter]
    fn set_sanitize(&mut self, v: bool) {
        self.params.sanitize = v;
    }

    /// controls whether or not Hs are removed before the molecule is returned
    #[getter(removeHs)]
    fn remove_hs(&self) -> bool {
        self.params.remove_hs
    }

    #[setter(removeHs)]
    fn set_remove_hs(&mut self, v: bool) {
        self.params.remove_hs = v;
    }

    /// ChemDraw format One of Auto, CDXML, CDX.  For data streams, Auto defaults to CDXML
    #[getter]
    fn format<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let name = match self.params.format {
            CdxmlFormat::Cdxml => "CDXML",
            CdxmlFormat::Cdx => "CDX",
            CdxmlFormat::Auto => "Auto",
        };
        py.import("rdoxide.Chem.rdmolfiles")?
            .getattr("CDXMLFormat")?
            .getattr(name)
    }

    #[setter]
    fn set_format(&mut self, v: &Bound<'_, PyAny>) -> PyResult<()> {
        self.params.format = format_from_obj(v)?;
        Ok(())
    }

    /// controls whether CDXML query constructs are preserved as query atoms and bonds
    #[getter(parseQueries)]
    fn parse_queries(&self) -> bool {
        self.params.parse_queries
    }

    #[setter(parseQueries)]
    fn set_parse_queries(&mut self, v: bool) {
        self.params.parse_queries = v;
    }

    /// controls whether unsupported CDXML query constructs raise instead of warning
    #[getter(strictQueryParsing)]
    fn strict_query_parsing(&self) -> bool {
        self.params.strict_query_parsing
    }

    #[setter(strictQueryParsing)]
    fn set_strict_query_parsing(&mut self, v: bool) {
        self.params.strict_query_parsing = v;
    }

    fn __setattr__(slf: &Bound<'_, Self>, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        crate::pyconv::set_declared_attribute(slf.as_any(), name, value)
    }
}

/// `(sanitize, removeHs)` or a `CDXMLParserParams` as the trailing arguments.
fn params_from_args(
    args: &Bound<'_, PyTuple>,
    default: CdxmlParserParams,
) -> PyResult<CdxmlParserParams> {
    let mut params = default;
    match args.len() {
        0 => {}
        1 if args
            .get_item(0)?
            .extract::<PyRef<'_, PyCdxmlParserParams>>()
            .is_ok() =>
        {
            params = args
                .get_item(0)?
                .extract::<PyRef<'_, PyCdxmlParserParams>>()?
                .params
                .clone();
        }
        _ => {
            params.sanitize = args.get_item(0)?.is_truthy()?;
            if args.len() > 1 {
                params.remove_hs = args.get_item(1)?.is_truthy()?;
            }
        }
    }
    Ok(params)
}

fn to_tuple<'py>(py: Python<'py>, mols: Vec<RWMol>) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(
        py,
        mols.into_iter()
            .map(|mol| Py::new(py, PyMol { mol }))
            .collect::<PyResult<Vec<_>>>()?,
    )
}

/// Construct a molecule from a cdxml file.
///
///   Note that the CDXML format is large and complex, the RDKit doesn't support
///   full functionality, just the base ones required for molecule and
///   reaction parsing.
///
///   ARGUMENTS:
///     - filename: the cdxml filename
///     - sanitize: if True, sanitize the molecules [default True]
///     - removeHs: if True, convert explicit Hs into implicit Hs. [default True]
///       (or a CDXMLParserParams object in place of both)
///
///   RETURNS:
///     a tuple of parsed Mol objects.
#[pyfunction]
#[pyo3(name = "MolsFromCDXMLFile", signature = (filename, *args))]
fn mols_from_cdxml_file_py<'py>(
    py: Python<'py>,
    filename: &str,
    args: &Bound<'py, PyTuple>,
) -> PyResult<Bound<'py, PyTuple>> {
    let params = params_from_args(args, CdxmlParserParams::default())?;
    match mols_from_cdxml_file(filename, &params) {
        Ok(mols) => to_tuple(py, mols),
        Err(RDKitError::BadFile(m)) => Err(PyIOError::new_err(m)),
        Err(RDKitError::FileParse(m)) => {
            rdoxide::rd_log_warning!("{m}");
            to_tuple(py, Vec::new())
        }
        Err(_) => to_tuple(py, Vec::new()),
    }
}

/// Construct a molecule from a cdxml string.
///
///   Note that the CDXML format is large and complex, the RDKit doesn't support
///   full functionality, just the base ones required for molecule and
///   reaction parsing.
///
///   ARGUMENTS:
///     - cdxml: the cdxml string
///     - sanitize: if True, sanitize the molecules [default True]
///     - removeHs: if True, convert explicit Hs into implicit Hs. [default True]
///       (or a CDXMLParserParams object in place of both)
///
///   RETURNS:
///     a tuple of parsed Mol objects.
#[pyfunction]
#[pyo3(name = "MolsFromCDXML", signature = (cdxml, *args))]
fn mols_from_cdxml_py<'py>(
    py: Python<'py>,
    cdxml: &Bound<'py, PyAny>,
    args: &Bound<'py, PyTuple>,
) -> PyResult<Bound<'py, PyTuple>> {
    let params = params_from_args(args, CdxmlParserParams::default())?;
    let text = text_arg(cdxml)?;
    to_tuple(py, mols_from_cdxml(&text, &params).py()?)
}

/// write a CDX or CDXML block from a molecule; needs the ChemDraw library,
/// so it always raises here (see HasChemDrawCDXSupport())
#[pyfunction]
#[pyo3(name = "MolToCDXMLBlock", signature = (mol, format = None))]
fn mol_to_cdxml_block_py(
    mol: &Bound<'_, PyMol>,
    format: Option<&Bound<'_, PyAny>>,
) -> PyResult<String> {
    let format = format.map_or(Ok(CdxmlFormat::Cdxml), format_from_obj)?;
    mol_to_cdxml_block(&mol.borrow().mol, format).py()
}

/// Returns whether ChemDraw CDX (binary) parsing is available (never, here).
#[pyfunction]
#[pyo3(name = "HasChemDrawCDXSupport")]
fn has_chem_draw_cdx_support_py() -> bool {
    has_chem_draw_cdx_support()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    add_enum(
        m,
        "CDXMLFormat",
        &[("CDXML", 0), ("CDX", 1), ("Auto", 2)],
        None,
    )?;
    m.add_class::<PyCdxmlParserParams>()?;
    m.add_function(wrap_pyfunction!(mols_from_cdxml_file_py, m)?)?;
    m.add_function(wrap_pyfunction!(mols_from_cdxml_py, m)?)?;
    m.add_function(wrap_pyfunction!(mol_to_cdxml_block_py, m)?)?;
    m.add_function(wrap_pyfunction!(has_chem_draw_cdx_support_py, m)?)?;
    Ok(())
}
