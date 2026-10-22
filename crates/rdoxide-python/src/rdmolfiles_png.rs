//! `rdoxide.Chem.rdmolfiles`: molecules in PNG metadata — port of the PNG
//! parts of upstream `Code/GraphMol/Wrap/rdmolfiles.cpp`.

use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyList, PyTuple};
use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::png;
use rdoxide::graphmol::pickle::get_default_pickle_properties;
use rdoxide::graphmol::smiles_parse::SmilesParserParams;
use rdoxide::graphmol::smiles_parse::cx_write::{RestoreBondDirOption, cx_fields};

use crate::errors::{IntoPyResult, to_py};
use crate::mol::PyMol;
use crate::pyconv::{set_declared_attribute, text_bytes};
use crate::rdmolfiles::PySmilesParserParams;
use crate::rdmolfiles_write::PySmilesWriteParams;

/// Parameters controlling metadata included in PNG images
#[pyclass(
    name = "PNGMetadataParams",
    module = "rdoxide.Chem.rdmolfiles",
    unsendable
)]
#[derive(Debug)]
#[allow(non_snake_case)]
pub struct PyPngMetadataParams {
    /// toggles inclusion of molecule pickle (default=True)
    #[pyo3(get, set)]
    includePkl: bool,
    /// toggles inclusion of molecule CXSMILES (default=True)
    #[pyo3(get, set)]
    includeSmiles: bool,
    /// toggles inclusion of molecule molblock (default=False)
    #[pyo3(get, set)]
    includeMol: bool,
    /// choose properties to be included in the pickle (default=rdkit.Chem.rdchem.PropertyPickleOptions.NoProps)
    #[pyo3(get, set)]
    propertyFlags: u32,
    /// choose SmilesWriteParams for the CXSMILES string (default=rdkit.Chem.rdmolfiles.SmilesWriteParams())
    #[pyo3(get)]
    smilesWriteParams: Py<PySmilesWriteParams>,
    /// choose CXSMILES fields to be included in the CXSMILES string (default=rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ALL)
    #[pyo3(get, set)]
    cxSmilesFlags: u32,
    /// choose what to do with bond dirs in the CXSMILES string (default=rdkit.Chem.rdmolfiles.RestoreBondDirOption.RestoreBondDirOptionClear)
    #[pyo3(get)]
    restoreBondDirs: Py<PyAny>,
}

fn restore_option_object(py: Python<'_>, value: i64) -> PyResult<Py<PyAny>> {
    Ok(py
        .import("rdoxide.Chem.rdmolfiles")?
        .getattr("RestoreBondDirOption")?
        .call1((value,))?
        .unbind())
}

#[pymethods]
#[allow(non_snake_case)]
impl PyPngMetadataParams {
    #[new]
    fn py_new(py: Python<'_>) -> PyResult<Self> {
        Ok(PyPngMetadataParams {
            includePkl: true,
            includeSmiles: true,
            includeMol: false,
            propertyFlags: get_default_pickle_properties(),
            smilesWriteParams: Py::new(py, PySmilesWriteParams::py_new())?,
            cxSmilesFlags: cx_fields::ALL,
            restoreBondDirs: restore_option_object(py, 1)?,
        })
    }

    #[setter]
    fn set_smilesWriteParams(
        &mut self,
        py: Python<'_>,
        value: PyRef<'_, PySmilesWriteParams>,
    ) -> PyResult<()> {
        self.smilesWriteParams = Py::new(py, value.clone())?;
        Ok(())
    }

    #[setter]
    fn set_restoreBondDirs(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        value.extract::<i64>()?;
        self.restoreBondDirs = value.clone().unbind();
        Ok(())
    }

    fn __setattr__(slf: &Bound<'_, Self>, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        set_declared_attribute(slf.as_any(), name, value)
    }
}

impl PyPngMetadataParams {
    fn to_core(&self, py: Python<'_>) -> PyResult<png::PngMetadataParams> {
        let restore = if self.restoreBondDirs.bind(py).extract::<i64>()? == 0 {
            RestoreBondDirOption::True
        } else {
            RestoreBondDirOption::Clear
        };
        Ok(png::PngMetadataParams {
            include_pkl: self.includePkl,
            include_smiles: self.includeSmiles,
            include_mol: self.includeMol,
            property_flags: self.propertyFlags,
            smiles_write_params: self.smilesWriteParams.borrow(py).to_core(),
            cx_smiles_flags: self.cxSmilesFlags,
            restore_bond_dirs: restore,
        })
    }
}

fn parser_params(params: Option<&Bound<'_, PyAny>>) -> PyResult<SmilesParserParams> {
    match params.filter(|p| !p.is_none()) {
        Some(p) => Ok(p.extract::<PyRef<'_, PySmilesParserParams>>()?.to_core()),
        None => Ok(SmilesParserParams::default()),
    }
}

/// `std::ifstream` + `BadFileException` for the readers
pub(crate) fn read_input(path: &str) -> PyResult<Vec<u8>> {
    std::fs::read(path).map_err(|_| to_py(RDKitError::BadFile(format!("Bad input file {path}"))))
}

/// The bytes of a file the writers read through a stream that is never
/// checked: a missing file is an empty PNG.
pub(crate) fn read_unchecked(path: &str) -> Vec<u8> {
    std::fs::read(path).unwrap_or_default()
}

fn py_bytes<'py>(py: Python<'py>, data: &[u8]) -> Bound<'py, PyBytes> {
    PyBytes::new(py, data)
}

/// Construct a molecule from metadata in a PNG string.
///
///      ARGUMENTS:
///
///        - png: the PNG string
///
///        - params: used to provide optional parameters for the metadata parsing
///
///      RETURNS:
///        a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromPNGString", signature = (png, params = None))]
fn mol_from_png_string(
    py: Python<'_>,
    png: &Bound<'_, PyAny>,
    params: Option<&Bound<'_, PyAny>>,
) -> PyResult<Option<PyMol>> {
    let params = parser_params(params)?;
    let data = text_bytes(png)?;
    finish(py.detach(|| png::png_to_mol(&data, &params)))
}

/// Construct a molecule from metadata in a PNG file.
///
///      ARGUMENTS:
///
///        - filename: the PNG filename
///
///        - params: used to provide optional parameters for the metadata parsing
///
///      RETURNS:
///        a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromPNGFile", signature = (filename, params = None))]
fn mol_from_png_file(
    py: Python<'_>,
    filename: &str,
    params: Option<&Bound<'_, PyAny>>,
) -> PyResult<Option<PyMol>> {
    let params = parser_params(params)?;
    let data = read_input(filename)?;
    finish(py.detach(|| png::png_to_mol(&data, &params)))
}

fn finish(res: rdoxide::Result<Option<rdoxide::graphmol::ROMol>>) -> PyResult<Option<PyMol>> {
    match res {
        Ok(mol) => Ok(mol.map(|mol| PyMol { mol })),
        Err(RDKitError::FileParse(m)) => {
            rdoxide::rd_log_warning!("{m}");
            Ok(None)
        }
        Err(_) => Ok(None),
    }
}

fn mols_tuple<'py>(
    py: Python<'py>,
    mols: Vec<rdoxide::graphmol::ROMol>,
) -> PyResult<Bound<'py, PyTuple>> {
    let items = mols
        .into_iter()
        .map(|mol| Py::new(py, PyMol { mol }))
        .collect::<PyResult<Vec<_>>>()?;
    PyTuple::new(py, items)
}

/// returns a tuple of molecules constructed from the PNG string
#[pyfunction]
#[pyo3(name = "MolsFromPNGString", signature = (png, tag = png::PKL_TAG, params = None))]
fn mols_from_png_string<'py>(
    py: Python<'py>,
    png: &Bound<'py, PyAny>,
    tag: &str,
    params: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let params = parser_params(params)?;
    let data = text_bytes(png)?;
    let mols = py.detach(|| png::png_to_mols(&data, tag, &params)).py()?;
    mols_tuple(py, mols)
}

/// returns a tuple of molecules constructed from the PNG file
#[pyfunction]
#[pyo3(name = "MolsFromPNGFile", signature = (filename, tag = png::PKL_TAG, params = None))]
fn mols_from_png_file<'py>(
    py: Python<'py>,
    filename: &str,
    tag: &str,
    params: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let params = parser_params(params)?;
    let data = read_input(filename)?;
    let mols = match py.detach(|| png::png_to_mols(&data, tag, &params)) {
        Ok(mols) => mols,
        Err(RDKitError::FileParse(m)) => {
            rdoxide::rd_log_warning!("{m}");
            Vec::new()
        }
        Err(_) => Vec::new(),
    };
    mols_tuple(py, mols)
}

fn metadata_from_dict(metadata: &Bound<'_, PyDict>) -> PyResult<png::PngMetadata> {
    metadata
        .iter()
        .map(|(k, v)| Ok((k.extract::<String>()?, v.extract::<String>()?.into_bytes())))
        .collect()
}

/// `addMolToPNG*` with either overload's arguments
fn mol_to_png(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    data: Vec<u8>,
    params: &png::PngMetadataParams,
) -> PyResult<Py<PyBytes>> {
    let mol = mol.borrow().mol.clone();
    let res = py
        .detach(|| png::add_mol_to_png(&mol, &data, params))
        .py()?;
    Ok(py_bytes(py, &res).unbind())
}

/// The third argument is `params` in one overload and `includePkl` in the other.
fn params_from_overloads(
    py: Python<'_>,
    third: Option<&Bound<'_, PyAny>>,
    include_smiles: bool,
    include_mol: bool,
) -> PyResult<png::PngMetadataParams> {
    match third {
        Some(p) if p.cast::<PyPngMetadataParams>().is_ok() => {
            p.cast::<PyPngMetadataParams>()?.borrow().to_core(py)
        }
        Some(v) => Ok(simple_params(v.is_truthy()?, include_smiles, include_mol)),
        None => Ok(simple_params(true, include_smiles, include_mol)),
    }
}

fn simple_params(
    include_pkl: bool,
    include_smiles: bool,
    include_mol: bool,
) -> png::PngMetadataParams {
    png::PngMetadataParams {
        include_pkl,
        include_smiles,
        include_mol,
        ..png::PngMetadataParams::default()
    }
}

/// Adds molecular metadata to PNG data read from a file.
///
///      ARGUMENTS:
///
///        - mol: the molecule
///
///        - filename: the PNG filename
///
///        - includePkl: include the RDKit's internal binary format in the output
///
///        - includeSmiles: include CXSmiles in the output
///
///        - includeMol: include CTAB (Mol block) in the output
///
///      RETURNS:
///        the updated PNG data
///
/// Alternatively: MolMetadataToPNGFile(mol, filename, params) with a PNGMetadataParams object.
#[pyfunction]
#[pyo3(name = "MolMetadataToPNGFile", signature = (mol, filename, includePkl = None, includeSmiles = true, includeMol = false))]
#[allow(non_snake_case)]
fn mol_metadata_to_png_file(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    filename: &Bound<'_, PyAny>,
    includePkl: Option<&Bound<'_, PyAny>>,
    includeSmiles: bool,
    includeMol: bool,
) -> PyResult<Py<PyBytes>> {
    let data = read_unchecked(&filename.extract::<String>()?);
    // the third argument is `params` in the other overload
    let params = params_from_overloads(py, includePkl, includeSmiles, includeMol)?;
    mol_to_png(py, mol, data, &params)
}

/// Adds molecular metadata to a PNG string.
///
///      ARGUMENTS:
///
///        - mol: the molecule
///
///        - png: the PNG string
///
///        - includePkl: include the RDKit's internal binary format in the output
///
///        - includeSmiles: include CXSmiles in the output
///
///        - includeMol: include CTAB (Mol block) in the output
///
///      RETURNS:
///        the updated PNG data
///
/// Alternatively: MolMetadataToPNGString(mol, png, params) with a PNGMetadataParams object.
#[pyfunction]
#[pyo3(name = "MolMetadataToPNGString", signature = (mol, png, includePkl = None, includeSmiles = true, includeMol = false))]
#[allow(non_snake_case)]
fn mol_metadata_to_png_string(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    png: &Bound<'_, PyAny>,
    includePkl: Option<&Bound<'_, PyAny>>,
    includeSmiles: bool,
    includeMol: bool,
) -> PyResult<Py<PyBytes>> {
    let data = png_bytes_arg(png)?;
    let params = params_from_overloads(py, includePkl, includeSmiles, includeMol)?;
    mol_to_png(py, mol, data, &params)
}

/// `python::extract<std::string>(png)`: bytes or str
fn png_bytes_arg(png: &Bound<'_, PyAny>) -> PyResult<Vec<u8>> {
    text_bytes(png)
}

/// Adds metadata to PNG data read from a file.
///
///      ARGUMENTS:
///
///        - metadata: a dict with the metadata to be added
///                    (keys and values should be strings)
///
///        - filename: the PNG filename
///
///      RETURNS:
///        the updated PNG data
#[pyfunction]
#[pyo3(name = "AddMetadataToPNGFile", signature = (metadata, filename))]
fn add_metadata_to_png_file(
    py: Python<'_>,
    metadata: &Bound<'_, PyDict>,
    filename: &Bound<'_, PyAny>,
) -> PyResult<Py<PyBytes>> {
    let metadata = metadata_from_dict(metadata)?;
    let data = read_unchecked(&filename.extract::<String>()?);
    let res = png::add_metadata_to_png(&data, &metadata, true).py()?;
    Ok(py_bytes(py, &res).unbind())
}

/// Adds metadata to a PNG string.
///
///      ARGUMENTS:
///
///        - metadata: a dict with the metadata to be added
///                    (keys and values should be strings)
///
///        - png: the PNG string
///
///      RETURNS:
///        the updated PNG data
#[pyfunction]
#[pyo3(name = "AddMetadataToPNGString", signature = (metadata, png))]
fn add_metadata_to_png_string(
    py: Python<'_>,
    metadata: &Bound<'_, PyDict>,
    png: &Bound<'_, PyAny>,
) -> PyResult<Py<PyBytes>> {
    let metadata = metadata_from_dict(metadata)?;
    let res = png::add_metadata_to_png(&png_bytes_arg(png)?, &metadata, true).py()?;
    Ok(py_bytes(py, &res).unbind())
}

fn metadata_object<'py>(
    py: Python<'py>,
    metadata: png::PngMetadata,
    as_list: bool,
) -> PyResult<Bound<'py, PyAny>> {
    if as_list {
        let list = PyList::empty(py);
        for (k, v) in metadata {
            list.append((k, py_bytes(py, &v)))?;
        }
        Ok(list.into_any())
    } else {
        let dict = PyDict::new(py);
        for (k, v) in metadata {
            dict.set_item(k, py_bytes(py, &v))?;
        }
        Ok(dict.into_any())
    }
}

/// Returns a dict with all metadata from the PNG file. Keys are strings, values are bytes. If asList is True, a list of (key, value) tuples is returned; this enables retrieving multiple values sharing the same key.
#[pyfunction]
#[pyo3(name = "MetadataFromPNGFile", signature = (filename, asList = false))]
#[allow(non_snake_case)]
fn metadata_from_png_file<'py>(
    py: Python<'py>,
    filename: &Bound<'py, PyAny>,
    asList: bool,
) -> PyResult<Bound<'py, PyAny>> {
    let path = filename.extract::<String>()?;
    let data = read_input(&path)?;
    metadata_object(py, png::png_to_metadata(&data).py()?, asList)
}

/// Returns a dict with all metadata from the PNG string. Keys are strings, values are bytes. If asList is True, a list of (key, value) tuples is returned; this enables retrieving multiple values sharing the same key.
#[pyfunction]
#[pyo3(name = "MetadataFromPNGString", signature = (png, asList = false))]
#[allow(non_snake_case)]
fn metadata_from_png_string<'py>(
    py: Python<'py>,
    png: &Bound<'py, PyAny>,
    asList: bool,
) -> PyResult<Bound<'py, PyAny>> {
    let data = png_bytes_arg(png)?;
    metadata_object(py, png::png_to_metadata(&data).py()?, asList)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPngMetadataParams>()?;
    m.add_function(wrap_pyfunction!(mol_from_png_string, m)?)?;
    m.add_function(wrap_pyfunction!(mol_from_png_file, m)?)?;
    m.add_function(wrap_pyfunction!(mols_from_png_string, m)?)?;
    m.add_function(wrap_pyfunction!(mols_from_png_file, m)?)?;
    m.add_function(wrap_pyfunction!(mol_metadata_to_png_file, m)?)?;
    m.add_function(wrap_pyfunction!(mol_metadata_to_png_string, m)?)?;
    m.add_function(wrap_pyfunction!(add_metadata_to_png_file, m)?)?;
    m.add_function(wrap_pyfunction!(add_metadata_to_png_string, m)?)?;
    m.add_function(wrap_pyfunction!(metadata_from_png_file, m)?)?;
    m.add_function(wrap_pyfunction!(metadata_from_png_string, m)?)?;
    Ok(())
}
