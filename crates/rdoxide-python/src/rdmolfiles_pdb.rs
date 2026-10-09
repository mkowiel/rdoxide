//! PDB reading and writing (upstream `Wrap/rdmolfiles.cpp` and
//! `Wrap/PDBWriter.cpp`): `MolFromPDBBlock`, `MolFromPDBFile`,
//! `MolToPDBBlock`, `MolToPDBFile` and `PDBWriter`.

use pyo3::prelude::*;
use rdoxide::RDKitError;
use rdoxide::graphmol::file_parsers::pdb_parser::{
    PdbParserParams, mol_from_pdb_data_stream, mol_from_pdb_file,
};
use rdoxide::graphmol::file_parsers::pdb_writer::{
    mol_to_pdb_block, mol_to_pdb_file, pdb_writer_record,
};

use crate::errors::IntoPyResult;
use crate::errors::{Check, Site};
use crate::mol::PyMol;
use crate::writer_sink::{WriterSink, close_sink, no_stream};

const WRITE: Site = Site {
    check: Check::Precondition,
    file: "Code/GraphMol/FileParsers/PDBWriter.cpp",
    line: 361,
    expr: "dp_ostream",
};
const FLUSH: Site = Site {
    check: Check::Precondition,
    file: "Code/GraphMol/FileParsers/MolWriters.h",
    line: 312,
    expr: "dp_ostream",
};

/// The upstream wrappers log `FileParseException`s as warnings and swallow
/// every other error, returning `None`.
fn finish(res: rdoxide::Result<Option<rdoxide::graphmol::RWMol>>) -> PyResult<Option<PyMol>> {
    match res {
        Ok(mol) => Ok(mol.map(|mol| PyMol { mol })),
        Err(RDKitError::BadFile(m)) => Err(pyo3::exceptions::PyIOError::new_err(m)),
        Err(RDKitError::FileParse(m)) => {
            rdoxide::rd_log_warning!("{m}");
            Ok(None)
        }
        Err(_) => Ok(None),
    }
}

fn params(
    sanitize: bool,
    remove_hs: bool,
    flavor: u32,
    proximity_bonding: bool,
) -> PdbParserParams {
    PdbParserParams {
        sanitize,
        remove_hs,
        flavor,
        proximity_bonding,
    }
}

/// Construct a molecule from a PDB file.
///
///   ARGUMENTS:
///
///     - pdbFileName: name of the file to read
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to true.
///
///     - removeHs: (optional) toggles removing hydrogens from the molecule.
///       This only make sense when sanitization is done.
///       Defaults to true.
///
///     - flavor: (optional)
///
///     - proximityBonding: (optional) toggles automatic proximity bonding
///
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromPDBFile", signature = (pdbFileName, sanitize = true, removeHs = true, flavor = 0, proximityBonding = true))]
#[allow(non_snake_case)]
fn py_mol_from_pdb_file(
    py: Python<'_>,
    pdbFileName: &str,
    sanitize: bool,
    removeHs: bool,
    flavor: u32,
    proximityBonding: bool,
) -> PyResult<Option<PyMol>> {
    let ps = params(sanitize, removeHs, flavor, proximityBonding);
    finish(py.detach(|| mol_from_pdb_file(pdbFileName, &ps)))
}

/// Construct a molecule from a PDB block.
///
///   ARGUMENTS:
///
///     - molBlock: string containing the PDB block
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to True.
///
///     - removeHs: (optional) toggles removing hydrogens from the molecule.
///       This only make sense when sanitization is done.
///       Defaults to true.
///
///     - flavor: (optional)
///
///     - proximityBonding: (optional) toggles automatic proximity bonding
///
///   RETURNS:
///
///     a Mol object, None on failure.
#[pyfunction]
#[pyo3(name = "MolFromPDBBlock", signature = (molBlock, sanitize = true, removeHs = true, flavor = 0, proximityBonding = true))]
#[allow(non_snake_case)]
fn py_mol_from_pdb_block(
    py: Python<'_>,
    molBlock: &Bound<'_, PyAny>,
    sanitize: bool,
    removeHs: bool,
    flavor: u32,
    proximityBonding: bool,
) -> PyResult<Option<PyMol>> {
    // pyObjectToString accepts str or bytes
    let text = match molBlock.extract::<String>() {
        Ok(s) => s,
        Err(_) => String::from_utf8_lossy(&molBlock.extract::<Vec<u8>>()?).into_owned(),
    };
    let ps = params(sanitize, removeHs, flavor, proximityBonding);
    finish(py.detach(|| mol_from_pdb_data_stream(&text, &ps)))
}

/// Returns a PDB block for a molecule
///   ARGUMENTS:
///
///     - mol: the molecule
///     - confId: (optional) selects which conformation to output (-1 = default)
///     - flavor: (optional)
///             - flavor & 1 : Write MODEL/ENDMDL lines around each record
///             - flavor & 2 : Don't write any CONECT records
///             - flavor & 4 : Write CONECT records in both directions
///             - flavor & 8 : Don't use multiple CONECTs to encode bond order
///             - flavor & 16 : Write MASTER record
///             - flavor & 32 : Write TER record
///
///   RETURNS:
///
///     a string
#[pyfunction]
#[pyo3(name = "MolToPDBBlock", signature = (mol, confId = -1, flavor = 0))]
#[allow(non_snake_case)]
fn py_mol_to_pdb_block(mol: PyRef<'_, PyMol>, confId: i32, flavor: u32) -> PyResult<String> {
    mol_to_pdb_block(&mol.mol, confId, flavor).py()
}

/// Writes a PDB file for a molecule
///   ARGUMENTS:
///
///     - mol: the molecule
///     - filename: name of the file to write
///     - confId: (optional) selects which conformation to output (-1 = default)
///     - flavor: (optional)
///             - flavor & 1 : Write MODEL/ENDMDL lines around each record
///             - flavor & 2 : Don't write any CONECT records
///             - flavor & 4 : Write CONECT records in both directions
///             - flavor & 8 : Don't use multiple CONECTs to encode bond order
///             - flavor & 16 : Write MASTER record
///             - flavor & 32 : Write TER record
#[pyfunction]
#[pyo3(name = "MolToPDBFile", signature = (mol, filename, confId = -1, flavor = 0))]
#[allow(non_snake_case)]
fn py_mol_to_pdb_file(
    mol: PyRef<'_, PyMol>,
    filename: &str,
    confId: i32,
    flavor: u32,
) -> PyResult<()> {
    mol_to_pdb_file(&mol.mol, filename, confId, flavor).py()
}

/// A class for writing molecules to PDB files.
#[pyclass(name = "PDBWriter", module = "rdkit.Chem.rdmolfiles", unsendable)]
pub struct PyPdbWriter {
    sink: Option<WriterSink>,
    flavor: u32,
    count: u32,
}

#[pymethods]
impl PyPdbWriter {
    /// Constructor.
    ///
    ///    ARGUMENTS:
    ///
    ///      - fileName: name of the output file. ('-' to write to stdout)
    ///      - flavor: (optional)
    #[new]
    #[pyo3(signature = (fileName, flavor = 0))]
    #[allow(non_snake_case)]
    fn new(fileName: &Bound<'_, PyAny>, flavor: u32) -> PyResult<Self> {
        Ok(PyPdbWriter {
            sink: Some(WriterSink::open(fileName)?),
            flavor,
            count: 0,
        })
    }

    fn __enter__(slf: Py<Self>) -> Py<Self> {
        slf
    }

    #[pyo3(signature = (exc_type, exc_val, traceback))]
    fn __exit__(
        &mut self,
        py: Python<'_>,
        exc_type: &Bound<'_, PyAny>,
        exc_val: &Bound<'_, PyAny>,
        traceback: &Bound<'_, PyAny>,
    ) -> PyResult<bool> {
        let _ = (exc_type, exc_val, traceback);
        self.close(py)?;
        Ok(false)
    }

    /// Writes a molecule to the output file.
    ///
    ///   ARGUMENTS:
    ///
    ///     - mol: the Mol to be written
    ///     - confId: (optional) ignored
    #[pyo3(signature = (mol, confId = -1))]
    #[allow(non_snake_case)]
    fn write(&mut self, py: Python<'_>, mol: PyRef<'_, PyMol>, confId: i32) -> PyResult<()> {
        let Some(sink) = self.sink.as_mut() else {
            return Err(no_stream(WRITE));
        };
        self.count += 1;
        let text = pdb_writer_record(&mol.mol, confId, self.flavor, self.count).py()?;
        sink.write(py, &text)
    }

    /// Flushes the output file (forces the disk file to be updated).
    fn flush(&mut self, py: Python<'_>) -> PyResult<()> {
        self.sink
            .as_mut()
            .ok_or_else(|| no_stream(FLUSH))?
            .flush(py)
    }

    /// Flushes the output file and closes it. The Writer cannot be used after this.
    fn close(&mut self, py: Python<'_>) -> PyResult<()> {
        close_sink(py, &mut self.sink)
    }

    /// Returns the number of molecules written so far.
    #[pyo3(name = "NumMols")]
    fn num_mols(&self) -> u32 {
        self.count
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_mol_from_pdb_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_from_pdb_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_pdb_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_pdb_file, m)?)?;
    m.add_class::<PyPdbWriter>()
}
