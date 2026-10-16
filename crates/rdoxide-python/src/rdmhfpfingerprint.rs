//! `rdoxide.Chem.rdMHFPFingerprint` (upstream
//! `GraphMol/Fingerprints/Wrap/MHFPWrapper.cpp`).
//!
//! The bulk methods encode the molecules on all available cores.

use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::graphmol::fingerprints::mhfp::{MhfpEncoder, ShinglingParams};
use rdoxide::graphmol::romol::ROMol;
use rdoxide::rdgeneral::threads::par_map;

use crate::datastructs::PyExplicitBitVect;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::vects::{PyVectJ, PyVectStr};

/// Number of threads the bulk methods use (all cores).
const ALL_CORES: i32 = 0;

/// Applies Boost.Python's overload semantics of the upstream wrapper: the
/// keyword defaults say `kekulize=False`, but they are only used for
/// arguments *before* the last one supplied; when neither `kekulize`,
/// `min_radius` nor `length` is given the C++ default (`true`) is used.
fn shingling(
    radius: u8,
    rings: bool,
    isomeric: bool,
    kekulize: Option<bool>,
    min_radius: Option<u8>,
    length: Option<u32>,
) -> (ShinglingParams, u32) {
    let kekulize_default = min_radius.is_none() && length.is_none();
    (
        ShinglingParams {
            radius: radius.into(),
            rings,
            isomeric,
            kekulize: kekulize.unwrap_or(kekulize_default),
            min_radius: min_radius.unwrap_or(1).into(),
        },
        length.unwrap_or(2048),
    )
}

/// Encodes every item of `items` in parallel, keeping the order; the first
/// error (in input order) is raised.
fn encode_bulk<T: Sync, R: Send>(
    py: Python<'_>,
    items: &[T],
    f: impl Fn(&T) -> rdoxide::Result<R> + Sync + Send,
) -> PyResult<Vec<R>> {
    py.detach(|| par_map(items, ALL_CORES, f))
        .into_iter()
        .map(|r| r.py())
        .collect()
}

fn vect_tuple(py: Python<'_>, rows: Vec<Vec<u32>>) -> PyResult<Bound<'_, PyTuple>> {
    let items = rows
        .into_iter()
        .map(|v| Py::new(py, PyVectJ { v }))
        .collect::<PyResult<Vec<_>>>()?;
    PyTuple::new(py, items)
}

fn bit_vect_tuple<'py>(
    py: Python<'py>,
    rows: Vec<rdoxide::datastructs::ExplicitBitVect>,
) -> PyResult<Bound<'py, PyTuple>> {
    let items = rows
        .into_iter()
        .map(|v| Py::new(py, PyExplicitBitVect::wrap(v)))
        .collect::<PyResult<Vec<_>>>()?;
    PyTuple::new(py, items)
}

fn mols_of(list: &Bound<'_, pyo3::types::PyList>) -> PyResult<Vec<ROMol>> {
    list.iter()
        .map(|m| Ok(m.extract::<PyRef<'_, PyMol>>()?.mol.clone()))
        .collect()
}

/// A MinHash encoder (MHFP and SECFP fingerprints).
#[pyclass(name = "MHFPEncoder", module = "rdoxide.Chem.rdMHFPFingerprint")]
pub struct PyMhfpEncoder {
    encoder: MhfpEncoder,
}

#[pymethods]
impl PyMhfpEncoder {
    #[new]
    #[pyo3(signature = (n_permutations = 2048, seed = 42))]
    fn new(n_permutations: usize, seed: u32) -> Self {
        PyMhfpEncoder {
            encoder: MhfpEncoder::new(n_permutations, seed),
        }
    }

    /// Creates a MHFP vector from a list of arbitrary strings.
    #[pyo3(name = "FromStringArray", signature = (vec))]
    fn min_hash_strings(&self, vec: Vec<String>) -> PyVectJ {
        PyVectJ {
            v: self.encoder.from_string_array(&vec),
        }
    }

    /// Creates a MHFP vector from a list of unsigned integers.
    #[pyo3(name = "FromArray", signature = (vec))]
    fn min_hash_ints(&self, vec: Vec<u32>) -> PyVectJ {
        PyVectJ {
            v: self.encoder.from_array(&vec),
        }
    }

    /// Creates a shingling (a list of circular n-grams / substructures) from a SMILES string.
    #[pyo3(name = "CreateShinglingFromSmiles", signature = (smiles, radius = 3, rings = true, isomeric = false, kekulize = None, min_radius = None))]
    fn create_shingling_from_smiles(
        &self,
        smiles: &str,
        radius: u8,
        rings: bool,
        isomeric: bool,
        kekulize: Option<bool>,
        min_radius: Option<u8>,
    ) -> PyResult<PyVectStr> {
        let (p, _) = shingling(radius, rings, isomeric, kekulize, min_radius, None);
        Ok(PyVectStr {
            v: self.encoder.create_shingling_from_smiles(smiles, &p).py()?,
        })
    }

    /// Creates a shingling (a list of circular n-grams / substructures) from a RDKit Mol instance.
    #[pyo3(name = "CreateShinglingFromMol", signature = (mol, radius = 3, rings = true, isomeric = false, kekulize = None, min_radius = None))]
    fn create_shingling_from_mol(
        &self,
        mol: PyRef<'_, PyMol>,
        radius: u8,
        rings: bool,
        isomeric: bool,
        kekulize: Option<bool>,
        min_radius: Option<u8>,
    ) -> PyResult<PyVectStr> {
        let (p, _) = shingling(radius, rings, isomeric, kekulize, min_radius, None);
        Ok(PyVectStr {
            v: self.encoder.create_shingling(&mol.mol, &p).py()?,
        })
    }

    /// Creates a MHFP vector from a SMILES string.
    #[pyo3(name = "EncodeSmiles", signature = (smiles, radius = 3, rings = true, isomeric = false, kekulize = None, min_radius = None))]
    fn encode_smiles(
        &self,
        smiles: &str,
        radius: u8,
        rings: bool,
        isomeric: bool,
        kekulize: Option<bool>,
        min_radius: Option<u8>,
    ) -> PyResult<PyVectJ> {
        let (p, _) = shingling(radius, rings, isomeric, kekulize, min_radius, None);
        Ok(PyVectJ {
            v: self.encoder.encode_smiles(smiles, &p).py()?,
        })
    }

    /// Creates a MHFP vector from an RDKit Mol instance.
    #[pyo3(name = "EncodeMol", signature = (mol, radius = 3, rings = true, isomeric = false, kekulize = None, min_radius = None))]
    fn encode_mol(
        &self,
        mol: PyRef<'_, PyMol>,
        radius: u8,
        rings: bool,
        isomeric: bool,
        kekulize: Option<bool>,
        min_radius: Option<u8>,
    ) -> PyResult<PyVectJ> {
        let (p, _) = shingling(radius, rings, isomeric, kekulize, min_radius, None);
        Ok(PyVectJ {
            v: self.encoder.encode(&mol.mol, &p).py()?,
        })
    }

    /// Creates a MHFP vector from a list of SMILES strings.
    #[pyo3(name = "EncodeSmilesBulk", signature = (smiles, radius = 3, rings = true, isomeric = false, kekulize = None, min_radius = None))]
    fn encode_smiles_bulk<'py>(
        &self,
        py: Python<'py>,
        smiles: Vec<String>,
        radius: u8,
        rings: bool,
        isomeric: bool,
        kekulize: Option<bool>,
        min_radius: Option<u8>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let (p, _) = shingling(radius, rings, isomeric, kekulize, min_radius, None);
        let rows = encode_bulk(py, &smiles, |s| self.encoder.encode_smiles(s, &p))?;
        vect_tuple(py, rows)
    }

    /// Creates a MHFP vector from a list of RDKit Mol instances.
    #[pyo3(name = "EncodeMolsBulk", signature = (mols, radius = 3, rings = true, isomeric = false, kekulize = None, min_radius = None))]
    fn encode_mols_bulk<'py>(
        &self,
        py: Python<'py>,
        mols: &Bound<'py, pyo3::types::PyList>,
        radius: u8,
        rings: bool,
        isomeric: bool,
        kekulize: Option<bool>,
        min_radius: Option<u8>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let (p, _) = shingling(radius, rings, isomeric, kekulize, min_radius, None);
        let mols = mols_of(mols)?;
        let rows = encode_bulk(py, &mols, |m| self.encoder.encode(m, &p))?;
        vect_tuple(py, rows)
    }

    /// Creates a SECFP binary vector from a SMILES string.
    #[pyo3(name = "EncodeSECFPSmiles", signature = (smiles, radius = 3, rings = true, isomeric = false, kekulize = None, min_radius = None, length = None))]
    fn encode_secfp_smiles(
        &self,
        smiles: &str,
        radius: u8,
        rings: bool,
        isomeric: bool,
        kekulize: Option<bool>,
        min_radius: Option<u8>,
        length: Option<u32>,
    ) -> PyResult<PyExplicitBitVect> {
        let (p, length) = shingling(radius, rings, isomeric, kekulize, min_radius, length);
        Ok(PyExplicitBitVect::wrap(
            self.encoder.encode_secfp_smiles(smiles, &p, length).py()?,
        ))
    }

    /// Creates a SECFP binary vector from an RDKit Mol instance.
    #[pyo3(name = "EncodeSECFPMol", signature = (mol, radius = 3, rings = true, isomeric = false, kekulize = None, min_radius = None, length = None))]
    fn encode_secfp_mol(
        &self,
        mol: PyRef<'_, PyMol>,
        radius: u8,
        rings: bool,
        isomeric: bool,
        kekulize: Option<bool>,
        min_radius: Option<u8>,
        length: Option<u32>,
    ) -> PyResult<PyExplicitBitVect> {
        let (p, length) = shingling(radius, rings, isomeric, kekulize, min_radius, length);
        Ok(PyExplicitBitVect::wrap(
            self.encoder.encode_secfp(&mol.mol, &p, length).py()?,
        ))
    }

    /// Creates a SECFP binary vector from a list of SMILES strings.
    #[pyo3(name = "EncodeSECFPSmilesBulk", signature = (smiles, radius = 3, rings = true, isomeric = false, kekulize = None, min_radius = None, length = None))]
    fn encode_secfp_smiles_bulk<'py>(
        &self,
        py: Python<'py>,
        smiles: Vec<String>,
        radius: u8,
        rings: bool,
        isomeric: bool,
        kekulize: Option<bool>,
        min_radius: Option<u8>,
        length: Option<u32>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let (p, length) = shingling(radius, rings, isomeric, kekulize, min_radius, length);
        let rows = encode_bulk(py, &smiles, |s| {
            self.encoder.encode_secfp_smiles(s, &p, length)
        })?;
        bit_vect_tuple(py, rows)
    }

    /// Creates a SECFP binary vector from a list of RDKit Mol instances.
    #[pyo3(name = "EncodeSECFPMolsBulk", signature = (mols, radius = 3, rings = true, isomeric = false, kekulize = None, min_radius = None, length = None))]
    fn encode_secfp_mols_bulk<'py>(
        &self,
        py: Python<'py>,
        mols: &Bound<'py, pyo3::types::PyList>,
        radius: u8,
        rings: bool,
        isomeric: bool,
        kekulize: Option<bool>,
        min_radius: Option<u8>,
        length: Option<u32>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let (p, length) = shingling(radius, rings, isomeric, kekulize, min_radius, length);
        let mols = mols_of(mols)?;
        let rows = encode_bulk(py, &mols, |m| self.encoder.encode_secfp(m, &p, length))?;
        bit_vect_tuple(py, rows)
    }

    #[staticmethod]
    #[pyo3(name = "Distance", signature = (a, b))]
    fn distance(a: Vec<u32>, b: Vec<u32>) -> f64 {
        MhfpEncoder::distance(&a, &b)
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMhfpEncoder>()
}
