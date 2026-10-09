//! Legacy path fingerprints of `rdkit.Chem.rdmolops` (`RDKFingerprint`,
//! `UnfoldedRDKFingerprintCountBased`, `PatternFingerprint`, `LayeredFingerprint`; upstream
//! `GraphMol/Wrap/MolOps.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use rdoxide::datastructs::BitVect;
use rdoxide::error::RDKitError;
use rdoxide::graphmol::fingerprints::layered::{
    LAYERED_FINGERPRINT_MOL_VERSION, LayeredFingerprintArgs, SUBSTRUCT_LAYERS,
    layered_fingerprint_mol,
};
use rdoxide::graphmol::fingerprints::legacy::{
    PathFpOutput, rdk_fingerprint_mol, unfolded_rdk_fingerprint_mol,
};
use rdoxide::graphmol::fingerprints::pattern::{
    PATTERN_FINGERPRINT_MOL_VERSION, pattern_fingerprint_mol,
};

use crate::datastructs::PyExplicitBitVect;
use crate::datastructs_siv::PyULongSIV;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::pyconv::python_object_to_vect;

fn requested(o: Option<&Bound<'_, PyAny>>) -> bool {
    o.is_some_and(|o| !o.is_none())
}

/// Copies the extra outputs into the caller's `atomBits` list / `bitInfo` dict.
fn fill_outputs(
    py: Python<'_>,
    out: &PathFpOutput,
    n_atoms: usize,
    atom_bits: Option<&Bound<'_, PyAny>>,
    bit_info: Option<&Bound<'_, PyAny>>,
) -> PyResult<()> {
    if let (Some(ab), Some(l)) = (&out.atom_bits, atom_bits) {
        let l = l.cast::<PyList>()?;
        for i in 0..n_atoms {
            let vals = ab.get(i).cloned().unwrap_or_default();
            l.append(PyList::new(py, vals)?)?;
        }
    }
    if let (Some(bi), Some(d)) = (&out.bit_info, bit_info) {
        let d = d.cast::<PyDict>()?;
        for (k, paths) in bi {
            if d.contains(k)? {
                continue;
            }
            let lst = PyList::empty(py);
            for p in paths {
                lst.append(PyList::new(py, p)?)?;
            }
            d.set_item(k, lst)?;
        }
    }
    Ok(())
}

/// Returns an RDKit topological fingerprint for a molecule
#[pyfunction]
#[pyo3(name = "RDKFingerprint", signature = (mol, minPath = 1, maxPath = 7, fpSize = 2048, nBitsPerHash = 2, useHs = true, tgtDensity = 0.0, minSize = 128, branchedPaths = true, useBondOrder = true, atomInvariants = None, fromAtoms = None, atomBits = None, bitInfo = None))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn rdk_fingerprint(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    minPath: u32,
    maxPath: u32,
    fpSize: u32,
    nBitsPerHash: u32,
    useHs: bool,
    tgtDensity: f64,
    minSize: u32,
    branchedPaths: bool,
    useBondOrder: bool,
    atomInvariants: Option<&Bound<'_, PyAny>>,
    fromAtoms: Option<&Bound<'_, PyAny>>,
    atomBits: Option<&Bound<'_, PyAny>>,
    bitInfo: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyExplicitBitVect> {
    let n = mol.mol.num_atoms();
    let invs = python_object_to_vect(atomInvariants, None)?;
    let from = python_object_to_vect(fromAtoms, Some(n))?;
    let mut out = PathFpOutput {
        atom_bits: requested(atomBits).then(Vec::new),
        bit_info: requested(bitInfo).then(Default::default),
    };
    let bv = rdk_fingerprint_mol(
        &mol.mol,
        minPath,
        maxPath,
        fpSize,
        nBitsPerHash,
        useHs,
        tgtDensity,
        minSize,
        branchedPaths,
        useBondOrder,
        invs.as_deref(),
        from.as_deref(),
        &mut out,
    )
    .py()?;
    fill_outputs(py, &out, n as usize, atomBits, bitInfo)?;
    Ok(PyExplicitBitVect::wrap(bv))
}

/// Returns an unfolded count-based version of the RDKit fingerprint for a molecule
#[pyfunction]
#[pyo3(name = "UnfoldedRDKFingerprintCountBased", signature = (mol, minPath = 1, maxPath = 7, useHs = true, branchedPaths = true, useBondOrder = true, atomInvariants = None, fromAtoms = None, atomBits = None, bitInfo = None))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn unfolded_rdk_fingerprint(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    minPath: u32,
    maxPath: u32,
    useHs: bool,
    branchedPaths: bool,
    useBondOrder: bool,
    atomInvariants: Option<&Bound<'_, PyAny>>,
    fromAtoms: Option<&Bound<'_, PyAny>>,
    atomBits: Option<&Bound<'_, PyAny>>,
    bitInfo: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyULongSIV> {
    let n = mol.mol.num_atoms();
    let invs = python_object_to_vect(atomInvariants, None)?;
    let from = python_object_to_vect(fromAtoms, Some(n))?;
    let mut out = PathFpOutput {
        atom_bits: requested(atomBits).then(Vec::new),
        bit_info: requested(bitInfo).then(Default::default),
    };
    let fp = unfolded_rdk_fingerprint_mol(
        &mol.mol,
        minPath,
        maxPath,
        useHs,
        branchedPaths,
        useBondOrder,
        invs.as_deref(),
        from.as_deref(),
        &mut out,
    )
    .py()?;
    fill_outputs(py, &out, n as usize, atomBits, bitInfo)?;
    Ok(PyULongSIV::wrap(fp))
}

/// `atomCounts` argument of the screening fingerprints: `None` when absent or
/// empty, otherwise the (unsigned) values, which must cover all atoms.
fn read_atom_counts(counts: Option<&Bound<'_, PyAny>>, n_atoms: u32) -> PyResult<Option<Vec<u32>>> {
    let Some(counts) = counts.filter(|c| c.is_truthy().unwrap_or(false)) else {
        return Ok(None);
    };
    if (counts.len()? as u32) < n_atoms {
        return Err(PyValueError::new_err(
            "atomCounts shorter than the number of atoms",
        ));
    }
    counts
        .try_iter()?
        .map(|item| item?.extract::<u32>())
        .collect::<PyResult<Vec<_>>>()
        .map(Some)
}

/// Returns a layered fingerprint for a molecule
///
/// NOTE: This function is experimental. The API or results may change from
///   release to release.
///
/// Layer definitions:
///    - 0x01: pure topology
///    - 0x02: bond order
///    - 0x04: atom types
///    - 0x08: presence of rings
///    - 0x10: ring sizes
///    - 0x20: aromaticity
#[pyfunction]
#[pyo3(name = "LayeredFingerprint", signature = (mol, layerFlags = 0xFFFF_FFFF, minPath = 1, maxPath = 7, fpSize = 2048, atomCounts = None, setOnlyBits = None, branchedPaths = true, fromAtoms = None))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn layered_fingerprint(
    mol: PyRef<'_, PyMol>,
    layerFlags: u32,
    minPath: u32,
    maxPath: u32,
    fpSize: u32,
    atomCounts: Option<&Bound<'_, PyAny>>,
    setOnlyBits: Option<PyRef<'_, PyExplicitBitVect>>,
    branchedPaths: bool,
    fromAtoms: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyExplicitBitVect> {
    let from = python_object_to_vect(fromAtoms, Some(mol.mol.num_atoms()))?;
    let mut counts = read_atom_counts(atomCounts, mol.mol.num_atoms())?;
    let args = LayeredFingerprintArgs {
        layer_flags: layerFlags,
        min_path: minPath,
        max_path: maxPath,
        fp_size: fpSize,
        branched_paths: branchedPaths,
        from_atoms: from.as_deref(),
    };
    let bv = layered_fingerprint_mol(
        &mol.mol,
        &args,
        counts.as_deref_mut(),
        setOnlyBits.as_deref().map(|b| &b.bv),
    )
    .py()?;
    if let (Some(counts), Some(list)) = (counts, atomCounts) {
        for (i, c) in counts.into_iter().enumerate() {
            list.set_item(i, c)?;
        }
    }
    Ok(PyExplicitBitVect::wrap(bv))
}

/// A fingerprint using SMARTS patterns
///
/// NOTE: This function is experimental. The API or results may change from
///   release to release.
#[pyfunction]
#[pyo3(name = "PatternFingerprint", signature = (mol, fpSize = 2048, atomCounts = None, setOnlyBits = None, tautomerFingerprints = false))]
#[allow(non_snake_case)]
fn pattern_fingerprint(
    mol: PyRef<'_, PyMol>,
    fpSize: u32,
    atomCounts: Option<&Bound<'_, PyAny>>,
    setOnlyBits: Option<PyRef<'_, PyExplicitBitVect>>,
    tautomerFingerprints: bool,
) -> PyResult<PyExplicitBitVect> {
    // upstream validates both arguments but neither changes the result
    read_atom_counts(atomCounts, mol.mol.num_atoms())?;
    if setOnlyBits.is_some_and(|b| b.bv.num_bits() != fpSize) {
        return Err(RDKitError::precondition("bad setOnlyBits size")).py();
    }
    let bv = pattern_fingerprint_mol(&mol.mol, fpSize, tautomerFingerprints).py()?;
    Ok(PyExplicitBitVect::wrap(bv))
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rdk_fingerprint, m)?)?;
    m.add_function(wrap_pyfunction!(pattern_fingerprint, m)?)?;
    m.add_function(wrap_pyfunction!(layered_fingerprint, m)?)?;
    m.add(
        "_LayeredFingerprint_version",
        LAYERED_FINGERPRINT_MOL_VERSION,
    )?;
    m.add("LayeredFingerprint_substructLayers", SUBSTRUCT_LAYERS)?;
    m.add(
        "_PatternFingerprint_version",
        PATTERN_FINGERPRINT_MOL_VERSION,
    )?;
    m.add_function(wrap_pyfunction!(unfolded_rdk_fingerprint, m)?)?;
    Ok(())
}
