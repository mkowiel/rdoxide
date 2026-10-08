//! Legacy path fingerprints of `rdkit.Chem.rdmolops` (`RDKFingerprint`,
//! `UnfoldedRDKFingerprintCountBased`; upstream `GraphMol/Wrap/MolOps.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use rdoxide::graphmol::fingerprints::legacy::{
    PathFpOutput, rdk_fingerprint_mol, unfolded_rdk_fingerprint_mol,
};

use crate::datastructs::PyExplicitBitVect;
use crate::datastructs_siv::PyULongSIV;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// `pythonObjectToVect`: `None` for falsy objects (None, 0, empty sequences).
fn opt_vect(o: Option<&Bound<'_, PyAny>>, max_v: Option<u32>) -> PyResult<Option<Vec<u32>>> {
    let Some(o) = o else { return Ok(None) };
    if !o.is_truthy()? {
        return Ok(None);
    }
    let mut v = Vec::new();
    for x in o.try_iter()? {
        let x: u32 = x?.extract()?;
        if let Some(m) = max_v
            && x >= m
        {
            return Err(PyValueError::new_err(
                "list element larger than allowed value",
            ));
        }
        v.push(x);
    }
    Ok(Some(v))
}

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
    let invs = opt_vect(atomInvariants, None)?;
    let from = opt_vect(fromAtoms, Some(n))?;
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
    let invs = opt_vect(atomInvariants, None)?;
    let from = opt_vect(fromAtoms, Some(n))?;
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

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rdk_fingerprint, m)?)?;
    m.add_function(wrap_pyfunction!(unfolded_rdk_fingerprint, m)?)?;
    Ok(())
}
