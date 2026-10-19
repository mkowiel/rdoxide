//! The Open3DALIGN functions of `rdoxide.Chem.rdMolAlign` (upstream
//! `MolAlign/Wrap/rdMolAlign.cpp`): `GetO3A`, `GetCrippenO3A`,
//! `GetO3AForProbeConfs`, `GetCrippenO3AForProbeConfs` and the `O3A` class.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use rdoxide::graphmol::descriptors::crippen::get_crippen_atom_contribs;
use rdoxide::graphmol::forcefield_helpers::mmff::MmffMolProperties;
use rdoxide::graphmol::mol_align_o3a::{self as o3a, O3A, O3AParams, O3AProps};
use rdoxide::graphmol::romol::ROMol;
use rdoxide::graphmol::substruct::MatchVect;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdforcefield::PyMmffMolProperties;
use crate::rdmolalign::{atom_map, doubles, result_tuple};

/// Open3DALIGN object
#[pyclass(name = "O3A", module = "rdoxide.Chem.rdMolAlign")]
pub struct PyO3A {
    o3a: O3A,
    prb: Py<PyMol>,
    reference: Py<PyMol>,
}

#[pymethods]
impl PyO3A {
    /// aligns probe molecule onto reference molecule
    #[pyo3(name = "Align")]
    fn align(&self, py: Python<'_>) -> PyResult<f64> {
        let reference = self.reference.bind(py);
        if self.prb.is(reference) {
            let copy = reference.borrow().mol.clone();
            return self
                .o3a
                .align(&mut self.prb.bind(py).borrow_mut().mol, &copy)
                .py();
        }
        let r = reference.borrow();
        self.o3a
            .align(&mut self.prb.bind(py).borrow_mut().mol, &r.mol)
            .py()
    }

    /// returns the transformation which aligns probe molecule onto reference molecule
    #[pyo3(name = "Trans")]
    fn trans<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let (p, r) = (self.prb.bind(py).borrow(), self.reference.bind(py).borrow());
        let (rmsd, trans) = self.o3a.trans(&p.mol, &r.mol).py()?;
        result_tuple(py, rmsd, &trans, None)
    }

    /// returns the O3AScore of the alignment
    #[pyo3(name = "Score")]
    fn score(&self) -> f64 {
        self.o3a.score
    }

    /// returns the AtomMap as found by Open3DALIGN
    #[pyo3(name = "Matches")]
    fn matches<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(py, self.o3a.matches.iter().map(|&(p, r)| vec![p, r]))
    }

    /// returns the weights of the atom-to-atom matches found by Open3DALIGN
    #[pyo3(name = "Weights")]
    fn weights(&self) -> Vec<f64> {
        self.o3a.weights.clone()
    }
}

/// Validated `constraintMap`/`constraintWeights` arguments.
fn constraints(
    prb: &ROMol,
    reference: &ROMol,
    map: Option<&Bound<'_, PyAny>>,
    weights: Option<&Bound<'_, PyAny>>,
) -> PyResult<(Option<MatchVect>, Option<Vec<f64>>)> {
    let map = match map {
        Some(m) => atom_map(m)?,
        None => None,
    };
    let Some(map_ref) = &map else {
        return Ok((None, None));
    };
    let weights = weights.map_or(Ok(None), doubles)?;
    if let Some(w) = &weights
        && w.len() != map_ref.len()
    {
        return Err(PyValueError::new_err(
            "The number of weights should match the number of constraints",
        ));
    }
    for &(p, r) in map_ref {
        if p < 0 || p >= prb.num_atoms() as i32 || r < 0 || r >= reference.num_atoms() as i32 {
            return Err(PyValueError::new_err("Constrained atom idx out of range"));
        }
        if prb.atom(p as u32).atomic_num() == 1 || reference.atom(r as u32).atomic_num() == 1 {
            return Err(PyValueError::new_err(
                "Constrained atoms must be heavy atoms",
            ));
        }
    }
    Ok((map, weights))
}

/// The MMFF properties passed in, or computed from the molecule.
fn mmff_props(
    props: Option<&Bound<'_, PyAny>>,
    mol: &ROMol,
    which: &str,
) -> PyResult<MmffMolProperties> {
    if let Some(p) = props.filter(|p| !p.is_none()) {
        return Ok(p.extract::<PyRef<'_, PyMmffMolProperties>>()?.props.clone());
    }
    let mut copy = mol.clone();
    let props = MmffMolProperties::new(&mut copy, "MMFF94").py()?;
    if !props.is_valid() {
        return Err(PyValueError::new_err(format!(
            "missing MMFF94 parameters for {which} molecule"
        )));
    }
    Ok(props)
}

/// Crippen logP contributions from the `(logP, MR)` list passed in (when it
/// has one entry per atom) or computed.
fn crippen_logp(contribs: Option<&Bound<'_, PyAny>>, mol: &mut ROMol) -> PyResult<Vec<f64>> {
    let n = mol.num_atoms() as usize;
    if let Some(c) = contribs.filter(|c| !c.is_none() && c.len().is_ok_and(|l| l == n)) {
        return c
            .try_iter()?
            .map(|t| t?.get_item(0)?.extract::<f64>())
            .collect();
    }
    Ok(get_crippen_atom_contribs(mol, true).py()?.logp)
}

#[allow(clippy::too_many_arguments)]
fn wrap(
    py: Python<'_>,
    result: O3A,
    prb: &Bound<'_, PyMol>,
    reference: &Bound<'_, PyMol>,
) -> PyResult<Py<PyO3A>> {
    Py::new(
        py,
        PyO3A {
            o3a: result,
            prb: prb.clone().unbind(),
            reference: reference.clone().unbind(),
        },
    )
}

/// Get an O3A object with atomMap and weights vectors to overlay
/// the probe molecule onto the reference molecule based on
/// MMFF atom types and charges
///
///   ARGUMENTS
///     - prbMol         molecule that is to be aligned
///     - refMol         molecule used as the reference for the alignment
///     - prbPyMMFFMolProperties   PyMMFFMolProperties object for the probe molecule as returned
///                      by MMFFGetMoleculeProperties()
///     - refPyMMFFMolProperties   PyMMFFMolProperties object for the reference molecule as returned
///                      by MMFFGetMoleculeProperties()
///     - prbCid         ID of the conformation in the probe to be used
///                      for the alignment (defaults to first conformation)
///     - refCid         ID of the conformation in the ref molecule to which
///                      the alignment is computed (defaults to first conformation)
///     - reflect        if true reflect the conformation of the probe molecule
///                      (defaults to false)
///     - maxIters       maximum number of iterations used in minimizing the RMSD
///                      (defaults to 50)
///     - options        least 2 significant bits encode accuracy
///                      (0: maximum, 3: minimum; default: 0)
///                      bit 3 triggers local optimization of the alignment
///                      (no global search); default: 0 (global search + local optimization)
///                      bit 4 triggers the use of the MMFF weights... see upstream docs
///     - constraintMap  a vector of pairs of atom IDs (probe AtomId, ref AtomId)
///                      which shall be used for the alignment (defaults to [])
///     - constraintWeights  optionally specify weights for each of the constraints
///                      (weights default to 100.0)
///
///   RETURNS
///     The O3A object
#[pyfunction]
#[pyo3(name = "GetO3A", signature = (prbMol, refMol, prbPyMMFFMolProperties = None, refPyMMFFMolProperties = None, prbCid = -1, refCid = -1, reflect = false, maxIters = 50, options = 0, constraintMap = None, constraintWeights = None))]
#[allow(clippy::too_many_arguments, non_snake_case)]
fn get_o3a(
    py: Python<'_>,
    prbMol: &Bound<'_, PyMol>,
    refMol: &Bound<'_, PyMol>,
    prbPyMMFFMolProperties: Option<&Bound<'_, PyAny>>,
    refPyMMFFMolProperties: Option<&Bound<'_, PyAny>>,
    prbCid: i32,
    refCid: i32,
    reflect: bool,
    maxIters: u32,
    options: u32,
    constraintMap: Option<&Bound<'_, PyAny>>,
    constraintWeights: Option<&Bound<'_, PyAny>>,
) -> PyResult<Py<PyO3A>> {
    let (prb, reference) = (prbMol.borrow().mol.clone(), refMol.borrow().mol.clone());
    let (map, weights) = constraints(&prb, &reference, constraintMap, constraintWeights)?;
    let prb_props = mmff_props(prbPyMMFFMolProperties, &prb, "probe")?;
    let ref_props = mmff_props(refPyMMFFMolProperties, &reference, "reference")?;
    let params = O3AParams {
        prb_cid: prbCid,
        ref_cid: refCid,
        reflect,
        max_iters: maxIters,
        options,
        constraint_map: map.as_ref(),
        constraint_weights: weights.as_deref(),
    };
    let props = O3AProps::Mmff {
        prb: &prb_props,
        reference: &ref_props,
    };
    let result = py
        .detach(|| O3A::new(&prb, &reference, &props, &params))
        .py()?;
    wrap(py, result, prbMol, refMol)
}

/// Get an O3A object with atomMap and weights vectors to overlay
/// the probe molecule onto the reference molecule based on
/// Crippen logP atom contributions
///
///   ARGUMENTS
///     - prbMol         molecule that is to be aligned
///     - refMol         molecule used as the reference for the alignment
///     - prbCrippenContribs  Crippen atom contributions for the probe molecule
///                      as a list of (logp, mr) tuples, as returned by
///                      _CalcCrippenContribs()
///     - refCrippenContribs  Crippen atom contributions for the reference molecule
///                      as a list of (logp, mr) tuples, as returned by
///                      _CalcCrippenContribs()
///     - prbCid, refCid, reflect, maxIters, options, constraintMap,
///       constraintWeights: as for GetO3A
///
///   RETURNS
///     The O3A object
#[pyfunction]
#[pyo3(name = "GetCrippenO3A", signature = (prbMol, refMol, prbCrippenContribs = None, refCrippenContribs = None, prbCid = -1, refCid = -1, reflect = false, maxIters = 50, options = 0, constraintMap = None, constraintWeights = None))]
#[allow(clippy::too_many_arguments, non_snake_case)]
fn get_crippen_o3a(
    py: Python<'_>,
    prbMol: &Bound<'_, PyMol>,
    refMol: &Bound<'_, PyMol>,
    prbCrippenContribs: Option<&Bound<'_, PyAny>>,
    refCrippenContribs: Option<&Bound<'_, PyAny>>,
    prbCid: i32,
    refCid: i32,
    reflect: bool,
    maxIters: u32,
    options: u32,
    constraintMap: Option<&Bound<'_, PyAny>>,
    constraintWeights: Option<&Bound<'_, PyAny>>,
) -> PyResult<Py<PyO3A>> {
    let (mut prb, mut reference) = (prbMol.borrow().mol.clone(), refMol.borrow().mol.clone());
    let (map, weights) = constraints(&prb, &reference, constraintMap, constraintWeights)?;
    let prb_logp = crippen_logp(prbCrippenContribs, &mut prb)?;
    let ref_logp = crippen_logp(refCrippenContribs, &mut reference)?;
    let params = O3AParams {
        prb_cid: prbCid,
        ref_cid: refCid,
        reflect,
        max_iters: maxIters,
        options,
        constraint_map: map.as_ref(),
        constraint_weights: weights.as_deref(),
    };
    let props = O3AProps::Crippen {
        prb: &prb_logp,
        reference: &ref_logp,
    };
    let result = py
        .detach(|| O3A::new(&prb, &reference, &props, &params))
        .py()?;
    wrap(py, result, prbMol, refMol)
}

fn wrap_all<'py>(
    py: Python<'py>,
    results: Vec<O3A>,
    prb: &Bound<'py, PyMol>,
    reference: &Bound<'py, PyMol>,
) -> PyResult<Bound<'py, PyTuple>> {
    let items = results
        .into_iter()
        .map(|r| wrap(py, r, prb, reference))
        .collect::<PyResult<Vec<_>>>()?;
    PyTuple::new(py, items)
}

/// Get a vector of O3A objects for the overlay of all
/// the probe molecule's conformations onto the reference molecule based on
/// MMFF atom types and charges
///
///   ARGUMENTS
///     - prbMol, refMol, numThreads, prbPyMMFFMolProperties, refPyMMFFMolProperties,
///       refCid, reflect, maxIters, options, constraintMap, constraintWeights:
///       as for GetO3A (numThreads = 0 uses all available threads)
///
///   RETURNS
///     A vector of O3A objects
#[pyfunction]
#[pyo3(name = "GetO3AForProbeConfs", signature = (prbMol, refMol, numThreads = 1, prbPyMMFFMolProperties = None, refPyMMFFMolProperties = None, refCid = -1, reflect = false, maxIters = 50, options = 0, constraintMap = None, constraintWeights = None))]
#[allow(clippy::too_many_arguments, non_snake_case)]
fn get_o3a_for_probe_confs<'py>(
    py: Python<'py>,
    prbMol: &Bound<'py, PyMol>,
    refMol: &Bound<'py, PyMol>,
    numThreads: i32,
    prbPyMMFFMolProperties: Option<&Bound<'py, PyAny>>,
    refPyMMFFMolProperties: Option<&Bound<'py, PyAny>>,
    refCid: i32,
    reflect: bool,
    maxIters: u32,
    options: u32,
    constraintMap: Option<&Bound<'py, PyAny>>,
    constraintWeights: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let (prb, reference) = (prbMol.borrow().mol.clone(), refMol.borrow().mol.clone());
    let (map, weights) = constraints(&prb, &reference, constraintMap, constraintWeights)?;
    let prb_props = mmff_props(prbPyMMFFMolProperties, &prb, "probe")?;
    let ref_props = mmff_props(refPyMMFFMolProperties, &reference, "reference")?;
    let params = O3AParams {
        prb_cid: -1,
        ref_cid: refCid,
        reflect,
        max_iters: maxIters,
        options,
        constraint_map: map.as_ref(),
        constraint_weights: weights.as_deref(),
    };
    let props = O3AProps::Mmff {
        prb: &prb_props,
        reference: &ref_props,
    };
    let results = py
        .detach(|| o3a::get_o3a_for_probe_confs(&prb, &reference, &props, numThreads, &params))
        .py()?;
    wrap_all(py, results, prbMol, refMol)
}

/// Get a vector of O3A objects for the overlay of all
/// the probe molecule's conformations onto the reference molecule based on
/// Crippen logP atom contributions
///
///   ARGUMENTS: as for GetO3AForProbeConfs, with Crippen contributions
///     (prbCrippenContribs, refCrippenContribs) in place of the MMFF properties
///
///   RETURNS
///     A vector of O3A objects
#[pyfunction]
#[pyo3(name = "GetCrippenO3AForProbeConfs", signature = (prbMol, refMol, numThreads = 1, prbCrippenContribs = None, refCrippenContribs = None, refCid = -1, reflect = false, maxIters = 50, options = 0, constraintMap = None, constraintWeights = None))]
#[allow(clippy::too_many_arguments, non_snake_case)]
fn get_crippen_o3a_for_probe_confs<'py>(
    py: Python<'py>,
    prbMol: &Bound<'py, PyMol>,
    refMol: &Bound<'py, PyMol>,
    numThreads: i32,
    prbCrippenContribs: Option<&Bound<'py, PyAny>>,
    refCrippenContribs: Option<&Bound<'py, PyAny>>,
    refCid: i32,
    reflect: bool,
    maxIters: u32,
    options: u32,
    constraintMap: Option<&Bound<'py, PyAny>>,
    constraintWeights: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let (mut prb, mut reference) = (prbMol.borrow().mol.clone(), refMol.borrow().mol.clone());
    let (map, weights) = constraints(&prb, &reference, constraintMap, constraintWeights)?;
    let prb_logp = crippen_logp(prbCrippenContribs, &mut prb)?;
    let ref_logp = crippen_logp(refCrippenContribs, &mut reference)?;
    let params = O3AParams {
        prb_cid: -1,
        ref_cid: refCid,
        reflect,
        max_iters: maxIters,
        options,
        constraint_map: map.as_ref(),
        constraint_weights: weights.as_deref(),
    };
    let props = O3AProps::Crippen {
        prb: &prb_logp,
        reference: &ref_logp,
    };
    let results = py
        .detach(|| o3a::get_o3a_for_probe_confs(&prb, &reference, &props, numThreads, &params))
        .py()?;
    wrap_all(py, results, prbMol, refMol)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyO3A>()?;
    m.add_function(wrap_pyfunction!(get_o3a, m)?)?;
    m.add_function(wrap_pyfunction!(get_crippen_o3a, m)?)?;
    m.add_function(wrap_pyfunction!(get_o3a_for_probe_confs, m)?)?;
    m.add_function(wrap_pyfunction!(get_crippen_o3a_for_probe_confs, m)?)?;
    Ok(())
}
