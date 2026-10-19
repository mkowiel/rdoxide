//! `rdoxide.Chem.rdMolAlign` (upstream `MolAlign/Wrap/rdMolAlign.cpp`):
//! alignment of molecules and conformers; the O3A functions are in
//! [`crate::rdmolalign_o3a`].

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rdoxide::geometry::Transform3D;
use rdoxide::graphmol::mol_align::{self as ma, BestAlignmentParams};
use rdoxide::graphmol::romol::ROMol;
use rdoxide::graphmol::substruct::MatchVect;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::pyconv::float64_array;
use crate::rdmolfiles_write::Args;

const DEFAULT_MAX_ITERS: u32 = 50;

/// `translateAtomMap`: a sequence of `(probe, reference)` pairs; empty is
/// "no map".
pub(crate) fn atom_map(obj: &Bound<'_, PyAny>) -> PyResult<Option<MatchVect>> {
    if obj.is_none() || obj.len()? == 0 {
        return Ok(None);
    }
    let mut res = MatchVect::new();
    for item in obj.try_iter()? {
        let item = item?;
        if item.len()? != 2 {
            return Err(PyValueError::new_err("Incorrect format for an atomMap"));
        }
        res.push((item.get_item(0)?.extract()?, item.get_item(1)?.extract()?));
    }
    Ok(Some(res))
}

/// `translateAtomMapSeq`: a sequence of atom maps.
fn atom_map_seq(obj: &Bound<'_, PyAny>) -> PyResult<Vec<MatchVect>> {
    obj.try_iter()?
        .map(|m| Ok(atom_map(&m?)?.unwrap_or_default()))
        .collect()
}

/// `translateDoubleSeq`: the values, or `None` for an empty sequence.
pub(crate) fn doubles(obj: &Bound<'_, PyAny>) -> PyResult<Option<Vec<f64>>> {
    if obj.is_none() || obj.len()? == 0 {
        return Ok(None);
    }
    obj.try_iter()?
        .map(|v| v?.extract())
        .collect::<PyResult<Vec<f64>>>()
        .map(Some)
}

/// `translateIntSeq`: the values, or `None` for an empty sequence.
fn unsigneds(obj: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Vec<u32>>> {
    match obj {
        Some(o) if !o.is_none() && o.len()? > 0 => o
            .try_iter()?
            .map(|v| v?.extract())
            .collect::<PyResult<Vec<u32>>>()
            .map(Some),
        _ => Ok(None),
    }
}

/// `generateRmsdTransMatchPyTuple`: `(rmsd, 4x4 matrix[, match])`.
pub(crate) fn result_tuple<'py>(
    py: Python<'py>,
    rmsd: f64,
    trans: &Transform3D,
    matched: Option<&MatchVect>,
) -> PyResult<Bound<'py, PyTuple>> {
    let matrix = float64_array(py, &trans.data)?.call_method1("reshape", ((4, 4),))?;
    let mut items = vec![rmsd.into_pyobject(py)?.into_any(), matrix];
    if let Some(m) = matched {
        items.push(PyTuple::new(py, m.iter().map(|&(p, r)| (p, r)))?.into_any());
    }
    PyTuple::new(py, items)
}

/// Runs `f` on a probe (mutable) and a reference molecule that may be the
/// same Python object.
fn with_mols<R>(
    prb: &Bound<'_, PyMol>,
    reference: &Bound<'_, PyMol>,
    f: impl FnOnce(&mut ROMol, &ROMol) -> PyResult<R>,
) -> PyResult<R> {
    if prb.is(reference) {
        let copy = prb.borrow().mol.clone();
        f(&mut prb.borrow_mut().mol, &copy)
    } else {
        let r = reference.borrow();
        let mut p = prb.borrow_mut();
        f(&mut p.mol, &r.mol)
    }
}

/// [`with_mols`] for functions that do not modify the probe.
fn with_mols_ref<R>(
    prb: &Bound<'_, PyMol>,
    reference: &Bound<'_, PyMol>,
    f: impl FnOnce(&ROMol, &ROMol) -> PyResult<R>,
) -> PyResult<R> {
    let (p, r) = (prb.borrow(), reference.borrow());
    f(&p.mol, &r.mol)
}

fn mol_arg<'a, 'py>(args: &Args<'a, 'py>, pos: usize, key: &str) -> PyResult<Bound<'py, PyMol>> {
    let name = args.name;
    args.get(pos, key)?
        .ok_or_else(|| {
            pyo3::exceptions::PyTypeError::new_err(format!(
                "{name}() missing required argument '{key}'"
            ))
        })?
        .cast_into::<PyMol>()
        .map_err(Into::into)
}

fn opt_flag(a: &Args<'_, '_>, pos: usize, key: &str, default: bool) -> PyResult<bool> {
    a.flag(pos, key, default)
}

/// Parameters controlling RMSD alignment
#[pyclass(
    name = "BestAlignmentParams",
    module = "rdoxide.Chem.rdMolAlign",
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
pub struct PyBestAlignmentParams {
    params: BestAlignmentParams,
}

impl PyBestAlignmentParams {
    /// The `pyBestAlignmentParams` constructor: map and weights must agree.
    fn from_parts(
        max_matches: i32,
        symmetrize: bool,
        num_threads: i32,
        map: Option<&Bound<'_, PyAny>>,
        weights: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<BestAlignmentParams> {
        let map = match map {
            Some(m) if !m.is_none() => atom_map_seq(m)?,
            _ => Vec::new(),
        };
        let weights = weights.map_or(Ok(None), doubles)?;
        if let (Some(w), Some(first)) = (&weights, map.first())
            && w.len() != first.len()
        {
            return Err(PyValueError::new_err(
                "Incorrect number of weights specified",
            ));
        }
        Ok(BestAlignmentParams {
            max_matches: max_matches as u32,
            symmetrize_conjugated_terminal_groups: symmetrize,
            ignore_hs: false,
            num_threads,
            map,
            weights,
        })
    }
}

#[pymethods]
impl PyBestAlignmentParams {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// maximum number of substructure matches to consider
    #[getter(maxMatches)]
    fn max_matches(&self) -> u32 {
        self.params.max_matches
    }
    #[setter(maxMatches)]
    fn set_max_matches(&mut self, v: i32) {
        self.params.max_matches = v as u32;
    }

    /// if true, conjugated terminal functional groups (like nitro or carboxylate) will be considered symmetrically.
    #[getter(symmetrizeConjugatedTerminalGroups)]
    fn symmetrize(&self) -> bool {
        self.params.symmetrize_conjugated_terminal_groups
    }
    #[setter(symmetrizeConjugatedTerminalGroups)]
    fn set_symmetrize(&mut self, v: bool) {
        self.params.symmetrize_conjugated_terminal_groups = v;
    }

    /// if true, hydrogens will be ignored in the alignment
    #[getter(ignoreHs)]
    fn ignore_hs(&self) -> bool {
        self.params.ignore_hs
    }
    #[setter(ignoreHs)]
    fn set_ignore_hs(&mut self, v: bool) {
        self.params.ignore_hs = v;
    }

    /// number of threads to use
    #[getter(numThreads)]
    fn num_threads(&self) -> i32 {
        self.params.num_threads
    }
    #[setter(numThreads)]
    fn set_num_threads(&mut self, v: i32) {
        self.params.num_threads = v;
    }

    /// the atom-atom mapping(s) used in the alignment
    #[getter]
    fn map<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let maps = self
            .params
            .map
            .iter()
            .map(|m| PyList::new(py, m.iter().map(|&(p, r)| (p, r))))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, maps)
    }
    #[setter]
    fn set_map(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        if !value.is_none() {
            self.params.map = atom_map_seq(value)?;
        }
        Ok(())
    }

    /// the weights used in the alignment
    #[getter]
    fn weights<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.params.weights.as_deref().unwrap_or_default())
    }
    #[setter]
    fn set_weights(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        if !value.is_none() {
            let w = doubles(value)?;
            if let (Some(w), Some(first)) = (&w, self.params.map.first())
                && w.len() != first.len()
            {
                return Err(PyValueError::new_err(
                    "Number of weights specified do not match number of atoms in the alignment map",
                ));
            }
            self.params.weights = w;
        }
        Ok(())
    }

    /// Only the documented attributes can be set (upstream's `safeSetattr`).
    fn __setattr__(slf: &Bound<'_, Self>, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        crate::pyconv::set_declared_attribute(slf.as_any(), name, value)
    }
}

/// The parameters of the `(…, params, …)` overloads, if `obj` holds them.
fn params_object(obj: Option<&Bound<'_, PyAny>>) -> Option<BestAlignmentParams> {
    obj?.cast::<PyBestAlignmentParams>()
        .ok()
        .map(|p| p.borrow().params.clone())
}

/// The legacy `(…, map, maxMatches, symmetrize…, weights, numThreads)`
/// argument group, starting at positional index `first`.
struct LegacyKeys {
    map: usize,
    max_matches: usize,
    symmetrize: usize,
    weights: usize,
    num_threads: usize,
}

fn legacy_params(a: &Args<'_, '_>, k: &LegacyKeys) -> PyResult<BestAlignmentParams> {
    let max_matches = a.int(k.max_matches, "maxMatches", 1_000_000)?;
    let symmetrize = a.flag(k.symmetrize, "symmetrizeConjugatedTerminalGroups", true)?;
    let num_threads = a.int(k.num_threads, "numThreads", 0)?;
    let map = a.get(k.map, "map")?;
    let weights = a.get(k.weights, "weights")?;
    PyBestAlignmentParams::from_parts(
        max_matches,
        symmetrize,
        num_threads,
        map.as_ref(),
        weights.as_ref(),
    )
}

/// Checks the weight count of the `GetAlignmentTransform`/`AlignMol` calls.
fn check_weights(
    prb: &Bound<'_, PyMol>,
    map: &Option<MatchVect>,
    weights: &Option<Vec<f64>>,
) -> PyResult<()> {
    let n_atoms = map
        .as_ref()
        .map_or_else(|| prb.borrow().mol.num_atoms() as usize, Vec::len);
    match weights {
        Some(w) if w.len() != n_atoms => Err(PyValueError::new_err(
            "Incorrect number of weights specified",
        )),
        _ => Ok(()),
    }
}

/// Compute the transformation required to align a molecule
///
/// The 3D transformation required to align the specied conformation in the probe molecule
/// to a specified conformation in the reference molecule is computed so that the root mean
/// squared distance between a specified set of atoms is minimized
///
/// ARGUMENTS
///  - prbMol    molecule that is to be aligned
///  - refMol    molecule used as the reference for the alignment
///  - prbCid    ID of the conformation in the probe to be used
///                   for the alignment (defaults to first conformation)
///  - refCid    ID of the conformation in the ref molecule to which
///                   the alignment is computed (defaults to first conformation)
///  - atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)
///                   used to compute the alignments. If this mapping is
///                   not specified an attempt is made to generate one by
///                   substructure matching
///  - weights   Optionally specify weights for each of the atom pairs
///  - reflect   if true reflect the conformation of the probe molecule
///  - maxIters  maximum number of iterations used in minimizing the RMSD
///
/// RETURNS
///  a tuple of (RMSD value, transform matrix)
#[pyfunction]
#[pyo3(name = "GetAlignmentTransform", signature = (prbMol, refMol, prbCid = -1, refCid = -1, atomMap = None, weights = None, reflect = false, maxIters = 50))]
#[allow(clippy::too_many_arguments, non_snake_case)]
fn get_alignment_transform<'py>(
    py: Python<'py>,
    prbMol: &Bound<'py, PyMol>,
    refMol: &Bound<'py, PyMol>,
    prbCid: i32,
    refCid: i32,
    atomMap: Option<&Bound<'py, PyAny>>,
    weights: Option<&Bound<'py, PyAny>>,
    reflect: bool,
    maxIters: u32,
) -> PyResult<Bound<'py, PyTuple>> {
    let map = atomMap.map_or(Ok(None), atom_map)?;
    let weights = weights.map_or(Ok(None), doubles)?;
    check_weights(prbMol, &map, &weights)?;
    let (rmsd, trans) = with_mols_ref(prbMol, refMol, |p, r| {
        py.detach(|| {
            ma::get_alignment_transform(
                p,
                r,
                prbCid,
                refCid,
                map.as_ref(),
                weights.as_deref(),
                reflect,
                maxIters,
            )
        })
        .py()
    })?;
    result_tuple(py, rmsd, &trans, None)
}

/// Compute the optimal RMS, transformation and atom map for aligning
/// two molecules, taking symmetry into account. Molecule coordinates
/// are left unaltered.
///
/// This function will attempt to align all permutations of matching atom
/// orders in both molecules, for some molecules it will lead to 'combinatorial
/// explosion' especially if hydrogens are present.
/// Use 'GetAlignmentTransform' to align molecules without changing the atom order.
///
/// Overloads:
///   GetBestAlignmentTransform(prbMol, refMol, prbCid=-1, refCid=-1, map=[], maxMatches=1000000,
///     symmetrizeConjugatedTerminalGroups=True, weights=[], reflect=False, maxIters=50, numThreads=1)
///   GetBestAlignmentTransform(prbMol, refMol, params, prbCid=-1, refCid=-1, reflect=False, maxIters=50)
///
/// RETURNS
///  a tuple of (RMSD value, best transform matrix, best atom map)
#[pyfunction]
#[pyo3(name = "GetBestAlignmentTransform", signature = (*args, **kwargs))]
fn get_best_alignment_transform<'py>(
    py: Python<'py>,
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let a = Args {
        name: "GetBestAlignmentTransform",
        args,
        kwargs,
    };
    let prb = mol_arg(&a, 0, "prbMol")?;
    let reference = mol_arg(&a, 1, "refMol")?;
    let (params, prb_cid, ref_cid, reflect, max_iters) =
        if let Some(p) = params_object(a.get(2, "params")?.as_ref()) {
            (
                p,
                a.int(3, "prbCid", -1)?,
                a.int(4, "refCid", -1)?,
                opt_flag(&a, 5, "reflect", false)?,
                a.int(6, "maxIters", DEFAULT_MAX_ITERS as i32)? as u32,
            )
        } else {
            let keys = LegacyKeys {
                map: 4,
                max_matches: 5,
                symmetrize: 6,
                weights: 7,
                num_threads: 10,
            };
            (
                legacy_params(&a, &keys)?,
                a.int(2, "prbCid", -1)?,
                a.int(3, "refCid", -1)?,
                opt_flag(&a, 8, "reflect", false)?,
                a.int(9, "maxIters", DEFAULT_MAX_ITERS as i32)? as u32,
            )
        };
    let (rmsd, trans, best) = with_mols_ref(&prb, &reference, |p, r| {
        py.detach(|| {
            ma::get_best_alignment_transform(p, r, &params, prb_cid, ref_cid, reflect, max_iters)
        })
        .py()
    })?;
    result_tuple(py, rmsd, &trans, Some(&best))
}

/// Optimally (minimum RMSD) align a molecule to another molecule
///
/// The 3D transformation required to align the specied conformation in the probe molecule
/// to a specified conformation in the reference molecule is computed so that the root mean
/// squared distance between a specified set of atoms is minimized. This transform is then
/// applied to the specified conformation in the probe molecule
///
/// ARGUMENTS
///  - prbMol    molecule that is to be aligned
///  - refMol    molecule used as the reference for the alignment
///  - prbCid    ID of the conformation in the probe to be used
///                   for the alignment (defaults to first conformation)
///  - refCid    ID of the conformation in the ref molecule to which
///                   the alignment is computed (defaults to first conformation)
///  - atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)
///                   used to compute the alignments. If this mapping is
///                   not specified an attempt is made to generate one by
///                   substructure matching
///  - weights   Optionally specify weights for each of the atom pairs
///  - reflect   if true reflect the conformation of the probe molecule
///  - maxIters  maximum number of iterations used in minimizing the RMSD
///
/// RETURNS
/// RMSD value
#[pyfunction]
#[pyo3(name = "AlignMol", signature = (prbMol, refMol, prbCid = -1, refCid = -1, atomMap = None, weights = None, reflect = false, maxIters = 50))]
#[allow(clippy::too_many_arguments, non_snake_case)]
fn align_mol(
    py: Python<'_>,
    prbMol: &Bound<'_, PyMol>,
    refMol: &Bound<'_, PyMol>,
    prbCid: i32,
    refCid: i32,
    atomMap: Option<&Bound<'_, PyAny>>,
    weights: Option<&Bound<'_, PyAny>>,
    reflect: bool,
    maxIters: u32,
) -> PyResult<f64> {
    let map = atomMap.map_or(Ok(None), atom_map)?;
    let weights = weights.map_or(Ok(None), doubles)?;
    check_weights(prbMol, &map, &weights)?;
    with_mols(prbMol, refMol, |p, r| {
        py.detach(|| {
            ma::align_mol(
                p,
                r,
                prbCid,
                refCid,
                map.as_ref(),
                weights.as_deref(),
                reflect,
                maxIters,
            )
        })
        .py()
    })
}

/// Returns the optimal RMS for aligning two molecules, taking
/// symmetry into account. As a side-effect, the probe molecule is
/// left in the aligned state.
///
/// Overloads:
///   GetBestRMS(prbMol, refMol, prbId=-1, refId=-1, map=None, maxMatches=1000000,
///     symmetrizeConjugatedTerminalGroups=True, weights=[], numThreads=1)
///   GetBestRMS(prbMol, refMol, params, prbId=-1, refId=-1)
#[pyfunction]
#[pyo3(name = "GetBestRMS", signature = (*args, **kwargs))]
fn get_best_rms(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<f64> {
    let a = Args {
        name: "GetBestRMS",
        args,
        kwargs,
    };
    let prb = mol_arg(&a, 0, "prbMol")?;
    let reference = mol_arg(&a, 1, "refMol")?;
    let (params, prb_id, ref_id) = if let Some(p) = params_object(a.get(2, "params")?.as_ref()) {
        (p, a.int(3, "prbId", -1)?, a.int(4, "refId", -1)?)
    } else {
        let keys = LegacyKeys {
            map: 4,
            max_matches: 5,
            symmetrize: 6,
            weights: 7,
            num_threads: 8,
        };
        (
            legacy_params(&a, &keys)?,
            a.int(2, "prbId", -1)?,
            a.int(3, "refId", -1)?,
        )
    };
    with_mols(&prb, &reference, |p, r| {
        py.detach(|| ma::get_best_rms(p, r, &params, prb_id, ref_id))
            .py()
    })
}

fn rms_tuple<'py>(py: Python<'py>, values: Vec<f64>) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(py, values)
}

/// Returns the symmetric distance matrix between the conformers of a molecule.
/// getBestRMS() is used to calculate the inter-conformer distances
///
/// Overloads:
///   GetAllConformerBestRMS(mol, numThreads=1, map=None, maxMatches=1000000,
///     symmetrizeConjugatedTerminalGroups=True, weights=[])
///   GetAllConformerBestRMS(mol, params)
///
/// RETURNS
///  a tuple with the lower triangle of the matrix: (1,0), (2,0), (2,1), …
#[pyfunction]
#[pyo3(name = "GetAllConformerBestRMS", signature = (*args, **kwargs))]
fn get_all_conformer_best_rms<'py>(
    py: Python<'py>,
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let a = Args {
        name: "GetAllConformerBestRMS",
        args,
        kwargs,
    };
    let mol = mol_arg(&a, 0, "mol")?;
    let params = match params_object(a.get(1, "params")?.as_ref()) {
        Some(p) => p,
        None => legacy_params(
            &a,
            &LegacyKeys {
                map: 2,
                max_matches: 3,
                symmetrize: 4,
                weights: 5,
                num_threads: 1,
            },
        )?,
    };
    let m = mol.borrow();
    let core = &m.mol;
    rms_tuple(
        py,
        py.detach(|| ma::get_all_conformer_best_rms(core, &params))
            .py()?,
    )
}

/// Returns the best RMS of every probe conformer against every reference
/// conformer (reference conformers in the outer loop).
#[pyfunction]
#[pyo3(name = "GetAllConformerBestRMSToRef", signature = (prbMol, refMol, params = None))]
#[allow(non_snake_case)]
fn get_all_conformer_best_rms_to_ref<'py>(
    py: Python<'py>,
    prbMol: &Bound<'py, PyMol>,
    refMol: &Bound<'py, PyMol>,
    params: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let params = match params.filter(|p| !p.is_none()) {
        None => BestAlignmentParams::default(),
        Some(p) => params_object(Some(p)).ok_or_else(|| {
            pyo3::exceptions::PyTypeError::new_err(
                "params must be a BestAlignmentParams object or None",
            )
        })?,
    };
    let values = with_mols_ref(prbMol, refMol, |p, r| {
        py.detach(|| ma::get_all_conformer_best_rms_to_ref(p, r, &params))
            .py()
    })?;
    rms_tuple(py, values)
}

/// Returns the RMS between two molecules, taking symmetry into account.
/// In contrast to getBestRMS, the RMS is computed 'in place', i.e.
/// probe molecules are not aligned to the reference ahead of the
/// RMS calculation.
#[pyfunction]
#[pyo3(name = "CalcRMS", signature = (prbMol, refMol, prbId = -1, refId = -1, map = None, maxMatches = 1000000, symmetrizeConjugatedTerminalGroups = true, weights = None))]
#[allow(clippy::too_many_arguments, non_snake_case)]
fn calc_rms(
    py: Python<'_>,
    prbMol: &Bound<'_, PyMol>,
    refMol: &Bound<'_, PyMol>,
    prbId: i32,
    refId: i32,
    map: Option<&Bound<'_, PyAny>>,
    maxMatches: i32,
    symmetrizeConjugatedTerminalGroups: bool,
    weights: Option<&Bound<'_, PyAny>>,
) -> PyResult<f64> {
    let maps = match map {
        Some(m) if !m.is_none() => atom_map_seq(m)?,
        _ => Vec::new(),
    };
    let weights = weights.map_or(Ok(None), doubles)?;
    with_mols_ref(prbMol, refMol, |p, r| {
        py.detach(|| {
            ma::calc_rms(
                p,
                r,
                prbId,
                refId,
                &maps,
                maxMatches as u32,
                symmetrizeConjugatedTerminalGroups,
                weights.as_deref(),
            )
        })
        .py()
    })
}

/// Align conformations in a molecule to each other
///
/// The first conformation in the molecule is used as the reference
///
/// ARGUMENTS
///  - mol          molecule of interest
///  - atomIds      List of atom ids to use a points for alignment - defaults to all atoms
///  - confIds      Ids of conformations to align - defaults to all conformers
///  - weights      Optionally specify weights for each of the atom pairs
///  - reflect      if true reflect the conformation of the probe molecule
///  - maxIters     maximum number of iterations used in minimizing the RMSD
///  - RMSlist      if provided, fills in the RMS values between the reference
///                 conformation and the other aligned conformations
#[pyfunction]
#[pyo3(name = "AlignMolConformers", signature = (mol, atomIds = None, confIds = None, weights = None, reflect = false, maxIters = 50, RMSlist = None))]
#[allow(clippy::too_many_arguments, non_snake_case)]
fn align_mol_conformers(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    atomIds: Option<&Bound<'_, PyAny>>,
    confIds: Option<&Bound<'_, PyAny>>,
    weights: Option<&Bound<'_, PyAny>>,
    reflect: bool,
    maxIters: u32,
    RMSlist: Option<&Bound<'_, PyAny>>,
) -> PyResult<()> {
    let weights = weights.map_or(Ok(None), doubles)?;
    let atom_ids = unsigneds(atomIds)?;
    let conf_ids = unsigneds(confIds)?;
    let mut m = mol.borrow_mut();
    let core = &mut m.mol;
    let rms = py
        .detach(|| {
            ma::align_mol_conformers(
                core,
                atom_ids.as_deref(),
                conf_ids.as_deref(),
                weights.as_deref(),
                reflect,
                maxIters,
            )
        })
        .py()?;
    if let Some(list) = RMSlist.filter(|l| !l.is_none()) {
        for v in rms {
            list.call_method1("append", (v,))?;
        }
    }
    Ok(())
}

/// Perform a random transformation on a molecule
///
/// ARGUMENTS
///  - mol    molecule that is to be transformed
///  - cid    ID of the conformation in the mol to be transformed
///           (defaults to first conformation)
///  - seed   seed used to initialize the random generator
///           (defaults to -1, that is no seeding)
#[pyfunction]
#[pyo3(name = "RandomTransform", signature = (mol, cid = -1, seed = -1))]
fn random_transform(mol: &Bound<'_, PyMol>, cid: i32, seed: i32) -> PyResult<()> {
    ma::random_transform(&mut mol.borrow_mut().mol, cid, seed).py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions to align a molecule to a second molecule",
    )?;
    m.add_class::<PyBestAlignmentParams>()?;
    m.add_function(wrap_pyfunction!(get_alignment_transform, m)?)?;
    m.add_function(wrap_pyfunction!(get_best_alignment_transform, m)?)?;
    m.add_function(wrap_pyfunction!(align_mol, m)?)?;
    m.add_function(wrap_pyfunction!(get_best_rms, m)?)?;
    m.add_function(wrap_pyfunction!(get_all_conformer_best_rms, m)?)?;
    m.add_function(wrap_pyfunction!(get_all_conformer_best_rms_to_ref, m)?)?;
    m.add_function(wrap_pyfunction!(calc_rms, m)?)?;
    m.add_function(wrap_pyfunction!(align_mol_conformers, m)?)?;
    m.add_function(wrap_pyfunction!(random_transform, m)?)?;
    crate::rdmolalign_o3a::init(m)
}
