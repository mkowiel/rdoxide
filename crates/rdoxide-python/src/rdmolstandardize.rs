//! `rdoxide.Chem.MolStandardize.rdMolStandardize` (upstream
//! `GraphMol/MolStandardize/nbWrap/rdMolStandardize.cpp`).

use pyo3::exceptions::{PyRuntimeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::graphmol::mol_standardize as ms;
use rdoxide::graphmol::mol_standardize::charge::{
    ChargeCorrection, Reionizer, Uncharger, charge_corrections,
};
use rdoxide::graphmol::mol_standardize::fragment::{FragmentRemover, LargestFragmentChooser};
use rdoxide::graphmol::mol_standardize::metal::{MetalDisconnector, MetalDisconnectorOptions};
use rdoxide::graphmol::mol_standardize::normalize::{Normalizer, normalizer_from_params};
use rdoxide::graphmol::mol_standardize::validate::{
    FeaturesValidation, Layout2DValidation, ValidationMethod, mol_vs_validation, validate_smiles,
};
use rdoxide::graphmol::mol_standardize::{
    CleanupParameters, default_cleanup_parameters, standardize_multiple_mols_in_place,
    update_cleanup_params_from_json,
};
use rdoxide::graphmol::smiles_parse::smarts_write::mol_to_smarts;
use rdoxide::graphmol::smiles_parse::write::SmilesWriteParams;
use rdoxide::graphmol::{Atom, RWMol};

use crate::atom::PyAtom;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Parameters controlling molecular standardization
#[pyclass(
    name = "CleanupParameters",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    skip_from_py_object
)]
#[derive(Clone, Debug, Default)]
pub struct PyCleanupParameters {
    pub(crate) params: CleanupParameters,
}

/// A read/write attribute of the parameters, as `def_rw` makes it.
macro_rules! attributes {
    ($($(#[$m:meta])* $py:literal, $get:ident, $set:ident, $field:ident: $ty:ty;)*) => {
        #[pymethods]
        impl PyCleanupParameters {
            $(
                $(#[$m])*
                #[getter($py)]
                fn $get(&self) -> $ty {
                    self.params.$field.clone()
                }
                #[setter($py)]
                fn $set(&mut self, value: $ty) {
                    self.params.$field = value;
                }
            )*
        }
    };
}

attributes! {
    /// file containing the normalization transformations
    "normalizationsFile", get_normalizations, set_normalizations, normalizations: String;
    /// file containing the acid and base definitions
    "acidbaseFile", get_acidbase_file, set_acidbase_file, acidbase_file: String;
    /// file containing the acid and base definitions
    "fragmentFile", get_fragment_file, set_fragment_file, fragment_file: String;
    /// file containing the tautomer transformations
    "tautomerTransformsFile", get_tautomer_transforms, set_tautomer_transforms, tautomer_transforms: String;
    /// maximum number of restarts
    "maxRestarts", get_max_restarts, set_max_restarts, max_restarts: i32;
    /// prefer organic fragments to inorganic ones when deciding what to keep
    "preferOrganic", get_prefer_organic, set_prefer_organic, prefer_organic: bool;
    /// apply atom-order dependent normalizations (like uncharging) in a canonical order
    "doCanonical", get_do_canonical, set_do_canonical, do_canonical: bool;
    /// maximum number of tautomers to generate (defaults to 1000)
    "maxTautomers", get_max_tautomers, set_max_tautomers, max_tautomers: i32;
    /// maximum number of transforms to apply during tautomer enumeration (defaults to 1000)
    "maxTransforms", get_max_transforms, set_max_transforms, max_transforms: i32;
    /// remove stereochemistry from sp3 centers involved in tautomerism (defaults to True)
    "tautomerRemoveSp3Stereo", get_tautomer_remove_sp3_stereo, set_tautomer_remove_sp3_stereo, tautomer_remove_sp3_stereo: bool;
    /// remove stereochemistry from double bonds involved in tautomerism (defaults to True)
    "tautomerRemoveBondStereo", get_tautomer_remove_bond_stereo, set_tautomer_remove_bond_stereo, tautomer_remove_bond_stereo: bool;
    /// remove isotopic Hs from centers involved in tautomerism (defaults to True)
    "tautomerRemoveIsotopicHs", get_tautomer_remove_isotopic_hs, set_tautomer_remove_isotopic_hs, tautomer_remove_isotopic_hs: bool;
    /// call AssignStereochemistry on all generated tautomers (defaults to True)
    "tautomerReassignStereo", get_tautomer_reassign_stereo, set_tautomer_reassign_stereo, tautomer_reassign_stereo: bool;
    /// Whether LargestFragmentChooser should use atom count as main criterion before MW (defaults to True)
    "largestFragmentChooserUseAtomCount", get_lfc_use_atom_count, set_lfc_use_atom_count, largest_fragment_chooser_use_atom_count: bool;
    /// whether LargestFragmentChooser should only count heavy atoms (defaults to False)
    "largestFragmentChooserCountHeavyAtomsOnly", get_lfc_count_heavy_atoms_only, set_lfc_count_heavy_atoms_only, largest_fragment_chooser_count_heavy_atoms_only: bool;
}

#[pymethods]
impl PyCleanupParameters {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// Only the documented attributes can be set (upstream's `safeSetattr`).
    fn __setattr__(slf: &Bound<'_, Self>, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        crate::pyconv::set_declared_attribute(slf.as_any(), name, value)
    }
}

/// updates the cleanup parameters from the provided JSON string
#[pyfunction]
#[pyo3(name = "UpdateParamsFromJSON")]
fn update_params_from_json(
    mut params: PyRefMut<'_, PyCleanupParameters>,
    json: &str,
) -> PyResult<()> {
    update_cleanup_params_from_json(&mut params.params, json).py()
}

// ------------------------------------------------------- function plumbing

fn params_or_default(params: Option<&PyCleanupParameters>) -> CleanupParameters {
    params.map_or_else(
        || default_cleanup_parameters().clone(),
        |p| p.params.clone(),
    )
}

fn none_error() -> PyErr {
    PyValueError::new_err("Molecule is None")
}

fn incompatible(name: &str) -> PyErr {
    PyTypeError::new_err(format!("{name}(): incompatible function arguments"))
}

/// A standardization step: (molecule, parameters, skipStandardize).
type InPlaceFn = fn(&mut RWMol, &CleanupParameters, bool) -> rdoxide::Result<()>;

/// The copying form, `X(mol, params=None[, skipStandardize=False])`.
fn run_copy(
    py: Python<'_>,
    mol: Option<PyRef<'_, PyMol>>,
    params: Option<PyRef<'_, PyCleanupParameters>>,
    skip_standardize: bool,
    func: InPlaceFn,
) -> PyResult<PyMol> {
    let mol = mol.ok_or_else(none_error)?;
    let params = params_or_default(params.as_deref());
    let mut res = mol.mol.clone();
    py.detach(|| func(&mut res, &params, skip_standardize))
        .py()?;
    Ok(PyMol { mol: res })
}

/// `XInPlace(mol, params=None[, skipStandardize])` and
/// `XInPlace(mols, numThreads, params=None[, skipStandardize])`.
fn run_in_place(
    name: &str,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
    with_skip: bool,
    func: InPlaceFn,
) -> PyResult<()> {
    let py = args.py();
    let first = match args.get_item(0) {
        Ok(a) => Some(a),
        Err(_) => kwargs
            .and_then(|k| {
                k.get_item("mol")
                    .transpose()
                    .or_else(|| k.get_item("mols").transpose())
            })
            .transpose()?,
    }
    .ok_or_else(|| incompatible(name))?;
    let multi = !first.is_instance_of::<PyMol>();
    let mut names = vec![if multi { "mols" } else { "mol" }];
    if multi {
        names.push("numThreads");
    }
    names.push("params");
    if with_skip {
        names.push("skipStandardize");
    }
    // bind the arguments by position, then by keyword
    let mut bound: Vec<Option<Bound<'_, PyAny>>> = vec![None; names.len()];
    if args.len() > names.len() {
        return Err(incompatible(name));
    }
    for (i, a) in args.iter().enumerate() {
        bound[i] = Some(a);
    }
    if let Some(kwargs) = kwargs {
        for (k, v) in kwargs.iter() {
            let k: String = k.extract()?;
            let i = names
                .iter()
                .position(|n| *n == k)
                .ok_or_else(|| incompatible(name))?;
            if bound[i].is_some() {
                return Err(incompatible(name));
            }
            bound[i] = Some(v);
        }
    }
    let arg = |n: &str| bound[names.iter().position(|x| *x == n).unwrap_or(0)].as_ref();
    let params = match arg("params").filter(|p| !p.is_none()) {
        Some(p) => Some(
            p.cast::<PyCleanupParameters>()
                .map_err(|_| incompatible(name))?
                .borrow()
                .params
                .clone(),
        ),
        None => None,
    };
    let params = params.unwrap_or_else(|| default_cleanup_parameters().clone());
    let skip = match arg("skipStandardize") {
        Some(s) if with_skip => s.extract::<bool>()?,
        _ => false,
    };
    if !multi {
        let mol = first.cast::<PyMol>().map_err(|_| incompatible(name))?;
        let mut mol = mol.borrow_mut();
        let mol = &mut mol.mol;
        return py.detach(|| func(mol, &params, skip)).py();
    }
    let num_threads: i32 = arg("numThreads")
        .ok_or_else(|| incompatible(name))?
        .extract()?;
    let seq: Vec<Bound<'_, PyAny>> = first
        .try_iter()
        .map_err(|_| incompatible(name))?
        .collect::<PyResult<_>>()?;
    let mut mols = Vec::with_capacity(seq.len());
    for (i, m) in seq.iter().enumerate() {
        if seq[..i].iter().any(|o| o.is(m)) {
            return Err(PyValueError::new_err("duplicate molecule in input list"));
        }
        mols.push(
            m.cast::<PyMol>()
                .map_err(|_| incompatible(name))?
                .borrow_mut(),
        );
    }
    let mut refs: Vec<&mut RWMol> = mols.iter_mut().map(|m| &mut m.mol).collect();
    py.detach(|| {
        standardize_multiple_mols_in_place(&mut refs, num_threads, |m| func(m, &params, skip))
    })
    .py()
}

/// The copying function `X(mol, params=None[, skipStandardize=False])`.
macro_rules! copy_function {
    ($(#[$doc:meta])* $copy:ident = $name:literal, skip = true, $func:expr) => {
        $(#[$doc])*
        #[pyfunction]
        #[pyo3(name = $name, signature = (mol, params = None, skipStandardize = false))]
        #[allow(non_snake_case)]
        fn $copy(
            py: Python<'_>,
            mol: Option<PyRef<'_, PyMol>>,
            params: Option<PyRef<'_, PyCleanupParameters>>,
            skipStandardize: bool,
        ) -> PyResult<PyMol> {
            run_copy(py, mol, params, skipStandardize, $func)
        }
    };
    ($(#[$doc:meta])* $copy:ident = $name:literal, skip = false, $func:expr) => {
        $(#[$doc])*
        #[pyfunction]
        #[pyo3(name = $name, signature = (mol, params = None))]
        fn $copy(
            py: Python<'_>,
            mol: Option<PyRef<'_, PyMol>>,
            params: Option<PyRef<'_, PyCleanupParameters>>,
        ) -> PyResult<PyMol> {
            run_copy(py, mol, params, false, $func)
        }
    };
}

/// The in-place function with its single- and multi-molecule overloads.
macro_rules! in_place_function {
    ($(#[$doc:meta])* $in_place:ident = $name:literal, skip = $skip:literal, $func:expr) => {
        $(#[$doc])*
        #[pyfunction]
        #[pyo3(name = $name, signature = (*args, **kwargs))]
        fn $in_place(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
            run_in_place($name, args, kwargs, $skip, $func)
        }
    };
}

copy_function! {
    /// Standardizes a molecule
    cleanup_py = "Cleanup", skip = false,
    |m, p, _| ms::cleanup_in_place(m, p)
}
in_place_function! {
    /// Standardizes a molecule in place
    cleanup_in_place_py = "CleanupInPlace", skip = false,
    |m, p, _| ms::cleanup_in_place(m, p)
}
copy_function! {
    /// Returns the tautomer parent of a given molecule. The fragment parent is
    /// the standardized canonical tautomer of the molecule
    tautomer_parent_py = "TautomerParent", skip = true, ms::tautomer_parent_in_place
}
in_place_function! {
    /// Generates the tautomer parent in place
    tautomer_parent_in_place_py = "TautomerParentInPlace", skip = true, ms::tautomer_parent_in_place
}
copy_function! {
    /// Returns the largest fragment after doing a cleanup
    fragment_parent_py = "FragmentParent", skip = true, ms::fragment_parent_in_place
}
in_place_function! {
    /// Generates the largest fragment in place
    fragment_parent_in_place_py = "FragmentParentInPlace", skip = true, ms::fragment_parent_in_place
}
copy_function! {
    /// Returns the stereo parent of the molecule
    stereo_parent_py = "StereoParent", skip = true, ms::stereo_parent_in_place
}
in_place_function! {
    /// Generates the stereo parent in place
    stereo_parent_in_place_py = "StereoParentInPlace", skip = true, ms::stereo_parent_in_place
}
copy_function! {
    /// removes all isotopes specifications from the given molecule
    isotope_parent_py = "IsotopeParent", skip = true, ms::isotope_parent_in_place
}
in_place_function! {
    /// Generates the isotope parent in place
    isotope_parent_in_place_py = "IsotopeParentInPlace", skip = true, ms::isotope_parent_in_place
}
copy_function! {
    /// Returns the uncharged version of the largest fragment
    charge_parent_py = "ChargeParent", skip = true, ms::charge_parent_in_place
}
in_place_function! {
    /// Generates the charge parent in place
    charge_parent_in_place_py = "ChargeParentInPlace", skip = true, ms::charge_parent_in_place
}
copy_function! {
    /// Returns the super parent. The super parent is the fragment, charge,
    /// isotope, stereo, and tautomer parent of the molecule.
    super_parent_py = "SuperParent", skip = true, ms::super_parent_in_place
}
in_place_function! {
    /// Generates the super parent in place
    super_parent_in_place_py = "SuperParentInPlace", skip = true, ms::super_parent_in_place
}

/// Convenience function for standardizing a SMILES
#[pyfunction]
#[pyo3(name = "StandardizeSmiles")]
fn standardize_smiles_py(py: Python<'_>, smiles: &str) -> PyResult<String> {
    py.detach(|| ms::standardize_smiles(smiles)).py()
}

copy_function! {
    /// Returns the canonical tautomer for the molecule
    canonical_tautomer_py = "CanonicalTautomer", skip = false,
    |m, p, _| { *m = ms::canonical_tautomer(m, p)?; Ok(()) }
}

copy_function! {
    /// Removes fragments from the molecule
    remove_fragments_py = "RemoveFragments", skip = false,
    |m, p, _| ms::remove_fragments_in_place(m, p)
}
in_place_function! {
    /// Removes fragments from the molecule, modifies the input molecule
    remove_fragments_in_place_py = "RemoveFragmentsInPlace", skip = false,
    |m, p, _| ms::remove_fragments_in_place(m, p)
}

copy_function! {
    /// Applies a series of standard transformations to correct functional
    /// groups and recombine charges
    normalize_py = "Normalize", skip = false,
    |m, p, _| { *m = ms::normalize(m, p)?; Ok(()) }
}
in_place_function! {
    /// Applies a series of standard transformations to correct functional
    /// groups and recombine charges, modifies the input molecule
    normalize_in_place_py = "NormalizeInPlace", skip = false,
    |m, p, _| ms::normalize_in_place(m, p)
}

copy_function! {
    /// Ensures the strongest acid groups are charged first
    reionize_py = "Reionize", skip = false,
    |m, p, _| { *m = ms::reionize(m, p)?; Ok(()) }
}
in_place_function! {
    /// Ensures the strongest acid groups are charged first, modifies the input molecule
    reionize_in_place_py = "ReionizeInPlace", skip = false,
    |m, p, _| ms::reionize_in_place(m, p)
}

fn add_standardize_functions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(cleanup_py, m)?)?;
    m.add_function(wrap_pyfunction!(cleanup_in_place_py, m)?)?;
    m.add_function(wrap_pyfunction!(standardize_smiles_py, m)?)?;
    m.add_function(wrap_pyfunction!(tautomer_parent_py, m)?)?;
    m.add_function(wrap_pyfunction!(tautomer_parent_in_place_py, m)?)?;
    m.add_function(wrap_pyfunction!(fragment_parent_py, m)?)?;
    m.add_function(wrap_pyfunction!(fragment_parent_in_place_py, m)?)?;
    m.add_function(wrap_pyfunction!(stereo_parent_py, m)?)?;
    m.add_function(wrap_pyfunction!(stereo_parent_in_place_py, m)?)?;
    m.add_function(wrap_pyfunction!(isotope_parent_py, m)?)?;
    m.add_function(wrap_pyfunction!(isotope_parent_in_place_py, m)?)?;
    m.add_function(wrap_pyfunction!(charge_parent_py, m)?)?;
    m.add_function(wrap_pyfunction!(charge_parent_in_place_py, m)?)?;
    m.add_function(wrap_pyfunction!(super_parent_py, m)?)?;
    m.add_function(wrap_pyfunction!(super_parent_in_place_py, m)?)?;
    m.add_function(wrap_pyfunction!(canonical_tautomer_py, m)?)?;
    m.add_function(wrap_pyfunction!(reionize_py, m)?)?;
    m.add_function(wrap_pyfunction!(reionize_in_place_py, m)?)?;
    m.add_function(wrap_pyfunction!(normalize_py, m)?)?;
    m.add_function(wrap_pyfunction!(normalize_in_place_py, m)?)?;
    m.add_function(wrap_pyfunction!(remove_fragments_py, m)?)?;
    m.add_function(wrap_pyfunction!(remove_fragments_in_place_py, m)?)?;
    Ok(())
}

// ------------------------------------------------------------------ Fragment

#[pyclass(
    name = "FragmentRemover",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize"
)]
#[derive(Debug)]
pub struct PyFragmentRemover {
    inner: FragmentRemover,
}

#[pymethods]
impl PyFragmentRemover {
    #[new]
    #[pyo3(signature = (fragmentFilename = "", leave_last = true, skip_if_all_match = false))]
    #[allow(non_snake_case)]
    fn new(fragmentFilename: &str, leave_last: bool, skip_if_all_match: bool) -> PyResult<Self> {
        let inner = if fragmentFilename.is_empty() {
            let mut remover = FragmentRemover::default();
            remover.leave_last = leave_last;
            remover.skip_if_all_match = skip_if_all_match;
            remover
        } else {
            FragmentRemover::from_file(fragmentFilename, leave_last, skip_if_all_match).py()?
        };
        Ok(PyFragmentRemover { inner })
    }

    fn remove(&self, py: Python<'_>, mol: PyRef<'_, PyMol>) -> PyResult<PyMol> {
        let mol = &mol.mol;
        let res = py.detach(|| self.inner.remove(mol)).py()?;
        Ok(PyMol { mol: res })
    }

    /// modifies the molecule in place
    #[pyo3(name = "removeInPlace")]
    fn remove_in_place(&self, py: Python<'_>, mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
        let mol = &mut mol.mol;
        py.detach(|| self.inner.remove_in_place(mol)).py()
    }
}

/// creates a FragmentRemover from a string containing parameter data
#[pyfunction]
#[pyo3(
    name = "FragmentRemoverFromData",
    signature = (fragmentData, leave_last = true, skip_if_all_match = false)
)]
#[allow(non_snake_case)]
fn fragment_remover_from_data(
    fragmentData: &str,
    leave_last: bool,
    skip_if_all_match: bool,
) -> PyResult<PyFragmentRemover> {
    let inner = FragmentRemover::from_text(fragmentData, leave_last, skip_if_all_match).py()?;
    Ok(PyFragmentRemover { inner })
}

#[pyclass(
    name = "LargestFragmentChooser",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize"
)]
#[derive(Debug)]
pub struct PyLargestFragmentChooser {
    inner: LargestFragmentChooser,
}

#[pymethods]
impl PyLargestFragmentChooser {
    /// `LargestFragmentChooser(preferOrganic=False)` or
    /// `LargestFragmentChooser(params)`.
    #[new]
    #[pyo3(signature = (preferOrganic = None, params = None))]
    #[allow(non_snake_case)]
    fn new(
        preferOrganic: Option<&Bound<'_, PyAny>>,
        params: Option<PyRef<'_, PyCleanupParameters>>,
    ) -> PyResult<Self> {
        let inner = match (preferOrganic, params) {
            (None, None) => LargestFragmentChooser::new(false),
            (None, Some(p)) => LargestFragmentChooser::from_params(&p.params),
            (Some(arg), None) => match arg.cast::<PyCleanupParameters>() {
                Ok(p) => LargestFragmentChooser::from_params(&p.borrow().params),
                Err(_) => LargestFragmentChooser::new(arg.extract::<bool>()?),
            },
            (Some(_), Some(_)) => {
                return Err(PyTypeError::new_err(
                    "__init__(): incompatible function arguments",
                ));
            }
        };
        Ok(PyLargestFragmentChooser { inner })
    }

    fn choose(&self, py: Python<'_>, mol: PyRef<'_, PyMol>) -> PyResult<PyMol> {
        let mol = &mol.mol;
        let res = py.detach(|| self.inner.choose(mol)).py()?;
        Ok(PyMol { mol: res })
    }

    /// modifies the molecule in place
    #[pyo3(name = "chooseInPlace")]
    fn choose_in_place(&self, py: Python<'_>, mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
        let mol = &mut mol.mol;
        py.detach(|| self.inner.choose_in_place(mol)).py()
    }
}

// ----------------------------------------------------------------- Normalize

#[pyclass(
    name = "Normalizer",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize"
)]
#[derive(Debug)]
pub struct PyNormalizer {
    inner: Normalizer,
}

#[pymethods]
impl PyNormalizer {
    /// `Normalizer()` or `Normalizer(normalizeFilename, maxRestarts)`.
    #[new]
    #[pyo3(signature = (normalizeFilename = None, maxRestarts = None))]
    #[allow(non_snake_case)]
    fn new(normalizeFilename: Option<String>, maxRestarts: Option<u32>) -> PyResult<Self> {
        let inner = match (normalizeFilename, maxRestarts) {
            (None, None) => Normalizer::default(),
            (Some(f), Some(n)) => Normalizer::from_file(&f, n).py()?,
            _ => return Err(incompatible("__init__")),
        };
        Ok(PyNormalizer { inner })
    }

    fn normalize(&self, py: Python<'_>, mol: PyRef<'_, PyMol>) -> PyResult<PyMol> {
        let mol = &mol.mol;
        let res = py.detach(|| self.inner.normalize(mol)).py()?;
        Ok(PyMol { mol: res })
    }

    /// modifies the input molecule
    #[pyo3(name = "normalizeInPlace")]
    fn normalize_in_place(&self, py: Python<'_>, mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
        let mol = &mut mol.mol;
        py.detach(|| self.inner.normalize_in_place(mol)).py()
    }
}

/// creates a Normalizer from a string containing normalization SMARTS
#[pyfunction]
#[pyo3(name = "NormalizerFromData")]
#[allow(non_snake_case)]
fn normalizer_from_data(
    paramData: &str,
    params: PyRef<'_, PyCleanupParameters>,
) -> PyResult<PyNormalizer> {
    let inner = Normalizer::from_text(paramData, params.params.max_restarts as u32).py()?;
    Ok(PyNormalizer { inner })
}

/// creates a Normalizer from CleanupParameters
#[pyfunction]
#[pyo3(name = "NormalizerFromParams")]
fn normalizer_from_params_py(params: PyRef<'_, PyCleanupParameters>) -> PyResult<PyNormalizer> {
    let inner = normalizer_from_params(&params.params).py()?;
    Ok(PyNormalizer { inner })
}

// -------------------------------------------------------------------- Charge

#[pyclass(
    name = "ChargeCorrection",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
pub struct PyChargeCorrection {
    inner: ChargeCorrection,
}

#[pymethods]
impl PyChargeCorrection {
    #[new]
    fn new(name: String, smarts: String, charge: i32) -> Self {
        PyChargeCorrection {
            inner: ChargeCorrection {
                name,
                smarts,
                charge,
            },
        }
    }

    #[getter(Name)]
    fn get_name(&self) -> String {
        self.inner.name.clone()
    }
    #[setter(Name)]
    fn set_name(&mut self, v: String) {
        self.inner.name = v;
    }
    #[getter(Smarts)]
    fn get_smarts(&self) -> String {
        self.inner.smarts.clone()
    }
    #[setter(Smarts)]
    fn set_smarts(&mut self, v: String) {
        self.inner.smarts = v;
    }
    #[getter(Charge)]
    fn get_charge(&self) -> i32 {
        self.inner.charge
    }
    #[setter(Charge)]
    fn set_charge(&mut self, v: i32) {
        self.inner.charge = v;
    }
}

fn extract_ccs(ccs: &Bound<'_, PyAny>) -> PyResult<Vec<ChargeCorrection>> {
    ccs.try_iter()?
        .map(|cc| Ok(cc?.cast::<PyChargeCorrection>()?.borrow().inner.clone()))
        .collect()
}

#[pyfunction]
#[pyo3(name = "CHARGE_CORRECTIONS")]
fn charge_corrections_py() -> Vec<PyChargeCorrection> {
    charge_corrections()
        .into_iter()
        .map(|inner| PyChargeCorrection { inner })
        .collect()
}

#[pyclass(
    name = "Reionizer",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize"
)]
#[derive(Debug)]
pub struct PyReionizer {
    inner: Reionizer,
}

#[pymethods]
impl PyReionizer {
    /// `Reionizer()`, `Reionizer(acidbaseFile)` or `Reionizer(acidbaseFile, ccs)`.
    #[new]
    #[pyo3(signature = (acidbaseFile = None, ccs = None))]
    #[allow(non_snake_case)]
    fn new(acidbaseFile: Option<String>, ccs: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let inner = match (acidbaseFile, ccs) {
            (None, None) => Reionizer::default(),
            (Some(f), ccs) => {
                let ccs = ccs.map(extract_ccs).transpose()?;
                Reionizer::from_file(&f, ccs).py()?
            }
            (None, Some(_)) => return Err(incompatible("__init__")),
        };
        Ok(PyReionizer { inner })
    }

    fn reionize(&self, py: Python<'_>, mol: PyRef<'_, PyMol>) -> PyResult<PyMol> {
        let mol = &mol.mol;
        let res = py.detach(|| self.inner.reionize(mol)).py()?;
        Ok(PyMol { mol: res })
    }

    /// modifies the input molecule
    #[pyo3(name = "reionizeInPlace")]
    fn reionize_in_place(&self, py: Python<'_>, mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
        let mol = &mut mol.mol;
        py.detach(|| self.inner.reionize_in_place(mol)).py()
    }
}

/// creates a reionizer from a string containing parameter data and a list of charge corrections
#[pyfunction]
#[pyo3(name = "ReionizerFromData", signature = (paramData, chargeCorrections = None))]
#[allow(non_snake_case)]
fn reionizer_from_data(
    paramData: &str,
    chargeCorrections: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyReionizer> {
    let ccs = match chargeCorrections {
        Some(c) => extract_ccs(c)?,
        None => Vec::new(),
    };
    let inner = Reionizer::from_text(paramData, Some(ccs)).py()?;
    Ok(PyReionizer { inner })
}

#[pyclass(
    name = "Uncharger",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize"
)]
#[derive(Debug)]
pub struct PyUncharger {
    inner: Uncharger,
}

#[pymethods]
impl PyUncharger {
    #[new]
    #[pyo3(signature = (canonicalOrder = true, force = false, protonationOnly = false))]
    #[allow(non_snake_case)]
    fn new(canonicalOrder: bool, force: bool, protonationOnly: bool) -> Self {
        PyUncharger {
            inner: Uncharger::new(canonicalOrder, force, protonationOnly),
        }
    }

    fn uncharge(&self, py: Python<'_>, mol: PyRef<'_, PyMol>) -> PyResult<PyMol> {
        let mol = &mol.mol;
        let res = py.detach(|| self.inner.uncharge(mol)).py()?;
        Ok(PyMol { mol: res })
    }

    /// modifies the input molecule
    #[pyo3(name = "unchargeInPlace")]
    fn uncharge_in_place(&self, py: Python<'_>, mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
        let mol = &mut mol.mol;
        py.detach(|| self.inner.uncharge_in_place(mol)).py()
    }
}

// --------------------------------------------------------------------- Metal

/// Metal Disconnector Options
#[pyclass(
    name = "MetalDisconnectorOptions",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
pub struct PyMetalDisconnectorOptions {
    inner: MetalDisconnectorOptions,
}

#[pymethods]
impl PyMetalDisconnectorOptions {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// Whether to split Grignard-type complexes. Default false.
    #[getter(splitGrignards)]
    fn get_split_grignards(&self) -> bool {
        self.inner.split_grignards
    }
    #[setter(splitGrignards)]
    fn set_split_grignards(&mut self, v: bool) {
        self.inner.split_grignards = v;
    }
    /// Whether to split metal-aromatic C bonds.  Default false.
    #[getter(splitAromaticC)]
    fn get_split_aromatic_c(&self) -> bool {
        self.inner.split_aromatic_c
    }
    #[setter(splitAromaticC)]
    fn set_split_aromatic_c(&mut self, v: bool) {
        self.inner.split_aromatic_c = v;
    }
    /// Whether to adjust charges on ligand atoms.  Default true.
    #[getter(adjustCharges)]
    fn get_adjust_charges(&self) -> bool {
        self.inner.adjust_charges
    }
    #[setter(adjustCharges)]
    fn set_adjust_charges(&mut self, v: bool) {
        self.inner.adjust_charges = v;
    }
    /// Whether to remove the dummy atoms representing haptic bonds.  Such
    /// dummies are bonded to the metal with a bond that has the
    /// MolFileBondEndPts prop set.  Default false.
    #[getter(removeHapticDummies)]
    fn get_remove_haptic_dummies(&self) -> bool {
        self.inner.remove_haptic_dummies
    }
    #[setter(removeHapticDummies)]
    fn set_remove_haptic_dummies(&mut self, v: bool) {
        self.inner.remove_haptic_dummies = v;
    }
}

/// a class to disconnect metals that are defined as covalently bonded to non-metals
#[pyclass(
    name = "MetalDisconnector",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize"
)]
#[derive(Debug)]
pub struct PyMetalDisconnector {
    inner: MetalDisconnector,
}

fn query_smarts(mol: &rdoxide::graphmol::ROMol) -> PyResult<String> {
    let mut mol = mol.clone();
    mol_to_smarts(&mut mol, &SmilesWriteParams::default()).py()
}

#[pymethods]
impl PyMetalDisconnector {
    #[new]
    #[pyo3(signature = (options = None))]
    fn new(options: Option<PyRef<'_, PyMetalDisconnectorOptions>>) -> Self {
        let options = options.map(|o| o.inner).unwrap_or_default();
        PyMetalDisconnector {
            inner: MetalDisconnector::new(options),
        }
    }

    /// SMARTS defining the metals to disconnect if attached to Nitrogen, Oxygen or Fluorine
    #[getter(MetalNof)]
    fn metal_nof(&self) -> PyResult<String> {
        query_smarts(self.inner.metal_nof())
    }

    /// SMARTS defining the metals to disconnect other inorganic elements
    #[getter(MetalNon)]
    fn metal_non(&self) -> PyResult<String> {
        query_smarts(self.inner.metal_non())
    }

    /// Set the query molecule defining the metals to disconnect from other inorganic elements.
    #[pyo3(name = "SetMetalNon")]
    fn set_metal_non(&mut self, mol: PyRef<'_, PyMol>) {
        self.inner.set_metal_non(mol.mol.clone());
    }

    /// Set the query molecule defining the metals to disconnect if attached to Nitrogen, Oxygen or Fluorine.
    #[pyo3(name = "SetMetalNof")]
    fn set_metal_nof(&mut self, mol: PyRef<'_, PyMol>) {
        self.inner.set_metal_nof(mol.mol.clone());
    }

    /// performs the disconnection
    #[pyo3(name = "Disconnect")]
    fn disconnect(&self, py: Python<'_>, mol: PyRef<'_, PyMol>) -> PyResult<PyMol> {
        let mol = &mol.mol;
        let res = py.detach(|| self.inner.disconnect(mol)).py()?;
        Ok(PyMol { mol: res })
    }

    /// performs the disconnection, modifies the input molecule
    #[pyo3(name = "DisconnectInPlace")]
    fn disconnect_in_place(&self, py: Python<'_>, mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
        let mol = &mut mol.mol;
        py.detach(|| self.inner.disconnect_in_place(mol)).py()
    }
}

/// Returns the molecule disconnected using the organometallics rules.
#[pyfunction]
#[pyo3(name = "DisconnectOrganometallics", signature = (mol, params = None))]
fn disconnect_organometallics_py(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    params: Option<PyRef<'_, PyMetalDisconnectorOptions>>,
) -> PyResult<PyMol> {
    let mdo = params.map(|p| p.inner);
    let mol = &mol.mol;
    let res = py
        .detach(|| ms::disconnect_organometallics(mol, mdo))
        .py()?;
    Ok(PyMol { mol: res })
}

/// Disconnects the molecule using the organometallics rules, modifies the input molecule
#[pyfunction]
#[pyo3(name = "DisconnectOrganometallicsInPlace", signature = (mol, params = None))]
fn disconnect_organometallics_in_place_py(
    py: Python<'_>,
    mol: Option<PyRefMut<'_, PyMol>>,
    params: Option<PyRef<'_, PyMetalDisconnectorOptions>>,
) -> PyResult<()> {
    let mut mol = mol.ok_or_else(none_error)?;
    let mdo = params.map(|p| p.inner);
    let mol = &mut mol.mol;
    py.detach(|| ms::disconnect_organometallics_in_place(mol, mdo))
        .py()
}

// ------------------------------------------------------------------ Validate

/// A member of a `MolVSValidation` built from Python objects.
#[derive(Debug)]
enum Member {
    Native(ValidationMethod),
    Python(Py<PyAny>),
}

#[pyclass(
    name = "ValidationMethod",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    subclass
)]
#[derive(Debug, Default)]
pub struct PyValidationMethod {
    inner: Option<ValidationMethod>,
    /// composite members when some of them are implemented in Python
    members: Option<Vec<Member>>,
}

impl PyValidationMethod {
    fn native(inner: ValidationMethod) -> Self {
        PyValidationMethod {
            inner: Some(inner),
            members: None,
        }
    }
}

#[pymethods]
impl PyValidationMethod {
    #[new]
    #[pyo3(signature = (*_args, **_kwargs))]
    fn new(_args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> Self {
        Self::default()
    }

    #[pyo3(signature = (mol, reportAllFailures = false))]
    #[allow(non_snake_case)]
    fn validate(
        &self,
        py: Python<'_>,
        mol: &Bound<'_, PyMol>,
        reportAllFailures: bool,
    ) -> PyResult<Vec<String>> {
        // validations may store ring information on the molecule, as
        // upstream's Layout2DValidation does
        if let Some(inner) = &self.inner {
            let mut m = mol.borrow_mut();
            let m = &mut m.mol;
            return py
                .detach(|| inner.validate_in_place(m, reportAllFailures))
                .py();
        }
        let Some(members) = &self.members else {
            return Err(PyRuntimeError::new_err(
                "ValidationMethod.validate(): tried to call a pure virtual function",
            ));
        };
        let mut errors = Vec::new();
        for m in members {
            let partial: Vec<String> = match m {
                Member::Native(v) => v
                    .validate_in_place(&mut mol.borrow_mut().mol, reportAllFailures)
                    .py()?,
                Member::Python(obj) => obj
                    .bind(py)
                    .call_method1("validate", (mol, reportAllFailures))?
                    .extract()?,
            };
            if !partial.is_empty() {
                errors.extend(partial);
                if !reportAllFailures {
                    break;
                }
            }
        }
        Ok(errors)
    }
}

fn sub<T: pyo3::PyClass<BaseType = PyValidationMethod>>(
    value: T,
    base: PyValidationMethod,
) -> PyClassInitializer<T> {
    PyClassInitializer::from(base).add_subclass(value)
}

#[rustfmt::skip]
macro_rules! simple_validation {
    ($name:ident, $py:literal, $value:expr) => {
        #[pyclass(
            name = $py,
            module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
            extends = PyValidationMethod
        )]
        #[derive(Debug)]
        pub struct $name;

        #[pymethods]
        impl $name {
            #[new]
            fn new() -> PyClassInitializer<Self> {
                sub($name, PyValidationMethod::native($value))
            }
        }
    };
}

simple_validation!(
    PyNoAtomValidation,
    "NoAtomValidation",
    ValidationMethod::NoAtom
);
simple_validation!(
    PyFragmentValidation,
    "FragmentValidation",
    ValidationMethod::Fragment
);
simple_validation!(
    PyNeutralValidation,
    "NeutralValidation",
    ValidationMethod::Neutral
);
simple_validation!(
    PyDisallowedRadicalValidation,
    "DisallowedRadicalValidation",
    ValidationMethod::DisallowedRadical
);
simple_validation!(
    PyStereoValidation,
    "StereoValidation",
    ValidationMethod::Stereo
);

/// Read/write attributes of a validation stored in the base class.
macro_rules! validation_attributes {
    ($name:ident, $variant:ident { $($py:literal, $get:ident, $set:ident, $field:ident: $ty:ty;)* }) => {
        #[pymethods]
        impl $name {
            $(
                #[getter($py)]
                fn $get(slf: PyRef<'_, Self>) -> $ty {
                    match &slf.as_super().inner {
                        Some(ValidationMethod::$variant(v)) => v.$field,
                        _ => Default::default(),
                    }
                }
                #[setter($py)]
                fn $set(mut slf: PyRefMut<'_, Self>, value: $ty) {
                    if let Some(ValidationMethod::$variant(v)) = &mut slf.as_super().inner {
                        v.$field = value;
                    }
                }
            )*
        }
    };
}

#[pyclass(
    name = "RDKitValidation",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    extends = PyValidationMethod
)]
#[derive(Debug)]
pub struct PyRDKitValidation;

#[pymethods]
impl PyRDKitValidation {
    #[new]
    #[pyo3(signature = (allowEmptyMolecules = false))]
    #[allow(non_snake_case)]
    fn new(allowEmptyMolecules: bool) -> PyClassInitializer<Self> {
        sub(
            PyRDKitValidation,
            PyValidationMethod::native(ValidationMethod::RDKit {
                allow_empty_molecules: allowEmptyMolecules,
            }),
        )
    }

    #[getter(allowEmptyMolecules)]
    fn get_allow_empty(slf: PyRef<'_, Self>) -> bool {
        matches!(
            slf.as_super().inner,
            Some(ValidationMethod::RDKit {
                allow_empty_molecules: true
            })
        )
    }
    #[setter(allowEmptyMolecules)]
    fn set_allow_empty(mut slf: PyRefMut<'_, Self>, value: bool) {
        slf.as_super().inner = Some(ValidationMethod::RDKit {
            allow_empty_molecules: value,
        });
    }
}

#[pyclass(
    name = "IsotopeValidation",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    extends = PyValidationMethod
)]
#[derive(Debug)]
pub struct PyIsotopeValidation;

#[pymethods]
impl PyIsotopeValidation {
    #[new]
    #[pyo3(signature = (strict = false))]
    fn new(strict: bool) -> PyClassInitializer<Self> {
        sub(
            PyIsotopeValidation,
            PyValidationMethod::native(ValidationMethod::Isotope { strict }),
        )
    }

    #[getter(strict)]
    fn get_strict(slf: PyRef<'_, Self>) -> bool {
        matches!(
            slf.as_super().inner,
            Some(ValidationMethod::Isotope { strict: true })
        )
    }
    #[setter(strict)]
    fn set_strict(mut slf: PyRefMut<'_, Self>, value: bool) {
        slf.as_super().inner = Some(ValidationMethod::Isotope { strict: value });
    }
}

#[pyclass(
    name = "Is2DValidation",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    extends = PyValidationMethod
)]
#[derive(Debug)]
pub struct PyIs2DValidation;

#[pymethods]
impl PyIs2DValidation {
    #[new]
    #[pyo3(signature = (threshold = 1e-3))]
    fn new(threshold: f64) -> PyClassInitializer<Self> {
        sub(
            PyIs2DValidation,
            PyValidationMethod::native(ValidationMethod::Is2D { threshold }),
        )
    }

    #[getter(threshold)]
    fn get_threshold(slf: PyRef<'_, Self>) -> f64 {
        match slf.as_super().inner {
            Some(ValidationMethod::Is2D { threshold }) => threshold,
            _ => 0.0,
        }
    }
    #[setter(threshold)]
    fn set_threshold(mut slf: PyRefMut<'_, Self>, value: f64) {
        slf.as_super().inner = Some(ValidationMethod::Is2D { threshold: value });
    }
}

#[pyclass(
    name = "FeaturesValidation",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    extends = PyValidationMethod
)]
#[derive(Debug)]
pub struct PyFeaturesValidation;

#[pymethods]
impl PyFeaturesValidation {
    #[new]
    #[pyo3(signature = (
        allowEnhancedStereo = false, allowAromaticBondType = false, allowDativeBondType = false,
        allowQueries = false, allowDummies = false, allowAtomAliases = false
    ))]
    #[allow(non_snake_case)]
    fn new(
        allowEnhancedStereo: bool,
        allowAromaticBondType: bool,
        allowDativeBondType: bool,
        allowQueries: bool,
        allowDummies: bool,
        allowAtomAliases: bool,
    ) -> PyClassInitializer<Self> {
        sub(
            PyFeaturesValidation,
            PyValidationMethod::native(ValidationMethod::Features(FeaturesValidation {
                allow_enhanced_stereo: allowEnhancedStereo,
                allow_aromatic_bond_type: allowAromaticBondType,
                allow_dative_bond_type: allowDativeBondType,
                allow_queries: allowQueries,
                allow_dummies: allowDummies,
                allow_atom_aliases: allowAtomAliases,
            })),
        )
    }
}

validation_attributes!(PyFeaturesValidation, Features {
    "allowEnhancedStereo", get_allow_enhanced_stereo, set_allow_enhanced_stereo, allow_enhanced_stereo: bool;
    "allowAromaticBondType", get_allow_aromatic_bond_type, set_allow_aromatic_bond_type, allow_aromatic_bond_type: bool;
    "allowDativeBondType", get_allow_dative_bond_type, set_allow_dative_bond_type, allow_dative_bond_type: bool;
    "allowQueries", get_allow_queries, set_allow_queries, allow_queries: bool;
    "allowDummies", get_allow_dummies, set_allow_dummies, allow_dummies: bool;
    "allowAtomAliases", get_allow_atom_aliases, set_allow_atom_aliases, allow_atom_aliases: bool;
});

#[pyclass(
    name = "Layout2DValidation",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    extends = PyValidationMethod
)]
#[derive(Debug)]
pub struct PyLayout2DValidation;

#[pymethods]
impl PyLayout2DValidation {
    /// The `minMedianBondLength` default is `false` (0.0) in the upstream
    /// wrapper.
    #[new]
    #[pyo3(signature = (
        clashLimit = 0.15, bondLengthLimit = 25.0, allowLongBondsInRings = true,
        allowAtomBondClashExemption = true, minMedianBondLength = 0.0
    ))]
    #[allow(non_snake_case)]
    fn new(
        clashLimit: f64,
        bondLengthLimit: f64,
        allowLongBondsInRings: bool,
        allowAtomBondClashExemption: bool,
        minMedianBondLength: f64,
    ) -> PyClassInitializer<Self> {
        sub(
            PyLayout2DValidation,
            PyValidationMethod::native(ValidationMethod::Layout2D(Layout2DValidation {
                clash_limit: clashLimit,
                bond_length_limit: bondLengthLimit,
                allow_long_bonds_in_rings: allowLongBondsInRings,
                allow_atom_bond_clash_exemption: allowAtomBondClashExemption,
                min_median_bond_length: minMedianBondLength,
            })),
        )
    }
}

validation_attributes!(PyLayout2DValidation, Layout2D {
    "clashLimit", get_clash_limit, set_clash_limit, clash_limit: f64;
    "bondLengthLimit", get_bond_length_limit, set_bond_length_limit, bond_length_limit: f64;
    "allowLongBondsInRings", get_allow_long_bonds, set_allow_long_bonds, allow_long_bonds_in_rings: bool;
    "allowAtomBondClashExemption", get_allow_clash_exemption, set_allow_clash_exemption, allow_atom_bond_clash_exemption: bool;
    "minMedianBondLength", get_min_median_bond_length, set_min_median_bond_length, min_median_bond_length: f64;
});

fn extract_atoms(py: Python<'_>, atoms: &Bound<'_, PyAny>, what: &str) -> PyResult<Vec<Atom>> {
    let mut res = Vec::new();
    for a in atoms.try_iter()? {
        let a = a?;
        let a = a.cast::<PyAtom>().map_err(|_| incompatible("__init__"))?;
        res.push(a.borrow().with(py, |at| at.copy_detached())?);
    }
    if res.is_empty() {
        return Err(PyValueError::new_err(format!("{what} must be non-empty")));
    }
    Ok(res)
}

#[pyclass(
    name = "AllowedAtomsValidation",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    extends = PyValidationMethod
)]
#[derive(Debug)]
pub struct PyAllowedAtomsValidation;

#[pymethods]
impl PyAllowedAtomsValidation {
    #[new]
    fn new(py: Python<'_>, atoms: &Bound<'_, PyAny>) -> PyResult<PyClassInitializer<Self>> {
        let atoms = extract_atoms(py, atoms, "allowedAtoms argument")?;
        Ok(sub(
            PyAllowedAtomsValidation,
            PyValidationMethod::native(ValidationMethod::AllowedAtoms(atoms)),
        ))
    }
}

#[pyclass(
    name = "DisallowedAtomsValidation",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    extends = PyValidationMethod
)]
#[derive(Debug)]
pub struct PyDisallowedAtomsValidation;

#[pymethods]
impl PyDisallowedAtomsValidation {
    #[new]
    fn new(py: Python<'_>, atoms: &Bound<'_, PyAny>) -> PyResult<PyClassInitializer<Self>> {
        let atoms = extract_atoms(py, atoms, "disallowedAtoms")?;
        Ok(sub(
            PyDisallowedAtomsValidation,
            PyValidationMethod::native(ValidationMethod::DisallowedAtoms(atoms)),
        ))
    }
}

#[pyclass(
    name = "MolVSValidation",
    module = "rdoxide.Chem.MolStandardize.rdMolStandardize",
    extends = PyValidationMethod
)]
#[derive(Debug)]
pub struct PyMolVSValidation;

#[pymethods]
impl PyMolVSValidation {
    /// `MolVSValidation()` or `MolVSValidation(validations)`; the
    /// validations are copied.
    #[new]
    #[pyo3(signature = (validations = None))]
    fn new(
        py: Python<'_>,
        validations: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<PyClassInitializer<Self>> {
        let Some(validations) = validations else {
            return Ok(sub(
                PyMolVSValidation,
                PyValidationMethod::native(mol_vs_validation()),
            ));
        };
        let mut members = Vec::new();
        for v in validations.try_iter()? {
            let v = v?;
            let base = v
                .cast::<PyValidationMethod>()
                .map_err(|_| incompatible("__init__"))?;
            let base = base.borrow();
            // a Python subclass overriding validate() is called back
            let overridden = v
                .get_type()
                .getattr("validate")?
                .ne(py.get_type::<PyValidationMethod>().getattr("validate")?)?;
            match (&base.inner, overridden) {
                (Some(inner), false) => members.push(Member::Native(inner.clone())),
                _ => members.push(Member::Python(v.clone().unbind())),
            }
        }
        if members.is_empty() {
            return Err(PyValueError::new_err(
                "validations argument must be non-empty",
            ));
        }
        if members.iter().all(|m| matches!(m, Member::Native(_))) {
            let natives = members
                .into_iter()
                .filter_map(|m| match m {
                    Member::Native(v) => Some(v),
                    Member::Python(_) => None,
                })
                .collect();
            return Ok(sub(
                PyMolVSValidation,
                PyValidationMethod::native(ValidationMethod::Composite(natives)),
            ));
        }
        Ok(sub(
            PyMolVSValidation,
            PyValidationMethod {
                inner: None,
                members: Some(members),
            },
        ))
    }
}

#[pyfunction]
#[pyo3(name = "ValidateSmiles")]
fn validate_smiles_py(py: Python<'_>, mol: &str) -> PyResult<Vec<String>> {
    py.detach(|| validate_smiles(mol)).py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyCleanupParameters>()?;
    m.add_function(wrap_pyfunction!(update_params_from_json, m)?)?;
    add_standardize_functions(m)?;
    m.add_class::<PyFragmentRemover>()?;
    m.add_function(wrap_pyfunction!(fragment_remover_from_data, m)?)?;
    m.add_class::<PyLargestFragmentChooser>()?;
    m.add_class::<PyNormalizer>()?;
    m.add_class::<PyChargeCorrection>()?;
    m.add_function(wrap_pyfunction!(charge_corrections_py, m)?)?;
    m.add_class::<PyReionizer>()?;
    m.add_function(wrap_pyfunction!(reionizer_from_data, m)?)?;
    m.add_class::<PyUncharger>()?;
    m.add_class::<PyValidationMethod>()?;
    m.add_class::<PyRDKitValidation>()?;
    m.add_class::<PyNoAtomValidation>()?;
    m.add_class::<PyFragmentValidation>()?;
    m.add_class::<PyNeutralValidation>()?;
    m.add_class::<PyIsotopeValidation>()?;
    m.add_class::<PyMolVSValidation>()?;
    m.add_class::<PyAllowedAtomsValidation>()?;
    m.add_class::<PyDisallowedAtomsValidation>()?;
    m.add_class::<PyFeaturesValidation>()?;
    m.add_class::<PyDisallowedRadicalValidation>()?;
    m.add_class::<PyIs2DValidation>()?;
    m.add_class::<PyLayout2DValidation>()?;
    m.add_class::<PyStereoValidation>()?;
    m.add_function(wrap_pyfunction!(validate_smiles_py, m)?)?;
    crate::rdmolstandardize_tautomer::init(m)?;
    crate::rdmolstandardize_pipeline::init(m)?;
    m.add_class::<PyMetalDisconnectorOptions>()?;
    m.add_class::<PyMetalDisconnector>()?;
    m.add_function(wrap_pyfunction!(disconnect_organometallics_py, m)?)?;
    m.add_function(wrap_pyfunction!(disconnect_organometallics_in_place_py, m)?)?;
    m.add_function(wrap_pyfunction!(normalizer_from_data, m)?)?;
    m.add_function(wrap_pyfunction!(normalizer_from_params_py, m)?)?;
    Ok(())
}
