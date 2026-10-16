//! `rdoxide.Chem.rdRGroupDecomposition` (upstream
//! `RGroupDecomposition/Wrap/rdRGroupComposition.cpp`): decomposition of
//! molecules into a core and its R groups.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rdoxide::graphmol::ROMol;
use rdoxide::graphmol::rgroup_decomposition::{
    RGroupColumns, RGroupDecomposition, RGroupDecompositionParameters, RGroupRows,
    check_for_timeout, relabel_mapped_dummies, rgroup_core_alignment, rgroup_labelling,
    rgroup_labels, rgroup_matching, rgroup_score,
};
use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};

use crate::enums::add_enum_exported;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::pyconv::set_declared_attribute;
use crate::substruct::PySubstructMatchParameters;

const MODULE: &str = "rdoxide.Chem.rdRGroupDecomposition";

/// RGroupDecompositionParameters controls how the RGroupDecomposition sets labelling and matches structures
///   OPTIONS:
///     - RGroupCoreAlignment: can be one of RGroupCoreAlignment.None_ or RGroupCoreAlignment.MCS
///                            If set to MCS, cores labels are mapped to each other using their
///                            Maximum common substructure overlap.
///     - RGroupLabels: optionally set where the rgroup labels to use are encoded.
///                      RGroupLabels.IsotopeLabels - labels are stored on isotopes
///                      RGroupLabels.AtomMapLabels - labels are stored on atommaps
///                      RGroupLabels.MDLRGroupLabels - labels are stored on MDL R-groups
///                      RGroupLabels.DummyAtomLabels - labels are stored on dummy atoms
///                      RGroupLabels.AtomIndexLabels - use the atom index as the label
///                      RGroupLabels.RelabelDuplicateLabels - fix any duplicate labels
///                      RGroupLabels.AutoDetect - auto detect the label [default]
///        Note: in all cases, any rgroups found on unlabelled atoms will be automatically
///               labelled.
///     - RGroupLabelling: choose where the rlabels are stored on the decomposition
///                         RGroupLabelling.AtomMap - store rgroups as atom maps (for smiles)
///                         RGroupLabelling.Isotope - store rgroups on the isotope
///                         RGroupLabelling.MDLRGroup - store rgroups as mdl rgroups (for molblocks)
///                        default: AtomMap | MDLRGroup
///     - onlyMatchAtRGroups: only allow rgroup decomposition at the specified rgroups
///     - removeAllHydrogenRGroups: remove all user-defined rgroups that only have hydrogens
///     - removeAllHydrogenRGroupsAndLabels: remove all user-defined rgroups that only have hydrogens, and also remove the corresponding labels from the core
///     - removeHydrogensPostMatch: remove all hydrogens from the output molecules
///     - allowNonTerminalRGroups: allow labelled Rgroups of degree 2 or more
///     - doTautomers: match all tautomers of a core against each input structure
///     - doEnumeration: expand input cores into enumerated mol bundles
///     - allowMultipleRGroupsOnUnlabelled: permit more than one rgroup to be attached to an unlabelled core atom
///     - allowMultipleCoresInSameMol: permit a core to match more than once in the same molecule if the sets of matched atoms are not equal (default=False)
#[pyclass(
    name = "RGroupDecompositionParameters",
    module = "rdoxide.Chem.rdRGroupDecomposition",
    skip_from_py_object
)]
#[derive(Debug)]
pub struct PyRGroupDecompositionParameters {
    p: RGroupDecompositionParameters,
    substruct: Option<Py<PySubstructMatchParameters>>,
}

impl PyRGroupDecompositionParameters {
    /// The parameters, with the substructure parameters as edited in Python.
    fn snapshot(&self, py: Python<'_>) -> RGroupDecompositionParameters {
        let mut p = self.p.clone();
        if let Some(s) = &self.substruct {
            p.substructmatch_params = s.borrow(py).base_params();
        }
        p
    }
}

macro_rules! params_props {
    ($( ($py:literal, $get:ident, $set:ident, $field:ident, $ty:ty) ),* $(,)?) => {
        #[pymethods]
        #[allow(non_snake_case)]
        impl PyRGroupDecompositionParameters {
            $(
                #[getter($py)]
                fn $get(&self) -> $ty {
                    self.p.$field
                }
                #[setter($py)]
                fn $set(&mut self, value: $ty) {
                    self.p.$field = value;
                }
            )*
        }
    };
}

params_props!(
    ("labels", get_labels, set_labels, labels, u32),
    (
        "matchingStrategy",
        get_matching,
        set_matching,
        matching_strategy,
        u32
    ),
    ("scoreMethod", get_score, set_score, score_method, u32),
    (
        "rgroupLabelling",
        get_labelling,
        set_labelling,
        rgroup_labelling,
        u32
    ),
    ("alignment", get_alignment, set_alignment, alignment, u32),
    ("chunkSize", get_chunk, set_chunk, chunk_size, u32),
    (
        "onlyMatchAtRGroups",
        get_only,
        set_only,
        only_match_at_rgroups,
        bool
    ),
    (
        "removeAllHydrogenRGroups",
        get_rah,
        set_rah,
        remove_all_hydrogen_rgroups,
        bool
    ),
    (
        "removeHydrogensPostMatch",
        get_rhpm,
        set_rhpm,
        remove_hydrogens_post_match,
        bool
    ),
    ("timeout", get_timeout, set_timeout, timeout, f64),
    (
        "gaPopulationSize",
        get_gaps,
        set_gaps,
        ga_population_size,
        i32
    ),
    (
        "gaMaximumOperations",
        get_gamo,
        set_gamo,
        ga_maximum_operations,
        i32
    ),
    (
        "gaNumberOperationsWithoutImprovement",
        get_ganowi,
        set_ganowi,
        ga_number_operations_without_improvement,
        i32
    ),
    ("gaRandomSeed", get_gars, set_gars, ga_random_seed, i32),
    ("gaNumberRuns", get_ganr, set_ganr, ga_number_runs, i32),
    ("gaParallelRuns", get_gapr, set_gapr, ga_parallel_runs, bool),
    (
        "allowNonTerminalRGroups",
        get_antr,
        set_antr,
        allow_non_terminal_rgroups,
        bool
    ),
    (
        "removeAllHydrogenRGroupsAndLabels",
        get_rahl,
        set_rahl,
        remove_all_hydrogen_rgroups_and_labels,
        bool
    ),
    (
        "allowMultipleRGroupsOnUnlabelled",
        get_amru,
        set_amru,
        allow_multiple_rgroups_on_unlabelled,
        bool
    ),
    (
        "allowMultipleCoresInSameMol",
        get_amc,
        set_amc,
        allow_multiple_cores_in_same_mol,
        bool
    ),
    ("doTautomers", get_dt, set_dt, do_tautomers, bool),
    ("doEnumeration", get_de, set_de, do_enumeration, bool),
    (
        "includeTargetMolInResults",
        get_itm,
        set_itm,
        include_target_mol_in_results,
        bool
    ),
);

#[pymethods]
impl PyRGroupDecompositionParameters {
    /// Constructor, takes no arguments
    #[new]
    fn py_new() -> Self {
        PyRGroupDecompositionParameters {
            p: RGroupDecompositionParameters::default(),
            substruct: None,
        }
    }

    /// parameters controlling the substructure matching of the cores
    #[getter(substructMatchParams)]
    #[allow(non_snake_case)]
    fn substruct_match_params(
        &mut self,
        py: Python<'_>,
    ) -> PyResult<Py<PySubstructMatchParameters>> {
        if self.substruct.is_none() {
            self.substruct = Some(Py::new(
                py,
                PySubstructMatchParameters::from_core(&self.p.substructmatch_params),
            )?);
        }
        Ok(self.substruct.as_ref().expect("set above").clone_ref(py))
    }

    fn __setattr__(slf: &Bound<'_, Self>, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        set_declared_attribute(slf.as_any(), name, value)
    }
}

fn mol_from_py(obj: &Bound<'_, PyAny>) -> PyResult<ROMol> {
    if obj.is_none() {
        return Err(PyValueError::new_err("reaction called with None reactants"));
    }
    Ok(obj.extract::<PyRef<'_, PyMol>>()?.mol.clone())
}

fn mols_from_py(objs: &Bound<'_, PyAny>) -> PyResult<Vec<ROMol>> {
    objs.try_iter()?.map(|o| mol_from_py(&o?)).collect()
}

fn params_of(
    py: Python<'_>,
    options: Option<&Bound<'_, PyRGroupDecompositionParameters>>,
) -> RGroupDecompositionParameters {
    options.map_or_else(RGroupDecompositionParameters::default, |p| {
        p.borrow().snapshot(py)
    })
}

fn cores_from_py(cores: &Bound<'_, PyAny>) -> PyResult<Vec<ROMol>> {
    match cores.extract::<PyRef<'_, PyMol>>() {
        Ok(m) => Ok(vec![m.mol.clone()]),
        Err(_) => mols_from_py(cores),
    }
}

fn output_mol(py: Python<'_>, mol: ROMol, as_smiles: bool) -> PyResult<Py<PyAny>> {
    if as_smiles {
        let mut mol = mol;
        let smiles = mol_to_smiles(&mut mol, &SmilesWriteParams::default()).py()?;
        Ok(smiles.into_pyobject(py)?.into_any().unbind())
    } else {
        Ok(Py::new(py, PyMol { mol })?.into_any())
    }
}

fn rows_to_py<'py>(
    py: Python<'py>,
    rows: RGroupRows,
    as_smiles: bool,
) -> PyResult<Bound<'py, PyList>> {
    let result = PyList::empty(py);
    for row in rows {
        let dict = PyDict::new(py);
        for (label, mol) in row {
            dict.set_item(label, output_mol(py, mol, as_smiles)?)?;
        }
        result.append(dict)?;
    }
    Ok(result)
}

fn columns_to_py<'py>(
    py: Python<'py>,
    columns: RGroupColumns,
    as_smiles: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    for (label, mols) in columns {
        let col = PyList::empty(py);
        for mol in mols {
            col.append(output_mol(py, mol, as_smiles)?)?;
        }
        result.set_item(label, col)?;
    }
    Ok(result)
}

/// Decomposes molecules into a core and its R groups
#[pyclass(
    name = "RGroupDecomposition",
    module = "rdoxide.Chem.rdRGroupDecomposition"
)]
pub struct PyRGroupDecomposition {
    decomp: RGroupDecomposition,
}

#[pymethods]
impl PyRGroupDecomposition {
    /// Construct from a molecule or sequence of molecules and a parameters object
    #[new]
    #[pyo3(signature = (cores, params = None))]
    fn py_new(
        py: Python<'_>,
        cores: &Bound<'_, PyAny>,
        params: Option<&Bound<'_, PyRGroupDecompositionParameters>>,
    ) -> PyResult<Self> {
        let cores = cores_from_py(cores)?;
        let params = params_of(py, params);
        let refs: Vec<&ROMol> = cores.iter().collect();
        Ok(PyRGroupDecomposition {
            decomp: RGroupDecomposition::new_multi(&refs, params).py()?,
        })
    }

    #[pyo3(name = "Add", signature = (mol))]
    fn add(&mut self, py: Python<'_>, mol: &Bound<'_, PyAny>) -> PyResult<i32> {
        let mol = mol_from_py(mol)?;
        let decomp = &mut self.decomp;
        py.detach(|| decomp.add(&mol)).py()
    }

    #[pyo3(name = "GetMatchingCoreIdx", signature = (mol, matches = None))]
    fn matching_core_idx(
        &self,
        py: Python<'_>,
        mol: &Bound<'_, PyAny>,
        matches: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<i32> {
        let mol = mol_from_py(mol)?;
        let decomp = &self.decomp;
        let (idx, found) = py.detach(|| decomp.matching_core_idx(&mol)).py()?;
        if let Some(matches) = matches.filter(|m| !m.is_none()) {
            if let Ok(list) = matches.cast::<PyList>() {
                for m in found {
                    let pairs = m
                        .into_iter()
                        .map(|p| PyTuple::new(py, [p.0, p.1]))
                        .collect::<PyResult<Vec<_>>>()?;
                    list.append(PyTuple::new(py, pairs)?)?;
                }
            }
        }
        Ok(idx)
    }

    /// Process the rgroups (must be done prior to GetRGroupsAsRows/Columns and GetRGroupLabels)
    #[pyo3(name = "Process")]
    fn process(&mut self, py: Python<'_>) -> bool {
        let decomp = &mut self.decomp;
        py.detach(|| decomp.process())
    }

    /// Process the rgroups and returns the score (must be done prior to GetRGroupsAsRows/Columns and GetRGroupLabels)
    #[pyo3(name = "ProcessAndScore")]
    fn process_and_score(&mut self, py: Python<'_>) -> (bool, f64) {
        let decomp = &mut self.decomp;
        let r = py.detach(|| decomp.process_and_score());
        (r.success, r.score)
    }

    /// Return the current list of found rgroups.
    /// Note, Process() should be called first
    #[pyo3(name = "GetRGroupLabels")]
    fn rgroup_labels(&self) -> PyResult<Vec<String>> {
        self.decomp.rgroup_labels().py()
    }

    /// Return the rgroups as rows (note: can be fed directly into a pandas datatable)
    ///   ARGUMENTS:
    ///    - asSmiles: if True return smiles strings, otherwise return molecules [default: False]
    ///     Row structure:
    ///        rows[idx] = {rgroup_label: molecule_or_smiles}
    #[pyo3(name = "GetRGroupsAsRows", signature = (asSmiles = false))]
    #[allow(non_snake_case)]
    fn rgroups_as_rows<'py>(
        &self,
        py: Python<'py>,
        asSmiles: bool,
    ) -> PyResult<Bound<'py, PyList>> {
        rows_to_py(py, self.decomp.rgroups_as_rows().py()?, asSmiles)
    }

    /// Return the rgroups as columns (note: can be fed directly into a pandas datatable)
    ///   ARGUMENTS:
    ///    - asSmiles: if True return smiles strings, otherwise return molecules [default: False]
    ///     Column structure:
    ///        columns[rgroup_label] = [ mols_or_smiles ]
    #[pyo3(name = "GetRGroupsAsColumns", signature = (asSmiles = false))]
    #[allow(non_snake_case)]
    fn rgroups_as_columns<'py>(
        &self,
        py: Python<'py>,
        asSmiles: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        columns_to_py(py, self.decomp.rgroups_as_columns().py()?, asSmiles)
    }
}

/// Decompose a collection of molecules into their Rgroups
///   ARGUMENTS:
///     - cores: a set of cores from most to least specific.
///              See RGroupDecompositionParameters for more details
///              on how the cores can be labelled
///     - mols: the molecules to be decomposed
///     - asSmiles: if True return smiles strings, otherwise return molecules [default: False]
///     - asRows: return the results as rows (default) otherwise return columns
///     - options: RGroupDecompositionParameters object that defines the parameters for the decomposition.
///              See RGroupDecompositionParameters for defaults
///
///   RETURNS: row_or_column_results, unmatched
///
///     Row structure:
///        rows[idx] = {rgroup_label: molecule_or_smiles}
///     Column structure:
///        columns[rgroup_label] = [ mols_or_smiles ]
///
///     unmatched is a vector of indices in the input mols that were not matched.
#[pyfunction]
#[pyo3(name = "RGroupDecompose", signature = (cores, mols, asSmiles = false, asRows = true, options = None))]
#[allow(non_snake_case)]
fn rgroup_decompose<'py>(
    py: Python<'py>,
    cores: &Bound<'py, PyAny>,
    mols: &Bound<'py, PyAny>,
    asSmiles: bool,
    asRows: bool,
    options: Option<&Bound<'py, PyRGroupDecompositionParameters>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let params = params_of(py, options);
    let timeout = params.timeout;
    let cores = cores_from_py(cores)?;
    let mols = mols_from_py(mols)?;
    let core_refs: Vec<&ROMol> = cores.iter().collect();
    let (decomp, unmatched) = py
        .detach(|| -> rdoxide::Result<_> {
            let t0 = std::time::Instant::now();
            let mut decomp = RGroupDecomposition::new_multi(&core_refs, params)?;
            let mut unmatched = Vec::new();
            for (i, mol) in mols.iter().enumerate() {
                if decomp.add(mol)? == -1 {
                    unmatched.push(i as u32);
                }
                check_for_timeout(t0, timeout)?;
            }
            decomp.process();
            Ok((decomp, unmatched))
        })
        .py()?;
    let result: Bound<'py, PyAny> = if asRows {
        rows_to_py(py, decomp.rgroups_as_rows().py()?, asSmiles)?.into_any()
    } else {
        columns_to_py(py, decomp.rgroups_as_columns().py()?, asSmiles)?.into_any()
    };
    PyTuple::new(py, [result, PyList::new(py, unmatched)?.into_any()])
}

/// Relabel dummy atoms bearing an R-group mapping (as
/// atom map number, isotope or MDLRGroup label) such that
/// they will be displayed by the rendering code as R# rather
/// than #*, *:#, #*:#, etc. By default, only the MDLRGroup label
/// is retained on output; this may be configured through the
/// outputLabels parameter.
/// In case there are multiple potential R-group mappings,
/// the priority on input is Atom map number > Isotope > MDLRGroup.
/// The inputLabels parameter allows to configure which mappings
/// are taken into consideration.
#[pyfunction]
#[pyo3(
    name = "RelabelMappedDummies",
    signature = (mol, inputLabels = rgroup_labelling::ATOM_MAP | rgroup_labelling::ISOTOPE | rgroup_labelling::MDL_RGROUP, outputLabels = rgroup_labelling::MDL_RGROUP)
)]
#[allow(non_snake_case)]
fn relabel_mapped_dummies_py(
    mol: &Bound<'_, PyMol>,
    inputLabels: u32,
    outputLabels: u32,
) -> PyResult<()> {
    relabel_mapped_dummies(&mut mol.borrow_mut().mol, inputLabels, outputLabels);
    Ok(())
}

fn members(names: &[(&'static str, u32)]) -> Vec<(&'static str, i64)> {
    names.iter().map(|&(n, v)| (n, i64::from(v))).collect()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.setattr(
        "__doc__",
        "Module containing RGroupDecomposition classes and functions.",
    )?;
    debug_assert_eq!(m.getattr("__name__")?.extract::<String>()?, MODULE);
    add_enum_exported(m, "RGroupLabels", &members(&rgroup_labels::NAMES), None)?;
    add_enum_exported(m, "RGroupMatching", &members(&rgroup_matching::NAMES), None)?;
    add_enum_exported(
        m,
        "RGroupLabelling",
        &members(&rgroup_labelling::NAMES),
        None,
    )?;
    let mut alignment = members(&rgroup_core_alignment::NAMES);
    alignment.insert(0, ("None", i64::from(rgroup_core_alignment::NO_ALIGNMENT)));
    add_enum_exported(m, "RGroupCoreAlignment", &alignment, None)?;
    add_enum_exported(m, "RGroupScore", &members(&rgroup_score::NAMES), None)?;
    m.add_class::<PyRGroupDecompositionParameters>()?;
    m.add_class::<PyRGroupDecomposition>()?;
    m.add_function(wrap_pyfunction!(rgroup_decompose, m)?)?;
    m.add_function(wrap_pyfunction!(relabel_mapped_dummies_py, m)?)?;
    Ok(())
}
