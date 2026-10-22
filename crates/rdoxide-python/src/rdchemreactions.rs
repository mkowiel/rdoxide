//! `rdoxide.Chem.rdChemReactions` (upstream `ChemReactions/Wrap/rdChemReactions.cpp`):
//! chemical reactions defined by templates, and running them on molecules.
//!
//! Template accessors (`GetReactantTemplate`, `GetReactants`, ...) return
//! copies where upstream returns references into the reaction. The MDL/Marvin/
//! CDXML/PNG readers and writers, fingerprints, reaction depiction and
//! enumeration are not provided yet.

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyBytes, PyDict, PyList, PyTuple};
use rdoxide::graphmol::chem_reactions::daylight::{
    ReactionSmartsParserParams, reaction_from_smarts, reaction_from_smiles,
};
use rdoxide::graphmol::chem_reactions::fingerprints::{
    FingerprintType, ReactionFingerprintParams, difference_fingerprint, structural_fingerprint,
};
use rdoxide::graphmol::chem_reactions::mdl::{
    reaction_from_rxn_block, reaction_from_rxn_file, reaction_from_rxn_mol, reaction_to_rxn_block,
    reaction_to_rxn_mol, reaction_to_v3k_rxn_block,
};
use rdoxide::graphmol::chem_reactions::pickler::{
    pickle_reaction, pickle_reaction_default, reaction_from_pickle,
};
use rdoxide::graphmol::chem_reactions::preprocess::{
    DEFAULT_PROP_NAME, ReactantLabels, add_recursive_queries_to_reaction, preprocess_reaction,
};
use rdoxide::graphmol::chem_reactions::sanitize::{
    chem_draw_rxn_adjust_params, default_rxn_adjust_params, flags as sf,
    match_only_at_rgroups_adjust_params, sanitize_rxn, sanitize_rxn_as_mols,
};
use rdoxide::graphmol::chem_reactions::writer::{
    reaction_to_cx_smarts, reaction_to_cx_smiles, reaction_to_smarts, reaction_to_smiles,
};
use rdoxide::graphmol::chem_reactions::{self as cr, ChemicalReaction, runner};
use rdoxide::graphmol::file_parsers::MolFileParserParams;
use rdoxide::graphmol::file_parsers::marvin::{
    MrvParserParams, mrv_block_is_reaction, mrv_file_is_reaction, reaction_from_mrv_block,
    reaction_from_mrv_file, reaction_to_mrv_block, reaction_to_mrv_file,
};
use rdoxide::graphmol::functional_group_hierarchy::flattened_functional_group_hierarchy;
use rdoxide::graphmol::romol::ROMol;
use rdoxide::graphmol::smiles_parse::cx_write::cx_fields;
use rdoxide::graphmol::smiles_parse::write::SmilesWriteParams;

use crate::datastructs::PyExplicitBitVect;
use crate::datastructs_siv::PyUIntSIV;
use crate::enums::{add_enum, add_enum_exported, to_py as enum_to_py};
use crate::errors::{IntoPyResult, to_py};
use crate::mol::PyMol;
use crate::rdmolfiles_write::PySmilesWriteParams;
use crate::rdmolops_adjust_query::PyAdjustQueryParameters;
use crate::substruct::PySubstructMatchParameters;

/// A class for storing and applying chemical reactions.
///
/// Sample Usage:
///   >>> from rdoxide import Chem
///   >>> from rdoxide.Chem import rdChemReactions
///   >>> rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
///   >>> rxn.GetNumReactantTemplates()
///   2
///   >>> rxn.GetNumProductTemplates()
///   1
///   >>> ps = rxn.RunReactants((Chem.MolFromSmiles('CC(=O)O'),Chem.MolFromSmiles('NC')))
///   >>> len(ps)
///   1
///   >>> Chem.MolToSmiles(ps[0][0])
///   'CNC(C)=O'
#[pyclass(
    name = "ChemicalReaction",
    module = "rdoxide.Chem.rdChemReactions",
    dict,
    skip_from_py_object
)]
pub struct PyChemicalReaction {
    pub(crate) rxn: ChemicalReaction,
    /// The object handed out by `GetSubstructParams`; it stays live, so edits
    /// to it must reach `rxn` before the reaction is used.
    params: Option<Py<PySubstructMatchParameters>>,
    /// Template molecules handed out to Python. Upstream returns references
    /// into the reaction, so edits to them (and changes the reaction makes
    /// to the templates) have to show up on both sides: [`Self::pull`] copies
    /// the Python-side state in before the reaction is used, [`Self::push`]
    /// copies the templates back out after it may have changed them.
    shared: Vec<SharedTemplate>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Role {
    Reactant,
    Product,
    Agent,
}

struct SharedTemplate {
    role: Role,
    idx: usize,
    mol: Py<PyMol>,
}

impl PyChemicalReaction {
    pub(crate) fn new(rxn: ChemicalReaction) -> Self {
        PyChemicalReaction {
            rxn,
            params: None,
            shared: Vec::new(),
        }
    }

    fn templates_mut(&mut self, role: Role) -> &mut Vec<ROMol> {
        match role {
            Role::Reactant => self.rxn.reactants_mut(),
            Role::Product => self.rxn.products_mut(),
            Role::Agent => self.rxn.agents_mut(),
        }
    }

    /// The Python object for template `idx` (the same object every time).
    fn shared_template(&mut self, py: Python<'_>, role: Role, idx: usize) -> PyResult<Py<PyMol>> {
        if let Some(s) = self.shared.iter().find(|s| s.role == role && s.idx == idx) {
            return Ok(s.mol.clone_ref(py));
        }
        let mol = Py::new(
            py,
            PyMol {
                mol: self.templates_mut(role)[idx].clone(),
            },
        )?;
        self.shared.push(SharedTemplate {
            role,
            idx,
            mol: mol.clone_ref(py),
        });
        Ok(mol)
    }

    /// All templates of a kind as a tuple of the shared molecules.
    fn shared_templates(&mut self, py: Python<'_>, role: Role) -> PyResult<Py<PyTuple>> {
        self.sync_params(py);
        let n = self.templates_mut(role).len();
        let mols = (0..n)
            .map(|i| self.shared_template(py, role, i))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyTuple::new(py, mols)?.unbind())
    }

    /// Copies edits made through the handed-out templates into the reaction.
    fn pull_templates(&mut self, py: Python<'_>) {
        let shared = std::mem::take(&mut self.shared);
        for s in &shared {
            if let (Ok(m), Some(slot)) = (
                s.mol.bind(py).try_borrow(),
                self.templates_mut(s.role).get_mut(s.idx),
            ) {
                *slot = m.mol.clone();
            }
        }
        self.shared = shared;
    }

    /// Copies the reaction's templates into the handed-out molecules.
    fn push(&mut self, py: Python<'_>) {
        let shared = std::mem::take(&mut self.shared);
        for s in &shared {
            let template = self.templates_mut(s.role).get(s.idx).cloned();
            if let (Some(t), Ok(mut m)) = (template, s.mol.bind(py).try_borrow_mut()) {
                m.mol = t;
            }
        }
        self.shared = shared;
    }

    /// The templates were restructured: handed-out molecules are detached.
    fn detach_templates(&mut self) {
        self.shared.clear();
    }

    pub fn with_dict<R>(
        &self,
        _py: Python<'_>,
        f: impl FnOnce(&rdoxide::rdgeneral::Dict) -> R,
    ) -> PyResult<R> {
        Ok(f(rdoxide::rdgeneral::RDProps::dict(&self.rxn)))
    }

    pub fn with_dict_mut<R>(
        &mut self,
        _py: Python<'_>,
        f: impl FnOnce(&mut rdoxide::rdgeneral::Dict) -> R,
    ) -> PyResult<R> {
        Ok(f(rdoxide::rdgeneral::RDProps::dict_mut(&mut self.rxn)))
    }

    /// Brings the reaction up to date with Python-side edits: the parameters
    /// from `GetSubstructParams` and the handed-out templates.
    fn sync_params(&mut self, py: Python<'_>) {
        if let Some(p) = &self.params {
            *self.rxn.substruct_params_mut() = p.borrow(py).base_params();
        }
        self.pull_templates(py);
    }

    /// Syncs the parameters and initializes the reaction if that has not happened yet.
    fn prepare(&mut self, py: Python<'_>) {
        self.sync_params(py);
        if !self.rxn.is_initialized() {
            self.rxn.init_reactant_matchers(false);
            self.push(py);
        }
    }
}

/// Borrows the reaction after bringing it up to date with Python-side edits.
pub(crate) fn synced<'py>(
    reaction: &Bound<'py, PyChemicalReaction>,
) -> PyRefMut<'py, PyChemicalReaction> {
    let mut guard = reaction.borrow_mut();
    guard.sync_params(reaction.py());
    guard
}

fn template_index(which: u32, n: usize) -> PyResult<usize> {
    if which as usize >= n {
        return Err(PyValueError::new_err("requested template index too high"));
    }
    Ok(which as usize)
}

pub(crate) fn mol_ref<'a>(obj: &'a Bound<'_, PyAny>) -> PyResult<PyRef<'a, PyMol>> {
    Ok(obj.extract::<PyRef<'_, PyMol>>()?)
}

pub(crate) fn products_to_py(py: Python<'_>, sets: Vec<Vec<ROMol>>) -> PyResult<Py<PyTuple>> {
    let tuples = sets
        .into_iter()
        .map(|set| {
            let mols = set
                .into_iter()
                .map(|mol| Py::new(py, PyMol { mol }))
                .collect::<PyResult<Vec<_>>>()?;
            PyTuple::new(py, mols)
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyTuple::new(py, tuples)?.unbind())
}

/// Moves the removed templates into the `targetList` Python list.
fn extend_list(py: Python<'_>, target: &Bound<'_, PyAny>, mols: Vec<ROMol>) -> PyResult<()> {
    let list = target.cast::<PyList>()?;
    for mol in mols {
        list.append(Py::new(py, PyMol { mol })?)?;
    }
    Ok(())
}

#[pymethods]
#[allow(non_snake_case)]
impl PyChemicalReaction {
    #[new]
    #[pyo3(signature = (other = None))]
    fn py_new(other: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        match other {
            None => Ok(PyChemicalReaction::new(ChemicalReaction::new())),
            Some(o) => {
                if let Ok(r) = o.cast::<PyChemicalReaction>() {
                    return Ok(PyChemicalReaction::new(synced(r).rxn.clone()));
                }
                let bytes: Vec<u8> = match o.cast::<PyBytes>() {
                    Ok(b) => b.as_bytes().to_vec(),
                    Err(_) => o.extract::<String>()?.into_bytes(),
                };
                Ok(PyChemicalReaction::new(reaction_from_pickle(&bytes).py()?))
            }
        }
    }

    /// Returns a binary string representation of the reaction.
    ///
    /// Alternatively: ToBinary(propertyFlags) to choose which properties are included.
    #[pyo3(name = "ToBinary", signature = (propertyFlags = None))]
    #[allow(non_snake_case)]
    fn binary<'py>(
        &mut self,
        py: Python<'py>,
        propertyFlags: Option<u32>,
    ) -> PyResult<Bound<'py, PyBytes>> {
        self.sync_params(py);
        let rxn = &self.rxn;
        let bytes = py
            .detach(|| match propertyFlags {
                Some(flags) => pickle_reaction(rxn, flags),
                None => pickle_reaction_default(rxn),
            })
            .py()?;
        Ok(PyBytes::new(py, &bytes))
    }

    /// `rdkit_pickle_suite`: `(type(self), (self.ToBinary(),), (self.__dict__,))`
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let py = slf.py();
        PyTuple::new(
            py,
            [
                slf.get_type().into_any(),
                PyTuple::new(py, [slf.call_method0("ToBinary")?])?.into_any(),
                PyTuple::new(py, [slf.getattr("__dict__")?])?.into_any(),
            ],
        )
    }

    fn __getstate__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(slf.py(), [slf.getattr("__dict__")?])
    }

    fn __setstate__(slf: &Bound<'_, Self>, state: &Bound<'_, PyAny>) -> PyResult<()> {
        let bad = || {
            PyValueError::new_err(format!(
                "expected 1-item tuple in call to __setstate__; got {state}"
            ))
        };
        let state = state.cast::<PyTuple>().map_err(|_| bad())?;
        if state.len() != 1 {
            return Err(bad());
        }
        slf.getattr("__dict__")?
            .call_method1("update", (state.get_item(0)?,))?;
        Ok(())
    }

    fn __getinitargs__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(slf.py(), [slf.call_method0("ToBinary")?])
    }

    /// returns the number of reactants this reaction expects
    #[pyo3(name = "GetNumReactantTemplates")]
    fn num_reactant_templates(&self) -> usize {
        self.rxn.num_reactant_templates()
    }

    /// returns the number of products this reaction generates
    #[pyo3(name = "GetNumProductTemplates")]
    fn num_product_templates(&self) -> usize {
        self.rxn.num_product_templates()
    }

    /// returns the number of agents this reaction expects
    #[pyo3(name = "GetNumAgentTemplates")]
    fn num_agent_templates(&self) -> usize {
        self.rxn.num_agent_templates()
    }

    /// adds a reactant (a Molecule) to the reaction
    #[pyo3(name = "AddReactantTemplate")]
    fn add_reactant_template(&mut self, py: Python<'_>, mol: &Bound<'_, PyMol>) -> usize {
        self.sync_params(py);
        self.rxn.add_reactant_template(mol.borrow().mol.clone())
    }

    /// adds a product (a Molecule)
    #[pyo3(name = "AddProductTemplate")]
    fn add_product_template(&mut self, py: Python<'_>, mol: &Bound<'_, PyMol>) -> usize {
        self.sync_params(py);
        self.rxn.add_product_template(mol.borrow().mol.clone())
    }

    /// adds a agent (a Molecule)
    #[pyo3(name = "AddAgentTemplate")]
    fn add_agent_template(&mut self, py: Python<'_>, mol: &Bound<'_, PyMol>) -> usize {
        self.sync_params(py);
        self.rxn.add_agent_template(mol.borrow().mol.clone())
    }

    /// Removes molecules with an atom mapping ratio below thresholdUnmappedAtoms from reactant templates to the agent templates or to a given targetList
    #[pyo3(name = "RemoveUnmappedReactantTemplates", signature = (thresholdUnmappedAtoms = 0.2, moveToAgentTemplates = true, targetList = None))]
    fn remove_unmapped_reactant_templates(
        &mut self,
        py: Python<'_>,
        thresholdUnmappedAtoms: f64,
        moveToAgentTemplates: bool,
        targetList: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        self.sync_params(py);
        let mut removed = Vec::new();
        let target = targetList.filter(|t| !t.is_none()).map(|_| &mut removed);
        self.detach_templates();
        self.rxn.remove_unmapped_reactant_templates(
            thresholdUnmappedAtoms,
            moveToAgentTemplates,
            target,
        );
        targetList.map_or(Ok(()), |t| extend_list(py, t, removed))
    }

    /// Removes molecules with an atom mapping ratio below thresholdUnmappedAtoms from product templates to the agent templates or to a given targetList
    #[pyo3(name = "RemoveUnmappedProductTemplates", signature = (thresholdUnmappedAtoms = 0.2, moveToAgentTemplates = true, targetList = None))]
    fn remove_unmapped_product_templates(
        &mut self,
        py: Python<'_>,
        thresholdUnmappedAtoms: f64,
        moveToAgentTemplates: bool,
        targetList: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        self.sync_params(py);
        let mut removed = Vec::new();
        let target = targetList.filter(|t| !t.is_none()).map(|_| &mut removed);
        self.detach_templates();
        self.rxn.remove_unmapped_product_templates(
            thresholdUnmappedAtoms,
            moveToAgentTemplates,
            target,
        );
        targetList.map_or(Ok(()), |t| extend_list(py, t, removed))
    }

    /// Removes agents from reaction. If targetList is provide the agents will be transferred to that list.
    #[pyo3(name = "RemoveAgentTemplates", signature = (targetList = None))]
    fn remove_agent_templates(
        &mut self,
        py: Python<'_>,
        targetList: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        self.sync_params(py);
        let mut removed = Vec::new();
        let target = targetList.filter(|t| !t.is_none()).map(|_| &mut removed);
        self.detach_templates();
        self.rxn.remove_agent_templates(target);
        targetList.map_or(Ok(()), |t| extend_list(py, t, removed))
    }

    /// apply the reaction to a sequence of reactant molecules and return the products as a tuple of tuples.  If maxProducts is not zero, stop the reaction when maxProducts have been generated [default=1000]. The product sets are generated on `numThreads` threads (0: all available); the result does not depend on the thread count.
    #[pyo3(name = "RunReactants", signature = (reactants, maxProducts = 1000, numThreads = 0))]
    fn run_reactants(
        &mut self,
        py: Python<'_>,
        reactants: &Bound<'_, PyAny>,
        maxProducts: u32,
        numThreads: i32,
    ) -> PyResult<Py<PyTuple>> {
        self.prepare(py);
        let mols = reactants
            .try_iter()?
            .map(|r| {
                let r = r?;
                if r.is_none() {
                    return Err(PyValueError::new_err("reaction called with None reactants"));
                }
                Ok(mol_ref(&r)?.mol.clone())
            })
            .collect::<PyResult<Vec<_>>>()?;
        let rxn = &self.rxn;
        let sets = py
            .detach(|| runner::run_reactants(rxn, &mols, maxProducts, numThreads))
            .py()?;
        products_to_py(py, sets)
    }

    /// apply the reaction to a single reactant
    #[pyo3(name = "RunReactant", signature = (reactant, reactionIdx, numThreads = 0))]
    fn run_reactant(
        &mut self,
        py: Python<'_>,
        reactant: &Bound<'_, PyMol>,
        reactionIdx: u32,
        numThreads: i32,
    ) -> PyResult<Py<PyTuple>> {
        self.prepare(py);
        let mol = reactant.borrow().mol.clone();
        let rxn = &self.rxn;
        let sets = py
            .detach(|| runner::run_reactant(rxn, &mol, reactionIdx as usize, numThreads))
            .py()?;
        products_to_py(py, sets)
    }

    /// apply the reaction to a single reactant in place. The reactant itself is modified. This can only be used for single reactant - single product reactions.
    #[pyo3(name = "RunReactantInPlace", signature = (reactant, removeUnmatchedAtoms = true))]
    fn run_reactant_in_place(
        &mut self,
        py: Python<'_>,
        reactant: &Bound<'_, PyMol>,
        removeUnmatchedAtoms: bool,
    ) -> PyResult<bool> {
        self.prepare(py);
        let mut mol = reactant.borrow_mut();
        let rxn = &self.rxn;
        let target = &mut mol.mol;
        py.detach(|| runner::run_reactant_in_place(rxn, target, removeUnmatchedAtoms))
            .py()
    }

    /// initializes the reaction so that it can be used
    #[pyo3(name = "Initialize", signature = (silent = false))]
    fn initialize(&mut self, py: Python<'_>, silent: bool) {
        self.sync_params(py);
        self.rxn.init_reactant_matchers(silent);
        self.push(py);
    }

    /// checks if the reaction is ready for use
    #[pyo3(name = "IsInitialized")]
    fn is_initialized(&self) -> bool {
        self.rxn.is_initialized()
    }

    /// checks the reaction for potential problems, returns (numWarnings,numErrors)
    #[pyo3(name = "Validate", signature = (silent = false))]
    fn validate(&mut self, py: Python<'_>, silent: bool) -> (u32, u32) {
        self.sync_params(py);
        let (_, warnings, errors) = self.rxn.validate(silent);
        self.push(py);
        (warnings, errors)
    }

    /// returns one of our product templates
    #[pyo3(name = "GetProductTemplate")]
    fn product_template(&mut self, py: Python<'_>, which: u32) -> PyResult<Py<PyMol>> {
        self.sync_params(py);
        let i = template_index(which, self.rxn.num_product_templates())?;
        self.shared_template(py, Role::Product, i)
    }

    /// returns one of our reactant templates
    #[pyo3(name = "GetReactantTemplate")]
    fn reactant_template(&mut self, py: Python<'_>, which: u32) -> PyResult<Py<PyMol>> {
        self.sync_params(py);
        let i = template_index(which, self.rxn.num_reactant_templates())?;
        self.shared_template(py, Role::Reactant, i)
    }

    /// returns one of our agent templates
    #[pyo3(name = "GetAgentTemplate")]
    fn agent_template(&mut self, py: Python<'_>, which: u32) -> PyResult<Py<PyMol>> {
        self.sync_params(py);
        let i = template_index(which, self.rxn.num_agent_templates())?;
        self.shared_template(py, Role::Agent, i)
    }

    /// EXPERT USER: indicates that the reaction can have implicit properties
    #[pyo3(name = "_setImplicitPropertiesFlag")]
    fn set_implicit_properties_flag(&mut self, val: bool) {
        self.rxn.set_implicit_properties_flag(val);
    }

    /// EXPERT USER: returns whether or not the reaction can have implicit properties
    #[pyo3(name = "_getImplicitPropertiesFlag")]
    fn implicit_properties_flag(&self) -> bool {
        self.rxn.implicit_properties_flag()
    }

    /// returns whether or not the molecule has a substructure match to one of the reactants.
    #[pyo3(name = "IsMoleculeReactant")]
    fn is_molecule_reactant(&mut self, py: Python<'_>, mol: &Bound<'_, PyMol>) -> PyResult<bool> {
        self.sync_params(py);
        let m = mol.borrow();
        Ok(!self
            .rxn
            .reactant_templates_matching(&m.mol, true)
            .py()?
            .is_empty())
    }

    /// returns whether or not the molecule has a substructure match to one of the products.
    #[pyo3(name = "IsMoleculeProduct")]
    fn is_molecule_product(&mut self, py: Python<'_>, mol: &Bound<'_, PyMol>) -> PyResult<bool> {
        self.sync_params(py);
        let m = mol.borrow();
        Ok(!self
            .rxn
            .product_templates_matching(&m.mol, true)
            .py()?
            .is_empty())
    }

    /// returns whether or not the molecule has a substructure match to one of the agents.
    #[pyo3(name = "IsMoleculeAgent")]
    fn is_molecule_agent(&mut self, py: Python<'_>, mol: &Bound<'_, PyMol>) -> PyResult<bool> {
        self.sync_params(py);
        let m = mol.borrow();
        Ok(self.rxn.agent_template_matching(&m.mol).py()?.is_some())
    }

    /// returns a sequence of sequences with the atoms that change in the reaction
    #[pyo3(name = "GetReactingAtoms", signature = (mappedAtomsOnly = false))]
    fn reacting_atoms(&mut self, py: Python<'_>, mappedAtomsOnly: bool) -> PyResult<Py<PyTuple>> {
        self.sync_params(py);
        let atoms = cr::get_reacting_atoms(&self.rxn, mappedAtomsOnly).py()?;
        let rows = atoms
            .into_iter()
            .map(|a| PyTuple::new(py, a))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyTuple::new(py, rows)?.unbind())
    }

    /// get the reactant templates
    #[pyo3(name = "GetReactants")]
    fn reactants(&mut self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        self.shared_templates(py, Role::Reactant)
    }

    /// get the product templates
    #[pyo3(name = "GetProducts")]
    fn products(&mut self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        self.shared_templates(py, Role::Product)
    }

    /// get the agent templates
    #[pyo3(name = "GetAgents")]
    fn agents(&mut self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        self.shared_templates(py, Role::Agent)
    }

    /// get the parameter object controlling the substructure matching
    #[pyo3(name = "GetSubstructParams")]
    fn substruct_params(&mut self, py: Python<'_>) -> PyResult<Py<PySubstructMatchParameters>> {
        if self.params.is_none() {
            self.params = Some(Py::new(
                py,
                PySubstructMatchParameters::from_core(self.rxn.substruct_params()),
            )?);
        }
        Ok(self
            .params
            .as_ref()
            .map(|p| p.clone_ref(py))
            .expect("set above"))
    }

    /// adds recursive queries and returns reactant labels
    #[pyo3(name = "AddRecursiveQueriesToReaction", signature = (queries = None, propName = DEFAULT_PROP_NAME, getLabels = false))]
    #[allow(non_snake_case)]
    fn add_recursive_queries(
        &mut self,
        py: Python<'_>,
        queries: Option<&Bound<'_, PyDict>>,
        propName: &str,
        getLabels: bool,
    ) -> PyResult<Option<Py<PyTuple>>> {
        let queries = query_map(queries)?;
        self.sync_params(py);
        let labels = add_recursive_queries_to_reaction(&mut self.rxn, &queries, propName).py()?;
        self.push(py);
        getLabels.then(|| labels_to_py(py, labels)).transpose()
    }

    /// Removes all computed properties from the reaction.
    #[pyo3(name = "ClearComputedProps")]
    fn clear_computed_props(&mut self) {
        rdoxide::rdgeneral::RDProps::clear_computed_props(&mut self.rxn);
    }
}

crate::props::impl_props!(
    PyChemicalReaction,
    computed = true,
    dict_defaults = (false, false),
    has_prop = i32
);

fn replacement_map(
    dict: Option<&Bound<'_, PyDict>>,
) -> PyResult<std::collections::BTreeMap<String, String>> {
    dict.map_or_else(
        || Ok(Default::default()),
        |d| {
            d.iter()
                .map(|(k, v)| Ok((k.extract()?, v.extract()?)))
                .collect()
        },
    )
}

fn parse_reaction(
    py: Python<'_>,
    text: &str,
    replacements: Option<&Bound<'_, PyDict>>,
    use_smiles: bool,
) -> PyResult<PyChemicalReaction> {
    let params = ReactionSmartsParserParams {
        replacements: replacement_map(replacements)?,
        ..ReactionSmartsParserParams::new()
    };
    let rxn = py
        .detach(|| {
            if use_smiles {
                reaction_from_smiles(text, &params)
            } else {
                reaction_from_smarts(text, &params)
            }
        })
        .py()?;
    Ok(PyChemicalReaction::new(rxn))
}

/// construct a ChemicalReaction from a reaction SMARTS string.
/// see the documentation for rdoxide.Chem.MolFromSmiles for an explanation
/// of the replacements argument.
#[pyfunction]
#[pyo3(name = "ReactionFromSmarts", signature = (SMARTS, replacements = None, useSmiles = false))]
#[allow(non_snake_case)]
fn py_reaction_from_smarts(
    py: Python<'_>,
    SMARTS: &str,
    replacements: Option<&Bound<'_, PyDict>>,
    useSmiles: bool,
) -> PyResult<PyChemicalReaction> {
    parse_reaction(py, SMARTS, replacements, useSmiles)
}

/// construct a ChemicalReaction from a reaction SMILES string.
/// see the documentation for rdoxide.Chem.MolFromSmiles for an explanation
/// of the replacements argument.
#[pyfunction]
#[pyo3(name = "ReactionFromSmiles", signature = (SMILES, replacements = None))]
#[allow(non_snake_case)]
fn py_reaction_from_smiles(
    py: Python<'_>,
    SMILES: &str,
    replacements: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyChemicalReaction> {
    parse_reaction(py, SMILES, replacements, true)
}

/// `ReactionToSmarts(reaction[, params])` / `ReactionToSmiles(reaction[, canonical | params])`
fn write_params(
    arg: Option<&Bound<'_, PyAny>>,
    allow_bool: bool,
    default_canonical: bool,
) -> PyResult<SmilesWriteParams> {
    match arg {
        None => Ok(SmilesWriteParams {
            canonical: default_canonical,
            ..SmilesWriteParams::default()
        }),
        Some(a) => {
            if let Ok(p) = a.extract::<PyRef<'_, PySmilesWriteParams>>() {
                Ok(p.to_core())
            } else if allow_bool && a.is_instance_of::<PyBool>() {
                Ok(SmilesWriteParams {
                    canonical: a.extract()?,
                    ..SmilesWriteParams::default()
                })
            } else {
                Err(PyValueError::new_err(
                    "params must be a SmilesWriteParams object",
                ))
            }
        }
    }
}

/// construct a reaction SMARTS string for a ChemicalReaction
#[pyfunction]
#[pyo3(name = "ReactionToSmarts", signature = (reaction, params = None))]
fn py_reaction_to_smarts(
    reaction: &Bound<'_, PyChemicalReaction>,
    params: Option<&Bound<'_, PyAny>>,
) -> PyResult<String> {
    let params = write_params(params, false, false)?;
    reaction_to_smarts(&synced(reaction).rxn, &params).py()
}

/// construct a reaction SMILES string for a ChemicalReaction
#[pyfunction]
#[pyo3(name = "ReactionToSmiles", signature = (reaction, params = None))]
fn py_reaction_to_smiles(
    reaction: &Bound<'_, PyChemicalReaction>,
    params: Option<&Bound<'_, PyAny>>,
) -> PyResult<String> {
    let params = write_params(params, true, true)?;
    reaction_to_smiles(&synced(reaction).rxn, &params).py()
}

/// construct a reaction CXSMARTS string for a ChemicalReaction
#[pyfunction]
#[pyo3(name = "ReactionToCXSmarts", signature = (reaction, params = None, flags = cx_fields::ALL))]
fn py_reaction_to_cx_smarts(
    py: Python<'_>,
    reaction: &Bound<'_, PyChemicalReaction>,
    params: Option<&Bound<'_, PyAny>>,
    flags: u32,
) -> PyResult<String> {
    let params = write_params(params, false, false)?;
    let rxn = synced(reaction).rxn.clone();
    py.detach(|| reaction_to_cx_smarts(&rxn, &params, flags))
        .py()
}

/// construct a reaction CXSMILES string for a ChemicalReaction
#[pyfunction]
#[pyo3(name = "ReactionToCXSmiles", signature = (reaction, params = None, flags = cx_fields::ALL))]
fn py_reaction_to_cx_smiles(
    py: Python<'_>,
    reaction: &Bound<'_, PyChemicalReaction>,
    params: Option<&Bound<'_, PyAny>>,
    flags: u32,
) -> PyResult<String> {
    let params = write_params(params, true, true)?;
    let rxn = synced(reaction).rxn.clone();
    py.detach(|| reaction_to_cx_smiles(&rxn, &params, flags))
        .py()
}

const FINGERPRINT_TYPE: &str = "FingerprintType";

fn fingerprint_type_members() -> [(&'static str, i64); 5] {
    [
        ("AtomPairFP", FingerprintType::AtomPair as i64),
        (
            "TopologicalTorsion",
            FingerprintType::TopologicalTorsion as i64,
        ),
        ("MorganFP", FingerprintType::Morgan as i64),
        ("RDKitFP", FingerprintType::RDKit as i64),
        ("PatternFP", FingerprintType::Pattern as i64),
    ]
}

/// A class for storing parameters to manipulate the calculation of
/// fingerprints of chemical reactions.
#[pyclass(
    name = "ReactionFingerprintParams",
    module = "rdoxide.Chem.rdChemReactions",
    from_py_object
)]
#[derive(Clone, Default)]
struct PyReactionFingerprintParams {
    params: ReactionFingerprintParams,
}

#[pymethods]
impl PyReactionFingerprintParams {
    /// Constructor: no arguments, or (includeAgents, bitRatioAgents,
    /// nonAgentWeight, agentWeight, fpSize, fpType)
    #[new]
    #[pyo3(signature = (*args))]
    fn new(args: &Bound<'_, PyTuple>) -> PyResult<Self> {
        match args.len() {
            0 => Ok(Self::default()),
            6 => {
                let fp_type: i64 = args.get_item(5)?.extract()?;
                Ok(PyReactionFingerprintParams {
                    params: ReactionFingerprintParams {
                        include_agents: args.get_item(0)?.extract()?,
                        bit_ratio_agents: args.get_item(1)?.extract()?,
                        non_agent_weight: args.get_item(2)?.extract()?,
                        agent_weight: args.get_item(3)?.extract()?,
                        fp_size: args.get_item(4)?.extract()?,
                        fp_type: fingerprint_type(fp_type)?,
                    },
                })
            }
            _ => Err(PyTypeError::new_err(
                "ReactionFingerprintParams() takes no arguments or (includeAgents, \
                 bitRatioAgents, nonAgentWeight, agentWeight, fpSize, fpType)",
            )),
        }
    }

    #[getter(fpSize)]
    fn fp_size(&self) -> u32 {
        self.params.fp_size
    }
    #[setter(fpSize)]
    fn set_fp_size(&mut self, v: u32) {
        self.params.fp_size = v;
    }
    #[getter(fpType)]
    fn fp_type(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        enum_to_py(py, FINGERPRINT_TYPE, self.params.fp_type as i64)
    }
    #[setter(fpType)]
    fn set_fp_type(&mut self, v: i64) -> PyResult<()> {
        self.params.fp_type = fingerprint_type(v)?;
        Ok(())
    }
    #[getter(bitRatioAgents)]
    fn bit_ratio_agents(&self) -> f64 {
        self.params.bit_ratio_agents
    }
    #[setter(bitRatioAgents)]
    fn set_bit_ratio_agents(&mut self, v: f64) {
        self.params.bit_ratio_agents = v;
    }
    #[getter(nonAgentWeight)]
    fn non_agent_weight(&self) -> u32 {
        self.params.non_agent_weight
    }
    #[setter(nonAgentWeight)]
    fn set_non_agent_weight(&mut self, v: u32) {
        self.params.non_agent_weight = v;
    }
    #[getter(agentWeight)]
    fn agent_weight(&self) -> i32 {
        self.params.agent_weight
    }
    #[setter(agentWeight)]
    fn set_agent_weight(&mut self, v: i32) {
        self.params.agent_weight = v;
    }
    #[getter(includeAgents)]
    fn include_agents(&self) -> bool {
        self.params.include_agents
    }
    #[setter(includeAgents)]
    fn set_include_agents(&mut self, v: bool) {
        self.params.include_agents = v;
    }
}

fn fingerprint_type(v: i64) -> PyResult<FingerprintType> {
    i32::try_from(v)
        .ok()
        .and_then(FingerprintType::from_i32)
        .ok_or_else(|| PyValueError::new_err("invalid FingerprintType"))
}

/// construct a difference fingerprint for a ChemicalReaction by subtracting
/// the reactant fingerprint from the product fingerprint
#[pyfunction]
#[pyo3(name = "CreateDifferenceFingerprintForReaction", signature = (reaction, ReactionFingerPrintParams = None))]
#[allow(non_snake_case)]
fn py_create_difference_fingerprint(
    py: Python<'_>,
    reaction: &Bound<'_, PyChemicalReaction>,
    ReactionFingerPrintParams: Option<PyReactionFingerprintParams>,
) -> PyResult<PyUIntSIV> {
    let params = ReactionFingerPrintParams
        .map_or_else(ReactionFingerprintParams::default_difference, |p| p.params);
    let rxn = synced(reaction).rxn.clone();
    py.detach(|| difference_fingerprint(&rxn, &params))
        .py()
        .map(PyUIntSIV::wrap)
}

/// construct a structural fingerprint for a ChemicalReaction by
/// concatenating the reactant fingerprint and the product fingerprint
#[pyfunction]
#[pyo3(name = "CreateStructuralFingerprintForReaction", signature = (reaction, ReactionFingerPrintParams = None))]
#[allow(non_snake_case)]
fn py_create_structural_fingerprint(
    py: Python<'_>,
    reaction: &Bound<'_, PyChemicalReaction>,
    ReactionFingerPrintParams: Option<PyReactionFingerprintParams>,
) -> PyResult<PyExplicitBitVect> {
    let params = ReactionFingerPrintParams
        .map_or_else(ReactionFingerprintParams::default_structural, |p| p.params);
    let rxn = synced(reaction).rxn.clone();
    py.detach(|| structural_fingerprint(&rxn, &params))
        .py()
        .map(PyExplicitBitVect::wrap)
}

/// construct a molecule for a ChemicalReaction with RXN role property set
#[pyfunction]
#[pyo3(name = "ReactionToMolecule")]
fn py_reaction_to_molecule(reaction: &Bound<'_, PyChemicalReaction>) -> PyResult<PyMol> {
    Ok(PyMol {
        mol: reaction_to_rxn_mol(&synced(reaction).rxn).py()?,
    })
}

fn rxn_parser_params(sanitize: bool, remove_hs: bool, strict_parsing: bool) -> MolFileParserParams {
    MolFileParserParams {
        sanitize,
        remove_hs,
        strict_parsing,
        ..MolFileParserParams::default()
    }
}

/// construct a ChemicalReaction from an MDL rxn file
#[pyfunction]
#[pyo3(name = "ReactionFromRxnFile", signature = (filename, sanitize = false, removeHs = false, strictParsing = true))]
#[allow(non_snake_case)]
fn py_reaction_from_rxn_file(
    py: Python<'_>,
    filename: &str,
    sanitize: bool,
    removeHs: bool,
    strictParsing: bool,
) -> PyResult<Option<PyChemicalReaction>> {
    let params = rxn_parser_params(sanitize, removeHs, strictParsing);
    let rxn = py
        .detach(|| reaction_from_rxn_file(filename, &params))
        .py()?;
    Ok(rxn.map(PyChemicalReaction::new))
}

/// construct a ChemicalReaction from a string in MDL rxn format
#[pyfunction]
#[pyo3(name = "ReactionFromRxnBlock", signature = (rxnblock, sanitize = false, removeHs = false, strictParsing = true))]
#[allow(non_snake_case)]
fn py_reaction_from_rxn_block(
    py: Python<'_>,
    rxnblock: &str,
    sanitize: bool,
    removeHs: bool,
    strictParsing: bool,
) -> PyResult<PyChemicalReaction> {
    let params = rxn_parser_params(sanitize, removeHs, strictParsing);
    let rxn = py
        .detach(|| reaction_from_rxn_block(rxnblock, &params))
        .py()?;
    Ok(PyChemicalReaction::new(rxn))
}

/// Construct reactions from the scheme steps of a CDXML file.
#[pyfunction]
#[pyo3(name = "ReactionsFromCDXMLFile", signature = (filename, sanitize = false, removeHs = false))]
#[allow(non_snake_case)]
fn py_reactions_from_cdxml_file<'py>(
    py: Python<'py>,
    filename: &str,
    sanitize: bool,
    removeHs: bool,
) -> PyResult<Bound<'py, PyTuple>> {
    reactions_tuple(
        py,
        rdoxide::graphmol::chem_reactions::cdxml::cdxml_file_to_chemical_reactions(
            filename, sanitize, removeHs,
        ),
    )
}

/// Construct reactions from the scheme steps of a CDXML block.
#[pyfunction]
#[pyo3(name = "ReactionsFromCDXMLBlock", signature = (cdxml, sanitize = false, removeHs = false))]
#[allow(non_snake_case)]
fn py_reactions_from_cdxml_block<'py>(
    py: Python<'py>,
    cdxml: &Bound<'py, PyAny>,
    sanitize: bool,
    removeHs: bool,
) -> PyResult<Bound<'py, PyTuple>> {
    let text = crate::rdmolfiles_sequence::text_arg(cdxml)?;
    reactions_tuple(
        py,
        rdoxide::graphmol::chem_reactions::cdxml::cdxml_to_chemical_reactions(
            &text, sanitize, removeHs,
        ),
    )
}

/// Upstream logs parse failures as warnings and returns what it has.
fn reactions_tuple<'py>(
    py: Python<'py>,
    res: rdoxide::Result<Vec<ChemicalReaction>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let rxns = match res {
        Ok(r) => r,
        Err(rdoxide::RDKitError::FileParse(m)) => {
            rdoxide::rd_log_warning!("{m}");
            Vec::new()
        }
        Err(_) => Vec::new(),
    };
    PyTuple::new(
        py,
        rxns.into_iter()
            .map(|r| Py::new(py, PyChemicalReaction::new(r)))
            .collect::<PyResult<Vec<_>>>()?,
    )
}

/// Compute 2D coordinates for a reaction. The templates are laid out side by side.
///
///   ARGUMENTS:
///     - reaction - the reaction of interest
///     - spacing - the amount of space left between components of the reaction
///     - canonOrient - orient the reactants and products in a canonical way
///     - updateProps - if set, properties such as conjugation and
///       hybridization will be calculated for the reactant and product
///       templates before generating coordinates. This should result in
///       better depictions, but can lead to errors in some cases.
///     - nFlipsPerSample - number of rotatable bonds that are flipped at
///       random at a time
///     - nSample - number of samples
///     - sampleSeed - seed for the random sampling process
///     - permuteDeg4Nodes - try permuting the drawing order of bonds around
///       atoms with four neighbors in order to improve the depiction
///     - bondLength - change the default bond length for depiction
#[pyfunction]
#[pyo3(name = "Compute2DCoordsForReaction", signature = (reaction, spacing = 1.0, updateProps = true, canonOrient = true, nFlipsPerSample = 0, nSample = 0, sampleSeed = 0, permuteDeg4Nodes = false, bondLength = -1.0))]
#[allow(clippy::too_many_arguments, non_snake_case)]
fn py_compute_2d_coords_for_reaction(
    py: Python<'_>,
    reaction: &Bound<'_, PyChemicalReaction>,
    spacing: f64,
    updateProps: bool,
    canonOrient: bool,
    nFlipsPerSample: u32,
    nSample: u32,
    sampleSeed: i32,
    permuteDeg4Nodes: bool,
    bondLength: f64,
) -> PyResult<()> {
    use rdoxide::graphmol::chem_reactions::depict::{
        ReactionDepictParams, compute_2d_coords_for_reaction,
    };
    use rdoxide::graphmol::depictor::depict_utils::{bond_len, set_bond_len};
    let mut guard = synced(reaction);
    let params = ReactionDepictParams {
        spacing,
        update_props: updateProps,
        canon_orient: canonOrient,
        n_flips_per_sample: nFlipsPerSample,
        n_samples: nSample,
        sample_seed: sampleSeed,
        permute_deg4_nodes: permuteDeg4Nodes,
    };
    let old = bond_len();
    if bondLength > 0.0 {
        set_bond_len(bondLength);
    }
    let res = compute_2d_coords_for_reaction(&mut guard.rxn, &params);
    if bondLength > 0.0 {
        set_bond_len(old);
    }
    guard.push(py);
    res.py()
}

/// construct a string in MDL rxn format for a ChemicalReaction
#[pyfunction]
#[pyo3(name = "ReactionToRxnBlock", signature = (reaction, separateAgents = false, forceV3000 = false))]
#[allow(non_snake_case)]
fn py_reaction_to_rxn_block(
    reaction: &Bound<'_, PyChemicalReaction>,
    separateAgents: bool,
    forceV3000: bool,
) -> PyResult<String> {
    reaction_to_rxn_block(&synced(reaction).rxn, separateAgents, forceV3000).py()
}

/// construct a string in MDL v3000 rxn format for a ChemicalReaction
#[pyfunction]
#[pyo3(name = "ReactionToV3KRxnBlock", signature = (reaction, separateAgents = false))]
#[allow(non_snake_case)]
fn py_reaction_to_v3k_rxn_block(
    reaction: &Bound<'_, PyChemicalReaction>,
    separateAgents: bool,
) -> PyResult<String> {
    reaction_to_v3k_rxn_block(&synced(reaction).rxn, separateAgents).py()
}

/// Maps the errors of a Marvin reaction reader the way upstream does: a
/// missing file is an `IOError`, everything else gives `None`.
fn mrv_reaction_result(
    res: rdoxide::Result<ChemicalReaction>,
) -> PyResult<Option<PyChemicalReaction>> {
    match res {
        Ok(rxn) => Ok(Some(PyChemicalReaction::new(rxn))),
        Err(e @ rdoxide::RDKitError::BadFile(_)) => Err(to_py(e)),
        Err(rdoxide::RDKitError::FileParse(m)) => {
            rdoxide::rd_log_warning!("{m}");
            Ok(None)
        }
        Err(_) => Ok(None),
    }
}

/// construct a ChemicalReaction from a Marvin (mrv) rxn file
#[pyfunction]
#[pyo3(name = "ReactionFromMrvFile", signature = (filename, sanitize = false, removeHs = false))]
#[allow(non_snake_case)]
fn py_reaction_from_mrv_file(
    py: Python<'_>,
    filename: &str,
    sanitize: bool,
    removeHs: bool,
) -> PyResult<Option<PyChemicalReaction>> {
    let params = MrvParserParams {
        sanitize,
        remove_hs: removeHs,
    };
    mrv_reaction_result(py.detach(|| reaction_from_mrv_file(filename, &params)))
}

/// construct a ChemicalReaction from a string in Marvin (mrv) format
#[pyfunction]
#[pyo3(name = "ReactionFromMrvBlock", signature = (rxnblock, sanitize = false, removeHs = false))]
#[allow(non_snake_case)]
fn py_reaction_from_mrv_block(
    py: Python<'_>,
    rxnblock: &Bound<'_, PyAny>,
    sanitize: bool,
    removeHs: bool,
) -> PyResult<Option<PyChemicalReaction>> {
    let text = match rxnblock.extract::<String>() {
        Ok(s) => s,
        Err(_) => String::from_utf8_lossy(&rxnblock.extract::<Vec<u8>>()?).into_owned(),
    };
    let params = MrvParserParams {
        sanitize,
        remove_hs: removeHs,
    };
    mrv_reaction_result(py.detach(|| reaction_from_mrv_block(&text, &params)))
}

/// returns whether or not an MRV file contains reaction data
#[pyfunction]
#[pyo3(name = "MrvFileIsReaction", signature = (filename))]
fn py_mrv_file_is_reaction(filename: &str) -> PyResult<bool> {
    mrv_file_is_reaction(filename).py()
}

/// returns whether or not an MRV block contains reaction data
#[pyfunction]
#[pyo3(name = "MrvBlockIsReaction", signature = (mrvData))]
#[allow(non_snake_case)]
fn py_mrv_block_is_reaction(mrvData: &str) -> PyResult<bool> {
    mrv_block_is_reaction(mrvData).py()
}

/// construct a string in Marvin (MRV) format for a ChemicalReaction
#[pyfunction]
#[pyo3(name = "ReactionToMrvBlock", signature = (reaction, prettyPrint = false))]
#[allow(non_snake_case)]
fn py_reaction_to_mrv_block(
    reaction: &Bound<'_, PyChemicalReaction>,
    prettyPrint: bool,
) -> PyResult<String> {
    reaction_to_mrv_block(&synced(reaction).rxn, prettyPrint).py()
}

/// write a Marvin (MRV) file for a ChemicalReaction
#[pyfunction]
#[pyo3(name = "ReactionToMrvFile", signature = (reaction, filename, prettyPrint = false))]
#[allow(non_snake_case)]
fn py_reaction_to_mrv_file(
    reaction: &Bound<'_, PyChemicalReaction>,
    filename: &str,
    prettyPrint: bool,
) -> PyResult<()> {
    reaction_to_mrv_file(&synced(reaction).rxn, filename, prettyPrint).py()
}

/// construct a ChemicalReaction from an molecule if the RXN role property of the molecule is set
#[pyfunction]
#[pyo3(name = "ReactionFromMolecule")]
fn py_reaction_from_molecule(mol: &Bound<'_, PyMol>) -> PyResult<PyChemicalReaction> {
    Ok(PyChemicalReaction::new(
        reaction_from_rxn_mol(&mol.borrow().mol).py()?,
    ))
}

/// tests if a molecule can be classified as an agent depending on the ratio of mapped atoms and a give threshold
#[pyfunction]
#[pyo3(name = "IsReactionTemplateMoleculeAgent")]
#[allow(non_snake_case)]
fn py_is_agent(molecule: &Bound<'_, PyMol>, agentThreshold: f64) -> bool {
    cr::is_reaction_template_molecule_agent(&molecule.borrow().mol, agentThreshold)
}

/// tests if a reaction obtains any atom mapping
#[pyfunction]
#[pyo3(name = "HasReactionAtomMapping")]
fn py_has_atom_mapping(rxn: &Bound<'_, PyChemicalReaction>) -> bool {
    cr::has_reaction_atom_mapping(&synced(rxn).rxn)
}

/// tests if the queryReaction is a substructure of a reaction
#[pyfunction]
#[pyo3(name = "HasReactionSubstructMatch", signature = (reaction, queryReaction, includeAgents = false))]
#[allow(non_snake_case)]
fn py_has_reaction_substruct_match(
    reaction: &Bound<'_, PyChemicalReaction>,
    queryReaction: &Bound<'_, PyChemicalReaction>,
    includeAgents: bool,
) -> PyResult<bool> {
    cr::has_reaction_substruct_match(
        &synced(reaction).rxn,
        &synced(queryReaction).rxn,
        includeAgents,
    )
    .py()
}

macro_rules! template_match_fn {
    ($(#[$doc:meta])* $py_name:literal, $fn_name:ident, $core:path) => {
        $(#[$doc])*
        #[pyfunction]
        #[pyo3(name = $py_name)]
        #[allow(non_snake_case)]
        fn $fn_name(
            reaction: &Bound<'_, PyChemicalReaction>,
            queryReaction: &Bound<'_, PyChemicalReaction>,
        ) -> PyResult<bool> {
                    $core(&synced(reaction).rxn, &synced(queryReaction).rxn).py()
        }
    };
}

template_match_fn!(
    /// tests if the agents of a queryReaction are the same as those of a reaction
    "HasAgentTemplateSubstructMatch",
    py_has_agent_match,
    cr::has_agent_template_substruct_match
);
template_match_fn!(
    /// tests if the products of a queryReaction are substructures of the products of a reaction
    "HasProductTemplateSubstructMatch",
    py_has_product_match,
    cr::has_product_template_substruct_match
);
template_match_fn!(
    /// tests if the reactants of a queryReaction are substructures of the reactants of a reaction
    "HasReactantTemplateSubstructMatch",
    py_has_reactant_match,
    cr::has_reactant_template_substruct_match
);

/// Caution: This is an expert-user function which will change a property (molInversionFlag) of your products.
/// This function is called by default using the RXN or SMARTS parser for reactions and should really only be called if reactions have been constructed some other way.
/// The function updates the stereochemistry of the product by considering 4 different cases: inversion, retention, removal, and introduction
#[pyfunction]
#[pyo3(name = "UpdateProductsStereochemistry")]
fn py_update_products_stereochem(reaction: &Bound<'_, PyChemicalReaction>) -> PyResult<()> {
    let mut guard = synced(reaction);
    let result = cr::update_products_stereochem(&mut guard.rxn).py();
    guard.push(reaction.py());
    result
}

/// reduce the product of a reaction to the side chains added by the reaction.
/// The output is a molecule with attached wildcards indicating where the product was attached.
/// The dummy atom has the same reaction-map number as the product atom (if available).
#[pyfunction]
#[pyo3(name = "ReduceProductToSideChains", signature = (product, addDummyAtoms = true))]
#[allow(non_snake_case)]
fn py_reduce_product_to_side_chains(
    product: &Bound<'_, PyMol>,
    addDummyAtoms: bool,
) -> PyResult<PyMol> {
    Ok(PyMol {
        mol: runner::reduce_product_to_side_chains(&product.borrow().mol, addDummyAtoms).py()?,
    })
}

/// Removes the mapping numbers from the molecules of a reaction
#[pyfunction]
#[pyo3(name = "RemoveMappingNumbersFromReactions")]
fn py_remove_mapping_numbers(reaction: &Bound<'_, PyChemicalReaction>) {
    let mut guard = synced(reaction);
    cr::remove_mapping_numbers_from_reactions(&mut guard.rxn);
    guard.push(reaction.py());
}

fn query_map(
    dict: Option<&Bound<'_, PyDict>>,
) -> PyResult<std::collections::BTreeMap<String, ROMol>> {
    dict.map_or_else(
        || Ok(Default::default()),
        |d| {
            d.iter()
                .map(|(k, v)| Ok((k.extract()?, mol_ref(&v)?.mol.clone())))
                .collect()
        },
    )
}

fn labels_to_py(py: Python<'_>, labels: ReactantLabels) -> PyResult<Py<PyTuple>> {
    let per_template = labels
        .into_iter()
        .map(|template| {
            let pairs = template
                .into_iter()
                .map(|(idx, name)| {
                    PyTuple::new(
                        py,
                        [
                            idx.into_pyobject(py)?.into_any(),
                            name.into_pyobject(py)?.into_any(),
                        ],
                    )
                })
                .collect::<PyResult<Vec<_>>>()?;
            PyTuple::new(py, pairs)
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyTuple::new(py, per_template)?.unbind())
}

#[doc = include_str!("docs/preprocess_reaction.txt")]
#[allow(clippy::doc_lazy_continuation, clippy::doc_overindented_list_items)]
#[pyfunction]
#[pyo3(name = "PreprocessReaction", signature = (reaction, queries = None, propName = DEFAULT_PROP_NAME))]
#[allow(non_snake_case)]
fn py_preprocess_reaction(
    py: Python<'_>,
    reaction: &Bound<'_, PyChemicalReaction>,
    queries: Option<&Bound<'_, PyDict>>,
    propName: &str,
) -> PyResult<Py<PyTuple>> {
    let mut queries = query_map(queries)?;
    if queries.is_empty() {
        queries = flattened_functional_group_hierarchy(true).clone();
    }
    let mut guard = synced(reaction);
    let n_reactants = guard.rxn.num_reactant_templates();
    let n_products = guard.rxn.num_product_templates();
    let (_, mut warnings, mut errors) = guard.rxn.validate(false);
    let mut labels = Vec::new();
    if errors == 0 {
        let result = preprocess_reaction(&mut guard.rxn, &queries, propName).py()?;
        (warnings, errors, labels) = (
            result.num_warnings,
            result.num_errors,
            result.reactant_labels,
        );
    }
    guard.push(py);
    PyTuple::new(
        py,
        [
            warnings.into_pyobject(py)?.into_any(),
            errors.into_pyobject(py)?.into_any(),
            n_reactants.into_pyobject(py)?.into_any(),
            n_products.into_pyobject(py)?.into_any(),
            labels_to_py(py, labels)?.into_bound(py).into_any(),
        ],
    )
    .map(Bound::unbind)
}

/// Returns the default adjustment parameters for reactant templates
#[pyfunction]
#[pyo3(name = "GetDefaultAdjustParams")]
fn py_default_adjust_params() -> PyAdjustQueryParameters {
    PyAdjustQueryParameters {
        params: default_rxn_adjust_params(),
    }
}

/// (deprecated, see MatchOnlyAtRgroupsAdjustParams)
///     Returns the chemdraw style adjustment parameters for reactant templates
#[pyfunction]
#[pyo3(name = "GetChemDrawRxnAdjustParams")]
fn py_chem_draw_adjust_params() -> PyAdjustQueryParameters {
    PyAdjustQueryParameters {
        params: chem_draw_rxn_adjust_params(),
    }
}

/// Only match at the specified rgroup locations in the reactant templates
#[pyfunction]
#[pyo3(name = "MatchOnlyAtRgroupsAdjustParams")]
fn py_match_only_at_rgroups_adjust_params() -> PyAdjustQueryParameters {
    PyAdjustQueryParameters {
        params: match_only_at_rgroups_adjust_params(),
    }
}

/// Does some sanitization of the reactant and product templates of a reaction.
///
///     - The reaction is modified in place.
///     - If sanitization fails, an exception will be thrown unless catchErrors is set
///
///   ARGUMENTS:
///
///     - rxn: the reaction to be modified
///     - sanitizeOps: (optional) reaction sanitization operations to be carried out
///       these should be constructed by or'ing together the
///       operations in rdoxide.Chem.rdChemReactions.SanitizeFlags
///     - optional adjustment parameters for changing the meaning of the substructure
///       matching done in the templates.  The default is
///       rdoxide.Chem.rdChemReactions.DefaultRxnAdjustParams which aromatizes
///       kekule structures if possible.
///     - catchErrors: (optional) if provided, instead of raising an exception
///       when sanitization fails (the default behavior), the
///       first operation that failed (as defined in rdoxide.Chem.rdChemReactions.SanitizeFlags)
///       is returned. Zero is returned on success.
///
///   The operations carried out by default are:
///     1) fixRGroups(): sets R group labels on mapped dummy atoms when possible
///     2) fixAtomMaps(): attempts to set atom maps on unmapped R groups
///     3) adjustTemplate(): calls adjustQueryProperties() on all reactant templates
///     4) fixHs(): merges explicit Hs in the reactant templates that don't map to heavy atoms
#[pyfunction]
#[pyo3(name = "SanitizeRxn", signature = (rxn, sanitizeOps = sf::SANITIZE_ALL, params = None, catchErrors = false))]
#[allow(non_snake_case)]
fn py_sanitize_rxn(
    py: Python<'_>,
    rxn: &Bound<'_, PyChemicalReaction>,
    sanitizeOps: u32,
    params: Option<&Bound<'_, PyAdjustQueryParameters>>,
    catchErrors: bool,
) -> PyResult<Py<PyAny>> {
    let params = params.map_or_else(default_rxn_adjust_params, |p| p.borrow().params.clone());
    let mut guard = synced(rxn);
    let result = sanitize_rxn(&mut guard.rxn, sanitizeOps, &params);
    guard.push(py);
    match result {
        Ok(()) => enum_to_py(py, "SanitizeFlags", i64::from(sf::SANITIZE_NONE)),
        Err((op, _)) if catchErrors => enum_to_py(py, "SanitizeFlags", i64::from(op)),
        Err((_, e)) => Err(crate::errors::to_py(e)),
    }
}

/// Does the usual molecular sanitization on each reactant, agent, and product of the reaction
#[pyfunction]
#[pyo3(name = "SanitizeRxnAsMols", signature = (rxn, sanitizeOps = 0xFFF_FFFF))]
#[allow(non_snake_case)]
fn py_sanitize_rxn_as_mols(rxn: &Bound<'_, PyChemicalReaction>, sanitizeOps: u32) -> PyResult<()> {
    let mut guard = synced(rxn);
    let result = sanitize_rxn_as_mols(&mut guard.rxn, sanitizeOps).py();
    guard.push(rxn.py());
    result
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing classes and functions for working with chemical reactions.",
    )?;
    m.add_class::<PyChemicalReaction>()?;
    m.add_function(wrap_pyfunction!(py_compute_2d_coords_for_reaction, m)?)?;
    m.add_function(wrap_pyfunction!(py_reactions_from_cdxml_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_reactions_from_cdxml_block, m)?)?;
    crate::rdchemreactions_enumerate::init(m)?;
    crate::rdchemreactions_png::init(m)?;
    add_enum(m, FINGERPRINT_TYPE, &fingerprint_type_members(), None)?;
    m.add_class::<PyReactionFingerprintParams>()?;
    m.add_function(wrap_pyfunction!(py_create_difference_fingerprint, m)?)?;
    m.add_function(wrap_pyfunction!(py_create_structural_fingerprint, m)?)?;
    add_enum_exported(
        m,
        "SanitizeFlags",
        &[
            ("SANITIZE_NONE", i64::from(sf::SANITIZE_NONE)),
            ("SANITIZE_ATOM_MAPS", i64::from(sf::SANITIZE_ATOM_MAPS)),
            (
                "SANITIZE_RGROUP_NAMES",
                i64::from(sf::SANITIZE_RGROUP_NAMES),
            ),
            (
                "SANITIZE_ADJUST_REACTANTS",
                i64::from(sf::SANITIZE_ADJUST_REACTANTS),
            ),
            ("SANITIZE_MERGEHS", i64::from(sf::SANITIZE_MERGEHS)),
            ("SANITIZE_ALL", i64::from(sf::SANITIZE_ALL)),
        ],
        None,
    )?;
    m.add_function(wrap_pyfunction!(py_preprocess_reaction, m)?)?;
    m.add_function(wrap_pyfunction!(py_default_adjust_params, m)?)?;
    m.add_function(wrap_pyfunction!(py_chem_draw_adjust_params, m)?)?;
    m.add_function(wrap_pyfunction!(py_match_only_at_rgroups_adjust_params, m)?)?;
    m.add_function(wrap_pyfunction!(py_sanitize_rxn, m)?)?;
    m.add_function(wrap_pyfunction!(py_sanitize_rxn_as_mols, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_from_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_from_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_to_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_to_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_to_cx_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_to_cx_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_to_molecule, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_from_molecule, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_from_rxn_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_from_rxn_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_to_rxn_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_to_v3k_rxn_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_from_mrv_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_from_mrv_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_mrv_file_is_reaction, m)?)?;
    m.add_function(wrap_pyfunction!(py_mrv_block_is_reaction, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_to_mrv_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_to_mrv_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_is_agent, m)?)?;
    m.add_function(wrap_pyfunction!(py_has_atom_mapping, m)?)?;
    m.add_function(wrap_pyfunction!(py_has_reaction_substruct_match, m)?)?;
    m.add_function(wrap_pyfunction!(py_has_agent_match, m)?)?;
    m.add_function(wrap_pyfunction!(py_has_product_match, m)?)?;
    m.add_function(wrap_pyfunction!(py_has_reactant_match, m)?)?;
    m.add_function(wrap_pyfunction!(py_update_products_stereochem, m)?)?;
    m.add_function(wrap_pyfunction!(py_reduce_product_to_side_chains, m)?)?;
    m.add_function(wrap_pyfunction!(py_remove_mapping_numbers, m)?)?;
    Ok(())
}
