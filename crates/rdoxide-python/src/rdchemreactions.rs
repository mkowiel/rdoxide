//! `rdkit.Chem.rdChemReactions` (upstream `ChemReactions/Wrap/rdChemReactions.cpp`):
//! chemical reactions defined by templates, and running them on molecules.
//!
//! Template accessors (`GetReactantTemplate`, `GetReactants`, ...) return
//! copies where upstream returns references into the reaction. The MDL/Marvin/
//! CDXML/PNG readers and writers, pickling, fingerprints, reaction depiction,
//! preprocessing and sanitization are not provided yet.

use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyList, PyTuple};
use rdoxide::graphmol::chem_reactions::daylight::{
    ReactionSmartsParserParams, reaction_from_smarts, reaction_from_smiles,
};
use rdoxide::graphmol::chem_reactions::mdl::{
    reaction_from_rxn_block, reaction_from_rxn_file, reaction_from_rxn_mol, reaction_to_rxn_block,
    reaction_to_rxn_mol, reaction_to_v3k_rxn_block,
};
use rdoxide::graphmol::chem_reactions::writer::{reaction_to_smarts, reaction_to_smiles};
use rdoxide::graphmol::chem_reactions::{self as cr, ChemicalReaction, runner};
use rdoxide::graphmol::file_parsers::MolFileParserParams;
use rdoxide::graphmol::romol::ROMol;
use rdoxide::graphmol::smiles_parse::write::SmilesWriteParams;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::rdmolfiles_write::PySmilesWriteParams;
use crate::substruct::PySubstructMatchParameters;

/// A class for storing and applying chemical reactions.
///
/// Sample Usage:
///   >>> from rdkit import Chem
///   >>> from rdkit.Chem import rdChemReactions
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
    module = "rdkit.Chem.rdChemReactions",
    skip_from_py_object
)]
pub struct PyChemicalReaction {
    rxn: ChemicalReaction,
    /// The object handed out by `GetSubstructParams`; it stays live, so edits
    /// to it must reach `rxn` before the reaction is used.
    params: Option<Py<PySubstructMatchParameters>>,
}

impl PyChemicalReaction {
    fn new(rxn: ChemicalReaction) -> Self {
        PyChemicalReaction { rxn, params: None }
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

    /// Copies the parameters edited through `GetSubstructParams` into the reaction.
    fn sync_params(&mut self, py: Python<'_>) {
        if let Some(p) = &self.params {
            *self.rxn.substruct_params_mut() = p.borrow(py).base_params();
        }
    }

    /// Syncs the parameters and initializes the reaction if that has not happened yet.
    fn prepare(&mut self, py: Python<'_>) {
        self.sync_params(py);
        if !self.rxn.is_initialized() {
            self.rxn.init_reactant_matchers(false);
        }
    }
}

fn template_index(which: u32, n: usize) -> PyResult<usize> {
    if which as usize >= n {
        return Err(PyValueError::new_err("requested template index too high"));
    }
    Ok(which as usize)
}

fn mol_ref<'a>(obj: &'a Bound<'_, PyAny>) -> PyResult<PyRef<'a, PyMol>> {
    Ok(obj.extract::<PyRef<'_, PyMol>>()?)
}

fn products_to_py(py: Python<'_>, sets: Vec<Vec<ROMol>>) -> PyResult<Py<PyTuple>> {
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

fn mols_to_py(py: Python<'_>, mols: &[ROMol]) -> PyResult<Py<PyTuple>> {
    let items = mols
        .iter()
        .map(|mol| Py::new(py, PyMol { mol: mol.clone() }))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyTuple::new(py, items)?.unbind())
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
                if let Ok(r) = o.extract::<PyRef<'_, PyChemicalReaction>>() {
                    return Ok(PyChemicalReaction::new(r.rxn.clone()));
                }
                Err(PyNotImplementedError::new_err(
                    "constructing a ChemicalReaction from a binary string is not supported yet",
                ))
            }
        }
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
    fn add_reactant_template(&mut self, mol: &Bound<'_, PyMol>) -> usize {
        self.rxn.add_reactant_template(mol.borrow().mol.clone())
    }

    /// adds a product (a Molecule)
    #[pyo3(name = "AddProductTemplate")]
    fn add_product_template(&mut self, mol: &Bound<'_, PyMol>) -> usize {
        self.rxn.add_product_template(mol.borrow().mol.clone())
    }

    /// adds a agent (a Molecule)
    #[pyo3(name = "AddAgentTemplate")]
    fn add_agent_template(&mut self, mol: &Bound<'_, PyMol>) -> usize {
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
        let mut removed = Vec::new();
        let target = targetList.filter(|t| !t.is_none()).map(|_| &mut removed);
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
        let mut removed = Vec::new();
        let target = targetList.filter(|t| !t.is_none()).map(|_| &mut removed);
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
        let mut removed = Vec::new();
        let target = targetList.filter(|t| !t.is_none()).map(|_| &mut removed);
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
    }

    /// checks if the reaction is ready for use
    #[pyo3(name = "IsInitialized")]
    fn is_initialized(&self) -> bool {
        self.rxn.is_initialized()
    }

    /// checks the reaction for potential problems, returns (numWarnings,numErrors)
    #[pyo3(name = "Validate", signature = (silent = false))]
    fn validate(&mut self, silent: bool) -> (u32, u32) {
        let (_, warnings, errors) = self.rxn.validate(silent);
        (warnings, errors)
    }

    /// returns one of our product templates
    #[pyo3(name = "GetProductTemplate")]
    fn product_template(&self, which: u32) -> PyResult<PyMol> {
        let i = template_index(which, self.rxn.num_product_templates())?;
        Ok(PyMol {
            mol: self.rxn.products()[i].clone(),
        })
    }

    /// returns one of our reactant templates
    #[pyo3(name = "GetReactantTemplate")]
    fn reactant_template(&self, which: u32) -> PyResult<PyMol> {
        let i = template_index(which, self.rxn.num_reactant_templates())?;
        Ok(PyMol {
            mol: self.rxn.reactants()[i].clone(),
        })
    }

    /// returns one of our agent templates
    #[pyo3(name = "GetAgentTemplate")]
    fn agent_template(&self, which: u32) -> PyResult<PyMol> {
        let i = template_index(which, self.rxn.num_agent_templates())?;
        Ok(PyMol {
            mol: self.rxn.agents()[i].clone(),
        })
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
    fn reacting_atoms(&self, py: Python<'_>, mappedAtomsOnly: bool) -> PyResult<Py<PyTuple>> {
        let atoms = cr::get_reacting_atoms(&self.rxn, mappedAtomsOnly).py()?;
        let rows = atoms
            .into_iter()
            .map(|a| PyTuple::new(py, a))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyTuple::new(py, rows)?.unbind())
    }

    /// get the reactant templates
    #[pyo3(name = "GetReactants")]
    fn reactants(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        mols_to_py(py, self.rxn.reactants())
    }

    /// get the product templates
    #[pyo3(name = "GetProducts")]
    fn products(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        mols_to_py(py, self.rxn.products())
    }

    /// get the agent templates
    #[pyo3(name = "GetAgents")]
    fn agents(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        mols_to_py(py, self.rxn.agents())
    }

    /// get the parameter object controlling the substructure matching
    #[pyo3(name = "GetSubstructParams")]
    fn substruct_params(&mut self, py: Python<'_>) -> PyResult<Py<PySubstructMatchParameters>> {
        if self.params.is_none() {
            self.params = Some(Py::new(py, PySubstructMatchParameters::new_default())?);
        }
        Ok(self
            .params
            .as_ref()
            .map(|p| p.clone_ref(py))
            .expect("set above"))
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
/// see the documentation for rdkit.Chem.MolFromSmiles for an explanation
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
/// see the documentation for rdkit.Chem.MolFromSmiles for an explanation
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
    reaction_to_smarts(&reaction.borrow().rxn, &params).py()
}

/// construct a reaction SMILES string for a ChemicalReaction
#[pyfunction]
#[pyo3(name = "ReactionToSmiles", signature = (reaction, params = None))]
fn py_reaction_to_smiles(
    reaction: &Bound<'_, PyChemicalReaction>,
    params: Option<&Bound<'_, PyAny>>,
) -> PyResult<String> {
    let params = write_params(params, true, true)?;
    reaction_to_smiles(&reaction.borrow().rxn, &params).py()
}

/// construct a molecule for a ChemicalReaction with RXN role property set
#[pyfunction]
#[pyo3(name = "ReactionToMolecule")]
fn py_reaction_to_molecule(reaction: &Bound<'_, PyChemicalReaction>) -> PyResult<PyMol> {
    Ok(PyMol {
        mol: reaction_to_rxn_mol(&reaction.borrow().rxn).py()?,
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

/// construct a string in MDL rxn format for a ChemicalReaction
#[pyfunction]
#[pyo3(name = "ReactionToRxnBlock", signature = (reaction, separateAgents = false, forceV3000 = false))]
#[allow(non_snake_case)]
fn py_reaction_to_rxn_block(
    reaction: &Bound<'_, PyChemicalReaction>,
    separateAgents: bool,
    forceV3000: bool,
) -> PyResult<String> {
    reaction_to_rxn_block(&reaction.borrow().rxn, separateAgents, forceV3000).py()
}

/// construct a string in MDL v3000 rxn format for a ChemicalReaction
#[pyfunction]
#[pyo3(name = "ReactionToV3KRxnBlock", signature = (reaction, separateAgents = false))]
#[allow(non_snake_case)]
fn py_reaction_to_v3k_rxn_block(
    reaction: &Bound<'_, PyChemicalReaction>,
    separateAgents: bool,
) -> PyResult<String> {
    reaction_to_v3k_rxn_block(&reaction.borrow().rxn, separateAgents).py()
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
    cr::has_reaction_atom_mapping(&rxn.borrow().rxn)
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
    let py = reaction.py();
    reaction.borrow_mut().sync_params(py);
    cr::has_reaction_substruct_match(
        &reaction.borrow().rxn,
        &queryReaction.borrow().rxn,
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
            let py = reaction.py();
            reaction.borrow_mut().sync_params(py);
            $core(&reaction.borrow().rxn, &queryReaction.borrow().rxn).py()
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
    cr::update_products_stereochem(&mut reaction.borrow_mut().rxn).py()
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
    cr::remove_mapping_numbers_from_reactions(&mut reaction.borrow_mut().rxn);
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing classes and functions for working with chemical reactions.",
    )?;
    m.add_class::<PyChemicalReaction>()?;
    m.add_function(wrap_pyfunction!(py_reaction_from_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_from_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_to_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_to_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_to_molecule, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_from_molecule, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_from_rxn_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_from_rxn_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_to_rxn_block, m)?)?;
    m.add_function(wrap_pyfunction!(py_reaction_to_v3k_rxn_block, m)?)?;
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
