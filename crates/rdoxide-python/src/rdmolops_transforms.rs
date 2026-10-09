//! Substructure transforms and stereo-bond marking of `rdkit.Chem.rdmolops`
//! (upstream `GraphMol/Wrap/MolOps.cpp`): `DeleteSubstructs`,
//! `ReplaceSubstructs`, `ReplaceSidechains`, `ReplaceCore`,
//! `MurckoDecompose` and `FindPotentialStereoBonds`.

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::graphmol::chem_transforms::{
    delete_substructs, murcko_decompose, replace_core, replace_core_with_match, replace_sidechains,
    replace_substructs,
};
use rdoxide::graphmol::chirality::legacy::find_potential_stereo_bonds;
use rdoxide::graphmol::substruct::MatchVect;

use crate::errors::{Check, IntoPyResult, Site, argument_error, boost_type_name, check_failed};
use crate::mol::PyMol;
use crate::rdmoldescriptors::bool_arg;
use crate::rdmolfiles_write::Args;

/// Removes atoms matching a substructure query from a molecule
///
///   ARGUMENTS:
///
///     - mol: the molecule to be modified
///
///     - query: the molecule to be used as a substructure query
///
///     - onlyFrags: (optional) if this toggle is set, atoms will only be removed if
///       the entire fragment in which they are found is matched by the query.
///       See below for examples.
///       Default value is 0 (remove the atoms whether or not the entire fragment matches)
///
///     - useChirality: (optional) match the substructure query using chirality
///
///   RETURNS: a new molecule with the substructure removed
///
///   NOTES:
///
///     - The original molecule is *not* modified.
///
///   EXAMPLES:
///
///    The following examples substitute SMILES/SMARTS strings for molecules, you'd have
///    to actually use molecules:
///
///     - DeleteSubstructs('CCOC','OC') -> 'CC'
///
///     - DeleteSubstructs('CCOC','OC',1) -> 'CCOC'
///
///     - DeleteSubstructs('CCOCCl.Cl','Cl',1) -> 'CCOCCl'
///
///     - DeleteSubstructs('CCOCCl.Cl','Cl') -> 'CCOC'
#[pyfunction]
#[pyo3(name = "DeleteSubstructs", signature = (mol, query, onlyFrags = None, useChirality = None))]
#[allow(non_snake_case)]
fn py_delete_substructs(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    query: PyRef<'_, PyMol>,
    onlyFrags: Option<&Bound<'_, PyAny>>,
    useChirality: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyMol> {
    let only_frags = bool_arg(onlyFrags, false)?;
    let use_chirality = bool_arg(useChirality, false)?;
    let (mol, query) = (&mol.mol, &query.mol);
    let res = py
        .detach(|| delete_substructs(mol, query, only_frags, use_chirality))
        .py()?;
    Ok(PyMol { mol: res })
}

/// Replaces atoms matching a substructure query in a molecule
///
///   ARGUMENTS:
///
///     - mol: the molecule to be modified
///
///     - query: the molecule to be used as a substructure query
///
///     - replacement: the molecule to be used as the replacement
///
///     - replaceAll: (optional) if this toggle is set, all substructures matching
///       the query will be replaced in a single result, otherwise each result will
///       contain a separate replacement.
///       Default value is False (return multiple replacements)
///     - replacementConnectionPoint: (optional) index of the atom in the replacement that
///       the bond should be made to.
///     - useChirality: (optional) match the substructure query using chirality
///
///   RETURNS: a tuple of new molecules with the substructures replaced removed
///
///   NOTES:
///
///     - The original molecule is *not* modified.
///     - A bond is only formed to the remaining atoms, if any, that were bonded
///       to the first atom in the substructure query. (For finer control over
///       substructure replacement, consider using ChemicalReaction.)
///
///   EXAMPLES:
///
///    The following examples substitute SMILES/SMARTS strings for molecules, you'd have
///    to actually use molecules:
///
///     - ReplaceSubstructs('CCOC','O[CH3]','NC') -> ('CCNC',)
///
///     - ReplaceSubstructs('COCCOC','O[CH3]','NC') -> ('COCCNC','CNCCOC')
///
///     - ReplaceSubstructs('COCCOC','O[CH3]','NC',True) -> ('CNCCNC',)
///
///     - ReplaceSubstructs('COCCOC','O[CH3]','CN',True,1) -> ('CNCCNC',)
///
///     - ReplaceSubstructs('CCOC','[CH3]O','NC') -> ('CC.CN',)
#[pyfunction]
#[pyo3(
    name = "ReplaceSubstructs",
    signature = (mol, query, replacement, replaceAll = None, replacementConnectionPoint = 0, useChirality = None)
)]
#[allow(non_snake_case)]
fn py_replace_substructs<'py>(
    py: Python<'py>,
    mol: PyRef<'py, PyMol>,
    query: PyRef<'py, PyMol>,
    replacement: PyRef<'py, PyMol>,
    replaceAll: Option<&Bound<'py, PyAny>>,
    replacementConnectionPoint: u32,
    useChirality: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let replace_all = bool_arg(replaceAll, false)?;
    let use_chirality = bool_arg(useChirality, false)?;
    if replacementConnectionPoint >= replacement.mol.num_atoms() {
        return Err(check_failed(
            Site {
                check: Check::Precondition,
                file: "Code/GraphMol/ChemTransforms/ChemTransforms.cpp",
                line: 186,
                expr: "replacementConnectionPoint < replacement.getNumAtoms()",
            },
            "bad replacementConnectionPoint",
        ));
    }
    let (mol, query, replacement) = (&mol.mol, &query.mol, &replacement.mol);
    let res = py
        .detach(|| {
            replace_substructs(
                mol,
                query,
                replacement,
                replace_all,
                replacementConnectionPoint,
                use_chirality,
            )
        })
        .py()?;
    let mols = res
        .into_iter()
        .map(|mol| Py::new(py, PyMol { mol }))
        .collect::<PyResult<Vec<_>>>()?;
    PyTuple::new(py, mols)
}

/// Postprocesses the results of a mol.GetSubstructMatches(core) call
/// where mol has explicit Hs and core bears terminal dummy atoms (i.e., R groups).
/// It returns the match with the largest number of non-hydrogen matches to
/// the terminal dummy atoms.
///
///   ARGUMENTS:
///
///     - mol: the molecule GetSubstructMatches was run on
///
///     - core: the molecule used as a substructure query
///   RETURNS: a new molecule with the sidechains removed
///
///   NOTES:
///
///     - The original molecule is *not* modified.
///
///   EXAMPLES:
///
///    The following examples substitute SMILES/SMARTS strings for molecules, you'd have
///    to actually use molecules:
///
///     - ReplaceSidechains('CCC1CCC1','C1CCC1') -> '[Xa]C1CCC1'
///
///     - ReplaceSidechains('CCC1CC1','C1CCC1') -> ''
///
///     - ReplaceSidechains('C1CC2C1CCC2','C1CCC1') -> '[Xa]C1CCC1[Xb]'
#[pyfunction]
#[pyo3(name = "ReplaceSidechains", signature = (mol, coreQuery, useChirality = None))]
#[allow(non_snake_case)]
fn py_replace_sidechains(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    coreQuery: PyRef<'_, PyMol>,
    useChirality: Option<&Bound<'_, PyAny>>,
) -> PyResult<Option<PyMol>> {
    let use_chirality = bool_arg(useChirality, false)?;
    let (mol, core) = (&mol.mol, &coreQuery.mol);
    let res = py
        .detach(|| replace_sidechains(mol, core, use_chirality))
        .py()?;
    Ok(res.map(|mol| PyMol { mol }))
}

/// Do a Murcko decomposition and return the scaffold
#[pyfunction]
#[pyo3(name = "MurckoDecompose", signature = (mol))]
fn py_murcko_decompose(py: Python<'_>, mut mol: PyRefMut<'_, PyMol>) -> PyResult<PyMol> {
    let core = &mut mol.mol;
    let res = py.detach(|| murcko_decompose(core)).py()?;
    Ok(PyMol { mol: res })
}

/// The `(core atom, molecule atom)` pairs of a `ReplaceCore` match given as
/// `(core_atom_idx, molecule_atom_idx)` pairs or, for a match from
/// `GetSubstructMatch`, as the molecule atom indices in core order.
fn match_vect_from_py(matches: &Bound<'_, PyAny>, num_core_atoms: u32) -> PyResult<MatchVect> {
    let length = matches.len()?;
    let mut res = MatchVect::with_capacity(length);
    for i in 0..length {
        let item = matches.get_item(i)?;
        // objects without a length count as a single molecule atom index
        let (v1, v2) = match item.len().unwrap_or(1) {
            1 => {
                if length != num_core_atoms as usize {
                    let entries = if num_core_atoms == 1 {
                        " entry"
                    } else {
                        " entries"
                    };
                    return Err(PyValueError::new_err(format!(
                        "When using input vector of type (molecule_atom_idx,...) supplied core requires {num_core_atoms}{entries}"
                    )));
                }
                (i as i32, item.extract()?)
            }
            2 => (item.get_item(0)?.extract()?, item.get_item(1)?.extract()?),
            _ => {
                return Err(PyValueError::new_err(
                    "Input not a vector of (core_atom_idx,molecule_atom_idx) or (molecule_atom_idx,...) entries",
                ));
            }
        };
        res.push((v1, v2));
    }
    Ok(res)
}

const REPLACE_CORE_QUERY_SIGNATURE: &str = "ReplaceCore(RDKit::ROMol mol, RDKit::ROMol coreQuery, bool replaceDummies=True, bool labelByIndex=False, bool requireDummyMatch=False, bool useChirality=False)";
const REPLACE_CORE_MATCH_SIGNATURE: &str = "ReplaceCore(RDKit::ROMol mol, RDKit::ROMol core, boost::python::api::object matches, bool replaceDummies=True, bool labelByIndex=False, bool requireDummyMatch=False)";

/// Boost.Python's bool converter takes objects whose type has `nb_bool`
/// (numbers), not sequences or `None`.
fn converts_to_bool(obj: &Bound<'_, PyAny>) -> bool {
    obj.is_instance_of::<pyo3::types::PyBool>()
        || obj.is_instance_of::<pyo3::types::PyInt>()
        || obj.is_instance_of::<pyo3::types::PyFloat>()
}

/// Removes the core of a molecule and labels the sidechains with dummy atoms
/// based on the matches indices given in the matching vector matches, or,
/// without matches, on the first match of the query core.
///
/// Calling:
///   ReplaceCore(mol,core,mol.GetSubstructMatch(core))
/// is the same as calling ReplaceCore(mol,core).
///
///   ARGUMENTS:
///
///     - mol: the molecule
///
///     - core: the molecule to be used as a substructure query for recognizing the core
///
///     - matches: a matching vector of the form (molecule_atom_idx, ...) or
///       ((core_atom_idx, molecule_atom_idx), ...) (only with a matches argument)
///
///     - replaceDummies: toggles replacement of atoms that match dummies in the query
///
///     - labelByIndex: toggles labeling the attachment point dummy atoms with
///       the index of the core atom they're attached to.
///
///     - requireDummyMatch: if the molecule has side chains that attach at points not
///       flagged with a dummy, it will be rejected (None is returned)
///
///     - useChirality: (only without matches) use chirality in the substructure matching
///
///   RETURNS: a new molecule with the core removed, or None when the core is
///   not found (or a required dummy match is missing)
///
///   NOTES:
///
///     - The original molecule is *not* modified.
#[pyfunction]
#[pyo3(name = "ReplaceCore", signature = (*args, **kwargs))]
fn py_replace_core(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Option<PyMol>> {
    let a = Args {
        name: "ReplaceCore",
        args,
        kwargs,
    };
    let kw_has = |k: &str| kwargs.is_some_and(|d| d.contains(k).unwrap_or(false));
    let third = a.get(2, "matches")?;
    // the query overload is registered last, so Boost tries it first
    let with_matches =
        kw_has("matches") || kw_has("core") || third.as_ref().is_some_and(|t| !converts_to_bool(t));
    let mol = a.get(0, "mol")?;
    let core = a.get(1, if with_matches { "core" } else { "coreQuery" })?;
    let molecules = mol
        .as_ref()
        .and_then(|m| m.cast::<PyMol>().ok())
        .zip(core.as_ref().and_then(|c| c.cast::<PyMol>().ok()));
    let (mol, core) = match molecules {
        Some(m) => m,
        None => {
            let got: Vec<String> = args.iter().map(|o| boost_type_name(&o)).collect();
            return Err(argument_error(
                py,
                "ReplaceCore",
                &got,
                &format!("{REPLACE_CORE_QUERY_SIGNATURE}\n    {REPLACE_CORE_MATCH_SIGNATURE}"),
            ));
        }
    };
    let (mol, core) = (&mol.borrow().mol, &core.borrow().mol);
    if with_matches {
        let matches = third.ok_or_else(|| {
            PyTypeError::new_err("ReplaceCore() missing required argument 'matches'")
        })?;
        let match_v = match_vect_from_py(&matches, core.num_atoms())?;
        let replace_dummies = a.flag(3, "replaceDummies", true)?;
        let label_by_index = a.flag(4, "labelByIndex", false)?;
        let require_dummy_match = a.flag(5, "requireDummyMatch", false)?;
        let res = py
            .detach(|| {
                replace_core_with_match(
                    mol,
                    core,
                    &match_v,
                    replace_dummies,
                    label_by_index,
                    require_dummy_match,
                )
            })
            .py()?;
        Ok(res.map(|mol| PyMol { mol }))
    } else {
        let replace_dummies = a.flag(2, "replaceDummies", true)?;
        let label_by_index = a.flag(3, "labelByIndex", false)?;
        let require_dummy_match = a.flag(4, "requireDummyMatch", false)?;
        let use_chirality = a.flag(5, "useChirality", false)?;
        let res = py
            .detach(|| {
                replace_core(
                    mol,
                    core,
                    replace_dummies,
                    label_by_index,
                    require_dummy_match,
                    use_chirality,
                )
            })
            .py()?;
        Ok(res.map(|mol| PyMol { mol }))
    }
}

/// Find bonds that could be cis/trans in a molecule and mark them as Bond.STEREOANY.
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
///     - cleanIt: (optional) if this option is set to true, any previous marking of _CIPCode
///                on the bond is cleared - otherwise it is left untouched
#[pyfunction]
#[pyo3(name = "FindPotentialStereoBonds", signature = (mol, cleanIt = None))]
#[allow(non_snake_case)]
fn py_find_potential_stereo_bonds(
    mut mol: PyRefMut<'_, PyMol>,
    cleanIt: Option<&Bound<'_, PyAny>>,
) -> PyResult<()> {
    let clean_it = bool_arg(cleanIt, false)?;
    find_potential_stereo_bonds(&mut mol.mol, clean_it).py()
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_delete_substructs, m)?)?;
    m.add_function(wrap_pyfunction!(py_replace_substructs, m)?)?;
    m.add_function(wrap_pyfunction!(py_replace_sidechains, m)?)?;
    m.add_function(wrap_pyfunction!(py_replace_core, m)?)?;
    m.add_function(wrap_pyfunction!(py_murcko_decompose, m)?)?;
    m.add_function(wrap_pyfunction!(py_find_potential_stereo_bonds, m)?)?;
    Ok(())
}
