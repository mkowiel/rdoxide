//! `rdoxide.Chem.rdmolfiles`: SMILES writing — port of the corresponding parts
//! of upstream `Code/GraphMol/Wrap/rdmolfiles.cpp`.

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rdoxide::graphmol::smiles_parse::cx_write::{
    RestoreBondDirOption, cx_fields, mol_fragment_to_cx_smiles, mol_to_cx_smiles,
};
use rdoxide::graphmol::smiles_parse::write::{
    SmilesWriteParams, mol_fragment_to_smiles, mol_to_random_smiles_vect, mol_to_smiles,
};

use rdoxide::graphmol::canon::canonicalize_enhanced_stereo;

use crate::enums::add_enum;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Parameters controlling SMILES writing
#[pyclass(
    name = "SmilesWriteParams",
    module = "rdoxide.Chem.rdmolfiles",
    get_all,
    set_all,
    skip_from_py_object
)]
#[derive(Clone)]
#[allow(non_snake_case)]
pub struct PySmilesWriteParams {
    /// include stereochemistry and isotope information
    doIsomericSmiles: bool,
    /// kekulize the molecule before generating the SMILES and output single/double bonds. NOTE that the output is not canonical and that this will thrown an exception if the molecule cannot be kekulized
    doKekule: bool,
    /// generate canonical SMILES
    canonical: bool,
    /// chiral centers are removed if they have duplicate sidechains
    cleanStereo: bool,
    /// include symbols for all bonds
    allBondsExplicit: bool,
    /// provide hydrogen counts for every atom
    allHsExplicit: bool,
    /// randomize the output order. The resulting SMILES is not canonical
    doRandom: bool,
    /// make sure the SMILES starts at the specified atom. The resulting SMILES is not canonical
    rootedAtAtom: i32,
    /// include the RDKit extension for dative bonds. Otherwise dative bonds will be written as single bonds
    includeDativeBonds: bool,
    /// ignore atom map numbers when canonicalizing the molecule
    ignoreAtomMapNumbers: bool,
}

impl PySmilesWriteParams {
    pub(crate) fn to_core(&self) -> SmilesWriteParams {
        SmilesWriteParams {
            do_isomeric_smiles: self.doIsomericSmiles,
            do_kekule: self.doKekule,
            canonical: self.canonical,
            clean_stereo: self.cleanStereo,
            all_bonds_explicit: self.allBondsExplicit,
            all_hs_explicit: self.allHsExplicit,
            do_random: self.doRandom,
            rooted_at_atom: self.rootedAtAtom,
            include_dative_bonds: self.includeDativeBonds,
            ignore_atom_map_numbers: self.ignoreAtomMapNumbers,
        }
    }
}

#[pymethods]
impl PySmilesWriteParams {
    #[new]
    fn py_new() -> Self {
        let d = SmilesWriteParams::default();
        PySmilesWriteParams {
            doIsomericSmiles: d.do_isomeric_smiles,
            doKekule: d.do_kekule,
            canonical: d.canonical,
            cleanStereo: d.clean_stereo,
            allBondsExplicit: d.all_bonds_explicit,
            allHsExplicit: d.all_hs_explicit,
            doRandom: d.do_random,
            rootedAtAtom: d.rooted_at_atom,
            includeDativeBonds: d.include_dative_bonds,
            ignoreAtomMapNumbers: d.ignore_atom_map_numbers,
        }
    }
}

/// Positional/keyword argument lookup for the overloaded functions.
pub(crate) struct Args<'a, 'py> {
    pub(crate) name: &'static str,
    pub(crate) args: &'a Bound<'py, PyTuple>,
    pub(crate) kwargs: Option<&'a Bound<'py, PyDict>>,
}

impl<'py> Args<'_, 'py> {
    pub(crate) fn get(&self, pos: usize, key: &str) -> PyResult<Option<Bound<'py, PyAny>>> {
        let kw = match self.kwargs {
            Some(k) => k.get_item(key)?,
            None => None,
        };
        let p = if pos < self.args.len() {
            Some(self.args.get_item(pos)?)
        } else {
            None
        };
        match (p, kw) {
            (Some(_), Some(_)) => Err(PyTypeError::new_err(format!(
                "{}() got multiple values for argument '{key}'",
                self.name
            ))),
            (a, b) => Ok(a.or(b)),
        }
    }
    pub(crate) fn flag(&self, pos: usize, key: &str, default: bool) -> PyResult<bool> {
        self.get(pos, key)?.map_or(Ok(default), |v| v.is_truthy())
    }
    pub(crate) fn int(&self, pos: usize, key: &str, default: i32) -> PyResult<i32> {
        self.get(pos, key)?.map_or(Ok(default), |v| v.extract())
    }
}

fn run_mol_to_smiles(
    py: Python<'_>,
    mol: &Bound<'_, PyAny>,
    params: SmilesWriteParams,
) -> PyResult<String> {
    let mol = mol.cast::<PyMol>()?;
    let mut m = mol.borrow_mut();
    let core = &mut m.mol;
    py.detach(|| mol_to_smiles(core, &params)).py()
}

/// Returns the canonical SMILES string for a molecule
///   ARGUMENTS:
///
///     - mol: the molecule
///     - isomericSmiles: (optional) include information about stereochemistry in
///       the SMILES.  Defaults to true.
///     - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds) in
///       the SMILES.  Defaults to false.
///     - rootedAtAtom: (optional) if non-negative, this forces the SMILES
///       to start at a particular atom. Defaults to -1.  If not -1, overrides
///       canonical setting.
///     - canonical: (optional) if false no attempt will be made to canonicalize
///       the molecule. Defaults to true.
///     - allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated
///       in the output SMILES. Defaults to false.
///     - allHsExplicit: (optional) if true, all H counts will be explicitly indicated
///       in the output SMILES. Defaults to false.
///     - doRandom: (optional) if true, randomize the traversal of the molecule graph,
///       so we can generate random smiles. Defaults to false.  If true, overrides
///       canonical setting.
///     - ignoreAtomMapNumbers (optional) if true, ignores any atom map numbers when
///       canonicalizing the molecule
///
///   RETURNS:
///
///     a string
///
/// Alternatively: MolToSmiles(mol, params) with a SmilesWriteParams object.
#[pyfunction]
#[pyo3(name = "MolToSmiles", signature = (*args, **kwargs))]
fn py_mol_to_smiles(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    let a = Args {
        name: "MolToSmiles",
        args,
        kwargs,
    };
    let mol = a
        .get(0, "mol")?
        .ok_or_else(|| PyTypeError::new_err("MolToSmiles() missing required argument 'mol'"))?;
    let second = match a.get(1, "params")? {
        Some(p) => p
            .cast::<PySmilesWriteParams>()
            .ok()
            .map(|p| p.borrow().to_core()),
        None => None,
    };
    if let Some(params) = second {
        return run_mol_to_smiles(py, &mol, params);
    }
    let params = SmilesWriteParams {
        do_isomeric_smiles: a.flag(1, "isomericSmiles", true)?,
        do_kekule: a.flag(2, "kekuleSmiles", false)?,
        rooted_at_atom: a.int(3, "rootedAtAtom", -1)?,
        canonical: a.flag(4, "canonical", true)?,
        all_bonds_explicit: a.flag(5, "allBondsExplicit", false)?,
        all_hs_explicit: a.flag(6, "allHsExplicit", false)?,
        do_random: a.flag(7, "doRandom", false)?,
        ignore_atom_map_numbers: a.flag(8, "ignoreAtomMapNumbers", false)?,
        ..SmilesWriteParams::default()
    };
    run_mol_to_smiles(py, &mol, params)
}

/// Returns the CXSMILES string for a molecule
///   ARGUMENTS:
///
///     - mol: the molecule
///     - isomericSmiles: (optional) include information about stereochemistry in
///       the SMILES.  Defaults to true.
///     - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds) in
///       the SMILES.  Defaults to false.
///     - rootedAtAtom: (optional) if non-negative, this forces the SMILES
///       to start at a particular atom. Defaults to -1.
///     - canonical: (optional) if false no attempt will be made to canonicalize
///       the molecule. Defaults to true.
///     - allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated
///       in the output SMILES. Defaults to false.
///     - allHsExplicit: (optional) if true, all H counts will be explicitly indicated
///       in the output SMILES. Defaults to false.
///     - doRandom: (optional) if true, randomizes the traversal of the molecule graph,
///       so we can generate random smiles. Defaults to false.  If true, overrides
///       canonical setting.
///
///   RETURNS:
///
///     a string
///
/// Alternatively: MolToCXSmiles(mol, params, flags=CXSmilesFields.CX_ALL,
/// restoreBondDirs=RestoreBondDirOption.RestoreBondDirOptionClear)
#[pyfunction]
#[pyo3(name = "MolToCXSmiles", signature = (*args, **kwargs))]
fn py_mol_to_cx_smiles(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    let a = Args {
        name: "MolToCXSmiles",
        args,
        kwargs,
    };
    let mol = a
        .get(0, "mol")?
        .ok_or_else(|| PyTypeError::new_err("MolToCXSmiles() missing required argument 'mol'"))?;
    let mol = mol.cast::<PyMol>()?.borrow().mol.clone();
    let (params, flags, restore) = match a.get(1, "params")? {
        Some(p) if p.cast::<PySmilesWriteParams>().is_ok() => {
            let params = p.cast::<PySmilesWriteParams>()?.borrow().to_core();
            let flags = a
                .get(2, "flags")?
                .map_or(Ok(cx_fields::ALL), |v| v.extract::<u32>())?;
            let restore = match a.get(3, "restoreBondDirs")? {
                Some(v) if v.extract::<i64>()? == 0 => RestoreBondDirOption::True,
                _ => RestoreBondDirOption::Clear,
            };
            (params, flags, restore)
        }
        _ => (
            SmilesWriteParams {
                do_isomeric_smiles: a.flag(1, "isomericSmiles", true)?,
                do_kekule: a.flag(2, "kekuleSmiles", false)?,
                rooted_at_atom: a.int(3, "rootedAtAtom", -1)?,
                canonical: a.flag(4, "canonical", true)?,
                all_bonds_explicit: a.flag(5, "allBondsExplicit", false)?,
                all_hs_explicit: a.flag(6, "allHsExplicit", false)?,
                do_random: a.flag(7, "doRandom", false)?,
                ..SmilesWriteParams::default()
            },
            cx_fields::ALL,
            RestoreBondDirOption::Clear,
        ),
    };
    py.detach(|| mol_to_cx_smiles(&mol, &params, flags, restore))
        .py()
}

/// `pythonObjectToVect<int>(obj, maxV)`
fn int_vect(obj: Option<Bound<'_, PyAny>>, max_v: usize) -> PyResult<Option<Vec<i32>>> {
    let Some(o) = obj.filter(|o| !o.is_none() && !o.extract::<i64>().is_ok_and(|v| v == 0)) else {
        return Ok(None);
    };
    let mut res = Vec::new();
    for item in o.try_iter()? {
        let v: i32 = item?.extract()?;
        if v < 0 || v as usize >= max_v {
            return Err(PyValueError::new_err(
                "list element larger than allowed value",
            ));
        }
        res.push(v);
    }
    Ok(Some(res))
}

fn str_vect(obj: Option<Bound<'_, PyAny>>) -> PyResult<Option<Vec<String>>> {
    let Some(o) = obj.filter(|o| !o.is_none() && !o.extract::<i64>().is_ok_and(|v| v == 0)) else {
        return Ok(None);
    };
    o.try_iter()?
        .map(|i| i?.extract::<String>())
        .collect::<PyResult<Vec<_>>>()
        .map(Some)
}

/// Shared body of `MolFragmentToSmiles` and `MolFragmentToCXSmiles`.
fn fragment_to_smiles(
    py: Python<'_>,
    name: &'static str,
    cx: bool,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    let a = Args { name, args, kwargs };
    let mol_obj = a
        .get(0, "mol")?
        .ok_or_else(|| PyTypeError::new_err(format!("{name}() missing required argument 'mol'")))?;
    let mol = mol_obj.cast::<PyMol>()?;
    let mut m = mol.borrow_mut();
    let (na, nb) = (m.mol.num_atoms() as usize, m.mol.num_bonds() as usize);
    let params_obj = a.get(1, "params")?;
    let explicit_params = params_obj
        .as_ref()
        .and_then(|p| p.cast::<PySmilesWriteParams>().ok())
        .map(|p| p.borrow().to_core());
    let (params, off) = match explicit_params {
        Some(p) => (p, 1),
        None => (
            SmilesWriteParams {
                do_isomeric_smiles: a.flag(5, "isomericSmiles", true)?,
                do_kekule: a.flag(6, "kekuleSmiles", false)?,
                rooted_at_atom: a.int(7, "rootedAtAtom", -1)?,
                canonical: a.flag(8, "canonical", true)?,
                all_bonds_explicit: a.flag(9, "allBondsExplicit", false)?,
                all_hs_explicit: a.flag(10, "allHsExplicit", false)?,
                ..SmilesWriteParams::default()
            },
            0,
        ),
    };
    let atoms = int_vect(a.get(1 + off, "atomsToUse")?, na)?
        .filter(|v| !v.is_empty())
        .ok_or_else(|| PyValueError::new_err("atomsToUse must not be empty"))?;
    let bonds = int_vect(a.get(2 + off, "bondsToUse")?, nb)?;
    let asym = str_vect(a.get(3 + off, "atomSymbols")?)?;
    let bsym = str_vect(a.get(4 + off, "bondSymbols")?)?;
    if asym.as_ref().is_some_and(|s| s.len() != na) {
        return Err(PyValueError::new_err(
            "length of atom symbol list != number of atoms",
        ));
    }
    if bsym.as_ref().is_some_and(|s| s.len() != nb) {
        return Err(PyValueError::new_err(
            "length of bond symbol list != number of bonds",
        ));
    }
    let core = &mut m.mol;
    let write = if cx {
        mol_fragment_to_cx_smiles
    } else {
        mol_fragment_to_smiles
    };
    py.detach(|| {
        write(
            core,
            &params,
            &atoms,
            bonds.as_deref(),
            asym.as_deref(),
            bsym.as_deref(),
        )
    })
    .py()
}

/// Returns the canonical SMILES string for a fragment of a molecule
///   ARGUMENTS:
///
///     - mol: the molecule
///     - atomsToUse : a list of atoms to include in the fragment
///     - bondsToUse : (optional) a list of bonds to include in the fragment
///       if not provided, all bonds between the atoms provided
///       will be included.
///     - atomSymbols : (optional) a list with the symbols to use for the atoms
///       in the SMILES. This should have be mol.GetNumAtoms() long.
///     - bondSymbols : (optional) a list with the symbols to use for the bonds
///       in the SMILES. This should have be mol.GetNumBonds() long.
///     - isomericSmiles: (optional) include information about stereochemistry in
///       the SMILES.  Defaults to true.
///     - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds) in
///       the SMILES.  Defaults to false.
///     - rootedAtAtom: (optional) if non-negative, this forces the SMILES
///       to start at a particular atom. Defaults to -1.  If not -1, over-rides
///       setting for canonical.
///     - canonical: (optional) if false no attempt will be made to canonicalize
///       the molecule. Defaults to true.
///     - allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated
///       in the output SMILES. Defaults to false.
///     - allHsExplicit: (optional) if true, all H counts will be explicitly indicated
///       in the output SMILES. Defaults to false.
///
///   RETURNS:
///
///     a string
///
/// Alternatively: MolFragmentToSmiles(mol, params, atomsToUse, bondsToUse=0,
/// atomSymbols=0, bondSymbols=0) with a SmilesWriteParams object.
#[pyfunction]
#[pyo3(name = "MolFragmentToSmiles", signature = (*args, **kwargs))]
fn py_mol_fragment_to_smiles(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    fragment_to_smiles(py, "MolFragmentToSmiles", false, args, kwargs)
}

/// Returns the CXSMILES string for a fragment of a molecule (same arguments
/// as MolFragmentToSmiles)
#[pyfunction]
#[pyo3(name = "MolFragmentToCXSmiles", signature = (*args, **kwargs))]
fn py_mol_fragment_to_cx_smiles(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    fragment_to_smiles(py, "MolFragmentToCXSmiles", true, args, kwargs)
}

/// returns a list of SMILES generated using the randomSmiles algorithm
#[pyfunction]
#[pyo3(name = "MolToRandomSmilesVect", signature = (mol, numSmiles, randomSeed = 0, isomericSmiles = true, kekuleSmiles = false, allBondsExplicit = false, allHsExplicit = false))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn py_mol_to_random_smiles_vect<'py>(
    py: Python<'py>,
    mut mol: PyRefMut<'_, PyMol>,
    numSmiles: u32,
    randomSeed: u32,
    isomericSmiles: bool,
    kekuleSmiles: bool,
    allBondsExplicit: bool,
    allHsExplicit: bool,
) -> PyResult<Bound<'py, PyList>> {
    let params = SmilesWriteParams {
        do_isomeric_smiles: isomericSmiles,
        do_kekule: kekuleSmiles,
        all_bonds_explicit: allBondsExplicit,
        all_hs_explicit: allHsExplicit,
        ..SmilesWriteParams::default()
    };
    let core = &mut mol.mol;
    let res = py
        .detach(|| mol_to_random_smiles_vect(core, numSmiles, randomSeed, &params))
        .py()?;
    PyList::new(py, res)
}

/// Rewrites the enhanced-stereo groups of a molecule into canonical form.
#[pyfunction]
#[pyo3(name = "CanonicalizeEnhancedStereo", signature = (mol))]
fn py_canonicalize_enhanced_stereo(mut mol: PyRefMut<'_, PyMol>) -> PyResult<()> {
    canonicalize_enhanced_stereo(&mut mol.mol, None).py()
}

/// Returns a SMARTS string for a molecule
///   ARGUMENTS:
///
///     - mol: the molecule
///     - isomericSmiles: (optional) include information about stereochemistry in
///       the SMARTS.  Defaults to true.
///     - rootedAtAtom: (optional) the atom index to start the SMARTS from.
///       Defaults to -1.
///
///   RETURNS:
///
///     a string
///
/// Alternatively: MolToSmarts(mol, params) with a SmilesWriteParams object.
#[pyfunction]
#[pyo3(name = "MolToSmarts", signature = (*args, **kwargs))]
fn py_mol_to_smarts(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    use rdoxide::graphmol::smiles_parse::smarts_write::mol_to_smarts;
    let a = Args {
        name: "MolToSmarts",
        args,
        kwargs,
    };
    let mol = a
        .get(0, "mol")?
        .ok_or_else(|| PyTypeError::new_err("MolToSmarts() missing required argument 'mol'"))?;
    let params = match a.get(1, "params")?.and_then(|p| {
        p.cast::<PySmilesWriteParams>()
            .ok()
            .map(|p| p.borrow().to_core())
    }) {
        Some(p) => p,
        None => SmilesWriteParams {
            do_isomeric_smiles: a.flag(1, "isomericSmiles", true)?,
            rooted_at_atom: a.int(2, "rootedAtAtom", -1)?,
            ..SmilesWriteParams::default()
        },
    };
    let mol = mol.cast::<PyMol>()?;
    let mut m = mol.borrow_mut();
    let core = &mut m.mol;
    py.detach(|| mol_to_smarts(core, &params)).py()
}

/// Returns a SMARTS string for a fragment of a molecule
///   ARGUMENTS:
///
///     - mol: the molecule
///     - atomsToUse: indices of atoms to include in the SMARTS string
///     - bondsToUse: indices of bonds to include in the SMARTS string (optional)
///     - isomericSmarts: (optional) include information about stereochemistry in
///       the SMARTS.  Defaults to true.
///
///   RETURNS:
///
///     a string
#[pyfunction]
#[pyo3(name = "MolFragmentToSmarts", signature = (mol, atomsToUse, bondsToUse = None, isomericSmarts = true))]
#[allow(non_snake_case)]
fn py_mol_fragment_to_smarts(
    py: Python<'_>,
    mut mol: PyRefMut<'_, PyMol>,
    atomsToUse: Bound<'_, PyAny>,
    bondsToUse: Option<Bound<'_, PyAny>>,
    isomericSmarts: bool,
) -> PyResult<String> {
    use rdoxide::graphmol::smiles_parse::smarts_write::mol_fragment_to_smarts;
    let (na, nb) = (mol.mol.num_atoms() as usize, mol.mol.num_bonds() as usize);
    let atoms = int_vect(Some(atomsToUse), na)?
        .filter(|v| !v.is_empty())
        .ok_or_else(|| PyValueError::new_err("atomsToUse argument must be non-empty"))?;
    let bonds = int_vect(bondsToUse, nb)?;
    let params = SmilesWriteParams {
        do_isomeric_smiles: isomericSmarts,
        ..SmilesWriteParams::default()
    };
    let core = &mut mol.mol;
    py.detach(|| mol_fragment_to_smarts(core, &params, &atoms, bonds.as_deref()))
        .py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_mol_to_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_fragment_to_smarts, m)?)?;
    m.add_class::<PySmilesWriteParams>()?;
    m.add_function(wrap_pyfunction!(py_mol_to_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_cx_smiles, m)?)?;
    add_enum(
        m,
        "CXSmilesFields",
        &[
            ("CX_NONE", i64::from(cx_fields::NONE)),
            ("CX_ATOM_LABELS", i64::from(cx_fields::ATOM_LABELS)),
            ("CX_MOLFILE_VALUES", i64::from(cx_fields::MOLFILE_VALUES)),
            ("CX_COORDS", i64::from(cx_fields::COORDS)),
            ("CX_RADICALS", i64::from(cx_fields::RADICALS)),
            ("CX_ATOM_PROPS", i64::from(cx_fields::ATOM_PROPS)),
            ("CX_LINKNODES", i64::from(cx_fields::LINKNODES)),
            ("CX_ENHANCEDSTEREO", i64::from(cx_fields::ENHANCEDSTEREO)),
            ("CX_SGROUPS", i64::from(cx_fields::SGROUPS)),
            ("CX_POLYMER", i64::from(cx_fields::POLYMER)),
            ("CX_BOND_CFG", i64::from(cx_fields::BOND_CFG)),
            (
                "CX_BOND_ATROPISOMER",
                i64::from(cx_fields::BOND_ATROPISOMER),
            ),
            (
                "CX_COORDINATE_BONDS",
                i64::from(cx_fields::COORDINATE_BONDS),
            ),
            ("CX_HYDROGEN_BONDS", i64::from(cx_fields::HYDROGEN_BONDS)),
            ("CX_ZERO_BONDS", i64::from(cx_fields::ZERO_BONDS)),
            ("CX_ALL", i64::from(cx_fields::ALL)),
            ("CX_ALL_BUT_COORDS", i64::from(cx_fields::ALL_BUT_COORDS)),
        ],
        None,
    )?;
    add_enum(
        m,
        "RestoreBondDirOption",
        &[
            ("RestoreBondDirOptionClear", 1),
            ("RestoreBondDirOptionTrue", 0),
        ],
        None,
    )?;
    m.add_function(wrap_pyfunction!(py_mol_fragment_to_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_fragment_to_cx_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(py_mol_to_random_smiles_vect, m)?)?;
    m.add_function(wrap_pyfunction!(py_canonicalize_enhanced_stereo, m)?)?;
    Ok(())
}
