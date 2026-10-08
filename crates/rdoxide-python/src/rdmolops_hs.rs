//! `rdkit.Chem.rdmolops`: hydrogen removal and stereo perception — port of
//! the corresponding parts of upstream `Code/GraphMol/Wrap/MolOps.cpp`.

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::graphmol::chirality::{bond_dirs, legacy};
use rdoxide::graphmol::mol_ops::add_hs::{AddHsParameters, add_hs};
use rdoxide::graphmol::mol_ops::{RemoveHsParameters, remove_hs};
use rdoxide::graphmol::romol::ROMol;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Parameters controlling which Hs are removed.
#[pyclass(
    name = "RemoveHsParameters",
    module = "rdkit.Chem.rdmolops",
    get_all,
    set_all,
    skip_from_py_object
)]
#[derive(Clone, Default)]
#[allow(non_snake_case)]
pub struct PyRemoveHsParameters {
    /// hydrogens that have no bonds
    removeDegreeZero: bool,
    /// hydrogens with two (or more) bonds
    removeHigherDegrees: bool,
    /// hydrogens with bonds only to other hydrogens
    removeOnlyHNeighbors: bool,
    /// hydrogens with non-default isotopes
    removeIsotopes: bool,
    /// hydrogens with non-default isotopes and store them in the _isotopicHs atom property such that AddHs() can add the same isotope at a later stage
    removeAndTrackIsotopes: bool,
    /// hydrogens with at least one dummy-atom neighbor
    removeDummyNeighbors: bool,
    /// hydrogens defining bond stereochemistry
    removeDefiningBondStereo: bool,
    /// hydrogens with wedged bonds to them
    removeWithWedgedBond: bool,
    /// hydrogens with queries defined
    removeWithQuery: bool,
    /// mapped hydrogens
    removeMapped: bool,
    /// hydrogens involved in SubstanceGroups
    removeInSGroups: bool,
    /// DEPRECATED
    removeNonimplicit: bool,
    /// hydrogens with formal charge -1
    removeHydrides: bool,
    /// hydrogens with neighbors that have non-tetrahedral stereochemistry
    removeNontetrahedralNeighbors: bool,
    /// display warning messages for some classes of removed Hs
    showWarnings: bool,
    /// DEPRECATED
    updateExplicitCount: bool,
}

impl PyRemoveHsParameters {
    fn from_core(p: &RemoveHsParameters) -> Self {
        Self {
            removeDegreeZero: p.remove_degree_zero,
            removeHigherDegrees: p.remove_higher_degrees,
            removeOnlyHNeighbors: p.remove_only_h_neighbors,
            removeIsotopes: p.remove_isotopes,
            removeAndTrackIsotopes: p.remove_and_track_isotopes,
            removeDummyNeighbors: p.remove_dummy_neighbors,
            removeDefiningBondStereo: p.remove_defining_bond_stereo,
            removeWithWedgedBond: p.remove_with_wedged_bond,
            removeWithQuery: p.remove_with_query,
            removeMapped: p.remove_mapped,
            removeInSGroups: p.remove_in_sgroups,
            removeNonimplicit: p.remove_nonimplicit,
            removeHydrides: p.remove_hydrides,
            removeNontetrahedralNeighbors: p.remove_nontetrahedral_neighbors,
            showWarnings: p.show_warnings,
            updateExplicitCount: p.update_explicit_count,
        }
    }

    fn to_core(&self) -> RemoveHsParameters {
        RemoveHsParameters {
            remove_degree_zero: self.removeDegreeZero,
            remove_higher_degrees: self.removeHigherDegrees,
            remove_only_h_neighbors: self.removeOnlyHNeighbors,
            remove_isotopes: self.removeIsotopes,
            remove_and_track_isotopes: self.removeAndTrackIsotopes,
            remove_dummy_neighbors: self.removeDummyNeighbors,
            remove_defining_bond_stereo: self.removeDefiningBondStereo,
            remove_with_wedged_bond: self.removeWithWedgedBond,
            remove_with_query: self.removeWithQuery,
            remove_mapped: self.removeMapped,
            remove_in_sgroups: self.removeInSGroups,
            show_warnings: self.showWarnings,
            remove_nonimplicit: self.removeNonimplicit,
            update_explicit_count: self.updateExplicitCount,
            remove_hydrides: self.removeHydrides,
            remove_nontetrahedral_neighbors: self.removeNontetrahedralNeighbors,
        }
    }
}

#[pymethods]
impl PyRemoveHsParameters {
    #[new]
    fn py_new() -> Self {
        Self::from_core(&RemoveHsParameters::default())
    }
}

fn run_remove_hs(mol: &ROMol, ps: &RemoveHsParameters, sanitize: bool) -> PyResult<PyMol> {
    let mut res = ROMol::copy_from(mol, false, -1);
    remove_hs(&mut res, ps, sanitize).py()?;
    Ok(PyMol { mol: res })
}

/// Removes any hydrogens from the graph of a molecule.
///
///   ARGUMENTS:
///
///     - mol: the molecule to be modified
///
///     - implicitOnly: (optional) if this toggle is set, only implicit Hs will
///       be removed from the graph.  Default value is 0 (remove implicit and explicit Hs).
///
///     - updateExplicitCount: (optional) if this toggle is set, the explicit H count on atoms with
///       Hs will be updated. Default value is 0 (do not update explicit H count).
///
///     - sanitize: (optional) if this toggle is set, the molecule will be sanitized after the Hs
///       are removed. Default value is 1 (do sanitize).
///
///   RETURNS: a new molecule with the Hs removed
///
///   NOTES:
///
///     - The original molecule is *not* modified.
///     - Hydrogens which aren't connected to a heavy atom will not be
///       removed.  This prevents molecules like [H][H] from having
///       all atoms removed.
///     - Labelled hydrogen (e.g. atoms with atomic number=1, but isotope > 1),
///       will not be removed.
///     - two coordinate Hs, like the central H in C[H-]C, will not be removed
///     - Hs connected to dummy atoms will not be removed
///     - Hs that are part of the definition of double bond Stereochemistry
///       will not be removed
///     - Hs that are not connected to anything else will not be removed
///
/// Alternatively: RemoveHs(mol, params, sanitize=True) with a RemoveHsParameters object.
#[pyfunction]
#[pyo3(name = "RemoveHs", signature = (*args, **kwargs))]
fn py_remove_hs(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMol> {
    let get = |pos: usize, name: &str| -> PyResult<Option<Bound<'_, PyAny>>> {
        let kw = match kwargs {
            Some(k) => k.get_item(name)?,
            None => None,
        };
        let p = if pos < args.len() {
            Some(args.get_item(pos)?)
        } else {
            None
        };
        match (p, kw) {
            (Some(_), Some(_)) => Err(PyTypeError::new_err(format!(
                "RemoveHs() got multiple values for argument '{name}'"
            ))),
            (a, b) => Ok(a.or(b)),
        }
    };
    let mol_obj = get(0, "mol")?
        .ok_or_else(|| PyTypeError::new_err("RemoveHs() missing required argument 'mol'"))?;
    let mol = mol_obj.cast::<PyMol>()?.borrow();
    let second = match kwargs.map(|k| k.get_item("params")).transpose()?.flatten() {
        Some(p) => Some(p),
        None => {
            if args.len() > 1 {
                Some(args.get_item(1)?)
            } else {
                None
            }
        }
    };
    if let Some(p) = second
        .as_ref()
        .and_then(|s| s.cast::<PyRemoveHsParameters>().ok())
    {
        let ps = p.borrow().to_core();
        let sanitize = match get(2, "sanitize")? {
            Some(s) => s.is_truthy()?,
            None => true,
        };
        return run_remove_hs(&mol.mol, &ps, sanitize);
    }
    let flag = |pos, name, default| -> PyResult<bool> {
        Ok(match get(pos, name)? {
            Some(v) => v.is_truthy()?,
            None => default,
        })
    };
    let implicit_only = flag(1, "implicitOnly", false)?;
    let update_explicit_count = flag(2, "updateExplicitCount", false)?;
    let sanitize = flag(3, "sanitize", true)?;
    let ps = RemoveHsParameters {
        remove_nonimplicit: !implicit_only,
        update_explicit_count,
        ..Default::default()
    };
    run_remove_hs(&mol.mol, &ps, sanitize)
}

/// Returns a copy of the molecule with all Hs removed.
#[pyfunction]
#[pyo3(name = "RemoveAllHs", signature = (mol, sanitize = true))]
fn py_remove_all_hs(mol: PyRef<'_, PyMol>, sanitize: bool) -> PyResult<PyMol> {
    let ps = RemoveHsParameters {
        remove_degree_zero: true,
        remove_higher_degrees: true,
        remove_only_h_neighbors: true,
        remove_isotopes: true,
        remove_dummy_neighbors: true,
        remove_defining_bond_stereo: true,
        remove_with_wedged_bond: true,
        remove_with_query: true,
        remove_nonimplicit: true,
        remove_in_sgroups: true,
        show_warnings: false,
        remove_hydrides: true,
        remove_nontetrahedral_neighbors: true,
        ..Default::default()
    };
    run_remove_hs(&mol.mol, &ps, sanitize)
}

/// Does the CIP stereochemistry assignment
///   for the molecule's atoms (R/S) and double bond (Z/E).
///   Chiral atoms will have a property '_CIPCode' indicating
///   their chiral code.
///
///   ARGUMENTS:
///
///     - mol: the molecule to use
///     - cleanIt: (optional) if provided, any existing values of the property `_CIPCode`
///         will be cleared, atoms with a chiral specifier that aren't
///       actually chiral (e.g. atoms with duplicate substituents or only 2 substituents,
///       etc.) will have their chiral code set to CHI_UNSPECIFIED. Bonds with
///       STEREOCIS/STEREOTRANS specified that have duplicate substituents based upon the CIP
///       atom ranks will be marked STEREONONE.
///     - force: (optional) causes the calculation to be repeated, even if it has already
///       been done
///     - flagPossibleStereoCenters (optional)   set the _ChiralityPossible property on
///       atoms that are possible stereocenters
#[pyfunction]
#[pyo3(name = "AssignStereochemistry", signature = (mol, cleanIt = false, force = false, flagPossibleStereoCenters = false))]
#[allow(non_snake_case)]
fn py_assign_stereochemistry(
    mut mol: PyRefMut<'_, PyMol>,
    cleanIt: bool,
    force: bool,
    flagPossibleStereoCenters: bool,
) -> PyResult<()> {
    legacy::assign_stereochemistry(&mut mol.mol, cleanIt, force, flagPossibleStereoCenters).py()
}

/// Computes the CIP ranks for the atoms in a molecule.
///   The ranks are stored as an atom property '_CIPRank' and returned as a tuple.
#[pyfunction]
#[pyo3(name = "ComputeAtomCIPRanks")]
fn py_compute_atom_cip_ranks(
    py: Python<'_>,
    mut mol: PyRefMut<'_, PyMol>,
) -> PyResult<Py<PyTuple>> {
    let mut ranks = Vec::new();
    legacy::assign_atom_cip_ranks(&mut mol.mol, &mut ranks).py()?;
    Ok(PyTuple::new(py, ranks)?.unbind())
}

/// Uses the stereo info on double bonds to set the directions of neighboring single bonds
///
///   ARGUMENTS:
///
///     - mol: the molecule to be modified
#[pyfunction]
#[pyo3(name = "SetDoubleBondNeighborDirections", signature = (mol, conf = None))]
fn py_set_double_bond_neighbor_directions(
    mut mol: PyRefMut<'_, PyMol>,
    conf: Option<&Bound<'_, PyAny>>,
) -> PyResult<()> {
    let conf = match conf {
        Some(c) if !c.is_none() => {
            let id: i64 = c.call_method0("GetId")?.extract()?;
            mol.mol
                .conformers()
                .iter()
                .position(|cf| i64::from(cf.id()) == id)
        }
        _ => None,
    };
    bond_dirs::set_double_bond_neighbor_directions(&mut mol.mol, conf).py()
}

/// Parameters controlling which Hs are added.
#[pyclass(
    name = "AddHsParameters",
    module = "rdkit.Chem.rdmolops",
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
pub struct PyAddHsParameters {
    /// only add explicit Hs
    #[pyo3(get, set, name = "explicitOnly")]
    explicit_only: bool,
    /// add coordinates to the Hs
    #[pyo3(get, set, name = "addCoords")]
    add_coords: bool,
    /// add residue info to the Hs
    #[pyo3(get, set, name = "addResidueInfo")]
    add_residue_info: bool,
    /// do not add Hs to query atoms or atoms with query bonds
    #[pyo3(get, set, name = "skipQueries")]
    skip_queries: bool,
}

#[pymethods]
impl PyAddHsParameters {
    #[new]
    fn new() -> Self {
        Self::default()
    }
}

/// Adds hydrogens to the graph of a molecule.
///
///   AddHs(mol, explicitOnly=False, addCoords=False, onlyOnAtoms=None, addResidueInfo=False)
///   or AddHs(mol, params, onlyOnAtoms=None) with an AddHsParameters object.
#[pyfunction]
#[pyo3(name = "AddHs", signature = (*args, **kwargs))]
fn py_add_hs(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMol> {
    let get = |pos: usize, name: &str| -> PyResult<Option<Bound<'_, PyAny>>> {
        let kw = match kwargs {
            Some(k) => k.get_item(name)?,
            None => None,
        };
        let p = if pos < args.len() {
            Some(args.get_item(pos)?)
        } else {
            None
        };
        match (p, kw) {
            (Some(_), Some(_)) => Err(PyTypeError::new_err(format!(
                "AddHs() got multiple values for argument '{name}'"
            ))),
            (a, b) => Ok(a.or(b)),
        }
    };
    let mol_obj = get(0, "mol")?
        .ok_or_else(|| PyTypeError::new_err("AddHs() missing required argument 'mol'"))?;
    let mol = mol_obj.cast::<PyMol>()?.borrow();
    let second = match kwargs.map(|k| k.get_item("params")).transpose()?.flatten() {
        Some(p) => Some(p),
        None => {
            if args.len() > 1 {
                Some(args.get_item(1)?)
            } else {
                None
            }
        }
    };
    let (params, only_on) = if let Some(p) = second
        .as_ref()
        .and_then(|s| s.cast::<PyAddHsParameters>().ok())
    {
        let p = p.borrow();
        (
            AddHsParameters {
                explicit_only: p.explicit_only,
                add_coords: p.add_coords,
                add_residue_info: p.add_residue_info,
                skip_queries: p.skip_queries,
            },
            get(2, "onlyOnAtoms")?,
        )
    } else {
        let flag = |pos, name| -> PyResult<bool> {
            Ok(match get(pos, name)? {
                Some(v) => v.is_truthy()?,
                None => false,
            })
        };
        (
            AddHsParameters {
                explicit_only: flag(1, "explicitOnly")?,
                add_coords: flag(2, "addCoords")?,
                add_residue_info: flag(4, "addResidueInfo")?,
                skip_queries: false,
            },
            get(3, "onlyOnAtoms")?,
        )
    };
    // pythonObjectToVect(onlyOnAtoms, nAtoms): falsy -> all atoms
    let only = match only_on {
        Some(o) if o.is_truthy()? => {
            let n = mol.mol.num_atoms();
            let mut v = Vec::new();
            for x in o.try_iter()? {
                let x: u32 = x?.extract()?;
                if x >= n {
                    return Err(pyo3::exceptions::PyValueError::new_err(
                        "list element larger than allowed value",
                    ));
                }
                v.push(x);
            }
            Some(v)
        }
        _ => None,
    };
    let mut res: ROMol = mol.mol.clone();
    add_hs(&mut res, &params, only.as_deref()).py()?;
    Ok(PyMol { mol: res })
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRemoveHsParameters>()?;
    m.add_class::<PyAddHsParameters>()?;
    m.add_function(wrap_pyfunction!(py_add_hs, m)?)?;
    m.add_function(wrap_pyfunction!(py_remove_hs, m)?)?;
    m.add_function(wrap_pyfunction!(py_remove_all_hs, m)?)?;
    m.add_function(wrap_pyfunction!(py_assign_stereochemistry, m)?)?;
    m.add_function(wrap_pyfunction!(py_compute_atom_cip_ranks, m)?)?;
    m.add_function(wrap_pyfunction!(py_set_double_bond_neighbor_directions, m)?)?;
    Ok(())
}
