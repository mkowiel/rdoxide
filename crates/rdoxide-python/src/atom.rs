//! `rdchem.Atom` — port of upstream `Code/GraphMol/Wrap/Atom.cpp`.

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::graphmol::atom as core_atom;
use rdoxide::graphmol::queries::{AtomQuery, CompositeQueryType};
use rdoxide::graphmol::query_ops::atom_expand_query;
use rdoxide::graphmol::{Atom, ROMol, ValenceType};
use rdoxide::rdgeneral::Dict;

use crate::bond::PyBond;
use crate::enums_impl::{PyChiralType, PyHybridizationType, PyValenceType};
use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::monomer_info::{self, PyMonomerInfo};
use crate::props::impl_props;

/// Where the atom lives.
#[derive(Debug)]
pub enum AtomHandle {
    /// A free-standing atom (`Chem.Atom(6)`).
    Owned(Atom),
    /// An atom inside a molecule.
    InMol { mol: Py<PyMol>, idx: u32 },
}

/// The class to store Atoms.
/// Note that, though it is possible to create one, having an Atom on its own
/// (i.e not associated with a molecule) is not particularly useful.
#[pyclass(name = "Atom", module = "rdoxide.Chem.rdchem", subclass)]
#[derive(Debug)]
pub struct PyAtom {
    pub handle: AtomHandle,
}

fn gone() -> PyErr {
    PyRuntimeError::new_err("atom is no longer part of its molecule")
}

fn no_owner() -> PyErr {
    crate::errors::to_py(rdoxide::RDKitError::precondition_at(
        "no owner",
        "dp_mol",
        "Code/GraphMol/Atom.h",
        158,
    ))
}

/// The class to store QueryAtoms.
/// These cannot currently be constructed directly from Python
#[pyclass(name = "QueryAtom", module = "rdoxide.Chem.rdchem", extends = PyAtom)]
pub struct PyQueryAtom;

impl PyAtom {
    /// The Python object for an atom: a `QueryAtom` when the atom carries a
    /// query (upstream returns the most derived class), else an `Atom`.
    pub fn from_handle(py: Python<'_>, handle: AtomHandle) -> PyResult<Py<PyAtom>> {
        let is_query = match &handle {
            AtomHandle::Owned(a) => a.has_query(),
            // the molecule may be borrowed mutably by the caller; the plain
            // class is the safe answer then
            AtomHandle::InMol { mol, idx } => mol.bind(py).try_borrow().is_ok_and(|m| {
                m.mol
                    .atoms()
                    .get(*idx as usize)
                    .is_some_and(Atom::has_query)
            }),
        };
        let init = PyClassInitializer::from(PyAtom { handle });
        if is_query {
            let obj = Py::new(py, init.add_subclass(PyQueryAtom))?;
            Ok(obj
                .into_bound(py)
                .into_any()
                .cast_into::<PyAtom>()?
                .unbind())
        } else {
            Py::new(py, init)
        }
    }

    pub fn in_mol(py: Python<'_>, mol: Py<PyMol>, idx: u32) -> PyResult<Py<PyAtom>> {
        Self::from_handle(py, AtomHandle::InMol { mol, idx })
    }

    /// A free-standing atom.
    pub fn owned(py: Python<'_>, atom: Atom) -> PyResult<Py<PyAtom>> {
        Self::from_handle(py, AtomHandle::Owned(atom))
    }

    /// Run `f` on the owning molecule and this atom's index (error for
    /// free-standing atoms).
    pub fn with_owner<R>(
        &self,
        f: impl FnOnce(&rdoxide::graphmol::romol::ROMol, u32) -> PyResult<R>,
    ) -> PyResult<R> {
        match &self.handle {
            AtomHandle::Owned(_) => Err(no_owner()),
            AtomHandle::InMol { mol, idx } => Python::attach(|py| {
                let m = mol.bind(py).borrow();
                if *idx >= m.mol.num_atoms() {
                    return Err(gone());
                }
                f(&m.mol, *idx)
            }),
        }
    }

    /// Run `f` on the atom.
    pub fn with<R>(&self, py: Python<'_>, f: impl FnOnce(&Atom) -> R) -> PyResult<R> {
        match &self.handle {
            AtomHandle::Owned(a) => Ok(f(a)),
            AtomHandle::InMol { mol, idx } => {
                let m = mol.bind(py).borrow();
                m.mol.atoms().get(*idx as usize).map(f).ok_or_else(gone)
            }
        }
    }

    /// Run `f` on the (mutable) atom.
    pub fn with_mut<R>(&mut self, py: Python<'_>, f: impl FnOnce(&mut Atom) -> R) -> PyResult<R> {
        match &mut self.handle {
            AtomHandle::Owned(a) => Ok(f(a)),
            AtomHandle::InMol { mol, idx } => {
                let mut m = mol.bind(py).borrow_mut();
                if *idx >= m.mol.num_atoms() {
                    return Err(gone());
                }
                Ok(f(m.mol.atom_mut(*idx)))
            }
        }
    }

    /// Run `f` with the owning molecule and atom index (`PRECONDITION(dp_mol)`).
    pub fn with_mol<R>(&self, py: Python<'_>, f: impl FnOnce(&ROMol, u32) -> R) -> PyResult<R> {
        match &self.handle {
            AtomHandle::Owned(_) => Err(no_owner()),
            AtomHandle::InMol { mol, idx } => {
                let m = mol.bind(py).borrow();
                if *idx >= m.mol.num_atoms() {
                    return Err(gone());
                }
                Ok(f(&m.mol, *idx))
            }
        }
    }

    fn with_mol_mut<R>(&self, py: Python<'_>, f: impl FnOnce(&mut ROMol, u32) -> R) -> PyResult<R> {
        match &self.handle {
            AtomHandle::Owned(_) => Err(no_owner()),
            AtomHandle::InMol { mol, idx } => {
                let mut m = mol.bind(py).borrow_mut();
                if *idx >= m.mol.num_atoms() {
                    return Err(gone());
                }
                Ok(f(&mut m.mol, *idx))
            }
        }
    }

    /// A copy of the atom's data.
    pub fn snapshot(&self, py: Python<'_>) -> PyResult<Atom> {
        self.with(py, Atom::clone)
    }

    pub fn with_dict<R>(&self, py: Python<'_>, f: impl FnOnce(&Dict) -> R) -> PyResult<R> {
        self.with(py, |a| f(rdoxide::rdgeneral::RDProps::dict(a)))
    }
    pub fn with_dict_mut<R>(
        &mut self,
        py: Python<'_>,
        f: impl FnOnce(&mut Dict) -> R,
    ) -> PyResult<R> {
        self.with_mut(py, |a| f(rdoxide::rdgeneral::RDProps::dict_mut(a)))
    }

    fn deprecation_warning(message: &str) {
        rdoxide::rd_log_warning!("DEPRECATION WARNING: {message}");
    }
}

impl_props!(PyAtom, computed = false, dict_defaults = (true, true));

#[derive(FromPyObject)]
enum AtomInit<'py> {
    Num(u32),
    Symbol(String),
    Atom(PyRef<'py, PyAtom>),
}

#[pymethods]
impl PyAtom {
    /// Sets an atomic property
    ///
    ///   ARGUMENTS:
    ///     - key: the name of the property to be set (an ExplicitBitVect).
    ///     - value: the property value (an ExplicitBitVect).
    #[pyo3(name = "SetExplicitBitVectProp", signature = (key, val))]
    fn set_explicit_bit_vect_prop(
        &mut self,
        py: Python<'_>,
        key: &str,
        val: PyRef<'_, crate::datastructs::PyExplicitBitVect>,
    ) -> PyResult<()> {
        let bv = val.bv.clone();
        self.with_dict_mut(py, |d| crate::props::set_prop(d, key, bv.into(), false))?
    }

    /// Returns the value of the property.
    ///
    ///   ARGUMENTS:
    ///     - key: the name of the property to return (a ExplicitBitVect).
    ///
    ///   RETURNS: an ExplicitBitVect
    ///
    ///   NOTE:
    ///     - If the property has not been set, a KeyError exception will be raised.
    #[pyo3(name = "GetExplicitBitVectProp", signature = (key))]
    fn get_explicit_bit_vect_prop(
        &self,
        py: Python<'_>,
        key: &str,
    ) -> PyResult<Py<crate::datastructs::PyExplicitBitVect>> {
        let found = self.with_dict(py, |d| d.get(key).cloned())?;
        match found {
            None => Err(pyo3::exceptions::PyKeyError::new_err(key.to_string())),
            Some(rdoxide::rdgeneral::RDValue::BitVect(bv)) => {
                Py::new(py, crate::datastructs::PyExplicitBitVect::wrap(bv))
            }
            Some(_) => Err(PyValueError::new_err(format!(
                "key `{key}` exists but does not result in unregistered C++ type reason: bad any cast"
            ))),
        }
    }

    #[new]
    #[pyo3(signature = (what))]
    fn py_new(py: Python<'_>, what: AtomInit<'_>) -> PyResult<Self> {
        let atom = match what {
            AtomInit::Num(n) => Atom::new(n),
            AtomInit::Symbol(s) => Atom::from_symbol(&s).py()?,
            AtomInit::Atom(a) => a.snapshot(py)?.copy_detached(),
        };
        Ok(PyAtom {
            handle: AtomHandle::Owned(atom),
        })
    }

    /// marker for unspecified int values
    #[classattr]
    #[allow(non_snake_case)]
    fn NOATOM() -> u32 {
        u32::MAX
    }

    /// Create a copy of the atom
    fn __copy__(&self, py: Python<'_>) -> PyResult<Py<PyAtom>> {
        Self::owned(py, self.snapshot(py)?.copy_detached())
    }

    /// Returns the atomic number.
    #[pyo3(name = "GetAtomicNum")]
    fn get_atomic_num(&self, py: Python<'_>) -> PyResult<u32> {
        self.with(py, Atom::atomic_num)
    }

    /// Sets the atomic number, takes an integer value as an argument
    #[pyo3(name = "SetAtomicNum")]
    #[allow(non_snake_case)]
    fn set_atomic_num(&mut self, py: Python<'_>, newNum: i32) -> PyResult<()> {
        self.with_mut(py, |a| a.set_atomic_num(newNum))
    }

    /// Returns the atomic symbol (a string)
    #[pyo3(name = "GetSymbol")]
    fn get_symbol(&self, py: Python<'_>) -> PyResult<String> {
        self.with(py, Atom::symbol)
    }

    /// Returns the atom's index (ordering in the molecule)
    #[pyo3(name = "GetIdx")]
    fn get_idx(&self, py: Python<'_>) -> PyResult<u32> {
        self.with(py, Atom::idx)
    }

    /// Returns the degree of the atom in the molecule.
    ///
    ///   The degree of an atom is defined to be its number of
    ///   directly-bonded neighbors.
    ///   The degree is independent of bond orders, but is dependent
    ///     on whether or not Hs are explicit in the graph.
    #[pyo3(name = "GetDegree")]
    fn get_degree(&self, py: Python<'_>) -> PyResult<u32> {
        match &self.handle {
            AtomHandle::Owned(_) => Ok(0),
            _ => self.with_mol(py, |m, i| m.atom_degree(i)),
        }
    }

    /// Returns the degree of the atom in the molecule including Hs.
    #[pyo3(name = "GetTotalDegree")]
    fn get_total_degree(&self, py: Python<'_>) -> PyResult<u32> {
        match &self.handle {
            AtomHandle::Owned(a) => a.num_implicit_hs().map(|h| h + a.num_explicit_hs()).py(),
            _ => self.with_mol(py, |m, i| m.total_degree(i))?.py(),
        }
    }

    /// Returns the total number of Hs (explicit and implicit) on the atom.
    ///
    ///   ARGUMENTS:
    ///
    ///     - includeNeighbors: (optional) toggles inclusion of neighboring H atoms in the sum.
    ///       Defaults to 0.
    #[pyo3(name = "GetTotalNumHs", signature = (includeNeighbors = false))]
    #[allow(non_snake_case)]
    fn get_total_num_hs(&self, py: Python<'_>, includeNeighbors: bool) -> PyResult<u32> {
        match &self.handle {
            AtomHandle::Owned(a) => a.num_implicit_hs().map(|h| h + a.num_explicit_hs()).py(),
            _ => self
                .with_mol(py, |m, i| m.total_num_hs(i, includeNeighbors))?
                .py(),
        }
    }

    /// Returns the total number of implicit Hs on the atom.
    #[pyo3(name = "GetNumImplicitHs")]
    fn get_num_implicit_hs(&self, py: Python<'_>) -> PyResult<u32> {
        self.with(py, Atom::num_implicit_hs)?.py()
    }

    /// DEPRECATED, please use GetValence(Chem.ValenceType.EXPLICIT) instead.
    /// Returns the explicit valence of the atom.
    #[pyo3(name = "GetExplicitValence")]
    fn get_explicit_valence(&self, py: Python<'_>) -> PyResult<u32> {
        Self::deprecation_warning("please use GetValence(Chem.ValenceType.EXPLICIT) instead");
        self.with(py, |a| a.valence(ValenceType::Explicit))?.py()
    }

    /// DEPRECATED, please use GetValence(Chem.ValenceType.IMPLICIT) instead.
    /// Returns the number of implicit Hs on the atom.
    #[pyo3(name = "GetImplicitValence")]
    fn get_implicit_valence(&self, py: Python<'_>) -> PyResult<u32> {
        Self::deprecation_warning("please use GetValence(Chem.ValenceType.IMPLICIT) instead");
        self.with(py, |a| a.valence(ValenceType::Implicit))?.py()
    }

    /// Returns the valence (Chem.ValenceType.EXPLICIT or Chem.ValenceType.IMPLICIT) of the atom.
    #[pyo3(name = "GetValence")]
    fn get_valence(&self, py: Python<'_>, which: &Bound<'_, PyAny>) -> PyResult<u32> {
        let which = PyValenceType::extract(which)?;
        self.with(py, |a| a.valence(which))?.py()
    }

    /// Returns the total valence (explicit + implicit) of the atom.
    #[pyo3(name = "GetTotalValence")]
    fn get_total_valence(&self, py: Python<'_>) -> PyResult<u32> {
        self.with(py, Atom::total_valence)?.py()
    }

    /// Returns whether the atom has a valence violation or not.
    #[pyo3(name = "HasValenceViolation")]
    fn has_valence_violation(&self, py: Python<'_>) -> PyResult<bool> {
        self.with_mol(py, core_atom::has_valence_violation)
    }

    #[pyo3(name = "GetFormalCharge")]
    fn get_formal_charge(&self, py: Python<'_>) -> PyResult<i32> {
        self.with(py, Atom::formal_charge)
    }

    #[pyo3(name = "SetFormalCharge")]
    fn set_formal_charge(&mut self, py: Python<'_>, what: i32) -> PyResult<()> {
        self.with_mut(py, |a| a.set_formal_charge(what))
    }

    /// Sets a marker on the atom that *disallows* implicit Hs.
    ///   This holds even if the atom would otherwise have implicit Hs added.
    #[pyo3(name = "SetNoImplicit")]
    fn set_no_implicit(&mut self, py: Python<'_>, what: bool) -> PyResult<()> {
        self.with_mut(py, |a| a.set_no_implicit(what))
    }

    /// Returns whether or not the atom is *allowed* to have implicit Hs.
    #[pyo3(name = "GetNoImplicit")]
    fn get_no_implicit(&self, py: Python<'_>) -> PyResult<bool> {
        self.with(py, Atom::no_implicit)
    }

    #[pyo3(name = "SetNumExplicitHs")]
    fn set_num_explicit_hs(&mut self, py: Python<'_>, what: u32) -> PyResult<()> {
        self.with_mut(py, |a| a.set_num_explicit_hs(what))
    }

    #[pyo3(name = "GetNumExplicitHs")]
    fn get_num_explicit_hs(&self, py: Python<'_>) -> PyResult<u32> {
        self.with(py, Atom::num_explicit_hs)
    }

    #[pyo3(name = "SetIsAromatic")]
    fn set_is_aromatic(&mut self, py: Python<'_>, what: bool) -> PyResult<()> {
        self.with_mut(py, |a| a.set_is_aromatic(what))
    }

    #[pyo3(name = "GetIsAromatic")]
    fn get_is_aromatic(&self, py: Python<'_>) -> PyResult<bool> {
        self.with(py, Atom::is_aromatic)
    }

    #[pyo3(name = "GetMass")]
    fn get_mass(&self, py: Python<'_>) -> PyResult<f64> {
        self.with(py, Atom::mass)
    }

    #[pyo3(name = "SetIsotope")]
    fn set_isotope(&mut self, py: Python<'_>, what: u32) -> PyResult<()> {
        self.with_mut(py, |a| a.set_isotope(what))
    }

    #[pyo3(name = "GetIsotope")]
    fn get_isotope(&self, py: Python<'_>) -> PyResult<u32> {
        self.with(py, Atom::isotope)
    }

    #[pyo3(name = "SetNumRadicalElectrons")]
    fn set_num_radical_electrons(&mut self, py: Python<'_>, num: u32) -> PyResult<()> {
        self.with_mut(py, |a| a.set_num_radical_electrons(num))
    }

    #[pyo3(name = "GetNumRadicalElectrons")]
    fn get_num_radical_electrons(&self, py: Python<'_>) -> PyResult<u32> {
        self.with(py, Atom::num_radical_electrons)
    }

    #[pyo3(name = "GetQueryType")]
    fn get_query_type(&self, py: Python<'_>) -> PyResult<String> {
        self.with(py, |a| a.query_type().to_string())
    }

    /// returns the SMARTS (or SMILES) string for an Atom
    #[pyo3(name = "GetSmarts", signature = (doKekule = false, allHsExplicit = false, isomericSmiles = true))]
    #[allow(non_snake_case)]
    fn get_smarts(
        &mut self,
        py: Python<'_>,
        doKekule: bool,
        allHsExplicit: bool,
        isomericSmiles: bool,
    ) -> PyResult<String> {
        use rdoxide::graphmol::smiles_parse::smarts_write::{get_atom_smarts, get_atom_smarts_for};
        use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, get_atom_smiles};
        let smiles_ps = SmilesWriteParams {
            do_isomeric_smiles: isomericSmiles,
            do_kekule: doKekule,
            all_hs_explicit: allHsExplicit,
            ..SmilesWriteParams::default()
        };
        let defaults = SmilesWriteParams::default();
        let has_query = self.with(py, Atom::has_query)?;
        match &self.handle {
            AtomHandle::InMol { .. } => self
                .with_mol_mut(py, |m, i| {
                    if has_query {
                        get_atom_smarts(m, i, &defaults)
                    } else {
                        get_atom_smiles(m, i, &smiles_ps)
                    }
                })?
                .py(),
            AtomHandle::Owned(_) => {
                if has_query {
                    self.with_mut(py, |a| get_atom_smarts_for(a, false, &defaults))?
                        .py()
                } else {
                    // upstream works on the free-standing atom; without a
                    // property cache this raises like upstream does
                    let a = self.with(py, Atom::copy_detached)?;
                    let mut tmp = ROMol::new();
                    tmp.add_atom(a, false);
                    get_atom_smiles(&tmp, 0, &smiles_ps).py()
                }
            }
        }
    }

    /// returns a text description of the query. Primarily intended for debugging purposes.
    #[pyo3(name = "DescribeQuery")]
    fn describe_query(&self, py: Python<'_>) -> PyResult<String> {
        self.with(py, |a| a.query().map_or(String::new(), |q| q.describe(0)))
    }

    #[pyo3(name = "SetChiralTag")]
    fn set_chiral_tag(&mut self, py: Python<'_>, what: &Bound<'_, PyAny>) -> PyResult<()> {
        let t = PyChiralType::extract(what)?;
        self.with_mut(py, |a| a.set_chiral_tag(t))
    }

    #[pyo3(name = "InvertChirality")]
    fn invert_chirality(&mut self, py: Python<'_>) -> PyResult<bool> {
        self.with_mut(py, Atom::invert_chirality)
    }

    #[pyo3(name = "GetChiralTag")]
    fn get_chiral_tag(&self, py: Python<'_>) -> PyResult<PyChiralType> {
        self.with(py, |a| PyChiralType(a.chiral_tag()))
    }

    #[pyo3(name = "SetHybridization")]
    fn set_hybridization(&mut self, py: Python<'_>, what: &Bound<'_, PyAny>) -> PyResult<()> {
        let h = PyHybridizationType::extract(what)?;
        self.with_mut(py, |a| a.set_hybridization(h))
    }

    #[pyo3(name = "GetHybridization")]
    fn get_hybridization(&self, py: Python<'_>) -> PyResult<PyHybridizationType> {
        self.with(py, |a| PyHybridizationType(a.hybridization()))
    }

    /// Returns whether or not this instance belongs to a molecule.
    #[pyo3(name = "HasOwningMol")]
    fn has_owning_mol(&self) -> bool {
        matches!(self.handle, AtomHandle::InMol { .. })
    }

    /// Returns the Mol that owns this atom.
    #[pyo3(name = "GetOwningMol")]
    fn get_owning_mol(&self, py: Python<'_>) -> PyResult<Py<PyMol>> {
        match &self.handle {
            AtomHandle::InMol { mol, .. } => Ok(mol.clone_ref(py)),
            AtomHandle::Owned(_) => Err(no_owner()),
        }
    }

    /// Returns a read-only sequence of the atom's neighbors
    #[pyo3(name = "GetNeighbors")]
    fn get_neighbors<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let AtomHandle::InMol { mol, .. } = &self.handle else {
            return Err(no_owner());
        };
        let nbrs: Vec<u32> = self.with_mol(py, |m, i| m.atom_neighbors(i).collect())?;
        let items = nbrs
            .into_iter()
            .map(|n| PyAtom::in_mol(py, mol.clone_ref(py), n))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, items)
    }

    /// Returns a read-only sequence of the atom's bonds
    #[pyo3(name = "GetBonds")]
    fn get_bonds<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let AtomHandle::InMol { mol, .. } = &self.handle else {
            return Err(no_owner());
        };
        let bonds: Vec<u32> = self.with_mol(py, |m, i| m.atom_bond_indices(i).collect())?;
        let items = bonds
            .into_iter()
            .map(|b| PyBond::in_mol(py, mol.clone_ref(py), b))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, items)
    }

    /// Returns whether or not this atom matches another Atom.
    ///
    ///   Each Atom (or query Atom) has a query function which is
    ///   used for this type of matching.
    ///
    ///   ARGUMENTS:
    ///     - other: the other Atom to which to compare
    #[pyo3(name = "Match")]
    fn py_match(&self, py: Python<'_>, what: PyRef<'_, PyAtom>) -> PyResult<bool> {
        let Some(query) = self.with(py, |a| a.query().cloned())? else {
            let other = what.snapshot(py)?;
            return self.with(py, |a| a.matches(&other));
        };
        // a query atom evaluates its query on the other atom in its context
        match &what.handle {
            AtomHandle::InMol { .. } => what.with_owner(|mol, idx| query.matches(mol, idx).py()),
            AtomHandle::Owned(atom) => {
                let mut single = ROMol::new();
                single.add_atom(atom.copy_detached(), false);
                query.matches(&single, 0).py()
            }
        }
    }

    /// Returns whether or not the atom is in a ring of a particular size.
    ///
    ///   ARGUMENTS:
    ///     - size: the ring size to look for
    #[pyo3(name = "IsInRingSize")]
    fn is_in_ring_size(&self, py: Python<'_>, size: i32) -> PyResult<bool> {
        self.ensure_rings(py)?;
        self.with_mol(py, |m, i| {
            m.ring_info().is_atom_in_ring_of_size(i, size as u32)
        })?
        .py()
    }

    /// Returns whether or not the atom is in a ring
    #[pyo3(name = "IsInRing")]
    fn is_in_ring(&self, py: Python<'_>) -> PyResult<bool> {
        self.ensure_rings(py)?;
        Ok(self
            .with_mol(py, |m, i| m.ring_info().num_atom_rings(i))?
            .py()?
            != 0)
    }

    /// Returns whether or not the atom has an associated query
    #[pyo3(name = "HasQuery")]
    fn has_query(&self, py: Python<'_>) -> PyResult<bool> {
        self.with(py, Atom::has_query)
    }

    /// Updates cached properties (implicit valence, etc.) on the atom.
    #[pyo3(name = "UpdatePropertyCache", signature = (strict = true))]
    fn update_property_cache(&mut self, py: Python<'_>, strict: bool) -> PyResult<()> {
        self.with_mol_mut(py, |m, i| m.update_atom_property_cache(i, strict))?
            .py()
    }

    /// Returns true or false depending on whether implicit and explicit valence of the atom have already been calculated.
    #[pyo3(name = "NeedsUpdatePropertyCache")]
    fn needs_update_property_cache(&self, py: Python<'_>) -> PyResult<bool> {
        self.with(py, Atom::needs_update_property_cache)
    }

    /// Clears implicit and explicit valence information from the atom.
    #[pyo3(name = "ClearPropertyCache")]
    fn clear_property_cache(&mut self, py: Python<'_>) -> PyResult<()> {
        self.with_mut(py, Atom::clear_property_cache)
    }

    /// Returns the atom's MonomerInfo object, if there is one.
    #[pyo3(name = "GetMonomerInfo")]
    fn get_monomer_info(slf: &Bound<'_, Self>) -> PyResult<Option<Py<PyMonomerInfo>>> {
        monomer_info::get_monomer_info(slf.py(), &slf.clone().unbind())
    }

    /// Returns the atom's MonomerInfo object, if there is one.
    #[pyo3(name = "GetPDBResidueInfo")]
    fn get_pdb_residue_info(slf: &Bound<'_, Self>) -> PyResult<Option<Py<PyMonomerInfo>>> {
        monomer_info::get_pdb_residue_info(slf.py(), &slf.clone().unbind())
    }

    /// Sets the atom's MonomerInfo object.
    #[pyo3(name = "SetMonomerInfo")]
    fn set_monomer_info(
        &mut self,
        py: Python<'_>,
        info: Option<&Bound<'_, PyMonomerInfo>>,
    ) -> PyResult<()> {
        monomer_info::set_monomer_info(py, self, info, false)
    }

    /// Sets the atom's MonomerInfo object.
    #[pyo3(name = "SetPDBResidueInfo")]
    fn set_pdb_residue_info(
        &mut self,
        py: Python<'_>,
        info: Option<&Bound<'_, PyMonomerInfo>>,
    ) -> PyResult<()> {
        monomer_info::set_monomer_info(py, self, info, true)
    }

    /// Gets the atoms map number, returns 0 if not set
    #[pyo3(name = "GetAtomMapNum")]
    fn get_atom_map_num(&self, py: Python<'_>) -> PyResult<i32> {
        self.with(py, Atom::atom_map_num)
    }

    /// Sets the atoms map number, a value of 0 clears the atom map
    #[pyo3(name = "SetAtomMapNum", signature = (mapno, strict = false))]
    fn set_atom_map_num(&mut self, py: Python<'_>, mapno: i32, strict: bool) -> PyResult<()> {
        self.with_mut(py, |a| a.set_atom_map_num(mapno, strict))?
            .py()
    }

    /// Sets the query on this atom (plain atoms have no query).
    #[pyo3(name = "SetQuery")]
    fn set_query(&mut self, other: Option<PyRef<'_, PyAtom>>) -> PyResult<()> {
        if other.is_none() {
            return Err(PyValueError::new_err("other Atom is null"));
        }
        Err(crate::errors::to_py(rdoxide::RDKitError::precondition(
            "plain atoms have no Query",
        )))
    }
}

impl PyAtom {
    /// Upstream calls `MolOps::findSSSR` when ring info is missing; that
    /// arrives with stage 3 (ring perception).
    fn ensure_rings(&self, py: Python<'_>) -> PyResult<()> {
        let ok = self.with_mol(py, |m, _| m.ring_info().is_sssr_or_better())?;
        if ok {
            return Ok(());
        }
        self.with_mol_mut(py, |m, _| crate::perception::find_sssr_if_needed(m))?
            .py()
    }
}

pub(crate) fn expand_how(how: Option<&Bound<'_, PyAny>>) -> PyResult<CompositeQueryType> {
    let how = how.map_or(Ok(0), |h| h.extract::<i64>())?;
    match how {
        0 => Ok(CompositeQueryType::And),
        1 => Ok(CompositeQueryType::Or),
        2 => Ok(CompositeQueryType::XOr),
        _ => Err(PyValueError::new_err("bad CompositeQueryType")),
    }
}

/// The query of `other`, if it has one; `other` must be a `QueryAtom`.
fn other_query(
    slf: &Bound<'_, PyQueryAtom>,
    other: &Bound<'_, PyAny>,
    method: &str,
    extra: &[Bound<'_, PyAny>],
) -> PyResult<Option<AtomQuery>> {
    let py = slf.py();
    if other.is_none() {
        return Err(PyValueError::new_err("other Atom is null"));
    }
    let Ok(other) = other.cast::<PyQueryAtom>() else {
        let mut got = vec![
            "QueryAtom".to_string(),
            crate::errors::boost_type_name(other),
        ];
        got.extend(extra.iter().map(crate::errors::boost_type_name));
        let expected = if method == "ExpandQuery" {
            "ExpandQuery(RDKit::QueryAtom* self, RDKit::QueryAtom const* other, Queries::CompositeQueryType how=rdoxide.Chem.rdchem.CompositeQueryType.COMPOSITE_AND, bool maintainOrder=True)"
        } else {
            "SetQuery(RDKit::QueryAtom* self, RDKit::QueryAtom const* other)"
        };
        return Err(crate::errors::argument_error_in(
            py,
            "QueryAtom",
            method,
            &got,
            expected,
        ));
    };
    let other: PyRef<'_, PyAtom> = other.as_super().borrow();
    other.with(py, |a| a.query().cloned())
}

#[pymethods]
impl PyQueryAtom {
    #[new]
    #[pyo3(signature = (*_args, **_kwargs))]
    fn py_new(
        _args: &Bound<'_, pyo3::types::PyTuple>,
        _kwargs: Option<&Bound<'_, pyo3::types::PyDict>>,
    ) -> PyResult<PyClassInitializer<Self>> {
        Err(PyRuntimeError::new_err(
            "This class cannot be instantiated from Python",
        ))
    }

    /// combines the query from other with ours
    #[pyo3(name = "ExpandQuery", signature = (other, how = None, maintainOrder = true))]
    #[allow(non_snake_case)]
    fn expand_query(
        slf: &Bound<'_, Self>,
        other: &Bound<'_, PyAny>,
        how: Option<&Bound<'_, PyAny>>,
        maintainOrder: bool,
    ) -> PyResult<()> {
        let py = slf.py();
        let extra: Vec<Bound<'_, PyAny>> = how.into_iter().cloned().collect();
        let how = expand_how(how)?;
        if let Some(query) = other_query(slf, other, "ExpandQuery", &extra)? {
            let mut this: PyRefMut<'_, PyAtom> = slf.as_super().borrow_mut();
            this.with_mut(py, |a| atom_expand_query(a, query, how, maintainOrder))?
                .py()?;
        }
        Ok(())
    }

    /// Replace our query with a copy of the other query
    #[pyo3(name = "SetQuery", signature = (other))]
    fn set_query(slf: &Bound<'_, Self>, other: &Bound<'_, PyAny>) -> PyResult<()> {
        let py = slf.py();
        if let Some(query) = other_query(slf, other, "SetQuery", &[])? {
            let mut this: PyRefMut<'_, PyAtom> = slf.as_super().borrow_mut();
            this.with_mut(py, |a| a.set_query(Some(query)))?;
        }
        Ok(())
    }
}
