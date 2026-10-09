//! `rdchem.Bond` — port of upstream `Code/GraphMol/Wrap/Bond.cpp`.

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use rdoxide::graphmol::queries::BondQuery;
use rdoxide::graphmol::query_ops::bond_expand_query;
use rdoxide::graphmol::{Atom, Bond, ROMol};
use rdoxide::rdgeneral::Dict;

use crate::atom::{AtomHandle, PyAtom};
use crate::enums_impl::{PyBondDir, PyBondStereo, PyBondType};
use crate::errors::IntoPyResult;
use crate::mol::PyMol;
use crate::props::impl_props;

#[derive(Debug)]
pub enum BondHandle {
    /// Bonds cannot be created from Python, but query bonds can.
    Owned(Bond),
    InMol {
        mol: Py<PyMol>,
        idx: u32,
    },
}

/// The class to store Bonds.
/// Note: unlike Atoms, is it currently impossible to construct Bonds from
/// Python.
#[pyclass(name = "Bond", module = "rdkit.Chem.rdchem", subclass)]
#[derive(Debug)]
pub struct PyBond {
    pub handle: BondHandle,
}

fn gone() -> PyErr {
    PyRuntimeError::new_err("bond is no longer part of its molecule")
}

fn no_owner() -> PyErr {
    crate::errors::to_py(rdoxide::RDKitError::precondition(
        "no owning molecule for bond",
    ))
}

/// The class to store QueryBonds.
/// These cannot currently be constructed directly from Python
#[pyclass(name = "QueryBond", module = "rdkit.Chem.rdchem", extends = PyBond)]
pub struct PyQueryBond;

impl PyBond {
    /// The Python object for a bond: a `QueryBond` when the bond carries a
    /// query, else a `Bond` (see `PyAtom::from_handle`).
    pub fn from_handle(py: Python<'_>, handle: BondHandle) -> PyResult<Py<PyBond>> {
        let is_query = match &handle {
            BondHandle::Owned(b) => b.has_query(),
            BondHandle::InMol { mol, idx } => mol.bind(py).try_borrow().is_ok_and(|m| {
                m.mol
                    .bonds()
                    .get(*idx as usize)
                    .is_some_and(Bond::has_query)
            }),
        };
        let init = PyClassInitializer::from(PyBond { handle });
        if is_query {
            let obj = Py::new(py, init.add_subclass(PyQueryBond))?;
            Ok(obj
                .into_bound(py)
                .into_any()
                .cast_into::<PyBond>()?
                .unbind())
        } else {
            Py::new(py, init)
        }
    }

    pub fn in_mol(py: Python<'_>, mol: Py<PyMol>, idx: u32) -> PyResult<Py<PyBond>> {
        Self::from_handle(py, BondHandle::InMol { mol, idx })
    }

    /// A free-standing bond.
    pub fn owned(py: Python<'_>, bond: Bond) -> PyResult<Py<PyBond>> {
        Self::from_handle(py, BondHandle::Owned(bond))
    }

    pub fn with<R>(&self, py: Python<'_>, f: impl FnOnce(&Bond) -> R) -> PyResult<R> {
        match &self.handle {
            BondHandle::Owned(b) => Ok(f(b)),
            BondHandle::InMol { mol, idx } => {
                let m = mol.bind(py).borrow();
                m.mol.bonds().get(*idx as usize).map(f).ok_or_else(gone)
            }
        }
    }

    pub fn with_mut<R>(&mut self, py: Python<'_>, f: impl FnOnce(&mut Bond) -> R) -> PyResult<R> {
        match &mut self.handle {
            BondHandle::Owned(b) => Ok(f(b)),
            BondHandle::InMol { mol, idx } => {
                let mut m = mol.bind(py).borrow_mut();
                if *idx >= m.mol.num_bonds() {
                    return Err(gone());
                }
                Ok(f(m.mol.bond_mut(*idx)))
            }
        }
    }

    fn with_mol<R>(&self, py: Python<'_>, f: impl FnOnce(&ROMol, u32) -> R) -> PyResult<R> {
        match &self.handle {
            BondHandle::Owned(_) => Err(no_owner()),
            BondHandle::InMol { mol, idx } => {
                let m = mol.bind(py).borrow();
                if *idx >= m.mol.num_bonds() {
                    return Err(gone());
                }
                Ok(f(&m.mol, *idx))
            }
        }
    }

    pub(crate) fn with_mol_mut<R>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut ROMol, u32) -> R,
    ) -> PyResult<R> {
        match &self.handle {
            BondHandle::Owned(_) => Err(no_owner()),
            BondHandle::InMol { mol, idx } => {
                let mut m = mol.bind(py).borrow_mut();
                if *idx >= m.mol.num_bonds() {
                    return Err(gone());
                }
                Ok(f(&mut m.mol, *idx))
            }
        }
    }

    fn owner(&self, py: Python<'_>) -> PyResult<Py<PyMol>> {
        match &self.handle {
            BondHandle::InMol { mol, .. } => Ok(mol.clone_ref(py)),
            BondHandle::Owned(_) => Err(no_owner()),
        }
    }

    pub fn snapshot(&self, py: Python<'_>) -> PyResult<Bond> {
        self.with(py, Bond::clone)
    }

    pub fn with_dict<R>(&self, py: Python<'_>, f: impl FnOnce(&Dict) -> R) -> PyResult<R> {
        self.with(py, |b| f(rdoxide::rdgeneral::RDProps::dict(b)))
    }
    pub fn with_dict_mut<R>(
        &mut self,
        py: Python<'_>,
        f: impl FnOnce(&mut Dict) -> R,
    ) -> PyResult<R> {
        self.with_mut(py, |b| f(rdoxide::rdgeneral::RDProps::dict_mut(b)))
    }

    fn atom_idx_of(py: Python<'_>, atom: &PyAtom) -> PyResult<u32> {
        atom.with(py, rdoxide::graphmol::Atom::idx)
    }
}

impl_props!(PyBond, computed = false, dict_defaults = (true, true));

#[pymethods]
impl PyBond {
    /// Returns whether or not this instance belongs to a molecule.
    #[pyo3(name = "HasOwningMol")]
    fn has_owning_mol(&self) -> bool {
        matches!(self.handle, BondHandle::InMol { .. })
    }

    /// Returns the Mol that owns this bond.
    #[pyo3(name = "GetOwningMol")]
    fn get_owning_mol(&self, py: Python<'_>) -> PyResult<Py<PyMol>> {
        self.owner(py)
    }

    /// Returns the type of the bond as a BondType
    #[pyo3(name = "GetBondType")]
    fn get_bond_type(&self, py: Python<'_>) -> PyResult<PyBondType> {
        self.with(py, |b| PyBondType(b.bond_type()))
    }

    /// Set the type of the bond as a BondType
    #[pyo3(name = "SetBondType")]
    #[allow(non_snake_case)]
    fn set_bond_type(&mut self, py: Python<'_>, bT: &Bound<'_, PyAny>) -> PyResult<()> {
        let t = PyBondType::extract(bT)?;
        self.with_mut(py, |b| b.set_bond_type(t))
    }

    /// Returns the type of the bond as a double (i.e. 1.0 for SINGLE, 1.5 for AROMATIC, 2.0 for DOUBLE)
    #[pyo3(name = "GetBondTypeAsDouble")]
    fn get_bond_type_as_double(&self, py: Python<'_>) -> PyResult<f64> {
        self.with(py, Bond::bond_type_as_double)?.py()
    }

    /// Returns the type of the bond as a BondDir
    #[pyo3(name = "GetBondDir")]
    fn get_bond_dir(&self, py: Python<'_>) -> PyResult<PyBondDir> {
        self.with(py, |b| PyBondDir(b.bond_dir()))
    }

    /// Set the type of the bond as a BondDir
    #[pyo3(name = "SetBondDir")]
    fn set_bond_dir(&mut self, py: Python<'_>, what: &Bound<'_, PyAny>) -> PyResult<()> {
        let d = PyBondDir::extract(what)?;
        self.with_mut(py, |b| b.set_bond_dir(d))
    }

    /// Returns the stereo configuration of the bond as a BondStereo
    #[pyo3(name = "GetStereo")]
    fn get_stereo(&self, py: Python<'_>) -> PyResult<PyBondStereo> {
        self.with(py, |b| PyBondStereo(b.stereo()))
    }

    /// Set the stereo configuration of the bond as a BondStereo
    #[pyo3(name = "SetStereo")]
    fn set_stereo(&mut self, py: Python<'_>, what: &Bound<'_, PyAny>) -> PyResult<()> {
        let s = PyBondStereo::extract(what)?;
        self.with_mut(py, |b| b.set_stereo(s))?.py()
    }

    /// Returns the indices of the atoms setting this bond's stereochemistry.
    #[pyo3(name = "GetStereoAtoms")]
    fn get_stereo_atoms(&self, py: Python<'_>) -> PyResult<Vec<i32>> {
        self.with(py, |b| b.stereo_atoms().to_vec())
    }

    /// Set the indices of the atoms setting this bond's stereochemistry.
    #[pyo3(name = "SetStereoAtoms")]
    #[allow(non_snake_case)]
    fn set_stereo_atoms(&mut self, py: Python<'_>, bgnIdx: u32, endIdx: u32) -> PyResult<()> {
        self.with_mol_mut(py, |m, i| m.set_bond_stereo_atoms(i, bgnIdx, endIdx))?
            .py()
    }

    /// Inverts the chirality of an atropisomeric bond.
    #[pyo3(name = "InvertChirality")]
    fn invert_chirality(&mut self, py: Python<'_>) -> PyResult<bool> {
        self.with_mut(py, Bond::invert_chirality)
    }

    /// Returns the contribution of the bond to the valence of an Atom.
    ///
    ///   ARGUMENTS:
    ///
    ///     - atom: the Atom to consider.
    #[pyo3(name = "GetValenceContrib")]
    fn get_valence_contrib(&self, py: Python<'_>, at: PyRef<'_, PyAtom>) -> PyResult<f64> {
        let idx = Self::atom_idx_of(py, &at)?;
        let same_mol = match (&self.handle, &at.handle) {
            (BondHandle::InMol { mol: m1, .. }, AtomHandle::InMol { mol: m2, .. }) => m1.is(m2),
            _ => false,
        };
        if !same_mol {
            return Ok(0.0);
        }
        self.with_mol(py, |m, b| m.bond(b).valence_contrib(m, idx))
    }

    /// Returns whether or not the bond is aromatic.
    #[pyo3(name = "GetIsAromatic")]
    fn get_is_aromatic(&self, py: Python<'_>) -> PyResult<bool> {
        self.with(py, Bond::is_aromatic)
    }

    #[pyo3(name = "SetIsAromatic")]
    fn set_is_aromatic(&mut self, py: Python<'_>, what: bool) -> PyResult<()> {
        self.with_mut(py, |b| b.set_is_aromatic(what))
    }

    /// Returns whether or not the bond is considered to be conjugated.
    #[pyo3(name = "GetIsConjugated")]
    fn get_is_conjugated(&self, py: Python<'_>) -> PyResult<bool> {
        self.with(py, Bond::is_conjugated)
    }

    #[pyo3(name = "SetIsConjugated")]
    fn set_is_conjugated(&mut self, py: Python<'_>, what: bool) -> PyResult<()> {
        self.with_mut(py, |b| b.set_is_conjugated(what))
    }

    /// Returns the bond's index (ordering in the molecule)
    #[pyo3(name = "GetIdx")]
    fn get_idx(&self, py: Python<'_>) -> PyResult<u32> {
        self.with(py, Bond::idx)
    }

    /// Returns the index of the bond's first atom.
    #[pyo3(name = "GetBeginAtomIdx")]
    fn get_begin_atom_idx(&self, py: Python<'_>) -> PyResult<u32> {
        self.with(py, Bond::begin_atom_idx)
    }

    /// Returns the index of the bond's first atom.
    #[pyo3(name = "GetEndAtomIdx")]
    fn get_end_atom_idx(&self, py: Python<'_>) -> PyResult<u32> {
        self.with(py, Bond::end_atom_idx)
    }

    /// Given the index of one of the bond's atoms, returns the
    /// index of the other.
    #[pyo3(name = "GetOtherAtomIdx")]
    #[allow(non_snake_case)]
    fn get_other_atom_idx(&self, py: Python<'_>, thisIdx: u32) -> PyResult<u32> {
        self.with(py, |b| b.other_atom_idx(thisIdx))?.py()
    }

    /// Returns the bond's first atom.
    #[pyo3(name = "GetBeginAtom")]
    fn get_begin_atom(&self, py: Python<'_>) -> PyResult<Py<PyAtom>> {
        let mol = self.owner(py)?;
        let i = self.with(py, Bond::begin_atom_idx)?;
        PyAtom::in_mol(py, mol, i)
    }

    /// Returns the bond's second atom.
    #[pyo3(name = "GetEndAtom")]
    fn get_end_atom(&self, py: Python<'_>) -> PyResult<Py<PyAtom>> {
        let mol = self.owner(py)?;
        let i = self.with(py, Bond::end_atom_idx)?;
        PyAtom::in_mol(py, mol, i)
    }

    /// Given one of the bond's atoms, returns the other one.
    #[pyo3(name = "GetOtherAtom")]
    fn get_other_atom(&self, py: Python<'_>, what: PyRef<'_, PyAtom>) -> PyResult<Py<PyAtom>> {
        let mol = self.owner(py)?;
        let idx = Self::atom_idx_of(py, &what)?;
        let other = self.with(py, |b| b.other_atom_idx(idx))?.py()?;
        PyAtom::in_mol(py, mol, other)
    }

    /// Returns whether or not this bond matches another Bond.
    ///
    ///   Each Bond (or query Bond) has a query function which is
    ///   used for this type of matching.
    ///
    ///   ARGUMENTS:
    ///     - other: the other Bond to which to compare
    #[pyo3(name = "Match")]
    fn py_match(&self, py: Python<'_>, what: PyRef<'_, PyBond>) -> PyResult<bool> {
        let Some(query) = self.with(py, |b| b.query().cloned())? else {
            let other = what.snapshot(py)?;
            return self.with(py, |b| b.matches(&other));
        };
        // a query bond evaluates its query on the other bond in its context
        match &what.handle {
            BondHandle::InMol { mol, idx } => {
                let m = mol.bind(py).borrow();
                query.matches(&m.mol, *idx).py()
            }
            BondHandle::Owned(bond) => {
                let mut single = ROMol::new();
                single.add_atom(Atom::new(0), false);
                single.add_atom(Atom::new(0), false);
                let mut copy = bond.clone();
                copy.set_begin_atom_idx_unchecked(0);
                copy.set_end_atom_idx_unchecked(1);
                let idx = single.add_bond_object(copy).py()?;
                query.matches(&single, idx).py()
            }
        }
    }

    /// Returns whether or not the bond is in a ring of a particular size.
    ///
    ///   ARGUMENTS:
    ///     - size: the ring size to look for
    #[pyo3(name = "IsInRingSize")]
    fn is_in_ring_size(&self, py: Python<'_>, size: i32) -> PyResult<bool> {
        self.ensure_rings(py)?;
        self.with_mol(py, |m, i| {
            m.ring_info().is_bond_in_ring_of_size(i, size as u32)
        })?
        .py()
    }

    /// Returns whether or not the bond is in a ring of any size.
    #[pyo3(name = "IsInRing")]
    fn is_in_ring(&self, py: Python<'_>) -> PyResult<bool> {
        self.ensure_rings(py)?;
        Ok(self
            .with_mol(py, |m, i| m.ring_info().num_bond_rings(i))?
            .py()?
            != 0)
    }

    /// Returns whether or not the bond has an associated query
    /// returns the SMARTS (or SMILES) string for a Bond
    #[pyo3(name = "GetSmarts", signature = (allBondsExplicit = false))]
    #[allow(non_snake_case)]
    fn get_smarts(&self, py: Python<'_>, allBondsExplicit: bool) -> PyResult<String> {
        use rdoxide::graphmol::smiles_parse::smarts_write::get_bond_smarts;
        use rdoxide::graphmol::smiles_parse::write::{SmilesWriteParams, get_bond_smiles};
        if self.with(py, Bond::has_query)? {
            return self
                .with(py, |b| {
                    get_bond_smarts(b, &SmilesWriteParams::default(), -1)
                })?
                .py();
        }
        let ps = SmilesWriteParams {
            do_kekule: false,
            all_bonds_explicit: allBondsExplicit,
            do_isomeric_smiles: false,
            ..SmilesWriteParams::default()
        };
        self.with_mol_mut(py, |m, i| get_bond_smiles(m, i, &ps, -1))?
            .py()
    }

    /// returns a text description of the query. Primarily intended for debugging purposes.
    #[pyo3(name = "DescribeQuery")]
    fn describe_query(&self, py: Python<'_>) -> PyResult<String> {
        self.with(py, |b| b.query().map_or(String::new(), |q| q.describe(0)))
    }

    #[pyo3(name = "HasQuery")]
    fn has_query(&self, py: Python<'_>) -> PyResult<bool> {
        self.with(py, Bond::has_query)
    }

    /// Sets the query on this bond (plain bonds have no query).
    #[pyo3(name = "SetQuery")]
    fn set_query(&mut self, other: Option<PyRef<'_, PyBond>>) -> PyResult<()> {
        if other.is_none() {
            return Err(PyValueError::new_err("other Bond is null"));
        }
        Err(crate::errors::to_py(rdoxide::RDKitError::precondition(
            "plain bonds have no Query",
        )))
    }
}

impl PyBond {
    fn ensure_rings(&self, py: Python<'_>) -> PyResult<()> {
        let ok = self.with_mol(py, |m, _| m.ring_info().is_sssr_or_better())?;
        if ok {
            return Ok(());
        }
        self.with_mol_mut(py, |m, _| crate::perception::find_sssr_if_needed(m))?
            .py()
    }
}

fn other_query(
    slf: &Bound<'_, PyQueryBond>,
    other: &Bound<'_, PyAny>,
    method: &str,
    extra: &[Bound<'_, PyAny>],
) -> PyResult<Option<BondQuery>> {
    let py = slf.py();
    let Ok(other) = other.cast::<PyQueryBond>() else {
        let mut got = vec![
            "QueryBond".to_string(),
            crate::errors::boost_type_name(other),
        ];
        got.extend(extra.iter().map(crate::errors::boost_type_name));
        let expected = if method == "ExpandQuery" {
            "ExpandQuery(RDKit::QueryBond* self, RDKit::QueryBond const* other, Queries::CompositeQueryType how=rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND, bool maintainOrder=True)"
        } else {
            "SetQuery(RDKit::QueryBond* self, RDKit::QueryBond const* other)"
        };
        return Err(crate::errors::argument_error_in(
            py,
            "QueryBond",
            method,
            &got,
            expected,
        ));
    };
    let other: PyRef<'_, PyBond> = other.as_super().borrow();
    other.with(py, |a| a.query().cloned())
}

#[pymethods]
impl PyQueryBond {
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
        let how = crate::atom::expand_how(how)?;
        if let Some(query) = other_query(slf, other, "ExpandQuery", &extra)? {
            let mut this: PyRefMut<'_, PyBond> = slf.as_super().borrow_mut();
            this.with_mut(py, |b| bond_expand_query(b, query, how, maintainOrder))?
                .py()?;
        }
        Ok(())
    }

    /// Replace our query with a copy of the other query
    #[pyo3(name = "SetQuery", signature = (other))]
    fn set_query(slf: &Bound<'_, Self>, other: &Bound<'_, PyAny>) -> PyResult<()> {
        let py = slf.py();
        if let Some(query) = other_query(slf, other, "SetQuery", &[])? {
            let mut this: PyRefMut<'_, PyBond> = slf.as_super().borrow_mut();
            this.with_mut(py, |b| b.set_query(Some(query)))?;
        }
        Ok(())
    }
}
