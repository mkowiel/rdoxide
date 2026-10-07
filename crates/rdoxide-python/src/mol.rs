//! `rdchem.Mol`, `rdchem.RWMol`, `rdchem.EditableMol` — port of upstream
//! `Code/GraphMol/Wrap/Mol.cpp` and `EditableMol.cpp`.

use pyo3::exceptions::{PyIndexError, PyRuntimeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::graphmol::{Atom, Bond, BondType, ROMol};
use rdoxide::rdgeneral::{Dict, RDProps};

use crate::atom::PyAtom;
use crate::bond::PyBond;
use crate::conformer::PyConformer;
use crate::enums_impl::PyBondType;
use crate::errors::IntoPyResult;
use crate::props::impl_props;
use crate::ringinfo::PyRingInfo;
use crate::stereo_group::PyStereoGroup;

/// The Molecule class.
///
///   In addition to the expected Atoms and Bonds, molecules contain:
///     - a collection of Atom and Bond bookmarks indexed with integers
///         that can be used to flag and retrieve particular Atoms or Bonds
///         using the {get|set}{Atom|Bond}Bookmark() methods.
///
///     - a set of string-valued properties. These can have arbitrary string
///         labels and can be set and retrieved using the {set|get}Prop() methods
///         Molecular properties can be tagged as being *computed*, in which case
///           they will be automatically cleared under certain circumstances (when the
///           molecule itself is modified, for example).
///         Molecules also have the concept of *private* properties, which are tagged
///           by beginning the property name with an underscore (_).
#[pyclass(name = "Mol", module = "rdkit.Chem.rdchem", subclass, weakref)]
#[derive(Debug)]
pub struct PyMol {
    pub mol: ROMol,
}

impl PyMol {
    pub fn with_dict<R>(&self, _py: Python<'_>, f: impl FnOnce(&Dict) -> R) -> PyResult<R> {
        Ok(f(self.mol.dict()))
    }
    pub fn with_dict_mut<R>(
        &mut self,
        _py: Python<'_>,
        f: impl FnOnce(&mut Dict) -> R,
    ) -> PyResult<R> {
        Ok(f(self.mol.dict_mut()))
    }
}

impl_props!(PyMol, computed = true, dict_defaults = (false, false));

/// Arguments accepted by the `Mol(...)` constructor.
#[derive(FromPyObject)]
#[allow(dead_code)] // pickle payloads are decoded once MolPickler lands (stage 6)
enum MolInit<'py> {
    Mol(PyRef<'py, PyMol>),
    Bytes(Vec<u8>),
    Str(String),
}

fn make_mol(other: Option<MolInit<'_>>, quick_copy: bool, conf_id: i32) -> PyResult<ROMol> {
    match other {
        None => Ok(ROMol::new()),
        Some(MolInit::Mol(m)) => Ok(ROMol::copy_from(&m.mol, quick_copy, conf_id)),
        Some(MolInit::Bytes(_)) | Some(MolInit::Str(_)) => {
            // MolPickler is ported in stage 6
            Err(PyRuntimeError::new_err(
                "rdoxide: molecule pickles are not supported yet",
            ))
        }
    }
}

#[pymethods]
impl PyMol {
    #[new]
    #[pyo3(signature = (*args, quickCopy = false, confId = -1))]
    #[allow(non_snake_case)]
    fn py_new(args: &Bound<'_, PyTuple>, quickCopy: bool, confId: i32) -> PyResult<Self> {
        Ok(PyMol {
            mol: Self::init_from_args(args, quickCopy, confId)?,
        })
    }

    fn __copy__(slf: &Bound<'_, Self>) -> PyResult<Py<PyAny>> {
        copy_same_type(slf)
    }

    #[pyo3(signature = (memo))]
    fn __deepcopy__(slf: &Bound<'_, Self>, memo: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        let _ = memo;
        let res = copy_same_type(slf)?;
        // copy the instance __dict__ of Python subclasses
        if let Ok(d) = slf.getattr("__dict__") {
            let copy = slf.py().import("copy")?;
            let nd = copy.call_method1("deepcopy", (d, memo))?;
            res.bind(slf.py())
                .getattr("__dict__")?
                .call_method1("update", (nd,))?;
        }
        Ok(res)
    }

    /// Returns the number of atoms in the molecule.
    ///
    ///   ARGUMENTS:
    ///     - onlyExplicit: (optional) include only explicit atoms (atoms in the molecular graph)
    ///                     defaults to 1.
    ///   NOTE: the onlyHeavy argument is deprecated
    #[pyo3(name = "GetNumAtoms", signature = (onlyHeavy = -1, onlyExplicit = true))]
    #[allow(non_snake_case)]
    fn get_num_atoms(&self, onlyHeavy: i32, onlyExplicit: bool) -> PyResult<u32> {
        if onlyHeavy > -1 {
            rdoxide::rd_log_warning!(
                "WARNING: the onlyHeavy argument to mol.GetNumAtoms() has been deprecated. Please use the onlyExplicit argument instead or mol.GetNumHeavyAtoms() if you want the heavy atom count."
            );
            if onlyHeavy != 0 {
                return Ok(self.mol.num_heavy_atoms());
            }
        }
        if onlyExplicit {
            Ok(self.mol.num_atoms())
        } else {
            self.mol.num_atoms_with_hs().py()
        }
    }

    /// Returns the number of heavy atoms (atomic number >1) in the molecule.
    #[pyo3(name = "GetNumHeavyAtoms")]
    fn get_num_heavy_atoms(&self) -> u32 {
        self.mol.num_heavy_atoms()
    }

    /// Returns a particular Atom.
    #[pyo3(name = "GetAtomWithIdx")]
    fn get_atom_with_idx(slf: &Bound<'_, Self>, idx: u32) -> PyResult<Py<PyAtom>> {
        slf.borrow().mol.get_atom(idx).py()?;
        PyAtom::in_mol(slf.py(), slf.clone().unbind(), idx)
    }

    /// Returns the number of Bonds in the molecule.
    #[pyo3(name = "GetNumBonds", signature = (onlyHeavy = true))]
    #[allow(non_snake_case)]
    fn get_num_bonds(&self, onlyHeavy: bool) -> PyResult<u32> {
        if onlyHeavy {
            Ok(self.mol.num_bonds())
        } else {
            self.mol.num_bonds_with_hs().py()
        }
    }

    /// Returns a particular Bond.
    #[pyo3(name = "GetBondWithIdx")]
    fn get_bond_with_idx(slf: &Bound<'_, Self>, idx: u32) -> PyResult<Py<PyBond>> {
        slf.borrow().mol.get_bond(idx).py()?;
        PyBond::in_mol(slf.py(), slf.clone().unbind(), idx)
    }

    /// Returns the bond between two atoms, if there is one.
    #[pyo3(name = "GetBondBetweenAtoms")]
    fn get_bond_between_atoms(
        slf: &Bound<'_, Self>,
        idx1: u32,
        idx2: u32,
    ) -> PyResult<Option<Py<PyBond>>> {
        let b = slf
            .borrow()
            .mol
            .get_bond_between_atoms(idx1, idx2)
            .py()?
            .map(|b| b.idx());
        b.map(|b| PyBond::in_mol(slf.py(), slf.clone().unbind(), b))
            .transpose()
    }

    /// Returns a read-only sequence containing all of the molecule's Atoms.
    #[pyo3(name = "GetAtoms")]
    fn get_atoms(slf: &Bound<'_, Self>) -> AtomSeq {
        let n = slf.borrow().mol.num_atoms();
        AtomSeq {
            mol: slf.clone().unbind(),
            indices: None,
            len_at_start: n,
            pos: 0,
        }
    }

    /// Returns a read-only sequence containing all of the molecule's Bonds.
    #[pyo3(name = "GetBonds")]
    fn get_bonds(slf: &Bound<'_, Self>) -> BondSeq {
        let n = slf.borrow().mol.num_bonds();
        BondSeq {
            mol: slf.clone().unbind(),
            indices: None,
            len_at_start: n,
            pos: 0,
        }
    }

    /// Returns a read-only sequence containing all of the molecule's aromatic Atoms.
    #[pyo3(name = "GetAromaticAtoms")]
    fn get_aromatic_atoms(slf: &Bound<'_, Self>) -> AtomSeq {
        let m = slf.borrow();
        let idx: Vec<u32> = m
            .mol
            .atoms()
            .iter()
            .filter(|a| a.is_aromatic())
            .map(Atom::idx)
            .collect();
        AtomSeq {
            mol: slf.clone().unbind(),
            len_at_start: m.mol.num_atoms(),
            indices: Some(idx),
            pos: 0,
        }
    }

    /// Returns the number of conformations on the molecule
    #[pyo3(name = "GetNumConformers")]
    fn get_num_conformers(&self) -> u32 {
        self.mol.num_conformers()
    }

    /// Adds a Conformer to the molecule and returns its ID.
    ///   ARGUMENTS:
    ///     - conf: the conformer to add (will be copied).
    ///     - assignId: (optional) toggles assignment of the conformer id (defaults to False)
    #[pyo3(name = "AddConformer", signature = (conf, assignId = false))]
    #[allow(non_snake_case)]
    fn add_conformer(
        &mut self,
        py: Python<'_>,
        conf: PyRef<'_, PyConformer>,
        assignId: bool,
    ) -> PyResult<u32> {
        let c = conf.snapshot(py)?;
        self.mol.add_conformer(c, assignId).py()
    }

    /// Get the conformer with a specified ID
    #[pyo3(name = "GetConformer", signature = (id = -1))]
    fn get_conformer(slf: &Bound<'_, Self>, id: i32) -> PyResult<Py<PyConformer>> {
        let cid = slf
            .borrow()
            .mol
            .conformer(id)
            .map_err(|_| PyValueError::new_err("Bad Conformer Id"))?
            .id();
        PyConformer::in_mol(slf.py(), slf.clone().unbind(), cid)
    }

    /// Returns a read-only sequence containing all of the molecule's Conformers.
    #[pyo3(name = "GetConformers")]
    fn get_conformers(slf: &Bound<'_, Self>) -> PyResult<Vec<Py<PyConformer>>> {
        let ids: Vec<u32> = slf
            .borrow()
            .mol
            .conformers()
            .iter()
            .map(|c| c.id())
            .collect();
        ids.into_iter()
            .map(|id| PyConformer::in_mol(slf.py(), slf.clone().unbind(), id))
            .collect()
    }

    /// Remove all the conformations on the molecule
    #[pyo3(name = "RemoveAllConformers")]
    fn remove_all_conformers(&mut self) {
        self.mol.clear_conformers();
    }

    /// Remove the conformer with the specified ID
    #[pyo3(name = "RemoveConformer")]
    fn remove_conformer(&mut self, id: u32) {
        self.mol.remove_conformer(id);
    }

    /// Queries whether or not the molecule has any queries
    #[pyo3(name = "HasQuery")]
    fn has_query(&self) -> bool {
        self.mol.has_query()
    }

    /// Sets the molecule's name (the `_Name` property).
    #[pyo3(name = "SetName")]
    fn set_name(&mut self, name: &str) {
        self.mol.set_name(name);
    }

    /// Returns the molecule's name (the `_Name` property).
    #[pyo3(name = "GetName")]
    fn get_name(&self) -> String {
        self.mol.name()
    }

    /// Removes all computed properties from the molecule.
    #[pyo3(name = "ClearComputedProps", signature = (includeRings = true))]
    #[allow(non_snake_case)]
    fn clear_computed_props(&mut self, includeRings: bool) {
        self.mol.clear_computed_props(includeRings);
    }

    /// Regenerates computed properties like implicit valence and ring information.
    #[pyo3(name = "UpdatePropertyCache", signature = (strict = true))]
    fn update_property_cache(&mut self, strict: bool) -> PyResult<()> {
        self.mol.update_property_cache(strict).py()
    }

    /// Returns true or false depending on whether implicit and explicit valence of the molecule have already been calculated.
    #[pyo3(name = "NeedsUpdatePropertyCache")]
    fn needs_update_property_cache(&self) -> bool {
        self.mol.needs_update_property_cache()
    }

    /// Clears implicit and explicit valence information from all atoms.
    #[pyo3(name = "ClearPropertyCache")]
    fn clear_property_cache(&mut self) {
        self.mol.clear_property_cache();
    }

    /// Returns a list of StereoGroups defining the relative stereochemistry of the atoms.
    #[pyo3(name = "GetStereoGroups")]
    fn get_stereo_groups(slf: &Bound<'_, Self>) -> Vec<PyStereoGroup> {
        let m = slf.borrow();
        m.mol
            .stereo_groups()
            .iter()
            .map(|g| PyStereoGroup {
                group: g.clone(),
                mol: Some(slf.clone().unbind()),
            })
            .collect()
    }

    /// Returns the number of molecule's RingInfo object.
    #[pyo3(name = "GetRingInfo")]
    fn get_ring_info(slf: &Bound<'_, Self>) -> PyRingInfo {
        PyRingInfo {
            mol: slf.clone().unbind(),
        }
    }

    /// Prints debugging information about the molecule.
    #[pyo3(name = "Debug", signature = (useStdout = true))]
    #[allow(non_snake_case)]
    fn debug(&self, py: Python<'_>, useStdout: bool) -> PyResult<()> {
        let text = crate::debug::debug_mol(&self.mol);
        if useStdout {
            py.import("sys")?
                .getattr("stdout")?
                .call_method1("write", (text,))?;
        } else {
            rdoxide::rdgeneral::log::log(
                rdoxide::rdgeneral::log::Level::Info,
                text.trim_end_matches('\n'),
            );
        }
        Ok(())
    }
}

impl PyMol {
    pub fn init_from_args(
        args: &Bound<'_, PyTuple>,
        quick_copy: bool,
        conf_id: i32,
    ) -> PyResult<ROMol> {
        match args.len() {
            0 => make_mol(None, quick_copy, conf_id),
            1..=3 => {
                let first: MolInit<'_> = args.get_item(0)?.extract()?;
                let qc = if args.len() > 1 {
                    args.get_item(1)?.extract::<bool>()?
                } else {
                    quick_copy
                };
                let cid = if args.len() > 2 {
                    args.get_item(2)?.extract::<i32>()?
                } else {
                    conf_id
                };
                make_mol(Some(first), qc, cid)
            }
            _ => Err(PyTypeError::new_err(
                "Mol() takes at most 3 positional arguments",
            )),
        }
    }
}

/// Copy a Mol (or RWMol / Python subclass) preserving its type.
fn copy_same_type(slf: &Bound<'_, PyMol>) -> PyResult<Py<PyAny>> {
    Ok(slf.get_type().call1((slf,))?.unbind())
}

/// Read-only sequence of atoms (`_ROAtomSeq`).
#[pyclass(name = "_ROAtomSeq", module = "rdkit.Chem.rdchem")]
#[derive(Debug)]
pub struct AtomSeq {
    pub(crate) mol: Py<PyMol>,
    /// explicit atom indices (None: all atoms)
    pub(crate) indices: Option<Vec<u32>>,
    pub(crate) len_at_start: u32,
    pub(crate) pos: usize,
}

impl AtomSeq {
    fn current_len(&self, py: Python<'_>) -> usize {
        match &self.indices {
            Some(v) => v.len(),
            None => self.mol.bind(py).borrow().mol.num_atoms() as usize,
        }
    }
    fn check_modified(&self, py: Python<'_>) -> PyResult<()> {
        if self.mol.bind(py).borrow().mol.num_atoms() != self.len_at_start {
            return Err(PyRuntimeError::new_err(
                "Sequence modified during iteration",
            ));
        }
        Ok(())
    }
    fn item(&self, py: Python<'_>, i: usize) -> PyResult<Py<PyAtom>> {
        let idx = match &self.indices {
            Some(v) => v[i],
            None => i as u32,
        };
        PyAtom::in_mol(py, self.mol.clone_ref(py), idx)
    }
}

#[pymethods]
impl AtomSeq {
    fn __len__(&self, py: Python<'_>) -> usize {
        self.current_len(py)
    }
    fn __getitem__(&self, py: Python<'_>, which: i64) -> PyResult<Py<PyAtom>> {
        let n = self.current_len(py) as i64;
        let i = if which < 0 { which + n } else { which };
        if i < 0 || i >= n {
            return Err(PyIndexError::new_err("End of sequence hit"));
        }
        self.check_modified(py)?;
        self.item(py, i as usize)
    }
    fn __iter__(mut slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf.pos = 0;
        slf
    }
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<PyAtom>>> {
        if self.pos >= self.current_len(py) {
            return Ok(None);
        }
        self.check_modified(py)?;
        let r = self.item(py, self.pos)?;
        self.pos += 1;
        Ok(Some(r))
    }
}

/// Read-only sequence of bonds (`_ROBondSeq`).
#[pyclass(name = "_ROBondSeq", module = "rdkit.Chem.rdchem")]
#[derive(Debug)]
pub struct BondSeq {
    pub(crate) mol: Py<PyMol>,
    pub(crate) indices: Option<Vec<u32>>,
    pub(crate) len_at_start: u32,
    pub(crate) pos: usize,
}

impl BondSeq {
    fn current_len(&self, py: Python<'_>) -> usize {
        match &self.indices {
            Some(v) => v.len(),
            None => self.mol.bind(py).borrow().mol.num_bonds() as usize,
        }
    }
    fn check_modified(&self, py: Python<'_>) -> PyResult<()> {
        if self.indices.is_none() && self.mol.bind(py).borrow().mol.num_bonds() != self.len_at_start
        {
            return Err(PyRuntimeError::new_err(
                "Sequence modified during iteration",
            ));
        }
        Ok(())
    }
    fn item(&self, py: Python<'_>, i: usize) -> PyResult<Py<PyBond>> {
        let idx = match &self.indices {
            Some(v) => v[i],
            None => i as u32,
        };
        PyBond::in_mol(py, self.mol.clone_ref(py), idx)
    }
}

#[pymethods]
impl BondSeq {
    fn __len__(&self, py: Python<'_>) -> usize {
        self.current_len(py)
    }
    fn __getitem__(&self, py: Python<'_>, which: i64) -> PyResult<Py<PyBond>> {
        let n = self.current_len(py) as i64;
        let i = if which < 0 { which + n } else { which };
        if i < 0 || i >= n {
            return Err(PyIndexError::new_err("End of sequence hit"));
        }
        self.check_modified(py)?;
        self.item(py, i as usize)
    }
    fn __iter__(mut slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf.pos = 0;
        slf
    }
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<PyBond>>> {
        if self.pos >= self.current_len(py) {
            return Ok(None);
        }
        self.check_modified(py)?;
        let r = self.item(py, self.pos)?;
        self.pos += 1;
        Ok(Some(r))
    }
}

/// The RW molecule class (read/write)
///
///   This class is a more-performant version of the EditableMolecule class in that
///   it is a 'live' molecule and shares the interface from the Mol class.
///   All changes are performed without the need to create a copy of the
///   molecule using GetMol() (this is still available, however).
///
///   n.b. Eventually this class may become a direct replacement for EditableMol
#[pyclass(name = "RWMol", module = "rdkit.Chem.rdchem", extends = PyMol, subclass)]
#[derive(Debug)]
pub struct PyRWMol;

fn rw<'a>(slf: &'a Bound<'a, PyRWMol>) -> PyRefMut<'a, PyMol> {
    slf.as_super().borrow_mut()
}

#[pymethods]
impl PyRWMol {
    #[new]
    #[pyo3(signature = (*args, quickCopy = false, confId = -1))]
    #[allow(non_snake_case)]
    fn py_new(
        args: &Bound<'_, PyTuple>,
        quickCopy: bool,
        confId: i32,
    ) -> PyResult<PyClassInitializer<Self>> {
        let base = PyMol {
            mol: PyMol::init_from_args(args, quickCopy, confId)?,
        };
        Ok(PyClassInitializer::from(base).add_subclass(PyRWMol))
    }

    fn __copy__(slf: &Bound<'_, Self>) -> PyResult<Py<PyAny>> {
        Ok(slf.get_type().call1((slf,))?.unbind())
    }

    #[pyo3(signature = (memo))]
    fn __deepcopy__(slf: &Bound<'_, Self>, memo: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        let _ = memo;
        Ok(slf.get_type().call1((slf,))?.unbind())
    }

    fn __enter__(slf: &Bound<'_, Self>) -> PyResult<Py<Self>> {
        rw(slf).mol.begin_batch_edit().py()?;
        Ok(slf.clone().unbind())
    }

    #[pyo3(signature = (exc_type, exc_val, traceback))]
    fn __exit__(
        slf: &Bound<'_, Self>,
        exc_type: &Bound<'_, PyAny>,
        exc_val: &Bound<'_, PyAny>,
        traceback: &Bound<'_, PyAny>,
    ) -> PyResult<bool> {
        let _ = (exc_val, traceback);
        let mut m = rw(slf);
        if exc_type.is_none() {
            m.mol.commit_batch_edit().py()?;
        } else {
            m.mol.rollback_batch_edit();
        }
        Ok(false)
    }

    /// Remove the specified atom from the molecule
    #[pyo3(name = "RemoveAtom")]
    fn remove_atom(slf: &Bound<'_, Self>, idx: u32) -> PyResult<()> {
        rw(slf).mol.remove_atom(idx).py()
    }

    /// Remove the specified bond from the molecule
    #[pyo3(name = "RemoveBond")]
    fn remove_bond(slf: &Bound<'_, Self>, idx1: u32, idx2: u32) -> PyResult<()> {
        rw(slf).mol.remove_bond(idx1, idx2).py()
    }

    /// add a bond, returns the new number of bonds
    #[pyo3(name = "AddBond", signature = (beginAtomIdx, endAtomIdx, order = None))]
    #[allow(non_snake_case)]
    fn add_bond(
        slf: &Bound<'_, Self>,
        beginAtomIdx: u32,
        endAtomIdx: u32,
        order: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<u32> {
        let bt = order
            .map(PyBondType::extract)
            .transpose()?
            .unwrap_or(BondType::Unspecified);
        rw(slf).mol.add_bond(beginAtomIdx, endAtomIdx, bt).py()
    }

    /// add an atom, returns the index of the newly added atom
    #[pyo3(name = "AddAtom", signature = (atom = None))]
    fn add_atom(
        slf: &Bound<'_, Self>,
        py: Python<'_>,
        atom: Option<PyRef<'_, PyAtom>>,
    ) -> PyResult<u32> {
        let a = match atom {
            Some(a) => a.snapshot(py)?,
            None => return Err(PyValueError::new_err("got None for atom")),
        };
        Ok(rw(slf).mol.add_atom(a.copy_detached(), true))
    }

    /// replaces the specified atom with the provided one
    /// If updateLabel is True, the new atom becomes the active atom
    /// If preserveProps is True preserve keep the existing props unless explicit set on the new atom
    #[pyo3(name = "ReplaceAtom", signature = (index, newAtom, updateLabel = false, preserveProps = false))]
    #[allow(non_snake_case)]
    fn replace_atom(
        slf: &Bound<'_, Self>,
        py: Python<'_>,
        index: u32,
        newAtom: PyRef<'_, PyAtom>,
        updateLabel: bool,
        preserveProps: bool,
    ) -> PyResult<()> {
        let _ = updateLabel;
        let a = newAtom.snapshot(py)?;
        rw(slf).mol.replace_atom(index, &a, preserveProps).py()
    }

    /// replaces the specified bond with the provided one.
    /// If preserveProps is True preserve keep the existing props unless explicit set on the new bond
    /// If keepSGroups is False, all Substance Groups referencing the bond will be dropped.
    #[pyo3(name = "ReplaceBond", signature = (index, newBond, preserveProps = false, keepSGroups = true))]
    #[allow(non_snake_case)]
    fn replace_bond(
        slf: &Bound<'_, Self>,
        py: Python<'_>,
        index: u32,
        newBond: PyRef<'_, PyBond>,
        preserveProps: bool,
        keepSGroups: bool,
    ) -> PyResult<()> {
        let b = newBond.snapshot(py)?;
        rw(slf)
            .mol
            .replace_bond(index, &b, preserveProps, keepSGroups)
            .py()
    }

    /// Returns a Mol (a normal molecule)
    #[pyo3(name = "GetMol")]
    fn get_mol(slf: &Bound<'_, Self>) -> PyMol {
        PyMol {
            mol: ROMol::copy_from(&slf.as_super().borrow().mol, false, -1),
        }
    }

    /// Set the stereo groups
    #[pyo3(name = "SetStereoGroups")]
    fn set_stereo_groups(
        slf: &Bound<'_, Self>,
        stereo_groups: Vec<PyRef<'_, PyStereoGroup>>,
    ) -> PyResult<()> {
        let groups = stereo_groups.iter().map(|g| g.group.clone()).collect();
        rw(slf).mol.set_stereo_groups(groups).py()
    }

    /// Insert (add) the given molecule into this molecule
    #[pyo3(name = "InsertMol")]
    fn insert_mol(slf: &Bound<'_, Self>, mol: &Bound<'_, PyMol>) -> PyResult<()> {
        let other = ROMol::copy_from(&mol.borrow().mol, false, -1);
        rw(slf).mol.insert_mol(&other).py()
    }

    /// starts batch editing
    #[pyo3(name = "BeginBatchEdit")]
    fn begin_batch_edit(slf: &Bound<'_, Self>) -> PyResult<()> {
        rw(slf).mol.begin_batch_edit().py()
    }

    /// cancels batch editing
    #[pyo3(name = "RollbackBatchEdit")]
    fn rollback_batch_edit(slf: &Bound<'_, Self>) {
        rw(slf).mol.rollback_batch_edit();
    }

    /// finishes batch editing and makes the actual changes
    #[pyo3(name = "CommitBatchEdit")]
    fn commit_batch_edit(slf: &Bound<'_, Self>) -> PyResult<()> {
        rw(slf).mol.commit_batch_edit().py()
    }
}

/// an editable molecule class
#[pyclass(name = "EditableMol", module = "rdkit.Chem.rdchem")]
#[derive(Debug)]
pub struct PyEditableMol {
    mol: Option<ROMol>,
}

impl PyEditableMol {
    fn m(&mut self) -> PyResult<&mut ROMol> {
        self.mol
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("Pre-condition Violation\n\tno molecule"))
    }
}

#[pymethods]
impl PyEditableMol {
    #[new]
    fn py_new(m: PyRef<'_, PyMol>) -> Self {
        PyEditableMol {
            mol: Some(ROMol::copy_from(&m.mol, false, -1)),
        }
    }

    /// Remove the specified atom from the molecule
    #[pyo3(name = "RemoveAtom")]
    fn remove_atom(&mut self, idx: u32) -> PyResult<()> {
        self.m()?.remove_atom(idx).py()
    }

    /// Remove the specified bond from the molecule
    #[pyo3(name = "RemoveBond")]
    fn remove_bond(&mut self, idx1: u32, idx2: u32) -> PyResult<()> {
        self.m()?.remove_bond(idx1, idx2).py()
    }

    /// add a bond, returns the index of the newly added bond
    #[pyo3(name = "AddBond", signature = (beginAtomIdx, endAtomIdx, order = None))]
    #[allow(non_snake_case)]
    fn add_bond(
        &mut self,
        beginAtomIdx: u32,
        endAtomIdx: u32,
        order: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<u32> {
        let bt = order
            .map(PyBondType::extract)
            .transpose()?
            .unwrap_or(BondType::Unspecified);
        self.m()?.add_bond(beginAtomIdx, endAtomIdx, bt).py()
    }

    /// add an atom, returns the index of the newly added atom
    #[pyo3(name = "AddAtom")]
    fn add_atom(&mut self, py: Python<'_>, atom: PyRef<'_, PyAtom>) -> PyResult<u32> {
        let a = atom.snapshot(py)?;
        Ok(self.m()?.add_atom(a.copy_detached(), true))
    }

    /// replaces the specified atom with the provided one
    #[pyo3(name = "ReplaceAtom", signature = (index, newAtom, updateLabel = false, preserveProps = false))]
    #[allow(non_snake_case)]
    fn replace_atom(
        &mut self,
        py: Python<'_>,
        index: u32,
        newAtom: PyRef<'_, PyAtom>,
        updateLabel: bool,
        preserveProps: bool,
    ) -> PyResult<()> {
        let _ = updateLabel;
        let a = newAtom.snapshot(py)?;
        self.m()?.replace_atom(index, &a, preserveProps).py()
    }

    /// replaces the specified bond with the provided one
    #[pyo3(name = "ReplaceBond", signature = (index, newBond, preserveProps = false))]
    #[allow(non_snake_case)]
    fn replace_bond(
        &mut self,
        py: Python<'_>,
        index: u32,
        newBond: PyRef<'_, PyBond>,
        preserveProps: bool,
    ) -> PyResult<()> {
        let b: Bond = newBond.snapshot(py)?;
        self.m()?.replace_bond(index, &b, preserveProps, true).py()
    }

    /// starts batch editing
    #[pyo3(name = "BeginBatchEdit")]
    fn begin_batch_edit(&mut self) -> PyResult<()> {
        self.m()?.begin_batch_edit().py()
    }

    /// cancels batch editing
    #[pyo3(name = "RollbackBatchEdit")]
    fn rollback_batch_edit(&mut self) -> PyResult<()> {
        self.m()?.rollback_batch_edit();
        Ok(())
    }

    /// finishes batch editing and makes the actual edits
    #[pyo3(name = "CommitBatchEdit")]
    fn commit_batch_edit(&mut self) -> PyResult<()> {
        self.m()?.commit_batch_edit().py()
    }

    /// Returns a Mol (a normal molecule)
    #[pyo3(name = "GetMol")]
    fn get_mol(&mut self) -> PyResult<PyMol> {
        let m = self.m()?;
        Ok(PyMol {
            mol: ROMol::copy_from(m, false, -1),
        })
    }
}
