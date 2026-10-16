//! `rdchem.SubstanceGroup` — port of upstream
//! `Code/GraphMol/Wrap/SubstanceGroup.cpp`.
//!
//! `GetMolSubstanceGroups` returns copies of the groups while
//! `GetMolSubstanceGroupWithIdx` and the `Create`/`Add` functions return
//! references into the molecule, as upstream does; both are `SubstanceGroup`
//! objects (see [`Handle`]).

use pyo3::exceptions::{PyIndexError, PyKeyError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyList, PySlice, PyTuple};
use rdoxide::RDKitError;
use rdoxide::geometry::Point3D;
use rdoxide::graphmol::ROMol;
use rdoxide::graphmol::substance_group::{
    AttachPoint, Bracket, CState, SubstanceGroup, add_substance_group,
};
use rdoxide::rdgeneral::{Dict, RDProps, RDValue};

use crate::errors::{IntoPyResult, no_converter, to_py};
use crate::geometry::PyPoint3D;
use crate::mol::PyMol;
use crate::props::impl_props;
use crate::vects::{PyVectJ, PyVectStr};

/// Where the group lives.
#[derive(Debug)]
enum Handle {
    /// A detached copy that still knows its molecule (what
    /// `GetMolSubstanceGroups` returns).
    Copy {
        mol: Py<PyMol>,
        group: SubstanceGroup,
    },
    /// A group inside the molecule.
    InMol { mol: Py<PyMol>, idx: usize },
}

/// A collection of atoms and bonds with associated properties
#[pyclass(name = "SubstanceGroup", module = "rdoxide.Chem.rdchem")]
#[derive(Debug)]
pub struct PySubstanceGroup {
    handle: Handle,
}

fn gone() -> PyErr {
    PyRuntimeError::new_err("SubstanceGroup is no longer part of its molecule")
}

impl PySubstanceGroup {
    fn copy_of(mol: Py<PyMol>, group: SubstanceGroup) -> Self {
        Self {
            handle: Handle::Copy { mol, group },
        }
    }

    fn in_mol(mol: Py<PyMol>, idx: usize) -> Self {
        Self {
            handle: Handle::InMol { mol, idx },
        }
    }

    fn mol(&self) -> &Py<PyMol> {
        match &self.handle {
            Handle::Copy { mol, .. } | Handle::InMol { mol, .. } => mol,
        }
    }

    fn read<R>(&self, py: Python<'_>, f: impl FnOnce(&ROMol, &SubstanceGroup) -> R) -> PyResult<R> {
        let m = self.mol().bind(py).borrow();
        let group = match &self.handle {
            Handle::Copy { group, .. } => group,
            Handle::InMol { idx, .. } => m.mol.substance_groups().get(*idx).ok_or_else(gone)?,
        };
        Ok(f(&m.mol, group))
    }

    /// Runs `f` on the group and the molecule it validates against. A group
    /// inside the molecule is edited as a copy, which is stored back when
    /// `f` succeeds.
    fn edit<R>(
        &mut self,
        py: Python<'_>,
        f: impl FnOnce(&ROMol, &mut SubstanceGroup) -> rdoxide::Result<R>,
    ) -> PyResult<R> {
        match &mut self.handle {
            Handle::Copy { mol, group } => f(&mol.bind(py).borrow().mol, group).py(),
            Handle::InMol { mol, idx } => {
                let mut m = mol.bind(py).borrow_mut();
                let mut group = m.mol.substance_groups().get(*idx).ok_or_else(gone)?.clone();
                let res = f(&m.mol, &mut group).py()?;
                m.mol.substance_groups_mut()[*idx] = group;
                Ok(res)
            }
        }
    }

    pub fn with_dict<R>(&self, py: Python<'_>, f: impl FnOnce(&Dict) -> R) -> PyResult<R> {
        self.read(py, |_, g| f(g.dict()))
    }

    pub fn with_dict_mut<R>(
        &mut self,
        py: Python<'_>,
        f: impl FnOnce(&mut Dict) -> R,
    ) -> PyResult<R> {
        match &mut self.handle {
            Handle::Copy { group, .. } => Ok(f(group.dict_mut())),
            Handle::InMol { mol, idx } => {
                let mut m = mol.bind(py).borrow_mut();
                let group = m
                    .mol
                    .substance_groups_mut()
                    .get_mut(*idx)
                    .ok_or_else(gone)?;
                Ok(f(group.dict_mut()))
            }
        }
    }
}

impl_props!(
    PySubstanceGroup,
    computed = true,
    dict_defaults = (true, true),
    has_prop = bool
);

/// Reads a Python iterable of unsigned integers (`pythonObjectToVect`).
fn index_vec(iterable: &Bound<'_, PyAny>) -> PyResult<Vec<u32>> {
    iterable
        .try_iter()?
        .map(|item| {
            let item = item?;
            item.extract::<u32>()
                .map_err(|_| no_converter("unsigned int", &item))
        })
        .collect()
}

fn point(item: &Bound<'_, PyAny>) -> PyResult<Point3D> {
    Ok(item.extract::<PyRef<'_, PyPoint3D>>()?.p)
}

#[pymethods]
impl PySubstanceGroup {
    /// returns the molecule owning this SubstanceGroup
    #[pyo3(name = "GetOwningMol")]
    fn get_owning_mol(&self, py: Python<'_>) -> Py<PyMol> {
        self.mol().clone_ref(py)
    }

    /// returns the index of this SubstanceGroup in the owning molecule's list.
    #[pyo3(name = "GetIndexInMol")]
    fn get_index_in_mol(&self, py: Python<'_>) -> PyResult<usize> {
        let groups = self.mol().bind(py).borrow().mol.substance_groups().len();
        if groups == 0 {
            return Err(to_py(RDKitError::invariant_at(
                "No SubstanceGroups found on owning molecule",
                "!sgroups.empty()",
                "Code/GraphMol/SubstanceGroup.cpp",
                48,
            )));
        }
        match &self.handle {
            Handle::InMol { idx, .. } if *idx < groups => Ok(*idx),
            _ => Err(PyRuntimeError::new_err(
                "Unable to find own index in owning mol SubstanceGroup collection\n",
            )),
        }
    }

    /// returns a list of the indices of the atoms in this SubstanceGroup
    #[pyo3(name = "GetAtoms")]
    fn get_atoms(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        PyVectJ::new_py(py, self.read(py, |_, g| g.atoms().to_vec())?)
    }

    /// returns a list of the indices of the parent atoms in this SubstanceGroup
    #[pyo3(name = "GetParentAtoms")]
    fn get_parent_atoms(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        PyVectJ::new_py(py, self.read(py, |_, g| g.parent_atoms().to_vec())?)
    }

    /// returns a list of the indices of the bonds in this SubstanceGroup
    #[pyo3(name = "GetBonds")]
    fn get_bonds(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        PyVectJ::new_py(py, self.read(py, |_, g| g.bonds().to_vec())?)
    }

    /// Set the list of the indices of the atoms in this SubstanceGroup.
    /// Note that this does not update properties, CStates or Attachment Points.
    #[pyo3(name = "SetAtoms")]
    fn set_atoms(&mut self, py: Python<'_>, iterable: &Bound<'_, PyAny>) -> PyResult<()> {
        let atoms = index_vec(iterable)?;
        self.edit(py, |mol, g| g.set_atoms(mol, atoms))
    }

    /// Set the list of the indices of the parent atoms in this SubstanceGroup.
    /// Note that this does not update properties, CStates or Attachment Points.
    #[pyo3(name = "SetParentAtoms")]
    fn set_parent_atoms(&mut self, py: Python<'_>, iterable: &Bound<'_, PyAny>) -> PyResult<()> {
        let atoms = index_vec(iterable)?;
        self.edit(py, |mol, g| g.set_parent_atoms(mol, atoms))
    }

    /// Set the list of the indices of the bonds in this SubstanceGroup.
    /// Note that this does not update properties, CStates or Attachment Points.
    #[pyo3(name = "SetBonds")]
    fn set_bonds(&mut self, py: Python<'_>, iterable: &Bound<'_, PyAny>) -> PyResult<()> {
        let bonds = index_vec(iterable)?;
        self.edit(py, |mol, g| g.set_bonds(mol, bonds))
    }

    #[pyo3(name = "AddAtomWithIdx")]
    fn add_atom_with_idx(&mut self, py: Python<'_>, idx: u32) -> PyResult<()> {
        self.edit(py, |mol, g| g.add_atom_with_idx(mol, idx))
    }

    #[pyo3(name = "AddBondWithIdx")]
    fn add_bond_with_idx(&mut self, py: Python<'_>, idx: u32) -> PyResult<()> {
        self.edit(py, |mol, g| g.add_bond_with_idx(mol, idx))
    }

    #[pyo3(name = "AddParentAtomWithIdx")]
    fn add_parent_atom_with_idx(&mut self, py: Python<'_>, idx: u32) -> PyResult<()> {
        self.edit(py, |mol, g| g.add_parent_atom_with_idx(mol, idx))
    }

    #[pyo3(name = "AddAtomWithBookmark")]
    fn add_atom_with_bookmark(&mut self, py: Python<'_>, mark: i32) -> PyResult<()> {
        self.edit(py, |mol, g| g.add_atom_with_bookmark(mol, mark))
    }

    #[pyo3(name = "AddParentAtomWithBookmark")]
    fn add_parent_atom_with_bookmark(&mut self, py: Python<'_>, mark: i32) -> PyResult<()> {
        self.edit(py, |mol, g| g.add_parent_atom_with_bookmark(mol, mark))
    }

    #[pyo3(name = "AddBondWithBookmark")]
    fn add_bond_with_bookmark(&mut self, py: Python<'_>, mark: i32) -> PyResult<()> {
        self.edit(py, |mol, g| g.add_bond_with_bookmark(mol, mark))
    }

    #[pyo3(name = "AddCState")]
    #[allow(non_snake_case)]
    fn add_cstate(
        &mut self,
        py: Python<'_>,
        bondIdx: u32,
        vector: PyRef<'_, PyPoint3D>,
    ) -> PyResult<()> {
        let vector = vector.p;
        self.edit(py, |mol, g| g.add_cstate(mol, bondIdx, vector))
    }

    #[pyo3(name = "GetCStates")]
    fn get_cstates<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let states = self.read(py, |_, g| g.cstates().to_vec())?;
        let items = states
            .into_iter()
            .map(|c| Py::new(py, PySubstanceGroupCState { state: c }))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, items)
    }

    #[pyo3(name = "AddAttachPoint")]
    #[allow(non_snake_case)]
    fn add_attach_point(
        &mut self,
        py: Python<'_>,
        aIdx: u32,
        lvIdx: i32,
        idStr: &str,
    ) -> PyResult<()> {
        self.edit(py, |_, g| {
            g.add_attach_point(aIdx, lvIdx, idStr);
            Ok(())
        })
    }

    #[pyo3(name = "GetAttachPoints")]
    fn get_attach_points<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let points = self.read(py, |_, g| g.attach_points().to_vec())?;
        let items = points
            .into_iter()
            .map(|p| Py::new(py, PySubstanceGroupAttach { point: p }))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, items)
    }

    #[pyo3(name = "AddBracket")]
    fn add_bracket(&mut self, py: Python<'_>, pts: &Bound<'_, PyAny>) -> PyResult<()> {
        let n = pts.len()?;
        if n != 2 && n != 3 {
            return Err(PyValueError::new_err("pts object have a length of 2 or 3"));
        }
        let mut bracket: Bracket = [Point3D::new(0.0, 0.0, 0.0); 3];
        for (slot, item) in bracket.iter_mut().zip(pts.try_iter()?) {
            *slot = point(&item?)?;
        }
        self.edit(py, |_, g| {
            g.add_bracket(bracket);
            Ok(())
        })
    }

    #[pyo3(name = "GetBrackets")]
    fn get_brackets<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let brackets = self.read(py, |_, g| g.brackets().to_vec())?;
        let items = brackets
            .into_iter()
            .map(|b| {
                let pts = b
                    .into_iter()
                    .map(|p| Py::new(py, PyPoint3D { p }))
                    .collect::<PyResult<Vec<_>>>()?;
                PyTuple::new(py, pts)
            })
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, items)
    }

    #[pyo3(name = "ClearBrackets")]
    fn clear_brackets(&mut self, py: Python<'_>) -> PyResult<()> {
        self.edit(py, |_, g| {
            g.clear_brackets();
            Ok(())
        })
    }

    #[pyo3(name = "ClearCStates")]
    fn clear_cstates(&mut self, py: Python<'_>) -> PyResult<()> {
        self.edit(py, |_, g| {
            g.clear_cstates();
            Ok(())
        })
    }

    #[pyo3(name = "ClearAttachPoints")]
    fn clear_attach_points(&mut self, py: Python<'_>) -> PyResult<()> {
        self.edit(py, |_, g| {
            g.clear_attach_points();
            Ok(())
        })
    }

    /// returns the value of a particular property
    #[pyo3(name = "GetUnsignedVectProp")]
    fn get_unsigned_vect_prop(&self, py: Python<'_>, key: &str) -> PyResult<Py<PyAny>> {
        let v = self.with_dict(py, |d| match d.get(key) {
            None => Err(PyKeyError::new_err(key.to_string())),
            Some(RDValue::VecUnsignedInt(v)) => Ok(v.clone()),
            Some(_) => Err(bad_vect_cast(key)),
        })??;
        PyVectJ::new_py(py, v)
    }

    /// returns the value of a particular property
    #[pyo3(name = "GetStringVectProp")]
    fn get_string_vect_prop(&self, py: Python<'_>, key: &str) -> PyResult<Py<PyAny>> {
        let v = self.with_dict(py, |d| match d.get(key) {
            None => Err(PyKeyError::new_err(key.to_string())),
            Some(RDValue::VecString(v)) => Ok(v.clone()),
            Some(_) => Err(bad_vect_cast(key)),
        })??;
        PyVectStr::new_py(py, v)
    }
}

fn bad_vect_cast(key: &str) -> PyErr {
    PyValueError::new_err(format!(
        "key `{key}` exists but does not result in unregistered C++ type reason: bad any cast"
    ))
}

/// CSTATE for a SubstanceGroup
#[pyclass(
    name = "SubstanceGroupCState",
    module = "rdoxide.Chem.rdchem",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
pub struct PySubstanceGroupCState {
    state: CState,
}

#[pymethods]
impl PySubstanceGroupCState {
    #[new]
    fn py_new() -> Self {
        Self {
            state: CState {
                bond_idx: 0,
                vector: Point3D::new(0.0, 0.0, 0.0),
            },
        }
    }

    #[getter(bondIdx)]
    fn bond_idx(&self) -> u32 {
        self.state.bond_idx
    }

    #[getter]
    fn vector(&self) -> PyPoint3D {
        PyPoint3D {
            p: self.state.vector,
        }
    }
}

/// AttachPoint for a SubstanceGroup
#[pyclass(
    name = "SubstanceGroupAttach",
    module = "rdoxide.Chem.rdchem",
    skip_from_py_object
)]
#[derive(Debug, Clone)]
pub struct PySubstanceGroupAttach {
    point: AttachPoint,
}

#[pymethods]
impl PySubstanceGroupAttach {
    #[new]
    fn py_new() -> Self {
        Self {
            point: AttachPoint {
                a_idx: 0,
                lv_idx: 0,
                id: String::new(),
            },
        }
    }

    /// attachment index
    #[getter(aIdx)]
    fn a_idx(&self) -> u32 {
        self.point.a_idx
    }

    /// leaving atom or index (0 for implied)
    #[getter(lvIdx)]
    fn lv_idx(&self) -> i32 {
        self.point.lv_idx
    }

    /// attachment id
    #[getter]
    fn id(&self) -> String {
        self.point.id.clone()
    }
}

/// The result of `GetMolSubstanceGroups`: a mutable list of group copies.
#[pyclass(name = "SubstanceGroup_VECT", module = "rdoxide.Chem.rdchem")]
#[derive(Debug, Default)]
pub struct PySubstanceGroupVect {
    items: Vec<Py<PySubstanceGroup>>,
}

impl PySubstanceGroupVect {
    fn index(&self, i: isize) -> PyResult<usize> {
        let n = self.items.len() as isize;
        let j = if i < 0 { i + n } else { i };
        if (0..n).contains(&j) {
            Ok(j as usize)
        } else {
            Err(PyIndexError::new_err("Index out of range"))
        }
    }
}

#[pymethods]
impl PySubstanceGroupVect {
    #[new]
    fn py_new() -> Self {
        Self::default()
    }

    fn __len__(&self) -> usize {
        self.items.len()
    }

    fn __getitem__(&self, py: Python<'_>, idx: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        if let Ok(slice) = idx.cast::<PySlice>() {
            let ind = slice.indices(self.items.len() as isize)?;
            let mut items = Vec::new();
            let mut i = ind.start;
            while (ind.step > 0 && i < ind.stop) || (ind.step < 0 && i > ind.stop) {
                items.push(self.items[i as usize].clone_ref(py));
                i += ind.step;
            }
            return Ok(Py::new(py, Self { items })?.into_any());
        }
        let i = self.index(idx.extract()?)?;
        Ok(self.items[i].clone_ref(py).into_any())
    }

    fn __setitem__(&mut self, idx: isize, value: Py<PySubstanceGroup>) -> PyResult<()> {
        let i = self.index(idx)?;
        self.items[i] = value;
        Ok(())
    }

    fn __delitem__(&mut self, idx: isize) -> PyResult<()> {
        let i = self.index(idx)?;
        self.items.remove(i);
        Ok(())
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyList::new(py, self.items.iter().map(|g| g.clone_ref(py)))?
            .into_any()
            .try_iter()
            .map(Bound::into_any)
    }

    fn append(&mut self, item: Py<PySubstanceGroup>) {
        self.items.push(item);
    }

    fn extend(&mut self, items: &Bound<'_, PyAny>) -> PyResult<()> {
        for item in items.try_iter()? {
            self.items.push(item?.extract()?);
        }
        Ok(())
    }
}

fn group_count(mol: &Bound<'_, PyMol>) -> usize {
    mol.borrow().mol.substance_groups().len()
}

/// Adds `group` to the molecule and returns a reference to it.
fn add_group(mol: &Bound<'_, PyMol>, group: SubstanceGroup) -> PySubstanceGroup {
    add_substance_group(&mut mol.borrow_mut().mol, group);
    PySubstanceGroup::in_mol(mol.clone().unbind(), group_count(mol) - 1)
}

/// returns a copy of the molecule's SubstanceGroups (if any)
#[pyfunction]
#[pyo3(name = "GetMolSubstanceGroups")]
fn get_mol_substance_groups(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
) -> PyResult<PySubstanceGroupVect> {
    let groups = mol.borrow().mol.substance_groups().to_vec();
    let items = groups
        .into_iter()
        .map(|g| Py::new(py, PySubstanceGroup::copy_of(mol.clone().unbind(), g)))
        .collect::<PyResult<_>>()?;
    Ok(PySubstanceGroupVect { items })
}

/// returns a particular SubstanceGroup from the molecule
#[pyfunction]
#[pyo3(name = "GetMolSubstanceGroupWithIdx")]
fn get_mol_substance_group_with_idx(
    mol: &Bound<'_, PyMol>,
    idx: u32,
) -> PyResult<PySubstanceGroup> {
    if idx as usize >= group_count(mol) {
        return Err(PyIndexError::new_err(idx.to_string()));
    }
    Ok(PySubstanceGroup::in_mol(mol.clone().unbind(), idx as usize))
}

/// removes all SubstanceGroups from a molecule (if any)
#[pyfunction]
#[pyo3(name = "ClearMolSubstanceGroups")]
fn clear_mol_substance_groups(mol: &Bound<'_, PyMol>) {
    mol.borrow_mut().mol.substance_groups_mut().clear();
}

/// creates a new SubstanceGroup associated with a molecule, returns the new SubstanceGroup
#[pyfunction]
#[pyo3(name = "CreateMolSubstanceGroup")]
fn create_mol_substance_group(mol: &Bound<'_, PyMol>, r#type: &str) -> PySubstanceGroup {
    add_group(mol, SubstanceGroup::new(r#type))
}

/// creates a new DATA SubstanceGroup associated with a molecule, returns the new SubstanceGroup
#[pyfunction]
#[pyo3(name = "CreateMolDataSubstanceGroup")]
#[allow(non_snake_case)]
fn create_mol_data_substance_group(
    mol: &Bound<'_, PyMol>,
    fieldName: &str,
    value: &str,
) -> PyResult<PySubstanceGroup> {
    let mut group = SubstanceGroup::new("DAT");
    group.set_prop("FIELDNAME", fieldName, false).py()?;
    group
        .set_prop("DATAFIELDS", vec![value.to_string()], false)
        .py()?;
    Ok(add_group(mol, group))
}

/// adds a copy of a SubstanceGroup to a molecule, returns the new SubstanceGroup
#[pyfunction]
#[pyo3(name = "AddMolSubstanceGroup")]
fn add_mol_substance_group(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    sgroup: &Bound<'_, PySubstanceGroup>,
) -> PyResult<PySubstanceGroup> {
    let group = sgroup.borrow().read(py, |_, g| g.clone())?;
    Ok(add_group(mol, group))
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySubstanceGroup>()?;
    m.add_class::<PySubstanceGroupCState>()?;
    m.add_class::<PySubstanceGroupAttach>()?;
    m.add_class::<PySubstanceGroupVect>()?;
    m.add_function(wrap_pyfunction!(get_mol_substance_groups, m)?)?;
    m.add_function(wrap_pyfunction!(get_mol_substance_group_with_idx, m)?)?;
    m.add_function(wrap_pyfunction!(clear_mol_substance_groups, m)?)?;
    m.add_function(wrap_pyfunction!(create_mol_substance_group, m)?)?;
    m.add_function(wrap_pyfunction!(create_mol_data_substance_group, m)?)?;
    m.add_function(wrap_pyfunction!(add_mol_substance_group, m)?)?;
    Ok(())
}
