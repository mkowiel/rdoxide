//! Sub-molecule extraction of `rdoxide.Chem.rdmolops` (upstream
//! `GraphMol/Wrap/MolOps.cpp`): `SubsetMethod`, `SubsetOptions`,
//! `SubsetInfo`, `UIntUIntMap`, `CopyMolSubset` and `PathToSubmol`.

use std::collections::BTreeMap;

use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rdoxide::graphmol::subgraphs::path_to_submol;
use rdoxide::graphmol::subset::{
    SubsetInfo, SubsetMethod, SubsetOptions, copy_mol_subset, copy_mol_subset_atoms_bonds,
};

use crate::enums::{add_enum, to_py};
use crate::errors::{IntoPyResult, argument_error, boost_type_name, expect_dict};
use crate::mol::PyMol;
use crate::pyconv::{extract_unsigned, python_object_to_vect};
use crate::rdmolfiles_write::Args;

/// An `unsigned int -> unsigned int` map (Boost's `map_indexing_suite`).
#[pyclass(
    name = "UIntUIntMap",
    module = "rdoxide.Chem.rdmolops",
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
pub struct PyUIntUIntMap {
    map: BTreeMap<u32, u32>,
}

#[pymethods]
impl PyUIntUIntMap {
    fn __len__(&self) -> usize {
        self.map.len()
    }

    fn __getitem__(&self, key: u32) -> PyResult<u32> {
        self.map
            .get(&key)
            .copied()
            .ok_or_else(|| PyKeyError::new_err("Invalid key"))
    }

    fn __setitem__(&mut self, key: u32, value: u32) {
        self.map.insert(key, value);
    }

    /// Deleting a missing key is not an error (it erases from a `std::map`).
    fn __delitem__(&mut self, key: u32) {
        self.map.remove(&key);
    }

    fn __contains__(&self, key: u32) -> bool {
        self.map.contains_key(&key)
    }

    /// Iterates over the `(key, value)` pairs.
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let pairs = PyList::new(py, self.map.iter().map(|(k, v)| (*k, *v)))?;
        Ok(pairs.try_iter()?.into_any())
    }
}

fn method_to_py(py: Python<'_>, method: SubsetMethod) -> PyResult<Py<PyAny>> {
    to_py(py, "SubsetMethod", method as i64)
}

/// Options controlling what `CopyMolSubset` copies.
#[pyclass(
    name = "SubsetOptions",
    module = "rdoxide.Chem.rdmolops",
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
pub struct PySubsetOptions {
    options: SubsetOptions,
}

#[pymethods]
impl PySubsetOptions {
    #[new]
    fn py_new() -> Self {
        Self::default()
    }

    /// Sanitize the resulting subset
    #[getter]
    fn sanitize(&self) -> bool {
        self.options.sanitize
    }
    #[setter]
    fn set_sanitize(&mut self, v: bool) {
        self.options.sanitize = v;
    }

    /// clear all computed props on the subsetted molecule
    #[getter(clearComputedProps)]
    fn clear_computed_props(&self) -> bool {
        self.options.clear_computed_props
    }
    #[setter(clearComputedProps)]
    fn set_clear_computed_props(&mut self, v: bool) {
        self.options.clear_computed_props = v;
    }

    /// Return the subset as a query
    #[getter(copyAsQuery)]
    fn copy_as_query(&self) -> bool {
        self.options.copy_as_query
    }
    #[setter(copyAsQuery)]
    fn set_copy_as_query(&mut self, v: bool) {
        self.options.copy_as_query = v;
    }

    /// Copy the active coordinates from the molecule
    #[getter(copyCoordinates)]
    fn copy_coordinates(&self) -> bool {
        self.options.copy_coordinates
    }
    #[setter(copyCoordinates)]
    fn set_copy_coordinates(&mut self, v: bool) {
        self.options.copy_coordinates = v;
    }

    /// What conformer idx to use for the coordinates default is -1
    #[getter(conformerIdx)]
    fn conformer_idx(&self) -> u32 {
        self.options.conformer_idx
    }
    #[setter(conformerIdx)]
    fn set_conformer_idx(&mut self, v: &Bound<'_, PyAny>) -> PyResult<()> {
        self.options.conformer_idx = extract_unsigned(v)?;
        Ok(())
    }

    /// Subsetting method to use
    #[getter]
    fn method(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        method_to_py(py, self.options.method)
    }
    #[setter]
    fn set_method(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let py = value.py();
        let class = method_to_py(py, SubsetMethod::default())?
            .bind(py)
            .get_type();
        self.options.method = if value.is_instance(&class)? {
            match value.extract::<i64>()? {
                0 => SubsetMethod::BondsBetweenAtoms,
                _ => SubsetMethod::Bonds,
            }
        } else {
            return Err(argument_error(
                py,
                "None",
                &["SubsetOptions".to_string(), boost_type_name(value)],
                "None(RDKit::SubsetOptions {lvalue}, RDKit::SubsetMethod)",
            ));
        };
        Ok(())
    }
}

/// Where the atoms and bonds of a subset came from.
#[pyclass(name = "SubsetInfo", module = "rdoxide.Chem.rdmolops")]
#[derive(Debug)]
pub struct PySubsetInfo {
    atom_mapping: Py<PyUIntUIntMap>,
    bond_mapping: Py<PyUIntUIntMap>,
}

#[pymethods]
impl PySubsetInfo {
    #[new]
    fn py_new(py: Python<'_>) -> PyResult<Self> {
        Ok(PySubsetInfo {
            atom_mapping: Py::new(py, PyUIntUIntMap::default())?,
            bond_mapping: Py::new(py, PyUIntUIntMap::default())?,
        })
    }

    /// mapping from the original atom index to the subset atom index
    #[getter(atomMapping)]
    fn atom_mapping(&self, py: Python<'_>) -> Py<PyUIntUIntMap> {
        self.atom_mapping.clone_ref(py)
    }
    #[setter(atomMapping)]
    fn set_atom_mapping(&mut self, v: PyRef<'_, PyUIntUIntMap>) -> PyResult<()> {
        self.atom_mapping = Py::new(v.py(), v.clone())?;
        Ok(())
    }

    ///  mapping from the original bond index to the subset bond index
    #[getter(bondMapping)]
    fn bond_mapping(&self, py: Python<'_>) -> Py<PyUIntUIntMap> {
        self.bond_mapping.clone_ref(py)
    }
    #[setter(bondMapping)]
    fn set_bond_mapping(&mut self, v: PyRef<'_, PyUIntUIntMap>) -> PyResult<()> {
        self.bond_mapping = Py::new(v.py(), v.clone())?;
        Ok(())
    }
}

impl PySubsetInfo {
    fn record(&self, py: Python<'_>, info: SubsetInfo) {
        self.atom_mapping.bind(py).borrow_mut().map = info.atom_mapping;
        self.bond_mapping.bind(py).borrow_mut().map = info.bond_mapping;
    }
}

const COPY_SIGNATURES: &str = concat!(
    "CopyMolSubset(RDKit::ROMol mol, boost::python::api::object path, RDKit::SubsetInfo {lvalue} subsetInfo, RDKit::SubsetOptions options=<rdoxide.Chem.rdmolops.SubsetOptions object at 0x0>)\n",
    "    CopyMolSubset(RDKit::ROMol mol, boost::python::api::object path, RDKit::SubsetOptions options=<rdoxide.Chem.rdmolops.SubsetOptions object at 0x0>)\n",
    "    CopyMolSubset(RDKit::ROMol mol, boost::python::api::object atomIndices, boost::python::api::object bondIndices, RDKit::SubsetInfo {lvalue} subsetInfo, RDKit::SubsetOptions options=<rdoxide.Chem.rdmolops.SubsetOptions object at 0x0>)\n",
    "    CopyMolSubset(RDKit::ROMol mol, boost::python::api::object atomIndices, boost::python::api::object bondIndices, RDKit::SubsetOptions options=<rdoxide.Chem.rdmolops.SubsetOptions object at 0x0>)"
);

/// A sequence of indices; `None`, 0 and empty sequences give an empty list.
fn indices(obj: Option<&Bound<'_, PyAny>>) -> PyResult<Vec<u32>> {
    Ok(python_object_to_vect(obj, None)?.unwrap_or_default())
}

/// Extract a subgraph from an ROMol. Bonds, atoms, substance groups and
/// stereo groups are only extracted to the subgraph if all participant entities
/// are contained within the given atoms and bonds.
///
/// ARGUMENTS:
///  - mol - starting mol
///  - atoms - indices atoms to extract
///  - bonds - indices bonds to extract
///  - subsetInfo - optional subsetInfo to record the atoms and bonds used
///  - options - optional subset options, note the method is ignored since all the atoms and bonds are specified
///
/// or, with a path instead of atoms and bonds:
///  - path - the indices of atoms or bonds to extract. If an index falls
///           outside of the acceptable indices, it is ignored.
///           Use SubsetMethod.BONDS to indicate a bond path and BONDS_BETWEEN_ATOMS
///           to indicate an atom path with any bond that includes atoms in the path.
#[pyfunction]
#[pyo3(name = "CopyMolSubset", signature = (*args, **kwargs))]
fn py_copy_mol_subset(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyMol> {
    let a = Args {
        name: "CopyMolSubset",
        args,
        kwargs,
    };
    let bad_call = || {
        let got: Vec<String> = args.iter().map(|o| boost_type_name(&o)).collect();
        argument_error(py, "CopyMolSubset", &got, COPY_SIGNATURES)
    };
    let Some(mol) = a.get(0, "mol")?.and_then(|m| m.cast_into::<PyMol>().ok()) else {
        return Err(bad_call());
    };
    let is_info = |o: &Bound<'_, PyAny>| o.cast::<PySubsetInfo>().is_ok();
    let is_options = |o: &Bound<'_, PyAny>| o.cast::<PySubsetOptions>().is_ok();

    // atoms + bonds, or a path: told apart by what follows the first list
    let first = a.get(1, "atomIndices")?.or(a.get(1, "path")?);
    if first.is_none() {
        return Err(bad_call());
    }
    let second = a.get(2, "bondIndices")?;
    let atoms_and_bonds = second
        .as_ref()
        .is_some_and(|s| !is_info(s) && !is_options(s));
    let base = if atoms_and_bonds { 3 } else { 2 };
    let info_arg = a.get(base, "subsetInfo")?;
    let (info, options_arg) = match info_arg {
        Some(i) if is_info(&i) => (Some(i), a.get(base + 1, "options")?),
        Some(o) if is_options(&o) => (None, Some(o)),
        Some(_) => return Err(bad_call()),
        None => (None, a.get(base + 1, "options")?),
    };
    let options_arg = options_arg.or(a.get(base, "options")?.filter(is_options));
    let options = match &options_arg {
        Some(o) => o
            .cast::<PySubsetOptions>()
            .map_err(|_| bad_call())?
            .borrow()
            .options
            .clone(),
        None => SubsetOptions::default(),
    };

    let m = mol.borrow();
    let (sub, subset_info) = if atoms_and_bonds {
        let atoms = indices(first.as_ref())?;
        let bonds = indices(second.as_ref())?;
        copy_mol_subset_atoms_bonds(&m.mol, &atoms, &bonds, &options).py()?
    } else {
        let path = indices(first.as_ref())?;
        copy_mol_subset(&m.mol, &path, &options).py()?
    };
    if let Some(info) = info {
        info.cast::<PySubsetInfo>()?
            .borrow()
            .record(py, subset_info);
    }
    Ok(PyMol { mol: sub })
}

/// Extracts the molecule made of the bonds of a path (and their atoms).
///
///   ARGUMENTS:
///
///     - mol: the molecule
///     - path: a sequence of bond indices
///     - useQuery: (optional) if True, the result atoms and bonds are queries
///     - atomMap: (optional) a dict that is filled with old -> new atom indices
#[pyfunction]
#[pyo3(name = "PathToSubmol", signature = (mol, path, useQuery = false, atomMap = None))]
#[allow(non_snake_case)]
fn py_path_to_submol(
    mol: PyRef<'_, PyMol>,
    path: &Bound<'_, PyAny>,
    useQuery: bool,
    atomMap: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyMol> {
    let mut bonds = Vec::new();
    for i in 0..path.len()? {
        bonds.push(path.get_item(i)?.extract::<u32>()? as i32);
    }
    let (sub, mapping) = path_to_submol(&mol.mol, &bonds, useQuery).py()?;
    if let Some(map) = atomMap.filter(|m| !m.is_none()) {
        // make sure the optional argument actually is a dictionary
        let map = expect_dict(map)?;
        map.clear();
        for (old, new) in mapping {
            map.set_item(old, new)?;
        }
    }
    Ok(PyMol { mol: sub })
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    add_enum(
        m,
        "SubsetMethod",
        &[("BONDS_BETWEEN_ATOMS", 0), ("BONDS", 1)],
        None,
    )?;
    m.add_class::<PyUIntUIntMap>()?;
    m.add_class::<PySubsetOptions>()?;
    m.add_class::<PySubsetInfo>()?;
    m.add_function(wrap_pyfunction!(py_copy_mol_subset, m)?)?;
    m.add_function(wrap_pyfunction!(py_path_to_submol, m)?)?;
    Ok(())
}
