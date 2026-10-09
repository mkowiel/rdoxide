//! `rdchem.AtomMonomerInfo` / `AtomPDBResidueInfo` — port of upstream
//! `Code/GraphMol/Wrap/MonomerInfo.cpp`.
//!
//! Upstream hands out an internal reference to the atom's info object, so
//! changes made through `atom.GetPDBResidueInfo()` are visible on the atom.
//! [`InfoHandle::OnAtom`] reproduces that; free-standing infos own their data.

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use rdoxide::graphmol::monomer_info::{AtomMonomerInfo, AtomMonomerType, PdbResidueData};

use crate::atom::PyAtom;
use crate::enums::py_enum;

py_enum!(PyAtomMonomerType, AtomMonomerType, "AtomMonomerType");

#[derive(Debug)]
pub enum InfoHandle {
    Owned(AtomMonomerInfo),
    OnAtom(Py<PyAtom>),
}

/// The class to store monomer information attached to Atoms
#[pyclass(name = "AtomMonomerInfo", module = "rdkit.Chem.rdchem", subclass)]
#[derive(Debug)]
pub struct PyMonomerInfo {
    handle: InfoHandle,
}

/// The class to store PDB residue information attached to Atoms
#[pyclass(name = "AtomPDBResidueInfo", module = "rdkit.Chem.rdchem", extends = PyMonomerInfo)]
pub struct PyPdbResidueInfo;

impl PyMonomerInfo {
    /// The Python object for `info` (the most derived class, as upstream).
    fn make(py: Python<'_>, handle: InfoHandle, is_pdb: bool) -> PyResult<Py<PyMonomerInfo>> {
        let init = PyClassInitializer::from(PyMonomerInfo { handle });
        if is_pdb {
            let obj = Py::new(py, init.add_subclass(PyPdbResidueInfo))?;
            Ok(obj
                .into_bound(py)
                .into_any()
                .cast_into::<PyMonomerInfo>()?
                .unbind())
        } else {
            Py::new(py, init)
        }
    }

    /// A view on the monomer info of `atom`, if it has one.
    pub fn view(py: Python<'_>, atom: &Py<PyAtom>) -> PyResult<Option<Py<PyMonomerInfo>>> {
        let is_pdb = atom.bind(py).borrow().with(py, |a| {
            a.monomer_info().map(AtomMonomerInfo::is_pdb_residue)
        })?;
        is_pdb
            .map(|p| Self::make(py, InfoHandle::OnAtom(atom.clone_ref(py)), p))
            .transpose()
    }

    /// A copy of the data (`AtomMonomerInfo::copy`).
    pub fn snapshot(&self, py: Python<'_>) -> PyResult<AtomMonomerInfo> {
        self.read(py, AtomMonomerInfo::clone)
    }

    fn read<R>(&self, py: Python<'_>, f: impl FnOnce(&AtomMonomerInfo) -> R) -> PyResult<R> {
        match &self.handle {
            InfoHandle::Owned(i) => Ok(f(i)),
            InfoHandle::OnAtom(atom) => atom
                .bind(py)
                .borrow()
                .with(py, |a| a.monomer_info().map(f))?
                .ok_or_else(Self::detached),
        }
    }

    fn write<R>(
        &mut self,
        py: Python<'_>,
        f: impl FnOnce(&mut AtomMonomerInfo) -> R,
    ) -> PyResult<R> {
        match &mut self.handle {
            InfoHandle::Owned(i) => Ok(f(i)),
            InfoHandle::OnAtom(atom) => atom
                .bind(py)
                .borrow_mut()
                .with_mut(py, |a| a.monomer_info_mut().map(f))?
                .ok_or_else(Self::detached),
        }
    }

    fn detached() -> PyErr {
        PyRuntimeError::new_err("the atom no longer has monomer info")
    }
}

/// `Get*`/`Set*` pairs for the fields of the base class.
macro_rules! base_fields {
    ($(($get:ident, $pyget:literal, $set:ident, $pyset:literal, $ty:ty, $field:ident)),* $(,)?) => {
        #[pymethods]
        impl PyMonomerInfo {
            $(
                #[pyo3(name = $pyget)]
                fn $get(&self, py: Python<'_>) -> PyResult<$ty> {
                    self.read(py, |i| i.$field.clone())
                }
                #[pyo3(name = $pyset)]
                fn $set(&mut self, py: Python<'_>, val: $ty) -> PyResult<()> {
                    self.write(py, |i| i.$field = val)
                }
            )*
        }
    };
}

/// `Get*`/`Set*` pairs for the fields of [`PdbResidueData`].
macro_rules! pdb_fields {
    ($(($get:ident, $pyget:literal, $set:ident, $pyset:literal, $ty:ty, $field:ident)),* $(,)?) => {
        #[pymethods]
        impl PyPdbResidueInfo {
            $(
                #[pyo3(name = $pyget)]
                fn $get(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<$ty> {
                    slf.as_super()
                        .read(py, |i| i.pdb.clone().unwrap_or_default().$field)
                }
                #[pyo3(name = $pyset)]
                fn $set(mut slf: PyRefMut<'_, Self>, py: Python<'_>, val: $ty) -> PyResult<()> {
                    slf.as_super().write(py, |i| i.pdb_mut().$field = val)
                }
            )*
        }
    };
}

base_fields! {
    (get_name, "GetName", set_name, "SetName", String, name),
    (get_residue_name, "GetResidueName", set_residue_name, "SetResidueName", String, residue_name),
    (get_residue_number, "GetResidueNumber", set_residue_number, "SetResidueNumber", i32, residue_number),
    (get_chain_id, "GetChainId", set_chain_id, "SetChainId", String, chain_id),
    (get_monomer_class, "GetMonomerClass", set_monomer_class, "SetMonomerClass", String, monomer_class),
}

pdb_fields! {
    (get_serial_number, "GetSerialNumber", set_serial_number, "SetSerialNumber", u32, serial_number),
    (get_alt_loc, "GetAltLoc", set_alt_loc, "SetAltLoc", String, alt_loc),
    (get_insertion_code, "GetInsertionCode", set_insertion_code, "SetInsertionCode", String, insertion_code),
    (get_occupancy, "GetOccupancy", set_occupancy, "SetOccupancy", f64, occupancy),
    (get_temp_factor, "GetTempFactor", set_temp_factor, "SetTempFactor", f64, temp_factor),
    (get_is_hetero_atom, "GetIsHeteroAtom", set_is_hetero_atom, "SetIsHeteroAtom", bool, is_hetero_atom),
    (get_secondary_structure, "GetSecondaryStructure", set_secondary_structure, "SetSecondaryStructure", u32, secondary_structure),
    (get_segment_number, "GetSegmentNumber", set_segment_number, "SetSegmentNumber", u32, segment_number),
}

#[pymethods]
impl PyMonomerInfo {
    #[new]
    #[pyo3(signature = (r#type=None, name="", residueName="", resNum=0, chainId="", monomerClass=""))]
    #[allow(non_snake_case)]
    fn new(
        r#type: Option<&Bound<'_, PyAny>>,
        name: &str,
        residueName: &str,
        resNum: i32,
        chainId: &str,
        monomerClass: &str,
    ) -> PyResult<Self> {
        let monomer_type = r#type
            .map(PyAtomMonomerType::extract)
            .transpose()?
            .unwrap_or_default();
        Ok(PyMonomerInfo {
            handle: InfoHandle::Owned(AtomMonomerInfo {
                monomer_type,
                name: name.to_string(),
                residue_name: residueName.to_string(),
                residue_number: resNum,
                chain_id: chainId.to_string(),
                monomer_class: monomerClass.to_string(),
                pdb: None,
            }),
        })
    }

    #[pyo3(name = "GetMonomerType")]
    fn get_monomer_type(&self, py: Python<'_>) -> PyResult<PyAtomMonomerType> {
        self.read(py, |i| PyAtomMonomerType(i.monomer_type))
    }

    #[pyo3(name = "SetMonomerType")]
    fn set_monomer_type(&mut self, py: Python<'_>, typ: &Bound<'_, PyAny>) -> PyResult<()> {
        let t = PyAtomMonomerType::extract(typ)?;
        self.write(py, |i| i.monomer_type = t)
    }
}

#[pymethods]
impl PyPdbResidueInfo {
    #[new]
    #[pyo3(signature = (atomName=None, serialNumber=1, altLoc="", residueName="", residueNumber=0,
        chainId="", insertionCode="", occupancy=1.0, tempFactor=0.0, isHeteroAtom=false,
        secondaryStructure=0, segmentNumber=0, monomerClass=""))]
    #[allow(non_snake_case, clippy::too_many_arguments)]
    fn new(
        atomName: Option<&str>,
        serialNumber: i32,
        altLoc: &str,
        residueName: &str,
        residueNumber: i32,
        chainId: &str,
        insertionCode: &str,
        occupancy: f64,
        tempFactor: f64,
        isHeteroAtom: bool,
        secondaryStructure: u32,
        segmentNumber: u32,
        monomerClass: &str,
    ) -> PyClassInitializer<Self> {
        // the no-argument constructor leaves every field at its default
        let info = match atomName {
            None => AtomMonomerInfo::new_pdb(""),
            Some(name) => AtomMonomerInfo {
                monomer_type: AtomMonomerType::PdbResidue,
                name: name.to_string(),
                residue_name: residueName.to_string(),
                residue_number: residueNumber,
                chain_id: chainId.to_string(),
                monomer_class: monomerClass.to_string(),
                pdb: Some(PdbResidueData {
                    serial_number: serialNumber as u32,
                    alt_loc: altLoc.to_string(),
                    insertion_code: insertionCode.to_string(),
                    occupancy,
                    temp_factor: tempFactor,
                    is_hetero_atom: isHeteroAtom,
                    secondary_structure: secondaryStructure,
                    segment_number: segmentNumber,
                }),
            },
        };
        PyClassInitializer::from(PyMonomerInfo {
            handle: InfoHandle::Owned(info),
        })
        .add_subclass(PyPdbResidueInfo)
    }
}

/// `Atom.GetMonomerInfo` etc.: shared by the `Atom` bindings.
pub fn get_monomer_info(py: Python<'_>, atom: &Py<PyAtom>) -> PyResult<Option<Py<PyMonomerInfo>>> {
    PyMonomerInfo::view(py, atom)
}

/// `GetPDBResidueInfo`: like [`get_monomer_info`] but rejects other types.
pub fn get_pdb_residue_info(
    py: Python<'_>,
    atom: &Py<PyAtom>,
) -> PyResult<Option<Py<PyMonomerInfo>>> {
    let is_pdb = atom
        .bind(py)
        .borrow()
        .with(py, |a| a.monomer_info().map(|i| i.monomer_type))?;
    match is_pdb {
        Some(t) if t != AtomMonomerType::PdbResidue => Err(not_pdb()),
        _ => PyMonomerInfo::view(py, atom),
    }
}

fn not_pdb() -> PyErr {
    PyValueError::new_err("MonomerInfo is not a PDB Residue")
}

/// `SetMonomerInfo` (`pdb_only` = `SetPDBResidueInfo`).
pub fn set_monomer_info(
    py: Python<'_>,
    atom: &mut PyAtom,
    info: Option<&Bound<'_, PyMonomerInfo>>,
    pdb_only: bool,
) -> PyResult<()> {
    let snapshot = info.map(|i| i.borrow().snapshot(py)).transpose()?;
    if let Some(i) = &snapshot {
        if pdb_only && i.monomer_type != AtomMonomerType::PdbResidue {
            return Err(not_pdb());
        }
    }
    let snapshot = snapshot.map(|mut i| {
        // a base-class object typed PDBRESIDUE is treated as one
        if i.monomer_type == AtomMonomerType::PdbResidue {
            i.pdb_mut();
        }
        i
    });
    atom.with_mut(py, |a| a.set_monomer_info(snapshot))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMonomerInfo>()?;
    m.add_class::<PyPdbResidueInfo>()?;
    Ok(())
}
