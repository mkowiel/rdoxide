//! `rdkit.Chem.rdchem.StereoInfo` with its enums and
//! `rdkit.Chem.rdmolops.FindPotentialStereo` (upstream
//! `GraphMol/Wrap/Chirality.cpp` and `ChiralityOps.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rdoxide::graphmol::chirality::find_stereo::{
    NOATOM, StereoDescriptor, StereoInfo, StereoSpecified, StereoType, find_potential_stereo,
};

use crate::enums::{add_enum, to_py};
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

const TYPE_ENUM: &str = "StereoType";
const SPECIFIED_ENUM: &str = "StereoSpecified";
const DESCRIPTOR_ENUM: &str = "StereoDescriptor";

const TYPES: [(&str, StereoType); 8] = [
    ("Unspecified", StereoType::Unspecified),
    ("Atom_Tetrahedral", StereoType::AtomTetrahedral),
    ("Atom_SquarePlanar", StereoType::AtomSquarePlanar),
    (
        "Atom_TrigonalBipyramidal",
        StereoType::AtomTrigonalBipyramidal,
    ),
    ("Atom_Octahedral", StereoType::AtomOctahedral),
    ("Bond_Double", StereoType::BondDouble),
    ("Bond_Cumulene_Even", StereoType::BondCumuleneEven),
    ("Bond_Atropisomer", StereoType::BondAtropisomer),
];
const SPECIFIED: [(&str, StereoSpecified); 3] = [
    ("Unspecified", StereoSpecified::Unspecified),
    ("Specified", StereoSpecified::Specified),
    ("Unknown", StereoSpecified::Unknown),
];
// "None" cannot be used as a Python name, hence "NoValue".
const DESCRIPTORS: [(&str, StereoDescriptor); 5] = [
    ("NoValue", StereoDescriptor::None),
    ("Tet_CW", StereoDescriptor::TetCw),
    ("Tet_CCW", StereoDescriptor::TetCcw),
    ("Bond_Cis", StereoDescriptor::BondCis),
    ("Bond_Trans", StereoDescriptor::BondTrans),
];

fn members<T: Copy>(
    table: &[(&'static str, T)],
    code: impl Fn(T) -> i64,
) -> Vec<(&'static str, i64)> {
    table.iter().map(|&(n, v)| (n, code(v))).collect()
}

/// Look `value` up in `table` by its integer code.
fn from_code<T: Copy>(
    table: &[(&'static str, T)],
    code: impl Fn(T) -> i64,
    value: i64,
    what: &str,
) -> PyResult<T> {
    table
        .iter()
        .map(|&(_, v)| v)
        .find(|&v| code(v) == value)
        .ok_or_else(|| PyValueError::new_err(format!("invalid {what} value: {value}")))
}

fn type_code(v: StereoType) -> i64 {
    v as i64
}
fn specified_code(v: StereoSpecified) -> i64 {
    v as i64
}
fn descriptor_code(v: StereoDescriptor) -> i64 {
    v as i64
}

/// Class describing stereochemistry
#[pyclass(name = "StereoInfo", module = "rdkit.Chem.rdchem", skip_from_py_object)]
#[derive(Clone, Default)]
pub struct PyStereoInfo {
    info: StereoInfo,
}

#[pymethods]
impl PyStereoInfo {
    #[classattr]
    #[pyo3(name = "NOATOM")]
    fn noatom() -> u32 {
        NOATOM
    }

    /// the type of stereo
    #[getter(r#type)]
    fn get_type(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(py, TYPE_ENUM, type_code(self.info.stereo_type))
    }
    #[setter(r#type)]
    fn set_type(&mut self, value: i64) -> PyResult<()> {
        self.info.stereo_type = from_code(&TYPES, type_code, value, TYPE_ENUM)?;
        Ok(())
    }

    /// whether or not it is specified
    #[getter]
    fn specified(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(py, SPECIFIED_ENUM, specified_code(self.info.specified))
    }
    #[setter(specified)]
    fn set_specified(&mut self, value: i64) -> PyResult<()> {
        self.info.specified = from_code(&SPECIFIED, specified_code, value, SPECIFIED_ENUM)?;
        Ok(())
    }

    /// index of the item the stereo concerns
    #[getter(centeredOn)]
    fn centered_on(&self) -> u32 {
        self.info.centered_on
    }
    #[setter(centeredOn)]
    fn set_centered_on(&mut self, value: u32) {
        self.info.centered_on = value;
    }

    /// stereo descriptor
    #[getter]
    fn descriptor(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(py, DESCRIPTOR_ENUM, descriptor_code(self.info.descriptor))
    }
    #[setter(descriptor)]
    fn set_descriptor(&mut self, value: i64) -> PyResult<()> {
        self.info.descriptor = from_code(&DESCRIPTORS, descriptor_code, value, DESCRIPTOR_ENUM)?;
        Ok(())
    }

    /// permutation index (used for non-tetrahedral chirality)
    #[getter]
    fn permutation(&self) -> u32 {
        self.info.permutation
    }
    #[setter(permutation)]
    fn set_permutation(&mut self, value: u32) {
        self.info.permutation = value;
    }

    /// indices of the atoms controlling the stereo
    #[getter(controllingAtoms)]
    fn controlling_atoms(&self) -> Vec<u32> {
        self.info.controlling_atoms.clone()
    }
}

/// find potential stereo elements in a molecule and returns them as StereoInfo objects
/// Note that this function is still somewhat experimental and the API
/// and results may change in a future release.
#[pyfunction]
#[pyo3(name = "FindPotentialStereo", signature = (mol, cleanIt = false, flagPossible = true))]
#[allow(non_snake_case)]
fn py_find_potential_stereo(
    mut mol: PyRefMut<'_, PyMol>,
    cleanIt: bool,
    flagPossible: bool,
) -> PyResult<Vec<PyStereoInfo>> {
    let infos = find_potential_stereo(&mut mol.mol, cleanIt, flagPossible).py()?;
    Ok(infos
        .into_iter()
        .map(|info| PyStereoInfo { info })
        .collect())
}

/// Registers the enums and `StereoInfo` in `rdkit.Chem.rdchem`.
pub fn add_chem(m: &Bound<'_, PyModule>) -> PyResult<()> {
    add_enum(m, TYPE_ENUM, &members(&TYPES, type_code), None)?;
    add_enum(
        m,
        SPECIFIED_ENUM,
        &members(&SPECIFIED, specified_code),
        None,
    )?;
    add_enum(
        m,
        DESCRIPTOR_ENUM,
        &members(&DESCRIPTORS, descriptor_code),
        None,
    )?;
    m.add_class::<PyStereoInfo>()
}

/// Registers `FindPotentialStereo` in `rdkit.Chem.rdmolops`.
pub fn add_ops(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_find_potential_stereo, m)?)
}
