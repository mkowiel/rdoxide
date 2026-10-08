//! `rdkit.Chem.rdMolDescriptors.BCUT2D` (upstream `Descriptors/BCUT.cpp` and
//! its wrappers in `Wrap/rdMolDescriptors.cpp`).

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::{PyList, PyString, PyTuple};
use rdoxide::graphmol::descriptors::bcut::{
    BCUT2D_VERSION, bcut2d, bcut2d_with_prop_name, bcut2d_with_props,
};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Implements BCUT descriptors From J. Chem. Inf. Comput. Sci., Vol. 39, No. 1, 1999
/// Diagonal elements are (currently) atomic mass, gasteiger charge, crippen logP and crippen MR
///
/// BCUT2D(mol) returns [mass eigen value high, mass eigen value low,
///                      gasteiger charge eigenvalue high, gasteiger charge low,
///                      crippen lowgp  eigenvalue high, crippen lowgp  low,
///                      crippen mr eigenvalue high, crippen mr low]
///
/// BCUT2D(mol, atom_props) returns a 2D BCUT (eigen value hi, eigenvalue low) given the
/// molecule and the specified atom props: a list or tuple of floats equal in size to the
/// number of atoms in mol, or the name of a double atom property.
#[pyfunction]
#[pyo3(name = "BCUT2D", signature = (mol, atom_props = None))]
fn py_bcut2d<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyMol>,
    atom_props: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyAny>> {
    let m = &mol.borrow().mol;
    let pair = |(hi, lo): (f64, f64)| -> PyResult<Bound<'py, PyAny>> {
        Ok(PyTuple::new(py, [hi, lo])?.into_any())
    };
    match atom_props {
        None => Ok(PyList::new(py, bcut2d(m).py()?)?.into_any()),
        Some(p) if p.is_instance_of::<PyString>() => {
            pair(bcut2d_with_prop_name(m, &p.extract::<String>()?).py()?)
        }
        Some(p) if p.is_instance_of::<PyList>() || p.is_instance_of::<PyTuple>() => {
            let props: Vec<f64> = p.extract()?;
            pair(bcut2d_with_props(m, &props).py()?)
        }
        Some(_) => Err(PyTypeError::new_err(
            "atom_props must be a list, a tuple or the name of an atom property",
        )),
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("_BCUT2D_version", BCUT2D_VERSION)?;
    m.add_function(wrap_pyfunction!(py_bcut2d, m)?)
}
