//! `rdoxide.Chem.rdMolProcessing` (upstream `MolProcessing/Wrap/rdMolProcessing.cpp`).

use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rdoxide::graphmol::mol_processing::{SupplierOptions, get_fingerprints_for_mols_in_file};

use crate::datastructs::PyExplicitBitVect;
use crate::errors::IntoPyResult;
use crate::rdfingerprintgenerator::PyFPGen;

/// Supplier Options
#[pyclass(
    name = "SupplierOptions",
    module = "rdoxide.Chem.rdMolProcessing",
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PySupplierOptions {
    options: SupplierOptions,
}

macro_rules! option_attributes {
    ($($py:literal => $get:ident / $set:ident : $field:ident : $ty:ty,)*) => {
        #[pymethods]
        impl PySupplierOptions {
            #[new]
            fn py_new() -> Self {
                Self::default()
            }
            $(
                #[getter($py)]
                fn $get(&self) -> $ty {
                    self.options.$field.clone()
                }
                #[setter($py)]
                fn $set(&mut self, value: $ty) {
                    self.options.$field = value;
                }
            )*
        }
    };
}

option_attributes! {
    "numThreads" => num_threads / set_num_threads : num_writer_threads : i32,
    "sanitize" => sanitize / set_sanitize : sanitize : bool,
    "removeHs" => remove_hs / set_remove_hs : remove_hs : bool,
    "strictParsing" => strict_parsing / set_strict_parsing : strict_parsing : bool,
    "delimiter" => delimiter / set_delimiter : delimiter : String,
    "smilesColumn" => smiles_column / set_smiles_column : smiles_column : i32,
    "nameColumn" => name_column / set_name_column : name_column : i32,
    "titleLine" => title_line / set_title_line : title_line : bool,
    "nameRecord" => name_record / set_name_record : name_record : String,
    "confId2D" => conf_id_2d / set_conf_id_2d : conf_id_2d : i32,
    "confId3D" => conf_id_3d / set_conf_id_3d : conf_id_3d : i32,
}

/// returns the fingerprints for the molecules in a file
#[pyfunction]
#[pyo3(name = "GetFingerprintsForMolsInFile", signature = (filename, generator = None, options = None))]
fn get_fingerprints_for_mols_in_file_py<'py>(
    py: Python<'py>,
    filename: &str,
    generator: Option<&Bound<'py, PyFPGen>>,
    options: Option<&Bound<'py, PySupplierOptions>>,
) -> PyResult<Bound<'py, PyTuple>> {
    let options = options.map_or_else(SupplierOptions::default, |o| o.borrow().options.clone());
    let generator = generator.map(|g| g.borrow());
    let gen_ref = generator.as_ref().map(|g| &g.gen_);
    let fps = py
        .detach(|| get_fingerprints_for_mols_in_file(filename, &options, gen_ref))
        .py()?;
    let items = fps
        .into_iter()
        .map(|fp| match fp {
            Some(fp) => Py::new(py, PyExplicitBitVect::wrap(fp)).map(Py::into_any),
            None => Ok(py.None()),
        })
        .collect::<PyResult<Vec<_>>>()?;
    PyTuple::new(py, items)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions for working with groups of molecules",
    )?;
    m.add_class::<PySupplierOptions>()?;
    m.add_function(wrap_pyfunction!(get_fingerprints_for_mols_in_file_py, m)?)?;
    Ok(())
}
