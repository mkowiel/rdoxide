//! Python bindings for rdoxide, API-compatible with RDKit's Boost.Python
//! wrappers. Built as the extension module `rdkit._rdoxide`; each RDKit
//! extension module (`rdkit.rdBase`, `rdkit.Chem.rdchem`, ...) is a submodule
//! here whose `__name__` is the upstream module name. The pure-Python shims in
//! `python/rdkit` replace themselves in `sys.modules` with these submodules.

use pyo3::prelude::*;
use pyo3::types::PyDict;

mod atom;
mod bond;
mod conformer;
mod debug;
mod enums;
mod enums_impl;
mod errors;
mod geometry;
mod mol;
mod perception;
mod props;
mod rdbase;
mod rdchem;
mod rdmolfiles;
mod rdmolfiles_molblock;
mod rdmolfiles_molwrite;
mod rdmolfiles_sdsupplier;
mod rdmolfiles_smarts;
mod rdmolfiles_write;
mod rdmolops;
mod rdmolops_hs;
mod ringinfo;
mod stereo_group;
mod substruct;

/// Create a submodule named `full_name` (e.g. `rdkit.Chem.rdchem`) and attach
/// it to `parent` under `attr`.
fn add_submodule<'py>(
    parent: &Bound<'py, PyModule>,
    attr: &str,
    full_name: &str,
    init: impl FnOnce(&Bound<'py, PyModule>) -> PyResult<()>,
) -> PyResult<()> {
    let py = parent.py();
    let m = PyModule::new(py, attr)?;
    m.setattr("__name__", full_name)?;
    init(&m)?;
    parent.add(attr, &m)?;
    // make `import rdkit.x.y` find it even before the shim runs
    let modules = py
        .import("sys")?
        .getattr("modules")?
        .cast_into::<PyDict>()?;
    modules.set_item(full_name, &m)?;
    Ok(())
}

#[pymodule]
fn _rdoxide(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", rdoxide::rdgeneral::RDOXIDE_VERSION)?;
    add_submodule(m, "rdBase", "rdkit.rdBase", rdbase::init)?;
    add_submodule(m, "rdGeometry", "rdkit.Geometry.rdGeometry", geometry::init)?;
    add_submodule(m, "rdchem", "rdkit.Chem.rdchem", rdchem::init)?;
    add_submodule(m, "rdmolfiles", "rdkit.Chem.rdmolfiles", rdmolfiles::init)?;
    add_submodule(m, "rdmolops", "rdkit.Chem.rdmolops", rdmolops::init)?;
    Ok(())
}
