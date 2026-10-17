//! `MolFromRDKitSVG` (upstream `Wrap/rdmolfiles.cpp`).

use pyo3::prelude::*;
use rdoxide::graphmol::file_parsers::svg_parser::rdkit_svg_to_mol;

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Construct a molecule from an RDKit-generated SVG string.
///
///   ARGUMENTS:
///
///     - svg: string containing the SVG data (must include molecule metadata)
///
///     - sanitize: (optional) toggles sanitization of the molecule.
///       Defaults to true.
///
///     - removeHs: (optional) toggles removing hydrogens from the molecule.
///       This only make sense when sanitization is done.
///       Defaults to true.
///
///   RETURNS:
///
///     a Mol object, None on failure.
///
///   NOTE: this functionality should be considered beta.
#[pyfunction]
#[pyo3(name = "MolFromRDKitSVG", signature = (svg, sanitize = true, removeHs = true))]
#[allow(non_snake_case)]
fn py_mol_from_rdkit_svg(
    py: Python<'_>,
    #[pyo3(from_py_with = crate::pyconv::text_string)] svg: String,
    #[pyo3(from_py_with = crate::pyconv::truthy)] sanitize: bool,
    #[pyo3(from_py_with = crate::pyconv::truthy)] removeHs: bool,
) -> PyResult<Option<PyMol>> {
    let res = py
        .detach(|| rdkit_svg_to_mol(&svg, sanitize, removeHs))
        .py()?;
    Ok(res.map(|mol| PyMol { mol }))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_mol_from_rdkit_svg, m)?)
}
