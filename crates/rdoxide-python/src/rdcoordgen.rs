//! `rdoxide.Chem.rdCoordGen` (upstream `External/CoordGen/Wrap/rdCoordGen.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use rdoxide::external::coordgen::add_coords::{CoordGenParams, add_coords};
use rdoxide::geometry::Point2D;
use std::collections::BTreeMap;

use crate::errors::IntoPyResult;
use crate::geometry::PyPoint2D;
use crate::mol::PyMol;

/// Parameters controlling coordinate generation
#[pyclass(name = "CoordGenParams", module = "rdoxide.Chem.rdCoordGen")]
#[derive(Debug)]
struct PyCoordGenParams {
    coord_map: BTreeMap<u32, Point2D>,
    template_mol: Option<Py<PyMol>>,
    #[pyo3(get, set, name = "coordgenScaling")]
    coordgen_scaling: f64,
    #[pyo3(get, set, name = "dbg_useConstrained")]
    dbg_use_constrained: bool,
    #[pyo3(get, set, name = "dbg_useFixed")]
    dbg_use_fixed: bool,
    #[pyo3(get, set, name = "templateFileDir")]
    template_file_dir: String,
    #[pyo3(get, set, name = "minimizerPrecision")]
    minimizer_precision: f32,
    #[pyo3(get, set, name = "treatNonterminalBondsToMetalAsZOBs")]
    treat_metal_zob: bool,
}

#[pymethods]
impl PyCoordGenParams {
    #[new]
    fn new() -> Self {
        let d = CoordGenParams::default();
        PyCoordGenParams {
            coord_map: BTreeMap::new(),
            template_mol: None,
            coordgen_scaling: d.coordgen_scaling,
            dbg_use_constrained: d.dbg_use_constrained,
            dbg_use_fixed: d.dbg_use_fixed,
            template_file_dir: String::new(),
            minimizer_precision: d.minimizer_precision,
            treat_metal_zob: false,
        }
    }

    /// expects a dictionary of Point2D objects with template coordinates
    #[pyo3(name = "SetCoordMap", signature = (coordMap))]
    #[allow(non_snake_case)]
    fn set_coord_map(&mut self, coordMap: &Bound<'_, PyDict>) -> PyResult<()> {
        self.coord_map.clear();
        for (k, v) in coordMap.iter() {
            let id: u32 = k.extract()?;
            let p: PyRef<'_, PyPoint2D> = v.extract()?;
            self.coord_map.insert(id, p.p);
        }
        Ok(())
    }

    /// sets a molecule to be used as the template
    #[pyo3(name = "SetTemplateMol", signature = (templ))]
    fn set_template_mol(&mut self, templ: Py<PyMol>) {
        self.template_mol = Some(templ);
    }

    /// highest quality (and slowest) precision setting
    #[getter(sketcherBestPrecision)]
    fn best(&self) -> f32 {
        CoordGenParams::SKETCHER_BEST_PRECISION
    }

    /// standard quality precision setting, the default for the coordgen project
    #[getter(sketcherStandardPrecision)]
    fn standard(&self) -> f32 {
        CoordGenParams::SKETCHER_STANDARD_PRECISION
    }

    /// faster precision setting
    #[getter(sketcherQuickPrecision)]
    fn quick(&self) -> f32 {
        CoordGenParams::SKETCHER_QUICK_PRECISION
    }

    /// "coarse" (fastest) precision setting, produces good-quality coordinates most of the time, this is the default
    /// setting for the RDKit
    #[getter(sketcherCoarsePrecision)]
    fn coarse(&self) -> f32 {
        CoordGenParams::SKETCHER_COARSE_PRECISION
    }
}

/// Add 2D coordinates.
/// ARGUMENTS:
///    - mol: molecule to modify
///    - params: (optional) parameters controlling the coordinate generation
#[pyfunction]
#[pyo3(name = "AddCoords", signature = (mol, params = None))]
fn add_coords_py(
    py: Python<'_>,
    mut mol: PyRefMut<'_, PyMol>,
    params: Option<PyRef<'_, PyCoordGenParams>>,
) -> PyResult<()> {
    let tpl = params
        .as_ref()
        .and_then(|p| p.template_mol.as_ref().map(|t| t.clone_ref(py)));
    let tpl_ref = tpl.as_ref().map(|t| t.borrow(py));
    let mut ps = CoordGenParams::default();
    if let Some(p) = &params {
        ps.coord_map = p.coord_map.clone();
        ps.coordgen_scaling = p.coordgen_scaling;
        ps.dbg_use_constrained = p.dbg_use_constrained;
        ps.dbg_use_fixed = p.dbg_use_fixed;
        ps.template_file_dir = p.template_file_dir.clone();
        ps.minimizer_precision = p.minimizer_precision;
        ps.treat_nonterminal_bonds_to_metal_as_zero_order = p.treat_metal_zob;
    }
    for &i in ps.coord_map.keys() {
        if i >= mol.mol.num_atoms() {
            return Err(PyValueError::new_err("atom index out of range"));
        }
    }
    ps.template_mol = tpl_ref.as_ref().map(|t| &t.mol);
    add_coords(&mut mol.mol, &ps).py()?;
    Ok(())
}

/// Sets the directory containing the templates.mae file (the templates are built in; accepted for compatibility)
#[pyfunction]
#[pyo3(name = "SetDefaultTemplateFileDir", signature = (dir))]
fn set_default_template_file_dir(dir: &str) {
    let _ = dir;
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing interface to the CoordGen library.",
    )?;
    m.add_class::<PyCoordGenParams>()?;
    m.add_function(wrap_pyfunction!(add_coords_py, m)?)?;
    m.add_function(wrap_pyfunction!(set_default_template_file_dir, m)?)?;
    Ok(())
}
