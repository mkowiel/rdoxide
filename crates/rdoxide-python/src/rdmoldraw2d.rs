//! `rdoxide.Chem.Draw.rdMolDraw2D` (upstream `MolDraw2D/Wrap/rdMolDraw2D.cpp`):
//! the options, the SVG drawer and the helper functions. Cairo output is not
//! ported.

use std::collections::BTreeMap;

use pyo3::exceptions::{PyKeyError, PyRuntimeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rdoxide::graphmol::mol_draw2d::contour::{
    ContourParams, contour_and_draw_gaussians, contour_and_draw_grid,
};
use rdoxide::graphmol::mol_draw2d::draw_mol::DrawMolHighlights;
use rdoxide::graphmol::mol_draw2d::drawer::{MolDraw2DSVG, MoleculeHighlights};
use rdoxide::graphmol::mol_draw2d::helpers::{
    DrawColour, LegendPosition, MolDrawOptions, MultiColourHighlightStyle, TextAlignType,
    avalon_palette, bw_palette, cdk_palette, default_palette,
};
use rdoxide::graphmol::mol_draw2d::options_json::update_mol_draw_options_from_json;
use rdoxide::graphmol::mol_draw2d::utils::{
    mean_bond_length, prepare_mol_for_drawing, set_acs1996_options,
};
use rdoxide::graphmol::{ROMol, RWMol};

use crate::enums::{add_enum_exported, to_py};
use crate::errors::IntoPyResult;
use crate::geometry::PyPoint2D;
use crate::mol::PyMol;
use crate::rdchemreactions::PyChemicalReaction;

// ----------------------------------------------------------------------
// conversions
// ----------------------------------------------------------------------

/// `pyTupleToDrawColour`: an (R, G, B[, A]) tuple with values in [0, 1]. The
/// values go through `float` like upstream.
fn py_colour(obj: &Bound<'_, PyAny>) -> PyResult<DrawColour> {
    let t = obj
        .cast::<PyTuple>()
        .map_err(|_| PyTypeError::new_err("colour must be a tuple"))?;
    let bad = || PyValueError::new_err("RGBA color value needs to be between 0 and 1.");
    let mut comps = [0.0f64, 0.0, 0.0, 1.0];
    for (i, comp) in comps.iter_mut().enumerate() {
        if i == 3 && t.len() <= 3 {
            break;
        }
        let v: f32 = t.get_item(i)?.extract()?;
        if !(0.0..=1.0).contains(&v) {
            return Err(bad());
        }
        *comp = f64::from(v);
    }
    Ok(DrawColour::with_alpha(
        comps[0], comps[1], comps[2], comps[3],
    ))
}

fn colour_tuple<'py>(py: Python<'py>, c: DrawColour) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(py, [c.r, c.g, c.b, c.a])
}

/// A colour for the colour dictionaries: no range check, and `float` precision.
fn dict_colour(obj: &Bound<'_, PyAny>) -> PyResult<DrawColour> {
    let t = obj
        .cast::<PyTuple>()
        .map_err(|_| PyTypeError::new_err("colour must be a tuple"))?;
    let c = |i: usize| -> PyResult<f64> { Ok(f64::from(t.get_item(i)?.extract::<f32>()?)) };
    let a = if t.len() > 3 { c(3)? } else { 1.0 };
    Ok(DrawColour::with_alpha(c(0)?, c(1)?, c(2)?, a))
}

fn truthy(obj: &Bound<'_, PyAny>) -> PyResult<bool> {
    obj.is_truthy()
}

/// An optional argument: `None` (or an empty container) means not given.
fn given<'a, 'py>(obj: Option<&'a Bound<'py, PyAny>>) -> PyResult<Option<&'a Bound<'py, PyAny>>> {
    match obj {
        Some(o) if !o.is_none() && truthy(o)? => Ok(Some(o)),
        _ => Ok(None),
    }
}

fn colour_map(obj: Option<&Bound<'_, PyAny>>) -> PyResult<Option<BTreeMap<i32, DrawColour>>> {
    let Some(o) = given(obj)? else {
        return Ok(None);
    };
    let d = o.cast::<PyDict>()?;
    let mut res = BTreeMap::new();
    for (k, v) in d.iter() {
        res.insert(k.extract()?, dict_colour(&v)?);
    }
    Ok(Some(res))
}

fn double_map(obj: Option<&Bound<'_, PyAny>>) -> PyResult<Option<BTreeMap<i32, f64>>> {
    let Some(o) = given(obj)? else {
        return Ok(None);
    };
    let d = o.cast::<PyDict>()?;
    let mut res = BTreeMap::new();
    for (k, v) in d.iter() {
        res.insert(k.extract()?, v.extract()?);
    }
    Ok(Some(res))
}

/// `pythonObjectToVect(obj, maxV)`: the elements must be below `max`.
fn int_list(obj: Option<&Bound<'_, PyAny>>, max: Option<i32>) -> PyResult<Option<Vec<i32>>> {
    let Some(o) = given(obj)? else {
        return Ok(None);
    };
    let mut res = Vec::new();
    for item in o.try_iter()? {
        let v: i32 = item?.extract()?;
        if max.is_some_and(|m| v >= m) {
            return Err(PyValueError::new_err(
                "list element larger than allowed value",
            ));
        }
        res.push(v);
    }
    Ok(Some(res))
}

/// Each element of a list of per-molecule values, checked against the number
/// of molecules.
fn per_mol<'py>(
    obj: Option<&Bound<'py, PyAny>>,
    n: usize,
    name: &str,
) -> PyResult<Option<Vec<Bound<'py, PyAny>>>> {
    let Some(o) = given(obj)? else {
        return Ok(None);
    };
    let items = o.try_iter()?.collect::<PyResult<Vec<_>>>()?;
    if items.len() != n {
        return Err(PyValueError::new_err(format!(
            "If {name} is provided it must be the same length as the molecule list."
        )));
    }
    Ok(Some(items))
}

fn highlights(
    mol: &ROMol,
    atoms: Option<&Bound<'_, PyAny>>,
    bonds: Option<&Bound<'_, PyAny>>,
    atom_colours: Option<&Bound<'_, PyAny>>,
    bond_colours: Option<&Bound<'_, PyAny>>,
    radii: Option<&Bound<'_, PyAny>>,
) -> PyResult<DrawMolHighlights> {
    Ok(DrawMolHighlights {
        atoms: int_list(atoms, Some(mol.num_atoms() as i32))?.unwrap_or_default(),
        bonds: int_list(bonds, Some(mol.num_bonds() as i32))?.unwrap_or_default(),
        atom_map: colour_map(atom_colours)?.unwrap_or_default(),
        bond_map: colour_map(bond_colours)?.unwrap_or_default(),
        radii: double_map(radii)?.unwrap_or_default(),
        ..Default::default()
    })
}

/// Binds positional and keyword arguments to `names`.
fn bind_args<'py>(
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
    names: &[&str],
    fname: &str,
) -> PyResult<Vec<Option<Bound<'py, PyAny>>>> {
    if args.len() > names.len() {
        return Err(PyTypeError::new_err(format!(
            "{fname}() takes at most {} arguments ({} given)",
            names.len() + 1,
            args.len() + 1
        )));
    }
    let mut res: Vec<Option<Bound<'py, PyAny>>> = vec![None; names.len()];
    for (i, a) in args.iter().enumerate() {
        res[i] = Some(a);
    }
    if let Some(kw) = kwargs {
        for (k, v) in kw.iter() {
            let key: String = k.extract()?;
            let Some(pos) = names.iter().position(|n| *n == key) else {
                return Err(PyTypeError::new_err(format!(
                    "{fname}() got an unexpected keyword argument '{key}'"
                )));
            };
            if res[pos].is_some() {
                return Err(PyTypeError::new_err(format!(
                    "{fname}() got multiple values for argument '{key}'"
                )));
            }
            res[pos] = Some(v);
        }
    }
    Ok(res)
}

fn opt_i32(v: &Option<Bound<'_, PyAny>>, default: i32) -> PyResult<i32> {
    v.as_ref().map_or(Ok(default), |o| o.extract())
}

fn opt_str(v: &Option<Bound<'_, PyAny>>) -> PyResult<String> {
    v.as_ref().map_or(Ok(String::new()), |o| o.extract())
}

fn mol_arg<'a>(obj: &'a Bound<'_, PyAny>) -> PyResult<PyRef<'a, PyMol>> {
    Ok(obj.extract()?)
}

// ----------------------------------------------------------------------
// MolDrawOptions
// ----------------------------------------------------------------------

/// Drawing options
#[pyclass(name = "MolDrawOptions", module = "rdoxide.Chem.Draw.rdMolDraw2D")]
#[derive(Debug)]
pub struct PyMolDrawOptions {
    /// The drawer the options belong to; they are standalone otherwise.
    owner: Option<Py<PyMolDraw2D>>,
    own: MolDrawOptions,
}

impl PyMolDrawOptions {
    fn read<R>(&self, py: Python<'_>, f: impl FnOnce(&MolDrawOptions) -> R) -> PyResult<R> {
        match &self.owner {
            Some(o) => Ok(f(o.bind(py).try_borrow()?.d.draw_options())),
            None => Ok(f(&self.own)),
        }
    }

    fn write<R>(
        &mut self,
        py: Python<'_>,
        f: impl FnOnce(&mut MolDrawOptions) -> R,
    ) -> PyResult<R> {
        match &self.owner {
            Some(o) => Ok(f(o.bind(py).try_borrow_mut()?.d.draw_options_mut())),
            None => Ok(f(&mut self.own)),
        }
    }
}

macro_rules! option_props {
    ($( ($py:literal, $get:ident, $set:ident, $field:ident, $ty:ty) ),* $(,)?) => {
        #[pymethods]
        impl PyMolDrawOptions {
            $(
                #[getter($py)]
                fn $get(&self, py: Python<'_>) -> PyResult<$ty> {
                    self.read(py, |o| o.$field.clone())
                }
                #[setter($py)]
                fn $set(&mut self, py: Python<'_>, value: $ty) -> PyResult<()> {
                    self.write(py, |o| o.$field = value)
                }
            )*
        }
    };
}

option_props!(
    ("dummiesAreAttachments", get_dummies_are_attachments, set_dummies_are_attachments, dummies_are_attachments, bool),
    ("circleAtoms", get_circle_atoms, set_circle_atoms, circle_atoms, bool),
    ("splitBonds", get_split_bonds, set_split_bonds, split_bonds, bool),
    ("atomLabelDeuteriumTritium", get_atom_label_dt, set_atom_label_dt, atom_label_deuterium_tritium, bool),
    ("continuousHighlight", get_continuous_highlight, set_continuous_highlight, continuous_highlight, bool),
    ("fillHighlights", get_fill_highlights, set_fill_highlights, fill_highlights, bool),
    ("highlightRadius", get_highlight_radius, set_highlight_radius, highlight_radius, f64),
    ("flagCloseContactsDist", get_flag_close_contacts_dist, set_flag_close_contacts_dist, flag_close_contacts_dist, i32),
    ("includeAtomTags", get_include_atom_tags, set_include_atom_tags, include_atom_tags, bool),
    ("clearBackground", get_clear_background, set_clear_background, clear_background, bool),
    ("legendFontSize", get_legend_font_size, set_legend_font_size, legend_font_size, i32),
    ("legendFraction", get_legend_fraction, set_legend_fraction, legend_fraction, f64),
    ("legendVerticalText", get_legend_vertical_text, set_legend_vertical_text, legend_vertical_text, bool),
    ("maxFontSize", get_max_font_size, set_max_font_size, max_font_size, i32),
    ("minFontSize", get_min_font_size, set_min_font_size, min_font_size, i32),
    ("fixedFontSize", get_fixed_font_size, set_fixed_font_size, fixed_font_size, i32),
    ("baseFontSize", get_base_font_size, set_base_font_size, base_font_size, f64),
    ("annotationFontScale", get_annotation_font_scale, set_annotation_font_scale, annotation_font_scale, f64),
    ("fontFile", get_font_file, set_font_file, font_file, String),
    ("multipleBondOffset", get_multiple_bond_offset, set_multiple_bond_offset, multiple_bond_offset, f64),
    ("padding", get_padding, set_padding, padding, f64),
    ("reagentPadding", get_reagent_padding, set_reagent_padding, component_padding, f64),
    ("bondLineWidth", get_bond_line_width, set_bond_line_width, bond_line_width, f64),
    ("scaleBondWidth", get_scale_bond_width, set_scale_bond_width, scale_bond_width, bool),
    ("scaleHighlightBondWidth", get_scale_hl_bond_width, set_scale_hl_bond_width, scale_highlight_bond_width, bool),
    ("highlightBondWidthMultiplier", get_hl_bond_width_mult, set_hl_bond_width_mult, highlight_bond_width_multiplier, i32),
    ("prepareMolsBeforeDrawing", get_prepare_mols, set_prepare_mols, prepare_mols_before_drawing, bool),
    ("fixedScale", get_fixed_scale, set_fixed_scale, fixed_scale, f64),
    ("fixedBondLength", get_fixed_bond_length, set_fixed_bond_length, fixed_bond_length, f64),
    ("rotate", get_rotate, set_rotate, rotate, f64),
    ("addStereoAnnotation", get_add_stereo_annotation, set_add_stereo_annotation, add_stereo_annotation, bool),
    ("showAllCIPCodes", get_show_all_cip_codes, set_show_all_cip_codes, show_all_cip_codes, bool),
    ("addAtomIndices", get_add_atom_indices, set_add_atom_indices, add_atom_indices, bool),
    ("addBondIndices", get_add_bond_indices, set_add_bond_indices, add_bond_indices, bool),
    ("isotopeLabels", get_isotope_labels, set_isotope_labels, isotope_labels, bool),
    ("dummyIsotopeLabels", get_dummy_isotope_labels, set_dummy_isotope_labels, dummy_isotope_labels, bool),
    ("atomHighlightsAreCircles", get_atom_hl_circles, set_atom_hl_circles, atom_highlights_are_circles, bool),
    ("centreMoleculesBeforeDrawing", get_centre_mols, set_centre_mols, centre_molecules_before_drawing, bool),
    ("additionalAtomLabelPadding", get_add_atom_label_padding, set_add_atom_label_padding, additional_atom_label_padding, f64),
    ("noAtomLabels", get_no_atom_labels, set_no_atom_labels, no_atom_labels, bool),
    ("explicitMethyl", get_explicit_methyl, set_explicit_methyl, explicit_methyl, bool),
    ("includeMetadata", get_include_metadata, set_include_metadata, include_metadata, bool),
    ("includeRadicals", get_include_radicals, set_include_radicals, include_radicals, bool),
    ("comicMode", get_comic_mode, set_comic_mode, comic_mode, bool),
    ("variableBondWidthMultiplier", get_var_bond_width_mult, set_var_bond_width_mult, variable_bond_width_multiplier, i32),
    ("variableAtomRadius", get_variable_atom_radius, set_variable_atom_radius, variable_atom_radius, f64),
    ("includeChiralFlagLabel", get_include_chiral_flag, set_include_chiral_flag, include_chiral_flag_label, bool),
    ("simplifiedStereoGroupLabel", get_simplified_sg_label, set_simplified_sg_label, simplified_stereo_group_label, bool),
    ("stereoGroupAndLabel", get_sg_and_label, set_sg_and_label, stereo_group_and_label, String),
    ("stereoGroupOrLabel", get_sg_or_label, set_sg_or_label, stereo_group_or_label, String),
    ("stereoGroupAbsLabel", get_sg_abs_label, set_sg_abs_label, stereo_group_abs_label, String),
    ("addStereoGroupAnnotation", get_add_sg_annotation, set_add_sg_annotation, add_stereo_group_annotation, bool),
    ("unspecifiedStereoIsUnknown", get_unspec_stereo_unknown, set_unspec_stereo_unknown, unspecified_stereo_is_unknown, bool),
    ("singleColourWedgeBonds", get_single_colour_wedges, set_single_colour_wedges, single_colour_wedge_bonds, bool),
    ("singleColourBonds", get_single_colour_bonds, set_single_colour_bonds, single_colour_bonds, bool),
    ("useMolBlockWedging", get_use_mb_wedging, set_use_mb_wedging, use_mol_block_wedging, bool),
    ("scalingFactor", get_scaling_factor, set_scaling_factor, scaling_factor, f64),
    ("drawMolsSameScale", get_draw_mols_same_scale, set_draw_mols_same_scale, draw_mols_same_scale, bool),
    ("useComplexQueryAtomSymbols", get_complex_query_symbols, set_complex_query_symbols, use_complex_query_atom_symbols, bool),
    ("bracketsAroundAtomLists", get_brackets_atom_lists, set_brackets_atom_lists, brackets_around_atom_lists, bool),
    ("standardColoursForHighlightedAtoms", get_std_colours_hl, set_std_colours_hl, standard_colours_for_highlighted_atoms, bool),
    ("drawingExtentsInclude", get_extents_include, set_extents_include, drawing_extents_include, u32),
    ("atomRegions", get_atom_regions, set_atom_regions, atom_regions, Vec<Vec<i32>>),
);

macro_rules! colour_options {
    ($( ($py:literal, $field:ident, $getm:literal, $setm:literal, $getp:ident, $setp:ident, $getf:ident, $setf:ident) ),* $(,)?) => {
        #[pymethods]
        impl PyMolDrawOptions {
            $(
                #[getter($py)]
                fn $getp<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
                    let c = self.read(py, |o| o.$field)?;
                    colour_tuple(py, c)
                }
                #[setter($py)]
                fn $setp(&mut self, py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<()> {
                    let c = py_colour(value)?;
                    self.write(py, |o| o.$field = c)
                }
                #[pyo3(name = $getm)]
                fn $getf<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
                    let c = self.read(py, |o| o.$field)?;
                    colour_tuple(py, c)
                }
                #[pyo3(name = $setm, signature = (tpl))]
                fn $setf(&mut self, py: Python<'_>, tpl: &Bound<'_, PyAny>) -> PyResult<()> {
                    let c = py_colour(tpl)?;
                    self.write(py, |o| o.$field = c)
                }
            )*
        }
    };
}

colour_options!(
    (
        "backgroundColour",
        background_colour,
        "getBackgroundColour",
        "setBackgroundColour",
        get_bg_colour_p,
        set_bg_colour_p,
        get_bg_colour,
        set_bg_colour
    ),
    (
        "queryColour",
        query_colour,
        "getQueryColour",
        "setQueryColour",
        get_qy_colour_p,
        set_qy_colour_p,
        get_qy_colour,
        set_qy_colour
    ),
    (
        "highlightColour",
        highlight_colour,
        "getHighlightColour",
        "setHighlightColour",
        get_hl_colour_p,
        set_hl_colour_p,
        get_hl_colour,
        set_hl_colour
    ),
    (
        "symbolColour",
        symbol_colour,
        "getSymbolColour",
        "setSymbolColour",
        get_sym_colour_p,
        set_sym_colour_p,
        get_sym_colour,
        set_sym_colour
    ),
    (
        "annotationColour",
        annotation_colour,
        "getAnnotationColour",
        "setAnnotationColour",
        get_ann_colour_p,
        set_ann_colour_p,
        get_ann_colour,
        set_ann_colour
    ),
    (
        "atomNoteColour",
        atom_note_colour,
        "getAtomNoteColour",
        "setAtomNoteColour",
        get_an_colour_p,
        set_an_colour_p,
        get_an_colour,
        set_an_colour
    ),
    (
        "bondNoteColour",
        bond_note_colour,
        "getBondNoteColour",
        "setBondNoteColour",
        get_bn_colour_p,
        set_bn_colour_p,
        get_bn_colour,
        set_bn_colour
    ),
    (
        "legendColour",
        legend_colour,
        "getLegendColour",
        "setLegendColour",
        get_leg_colour_p,
        set_leg_colour_p,
        get_leg_colour,
        set_leg_colour
    ),
    (
        "variableAttachmentColour",
        variable_attachment_colour,
        "getVariableAttachmentColour",
        "setVariableAttachmentColour",
        get_va_colour_p,
        set_va_colour_p,
        get_va_colour,
        set_va_colour
    ),
);

#[pymethods]
#[allow(non_snake_case)]
impl PyMolDrawOptions {
    #[new]
    fn new() -> Self {
        PyMolDrawOptions {
            owner: None,
            own: MolDrawOptions::default(),
        }
    }

    #[getter(legendPosition)]
    fn get_legend_position(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let v = self.read(py, |o| o.legend_position)?;
        to_py(py, "LegendPosition", v as i64)
    }

    #[setter(legendPosition)]
    fn set_legend_position(&mut self, py: Python<'_>, value: i64) -> PyResult<()> {
        let v = match value {
            0 => LegendPosition::Bottom,
            1 => LegendPosition::Top,
            2 => LegendPosition::Left,
            3 => LegendPosition::Right,
            _ => return Err(PyValueError::new_err("invalid LegendPosition value")),
        };
        self.write(py, |o| o.legend_position = v)
    }

    #[getter(multiColourHighlightStyle)]
    fn get_mch_style(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let v = self.read(py, |o| o.multi_colour_highlight_style)?;
        to_py(py, "MultiColourHighlightStyle", v as i64)
    }

    #[setter(multiColourHighlightStyle)]
    fn set_mch_style(&mut self, py: Python<'_>, value: i64) -> PyResult<()> {
        let v = match value {
            0 => MultiColourHighlightStyle::CircleAndLine,
            1 => MultiColourHighlightStyle::Lasso,
            _ => {
                return Err(PyValueError::new_err(
                    "invalid MultiColourHighlightStyle value",
                ));
            }
        };
        self.write(py, |o| o.multi_colour_highlight_style = v)
    }

    /// maps indices to atom labels
    #[getter(atomLabels)]
    fn get_atom_labels(slf: Bound<'_, Self>) -> PyIntStringMap {
        PyIntStringMap {
            options: slf.unbind(),
        }
    }

    #[setter(atomLabels)]
    fn set_atom_labels(&mut self, py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let map: BTreeMap<i32, String> = value.extract()?;
        self.write(py, |o| o.atom_labels = map)
    }

    /// use the default colour palette for atoms and bonds
    fn useDefaultAtomPalette(&mut self, py: Python<'_>) -> PyResult<()> {
        self.write(py, |o| o.atom_colour_palette = default_palette())
    }

    /// use a black and white palette for atoms and bonds
    fn useBWAtomPalette(&mut self, py: Python<'_>) -> PyResult<()> {
        self.write(py, |o| o.atom_colour_palette = bw_palette())
    }

    /// use the Avalon renderer palette for atoms and bonds
    fn useAvalonAtomPalette(&mut self, py: Python<'_>) -> PyResult<()> {
        self.write(py, |o| o.atom_colour_palette = avalon_palette())
    }

    /// use the CDK palette for atoms and bonds
    fn useCDKAtomPalette(&mut self, py: Python<'_>) -> PyResult<()> {
        self.write(py, |o| o.atom_colour_palette = cdk_palette())
    }

    /// updates the palette for atoms and bonds from a dictionary mapping ints to 3-tuples
    #[pyo3(signature = (cmap))]
    fn updateAtomPalette(&mut self, py: Python<'_>, cmap: &Bound<'_, PyAny>) -> PyResult<()> {
        let map = colour_map_required(cmap)?;
        self.write(py, |o| o.atom_colour_palette.extend(map))
    }

    /// sets the palette for atoms and bonds from a dictionary mapping ints to 3-tuples
    #[pyo3(signature = (cmap))]
    fn setAtomPalette(&mut self, py: Python<'_>, cmap: &Bound<'_, PyAny>) -> PyResult<()> {
        let map = colour_map_required(cmap)?;
        self.write(py, |o| o.atom_colour_palette = map.into_iter().collect())
    }

    /// returns the current atom palette as a dictionary mapping ints to 4-tuples
    fn getAtomPalette<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let palette = self.read(py, |o| o.atom_colour_palette.clone())?;
        let res = PyDict::new(py);
        for (k, v) in palette {
            res.set_item(k, colour_tuple(py, v)?)?;
        }
        Ok(res)
    }
}

fn colour_map_required(obj: &Bound<'_, PyAny>) -> PyResult<BTreeMap<i32, DrawColour>> {
    let d = obj.cast::<PyDict>()?;
    let mut res = BTreeMap::new();
    for (k, v) in d.iter() {
        res.insert(k.extract()?, dict_colour(&v)?);
    }
    Ok(res)
}

/// A view of `MolDrawOptions.atomLabels` that writes through.
#[pyclass(name = "IntStringMap", module = "rdoxide.Chem.Draw.rdMolDraw2D")]
#[derive(Debug)]
pub struct PyIntStringMap {
    options: Py<PyMolDrawOptions>,
}

impl PyIntStringMap {
    fn snapshot(&self, py: Python<'_>) -> PyResult<BTreeMap<i32, String>> {
        self.options
            .bind(py)
            .try_borrow()?
            .read(py, |o| o.atom_labels.clone())
    }
}

#[pymethods]
impl PyIntStringMap {
    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        Ok(self.snapshot(py)?.len())
    }

    fn __getitem__(&self, py: Python<'_>, key: i32) -> PyResult<String> {
        self.snapshot(py)?
            .get(&key)
            .cloned()
            .ok_or_else(|| PyKeyError::new_err(key))
    }

    fn __setitem__(&self, py: Python<'_>, key: i32, value: String) -> PyResult<()> {
        self.options.bind(py).try_borrow_mut()?.write(py, |o| {
            o.atom_labels.insert(key, value);
        })
    }

    fn __delitem__(&self, py: Python<'_>, key: i32) -> PyResult<()> {
        let removed = self
            .options
            .bind(py)
            .try_borrow_mut()?
            .write(py, |o| o.atom_labels.remove(&key))?;
        removed.map(|_| ()).ok_or_else(|| PyKeyError::new_err(key))
    }

    fn __contains__(&self, py: Python<'_>, key: i32) -> PyResult<bool> {
        Ok(self.snapshot(py)?.contains_key(&key))
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        PyList::new(py, self.snapshot(py)?.keys().copied())?
            .into_any()
            .try_iter()
            .map(|i| i.into_any())
    }

    fn keys(&self, py: Python<'_>) -> PyResult<Vec<i32>> {
        Ok(self.snapshot(py)?.keys().copied().collect())
    }

    fn values(&self, py: Python<'_>) -> PyResult<Vec<String>> {
        Ok(self.snapshot(py)?.values().cloned().collect())
    }

    fn items(&self, py: Python<'_>) -> PyResult<Vec<(i32, String)>> {
        Ok(self.snapshot(py)?.into_iter().collect())
    }

    #[pyo3(signature = (key, default = None))]
    fn get(&self, py: Python<'_>, key: i32, default: Option<String>) -> PyResult<Option<String>> {
        Ok(self.snapshot(py)?.get(&key).cloned().or(default))
    }
}

// ----------------------------------------------------------------------
// the drawers
// ----------------------------------------------------------------------

/// Base class for the molecule drawers.
#[pyclass(name = "MolDraw2D", module = "rdoxide.Chem.Draw.rdMolDraw2D", subclass)]
#[derive(Debug)]
pub struct PyMolDraw2D {
    d: MolDraw2DSVG,
}

fn align_from_int(align: i32) -> PyResult<TextAlignType> {
    match align {
        0 => Ok(TextAlignType::Middle),
        1 => Ok(TextAlignType::Start),
        2 => Ok(TextAlignType::End),
        _ => Err(PyValueError::new_err("align must be 0, 1, or 2")),
    }
}

#[pymethods]
#[allow(non_snake_case)]
impl PyMolDraw2D {
    /// change the default font size. The units are, roughly, pixels.
    #[pyo3(signature = (new_size))]
    fn SetFontSize(&mut self, new_size: f64) {
        self.d.set_font_size(new_size);
    }

    /// get the default font size. The units are, roughly, pixels.
    fn FontSize(&self) -> f64 {
        self.d.font_size()
    }

    /// renders a molecule
    #[pyo3(signature = (mol, *args, **kwargs))]
    fn DrawMolecule(
        &mut self,
        py: Python<'_>,
        mol: &Bound<'_, PyAny>,
        args: &Bound<'_, PyTuple>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<()> {
        let m = mol_arg(mol)?;
        // the overload with bonds has them as the second highlight argument
        let with_bonds = kwargs.is_some_and(|k| {
            k.contains("highlightBonds").unwrap_or(false)
                || k.contains("highlightBondColors").unwrap_or(false)
        }) || args
            .get_item(1)
            .is_ok_and(|a| !a.is_none() && !a.is_instance_of::<PyDict>());
        let (hl, conf_id, legend) = if with_bonds {
            let a = bind_args(
                args,
                kwargs,
                &[
                    "highlightAtoms",
                    "highlightBonds",
                    "highlightAtomColors",
                    "highlightBondColors",
                    "highlightAtomRadii",
                    "confId",
                    "legend",
                ],
                "DrawMolecule",
            )?;
            (
                highlights(
                    &m.mol,
                    a[0].as_ref(),
                    a[1].as_ref(),
                    a[2].as_ref(),
                    a[3].as_ref(),
                    a[4].as_ref(),
                )?,
                opt_i32(&a[5], -1)?,
                opt_str(&a[6])?,
            )
        } else {
            let a = bind_args(
                args,
                kwargs,
                &[
                    "highlightAtoms",
                    "highlightAtomColors",
                    "highlightAtomRadii",
                    "confId",
                    "legend",
                ],
                "DrawMolecule",
            )?;
            let atoms = int_list(a[0].as_ref(), Some(m.mol.num_atoms() as i32))?;
            let bonds = atoms
                .as_ref()
                .map(|at| {
                    rdoxide::graphmol::mol_draw2d::draw_mol_utils::get_bond_highlights_for_atoms(
                        &m.mol, at,
                    )
                })
                .unwrap_or_default();
            (
                DrawMolHighlights {
                    atoms: atoms.unwrap_or_default(),
                    bonds,
                    atom_map: colour_map(a[1].as_ref())?.unwrap_or_default(),
                    radii: double_map(a[2].as_ref())?.unwrap_or_default(),
                    ..Default::default()
                },
                opt_i32(&a[3], -1)?,
                opt_str(&a[4])?,
            )
        };
        self.d
            .draw_molecule_with(&m.mol, &legend, hl, conf_id)
            .py()?;
        let _ = py;
        Ok(())
    }

    /// returns the width and height required to draw a molecule at the current size
    #[pyo3(signature = (mol, highlightAtoms = None, highlightBonds = None, highlightAtomColors = None, highlightBondColors = None, highlightAtomRadii = None, confId = -1, legend = ""))]
    fn GetMolSize(
        &mut self,
        mol: PyRef<'_, PyMol>,
        highlightAtoms: Option<&Bound<'_, PyAny>>,
        highlightBonds: Option<&Bound<'_, PyAny>>,
        highlightAtomColors: Option<&Bound<'_, PyAny>>,
        highlightBondColors: Option<&Bound<'_, PyAny>>,
        highlightAtomRadii: Option<&Bound<'_, PyAny>>,
        confId: i32,
        legend: &str,
    ) -> PyResult<(i32, i32)> {
        let hl = highlights(
            &mol.mol,
            highlightAtoms,
            highlightBonds,
            highlightAtomColors,
            highlightBondColors,
            highlightAtomRadii,
        )?;
        self.d.get_mol_size(&mol.mol, legend, hl, confId).py()
    }

    /// renders a molecule with multiple highlight colours
    #[pyo3(signature = (mol, legend, highlight_atom_map, highlight_bond_map, highlight_radii, highlight_linewidth_multipliers, confId = -1))]
    fn DrawMoleculeWithHighlights(
        &mut self,
        mol: PyRef<'_, PyMol>,
        legend: &str,
        highlight_atom_map: &Bound<'_, PyAny>,
        highlight_bond_map: &Bound<'_, PyAny>,
        highlight_radii: &Bound<'_, PyAny>,
        highlight_linewidth_multipliers: &Bound<'_, PyAny>,
        confId: i32,
    ) -> PyResult<()> {
        let colour_lists = |obj: &Bound<'_, PyAny>| -> PyResult<BTreeMap<i32, Vec<DrawColour>>> {
            let mut res = BTreeMap::new();
            if given(Some(obj))?.is_some() {
                for (k, v) in obj.cast::<PyDict>()?.iter() {
                    let cols = v
                        .try_iter()?
                        .map(|c| py_colour(&c?))
                        .collect::<PyResult<Vec<_>>>()?;
                    res.insert(k.extract()?, cols);
                }
            }
            Ok(res)
        };
        let atoms = colour_lists(highlight_atom_map)?;
        let bonds = colour_lists(highlight_bond_map)?;
        let radii = double_map(Some(highlight_radii))?.unwrap_or_default();
        let mults: BTreeMap<i32, i32> = if given(Some(highlight_linewidth_multipliers))?.is_some() {
            highlight_linewidth_multipliers.extract()?
        } else {
            BTreeMap::new()
        };
        self.d
            .draw_molecule_with_highlights(&mol.mol, legend, atoms, bonds, radii, mults, confId)
            .py()
    }

    /// renders multiple molecules
    #[pyo3(signature = (mols, highlightAtoms = None, highlightBonds = None, highlightAtomColors = None, highlightBondColors = None, highlightAtomRadii = None, confIds = None, legends = None))]
    fn DrawMolecules(
        &mut self,
        mols: &Bound<'_, PyAny>,
        highlightAtoms: Option<&Bound<'_, PyAny>>,
        highlightBonds: Option<&Bound<'_, PyAny>>,
        highlightAtomColors: Option<&Bound<'_, PyAny>>,
        highlightBondColors: Option<&Bound<'_, PyAny>>,
        highlightAtomRadii: Option<&Bound<'_, PyAny>>,
        confIds: Option<&Bound<'_, PyAny>>,
        legends: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        let mut holders: Vec<Option<PyRef<'_, PyMol>>> = Vec::new();
        for item in mols.try_iter()? {
            let item = item?;
            holders.push(if item.is_none() {
                None
            } else {
                Some(item.extract()?)
            });
        }
        if holders.is_empty() {
            return Ok(());
        }
        let n = holders.len();
        let atoms = per_mol(highlightAtoms, n, "highlightAtoms")?;
        let bonds = per_mol(highlightBonds, n, "highlightBonds")?;
        let atom_maps = per_mol(highlightAtomColors, n, "highlightAtomMap")?;
        let bond_maps = per_mol(highlightBondColors, n, "highlightBondMap")?;
        let radii = per_mol(highlightAtomRadii, n, "highlightAtomRadii")?;
        let any_hl = atoms.is_some()
            || bonds.is_some()
            || atom_maps.is_some()
            || bond_maps.is_some()
            || radii.is_some();
        let mut hls = Vec::with_capacity(n);
        for i in 0..n {
            hls.push(MoleculeHighlights {
                atoms: match &atoms {
                    Some(a) => int_list(Some(&a[i]), None)?.or(Some(Vec::new())),
                    None => None,
                },
                bonds: match &bonds {
                    Some(a) => int_list(Some(&a[i]), None)?.or(Some(Vec::new())),
                    None => None,
                },
                atom_map: match &atom_maps {
                    Some(a) => Some(colour_map(Some(&a[i]))?.unwrap_or_default()),
                    None => None,
                },
                bond_map: match &bond_maps {
                    Some(a) => Some(colour_map(Some(&a[i]))?.unwrap_or_default()),
                    None => None,
                },
                radii: match &radii {
                    Some(a) => Some(double_map(Some(&a[i]))?.unwrap_or_default()),
                    None => None,
                },
            });
        }
        let conf_ids = match given(confIds)? {
            Some(c) => Some(
                c.try_iter()?
                    .map(|v| v?.extract())
                    .collect::<PyResult<Vec<i32>>>()?,
            ),
            None => None,
        };
        let legends = match given(legends)? {
            Some(l) => Some(
                l.try_iter()?
                    .map(|v| v?.extract())
                    .collect::<PyResult<Vec<String>>>()?,
            ),
            None => None,
        };
        let refs: Vec<Option<&ROMol>> =
            holders.iter().map(|h| h.as_ref().map(|m| &m.mol)).collect();
        self.d
            .draw_molecules(
                &refs,
                legends.as_deref(),
                if any_hl { Some(&hls) } else { None },
                conf_ids.as_deref(),
            )
            .py()
    }

    /// renders a reaction
    #[pyo3(signature = (rxn, highlightByReactant = false, highlightColorsReactants = None, confIds = None))]
    fn DrawReaction(
        &mut self,
        rxn: PyRef<'_, PyChemicalReaction>,
        highlightByReactant: bool,
        highlightColorsReactants: Option<&Bound<'_, PyAny>>,
        confIds: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        let colours = match given(highlightColorsReactants)? {
            Some(c) => Some(
                c.try_iter()?
                    .map(|t| py_colour(&t?))
                    .collect::<PyResult<Vec<_>>>()?,
            ),
            None => None,
        };
        let conf_ids = int_list(confIds, None)?;
        self.d
            .draw_reaction(
                &rxn.rxn,
                highlightByReactant,
                colours.as_deref(),
                conf_ids.as_deref(),
            )
            .py()
    }

    /// get the width of the drawing canvas
    fn Width(&self) -> i32 {
        self.d.width()
    }

    /// get the height of the drawing canvas
    fn Height(&self) -> i32 {
        self.d.height()
    }

    /// set the offset (in drawing coordinates) for the drawing
    #[pyo3(signature = (x, y))]
    fn SetOffset(&mut self, x: i32, y: i32) {
        self.d.set_offset(x, y);
    }

    /// returns the offset (in drawing coordinates) for the drawing
    fn Offset(&self) -> PyPoint2D {
        PyPoint2D { p: self.d.offset() }
    }

    /// uses the values provided to set the drawing scaling
    #[pyo3(signature = (width, height, minv, maxv, mol = None))]
    fn SetScale(
        &mut self,
        width: i32,
        height: i32,
        minv: PyRef<'_, PyPoint2D>,
        maxv: PyRef<'_, PyPoint2D>,
        mol: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        let holder = match mol {
            Some(m) if !m.is_none() => Some(mol_arg(m)?),
            _ => None,
        };
        self.d
            .set_scale_for_range(
                width,
                height,
                minv.p,
                maxv.p,
                holder.as_ref().map(|m| &m.mol),
            )
            .py()
    }

    /// returns whether or not FlexiMode is being used
    fn FlexiMode(&self) -> bool {
        self.d.flexi_mode()
    }

    /// when FlexiMode is set, molecules will always been drawn with the default values for bond length, font size, etc.
    #[pyo3(signature = (mode))]
    fn SetFlexiMode(&mut self, mode: bool) {
        self.d.set_flexi_mode(mode);
    }

    /// set the line width being used
    #[pyo3(signature = (width))]
    fn SetLineWidth(&mut self, width: f64) {
        self.d.set_line_width(width);
    }

    /// set the color being used fr drawing and filling
    #[pyo3(signature = (tpl))]
    fn SetColour(&mut self, tpl: &Bound<'_, PyAny>) -> PyResult<()> {
        self.d.set_colour(py_colour(tpl)?);
        Ok(())
    }

    /// returns the line width being used
    fn LineWidth(&self) -> f64 {
        self.d.line_width()
    }

    /// sets whether or not polygons are filled
    #[pyo3(signature = (val))]
    fn SetFillPolys(&mut self, val: bool) {
        self.d.set_fill_polys(val);
    }

    /// returns whether or not polygons are being filled
    fn FillPolys(&self) -> bool {
        self.d.fill_polys()
    }

    /// draws a line with the current drawing style. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
    #[pyo3(signature = (cds1, cds2, rawCoords = false))]
    fn DrawLine(
        &mut self,
        cds1: PyRef<'_, PyPoint2D>,
        cds2: PyRef<'_, PyPoint2D>,
        rawCoords: bool,
    ) -> PyResult<()> {
        self.d.draw_line(cds1.p, cds2.p, rawCoords).py()
    }

    /// draws an arrow with the current drawing style. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels. If asPolygon is true the head of the arrow will be drawn as a triangle, otherwise two lines are used. The fraction of the arrow length to use for the head is given by frac. The angle of the arrowhead (the angle between the main line and each arrowhead line) is given by angle. The color is a tuple of 3 floats (0-1) in red, green, blue (RGB) order.
    #[pyo3(signature = (cds1, cds2, asPolygon = false, frac = 0.05, angle = std::f64::consts::PI / 6.0, color = None, rawCoords = false))]
    fn DrawArrow(
        &mut self,
        cds1: PyRef<'_, PyPoint2D>,
        cds2: PyRef<'_, PyPoint2D>,
        asPolygon: bool,
        frac: f64,
        angle: f64,
        color: Option<&Bound<'_, PyAny>>,
        rawCoords: bool,
    ) -> PyResult<()> {
        let col = match color {
            Some(c) if !c.is_none() => py_colour(c)?,
            _ => DrawColour::new(0.0, 0.0, 0.0),
        };
        self.d
            .draw_arrow(cds1.p, cds2.p, asPolygon, frac, angle, col, rawCoords)
            .py()
    }

    /// draws a triangle with the current drawing style. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
    #[pyo3(signature = (cds1, cds2, cds3, rawCoords = false))]
    fn DrawTriangle(
        &mut self,
        cds1: PyRef<'_, PyPoint2D>,
        cds2: PyRef<'_, PyPoint2D>,
        cds3: PyRef<'_, PyPoint2D>,
        rawCoords: bool,
    ) -> PyResult<()> {
        self.d.draw_triangle(cds1.p, cds2.p, cds3.p, rawCoords).py()
    }

    /// draws a polygon with the current drawing style. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
    #[pyo3(signature = (cds, rawCoords = false))]
    fn DrawPolygon(&mut self, cds: &Bound<'_, PyAny>, rawCoords: bool) -> PyResult<()> {
        let pts = cds
            .try_iter()?
            .map(|p| Ok(p?.extract::<PyRef<'_, PyPoint2D>>()?.p))
            .collect::<PyResult<Vec<_>>>()?;
        if pts.is_empty() {
            return Err(PyValueError::new_err("cds argument must be non-empty"));
        }
        self.d.draw_polygon(&pts, rawCoords).py()
    }

    /// draws a triangle with the current drawing style in the rectangle defined by the two points. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
    #[pyo3(signature = (cds1, cds2, rawCoords = false))]
    fn DrawEllipse(
        &mut self,
        cds1: PyRef<'_, PyPoint2D>,
        cds2: PyRef<'_, PyPoint2D>,
        rawCoords: bool,
    ) -> PyResult<()> {
        self.d.draw_ellipse(cds1.p, cds2.p, rawCoords).py()
    }

    /// draws a rectangle with the current drawing style in the rectangle defined by the two points. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
    #[pyo3(signature = (cds1, cds2, rawCoords = false))]
    fn DrawRect(
        &mut self,
        cds1: PyRef<'_, PyPoint2D>,
        cds2: PyRef<'_, PyPoint2D>,
        rawCoords: bool,
    ) -> PyResult<()> {
        self.d.draw_rect(cds1.p, cds2.p, rawCoords).py()
    }

    /// draws an arc with the current drawing style. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.The angles are in degrees; angle2 should be > angle1.
    #[pyo3(signature = (center, radius, angle1, angle2, rawCoords = false))]
    fn DrawArc(
        &mut self,
        center: PyRef<'_, PyPoint2D>,
        radius: f64,
        angle1: f64,
        angle2: f64,
        rawCoords: bool,
    ) -> PyResult<()> {
        self.d
            .draw_arc(center.p, radius, radius, angle1, angle2, rawCoords)
            .py()
    }

    /// draw a line indicating the presence of an attachment point (normally a squiggle line perpendicular to a bond). The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
    #[pyo3(signature = (cds1, cds2, color, len = 1.0, nSegments = 16, rawCoords = false))]
    fn DrawAttachmentLine(
        &mut self,
        cds1: PyRef<'_, PyPoint2D>,
        cds2: PyRef<'_, PyPoint2D>,
        color: &Bound<'_, PyAny>,
        len: f64,
        nSegments: u32,
        rawCoords: bool,
    ) -> PyResult<()> {
        self.d
            .draw_attachment_line(cds1.p, cds2.p, py_colour(color)?, len, nSegments, rawCoords)
            .py()
    }

    /// draw a line indicating the presence of an attachment point (normally a squiggle line perpendicular to a bond). The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
    #[pyo3(signature = (cds1, cds2, color1, color2, nSegments = 16, vertOffset = 0.05, rawCoords = false))]
    fn DrawWavyLine(
        &mut self,
        cds1: PyRef<'_, PyPoint2D>,
        cds2: PyRef<'_, PyPoint2D>,
        color1: &Bound<'_, PyAny>,
        color2: &Bound<'_, PyAny>,
        nSegments: u32,
        vertOffset: f64,
        rawCoords: bool,
    ) -> PyResult<()> {
        self.d
            .draw_wavy_line(
                cds1.p,
                cds2.p,
                py_colour(color1)?,
                py_colour(color2)?,
                nSegments,
                vertOffset,
                rawCoords,
            )
            .py()
    }

    /// add text to the canvas. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels. With an align argument (0 (=MIDDLE), 1 (=START), or 2 (=END)) the text is aligned.
    #[pyo3(signature = (string, pos, *args, **kwargs))]
    fn DrawString(
        &mut self,
        string: &str,
        pos: PyRef<'_, PyPoint2D>,
        args: &Bound<'_, PyTuple>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<()> {
        let aligned = kwargs.is_some_and(|k| k.contains("align").unwrap_or(false))
            || args
                .get_item(0)
                .is_ok_and(|a| !a.is_instance_of::<pyo3::types::PyBool>());
        if aligned {
            let a = bind_args(args, kwargs, &["align", "rawCoords"], "DrawString")?;
            let align = a[0]
                .as_ref()
                .ok_or_else(|| PyTypeError::new_err("DrawString() missing argument 'align'"))?
                .extract()?;
            let raw = a[1].as_ref().map_or(Ok(false), |o| o.extract())?;
            self.d
                .draw_string_aligned(string, pos.p, align_from_int(align)?, raw)
                .py()
        } else {
            let a = bind_args(args, kwargs, &["rawCoords"], "DrawString")?;
            let raw = a[0].as_ref().map_or(Ok(false), |o| o.extract())?;
            self.d.draw_string(string, pos.p, raw).py()
        }
    }

    /// get the coordinates in drawing space for a particular point in molecule space, or for a particular atom
    #[pyo3(signature = (point))]
    fn GetDrawCoords(&self, point: &Bound<'_, PyAny>) -> PyResult<PyPoint2D> {
        if let Ok(p) = point.extract::<PyRef<'_, PyPoint2D>>() {
            Ok(PyPoint2D {
                p: self.d.get_draw_coords(p.p).py()?,
            })
        } else {
            let idx: usize = point.extract()?;
            Ok(PyPoint2D {
                p: self.d.get_draw_coords_atom(idx).py()?,
            })
        }
    }

    /// clears the drawing by filling it with the background color
    fn ClearDrawing(&mut self) {
        self.d.clear_drawing();
    }

    /// Returns a modifiable version of the current drawing options
    fn drawOptions(slf: Bound<'_, Self>) -> PyMolDrawOptions {
        PyMolDrawOptions {
            owner: Some(slf.unbind()),
            own: MolDrawOptions::default(),
        }
    }

    /// Copies the drawing options passed in over our drawing options
    #[pyo3(signature = (opts))]
    fn SetDrawOptions(slf: &Bound<'_, Self>, opts: &Bound<'_, PyMolDrawOptions>) -> PyResult<()> {
        // the options may be this drawer's own, so copy them before borrowing
        // the drawer mutably
        let copy = opts.try_borrow()?.read(slf.py(), |o| o.clone())?;
        *slf.try_borrow_mut()?.d.draw_options_mut() = copy;
        Ok(())
    }
}

/// SVG molecule drawer
#[pyclass(
    name = "MolDraw2DSVG",
    module = "rdoxide.Chem.Draw.rdMolDraw2D",
    extends = PyMolDraw2D
)]
#[derive(Debug)]
pub struct PyMolDraw2DSVG;

#[pymethods]
#[allow(non_snake_case)]
impl PyMolDraw2DSVG {
    #[new]
    #[pyo3(signature = (width, height, panelWidth = -1, panelHeight = -1, noFreetype = false))]
    fn new(
        width: i32,
        height: i32,
        panelWidth: i32,
        panelHeight: i32,
        noFreetype: bool,
    ) -> PyResult<PyClassInitializer<Self>> {
        let d = MolDraw2DSVG::with_text(
            width,
            height,
            panelWidth,
            panelHeight,
            MolDrawOptions::default(),
            noFreetype,
        )
        .py()?;
        Ok(PyClassInitializer::from(PyMolDraw2D { d }).add_subclass(PyMolDraw2DSVG))
    }

    /// add the last bits of SVG to finish the drawing
    fn FinishDrawing(mut slf: PyRefMut<'_, Self>) {
        slf.as_super().d.finish_drawing();
    }

    /// add RDKit-specific information to the bottom of the drawing
    #[pyo3(signature = (mol, confId = -1))]
    fn AddMoleculeMetadata(
        mut slf: PyRefMut<'_, Self>,
        mol: PyRef<'_, PyMol>,
        confId: i32,
    ) -> PyResult<()> {
        slf.as_super()
            .d
            .add_molecule_metadata(&mol.mol, confId)
            .py()
    }

    /// allow atom selection in the SVG
    #[pyo3(signature = (mol, radius = 0.2, events = None))]
    fn TagAtoms(
        mut slf: PyRefMut<'_, Self>,
        mol: PyRef<'_, PyMol>,
        radius: f64,
        events: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        let mut map = BTreeMap::new();
        if let Some(e) = given(events)? {
            for (k, v) in e.cast::<PyDict>()?.iter() {
                map.insert(k.extract::<String>()?, v.extract::<String>()?);
            }
        }
        slf.as_super().d.tag_atoms(&mol.mol, radius, &map).py()
    }

    /// return the SVG
    fn GetDrawingText(mut slf: PyRefMut<'_, Self>) -> String {
        slf.as_super().d.get_drawing_text()
    }
}

// ----------------------------------------------------------------------
// functions
// ----------------------------------------------------------------------

/// Does some cleanup operations on the molecule to prepare it to draw nicely.
/// The operations include: kekulization, addition of chiral Hs (so that we can draw
/// wedges to them), wedging of bonds at chiral centers, and generation of a 2D
/// conformation if the molecule does not already have a conformation
///
/// Returns a modified copy of the molecule.
#[pyfunction]
#[pyo3(name = "PrepareMolForDrawing", signature = (mol, kekulize = true, addChiralHs = true, wedgeBonds = true, forceCoords = false, wavyBonds = false))]
#[allow(non_snake_case)]
fn prepare_mol_for_drawing_py(
    py: Python<'_>,
    mol: Option<PyRef<'_, PyMol>>,
    kekulize: bool,
    addChiralHs: bool,
    wedgeBonds: bool,
    forceCoords: bool,
    wavyBonds: bool,
) -> PyResult<Py<PyMol>> {
    let mol = mol.ok_or_else(|| {
        PyRuntimeError::new_err("Pre-condition Violation\n\tmolecule must not be None")
    })?;
    let mut copy: RWMol = mol.mol.clone();
    prepare_mol_for_drawing(
        &mut copy,
        kekulize,
        addChiralHs,
        wedgeBonds,
        forceCoords,
        wavyBonds,
    )
    .py()?;
    Py::new(py, PyMol { mol: copy })
}

/// Preps a molecule for drawing and actually draws it
#[pyfunction]
#[pyo3(name = "PrepareAndDrawMolecule", signature = (drawer, mol, legend = "", highlightAtoms = None, highlightBonds = None, highlightAtomColors = None, highlightBondColors = None, highlightAtomRadii = None, confId = -1, kekulize = true))]
#[allow(non_snake_case)]
fn prepare_and_draw_molecule_py(
    mut drawer: PyRefMut<'_, PyMolDraw2D>,
    mol: PyRef<'_, PyMol>,
    legend: &str,
    highlightAtoms: Option<&Bound<'_, PyAny>>,
    highlightBonds: Option<&Bound<'_, PyAny>>,
    highlightAtomColors: Option<&Bound<'_, PyAny>>,
    highlightBondColors: Option<&Bound<'_, PyAny>>,
    highlightAtomRadii: Option<&Bound<'_, PyAny>>,
    confId: i32,
    kekulize: bool,
) -> PyResult<()> {
    let hl = highlights(
        &mol.mol,
        highlightAtoms,
        highlightBonds,
        highlightAtomColors,
        highlightBondColors,
        highlightAtomRadii,
    )?;
    drawer
        .d
        .prepare_and_draw_molecule(
            &mol.mol, legend, hl, confId, kekulize, true, true, false, false,
        )
        .py()
}

/// Draws molecule in ACS 1996 mode.
#[pyfunction]
#[pyo3(name = "DrawMoleculeACS1996", signature = (drawer, mol, legend = "", highlightAtoms = None, highlightBonds = None, highlightAtomColors = None, highlightBondColors = None, highlightAtomRadii = None, confId = -1))]
#[allow(non_snake_case)]
fn draw_molecule_acs1996_py(
    mut drawer: PyRefMut<'_, PyMolDraw2D>,
    mol: PyRef<'_, PyMol>,
    legend: &str,
    highlightAtoms: Option<&Bound<'_, PyAny>>,
    highlightBonds: Option<&Bound<'_, PyAny>>,
    highlightAtomColors: Option<&Bound<'_, PyAny>>,
    highlightBondColors: Option<&Bound<'_, PyAny>>,
    highlightAtomRadii: Option<&Bound<'_, PyAny>>,
    confId: i32,
) -> PyResult<()> {
    let hl = highlights(
        &mol.mol,
        highlightAtoms,
        highlightBonds,
        highlightAtomColors,
        highlightBondColors,
        highlightAtomRadii,
    )?;
    drawer.d.draw_mol_acs1996(&mol.mol, legend, hl, confId).py()
}

/// Returns svg for a molecule
#[pyfunction]
#[pyo3(name = "MolToSVG", signature = (mol, width = 300, height = 300, highlightAtoms = None, kekulize = true, lineWidthMult = 1, includeAtomCircles = true, confId = -1))]
#[allow(non_snake_case)]
fn mol_to_svg(
    mol: PyRef<'_, PyMol>,
    width: i32,
    height: i32,
    highlightAtoms: Option<&Bound<'_, PyAny>>,
    kekulize: bool,
    lineWidthMult: u32,
    includeAtomCircles: bool,
    confId: i32,
) -> PyResult<String> {
    let _ = kekulize;
    let atoms = int_list(highlightAtoms, Some(mol.mol.num_atoms() as i32))?;
    let mut drawer = MolDraw2DSVG::new(width, height, -1, -1).py()?;
    let lw = drawer.line_width();
    drawer.set_line_width(lw * f64::from(lineWidthMult));
    drawer.draw_options_mut().circle_atoms = includeAtomCircles;
    drawer.draw_options_mut().prepare_mols_before_drawing = false;
    match atoms {
        Some(a) => drawer.draw_molecule_atoms(&mol.mol, "", &a, None, None, confId),
        None => drawer.draw_molecule(&mol.mol, "", confId),
    }
    .py()?;
    drawer.finish_drawing();
    Ok(drawer.get_drawing_text())
}

/// Returns ACS 1996 mode svg for a molecule
#[pyfunction]
#[pyo3(name = "MolToACS1996SVG", signature = (mol, legend = "", highlightAtoms = None, highlightBonds = None, highlightAtomColors = None, highlightBondColors = None, highlightAtomRadii = None, confId = -1))]
#[allow(non_snake_case)]
fn mol_to_acs1996_svg(
    mol: PyRef<'_, PyMol>,
    legend: &str,
    highlightAtoms: Option<&Bound<'_, PyAny>>,
    highlightBonds: Option<&Bound<'_, PyAny>>,
    highlightAtomColors: Option<&Bound<'_, PyAny>>,
    highlightBondColors: Option<&Bound<'_, PyAny>>,
    highlightAtomRadii: Option<&Bound<'_, PyAny>>,
    confId: i32,
) -> PyResult<String> {
    let hl = highlights(
        &mol.mol,
        highlightAtoms,
        highlightBonds,
        highlightAtomColors,
        highlightBondColors,
        highlightAtomRadii,
    )?;
    let mut drawer = MolDraw2DSVG::new(-1, -1, -1, -1).py()?;
    drawer.draw_mol_acs1996(&mol.mol, legend, hl, confId).py()?;
    drawer.finish_drawing();
    Ok(drawer.get_drawing_text())
}

/// Set the draw options to produce something as close as possible to
/// the ACS 1996 guidelines as described at
/// https://en.wikipedia.org/wiki/Wikipedia:Manual_of_Style/Chemistry/Structure_drawing
///
///  - MolDrawOptions opt - the options what will be changed
///  - float meanBondLength - mean bond length of the molecule
///
///  Works best if the MolDraw2D object is created with width and height -1 (a
///  flexiCanvas).
///  The mean bond length may be calculated with MeanBondLength.
///  It is used to calculate the offset for the lines in multiple bonds.
///
///  Options changed are:
///    bondLineWidth = 0.6
///    scaleBondWidth = false
///    scalingFactor = 14.4 / meanBondLen
///    multipleBondOffset = 0.18
///    highlightBondWidthMultiplier = 32
///    setMonochromeMode - black and white
///    fixedFontSize = 10
///    additionalAtomLabelPadding = 0.066
///    fontFile - if it isn't set already, then if RDBASE is set and the file
///               exists, uses $RDBASE/Data/Fonts/FreeSans.ttf.  Otherwise uses
///               BuiltinRobotoRegular.
#[pyfunction]
#[pyo3(name = "SetACS1996Mode", signature = (drawOptions, meanBondLength))]
#[allow(non_snake_case)]
fn set_acs1996_mode(
    py: Python<'_>,
    mut drawOptions: PyRefMut<'_, PyMolDrawOptions>,
    meanBondLength: f64,
) -> PyResult<()> {
    drawOptions
        .write(py, |o| set_acs1996_options(o, meanBondLength))?
        .py()
}

/// Calculate the mean bond length for the molecule.
#[pyfunction]
#[pyo3(name = "MeanBondLength", signature = (mol, confId = -1))]
#[allow(non_snake_case)]
fn mean_bond_length_py(mol: PyRef<'_, PyMol>, confId: i32) -> PyResult<f64> {
    mean_bond_length(&mol.mol, confId).py()
}

/// set dark mode for a MolDrawOptions or MolDraw2D object
#[pyfunction]
#[pyo3(name = "SetDarkMode", signature = (d2d))]
fn set_dark_mode(py: Python<'_>, d2d: &Bound<'_, PyAny>) -> PyResult<()> {
    if let Ok(mut o) = d2d.extract::<PyRefMut<'_, PyMolDrawOptions>>() {
        o.write(py, |o| o.set_dark_mode())
    } else {
        let mut d: PyRefMut<'_, PyMolDraw2D> = d2d.extract()?;
        d.d.draw_options_mut().set_dark_mode();
        Ok(())
    }
}

/// set monochrome mode for a MolDrawOptions or MolDraw2D object
#[pyfunction]
#[pyo3(name = "SetMonochromeMode", signature = (options, fgColour, bgColour))]
#[allow(non_snake_case)]
fn set_monochrome_mode(
    py: Python<'_>,
    options: &Bound<'_, PyAny>,
    fgColour: &Bound<'_, PyAny>,
    bgColour: &Bound<'_, PyAny>,
) -> PyResult<()> {
    let fg = py_colour(fgColour)?;
    let bg = py_colour(bgColour)?;
    if let Ok(mut o) = options.extract::<PyRefMut<'_, PyMolDrawOptions>>() {
        o.write(py, |o| o.set_monochrome_mode(fg, bg))
    } else {
        let mut d: PyRefMut<'_, PyMolDraw2D> = options.extract()?;
        d.d.draw_options_mut().set_monochrome_mode(fg, bg);
        Ok(())
    }
}

/// Updates a MolDrawOptions object from a JSON string
#[pyfunction]
#[pyo3(name = "UpdateMolDrawOptionsFromJSON", signature = (opts, json))]
fn update_mol_draw_options_from_json_py(
    py: Python<'_>,
    mut opts: PyRefMut<'_, PyMolDrawOptions>,
    json: &str,
) -> PyResult<()> {
    opts.write(py, |o| update_mol_draw_options_from_json(o, json))?
        .py()
}

/// Updates the drawing options of a MolDraw2D object from a JSON string
#[pyfunction]
#[pyo3(name = "UpdateDrawerParamsFromJSON", signature = (drawer, json))]
fn update_drawer_params_from_json_py(
    mut drawer: PyRefMut<'_, PyMolDraw2D>,
    json: &str,
) -> PyResult<()> {
    update_mol_draw_options_from_json(drawer.d.draw_options_mut(), json).py()
}

// ----------------------------------------------------------------------
// contouring
// ----------------------------------------------------------------------

/// Parameters for drawing contours
#[pyclass(
    name = "ContourParams",
    module = "rdoxide.Chem.Draw.rdMolDraw2D",
    skip_from_py_object
)]
#[derive(Debug, Clone, Default)]
pub struct PyContourParams {
    p: ContourParams,
}

macro_rules! contour_props {
    ($( ($py:literal, $get:ident, $set:ident, $field:ident, $ty:ty) ),* $(,)?) => {
        #[pymethods]
        impl PyContourParams {
            $(
                #[getter($py)]
                fn $get(&self) -> $ty {
                    self.p.$field
                }
                #[setter($py)]
                fn $set(&mut self, value: $ty) {
                    self.p.$field = value;
                }
            )*
        }
    };
}

contour_props!(
    ("setScale", get_set_scale, set_set_scale, set_scale, bool),
    (
        "dashNegative",
        get_dash_negative,
        set_dash_negative,
        dash_negative,
        bool
    ),
    ("fillGrid", get_fill_grid, set_fill_grid, fill_grid, bool),
    (
        "gridResolution",
        get_grid_resolution,
        set_grid_resolution,
        grid_resolution,
        f64
    ),
    (
        "contourWidth",
        get_contour_width,
        set_contour_width,
        contour_width,
        f64
    ),
    (
        "extraGridPadding",
        get_extra_grid_padding,
        set_extra_grid_padding,
        extra_grid_padding,
        f64
    ),
    (
        "drawAsLines",
        get_draw_as_lines,
        set_draw_as_lines,
        draw_as_lines,
        bool
    ),
    (
        "coordScaleForQuantization",
        get_coord_scale,
        set_coord_scale,
        coord_scale_for_quantization,
        f64
    ),
    (
        "isovalScaleForQuantization",
        get_isoval_scale,
        set_isoval_scale,
        isoval_scale_for_quantization,
        f64
    ),
    (
        "useFillThreshold",
        get_use_fill_threshold,
        set_use_fill_threshold,
        use_fill_threshold,
        bool
    ),
    (
        "fillThreshold",
        get_fill_threshold,
        set_fill_threshold,
        fill_threshold,
        f64
    ),
    (
        "fillThresholdIsFraction",
        get_fill_threshold_is_fraction,
        set_fill_threshold_is_fraction,
        fill_threshold_is_fraction,
        bool
    ),
);

#[pymethods]
#[allow(non_snake_case)]
impl PyContourParams {
    #[new]
    fn new() -> Self {
        PyContourParams::default()
    }

    /// the color map to use when filling the grid
    #[getter(colourMap)]
    fn get_colour_map<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let items = self
            .p
            .colour_map
            .iter()
            .map(|c| colour_tuple(py, *c))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, items)
    }

    #[setter(colourMap)]
    fn set_colour_map(&mut self, colours: &Bound<'_, PyAny>) -> PyResult<()> {
        self.p.colour_map = colours
            .try_iter()?
            .map(|c| py_colour(&c?))
            .collect::<PyResult<_>>()?;
        Ok(())
    }

    #[pyo3(signature = (colours))]
    fn setColourMap(&mut self, colours: &Bound<'_, PyAny>) -> PyResult<()> {
        self.set_colour_map(colours)
    }

    /// the color to use for drawing the contours
    #[getter(contourColour)]
    fn get_contour_colour<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        colour_tuple(py, self.p.contour_colour)
    }

    #[setter(contourColour)]
    fn set_contour_colour(&mut self, tpl: &Bound<'_, PyAny>) -> PyResult<()> {
        self.p.contour_colour = py_colour(tpl)?;
        Ok(())
    }

    #[pyo3(signature = (colour))]
    fn setContourColour(&mut self, colour: &Bound<'_, PyAny>) -> PyResult<()> {
        self.set_contour_colour(colour)
    }
}

fn float_list(obj: Option<&Bound<'_, PyAny>>, name: &str) -> PyResult<Vec<f64>> {
    match given(obj)? {
        Some(o) => o.try_iter()?.map(|v| v?.extract()).collect(),
        None => Err(PyValueError::new_err(format!(
            "{name} argument must be non-empty"
        ))),
    }
}

fn optional_mol<'a>(obj: Option<&'a Bound<'_, PyAny>>) -> PyResult<Option<PyRef<'a, PyMol>>> {
    match obj {
        Some(m) if !m.is_none() => Ok(Some(mol_arg(m)?)),
        _ => Ok(None),
    }
}

/// Generates and draws contours for a set of gaussians
///
///   - drawer: the MolDraw2D object to use
///   - locs: locations of the gaussians
///   - heights: the heights (or weights) of the gaussians
///   - widths: the standard deviations of the gaussians
///   - nContours: the number of contours to draw
///   - levels: the contours to use
///   - ps: additional parameters controlling the contouring.
///   - mol: molecule used to help set scale.
///
///   The values are calculated on a grid with spacing params.gridResolution.
///   If params.setScale  is set, the grid size will be calculated based on the
///   locations of the gaussians and params.extraGridPadding. Otherwise the current
///   size of the viewport will be used.
///
///   If the levels argument is empty, the contour levels will be determined
///   automatically from the max and min values on the grid and levels will
///   be updated to include the contour levels.
///
///   If params.fillGrid is set, the data on the grid will also be drawn using
///   the color scheme in params.colourMap
///
///   If mol is not 0, uses the molecule to help set the scale, assuming that
///   it will be drawn over the plot, so needs to fit on it.
#[pyfunction]
#[pyo3(name = "ContourAndDrawGaussians", signature = (drawer, locs, heights, widths, nContours = 10, levels = None, params = None, mol = None))]
#[allow(non_snake_case)]
fn contour_and_draw_gaussians_py(
    mut drawer: PyRefMut<'_, PyMolDraw2D>,
    locs: &Bound<'_, PyAny>,
    heights: &Bound<'_, PyAny>,
    widths: &Bound<'_, PyAny>,
    nContours: usize,
    levels: Option<&Bound<'_, PyAny>>,
    params: Option<PyRef<'_, PyContourParams>>,
    mol: Option<&Bound<'_, PyAny>>,
) -> PyResult<()> {
    let locs = match given(Some(locs))? {
        Some(l) => l
            .try_iter()?
            .map(|p| Ok(p?.extract::<PyRef<'_, PyPoint2D>>()?.p))
            .collect::<PyResult<Vec<_>>>()?,
        None => return Err(PyValueError::new_err("locs argument must be non-empty")),
    };
    let heights = float_list(Some(heights), "heights")?;
    let widths = float_list(Some(widths), "widths")?;
    let mut levels = match given(levels)? {
        Some(l) => float_list(Some(l), "levels")?,
        None => Vec::new(),
    };
    let p = params.map_or_else(ContourParams::default, |p| p.p.clone());
    let mol = optional_mol(mol)?;
    contour_and_draw_gaussians(
        &mut drawer.d,
        &locs,
        &heights,
        &widths,
        nContours,
        &mut levels,
        &p,
        mol.as_ref().map(|m| &m.mol),
    )
    .py()
}

/// Generates and draws contours for data on a grid
///
///   - drawer: the MolDraw2D object to use
///   - data: numpy array with the data to be contoured
///   - xcoords: the x coordinates of the grid
///   - ycoords: the y coordinates of the grid
///   - nContours: the number of contours to draw
///   - levels: the contours to use
///   - ps: additional parameters controlling the contouring
///   - mol: molecule used to help set scale.
///
///   The values are calculated on a grid with spacing params.gridResolution.
///   If params.setScale  is set, the grid size will be calculated based on the
///   locations of the gaussians and params.extraGridPadding. Otherwise the current
///   size of the viewport will be used.
///
///   If the levels argument is empty, the contour levels will be determined
///   automatically from the max and min values on the grid and levels will
///   be updated to include the contour levels.
///
///   If params.fillGrid is set, the data on the grid will also be drawn using
///   the color scheme in params.colourMap
///
///   If mol is not 0, uses the molecule to help set the scale, assuming that
///   it will be drawn over the plot, so needs to fit on it.
#[pyfunction]
#[pyo3(name = "ContourAndDrawGrid", signature = (drawer, data, xcoords, ycoords, nContours = 10, levels = None, params = None, mol = None))]
#[allow(non_snake_case)]
fn contour_and_draw_grid_py(
    py: Python<'_>,
    mut drawer: PyRefMut<'_, PyMolDraw2D>,
    data: &Bound<'_, PyAny>,
    xcoords: &Bound<'_, PyAny>,
    ycoords: &Bound<'_, PyAny>,
    nContours: usize,
    levels: Option<&Bound<'_, PyAny>>,
    params: Option<PyRef<'_, PyContourParams>>,
    mol: Option<&Bound<'_, PyAny>>,
) -> PyResult<()> {
    let np = py.import("numpy")?;
    if !data.is_instance(&np.getattr("ndarray")?)? {
        return Err(PyValueError::new_err("data argument must be a numpy array"));
    }
    let arr = np
        .call_method1("ascontiguousarray", (data, "float64"))
        .map_err(|_| PyValueError::new_err("could not convert data argument"))?;
    let shape: Vec<usize> = arr.getattr("shape")?.extract()?;
    let &[n_x, n_y] = shape.as_slice() else {
        return Err(PyValueError::new_err("could not convert data argument"));
    };
    let values = crate::pyconv::float64_flat(&arr)?;
    let xs = float_list(Some(xcoords), "xcoords")?;
    let ys = float_list(Some(ycoords), "ycoords")?;
    let mut levels = match given(levels)? {
        Some(l) => float_list(Some(l), "levels")?,
        None => Vec::new(),
    };
    if n_x != xs.len() {
        return Err(PyValueError::new_err(
            "data array and xcoords sizes do not match.\nDid you forget to call np.transpose() on the array?",
        ));
    }
    if n_y != ys.len() {
        return Err(PyValueError::new_err(
            "data array and ycoords sizes do not match",
        ));
    }
    let p = params.map_or_else(ContourParams::default, |p| p.p.clone());
    let mol = optional_mol(mol)?;
    contour_and_draw_grid(
        &mut drawer.d,
        &values,
        &xs,
        &ys,
        nContours,
        &mut levels,
        &p,
        mol.as_ref().map(|m| &m.mol),
    )
    .py()
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing a C++ implementation of 2D molecule drawing",
    )?;
    add_enum_exported(
        m,
        "MultiColourHighlightStyle",
        &[("CircleAndLine", 0), ("Lasso", 1)],
        None,
    )?;
    add_enum_exported(
        m,
        "LegendPosition",
        &[("Bottom", 0), ("Top", 1), ("Left", 2), ("Right", 3)],
        None,
    )?;
    add_enum_exported(
        m,
        "DrawElement",
        &[
            ("NONE", 0),
            ("PRESHAPES", 1),
            ("BONDS", 2),
            ("ATOMLABELS", 4),
            ("HIGHLIGHTS", 8),
            ("ANNOTATIONS", 16),
            ("RADICALS", 32),
            ("POSTSHAPES", 64),
            ("ALL", 0x7fff_ffff),
        ],
        None,
    )?;
    m.add_class::<PyIntStringMap>()?;
    m.add_class::<PyMolDrawOptions>()?;
    m.add_class::<PyMolDraw2D>()?;
    m.add_class::<PyMolDraw2DSVG>()?;
    m.add_class::<PyContourParams>()?;
    m.add_function(wrap_pyfunction!(contour_and_draw_gaussians_py, m)?)?;
    m.add_function(wrap_pyfunction!(contour_and_draw_grid_py, m)?)?;
    m.add_function(wrap_pyfunction!(prepare_mol_for_drawing_py, m)?)?;
    m.add_function(wrap_pyfunction!(prepare_and_draw_molecule_py, m)?)?;
    m.add_function(wrap_pyfunction!(draw_molecule_acs1996_py, m)?)?;
    m.add_function(wrap_pyfunction!(mol_to_svg, m)?)?;
    m.add_function(wrap_pyfunction!(mol_to_acs1996_svg, m)?)?;
    m.add_function(wrap_pyfunction!(set_acs1996_mode, m)?)?;
    m.add_function(wrap_pyfunction!(mean_bond_length_py, m)?)?;
    m.add_function(wrap_pyfunction!(set_dark_mode, m)?)?;
    m.add_function(wrap_pyfunction!(update_mol_draw_options_from_json_py, m)?)?;
    m.add_function(wrap_pyfunction!(update_drawer_params_from_json_py, m)?)?;
    m.add_function(wrap_pyfunction!(set_monochrome_mode, m)?)?;
    Ok(())
}
