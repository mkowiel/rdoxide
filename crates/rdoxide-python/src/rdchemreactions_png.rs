//! Reactions in PNG metadata in `rdoxide.Chem.rdChemReactions` (upstream
//! `ChemReactions/Wrap/rdChemReactions.cpp`, the `RDK_USE_BOOST_IOSTREAMS`
//! block, and `PNGParser.cpp`).

use pyo3::prelude::*;
use pyo3::types::PyBytes;
use rdoxide::graphmol::chem_reactions::png::{
    ReactionPngOptions, add_chemical_reaction_to_png, png_to_chemical_reaction,
};

use crate::errors::IntoPyResult;
use crate::pyconv::text_bytes;
use crate::rdchemreactions::{PyChemicalReaction, synced};
use crate::rdmolfiles_png::{read_input, read_unchecked};

/// construct a ChemicalReaction from metadata in a PNG file
#[pyfunction]
#[pyo3(name = "ReactionFromPNGFile", signature = (fname))]
fn reaction_from_png_file(py: Python<'_>, fname: &str) -> PyResult<PyChemicalReaction> {
    let data = read_input(fname)?;
    let rxn = py.detach(|| png_to_chemical_reaction(&data)).py()?;
    Ok(PyChemicalReaction::new(rxn))
}

/// construct a ChemicalReaction from an string with PNG data
#[pyfunction]
#[pyo3(name = "ReactionFromPNGString", signature = (data))]
fn reaction_from_png_string(
    py: Python<'_>,
    data: &Bound<'_, PyAny>,
) -> PyResult<PyChemicalReaction> {
    let data = text_bytes(data)?;
    let rxn = py.detach(|| png_to_chemical_reaction(&data)).py()?;
    Ok(PyChemicalReaction::new(rxn))
}

fn add_to_png<'py>(
    py: Python<'py>,
    reaction: &Bound<'py, PyChemicalReaction>,
    data: Vec<u8>,
    opts: ReactionPngOptions,
) -> PyResult<Bound<'py, PyBytes>> {
    let rxn = synced(reaction).rxn.clone();
    let out = py
        .detach(|| add_chemical_reaction_to_png(&rxn, &data, &opts))
        .py()?;
    Ok(PyBytes::new(py, &out))
}

/// Reads the contents of a PNG file and adds metadata about a reaction to it.
/// The modified file contents are returned.
#[pyfunction]
#[pyo3(name = "ReactionMetadataToPNGFile", signature = (mol, filename, includePkl = true, includeSmiles = true, includeSmarts = false, includeMol = false))]
#[allow(non_snake_case)]
fn reaction_metadata_to_png_file<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyChemicalReaction>,
    filename: &str,
    includePkl: bool,
    includeSmiles: bool,
    includeSmarts: bool,
    includeMol: bool,
) -> PyResult<Bound<'py, PyBytes>> {
    let opts = ReactionPngOptions {
        include_pkl: includePkl,
        include_smiles: includeSmiles,
        include_smarts: includeSmarts,
        include_rxn: includeMol,
    };
    add_to_png(py, mol, read_unchecked(filename), opts)
}

/// Adds metadata about a reaction to the PNG string passed in.
/// The modified string is returned.
#[pyfunction]
#[pyo3(name = "ReactionMetadataToPNGString", signature = (mol, pngdata, includePkl = true, includeSmiles = true, includeSmarts = false, includeRxn = false))]
#[allow(non_snake_case)]
fn reaction_metadata_to_png_string<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyChemicalReaction>,
    pngdata: &Bound<'py, PyAny>,
    includePkl: bool,
    includeSmiles: bool,
    includeSmarts: bool,
    includeRxn: bool,
) -> PyResult<Bound<'py, PyBytes>> {
    let opts = ReactionPngOptions {
        include_pkl: includePkl,
        include_smiles: includeSmiles,
        include_smarts: includeSmarts,
        include_rxn: includeRxn,
    };
    add_to_png(py, mol, text_bytes(pngdata)?, opts)
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(reaction_from_png_file, m)?)?;
    m.add_function(wrap_pyfunction!(reaction_from_png_string, m)?)?;
    m.add_function(wrap_pyfunction!(reaction_metadata_to_png_file, m)?)?;
    m.add_function(wrap_pyfunction!(reaction_metadata_to_png_string, m)?)?;
    Ok(())
}
