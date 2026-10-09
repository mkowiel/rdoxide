//! `rdkit.Chem.rdDetermineBonds` (upstream
//! `GraphMol/DetermineBonds/Wrap/rdDetermineBonds.cpp`).

use pyo3::prelude::*;
use rdoxide::graphmol::determine_bonds::{
    HUECKEL_SUPPORT, determine_bond_orders, determine_bonds, determine_connectivity,
};

use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Assigns atomic connectivity to a molecule using atomic coordinates,
/// disregarding pre-existing bonds
///
/// Args:
///    mol : the molecule of interest; it must have a 3D conformer
///    useHueckel : (optional) if this is  \c true, extended Hueckel theory
///        will be used to determine connectivity rather than the van der Waals
///        or connect-the-dots methods
///    charge : (optional) the charge of the molecule; it must be provided if
///        the Hueckel method is used and charge is non-zero
///    covFactor : (optional) the factor with which to multiply each covalent
///        radius if the van der Waals method is used
///    useVdw: (optional) if this is false, the connect-the-dots method
///        will be used instead of the van der Waals method
#[pyfunction]
#[pyo3(name = "DetermineConnectivity", signature = (mol, useHueckel = false, charge = 0, covFactor = 1.3, useVdw = false))]
#[allow(non_snake_case)]
fn py_determine_connectivity(
    mut mol: PyRefMut<'_, PyMol>,
    useHueckel: bool,
    charge: i32,
    covFactor: f64,
    useVdw: bool,
) -> PyResult<()> {
    determine_connectivity(&mut mol.mol, useHueckel, charge, covFactor, useVdw).py()
}

/// Assigns bond orders to a molecule with known connectivity, using atomic coordinates
///
/// Args:
///    mol : the molecule of interest; it must have a 3D conformer
///    charge : (optional) the charge of the molecule
///    allowChargedFragments : (optional) if this is true, formal charges
///        will be placed on atoms according to their valency; otherwise, radical
///        electrons will be placed on the atoms
///    embedChiral : (optional) if this is true,
///       chirality information will be embedded into the molecule; the function calls
///       sanitizeMol() when this is true
///    useAtomMap : (optional) if this is true, an atom map will be created for the
///       molecule
///    maxIterations: (optional) maximum number of iterations to run in the bond order
///    determination algorithm, after which a MaxFindBondOrdersItersExceeded
///    exception will be thrown. Defaults to 0 (no limit)
#[pyfunction]
#[pyo3(name = "DetermineBondOrders", signature = (mol, charge = 0, allowChargedFragments = true, embedChiral = true, useAtomMap = false, maxIterations = 0))]
#[allow(non_snake_case)]
fn py_determine_bond_orders(
    py: Python<'_>,
    mut mol: PyRefMut<'_, PyMol>,
    charge: i32,
    allowChargedFragments: bool,
    embedChiral: bool,
    useAtomMap: bool,
    maxIterations: usize,
) -> PyResult<()> {
    let m = &mut mol.mol;
    py.detach(|| {
        determine_bond_orders(
            m,
            charge,
            allowChargedFragments,
            embedChiral,
            useAtomMap,
            maxIterations,
        )
    })
    .py()
}

/// Assigns atomic connectivity and bond orders to a molecule using atomic
/// coordinates, disregarding pre-existing bonds; it is recommended to sanitize
/// the molecule after calling this function if embedChiral is not set to true
///
/// Args:
///    mol : the molecule of interest; it must have a 3D conformer
///    useHueckel : (optional) if this is  \c true, extended Hueckel theory
///        will be used to determine connectivity rather than the van der Waals
///        or connect-the-dots methods
///    charge : (optional) the charge of the molecule; it must be provided if
///        the Hueckel method is used and charge is non-zero
///    covFactor : (optional) the factor with which to multiply each covalent
///        radius if the van der Waals method is used
///    allowChargedFragments : (optional) if this is true, formal charges
///        will be placed on atoms according to their valency; otherwise, radical
///        electrons will be placed on the atoms
///    embedChiral : (optional) if this is true,
///       chirality information will be embedded into the molecule; the function calls
///       sanitizeMol() when this is true
///    useAtomMap : (optional) if this is true, an atom map will be created for the
///       molecule
///    useVdw: (optional) if this is false, the connect-the-dots method
///        will be used instead of the van der Waals method
///    maxIterations: (optional) maximum number of iterations to run in the bond order
///    determination algorithm, after which a MaxFindBondOrdersItersExceeded
///    exception will be thrown. Defaults to 0 (no limit)
#[pyfunction]
#[pyo3(name = "DetermineBonds", signature = (mol, useHueckel = false, charge = 0, covFactor = 1.3, allowChargedFragments = true, embedChiral = true, useAtomMap = false, useVdw = false, maxIterations = 0))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn py_determine_bonds(
    py: Python<'_>,
    mut mol: PyRefMut<'_, PyMol>,
    useHueckel: bool,
    charge: i32,
    covFactor: f64,
    allowChargedFragments: bool,
    embedChiral: bool,
    useAtomMap: bool,
    useVdw: bool,
    maxIterations: usize,
) -> PyResult<()> {
    let m = &mut mol.mol;
    py.detach(|| {
        determine_bonds(
            m,
            useHueckel,
            charge,
            covFactor,
            allowChargedFragments,
            embedChiral,
            useAtomMap,
            useVdw,
            maxIterations,
        )
    })
    .py()
}

/// Whether extended Hueckel connectivity (YAeHMOP) support is available.
#[pyfunction]
#[pyo3(name = "hueckelEnabled")]
fn py_hueckel_enabled() -> bool {
    HUECKEL_SUPPORT
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing a C++ implementation of the xyz2mol algorithm. This is based on xyz2mol: https://github.com/jensengroup/xyz2mol",
    )?;
    m.add_function(wrap_pyfunction!(py_determine_connectivity, m)?)?;
    m.add_function(wrap_pyfunction!(py_determine_bond_orders, m)?)?;
    m.add_function(wrap_pyfunction!(py_determine_bonds, m)?)?;
    m.add_function(wrap_pyfunction!(py_hueckel_enabled, m)?)?;
    Ok(())
}
