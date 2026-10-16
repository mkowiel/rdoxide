//! `AdjustQueryProperties`, `AdjustQueryParameters` and `AdjustQueryWhichFlags`
//! of `rdoxide.Chem.rdmolops` (upstream `GraphMol/Wrap/MolOps.cpp`).

use pyo3::prelude::*;
use rdoxide::graphmol::adjust_query::{AdjustQueryParameters, adjust_query_properties};
use rdoxide::graphmol::generic_groups::{
    adjust_query_properties_with_generic_groups, convert_generic_queries_to_substance_groups,
    set_generic_queries_from_properties,
};

use crate::enums::add_enum_exported;
use crate::errors::{IntoPyResult, argument_error, boost_type_name, no_converter};
use crate::mol::PyMol;

const WHICH_FLAGS: [(&str, i64); 7] = [
    ("ADJUST_IGNORENONE", 0x0),
    ("ADJUST_IGNORECHAINS", 0x1),
    ("ADJUST_IGNORERINGS", 0x4),
    ("ADJUST_IGNOREDUMMIES", 0x2),
    ("ADJUST_IGNORENONDUMMIES", 0x8),
    ("ADJUST_IGNOREMAPPED", 0x10),
    ("ADJUST_IGNOREALL", 0xFFF_FFFF),
];

/// Parameters controlling which components of the query atoms/bonds are adjusted.
///
/// Note that some of the options here are either directly contradictory or make
///   no sense when combined with each other. We generally assume that client code
///   is doing something sensible and don't attempt to detect possible conflicts or
///   problems.
///
/// A note on the flags controlling which atoms/bonds are modified:
///    These generally limit the set of atoms/bonds to be modified.
///    For example:
///        - ADJUST_IGNORERINGS atoms/bonds in rings will not be modified.
///        - ADJUST_IGNORENONE causes all atoms/bonds to be modified
///        - ADJUST_IGNOREALL no atoms/bonds will be modified
///    Some of the options obviously make no sense for bonds
#[pyclass(
    name = "AdjustQueryParameters",
    module = "rdoxide.Chem.rdmolops",
    dict,
    skip_from_py_object
)]
#[derive(Clone, Default)]
pub struct PyAdjustQueryParameters {
    pub(crate) params: AdjustQueryParameters,
}

/// A read/write attribute of the parameters, as `def_readwrite` makes it.
macro_rules! attributes {
    ($($(#[$m:meta])* $py:literal, $get:ident, $set:ident, $field:ident: $ty:ty;)*) => {
        #[pymethods]
        impl PyAdjustQueryParameters {
            $(
                $(#[$m])*
                #[getter($py)]
                fn $get(&self) -> $ty {
                    self.params.$field
                }
                #[setter($py)]
                fn $set(&mut self, value: $ty) {
                    self.params.$field = value;
                }
            )*
        }
    };
}

attributes! {
    /// add degree queries
    "adjustDegree", get_adjust_degree, set_adjust_degree, adjust_degree: bool;
    /// controls which atoms have their degree queries changed
    "adjustDegreeFlags", get_adjust_degree_flags, set_adjust_degree_flags, adjust_degree_flags: u32;
    /// adjust the heavy-atom degree
    "adjustHeavyDegree", get_adjust_heavy_degree, set_adjust_heavy_degree, adjust_heavy_degree: bool;
    /// controls which atoms have their heavy-atom degree queries changed
    "adjustHeavyDegreeFlags", get_adjust_heavy_degree_flags, set_adjust_heavy_degree_flags, adjust_heavy_degree_flags: u32;
    /// add ring-count queries
    "adjustRingCount", get_adjust_ring_count, set_adjust_ring_count, adjust_ring_count: bool;
    /// controls which atoms have ring-count queries added
    "adjustRingCountFlags", get_adjust_ring_count_flags, set_adjust_ring_count_flags, adjust_ring_count_flags: u32;
    /// convert dummy atoms without isotope labels to any-atom queries
    "makeDummiesQueries", get_make_dummies_queries, set_make_dummies_queries, make_dummies_queries: bool;
    /// perceive and set aromaticity
    "aromatizeIfPossible", get_aromatize_if_possible, set_aromatize_if_possible, aromatize_if_possible: bool;
    /// converts bonds to generic queries (any bonds)
    "makeBondsGeneric", get_make_bonds_generic, set_make_bonds_generic, make_bonds_generic: bool;
    /// controls which bonds are converted to generic queries
    "makeBondsGenericFlags", get_make_bonds_generic_flags, set_make_bonds_generic_flags, make_bonds_generic_flags: u32;
    /// convert atoms to generic queries (any atoms)
    "makeAtomsGeneric", get_make_atoms_generic, set_make_atoms_generic, make_atoms_generic: bool;
    /// controls which atoms are converted to generic queries
    "makeAtomsGenericFlags", get_make_atoms_generic_flags, set_make_atoms_generic_flags, make_atoms_generic_flags: u32;
    /// add ring-chain queries to atoms
    "adjustRingChain", get_adjust_ring_chain, set_adjust_ring_chain, adjust_ring_chain: bool;
    /// controls which atoms have ring-chain queries added
    "adjustRingChainFlags", get_adjust_ring_chain_flags, set_adjust_ring_chain_flags, adjust_ring_chain_flags: u32;
    /// if this is set sterochemistry information will be removed from double bonds that do not have the stereoCare property set
    "useStereoCareForBonds", get_use_stereo_care_for_bonds, set_use_stereo_care_for_bonds, use_stereo_care_for_bonds: bool;
    /// set bond queries in conjugated five-rings to SINGLE|DOUBLE|AROMATIC
    "adjustConjugatedFiveRings", get_adjust_conjugated_five_rings, set_adjust_conjugated_five_rings, adjust_conjugated_five_rings: bool;
    /// uses the 5-ring aromaticity behavior of the (former) MDL software as documented in the Chemical Representation Guide
    "setMDLFiveRingAromaticity", get_set_mdl_five_ring_aromaticity, set_set_mdl_five_ring_aromaticity, set_mdl_five_ring_aromaticity: bool;
    /// set single bonds bewteen aromatic or conjugated atoms and degree-one neighbors to SINGLE|AROMATIC
    "adjustSingleBondsToDegreeOneNeighbors", get_adjust_single_bonds_to_degree_one, set_adjust_single_bonds_to_degree_one, adjust_single_bonds_to_degree_one_neighbors: bool;
    /// sets non-ring single bonds between two aromatic or conjugated atoms to SINGLE|AROMATIC
    "adjustSingleBondsBetweenAromaticAtoms", get_adjust_single_bonds_between_aromatic, set_adjust_single_bonds_between_aromatic, adjust_single_bonds_between_aromatic_atoms: bool;
}

#[pymethods]
impl PyAdjustQueryParameters {
    #[new]
    fn py_new() -> Self {
        Self::default()
    }

    /// Returns an AdjustQueryParameters object with all parameters set to false
    #[staticmethod]
    #[pyo3(name = "NoAdjustments")]
    fn no_adjustments() -> Self {
        PyAdjustQueryParameters {
            params: AdjustQueryParameters::no_adjustments(),
        }
    }
}

/// Returns a new molecule where the query properties of atoms have been modified.
#[pyfunction]
#[pyo3(name = "AdjustQueryProperties", signature = (mol, params = None))]
fn py_adjust_query_properties(
    py: Python<'_>,
    mol: &Bound<'_, PyAny>,
    params: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyMol> {
    let bad_call = || {
        let mut got = vec![boost_type_name(mol)];
        got.extend(params.map(boost_type_name));
        argument_error(
            py,
            "AdjustQueryProperties",
            &got,
            "AdjustQueryProperties(RDKit::ROMol mol, boost::python::api::object params=None)",
        )
    };
    let mol = mol.cast::<PyMol>().map_err(|_| bad_call())?;
    let params = match params.filter(|p| !p.is_none()) {
        None => AdjustQueryParameters::default(),
        Some(p) => p
            .cast::<PyAdjustQueryParameters>()
            .map_err(|_| no_converter("RDKit::MolOps::AdjustQueryParameters", p))?
            .borrow()
            .params
            .clone(),
    };
    let core = mol.borrow().mol.clone();
    let adjusted = py.detach(|| adjust_query_properties(&core, &params)).py()?;
    Ok(PyMol { mol: adjusted })
}

/// Returns a new molecule where the query properties of atoms have been modified and generic group queries have been prepared.
#[pyfunction]
#[pyo3(name = "AdjustQueryPropertiesWithGenericGroups", signature = (mol, params = None))]
fn py_adjust_query_properties_with_generic_groups(
    py: Python<'_>,
    mol: &Bound<'_, PyMol>,
    params: Option<&Bound<'_, PyAdjustQueryParameters>>,
) -> PyResult<PyMol> {
    let params = params.map_or_else(AdjustQueryParameters::default, |p| {
        p.borrow().params.clone()
    });
    let core = mol.borrow().mol.clone();
    let adjusted = py
        .detach(|| adjust_query_properties_with_generic_groups(&core, &params))
        .py()?;
    Ok(PyMol { mol: adjusted })
}

/// Sets the generic query tags of the atoms from their atom labels and/or SUP substance groups.
#[pyfunction]
#[pyo3(name = "SetGenericQueriesFromProperties", signature = (mol, useAtomLabels = true, useSGroups = true))]
#[allow(non_snake_case)]
fn py_set_generic_queries_from_properties(
    mol: &Bound<'_, PyMol>,
    useAtomLabels: bool,
    useSGroups: bool,
) {
    set_generic_queries_from_properties(&mut mol.borrow_mut().mol, useAtomLabels, useSGroups);
}

/// Converts the generic query tags of the atoms into SUP substance groups.
#[pyfunction]
#[pyo3(name = "ConvertGenericQueriesToSubstanceGroups")]
fn py_convert_generic_queries_to_substance_groups(mol: &Bound<'_, PyMol>) -> PyResult<()> {
    convert_generic_queries_to_substance_groups(&mut mol.borrow_mut().mol).py()
}

pub fn add(m: &Bound<'_, PyModule>) -> PyResult<()> {
    add_enum_exported(m, "AdjustQueryWhichFlags", &WHICH_FLAGS, None)?;
    m.add_class::<PyAdjustQueryParameters>()?;
    m.add_function(wrap_pyfunction!(py_adjust_query_properties, m)?)?;
    m.add_function(wrap_pyfunction!(
        py_adjust_query_properties_with_generic_groups,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(py_set_generic_queries_from_properties, m)?)?;
    m.add_function(wrap_pyfunction!(
        py_convert_generic_queries_to_substance_groups,
        m
    )?)?;
    Ok(())
}
