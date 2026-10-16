//! `rdoxide.Chem.Scaffolds.rdScaffoldNetwork` (upstream
//! `GraphMol/ScaffoldNetwork/Wrap/rdScaffoldNetwork.cpp`). Pickling of
//! networks is not supported (upstream needs Boost serialization).

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use rdoxide::graphmol::romol::ROMol;
use rdoxide::graphmol::scaffold_network::{
    self, EdgeType, NetworkEdge, ScaffoldNetwork, ScaffoldNetworkParams,
};

use crate::enums::add_enum;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

/// Scaffold network parameters
#[pyclass(
    name = "ScaffoldNetworkParams",
    module = "rdoxide.Chem.Scaffolds.rdScaffoldNetwork",
    skip_from_py_object
)]
#[derive(Clone)]
struct PyScaffoldNetworkParams {
    params: ScaffoldNetworkParams,
}

macro_rules! bool_attrs {
    ($($(#[$doc:meta])* $py:literal, $get:ident, $set:ident, $field:ident;)*) => {
        #[pymethods]
        impl PyScaffoldNetworkParams {
            $(
                $(#[$doc])*
                #[getter($py)]
                fn $get(&self) -> bool {
                    self.params.$field
                }
                #[setter($py)]
                fn $set(&mut self, v: bool) {
                    self.params.$field = v;
                }
            )*
        }
    };
}

bool_attrs! {
    /// include scaffolds with all atoms replaced by dummies
    "includeGenericScaffolds", get_igs, set_igs, include_generic_scaffolds;
    /// include scaffolds with all bonds replaced by single bonds
    "includeGenericBondScaffolds", get_igbs, set_igbs, include_generic_bond_scaffolds;
    /// remove attachment points from scaffolds and include the result
    "includeScaffoldsWithoutAttachments", get_iswo, set_iswo, include_scaffolds_without_attachments;
    /// Include the version of the scaffold with attachment points
    "includeScaffoldsWithAttachments", get_isw, set_isw, include_scaffolds_with_attachments;
    /// Include molecules names of the input molecules
    "includeNames", get_in, set_in, include_names;
    /// keep only the first fragment from the bond breaking rule
    "keepOnlyFirstFragment", get_kof, set_kof, keep_only_first_fragment;
    /// Do a pruning/flattening step before starting fragmenting
    "pruneBeforeFragmenting", get_pbf, set_pbf, prune_before_fragmenting;
    /// remove isotopes when flattening
    "flattenIsotopes", get_fi, set_fi, flatten_isotopes;
    /// remove chirality and bond stereo when flattening
    "flattenChirality", get_fc, set_fc, flatten_chirality;
    /// keep only the largest fragment when doing flattening
    "flattenKeepLargest", get_fkl, set_fkl, flatten_keep_largest;
    /// keep track of the number of molecules each scaffold was found in
    "collectMolCounts", get_cmc, set_cmc, collect_mol_counts;
}

#[pymethods]
impl PyScaffoldNetworkParams {
    /// Constructor taking a list of Reaction SMARTS for the fragmentation reactions
    #[new]
    #[pyo3(signature = (bondBreakerSmartsList = None))]
    #[allow(non_snake_case)]
    fn new(bondBreakerSmartsList: Option<Vec<String>>) -> PyResult<Self> {
        let params = match bondBreakerSmartsList {
            Some(list) => ScaffoldNetworkParams::new(&list).py()?,
            None => ScaffoldNetworkParams::default(),
        };
        Ok(PyScaffoldNetworkParams { params })
    }

    fn __setattr__(slf: &Bound<'_, Self>, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        crate::pyconv::set_declared_attribute(slf.as_any(), name, value)
    }
}

/// A scaffold network edge
#[pyclass(
    name = "NetworkEdge",
    module = "rdoxide.Chem.Scaffolds.rdScaffoldNetwork",
    skip_from_py_object
)]
#[derive(Clone)]
struct PyNetworkEdge {
    edge: NetworkEdge,
}

#[pymethods]
impl PyNetworkEdge {
    /// index of the begin node in node list
    #[getter(beginIdx)]
    fn begin_idx(&self) -> usize {
        self.edge.begin_idx
    }

    /// index of the end node in node list
    #[getter(endIdx)]
    fn end_idx(&self) -> usize {
        self.edge.end_idx
    }

    /// type of the edge
    #[getter(r#type)]
    fn edge_type(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let cls = py
            .import("rdoxide.Chem.Scaffolds.rdScaffoldNetwork")?
            .getattr("EdgeType")?;
        Ok(cls.call1((self.edge.edge_type as i64,))?.unbind())
    }

    fn __str__(&self) -> String {
        self.edge.to_string()
    }
}

/// A scaffold network
#[pyclass(
    name = "ScaffoldNetwork",
    module = "rdoxide.Chem.Scaffolds.rdScaffoldNetwork"
)]
#[derive(Default)]
struct PyScaffoldNetwork {
    network: ScaffoldNetwork,
}

#[pymethods]
impl PyScaffoldNetwork {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// the sequence of SMILES defining the nodes
    #[getter]
    fn nodes(&self) -> Vec<String> {
        self.network.nodes.clone()
    }

    /// the number of times each node was encountered while building the network.
    #[getter]
    fn counts(&self) -> Vec<u32> {
        self.network.counts.clone()
    }

    /// the number of moleclues each node was found in.
    #[getter(molCounts)]
    fn mol_counts(&self) -> Vec<u32> {
        self.network.mol_counts.clone()
    }

    /// the sequence of network edges
    #[getter]
    fn edges(&self) -> Vec<PyNetworkEdge> {
        self.network
            .edges
            .iter()
            .map(|&edge| PyNetworkEdge { edge })
            .collect()
    }

    fn __reduce__(&self) -> PyResult<()> {
        Err(PyRuntimeError::new_err(
            "Pickling of ScaffoldNetwork instances is not enabled",
        ))
    }
}

/// The molecules of a Python sequence; `None` entries are errors.
fn molecules(mols: &Bound<'_, PyAny>) -> PyResult<Vec<ROMol>> {
    let mut res = Vec::new();
    for item in mols.try_iter()? {
        let item = item?;
        if item.is_none() {
            return Err(PyValueError::new_err(
                "updateScaffoldNetwork called with null molecule",
            ));
        }
        res.push(item.extract::<PyRef<'_, PyMol>>()?.mol.clone());
    }
    Ok(res)
}

/// create (and return) a new network from a sequence of molecules
#[pyfunction]
#[pyo3(name = "CreateScaffoldNetwork", signature = (mols, params))]
fn create_scaffold_network(
    py: Python<'_>,
    mols: &Bound<'_, PyAny>,
    params: PyRef<'_, PyScaffoldNetworkParams>,
) -> PyResult<PyScaffoldNetwork> {
    let mols = molecules(mols)?;
    let params = &params.params;
    let network = py
        .detach(|| scaffold_network::create_scaffold_network(&mols, params))
        .py()?;
    Ok(PyScaffoldNetwork { network })
}

/// update an existing network by adding molecules
#[pyfunction]
#[pyo3(name = "UpdateScaffoldNetwork", signature = (mols, network, params))]
fn update_scaffold_network(
    py: Python<'_>,
    mols: &Bound<'_, PyAny>,
    mut network: PyRefMut<'_, PyScaffoldNetwork>,
    params: PyRef<'_, PyScaffoldNetworkParams>,
) -> PyResult<()> {
    let mols = molecules(mols)?;
    let params = &params.params;
    let network = &mut network.network;
    py.detach(|| scaffold_network::update_scaffold_network(&mols, network, params))
        .py()
}

/// Returns parameters for generating scaffolds using BRICS fragmentation rules
#[pyfunction]
#[pyo3(name = "BRICSScaffoldParams")]
fn brics_scaffold_params() -> PyResult<PyScaffoldNetworkParams> {
    Ok(PyScaffoldNetworkParams {
        params: scaffold_network::get_brics_network_params().py()?,
    })
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions for creating a Scaffold Network",
    )?;
    add_enum(
        m,
        "EdgeType",
        &[
            ("Fragment", EdgeType::Fragment as i64),
            ("Generic", EdgeType::Generic as i64),
            ("GenericBond", EdgeType::GenericBond as i64),
            ("RemoveAttachment", EdgeType::RemoveAttachment as i64),
            ("Initialize", EdgeType::Initialize as i64),
        ],
        None,
    )?;
    m.add_class::<PyScaffoldNetworkParams>()?;
    m.add_class::<PyNetworkEdge>()?;
    m.add_class::<PyScaffoldNetwork>()?;
    m.add_function(wrap_pyfunction!(create_scaffold_network, m)?)?;
    m.add_function(wrap_pyfunction!(update_scaffold_network, m)?)?;
    m.add_function(wrap_pyfunction!(brics_scaffold_params, m)?)?;
    Ok(())
}
