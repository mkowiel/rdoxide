//! `rdoxide.Chem.rdFreeSASA` (upstream `External/FreeSASA/Wrap/rdFreeSASA.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAnyMethods;
use rdoxide::external::freesasa::Algorithm;
use rdoxide::external::freesasa::classifier::Classifier;
use rdoxide::external::freesasa::rdfreesasa::{
    SasaOpts, calc_sasa, classify_atoms, make_apolar_atom_query, make_polar_atom_query,
};
use rdoxide::graphmol::atom::Atom;

use crate::atom::PyAtom;
use crate::enums::{add_enum_exported, to_py};
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

fn algorithm_from(v: i64) -> PyResult<Algorithm> {
    match v {
        0 => Ok(Algorithm::LeeRichards),
        1 => Ok(Algorithm::ShrakeRupley),
        _ => Err(PyValueError::new_err("invalid SASAAlgorithm value")),
    }
}

fn classifier_from(v: i64) -> PyResult<Classifier> {
    match v {
        0 => Ok(Classifier::Protor),
        1 => Ok(Classifier::Naccess),
        2 => Ok(Classifier::Oons),
        _ => Err(PyValueError::new_err("invalid SASAClassifier value")),
    }
}

/// FreeSASA options: algorithm, classifier and probe radius
#[pyclass(name = "SASAOpts", module = "rdoxide.Chem.rdFreeSASA")]
#[derive(Debug, Clone)]
struct PySasaOpts {
    opts: SasaOpts,
}

#[pymethods]
impl PySasaOpts {
    /// Constructor takes no arguments, or the algorithm, the classifier and optionally the probe radius
    #[new]
    #[pyo3(signature = (alg = None, cls = None, pr = None))]
    fn new(
        alg: Option<&Bound<'_, PyAny>>,
        cls: Option<&Bound<'_, PyAny>>,
        pr: Option<f64>,
    ) -> PyResult<Self> {
        let mut opts = SasaOpts::default();
        if let Some(a) = alg {
            opts.algorithm = algorithm_from(a.extract()?)?;
        }
        if let Some(c) = cls {
            opts.classifier = classifier_from(c.extract()?)?;
        }
        if let Some(p) = pr {
            opts.probe_radius = p;
        }
        Ok(PySasaOpts { opts })
    }

    #[getter]
    fn algorithm(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(py, "SASAAlgorithm", self.opts.algorithm as i64)
    }

    #[setter(algorithm)]
    fn set_algorithm(&mut self, v: &Bound<'_, PyAny>) -> PyResult<()> {
        self.opts.algorithm = algorithm_from(v.extract()?)?;
        Ok(())
    }

    #[getter]
    fn classifier(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_py(py, "SASAClassifier", self.opts.classifier as i64)
    }

    #[setter(classifier)]
    fn set_classifier(&mut self, v: &Bound<'_, PyAny>) -> PyResult<()> {
        self.opts.classifier = classifier_from(v.extract()?)?;
        Ok(())
    }

    #[getter(probeRadius)]
    fn probe_radius(&self) -> f64 {
        self.opts.probe_radius
    }

    #[setter(probeRadius)]
    fn set_probe_radius(&mut self, v: f64) {
        self.opts.probe_radius = v;
    }
}

/// Classify the atoms in the molecule returning their radii if possible.
/// ARGUMENTS:
///    - mol: molecule to classify
///    - options: FreeSASA options class specifying the classification method.
///               Current classifiers are Protor, NACCESS and OONS
///               classification is stored as atom property 'SASAClass' for the integer value
///                and 'SASAClassName' for the string name of the class, Polar, APolar...
///
/// RETURNS:
///   list of radii where radii[atom.GetIdx()] is the radii of the atom.
///   If classification fails, NONE is returned
#[pyfunction]
#[pyo3(name = "classifyAtoms", signature = (mol, options = None))]
fn classify_atoms_py(
    mut mol: PyRefMut<'_, PyMol>,
    options: Option<PyRef<'_, PySasaOpts>>,
) -> PyResult<Vec<f64>> {
    let opts = options.map_or_else(SasaOpts::default, |o| o.opts);
    let (ok, radii) = classify_atoms(&mut mol.mol, &opts).py()?;
    Ok(if ok { radii } else { Vec::new() })
}

/// Compute the Solvent Accessible Surface Area using the FreeSASA library
/// ARGUMENTS:
///   - mol: The molecule to compute.
///   - radii:  A list of atom raddii where radii[atom.GetIdx()] is the radius of the atom
///             These can be passed in or calculated with classifyAtoms for some proteins
///   - confIdx: Specify the conformer to use for the 3D geometry  [default -1]
///   - query: Pass along a query atom to compute the SASA for a subset of atoms.
///            precanned query atoms can be made with MakeFreeSasaPolarAtomQuery and
///            MakeFreeSasaAPolarAtomQuery for classified polar and apolar atoms respectively.
///   - opts: a SASAOpts class specifying the algorithm to use
///
/// RETURNS:
/// The computed solvent accessible surface area.
#[pyfunction]
#[pyo3(name = "CalcSASA", signature = (mol, radii, confIdx = -1, query = None, opts = None))]
#[allow(non_snake_case)]
fn calc_sasa_py(
    py: Python<'_>,
    mut mol: PyRefMut<'_, PyMol>,
    radii: Vec<f64>,
    confIdx: i32,
    query: Option<PyRef<'_, PyAtom>>,
    opts: Option<PyRef<'_, PySasaOpts>>,
) -> PyResult<f64> {
    let opts = opts.map_or_else(SasaOpts::default, |o| o.opts);
    let q = match &query {
        Some(a) => Some(
            a.with(py, |atom| atom.query().cloned())?
                .ok_or_else(|| PyValueError::new_err("Query is not a query atom!"))?,
        ),
        None => None,
    };
    calc_sasa(&mut mol.mol, &radii, confIdx, q.as_ref(), &opts).py()
}

fn query_atom(py: Python<'_>, q: rdoxide::graphmol::queries::AtomQuery) -> PyResult<Py<PyAtom>> {
    let mut atom = Atom::new(0);
    atom.set_query(Some(q));
    PyAtom::owned(py, atom)
}

/// Returns an APolar atom query for use with CalcSASA.  An apolar atom has the SASAClass
/// and SASAClassName set to the APOLAR class.  (see classifyAtoms)
#[pyfunction]
#[pyo3(name = "MakeFreeSasaAPolarAtomQuery")]
fn make_apolar(py: Python<'_>) -> PyResult<Py<PyAtom>> {
    query_atom(py, make_apolar_atom_query())
}

/// Returns a polar atom query for use with CalcSASA.  An polar atom has the SASAClass
/// and SASAClassName set to the POLAR class.  (see classifyAtoms)
#[pyfunction]
#[pyo3(name = "MakeFreeSasaPolarAtomQuery")]
fn make_polar(py: Python<'_>) -> PyResult<Py<PyAtom>> {
    query_atom(py, make_polar_atom_query())
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing rdFreeSASA classes and functions.",
    )?;
    add_enum_exported(
        m,
        "SASAAlgorithm",
        &[("LeeRichards", 0), ("ShrakeRupley", 1)],
        None,
    )?;
    add_enum_exported(
        m,
        "SASAClassifier",
        &[("Protor", 0), ("NACCESS", 1), ("OONS", 2)],
        None,
    )?;
    add_enum_exported(
        m,
        "SASAClass",
        &[("Unclassified", 0), ("APolar", 1), ("Polar", 2)],
        None,
    )?;
    m.add_class::<PySasaOpts>()?;
    m.add_function(wrap_pyfunction!(classify_atoms_py, m)?)?;
    m.add_function(wrap_pyfunction!(calc_sasa_py, m)?)?;
    m.add_function(wrap_pyfunction!(make_apolar, m)?)?;
    m.add_function(wrap_pyfunction!(make_polar, m)?)?;
    Ok(())
}
