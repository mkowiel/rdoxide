//! `rdoxide.Chem.rdFingerprintGenerator` (upstream
//! `GraphMol/Fingerprints/Wrap/{FingerprintGeneratorWrapper,MorganWrapper}.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rdoxide::graphmol::fingerprints::atom_pair::AtomPairArguments;
use rdoxide::graphmol::fingerprints::json::{generator_from_json, generator_to_json};
use rdoxide::graphmol::fingerprints::morgan::MorganArguments;
use rdoxide::graphmol::fingerprints::rdkit_fp::RDKitFPArguments;
use rdoxide::graphmol::fingerprints::torsion::TopologicalTorsionArguments;
use rdoxide::graphmol::fingerprints::{
    AdditionalOutput, AtomInvariantsGenerator, BondInvariantsGenerator, EnvGenerator,
    FingerprintArguments, FingerprintFuncArguments, FingerprintGenerator,
};
use rdoxide::graphmol::romol::ROMol;
use rdoxide::rdgeneral::threads::par_map;

use crate::datastructs::{PyExplicitBitVect, PySparseBitVect};
use crate::datastructs_siv::{PyUIntSIV, PyULongSIV};
use crate::enums::add_enum_exported;
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

#[pyclass(
    name = "AtomInvariantsGenerator",
    module = "rdoxide.Chem.rdFingerprintGenerator"
)]
#[derive(Debug)]
pub struct PyAtomInvGen {
    gen_: AtomInvariantsGenerator,
}

#[pyclass(
    name = "BondInvariantsGenerator",
    module = "rdoxide.Chem.rdFingerprintGenerator"
)]
#[derive(Debug)]
pub struct PyBondInvGen {
    gen_: BondInvariantsGenerator,
}

/// `AdditionalOutput`
#[pyclass(
    name = "AdditionalOutput",
    module = "rdoxide.Chem.rdFingerprintGenerator"
)]
#[derive(Debug, Default)]
pub struct PyAdditionalOutput {
    ao: AdditionalOutput,
}

#[pymethods]
impl PyAdditionalOutput {
    #[new]
    fn new() -> Self {
        Self::default()
    }
    #[pyo3(name = "AllocateAtomToBits")]
    fn allocate_atom_to_bits(&mut self) {
        self.ao.atom_to_bits = Some(Vec::new());
    }
    #[pyo3(name = "AllocateBitInfoMap")]
    fn allocate_bit_info_map(&mut self) {
        self.ao.bit_info_map = Some(Default::default());
    }
    #[pyo3(name = "AllocateBitPaths")]
    fn allocate_bit_paths(&mut self) {
        self.ao.bit_paths = Some(Default::default());
    }
    #[pyo3(name = "AllocateAtomCounts")]
    fn allocate_atom_counts(&mut self) {
        self.ao.atom_counts = Some(Vec::new());
    }
    #[pyo3(name = "AllocateAtomsPerBit")]
    fn allocate_atoms_per_bit(&mut self) {
        self.ao.atoms_per_bit = Some(Default::default());
    }
    #[pyo3(name = "CollectAtomToBits")]
    fn collect_atom_to_bits(&mut self) {
        self.allocate_atom_to_bits();
    }
    #[pyo3(name = "CollectBitInfoMap")]
    fn collect_bit_info_map(&mut self) {
        self.allocate_bit_info_map();
    }
    #[pyo3(name = "CollectBitPaths")]
    fn collect_bit_paths(&mut self) {
        self.allocate_bit_paths();
    }
    #[pyo3(name = "CollectAtomCounts")]
    fn collect_atom_counts(&mut self) {
        self.allocate_atom_counts();
    }
    #[pyo3(name = "CollectAtomsPerBit")]
    fn collect_atoms_per_bit(&mut self) {
        self.allocate_atoms_per_bit();
    }
    #[pyo3(name = "GetAtomToBits")]
    fn get_atom_to_bits(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let Some(v) = &self.ao.atom_to_bits else {
            return Ok(py.None());
        };
        let items = v
            .iter()
            .map(|l| PyTuple::new(py, l))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyTuple::new(py, items)?.into_any().unbind())
    }
    #[pyo3(name = "GetBitInfoMap")]
    fn get_bit_info_map(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let Some(m) = &self.ao.bit_info_map else {
            return Ok(py.None());
        };
        let d = PyDict::new(py);
        for (k, v) in m {
            d.set_item(k, PyTuple::new(py, v)?)?;
        }
        Ok(d.into_any().unbind())
    }
    #[pyo3(name = "GetBitPaths")]
    fn get_bit_paths(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let Some(m) = &self.ao.bit_paths else {
            return Ok(py.None());
        };
        nested_map(py, m)
    }
    #[pyo3(name = "GetAtomCounts")]
    fn get_atom_counts(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let Some(v) = &self.ao.atom_counts else {
            return Ok(py.None());
        };
        Ok(PyTuple::new(py, v)?.into_any().unbind())
    }
    #[pyo3(name = "GetAtomsPerBit")]
    fn get_atoms_per_bit(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let Some(m) = &self.ao.atoms_per_bit else {
            return Ok(py.None());
        };
        nested_map(py, m)
    }
}

fn nested_map<T: Copy + for<'py> IntoPyObject<'py>>(
    py: Python<'_>,
    m: &std::collections::BTreeMap<u64, Vec<Vec<T>>>,
) -> PyResult<Py<PyAny>> {
    let d = PyDict::new(py);
    for (k, v) in m {
        let inner = v
            .iter()
            .map(|l| PyTuple::new(py, l.iter().copied()))
            .collect::<PyResult<Vec<_>>>()?;
        d.set_item(k, PyTuple::new(py, inner)?)?;
    }
    Ok(d.into_any().unbind())
}

/// `FingerprintGenerator64`
#[pyclass(
    name = "FingerprintGenerator64",
    module = "rdoxide.Chem.rdFingerprintGenerator"
)]
#[derive(Debug)]
pub struct PyFPGen {
    pub(crate) gen_: FingerprintGenerator,
}

/// Owned versions of the optional per-call arguments (`convertPyArguments`).
struct CallArgs {
    from_atoms: Option<Vec<u32>>,
    ignore_atoms: Option<Vec<u32>>,
    atom_invs: Option<Vec<u32>>,
    bond_invs: Option<Vec<u32>>,
    conf_id: i32,
}

fn opt_list(o: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Vec<u32>>> {
    match o {
        Some(o) if !o.is_none() => {
            let v: Vec<u32> = o.extract()?;
            Ok(if v.is_empty() { None } else { Some(v) })
        }
        _ => Ok(None),
    }
}

impl CallArgs {
    fn new(
        from_atoms: Option<&Bound<'_, PyAny>>,
        ignore_atoms: Option<&Bound<'_, PyAny>>,
        conf_id: i32,
        atom_invs: Option<&Bound<'_, PyAny>>,
        bond_invs: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Ok(CallArgs {
            from_atoms: opt_list(from_atoms)?,
            ignore_atoms: opt_list(ignore_atoms)?,
            atom_invs: opt_list(atom_invs)?,
            bond_invs: opt_list(bond_invs)?,
            conf_id,
        })
    }
    fn func_args<'a>(
        &'a self,
        ao: Option<&'a mut AdditionalOutput>,
    ) -> FingerprintFuncArguments<'a> {
        FingerprintFuncArguments {
            from_atoms: self.from_atoms.as_deref(),
            ignore_atoms: self.ignore_atoms.as_deref(),
            conf_id: self.conf_id,
            additional_output: ao,
            custom_atom_invariants: self.atom_invs.as_deref(),
            custom_bond_invariants: self.bond_invs.as_deref(),
        }
    }
}

/// Runs `f` with the call arguments and the optional AdditionalOutput.
fn with_args<R>(
    additional_output: Option<&Bound<'_, PyAny>>,
    args: CallArgs,
    f: impl FnOnce(&mut FingerprintFuncArguments<'_>) -> rdoxide::Result<R>,
) -> PyResult<R> {
    match additional_output {
        Some(o) if !o.is_none() => {
            let mut ao = o.extract::<PyRefMut<'_, PyAdditionalOutput>>()?;
            let mut fa = args.func_args(Some(&mut ao.ao));
            f(&mut fa).py()
        }
        _ => {
            let mut fa = args.func_args(None);
            f(&mut fa).py()
        }
    }
}

fn numpy_array(py: Python<'_>, vals: Vec<u32>, dtype: &str) -> PyResult<Py<PyAny>> {
    let np = py.import("numpy")?;
    Ok(np.call_method1("array", (vals, dtype))?.unbind())
}

#[pymethods]
impl PyFPGen {
    #[pyo3(name = "GetFingerprint", signature = (mol, fromAtoms = None, ignoreAtoms = None, confId = -1, customAtomInvariants = None, customBondInvariants = None, additionalOutput = None))]
    #[allow(non_snake_case, clippy::too_many_arguments)]
    fn get_fingerprint(
        &self,
        mol: PyRef<'_, PyMol>,
        fromAtoms: Option<&Bound<'_, PyAny>>,
        ignoreAtoms: Option<&Bound<'_, PyAny>>,
        confId: i32,
        customAtomInvariants: Option<&Bound<'_, PyAny>>,
        customBondInvariants: Option<&Bound<'_, PyAny>>,
        additionalOutput: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<PyExplicitBitVect> {
        let args = CallArgs::new(
            fromAtoms,
            ignoreAtoms,
            confId,
            customAtomInvariants,
            customBondInvariants,
        )?;
        let r = with_args(additionalOutput, args, |fa| {
            self.gen_.fingerprint(&mol.mol, fa)
        })?;
        Ok(PyExplicitBitVect::wrap(r))
    }
    #[pyo3(name = "GetSparseFingerprint", signature = (mol, fromAtoms = None, ignoreAtoms = None, confId = -1, customAtomInvariants = None, customBondInvariants = None, additionalOutput = None))]
    #[allow(non_snake_case, clippy::too_many_arguments)]
    fn get_sparse_fingerprint(
        &self,
        mol: PyRef<'_, PyMol>,
        fromAtoms: Option<&Bound<'_, PyAny>>,
        ignoreAtoms: Option<&Bound<'_, PyAny>>,
        confId: i32,
        customAtomInvariants: Option<&Bound<'_, PyAny>>,
        customBondInvariants: Option<&Bound<'_, PyAny>>,
        additionalOutput: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<PySparseBitVect> {
        let args = CallArgs::new(
            fromAtoms,
            ignoreAtoms,
            confId,
            customAtomInvariants,
            customBondInvariants,
        )?;
        let r = with_args(additionalOutput, args, |fa| {
            self.gen_.sparse_fingerprint(&mol.mol, fa)
        })?;
        Ok(PySparseBitVect::wrap(r))
    }
    #[pyo3(name = "GetCountFingerprint", signature = (mol, fromAtoms = None, ignoreAtoms = None, confId = -1, customAtomInvariants = None, customBondInvariants = None, additionalOutput = None))]
    #[allow(non_snake_case, clippy::too_many_arguments)]
    fn get_count_fingerprint(
        &self,
        mol: PyRef<'_, PyMol>,
        fromAtoms: Option<&Bound<'_, PyAny>>,
        ignoreAtoms: Option<&Bound<'_, PyAny>>,
        confId: i32,
        customAtomInvariants: Option<&Bound<'_, PyAny>>,
        customBondInvariants: Option<&Bound<'_, PyAny>>,
        additionalOutput: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<PyUIntSIV> {
        let args = CallArgs::new(
            fromAtoms,
            ignoreAtoms,
            confId,
            customAtomInvariants,
            customBondInvariants,
        )?;
        let r = with_args(additionalOutput, args, |fa| {
            self.gen_.count_fingerprint(&mol.mol, fa)
        })?;
        Ok(PyUIntSIV::wrap(r))
    }
    #[pyo3(name = "GetSparseCountFingerprint", signature = (mol, fromAtoms = None, ignoreAtoms = None, confId = -1, customAtomInvariants = None, customBondInvariants = None, additionalOutput = None))]
    #[allow(non_snake_case, clippy::too_many_arguments)]
    fn get_sparse_count_fingerprint(
        &self,
        mol: PyRef<'_, PyMol>,
        fromAtoms: Option<&Bound<'_, PyAny>>,
        ignoreAtoms: Option<&Bound<'_, PyAny>>,
        confId: i32,
        customAtomInvariants: Option<&Bound<'_, PyAny>>,
        customBondInvariants: Option<&Bound<'_, PyAny>>,
        additionalOutput: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<PyULongSIV> {
        let args = CallArgs::new(
            fromAtoms,
            ignoreAtoms,
            confId,
            customAtomInvariants,
            customBondInvariants,
        )?;
        let r = with_args(additionalOutput, args, |fa| {
            self.gen_.sparse_count_fingerprint(&mol.mol, fa)
        })?;
        Ok(PyULongSIV::wrap(r))
    }

    #[pyo3(name = "GetFingerprintAsNumPy", signature = (mol, fromAtoms = None, ignoreAtoms = None, confId = -1, customAtomInvariants = None, customBondInvariants = None, additionalOutput = None))]
    #[allow(non_snake_case, clippy::too_many_arguments)]
    fn get_fingerprint_as_numpy(
        &self,
        py: Python<'_>,
        mol: PyRef<'_, PyMol>,
        fromAtoms: Option<&Bound<'_, PyAny>>,
        ignoreAtoms: Option<&Bound<'_, PyAny>>,
        confId: i32,
        customAtomInvariants: Option<&Bound<'_, PyAny>>,
        customBondInvariants: Option<&Bound<'_, PyAny>>,
        additionalOutput: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        use rdoxide::datastructs::BitVect;
        let args = CallArgs::new(
            fromAtoms,
            ignoreAtoms,
            confId,
            customAtomInvariants,
            customBondInvariants,
        )?;
        let bv = with_args(additionalOutput, args, |fa| {
            self.gen_.fingerprint(&mol.mol, fa)
        })?;
        let mut vals = vec![0u32; bv.num_bits() as usize];
        for b in bv.on_bits() {
            vals[b as usize] = 1;
        }
        numpy_array(py, vals, "uint8")
    }

    #[pyo3(name = "GetCountFingerprintAsNumPy", signature = (mol, fromAtoms = None, ignoreAtoms = None, confId = -1, customAtomInvariants = None, customBondInvariants = None, additionalOutput = None))]
    #[allow(non_snake_case, clippy::too_many_arguments)]
    fn get_count_fingerprint_as_numpy(
        &self,
        py: Python<'_>,
        mol: PyRef<'_, PyMol>,
        fromAtoms: Option<&Bound<'_, PyAny>>,
        ignoreAtoms: Option<&Bound<'_, PyAny>>,
        confId: i32,
        customAtomInvariants: Option<&Bound<'_, PyAny>>,
        customBondInvariants: Option<&Bound<'_, PyAny>>,
        additionalOutput: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let args = CallArgs::new(
            fromAtoms,
            ignoreAtoms,
            confId,
            customAtomInvariants,
            customBondInvariants,
        )?;
        let fp = with_args(additionalOutput, args, |fa| {
            self.gen_.count_fingerprint(&mol.mol, fa)
        })?;
        let mut vals = vec![0u32; fp.length() as usize];
        for (&k, &v) in fp.nonzero_elements() {
            vals[k as usize] = v as u32;
        }
        numpy_array(py, vals, "uint32")
    }

    #[pyo3(name = "GetFingerprints", signature = (mols, numThreads = 1))]
    #[allow(non_snake_case)]
    fn get_fingerprints<'py>(
        &self,
        py: Python<'py>,
        mols: &Bound<'py, PyAny>,
        numThreads: i32,
    ) -> PyResult<Bound<'py, PyTuple>> {
        self.bulk(py, mols, numThreads, |g, m| {
            Ok(PyExplicitBitVect::wrap(
                g.fingerprint(m, &mut Default::default())?,
            ))
        })
    }
    #[pyo3(name = "GetCountFingerprints", signature = (mols, numThreads = 1))]
    #[allow(non_snake_case)]
    fn get_count_fingerprints<'py>(
        &self,
        py: Python<'py>,
        mols: &Bound<'py, PyAny>,
        numThreads: i32,
    ) -> PyResult<Bound<'py, PyTuple>> {
        self.bulk(py, mols, numThreads, |g, m| {
            Ok(PyUIntSIV::wrap(
                g.count_fingerprint(m, &mut Default::default())?,
            ))
        })
    }
    #[pyo3(name = "GetSparseFingerprints", signature = (mols, numThreads = 1))]
    #[allow(non_snake_case)]
    fn get_sparse_fingerprints<'py>(
        &self,
        py: Python<'py>,
        mols: &Bound<'py, PyAny>,
        numThreads: i32,
    ) -> PyResult<Bound<'py, PyTuple>> {
        self.bulk(py, mols, numThreads, |g, m| {
            Ok(PySparseBitVect::wrap(
                g.sparse_fingerprint(m, &mut Default::default())?,
            ))
        })
    }
    #[pyo3(name = "GetSparseCountFingerprints", signature = (mols, numThreads = 1))]
    #[allow(non_snake_case)]
    fn get_sparse_count_fingerprints<'py>(
        &self,
        py: Python<'py>,
        mols: &Bound<'py, PyAny>,
        numThreads: i32,
    ) -> PyResult<Bound<'py, PyTuple>> {
        self.bulk(py, mols, numThreads, |g, m| {
            Ok(PyULongSIV::wrap(
                g.sparse_count_fingerprint(m, &mut Default::default())?,
            ))
        })
    }

    /// Serialize the generator to a JSON string.
    #[pyo3(name = "ToJSON")]
    fn to_json(&self) -> PyResult<String> {
        generator_to_json(&self.gen_).py()
    }

    #[pyo3(name = "GetInfoString")]
    fn get_info_string(&self) -> String {
        self.gen_.info_string()
    }

    /// Returns the fingerprint options object (changes to it modify the generator)
    #[pyo3(name = "GetOptions")]
    fn get_options(slf: Bound<'_, Self>) -> PyResult<Py<PyAny>> {
        let py = slf.py();
        let base = PyFingerprintOptions {
            owner: slf.unbind(),
        };
        let kind = match base.owner.borrow(py).gen_.env {
            EnvGenerator::Morgan(_) => 0,
            EnvGenerator::AtomPair(_) => 1,
            EnvGenerator::TopologicalTorsion(_) => 2,
            EnvGenerator::RDKitFP(_) => 3,
        };
        let init = PyClassInitializer::from(base);
        match kind {
            0 => Ok(Py::new(py, init.add_subclass(PyMorganOptions))?.into_any()),
            1 => Ok(Py::new(py, init.add_subclass(PyAtomPairOptions))?.into_any()),
            2 => Ok(Py::new(py, init.add_subclass(PyTorsionOptions))?.into_any()),
            _ => Ok(Py::new(py, init.add_subclass(PyRDKitFPOptions))?.into_any()),
        }
    }
}

impl PyFPGen {
    fn bulk<'py, T>(
        &self,
        py: Python<'py>,
        mols: &Bound<'py, PyAny>,
        num_threads: i32,
        f: impl Fn(&FingerprintGenerator, &ROMol) -> rdoxide::error::Result<T> + Sync,
    ) -> PyResult<Bound<'py, PyTuple>>
    where
        T: pyo3::PyClass + Into<PyClassInitializer<T>> + Send,
    {
        bulk_map(py, &self.gen_, mols, num_threads, f)
    }
}

/// Runs `f` over every molecule of `mols` (`None` entries stay `None`) on
/// `num_threads` threads with the GIL released. Results keep the input order
/// and the first error in that order is raised.
fn bulk_map<'py, T>(
    py: Python<'py>,
    gen_: &FingerprintGenerator,
    mols: &Bound<'py, PyAny>,
    num_threads: i32,
    f: impl Fn(&FingerprintGenerator, &ROMol) -> rdoxide::error::Result<T> + Sync,
) -> PyResult<Bound<'py, PyTuple>>
where
    T: pyo3::PyClass + Into<PyClassInitializer<T>> + Send,
{
    let mut borrowed = Vec::new();
    for m in mols.try_iter()? {
        let m = m?;
        borrowed.push(if m.is_none() {
            None
        } else {
            Some(m.extract::<PyRef<'_, PyMol>>()?)
        });
    }
    let refs: Vec<Option<&ROMol>> = borrowed
        .iter()
        .map(|m| m.as_ref().map(|m| &m.mol))
        .collect();
    let results = py.detach(|| par_map(&refs, num_threads, |m| m.map(|m| f(gen_, m)).transpose()));
    let mut out = Vec::with_capacity(results.len());
    for r in results {
        out.push(match r.py()? {
            Some(v) => Py::new(py, v)?.into_any(),
            None => py.None(),
        });
    }
    PyTuple::new(py, out)
}

/// The generator `FPType` stands for, with the upstream defaults of the
/// matching `Get*Generator` function (Morgan uses radius 2).
fn default_generator(py: Python<'_>, kind: i64) -> PyResult<FingerprintGenerator> {
    let module = py.import("rdoxide.Chem.rdFingerprintGenerator")?;
    let (factory, args) = match kind {
        FP_RDKIT => ("GetRDKitFPGenerator", PyTuple::empty(py)),
        FP_MORGAN => ("GetMorganGenerator", PyTuple::new(py, [2u32])?),
        FP_ATOM_PAIR => ("GetAtomPairGenerator", PyTuple::empty(py)),
        FP_TORSION => ("GetTopologicalTorsionGenerator", PyTuple::empty(py)),
        _ => {
            return Err(PyValueError::new_err(
                "Fingerprint type not implemented for bulk generation",
            ));
        }
    };
    let gen_ = module.getattr(factory)?.call1(args)?;
    Ok(gen_.extract::<PyRef<'_, PyFPGen>>()?.gen_.clone())
}

const FP_RDKIT: i64 = 1;
const FP_MORGAN: i64 = 2;
const FP_ATOM_PAIR: i64 = 3;
const FP_TORSION: i64 = 4;

macro_rules! bulk_function {
    ($(#[$doc:meta])* $fn_name:ident, $py_name:literal, $method:ident, $wrap:ty) => {
        $(#[$doc])*
        #[pyfunction]
        #[pyo3(name = $py_name, signature = (molecules = None, fpType = None, numThreads = 0))]
        #[allow(non_snake_case)]
        fn $fn_name<'py>(
            py: Python<'py>,
            molecules: Option<&Bound<'py, PyAny>>,
            fpType: Option<&Bound<'py, PyAny>>,
            numThreads: i32,
        ) -> PyResult<Bound<'py, PyList>> {
            let kind = match fpType {
                Some(t) => t
                    .extract::<i64>()
                    .or_else(|_| t.getattr("real")?.extract::<i64>())?,
                None => FP_MORGAN,
            };
            let gen_ = default_generator(py, kind)?;
            let none = PyList::empty(py);
            let mols = molecules.unwrap_or(none.as_any());
            let res = bulk_map(py, &gen_, mols, numThreads, |g, m| {
                Ok(<$wrap>::wrap(g.$method(m, &mut Default::default())?))
            })?;
            PyList::new(py, res)
        }
    };
}

bulk_function!(
    /// Generates sparse count fingerprints for a list of molecules.
    py_sparse_count_fps, "GetSparseCountFPs", sparse_count_fingerprint, PyULongSIV
);
bulk_function!(
    /// Generates sparse bit fingerprints for a list of molecules.
    py_sparse_fps, "GetSparseFPs", sparse_fingerprint, PySparseBitVect
);
bulk_function!(
    /// Generates count fingerprints for a list of molecules.
    py_count_fps, "GetCountFPs", count_fingerprint, PyUIntSIV
);
bulk_function!(
    /// Generates fingerprints for a list of molecules.
    py_fps, "GetFPs", fingerprint, PyExplicitBitVect
);

/// Deserialize a FingerprintGenerator from a JSON string
#[pyfunction]
#[pyo3(name = "FingerprintGeneratorFromJSON", signature = (jsonString))]
#[allow(non_snake_case)]
fn fingerprint_generator_from_json(jsonString: &str) -> PyResult<PyFPGen> {
    Ok(PyFPGen {
        gen_: generator_from_json(jsonString).py()?,
    })
}

/// `FingerprintOptions`: a view on the options of its generator.
#[pyclass(
    name = "FingerprintOptions",
    module = "rdoxide.Chem.rdFingerprintGenerator",
    subclass
)]
#[derive(Debug)]
pub struct PyFingerprintOptions {
    owner: Py<PyFPGen>,
}

impl PyFingerprintOptions {
    fn with<R>(&self, py: Python<'_>, f: impl FnOnce(&mut FingerprintGenerator) -> R) -> R {
        f(&mut self.owner.borrow_mut(py).gen_)
    }
}

#[pymethods]
impl PyFingerprintOptions {
    /// use count simulation
    #[getter(countSimulation)]
    fn get_count_simulation(&self, py: Python<'_>) -> bool {
        self.with(py, |g| g.args.count_simulation)
    }
    #[setter(countSimulation)]
    fn set_count_simulation(&self, py: Python<'_>, v: bool) {
        self.with(py, |g| g.args.count_simulation = v)
    }
    /// include chirality in atom invariants (not for all fingerprints)
    #[getter(includeChirality)]
    fn get_include_chirality(&self, py: Python<'_>) -> bool {
        self.with(py, |g| g.args.include_chirality)
    }
    #[setter(includeChirality)]
    fn set_include_chirality(&self, py: Python<'_>, v: bool) {
        self.with(py, |g| g.args.include_chirality = v)
    }
    /// size of the fingerprints created
    #[getter(fpSize)]
    fn get_fp_size(&self, py: Python<'_>) -> u32 {
        self.with(py, |g| g.args.fp_size)
    }
    #[setter(fpSize)]
    fn set_fp_size(&self, py: Python<'_>, v: u32) {
        self.with(py, |g| g.args.fp_size = v)
    }
    /// number of bits to set for each feature
    #[getter(numBitsPerFeature)]
    fn get_num_bits_per_feature(&self, py: Python<'_>) -> u32 {
        self.with(py, |g| g.args.num_bits_per_feature)
    }
    #[setter(numBitsPerFeature)]
    fn set_num_bits_per_feature(&self, py: Python<'_>, v: u32) {
        self.with(py, |g| g.args.num_bits_per_feature = v)
    }
    /// set the bins for the count bounds
    #[pyo3(name = "SetCountBounds", signature = (bounds))]
    fn set_count_bounds(&self, py: Python<'_>, bounds: Vec<u32>) {
        self.with(py, |g| g.args.count_bounds = bounds)
    }
}

/// `MorganFingerprintOptions`
#[pyclass(name = "MorganFingerprintOptions", module = "rdoxide.Chem.rdFingerprintGenerator", extends = PyFingerprintOptions)]
#[derive(Debug)]
pub struct PyMorganOptions;

fn morgan_args(g: &mut FingerprintGenerator) -> &mut MorganArguments {
    match &mut g.env {
        EnvGenerator::Morgan(m) => m,
        _ => unreachable!("MorganFingerprintOptions of a non-Morgan generator"),
    }
}

fn atom_pair_args(g: &mut FingerprintGenerator) -> &mut AtomPairArguments {
    match &mut g.env {
        EnvGenerator::AtomPair(a) => a,
        _ => unreachable!("AtomPairFingerprintOptions of a non-atom-pair generator"),
    }
}

/// `AtomPairFingerprintOptions`
#[pyclass(name = "AtomPairFingerprintOptions", module = "rdoxide.Chem.rdFingerprintGenerator", extends = PyFingerprintOptions)]
#[derive(Debug)]
pub struct PyAtomPairOptions;

#[pymethods]
impl PyAtomPairOptions {
    /// use 2D distances
    #[getter(use2D)]
    fn get_use_2d(slf: PyRef<'_, Self>) -> bool {
        slf.as_super().with(slf.py(), |g| atom_pair_args(g).use_2d)
    }
    #[setter(use2D)]
    fn set_use_2d(slf: PyRef<'_, Self>, v: bool) {
        slf.as_super()
            .with(slf.py(), |g| atom_pair_args(g).use_2d = v)
    }
    /// minimum distance to be included
    #[getter(minDistance)]
    fn get_min_distance(slf: PyRef<'_, Self>) -> u32 {
        slf.as_super()
            .with(slf.py(), |g| atom_pair_args(g).min_distance)
    }
    #[setter(minDistance)]
    fn set_min_distance(slf: PyRef<'_, Self>, v: u32) {
        slf.as_super()
            .with(slf.py(), |g| atom_pair_args(g).min_distance = v)
    }
    /// maximum distance to be included
    #[getter(maxDistance)]
    fn get_max_distance(slf: PyRef<'_, Self>) -> u32 {
        slf.as_super()
            .with(slf.py(), |g| atom_pair_args(g).max_distance)
    }
    #[setter(maxDistance)]
    fn set_max_distance(slf: PyRef<'_, Self>, v: u32) {
        slf.as_super()
            .with(slf.py(), |g| atom_pair_args(g).max_distance = v)
    }
}

fn torsion_args(g: &mut FingerprintGenerator) -> &mut TopologicalTorsionArguments {
    match &mut g.env {
        EnvGenerator::TopologicalTorsion(t) => t,
        _ => unreachable!("TopologicalTorsionFingerprintOptions of another generator"),
    }
}

/// `TopologicalTorsionFingerprintOptions`
#[pyclass(name = "TopologicalTorsionFingerprintOptions", module = "rdoxide.Chem.rdFingerprintGenerator", extends = PyFingerprintOptions)]
#[derive(Debug)]
pub struct PyTorsionOptions;

#[pymethods]
impl PyTorsionOptions {
    /// number of atoms to be included in the paths
    #[getter(torsionAtomCount)]
    fn get_torsion_atom_count(slf: PyRef<'_, Self>) -> u32 {
        slf.as_super()
            .with(slf.py(), |g| torsion_args(g).torsion_atom_count)
    }
    #[setter(torsionAtomCount)]
    fn set_torsion_atom_count(slf: PyRef<'_, Self>, v: u32) {
        slf.as_super()
            .with(slf.py(), |g| torsion_args(g).torsion_atom_count = v)
    }
    /// whether or not to only include paths which are the shortest path between the start and end atoms
    #[getter(onlyShortestPaths)]
    fn get_only_shortest_paths(slf: PyRef<'_, Self>) -> bool {
        slf.as_super()
            .with(slf.py(), |g| torsion_args(g).only_shortest_paths)
    }
    #[setter(onlyShortestPaths)]
    fn set_only_shortest_paths(slf: PyRef<'_, Self>, v: bool) {
        slf.as_super()
            .with(slf.py(), |g| torsion_args(g).only_shortest_paths = v)
    }
}

/// Get a topological torsion fingerprint generator
#[pyfunction]
#[pyo3(name = "GetTopologicalTorsionGenerator", signature = (includeChirality = false, torsionAtomCount = 4, countSimulation = true, countBounds = None, fpSize = 2048, atomInvariantsGenerator = None))]
#[allow(non_snake_case)]
fn get_topological_torsion_generator(
    includeChirality: bool,
    torsionAtomCount: u32,
    countSimulation: bool,
    countBounds: Option<&Bound<'_, PyAny>>,
    fpSize: u32,
    atomInvariantsGenerator: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyFPGen> {
    let count_bounds = count_bounds_arg(countBounds)?;
    if countSimulation && count_bounds.is_empty() {
        return Err(rdoxide::error::RDKitError::precondition(
            "bad count bounds provided",
        ))
        .py();
    }
    let atom_inv_gen = match atomInvariantsGenerator {
        Some(o) if !o.is_none() => o.extract::<PyRef<'_, PyAtomInvGen>>()?.gen_.clone(),
        _ => AtomInvariantsGenerator::AtomPair {
            include_chirality: includeChirality,
            topological_torsion_correction: true,
        },
    };
    Ok(PyFPGen {
        gen_: FingerprintGenerator {
            env: EnvGenerator::TopologicalTorsion(TopologicalTorsionArguments {
                torsion_atom_count: torsionAtomCount,
                only_shortest_paths: false,
            }),
            args: FingerprintArguments {
                count_simulation: countSimulation,
                include_chirality: includeChirality,
                count_bounds,
                fp_size: fpSize,
                num_bits_per_feature: 1,
            },
            atom_inv_gen: Some(atom_inv_gen),
            bond_inv_gen: None,
            output_64: true,
        },
    })
}

fn rdkit_fp_args(g: &mut FingerprintGenerator) -> &mut RDKitFPArguments {
    match &mut g.env {
        EnvGenerator::RDKitFP(r) => r,
        _ => unreachable!("RDKitFingerprintOptions of another generator"),
    }
}

/// `RDKitFingerprintOptions`
#[pyclass(name = "RDKitFingerprintOptions", module = "rdoxide.Chem.rdFingerprintGenerator", extends = PyFingerprintOptions)]
#[derive(Debug)]
pub struct PyRDKitFPOptions;

#[pymethods]
impl PyRDKitFPOptions {
    /// minimum path length (in bonds) to be included
    #[getter(minPath)]
    fn get_min_path(slf: PyRef<'_, Self>) -> u32 {
        slf.as_super().with(slf.py(), |g| rdkit_fp_args(g).min_path)
    }
    #[setter(minPath)]
    fn set_min_path(slf: PyRef<'_, Self>, v: u32) {
        slf.as_super()
            .with(slf.py(), |g| rdkit_fp_args(g).min_path = v)
    }
    /// maximum path length (in bonds) to be included
    #[getter(maxPath)]
    fn get_max_path(slf: PyRef<'_, Self>) -> u32 {
        slf.as_super().with(slf.py(), |g| rdkit_fp_args(g).max_path)
    }
    #[setter(maxPath)]
    fn set_max_path(slf: PyRef<'_, Self>, v: u32) {
        slf.as_super()
            .with(slf.py(), |g| rdkit_fp_args(g).max_path = v)
    }
    /// use explicit Hs in the paths (if molecule has explicit Hs)
    #[getter(useHs)]
    fn get_use_hs(slf: PyRef<'_, Self>) -> bool {
        slf.as_super().with(slf.py(), |g| rdkit_fp_args(g).use_hs)
    }
    #[setter(useHs)]
    fn set_use_hs(slf: PyRef<'_, Self>, v: bool) {
        slf.as_super()
            .with(slf.py(), |g| rdkit_fp_args(g).use_hs = v)
    }
    /// generate branched subgraphs, not just linear ones
    #[getter(branchedPaths)]
    fn get_branched_paths(slf: PyRef<'_, Self>) -> bool {
        slf.as_super()
            .with(slf.py(), |g| rdkit_fp_args(g).branched_paths)
    }
    #[setter(branchedPaths)]
    fn set_branched_paths(slf: PyRef<'_, Self>, v: bool) {
        slf.as_super()
            .with(slf.py(), |g| rdkit_fp_args(g).branched_paths = v)
    }
    /// include bond orders in the path hashes
    #[getter(useBondOrder)]
    fn get_use_bond_order(slf: PyRef<'_, Self>) -> bool {
        slf.as_super()
            .with(slf.py(), |g| rdkit_fp_args(g).use_bond_order)
    }
    #[setter(useBondOrder)]
    fn set_use_bond_order(slf: PyRef<'_, Self>, v: bool) {
        slf.as_super()
            .with(slf.py(), |g| rdkit_fp_args(g).use_bond_order = v)
    }
}

/// Get an RDKit fingerprint generator
#[pyfunction]
#[pyo3(name = "GetRDKitFPGenerator", signature = (minPath = 1, maxPath = 7, useHs = true, branchedPaths = true, useBondOrder = true, countSimulation = false, countBounds = None, fpSize = 2048, numBitsPerFeature = 2, atomInvariantsGenerator = None))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn get_rdkit_fp_generator(
    minPath: u32,
    maxPath: u32,
    useHs: bool,
    branchedPaths: bool,
    useBondOrder: bool,
    countSimulation: bool,
    countBounds: Option<&Bound<'_, PyAny>>,
    fpSize: u32,
    numBitsPerFeature: u32,
    atomInvariantsGenerator: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyFPGen> {
    let count_bounds = count_bounds_arg(countBounds)?;
    if countSimulation && count_bounds.is_empty() {
        return Err(rdoxide::error::RDKitError::precondition(
            "bad count bounds provided",
        ))
        .py();
    }
    if numBitsPerFeature == 0 {
        return Err(rdoxide::error::RDKitError::precondition(
            "numBitsPerFeature must be >0",
        ))
        .py();
    }
    if minPath == 0 {
        return Err(rdoxide::error::RDKitError::precondition("minPath==0")).py();
    }
    if maxPath < minPath {
        return Err(rdoxide::error::RDKitError::precondition("maxPath<minPath")).py();
    }
    let atom_inv_gen = match atomInvariantsGenerator {
        Some(o) if !o.is_none() => o.extract::<PyRef<'_, PyAtomInvGen>>()?.gen_.clone(),
        _ => AtomInvariantsGenerator::RDKitFP,
    };
    Ok(PyFPGen {
        gen_: FingerprintGenerator {
            env: EnvGenerator::RDKitFP(RDKitFPArguments {
                min_path: minPath,
                max_path: maxPath,
                use_hs: useHs,
                branched_paths: branchedPaths,
                use_bond_order: useBondOrder,
            }),
            args: FingerprintArguments {
                count_simulation: countSimulation,
                include_chirality: false,
                count_bounds,
                fp_size: fpSize,
                num_bits_per_feature: numBitsPerFeature,
            },
            atom_inv_gen: Some(atom_inv_gen),
            bond_inv_gen: None,
            output_64: true,
        },
    })
}

/// Get an RDKit atom invariants generator
#[pyfunction]
#[pyo3(name = "GetRDKitAtomInvGen")]
fn get_rdkit_atom_inv_gen() -> PyAtomInvGen {
    PyAtomInvGen {
        gen_: AtomInvariantsGenerator::RDKitFP,
    }
}

fn count_bounds_arg(o: Option<&Bound<'_, PyAny>>) -> PyResult<Vec<u32>> {
    match o {
        Some(o) if !o.is_none() => o.extract::<Vec<u32>>(),
        _ => Ok(vec![1, 2, 4, 8]),
    }
}

/// Get an atom pair fingerprint generator
#[pyfunction]
#[pyo3(name = "GetAtomPairGenerator", signature = (minDistance = 1, maxDistance = 30, includeChirality = false, use2D = true, countSimulation = true, countBounds = None, fpSize = 2048, atomInvariantsGenerator = None))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn get_atom_pair_generator(
    minDistance: u32,
    maxDistance: u32,
    includeChirality: bool,
    use2D: bool,
    countSimulation: bool,
    countBounds: Option<&Bound<'_, PyAny>>,
    fpSize: u32,
    atomInvariantsGenerator: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyFPGen> {
    let count_bounds = count_bounds_arg(countBounds)?;
    if countSimulation && count_bounds.is_empty() {
        return Err(rdoxide::error::RDKitError::precondition(
            "bad count bounds provided",
        ))
        .py();
    }
    if minDistance > maxDistance {
        return Err(rdoxide::error::RDKitError::precondition(
            "bad distances provided",
        ))
        .py();
    }
    let atom_inv_gen = match atomInvariantsGenerator {
        Some(o) if !o.is_none() => o.extract::<PyRef<'_, PyAtomInvGen>>()?.gen_.clone(),
        _ => AtomInvariantsGenerator::AtomPair {
            include_chirality: includeChirality,
            topological_torsion_correction: false,
        },
    };
    Ok(PyFPGen {
        gen_: FingerprintGenerator {
            env: EnvGenerator::AtomPair(AtomPairArguments {
                use_2d: use2D,
                min_distance: minDistance,
                max_distance: maxDistance,
            }),
            args: FingerprintArguments {
                count_simulation: countSimulation,
                include_chirality: includeChirality,
                count_bounds,
                fp_size: fpSize,
                num_bits_per_feature: 1,
            },
            atom_inv_gen: Some(atom_inv_gen),
            bond_inv_gen: None,
            output_64: true,
        },
    })
}

/// Get an atom pair atom-invariant generator
#[pyfunction]
#[pyo3(name = "GetAtomPairAtomInvGen", signature = (includeChirality = false))]
#[allow(non_snake_case)]
fn get_atom_pair_atom_inv_gen(includeChirality: bool) -> PyAtomInvGen {
    PyAtomInvGen {
        gen_: AtomInvariantsGenerator::AtomPair {
            include_chirality: includeChirality,
            topological_torsion_correction: false,
        },
    }
}

#[pymethods]
impl PyMorganOptions {
    /// use include atoms which have nonzero invariants
    #[getter(onlyNonzeroInvariants)]
    fn get_only_nonzero(slf: PyRef<'_, Self>) -> bool {
        slf.as_super()
            .with(slf.py(), |g| morgan_args(g).only_nonzero_invariants)
    }
    #[setter(onlyNonzeroInvariants)]
    fn set_only_nonzero(slf: PyRef<'_, Self>, v: bool) {
        slf.as_super()
            .with(slf.py(), |g| morgan_args(g).only_nonzero_invariants = v)
    }
    /// the radius of the fingerprints to generate
    #[getter(radius)]
    fn get_radius(slf: PyRef<'_, Self>) -> u32 {
        slf.as_super().with(slf.py(), |g| morgan_args(g).radius)
    }
    #[setter(radius)]
    fn set_radius(slf: PyRef<'_, Self>, v: u32) {
        slf.as_super().with(slf.py(), |g| morgan_args(g).radius = v)
    }
    /// include redundant environments in the fingerprint
    #[getter(includeRedundantEnvironments)]
    fn get_redundant(slf: PyRef<'_, Self>) -> bool {
        slf.as_super()
            .with(slf.py(), |g| morgan_args(g).include_redundant_environments)
    }
    #[setter(includeRedundantEnvironments)]
    fn set_redundant(slf: PyRef<'_, Self>, v: bool) {
        slf.as_super().with(slf.py(), |g| {
            morgan_args(g).include_redundant_environments = v
        })
    }
}

/// Get a morgan fingerprint generator
#[pyfunction]
#[pyo3(name = "GetMorganGenerator", signature = (radius = 3, countSimulation = false, includeChirality = false, useBondTypes = true, onlyNonzeroInvariants = false, includeRingMembership = true, countBounds = None, fpSize = 2048, atomInvariantsGenerator = None, bondInvariantsGenerator = None, includeRedundantEnvironments = false))]
#[allow(non_snake_case, clippy::too_many_arguments, unused_variables)]
fn get_morgan_generator(
    radius: u32,
    countSimulation: bool,
    includeChirality: bool,
    useBondTypes: bool,
    onlyNonzeroInvariants: bool,
    includeRingMembership: bool,
    countBounds: Option<&Bound<'_, PyAny>>,
    fpSize: u32,
    atomInvariantsGenerator: Option<&Bound<'_, PyAny>>,
    bondInvariantsGenerator: Option<&Bound<'_, PyAny>>,
    includeRedundantEnvironments: bool,
) -> PyResult<PyFPGen> {
    let count_bounds = match countBounds {
        Some(o) if !o.is_none() => o.extract::<Vec<u32>>()?,
        _ => vec![1, 2, 4, 8],
    };
    if countSimulation && count_bounds.is_empty() {
        return Err(rdoxide::error::RDKitError::precondition(
            "bad count bounds provided",
        ))
        .py();
    }
    // upstream ignores includeRingMembership here
    let atom_inv_gen = match atomInvariantsGenerator {
        Some(o) if !o.is_none() => o.extract::<PyRef<'_, PyAtomInvGen>>()?.gen_.clone(),
        _ => AtomInvariantsGenerator::Morgan {
            include_ring_membership: true,
        },
    };
    let bond_inv_gen = match bondInvariantsGenerator {
        Some(o) if !o.is_none() => o.extract::<PyRef<'_, PyBondInvGen>>()?.gen_.clone(),
        _ => BondInvariantsGenerator::Morgan {
            use_bond_types: useBondTypes,
            use_chirality: includeChirality,
        },
    };
    Ok(PyFPGen {
        gen_: FingerprintGenerator {
            env: EnvGenerator::Morgan(MorganArguments {
                radius,
                only_nonzero_invariants: onlyNonzeroInvariants,
                include_redundant_environments: includeRedundantEnvironments,
            }),
            args: FingerprintArguments {
                count_simulation: countSimulation,
                include_chirality: includeChirality,
                count_bounds,
                fp_size: fpSize,
                num_bits_per_feature: 1,
            },
            atom_inv_gen: Some(atom_inv_gen),
            bond_inv_gen: Some(bond_inv_gen),
            output_64: true,
        },
    })
}

/// Get a morgan atom invariants generator
#[pyfunction]
#[pyo3(name = "GetMorganAtomInvGen", signature = (includeRingMembership = true))]
#[allow(non_snake_case)]
fn get_morgan_atom_inv_gen(includeRingMembership: bool) -> PyAtomInvGen {
    PyAtomInvGen {
        gen_: AtomInvariantsGenerator::Morgan {
            include_ring_membership: includeRingMembership,
        },
    }
}

/// Get a morgan feature atom invariants generator
#[pyfunction]
#[pyo3(name = "GetMorganFeatureAtomInvGen", signature = (patterns = None))]
fn get_morgan_feature_atom_inv_gen(patterns: Option<&Bound<'_, PyAny>>) -> PyResult<PyAtomInvGen> {
    let pats = match patterns {
        // upstream silently ignores anything that is not a sequence of molecules
        Some(p) if !p.is_none() => p
            .try_iter()
            .and_then(|it| {
                it.map(|m| Ok(m?.extract::<PyRef<'_, PyMol>>()?.mol.clone()))
                    .collect::<PyResult<Vec<_>>>()
            })
            .ok(),
        _ => None,
    };
    Ok(PyAtomInvGen {
        gen_: AtomInvariantsGenerator::MorganFeature { patterns: pats },
    })
}

/// Get a morgan bond invariants generator
#[pyfunction]
#[pyo3(name = "GetMorganBondInvGen", signature = (useBondTypes = true, useChirality = false))]
#[allow(non_snake_case)]
fn get_morgan_bond_inv_gen(useBondTypes: bool, useChirality: bool) -> PyBondInvGen {
    PyBondInvGen {
        gen_: BondInvariantsGenerator::Morgan {
            use_bond_types: useBondTypes,
            use_chirality: useChirality,
        },
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyAtomInvGen>()?;
    m.add_class::<PyBondInvGen>()?;
    m.add_class::<PyAdditionalOutput>()?;
    m.add_class::<PyFPGen>()?;
    m.add_class::<PyFingerprintOptions>()?;
    m.add_class::<PyMorganOptions>()?;
    m.add_class::<PyAtomPairOptions>()?;
    m.add_class::<PyTorsionOptions>()?;
    m.add_class::<PyRDKitFPOptions>()?;
    add_enum_exported(
        m,
        "FPType",
        &[
            ("RDKitFP", FP_RDKIT),
            ("MorganFP", FP_MORGAN),
            ("AtomPairFP", FP_ATOM_PAIR),
            ("TopologicalTorsionFP", FP_TORSION),
        ],
        None,
    )?;
    m.add_function(wrap_pyfunction!(fingerprint_generator_from_json, m)?)?;
    m.add_function(wrap_pyfunction!(py_sparse_count_fps, m)?)?;
    m.add_function(wrap_pyfunction!(py_sparse_fps, m)?)?;
    m.add_function(wrap_pyfunction!(py_count_fps, m)?)?;
    m.add_function(wrap_pyfunction!(py_fps, m)?)?;
    m.add_function(wrap_pyfunction!(get_rdkit_fp_generator, m)?)?;
    m.add_function(wrap_pyfunction!(get_rdkit_atom_inv_gen, m)?)?;
    m.add_function(wrap_pyfunction!(get_topological_torsion_generator, m)?)?;
    m.add_function(wrap_pyfunction!(get_atom_pair_generator, m)?)?;
    m.add_function(wrap_pyfunction!(get_atom_pair_atom_inv_gen, m)?)?;
    m.add_function(wrap_pyfunction!(get_morgan_generator, m)?)?;
    m.add_function(wrap_pyfunction!(get_morgan_atom_inv_gen, m)?)?;
    m.add_function(wrap_pyfunction!(get_morgan_bond_inv_gen, m)?)?;
    m.add_function(wrap_pyfunction!(get_morgan_feature_atom_inv_gen, m)?)?;
    Ok(())
}
