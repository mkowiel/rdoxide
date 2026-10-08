//! `rdkit.Chem.rdMolDescriptors` (upstream
//! `GraphMol/Descriptors/Wrap/rdMolDescriptors.cpp`): the legacy fingerprint
//! functions so far.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rdoxide::graphmol::fingerprints::atom_pair::{CODE_SIZE, get_atom_code, get_atom_pair_code};
use rdoxide::graphmol::fingerprints::legacy::{self, BitInfoMap, MorganKind, MorganResult};
use rdoxide::graphmol::fingerprints::morgan::{
    get_connectivity_invariants, get_feature_invariants,
};

use crate::atom::PyAtom;
use crate::datastructs::PyExplicitBitVect;
use crate::datastructs_siv::{PyIntSIV, PyLongSIV, PyUIntSIV};
use crate::errors::IntoPyResult;
use crate::mol::PyMol;

fn deprecation_warning(msg: &str) {
    rdoxide::rd_log_warning!("DEPRECATION WARNING: {msg}");
}

/// `pythonObjectToVect(obj, maxV)`: `None` for falsy objects.
fn opt_vect(o: Option<&Bound<'_, PyAny>>, max_v: u32) -> PyResult<Option<Vec<u32>>> {
    let Some(o) = o else { return Ok(None) };
    if !o.is_truthy()? {
        return Ok(None);
    }
    let mut v = Vec::new();
    for x in o.try_iter()? {
        let x: u32 = x?.extract()?;
        if x >= max_v {
            return Err(PyValueError::new_err(
                "list element larger than allowed value",
            ));
        }
        v.push(x);
    }
    Ok(Some(v))
}

/// The Morgan helpers' handling of `invariants` / `useFeatures`.
fn morgan_invariants(
    mol: &PyMol,
    invariants: Option<&Bound<'_, PyAny>>,
    use_features: bool,
) -> PyResult<Option<Vec<u32>>> {
    if let Some(inv) = invariants
        && inv.is_truthy()?
    {
        let n = inv.len()?;
        if n == 0 {
            return Ok(None);
        }
        if n != mol.mol.num_atoms() as usize {
            return Err(PyValueError::new_err(
                "length of invariant vector != number of atoms",
            ));
        }
        return Ok(Some(
            (0..n)
                .map(|i| inv.get_item(i)?.extract::<u32>())
                .collect::<PyResult<_>>()?,
        ));
    }
    if use_features {
        return Ok(Some(get_feature_invariants(&mol.mol, None).py()?));
    }
    Ok(None)
}

fn fill_bit_info(py: Python<'_>, d: Option<&Bound<'_, PyAny>>, bi: &BitInfoMap) -> PyResult<()> {
    if let Some(d) = d {
        let d = d.cast::<PyDict>()?;
        d.clear();
        for (k, v) in bi {
            let items = v
                .iter()
                .map(|&(a, r)| PyTuple::new(py, [a, r]))
                .collect::<PyResult<Vec<_>>>()?;
            d.set_item(k, PyTuple::new(py, items)?)?;
        }
    }
    Ok(())
}

fn requested<'a, 'py>(o: Option<&'a Bound<'py, PyAny>>) -> Option<&'a Bound<'py, PyAny>> {
    o.filter(|o| !o.is_none())
}

#[allow(clippy::too_many_arguments)]
fn morgan_helper(
    py: Python<'_>,
    mol: &PyMol,
    radius: u32,
    kind: MorganKind,
    invariants: Option<&Bound<'_, PyAny>>,
    from_atoms: Option<&Bound<'_, PyAny>>,
    use_chirality: bool,
    use_bond_types: bool,
    use_features: bool,
    bit_info: Option<&Bound<'_, PyAny>>,
    include_redundant: bool,
) -> PyResult<MorganResult> {
    deprecation_warning("please use MorganGenerator");
    let invs = morgan_invariants(mol, invariants, use_features)?;
    let froms = match kind {
        // GetMorganFingerprintAsBitVect range-checks fromAtoms
        MorganKind::BitVect(_) => opt_vect(from_atoms, mol.mol.num_atoms())?,
        _ => opt_vect(from_atoms, u32::MAX)?,
    };
    let bit_info = requested(bit_info);
    if let Some(b) = bit_info {
        b.cast::<PyDict>()?;
    }
    let mut bi = BitInfoMap::new();
    let res = legacy::morgan_fingerprint(
        &mol.mol,
        radius,
        kind,
        invs.as_deref(),
        froms.as_deref(),
        use_chirality,
        use_bond_types,
        bit_info.is_some().then_some(&mut bi),
        include_redundant,
    )
    .py()?;
    fill_bit_info(py, bit_info, &bi)?;
    Ok(res)
}

fn counts(r: MorganResult) -> PyUIntSIV {
    match r {
        MorganResult::Counts(c) => PyUIntSIV::wrap(c),
        MorganResult::Bits(_) => unreachable!("count result expected"),
    }
}

/// Returns a Morgan fingerprint for a molecule
#[pyfunction]
#[pyo3(name = "GetMorganFingerprint", signature = (mol, radius, invariants = None, fromAtoms = None, useChirality = false, useBondTypes = true, useFeatures = false, useCounts = true, bitInfo = None, includeRedundantEnvironments = false))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn get_morgan_fingerprint(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    radius: u32,
    invariants: Option<&Bound<'_, PyAny>>,
    fromAtoms: Option<&Bound<'_, PyAny>>,
    useChirality: bool,
    useBondTypes: bool,
    useFeatures: bool,
    useCounts: bool,
    bitInfo: Option<&Bound<'_, PyAny>>,
    includeRedundantEnvironments: bool,
) -> PyResult<PyUIntSIV> {
    let r = morgan_helper(
        py,
        &mol,
        radius,
        MorganKind::Unfolded {
            use_counts: useCounts,
        },
        invariants,
        fromAtoms,
        useChirality,
        useBondTypes,
        useFeatures,
        bitInfo,
        includeRedundantEnvironments,
    )?;
    Ok(counts(r))
}

/// Returns a hashed Morgan fingerprint for a molecule
#[pyfunction]
#[pyo3(name = "GetHashedMorganFingerprint", signature = (mol, radius, nBits = 2048, invariants = None, fromAtoms = None, useChirality = false, useBondTypes = true, useFeatures = false, bitInfo = None, includeRedundantEnvironments = false))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn get_hashed_morgan_fingerprint(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    radius: u32,
    nBits: i64,
    invariants: Option<&Bound<'_, PyAny>>,
    fromAtoms: Option<&Bound<'_, PyAny>>,
    useChirality: bool,
    useBondTypes: bool,
    useFeatures: bool,
    bitInfo: Option<&Bound<'_, PyAny>>,
    includeRedundantEnvironments: bool,
) -> PyResult<PyUIntSIV> {
    let kind = if nBits < 0 {
        MorganKind::Unfolded { use_counts: true }
    } else {
        MorganKind::Hashed(nBits as u32)
    };
    let r = morgan_helper(
        py,
        &mol,
        radius,
        kind,
        invariants,
        fromAtoms,
        useChirality,
        useBondTypes,
        useFeatures,
        bitInfo,
        includeRedundantEnvironments,
    )?;
    Ok(counts(r))
}

/// Returns a Morgan fingerprint for a molecule as a bit vector
#[pyfunction]
#[pyo3(name = "GetMorganFingerprintAsBitVect", signature = (mol, radius, nBits = 2048, invariants = None, fromAtoms = None, useChirality = false, useBondTypes = true, useFeatures = false, bitInfo = None, includeRedundantEnvironments = false))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn get_morgan_fingerprint_bv(
    py: Python<'_>,
    mol: PyRef<'_, PyMol>,
    radius: u32,
    nBits: u32,
    invariants: Option<&Bound<'_, PyAny>>,
    fromAtoms: Option<&Bound<'_, PyAny>>,
    useChirality: bool,
    useBondTypes: bool,
    useFeatures: bool,
    bitInfo: Option<&Bound<'_, PyAny>>,
    includeRedundantEnvironments: bool,
) -> PyResult<PyExplicitBitVect> {
    match morgan_helper(
        py,
        &mol,
        radius,
        MorganKind::BitVect(nBits),
        invariants,
        fromAtoms,
        useChirality,
        useBondTypes,
        useFeatures,
        bitInfo,
        includeRedundantEnvironments,
    )? {
        MorganResult::Bits(b) => Ok(PyExplicitBitVect::wrap(b)),
        MorganResult::Counts(_) => unreachable!("bit vector result expected"),
    }
}

/// Returns connectivity invariants (ECFP-like) for a molecule.
#[pyfunction]
#[pyo3(name = "GetConnectivityInvariants", signature = (mol, includeRingMembership = true))]
#[allow(non_snake_case)]
fn get_connectivity_invariants_py<'py>(
    py: Python<'py>,
    mol: PyRef<'_, PyMol>,
    includeRingMembership: bool,
) -> PyResult<Bound<'py, PyList>> {
    PyList::new(
        py,
        get_connectivity_invariants(&mol.mol, includeRingMembership).py()?,
    )
}

/// Returns feature invariants (FCFP-like) for a molecule.
#[pyfunction]
#[pyo3(name = "GetFeatureInvariants", signature = (mol))]
fn get_feature_invariants_py<'py>(
    py: Python<'py>,
    mol: PyRef<'_, PyMol>,
) -> PyResult<Bound<'py, PyList>> {
    PyList::new(py, get_feature_invariants(&mol.mol, None).py()?)
}

/// Returns the atom code (hash) for an atom
#[pyfunction]
#[pyo3(name = "GetAtomPairAtomCode", signature = (atom, branchSubtract = 0, includeChirality = false))]
#[allow(non_snake_case)]
fn get_atom_pair_atom_code(
    atom: PyRef<'_, PyAtom>,
    branchSubtract: u32,
    includeChirality: bool,
) -> PyResult<u32> {
    atom.with_owner(|mol, idx| get_atom_code(mol, idx, branchSubtract, includeChirality).py())
}

/// Returns the atom-pair code (hash) for a pair of atoms separated by a certain number of bonds
#[pyfunction]
#[pyo3(name = "GetAtomPairCode", signature = (atom1Code, atom2Code, distance, includeChirality = false))]
#[allow(non_snake_case)]
fn get_atom_pair_code_py(
    atom1Code: u32,
    atom2Code: u32,
    distance: u32,
    includeChirality: bool,
) -> PyResult<u32> {
    get_atom_pair_code(atom1Code, atom2Code, distance, includeChirality).py()
}

struct ApArgs {
    from: Option<Vec<u32>>,
    ignore: Option<Vec<u32>>,
    invs: Option<Vec<u32>>,
}

fn ap_args(
    mol: &PyMol,
    from: Option<&Bound<'_, PyAny>>,
    ignore: Option<&Bound<'_, PyAny>>,
    invs: Option<&Bound<'_, PyAny>>,
) -> PyResult<ApArgs> {
    let n = mol.mol.num_atoms();
    Ok(ApArgs {
        from: opt_vect(from, n)?,
        ignore: opt_vect(ignore, n)?,
        invs: opt_vect(invs, 1 << CODE_SIZE)?,
    })
}

/// Returns the atom-pair fingerprint for a molecule as an IntSparseIntVect
#[pyfunction]
#[pyo3(name = "GetAtomPairFingerprint", signature = (mol, minLength = 1, maxLength = 30, fromAtoms = None, ignoreAtoms = None, atomInvariants = None, includeChirality = false, use2D = true, confId = -1))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn get_atom_pair_fingerprint(
    mol: PyRef<'_, PyMol>,
    minLength: u32,
    maxLength: u32,
    fromAtoms: Option<&Bound<'_, PyAny>>,
    ignoreAtoms: Option<&Bound<'_, PyAny>>,
    atomInvariants: Option<&Bound<'_, PyAny>>,
    includeChirality: bool,
    use2D: bool,
    confId: i32,
) -> PyResult<PyIntSIV> {
    let a = ap_args(&mol, fromAtoms, ignoreAtoms, atomInvariants)?;
    deprecation_warning("please use AtomPairGenerator");
    let r = legacy::get_atom_pair_fingerprint(
        &mol.mol,
        minLength,
        maxLength,
        a.from.as_deref(),
        a.ignore.as_deref(),
        a.invs.as_deref(),
        includeChirality,
        use2D,
        confId,
    )
    .py()?;
    Ok(PyIntSIV::wrap(r))
}

/// Returns the hashed atom-pair fingerprint for a molecule as an IntSparseIntVect
#[pyfunction]
#[pyo3(name = "GetHashedAtomPairFingerprint", signature = (mol, nBits = 2048, minLength = 1, maxLength = 30, fromAtoms = None, ignoreAtoms = None, atomInvariants = None, includeChirality = false, use2D = true, confId = -1))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn get_hashed_atom_pair_fingerprint(
    mol: PyRef<'_, PyMol>,
    nBits: u32,
    minLength: u32,
    maxLength: u32,
    fromAtoms: Option<&Bound<'_, PyAny>>,
    ignoreAtoms: Option<&Bound<'_, PyAny>>,
    atomInvariants: Option<&Bound<'_, PyAny>>,
    includeChirality: bool,
    use2D: bool,
    confId: i32,
) -> PyResult<PyIntSIV> {
    let a = ap_args(&mol, fromAtoms, ignoreAtoms, atomInvariants)?;
    deprecation_warning("please use AtomPairGenerator");
    let r = legacy::get_hashed_atom_pair_fingerprint(
        &mol.mol,
        nBits,
        minLength,
        maxLength,
        a.from.as_deref(),
        a.ignore.as_deref(),
        a.invs.as_deref(),
        includeChirality,
        use2D,
        confId,
    )
    .py()?;
    Ok(PyIntSIV::wrap(r))
}

/// Returns the atom-pair fingerprint for a molecule as an ExplicitBitVect
#[pyfunction]
#[pyo3(name = "GetHashedAtomPairFingerprintAsBitVect", signature = (mol, nBits = 2048, minLength = 1, maxLength = 30, fromAtoms = None, ignoreAtoms = None, atomInvariants = None, nBitsPerEntry = 4, includeChirality = false, use2D = true, confId = -1))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn get_hashed_atom_pair_fingerprint_bv(
    mol: PyRef<'_, PyMol>,
    nBits: u32,
    minLength: u32,
    maxLength: u32,
    fromAtoms: Option<&Bound<'_, PyAny>>,
    ignoreAtoms: Option<&Bound<'_, PyAny>>,
    atomInvariants: Option<&Bound<'_, PyAny>>,
    nBitsPerEntry: u32,
    includeChirality: bool,
    use2D: bool,
    confId: i32,
) -> PyResult<PyExplicitBitVect> {
    let a = ap_args(&mol, fromAtoms, ignoreAtoms, atomInvariants)?;
    if minLength <= maxLength {
        deprecation_warning("please use AtomPairGenerator");
    }
    let r = legacy::get_hashed_atom_pair_fingerprint_as_bit_vect(
        &mol.mol,
        nBits,
        minLength,
        maxLength,
        a.from.as_deref(),
        a.ignore.as_deref(),
        a.invs.as_deref(),
        nBitsPerEntry,
        includeChirality,
        use2D,
        confId,
    )
    .py()?;
    Ok(PyExplicitBitVect::wrap(r))
}

/// Returns the topological-torsion fingerprint for a molecule as a LongIntSparseIntVect
#[pyfunction]
#[pyo3(name = "GetTopologicalTorsionFingerprint", signature = (mol, targetSize = 4, fromAtoms = None, ignoreAtoms = None, atomInvariants = None, includeChirality = false))]
#[allow(non_snake_case)]
fn get_topological_torsion_fingerprint(
    mol: PyRef<'_, PyMol>,
    targetSize: u32,
    fromAtoms: Option<&Bound<'_, PyAny>>,
    ignoreAtoms: Option<&Bound<'_, PyAny>>,
    atomInvariants: Option<&Bound<'_, PyAny>>,
    includeChirality: bool,
) -> PyResult<PyLongSIV> {
    let a = ap_args(&mol, fromAtoms, ignoreAtoms, atomInvariants)?;
    if targetSize * CODE_SIZE > 64 {
        return Err(PyValueError::new_err(format!(
            "Maximum supported topological torsion path length is {}\n",
            64 / CODE_SIZE
        )));
    }
    deprecation_warning("please use TopologicalTorsionGenerator");
    let r = legacy::get_topological_torsion_fingerprint(
        &mol.mol,
        targetSize,
        a.from.as_deref(),
        a.ignore.as_deref(),
        a.invs.as_deref(),
        includeChirality,
    )
    .py()?;
    Ok(PyLongSIV::wrap(r))
}

/// Returns the hashed topological-torsion fingerprint for a molecule as a LongIntSparseIntVect
#[pyfunction]
#[pyo3(name = "GetHashedTopologicalTorsionFingerprint", signature = (mol, nBits = 2048, targetSize = 4, fromAtoms = None, ignoreAtoms = None, atomInvariants = None, includeChirality = false))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn get_hashed_topological_torsion_fingerprint(
    mol: PyRef<'_, PyMol>,
    nBits: u32,
    targetSize: u32,
    fromAtoms: Option<&Bound<'_, PyAny>>,
    ignoreAtoms: Option<&Bound<'_, PyAny>>,
    atomInvariants: Option<&Bound<'_, PyAny>>,
    includeChirality: bool,
) -> PyResult<PyLongSIV> {
    let a = ap_args(&mol, fromAtoms, ignoreAtoms, atomInvariants)?;
    deprecation_warning("please use TopologicalTorsionGenerator");
    let r = legacy::get_hashed_topological_torsion_fingerprint(
        &mol.mol,
        nBits,
        targetSize,
        a.from.as_deref(),
        a.ignore.as_deref(),
        a.invs.as_deref(),
        includeChirality,
    )
    .py()?;
    Ok(PyLongSIV::wrap(r))
}

/// Returns the topological-torsion fingerprint for a molecule as an ExplicitBitVect
#[pyfunction]
#[pyo3(name = "GetHashedTopologicalTorsionFingerprintAsBitVect", signature = (mol, nBits = 2048, targetSize = 4, fromAtoms = None, ignoreAtoms = None, atomInvariants = None, nBitsPerEntry = 4, includeChirality = false))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn get_hashed_topological_torsion_fingerprint_bv(
    mol: PyRef<'_, PyMol>,
    nBits: u32,
    targetSize: u32,
    fromAtoms: Option<&Bound<'_, PyAny>>,
    ignoreAtoms: Option<&Bound<'_, PyAny>>,
    atomInvariants: Option<&Bound<'_, PyAny>>,
    nBitsPerEntry: u32,
    includeChirality: bool,
) -> PyResult<PyExplicitBitVect> {
    let a = ap_args(&mol, fromAtoms, ignoreAtoms, atomInvariants)?;
    deprecation_warning("please use TopologicalTorsionGenerator");
    let r = legacy::get_hashed_topological_torsion_fingerprint_as_bit_vect(
        &mol.mol,
        nBits,
        targetSize,
        a.from.as_deref(),
        a.ignore.as_deref(),
        a.invs.as_deref(),
        nBitsPerEntry,
        includeChirality,
    )
    .py()?;
    Ok(PyExplicitBitVect::wrap(r))
}

/// Returns the MACCS keys for a molecule as an ExplicitBitVect
#[pyfunction]
#[pyo3(name = "GetMACCSKeysFingerprint", signature = (mol))]
fn get_maccs_keys_fingerprint(mol: PyRef<'_, PyMol>) -> PyResult<PyExplicitBitVect> {
    Ok(PyExplicitBitVect::wrap(
        rdoxide::graphmol::fingerprints::maccs::get_fingerprint_as_bit_vect(&mol.mol).py()?,
    ))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions to compute molecular descriptors",
    )?;
    m.add("_MorganFingerprint_version", "1.0.0")?;
    m.add_function(wrap_pyfunction!(get_morgan_fingerprint, m)?)?;
    m.add_function(wrap_pyfunction!(get_maccs_keys_fingerprint, m)?)?;
    m.add_function(wrap_pyfunction!(get_hashed_morgan_fingerprint, m)?)?;
    m.add_function(wrap_pyfunction!(get_morgan_fingerprint_bv, m)?)?;
    m.add_function(wrap_pyfunction!(get_connectivity_invariants_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_feature_invariants_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_atom_pair_atom_code, m)?)?;
    m.add_function(wrap_pyfunction!(get_atom_pair_code_py, m)?)?;
    m.add_function(wrap_pyfunction!(get_atom_pair_fingerprint, m)?)?;
    m.add_function(wrap_pyfunction!(get_hashed_atom_pair_fingerprint, m)?)?;
    m.add_function(wrap_pyfunction!(get_hashed_atom_pair_fingerprint_bv, m)?)?;
    m.add_function(wrap_pyfunction!(get_topological_torsion_fingerprint, m)?)?;
    m.add_function(wrap_pyfunction!(
        get_hashed_topological_torsion_fingerprint,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(
        get_hashed_topological_torsion_fingerprint_bv,
        m
    )?)?;
    Ok(())
}
