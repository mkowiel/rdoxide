//! `rdkit.Chem.rdMolDescriptors` (upstream
//! `GraphMol/Descriptors/Wrap/rdMolDescriptors.cpp`): the legacy fingerprint
//! functions so far.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rdoxide::graphmol::descriptors::crippen::{
    CRIPPEN_VERSION, calc_crippen_descriptors, get_crippen_atom_contribs,
};
use rdoxide::graphmol::descriptors::lipinski::{self as lip, NumRotatableBondsOptions};
use rdoxide::graphmol::fingerprints::atom_pair::{CODE_SIZE, get_atom_code, get_atom_pair_code};
use rdoxide::graphmol::fingerprints::legacy::{self, BitInfoMap, MorganKind, MorganResult};
use rdoxide::graphmol::fingerprints::morgan::{
    get_connectivity_invariants, get_feature_invariants,
};

use crate::atom::PyAtom;
use crate::datastructs::PyExplicitBitVect;
use crate::datastructs_siv::{PyIntSIV, PyLongSIV, PyUIntSIV};
use crate::enums::add_enum_exported;
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

/// `_CalcCrippenContribs(mol, force=False, atomTypes=[], atomTypeLabels=[])`: per-atom
/// (logP, MR) contributions. Non-empty `atomTypes` / `atomTypeLabels` lists (of the
/// molecule's length) receive the matched parameter index and label of each atom.
#[pyfunction]
#[pyo3(
    name = "_CalcCrippenContribs",
    signature = (mol, force = None, atomTypes = None, atomTypeLabels = None)
)]
#[allow(non_snake_case)]
fn calc_crippen_contribs<'py>(
    py: Python<'py>,
    mol: &Bound<'py, PyMol>,
    force: Option<&Bound<'py, PyAny>>,
    atomTypes: Option<Bound<'py, PyList>>,
    atomTypeLabels: Option<Bound<'py, PyList>>,
) -> PyResult<Bound<'py, PyList>> {
    let force = bool_arg(force, false)?;
    let mut m = mol.borrow_mut();
    let n = m.mol.num_atoms() as usize;
    let n_types = atomTypes.as_ref().map_or(0, |l| l.len());
    if n_types != 0 && n_types != n {
        return Err(PyValueError::new_err(
            "if atomTypes vector is provided, it must be as long as the number of atoms",
        ));
    }
    let n_labels = atomTypeLabels.as_ref().map_or(0, |l| l.len());
    if n_labels != 0 && n_labels != n {
        return Err(PyValueError::new_err(
            "if atomTypeLabels vector is provided, it must be as long as the number of atoms",
        ));
    }
    let c = get_crippen_atom_contribs(&mut m.mol, force).py()?;
    let contribs = PyList::empty(py);
    for i in 0..n {
        contribs.append((c.logp[i], c.mr[i]))?;
    }
    if let Some(types) = &atomTypes
        && n_types != 0
    {
        for (i, t) in c.atom_types.iter().enumerate() {
            types.set_item(i, *t)?;
        }
    }
    if let Some(labels) = &atomTypeLabels
        && n_labels != 0
    {
        for (i, l) in c.atom_type_labels.iter().enumerate() {
            labels.set_item(i, l)?;
        }
    }
    Ok(contribs)
}

/// A bool-typed argument with a default. Upstream's Boost.Python bool converter
/// accepts any object, converting it by truthiness (`includeHs=1` is `True`).
pub(crate) fn bool_arg(v: Option<&Bound<'_, PyAny>>, default: bool) -> PyResult<bool> {
    v.map_or(Ok(default), |v| v.is_truthy())
}

/// `CalcCrippenDescriptors(mol, includeHs=True, force=False)`: the Wildman-Crippen
/// (logP, MR) pair of a molecule.
#[pyfunction]
#[pyo3(
    name = "CalcCrippenDescriptors",
    signature = (mol, includeHs = None, force = None)
)]
#[allow(non_snake_case)]
fn calc_crippen_descriptors_py(
    mol: &Bound<'_, PyMol>,
    includeHs: Option<&Bound<'_, PyAny>>,
    force: Option<&Bound<'_, PyAny>>,
) -> PyResult<(f64, f64)> {
    let include_hs = bool_arg(includeHs, true)?;
    let force = bool_arg(force, false)?;
    let mut m = mol.borrow_mut();
    calc_crippen_descriptors(&mut m.mol, include_hs, force).py()
}

/// Defines `(mol) -> count` functions over a shared molecule borrow.
macro_rules! mol_count_fn {
    ($fname:ident, $pyname:literal, $ret:ty, |$m:ident| $body:expr) => {
        #[pyfunction]
        #[pyo3(name = $pyname, signature = (mol))]
        #[allow(non_snake_case)]
        fn $fname(mol: &Bound<'_, PyMol>) -> PyResult<$ret> {
            let mol_ref = mol.borrow();
            let $m = &mol_ref.mol;
            $body
        }
    };
}

mol_count_fn!(calc_num_lipinski_hba, "CalcNumLipinskiHBA", u32, |m| Ok(
    lip::calc_lipinski_hba(m)
));
mol_count_fn!(calc_num_lipinski_hbd, "CalcNumLipinskiHBD", u32, |m| {
    lip::calc_lipinski_hbd(m).py()
});
mol_count_fn!(calc_num_hbd, "CalcNumHBD", u32, |m| lip::calc_num_hbd(m)
    .py());
mol_count_fn!(calc_num_hba, "CalcNumHBA", u32, |m| lip::calc_num_hba(m)
    .py());
mol_count_fn!(calc_num_heteroatoms, "CalcNumHeteroatoms", u32, |m| {
    lip::calc_num_heteroatoms(m).py()
});
mol_count_fn!(calc_num_amide_bonds, "CalcNumAmideBonds", u32, |m| {
    lip::calc_num_amide_bonds(m).py()
});
mol_count_fn!(calc_num_rings, "CalcNumRings", u32, |m| {
    lip::calc_num_rings(m).py()
});
mol_count_fn!(calc_fraction_csp3, "CalcFractionCSP3", f64, |m| {
    lip::calc_fraction_csp3(m).py()
});
mol_count_fn!(calc_num_heterocycles, "CalcNumHeterocycles", u32, |m| Ok(
    lip::calc_num_heterocycles(m)
));
mol_count_fn!(
    calc_num_aromatic_rings,
    "CalcNumAromaticRings",
    u32,
    |m| Ok(lip::calc_num_aromatic_rings(m))
);
mol_count_fn!(
    calc_num_saturated_rings,
    "CalcNumSaturatedRings",
    u32,
    |m| Ok(lip::calc_num_saturated_rings(m))
);
mol_count_fn!(
    calc_num_aliphatic_rings,
    "CalcNumAliphaticRings",
    u32,
    |m| Ok(lip::calc_num_aliphatic_rings(m))
);
mol_count_fn!(
    calc_num_aromatic_heterocycles,
    "CalcNumAromaticHeterocycles",
    u32,
    |m| Ok(lip::calc_num_aromatic_heterocycles(m))
);
mol_count_fn!(
    calc_num_aromatic_carbocycles,
    "CalcNumAromaticCarbocycles",
    u32,
    |m| Ok(lip::calc_num_aromatic_carbocycles(m))
);
mol_count_fn!(
    calc_num_aliphatic_heterocycles,
    "CalcNumAliphaticHeterocycles",
    u32,
    |m| Ok(lip::calc_num_aliphatic_heterocycles(m))
);
mol_count_fn!(
    calc_num_aliphatic_carbocycles,
    "CalcNumAliphaticCarbocycles",
    u32,
    |m| Ok(lip::calc_num_aliphatic_carbocycles(m))
);
mol_count_fn!(
    calc_num_saturated_heterocycles,
    "CalcNumSaturatedHeterocycles",
    u32,
    |m| Ok(lip::calc_num_saturated_heterocycles(m))
);
mol_count_fn!(
    calc_num_saturated_carbocycles,
    "CalcNumSaturatedCarbocycles",
    u32,
    |m| Ok(lip::calc_num_saturated_carbocycles(m))
);
mol_count_fn!(
    calc_num_atom_stereo_centers,
    "CalcNumAtomStereoCenters",
    u32,
    |m| lip::num_atom_stereo_centers(m).py()
);
mol_count_fn!(
    calc_num_unspecified_atom_stereo_centers,
    "CalcNumUnspecifiedAtomStereoCenters",
    u32,
    |m| lip::num_unspecified_atom_stereo_centers(m).py()
);

/// `CalcNumRotatableBonds(mol, strict=Default)`. Upstream tries the bool overload
/// first, so a plain value is converted by truthiness (`True` is `Strict`); only
/// `NumRotatableBondsOptions` members select the other definitions.
#[pyfunction]
#[pyo3(name = "CalcNumRotatableBonds", signature = (mol, strict = None))]
fn calc_num_rotatable_bonds_py(
    mol: &Bound<'_, PyMol>,
    strict: Option<&Bound<'_, PyAny>>,
) -> PyResult<u32> {
    let opt = match strict {
        None => NumRotatableBondsOptions::Default,
        Some(v) if v.get_type().name()?.to_cow()? != "NumRotatableBondsOptions" => {
            if v.is_truthy()? {
                NumRotatableBondsOptions::Strict
            } else {
                NumRotatableBondsOptions::NonStrict
            }
        }
        Some(v) => match v.extract::<i32>()? {
            -1 => NumRotatableBondsOptions::Default,
            0 => NumRotatableBondsOptions::NonStrict,
            1 => NumRotatableBondsOptions::Strict,
            2 => NumRotatableBondsOptions::StrictLinkages,
            // any other value takes the final branch of upstream's definition switch
            _ => NumRotatableBondsOptions::StrictLinkages,
        },
    };
    lip::calc_num_rotatable_bonds(&mol.borrow().mol, opt).py()
}

/// `CalcNumSpiroAtoms(mol, atoms=None)`: appends the spiro atoms to `atoms` when given.
#[pyfunction]
#[pyo3(name = "CalcNumSpiroAtoms", signature = (mol, atoms = None))]
fn calc_num_spiro_atoms_py(
    mol: &Bound<'_, PyMol>,
    atoms: Option<Bound<'_, PyList>>,
) -> PyResult<u32> {
    let found = lip::spiro_atoms(&mut mol.borrow_mut().mol).py()?;
    if let Some(list) = atoms {
        for a in &found {
            list.append(*a)?;
        }
    }
    Ok(found.len() as u32)
}

/// `CalcNumBridgeheadAtoms(mol, atoms=None)`: appends the bridgehead atoms to `atoms` when given.
#[pyfunction]
#[pyo3(name = "CalcNumBridgeheadAtoms", signature = (mol, atoms = None))]
fn calc_num_bridgehead_atoms_py(
    mol: &Bound<'_, PyMol>,
    atoms: Option<Bound<'_, PyList>>,
) -> PyResult<u32> {
    let found = lip::bridgehead_atoms(&mut mol.borrow_mut().mol).py()?;
    if let Some(list) = atoms {
        for a in &found {
            list.append(*a)?;
        }
    }
    Ok(found.len() as u32)
}

/// Version strings of the Lipinski-style descriptors (`_<Name>_version`).
const LIPINSKI_VERSIONS: &[(&str, &str)] = &[
    ("_CalcNumLipinskiHBA_version", lip::LIPINSKI_HBA_VERSION),
    ("_CalcNumLipinskiHBD_version", lip::LIPINSKI_HBD_VERSION),
    (
        "_CalcNumRotatableBonds_version",
        lip::NUM_ROTATABLE_BONDS_VERSION,
    ),
    ("_CalcNumHBD_version", lip::NUM_HBD_VERSION),
    ("_CalcNumHBA_version", lip::NUM_HBA_VERSION),
    ("_CalcNumHeteroatoms_version", lip::NUM_HETEROATOMS_VERSION),
    ("_CalcNumAmideBonds_version", lip::NUM_AMIDE_BONDS_VERSION),
    ("_CalcFractionCSP3_version", lip::FRACTION_CSP3_VERSION),
    ("_CalcNumRings_version", lip::NUM_RINGS_VERSION),
    (
        "_CalcNumHeterocycles_version",
        lip::NUM_HETEROCYCLES_VERSION,
    ),
    (
        "_CalcNumAromaticRings_version",
        lip::NUM_AROMATIC_RINGS_VERSION,
    ),
    (
        "_CalcNumSaturatedRings_version",
        lip::NUM_SATURATED_RINGS_VERSION,
    ),
    (
        "_CalcNumAliphaticRings_version",
        lip::NUM_ALIPHATIC_RINGS_VERSION,
    ),
    (
        "_CalcNumAromaticHeterocycles_version",
        lip::NUM_AROMATIC_HETEROCYCLES_VERSION,
    ),
    (
        "_CalcNumAromaticCarbocycles_version",
        lip::NUM_AROMATIC_CARBOCYCLES_VERSION,
    ),
    (
        "_CalcNumSaturatedHeterocycles_version",
        lip::NUM_SATURATED_HETEROCYCLES_VERSION,
    ),
    (
        "_CalcNumSaturatedCarbocycles_version",
        lip::NUM_SATURATED_CARBOCYCLES_VERSION,
    ),
    (
        "_CalcNumAliphaticHeterocycles_version",
        lip::NUM_ALIPHATIC_HETEROCYCLES_VERSION,
    ),
    (
        "_CalcNumAliphaticCarbocycles_version",
        lip::NUM_ALIPHATIC_CARBOCYCLES_VERSION,
    ),
    ("_CalcNumSpiroAtoms_version", lip::NUM_SPIRO_ATOMS_VERSION),
    (
        "_CalcNumBridgeheadAtoms_version",
        lip::NUM_BRIDGEHEAD_ATOMS_VERSION,
    ),
    (
        "_CalcNumAtomStereoCenters_version",
        lip::NUM_ATOM_STEREO_CENTERS_VERSION,
    ),
    (
        "_CalcNumUnspecifiedAtomStereoCenters_version",
        lip::NUM_UNSPECIFIED_ATOM_STEREO_CENTERS_VERSION,
    ),
];

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions to compute molecular descriptors",
    )?;
    m.add("_MorganFingerprint_version", "1.0.0")?;
    m.add_function(wrap_pyfunction!(get_morgan_fingerprint, m)?)?;
    m.add_function(wrap_pyfunction!(get_maccs_keys_fingerprint, m)?)?;
    m.add_function(wrap_pyfunction!(calc_crippen_contribs, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_lipinski_hba, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_lipinski_hbd, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_rotatable_bonds_py, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_hbd, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_hba, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_heteroatoms, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_amide_bonds, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_rings, m)?)?;
    crate::rdmoldescriptors_props::init(m)?;
    crate::rdmoldescriptors_molsurf::init(m)?;
    crate::rdmoldescriptors_connectivity::init(m)?;
    crate::rdmoldescriptors_bcut::init(m)?;
    crate::rdmoldescriptors_autocorr::init(m)?;
    crate::rdmoldescriptors_profiles::init(m)?;
    crate::rdmoldescriptors_usr::init(m)?;
    crate::rdmoldescriptors_3d::init(m)?;
    crate::rdmoldescriptors_property::init(m)?;
    m.add_function(wrap_pyfunction!(calc_fraction_csp3, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_heterocycles, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_aromatic_rings, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_saturated_rings, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_aliphatic_rings, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_aromatic_heterocycles, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_aromatic_carbocycles, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_aliphatic_heterocycles, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_aliphatic_carbocycles, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_saturated_heterocycles, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_saturated_carbocycles, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_spiro_atoms_py, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_bridgehead_atoms_py, m)?)?;
    m.add_function(wrap_pyfunction!(calc_num_atom_stereo_centers, m)?)?;
    m.add_function(wrap_pyfunction!(
        calc_num_unspecified_atom_stereo_centers,
        m
    )?)?;
    for (name, version) in LIPINSKI_VERSIONS {
        m.add(*name, *version)?;
    }
    add_enum_exported(
        m,
        "NumRotatableBondsOptions",
        &[
            ("Default", -1),
            ("NonStrict", 0),
            ("Strict", 1),
            ("StrictLinkages", 2),
        ],
        None,
    )?;
    m.add_function(wrap_pyfunction!(calc_crippen_descriptors_py, m)?)?;
    m.add("_CalcCrippenDescriptors_version", CRIPPEN_VERSION)?;
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
