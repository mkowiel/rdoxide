//! `rdoxide.Chem.rdMIF`: molecular interaction fields (upstream
//! `GraphMol/MolInteractionFields/Wrap/rdMIF.cpp`).

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rdoxide::geometry::Point3D;
use rdoxide::graphmol::mol_interaction_fields::{
    Coulomb, CoulombDielectric, Descriptor, HBond, Hydrophilic, VdWaals, calculate_descriptors,
    construct_grid, read_from_cube_file, write_to_cube_file,
};

use crate::errors::IntoPyResult;
use crate::geometry_grids::PyUniformRealValueGrid3D;
use crate::mol::PyMol;

/// The arguments of a call by name, in the order of `names`; a name may be
/// given by position or keyword.
pub(crate) fn gather<'py>(
    args: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
    names: &[&str],
) -> PyResult<Vec<Option<Bound<'py, PyAny>>>> {
    if args.len() > names.len() {
        return Err(pyo3::exceptions::PyTypeError::new_err("too many arguments"));
    }
    let mut out: Vec<Option<Bound<'py, PyAny>>> = Vec::new();
    for (i, name) in names.iter().enumerate() {
        out.push(if i < args.len() {
            Some(args.get_item(i)?)
        } else if let Some(k) = kwargs {
            k.get_item(name)?
        } else {
            None
        });
    }
    if let Some(k) = kwargs {
        for key in k.keys() {
            let key: String = key.extract()?;
            if !names.contains(&key.as_str()) {
                return Err(pyo3::exceptions::PyTypeError::new_err(format!(
                    "unexpected keyword argument '{key}'"
                )));
            }
        }
    }
    Ok(out)
}

fn get<T: for<'a, 'py> FromPyObject<'a, 'py>>(
    v: &[Option<Bound<'_, PyAny>>],
    i: usize,
    default: T,
) -> PyResult<T>
where
    for<'a, 'py> <T as FromPyObject<'a, 'py>>::Error: Into<PyErr>,
{
    match &v[i] {
        Some(o) => o.extract::<T>().map_err(Into::into),
        None => Ok(default),
    }
}

/// `extractChargesAndPositions`
fn charges_and_positions(
    charges: &Bound<'_, PyAny>,
    positions: &Bound<'_, PyAny>,
) -> PyResult<(Vec<f64>, Vec<Point3D>)> {
    let seq_err = |what: &str| PyValueError::new_err(format!("{what} argument must be a sequence"));
    let ps: Vec<Bound<'_, PyAny>> = positions
        .try_iter()
        .map_err(|_| seq_err("positions"))?
        .collect::<PyResult<_>>()?;
    let cs: Vec<Bound<'_, PyAny>> = charges
        .try_iter()
        .map_err(|_| seq_err("charges"))?
        .collect::<PyResult<_>>()?;
    if ps.len() != cs.len() {
        return Err(PyValueError::new_err(
            "positions and charges must have the same length",
        ));
    }
    let mut pos = Vec::new();
    let mut ch = Vec::new();
    for (p, c) in ps.iter().zip(&cs) {
        let xyz: Vec<f64> = p
            .try_iter()
            .ok()
            .and_then(|it| {
                it.map(|v| v.and_then(|v| v.extract::<f64>()))
                    .collect::<PyResult<_>>()
                    .ok()
            })
            .filter(|v: &Vec<f64>| v.len() == 3)
            .ok_or_else(|| {
                PyValueError::new_err("all elements in positions argument must be x,y,z sequences")
            })?;
        pos.push(Point3D::new(xyz[0], xyz[1], xyz[2]));
        ch.push(c.extract::<f64>()?);
    }
    Ok((ch, pos))
}

/// Electrostatic interaction (Coulomb energy) between probe and molecule in
/// vacuum (no dielectric).
#[pyclass(name = "Coulomb", module = "rdoxide.Chem.rdMIF")]
pub struct PyCoulomb {
    d: Coulomb,
}

#[pymethods]
impl PyCoulomb {
    /// `Coulomb(mol, confId=-1, probeCharge=1.0, absVal=False,
    /// chargeKey='_GasteigerCharge', softcoreParam=0.0, cutoff=1.0)` or
    /// `Coulomb(charges, positions, probeCharge=1.0, absVal=False,
    /// softcoreParam=0.0, cutoff=1.0)`
    #[new]
    #[pyo3(signature = (*args, **kwargs))]
    fn py_new(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let first = args
            .get_item(0)
            .ok()
            .or(kwargs.and_then(|k| k.get_item("mol").ok().flatten()));
        if let Some(m) = first
            .as_ref()
            .and_then(|m| m.extract::<PyRef<'_, PyMol>>().ok())
        {
            let v = gather(
                args,
                kwargs,
                &[
                    "mol",
                    "confId",
                    "probeCharge",
                    "absVal",
                    "chargeKey",
                    "softcoreParam",
                    "cutoff",
                ],
            )?;
            let key: String = get(&v, 4, "_GasteigerCharge".to_string())?;
            return Ok(Self {
                d: Coulomb::from_mol(
                    &m.mol,
                    get(&v, 1, -1)?,
                    get(&v, 2, 1.0)?,
                    get(&v, 3, false)?,
                    &key,
                    get(&v, 5, 0.0)?,
                    get(&v, 6, 1.0)?,
                )
                .py()?,
            });
        }
        let v = gather(
            args,
            kwargs,
            &[
                "charges",
                "positions",
                "probeCharge",
                "absVal",
                "softcoreParam",
                "cutoff",
            ],
        )?;
        let (Some(c), Some(p)) = (&v[0], &v[1]) else {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "charges and positions are required",
            ));
        };
        let (ch, pos) = charges_and_positions(c, p)?;
        Ok(Self {
            d: Coulomb::new(
                &ch,
                &pos,
                get(&v, 2, 1.0)?,
                get(&v, 3, false)?,
                get(&v, 4, 0.0)?,
                get(&v, 5, 1.0)?,
            )
            .py()?,
        })
    }

    /// Calculates the electrostatic interaction between probe and molecule.
    #[pyo3(signature = (x, y, z, threshold))]
    fn __call__(&self, x: f64, y: f64, z: f64, threshold: f64) -> f64 {
        self.d.value(x, y, z, threshold)
    }
}

/// Electrostatic interaction with a distance-dependent dielectric (same
/// energy term as used in GRID MIFs).
#[pyclass(name = "CoulombDielectric", module = "rdoxide.Chem.rdMIF")]
pub struct PyCoulombDielectric {
    d: CoulombDielectric,
}

#[pymethods]
impl PyCoulombDielectric {
    /// `CoulombDielectric(mol, confId=-1, probeCharge=1.0, absVal=False,
    /// chargeKey='_GasteigerCharge', softcoreParam=0.0, cutoff=1.0,
    /// epsilon=80.0, xi=4.0)` or the same with `(charges, positions, ...)`
    #[new]
    #[pyo3(signature = (*args, **kwargs))]
    fn py_new(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let first = args
            .get_item(0)
            .ok()
            .or(kwargs.and_then(|k| k.get_item("mol").ok().flatten()));
        if let Some(m) = first
            .as_ref()
            .and_then(|m| m.extract::<PyRef<'_, PyMol>>().ok())
        {
            let v = gather(
                args,
                kwargs,
                &[
                    "mol",
                    "confId",
                    "probeCharge",
                    "absVal",
                    "chargeKey",
                    "softcoreParam",
                    "cutoff",
                    "epsilon",
                    "xi",
                ],
            )?;
            let key: String = get(&v, 4, "_GasteigerCharge".to_string())?;
            return Ok(Self {
                d: CoulombDielectric::from_mol(
                    &m.mol,
                    get(&v, 1, -1)?,
                    get(&v, 2, 1.0)?,
                    get(&v, 3, false)?,
                    &key,
                    get(&v, 5, 0.0)?,
                    get(&v, 6, 1.0)?,
                    get(&v, 7, 80.0)?,
                    get(&v, 8, 4.0)?,
                )
                .py()?,
            });
        }
        let v = gather(
            args,
            kwargs,
            &[
                "charges",
                "positions",
                "probeCharge",
                "absVal",
                "softcoreParam",
                "cutoff",
                "epsilon",
                "xi",
            ],
        )?;
        let (Some(c), Some(p)) = (&v[0], &v[1]) else {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "charges and positions are required",
            ));
        };
        let (ch, pos) = charges_and_positions(c, p)?;
        Ok(Self {
            d: CoulombDielectric::new(
                &ch,
                &pos,
                get(&v, 2, 1.0)?,
                get(&v, 3, false)?,
                get(&v, 4, 0.0)?,
                get(&v, 5, 1.0)?,
                get(&v, 6, 80.0)?,
                get(&v, 7, 4.0)?,
            )
            .py()?,
        })
    }

    /// Calculates the electrostatic interaction between probe and molecule.
    #[pyo3(signature = (x, y, z, threshold))]
    fn __call__(&self, x: f64, y: f64, z: f64, threshold: f64) -> f64 {
        self.d.value(x, y, z, threshold)
    }
}

/// van der Waals interactions between a molecule and a probe at a
/// gridpoint, based on the MMFF forcefield.
#[pyclass(name = "MMFFVdWaals", module = "rdoxide.Chem.rdMIF")]
pub struct PyMMFFVdWaals {
    d: VdWaals,
}

#[pymethods]
impl PyMMFFVdWaals {
    #[new]
    #[pyo3(signature = (mol, confId = -1, probeAtomType = 6, scaling = false, cutoff = 1.0))]
    #[allow(non_snake_case)]
    fn py_new(
        mol: PyRef<'_, PyMol>,
        confId: i32,
        probeAtomType: u32,
        scaling: bool,
        cutoff: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            d: VdWaals::mmff(&mol.mol, confId, probeAtomType, scaling, cutoff).py()?,
        })
    }

    /// Calculates the van der Waals interaction between molecule and probe.
    #[pyo3(signature = (x, y, z, threshold))]
    fn __call__(&self, x: f64, y: f64, z: f64, threshold: f64) -> f64 {
        self.d.value(x, y, z, threshold)
    }
}

/// van der Waals interactions between a molecule and a probe at a
/// gridpoint, based on the UFF forcefield.
#[pyclass(name = "UFFVdWaals", module = "rdoxide.Chem.rdMIF")]
pub struct PyUFFVdWaals {
    d: VdWaals,
}

#[pymethods]
impl PyUFFVdWaals {
    #[new]
    #[pyo3(signature = (mol, confId = -1, probeAtomType = "O_3", cutoff = 1.0))]
    #[allow(non_snake_case)]
    fn py_new(
        mol: PyRef<'_, PyMol>,
        confId: i32,
        probeAtomType: &str,
        cutoff: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            d: VdWaals::uff(&mol.mol, confId, probeAtomType, cutoff).py()?,
        })
    }

    /// Calculates the van der Waals interaction between molecule and probe.
    #[pyo3(signature = (x, y, z, threshold))]
    fn __call__(&self, x: f64, y: f64, z: f64, threshold: f64) -> f64 {
        self.d.value(x, y, z, threshold)
    }
}

/// Constructs a UniformRealValueGrid3D (3D grid with real values at
/// gridpoints) fitting to a molecule.
#[pyfunction]
#[pyo3(name = "ConstructGrid", signature = (mol, confId = -1, margin = 5.0, spacing = 0.5))]
#[allow(non_snake_case)]
fn py_construct_grid(
    mol: PyRef<'_, PyMol>,
    confId: i32,
    margin: f64,
    spacing: f64,
) -> PyResult<PyUniformRealValueGrid3D> {
    Ok(PyUniformRealValueGrid3D {
        g: construct_grid(&mol.mol, confId, margin, spacing).py()?,
    })
}

/// Hydrogen bonding interaction between a molecule and a probe at a
/// gridpoint (the probe is "OH", "O", "NH" or "N").
#[pyclass(name = "HBond", module = "rdoxide.Chem.rdMIF")]
pub struct PyHBond {
    d: HBond,
}

#[pymethods]
impl PyHBond {
    #[new]
    #[pyo3(signature = (mol, confId = -1, probeAtomType = "OH", fixed = true, cutoff = 1.0))]
    #[allow(non_snake_case)]
    fn py_new(
        mol: PyRef<'_, PyMol>,
        confId: i32,
        probeAtomType: &str,
        fixed: bool,
        cutoff: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            d: HBond::new(&mol.mol, confId, probeAtomType, fixed, cutoff).py()?,
        })
    }

    #[pyo3(name = "__call__")]
    fn call(&self, x: f64, y: f64, z: f64, threshold: f64) -> f64 {
        self.d.value(x, y, z, threshold)
    }
}

/// Hydrophilic field: the minimum of the hydrogen bond interactions with a
/// water oxygen and a water hydrogen probe.
#[pyclass(name = "Hydrophilic", module = "rdoxide.Chem.rdMIF")]
pub struct PyHydrophilic {
    d: Hydrophilic,
}

#[pymethods]
impl PyHydrophilic {
    #[new]
    #[pyo3(signature = (mol, confId = -1, fixed = true, cutoff = 1.0))]
    #[allow(non_snake_case)]
    fn py_new(mol: PyRef<'_, PyMol>, confId: i32, fixed: bool, cutoff: f64) -> PyResult<Self> {
        Ok(Self {
            d: Hydrophilic::new(&mol.mol, confId, fixed, cutoff).py()?,
        })
    }

    #[pyo3(name = "__call__")]
    fn call(&self, x: f64, y: f64, z: f64, threshold: f64) -> f64 {
        self.d.value(x, y, z, threshold)
    }
}

fn run<D: Descriptor>(
    py: Python<'_>,
    grid: &Bound<'_, PyUniformRealValueGrid3D>,
    d: &D,
    threshold: f64,
) -> PyResult<()> {
    let mut grid = grid.borrow_mut();
    let g = &mut grid.g;
    py.detach(|| calculate_descriptors(g, d, threshold, 0)).py()
}

/// Calculates descriptors (to be specified as parameter) of a molecule at
/// every gridpoint of a grid.
#[pyfunction]
#[pyo3(name = "CalculateDescriptors", signature = (grid, descriptor, threshold = -1.0))]
fn py_calculate_descriptors(
    py: Python<'_>,
    grid: &Bound<'_, PyUniformRealValueGrid3D>,
    descriptor: &Bound<'_, PyAny>,
    threshold: f64,
) -> PyResult<()> {
    if let Ok(d) = descriptor.extract::<PyRef<'_, PyCoulomb>>() {
        run(py, grid, &d.d, threshold)
    } else if let Ok(d) = descriptor.extract::<PyRef<'_, PyCoulombDielectric>>() {
        run(py, grid, &d.d, threshold)
    } else if let Ok(d) = descriptor.extract::<PyRef<'_, PyMMFFVdWaals>>() {
        run(py, grid, &d.d, threshold)
    } else if let Ok(d) = descriptor.extract::<PyRef<'_, PyUFFVdWaals>>() {
        run(py, grid, &d.d, threshold)
    } else if let Ok(d) = descriptor.extract::<PyRef<'_, PyHBond>>() {
        run(py, grid, &d.d, threshold)
    } else if let Ok(d) = descriptor.extract::<PyRef<'_, PyHydrophilic>>() {
        run(py, grid, &d.d, threshold)
    } else {
        Err(pyo3::exceptions::PyTypeError::new_err(
            "descriptor must be one of the descriptor classes of rdMIF",
        ))
    }
}

/// Writes Grid to a file in Gaussian CUBE format.
#[pyfunction]
#[pyo3(name = "WriteToCubeFile", signature = (grid, filename, mol = None, confId = -1))]
#[allow(non_snake_case)]
fn py_write_to_cube_file(
    grid: PyRef<'_, PyUniformRealValueGrid3D>,
    filename: &str,
    mol: Option<PyRef<'_, PyMol>>,
    confId: i32,
) -> PyResult<()> {
    write_to_cube_file(&grid.g, filename, mol.as_ref().map(|m| &m.mol), confId).py()
}

/// Reads Grid from a file in Gaussian CUBE format.
///
/// RETURNS: a tuple where the first element is the grid and the second
/// element is the molecule object associated to the grid (only atoms and
/// coordinates, no bonds; None if no molecule was associated to the grid)
#[pyfunction]
#[pyo3(name = "ReadFromCubeFile", signature = (filename))]
fn py_read_from_cube_file(filename: &str) -> PyResult<(PyUniformRealValueGrid3D, Option<PyMol>)> {
    let (g, mol) = read_from_cube_file(filename).py()?;
    Ok((PyUniformRealValueGrid3D { g }, mol.map(|mol| PyMol { mol })))
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Module containing functions for calculating molecular interaction fields (MIFs)\n  NOTE: This functionality is experimental and the API and/or results may change in future releases.",
    )?;
    m.add_class::<PyCoulomb>()?;
    m.add_class::<PyCoulombDielectric>()?;
    m.add_class::<PyMMFFVdWaals>()?;
    m.add_class::<PyUFFVdWaals>()?;
    m.add_class::<PyHBond>()?;
    m.add_class::<PyHydrophilic>()?;
    m.add_function(wrap_pyfunction!(py_construct_grid, m)?)?;
    m.add_function(wrap_pyfunction!(py_calculate_descriptors, m)?)?;
    m.add_function(wrap_pyfunction!(py_write_to_cube_file, m)?)?;
    m.add_function(wrap_pyfunction!(py_read_from_cube_file, m)?)?;
    Ok(())
}
