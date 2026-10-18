//! `rdMolDescriptors.DoubleCubicLatticeVolume` (upstream `Descriptors/DCLV.cpp`
//! and its wrapper in `Wrap/rdMolDescriptors.cpp`).

use std::collections::BTreeMap;

use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use rdoxide::graphmol::descriptors::dclv::DoubleCubicLatticeVolume;

use crate::errors::IntoPyResult;
use crate::geometry::PyPoint3D;
use crate::mol::PyMol;

/// Class for the Double Cubic Lattice Volume method
///
/// Calculation of the Shrake and Rupley surface area and volume using the
/// Double Cubic Lattice Method (Eisenhaber et al., J. Comput. Chem. 16, 273 (1995)).
///
///   ARGUMENTS:
///    - mol: input molecule or protein
///    - radii: (optional) radii for the atoms of mol, default: van der Waals radii
///    - isProtein: flag to indicate if the input is a protein (default=False, free ligand).
///    - includeLigand: flag to include or exclude a bound ligand when input is a protein (default=True)
///    - probeRadius: radius of the solvent probe (default=1.4)
///    - confId: conformer ID to consider (default=-1)
#[pyclass(
    name = "DoubleCubicLatticeVolume",
    module = "rdoxide.Chem.rdMolDescriptors"
)]
pub struct PyDoubleCubicLatticeVolume {
    inner: DoubleCubicLatticeVolume,
}

/// The arguments of `DoubleCubicLatticeVolume(mol, [radii,] isProtein=False,
/// includeLigand=True, probeRadius=1.4, confId=-1)`, positional or by name.
struct Arguments {
    radii: Vec<f64>,
    is_protein: bool,
    include_ligand: bool,
    probe_radius: f64,
    conf_id: i32,
}

impl Arguments {
    fn parse(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let mut rest: Vec<Bound<'_, PyAny>> = args.iter().skip(1).collect();
        let mut radii = Vec::new();
        // the radii come first when given by position
        if let Some(first) = rest.first() {
            if first.is_instance_of::<PyList>() || first.is_instance_of::<PyTuple>() {
                radii = first.extract()?;
                rest.remove(0);
            }
        }
        let get = |position: usize, name: &str| -> PyResult<Option<Bound<'_, PyAny>>> {
            if let Some(v) = rest.get(position) {
                return Ok(Some(v.clone()));
            }
            kwargs.map_or(Ok(None), |k| k.get_item(name))
        };
        if let Some(r) = kwargs.map_or(Ok(None), |k| k.get_item("radii"))? {
            radii = r.extract()?;
        }
        Ok(Arguments {
            radii,
            is_protein: get(0, "isProtein")?.map_or(Ok(false), |v| v.is_truthy())?,
            include_ligand: get(1, "includeLigand")?.map_or(Ok(true), |v| v.is_truthy())?,
            probe_radius: get(2, "probeRadius")?.map_or(Ok(1.4), |v| v.extract())?,
            conf_id: get(3, "confId")?.map_or(Ok(-1), |v| v.extract())?,
        })
    }
}

/// The atoms named by a sequence of indices as a flag per atom.
fn atom_flags(indices: &Bound<'_, PyAny>, num_atoms: u32) -> PyResult<Vec<bool>> {
    let indices: Vec<u32> = indices.extract()?;
    if indices.is_empty() {
        return Err(PyValueError::new_err(
            "No atom indices supplied for Partial Surface Area",
        ));
    }
    let mut flags = vec![false; num_atoms as usize];
    for i in indices {
        *flags
            .get_mut(i as usize)
            .ok_or_else(|| PyIndexError::new_err("atom index out of range"))? = true;
    }
    Ok(flags)
}

#[pymethods]
impl PyDoubleCubicLatticeVolume {
    #[new]
    #[pyo3(signature = (*args, **kwargs))]
    fn new(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let mol = match args.iter().next() {
            Some(m) => m.extract::<PyRef<'_, PyMol>>()?,
            None => match kwargs.map_or(Ok(None), |k| k.get_item("mol"))? {
                Some(m) => m.extract::<PyRef<'_, PyMol>>()?,
                None => return Err(PyTypeError::new_err("missing the molecule argument")),
            },
        };
        let a = Arguments::parse(args, kwargs)?;
        let inner = DoubleCubicLatticeVolume::new(
            &mol.mol,
            a.radii,
            a.is_protein,
            a.include_ligand,
            a.probe_radius,
            a.conf_id,
        )
        .py()?;
        Ok(PyDoubleCubicLatticeVolume { inner })
    }

    /// Get the Surface Area of the Molecule or Protein
    #[pyo3(name = "GetSurfaceArea")]
    fn get_surface_area(&mut self, py: Python<'_>) -> PyResult<f64> {
        let inner = &mut self.inner;
        py.detach(|| inner.surface_area()).py()
    }

    /// Get the surface area of atom with atom_idx
    #[pyo3(name = "GetAtomSurfaceArea", signature = (atom_idx))]
    fn get_atom_surface_area(&mut self, py: Python<'_>, atom_idx: u32) -> PyResult<f64> {
        let inner = &mut self.inner;
        py.detach(|| inner.atom_surface_area(atom_idx)).py()
    }

    /// Get the Polar Surface Area of the Molecule or Protein
    #[pyo3(name = "GetPolarSurfaceArea", signature = (includeSandP = false, includeHs = false))]
    #[allow(non_snake_case)]
    fn get_polar_surface_area(
        &mut self,
        py: Python<'_>,
        includeSandP: bool,
        includeHs: bool,
    ) -> PyResult<f64> {
        let inner = &mut self.inner;
        py.detach(|| inner.polar_surface_area(includeSandP, includeHs))
            .py()
    }

    /// Get the Partial Surface Area of the Molecule or Protein for specified subset of atoms
    #[pyo3(name = "GetPartialSurfaceArea", signature = (atomIndices))]
    #[allow(non_snake_case)]
    fn get_partial_surface_area(
        &mut self,
        py: Python<'_>,
        atomIndices: &Bound<'_, PyAny>,
    ) -> PyResult<f64> {
        let flags = atom_flags(atomIndices, self.inner.num_atoms())?;
        let inner = &mut self.inner;
        py.detach(|| inner.partial_surface_area(&flags)).py()
    }

    /// Get the set of points representing the surface. If allPoints is True, returns all surface points; otherwise, returns the standard surface points.
    #[pyo3(name = "GetSurfacePoints", signature = (allPoints = false))]
    #[allow(non_snake_case)]
    fn get_surface_points(
        &mut self,
        py: Python<'_>,
        allPoints: bool,
    ) -> PyResult<BTreeMap<u32, Vec<PyPoint3D>>> {
        let inner = &mut self.inner;
        let points = py
            .detach(|| inner.surface_points(allPoints).cloned())
            .py()?;
        Ok(points
            .into_iter()
            .map(|(i, v)| (i, v.into_iter().map(|p| PyPoint3D { p }).collect()))
            .collect())
    }

    /// Get the Total Volume of the Molecule or Protein
    #[pyo3(name = "GetVolume")]
    fn get_volume(&mut self, py: Python<'_>) -> PyResult<f64> {
        let inner = &mut self.inner;
        py.detach(|| inner.volume()).py()
    }

    /// Get the van der Waals Volume of the Molecule or Protein
    #[pyo3(name = "GetVDWVolume")]
    fn get_vdw_volume(&mut self, py: Python<'_>) -> PyResult<f64> {
        let inner = &mut self.inner;
        py.detach(|| inner.vdw_volume()).py()
    }

    /// Get the volume atom of atom_idx with volume for specified Probe Radius
    #[pyo3(name = "GetAtomVolume", signature = (atomIdx, solventRadius))]
    #[allow(non_snake_case)]
    fn get_atom_volume(
        &mut self,
        py: Python<'_>,
        atomIdx: u32,
        solventRadius: f64,
    ) -> PyResult<f64> {
        let inner = &mut self.inner;
        py.detach(|| inner.atom_volume(atomIdx, solventRadius)).py()
    }

    /// Get the Polar Volume of the Molecule or Protein
    #[pyo3(name = "GetPolarVolume", signature = (includeSandP = false, includeHs = false))]
    #[allow(non_snake_case)]
    fn get_polar_volume(
        &mut self,
        py: Python<'_>,
        includeSandP: bool,
        includeHs: bool,
    ) -> PyResult<f64> {
        let inner = &mut self.inner;
        py.detach(|| inner.polar_volume(includeSandP, includeHs))
            .py()
    }

    /// Get the Partial Volume of the Molecule or Protein for specified subset of atoms
    #[pyo3(name = "GetPartialVolume", signature = (atomIdx))]
    #[allow(non_snake_case)]
    fn get_partial_volume(&mut self, py: Python<'_>, atomIdx: &Bound<'_, PyAny>) -> PyResult<f64> {
        let flags = atom_flags(atomIdx, self.inner.num_atoms())?;
        let inner = &mut self.inner;
        py.detach(|| inner.partial_volume(&flags)).py()
    }

    /// Get the Compactness of the Protein
    #[pyo3(name = "GetCompactness")]
    fn get_compactness(&mut self, py: Python<'_>) -> PyResult<f64> {
        let inner = &mut self.inner;
        py.detach(|| inner.compactness()).py()
    }

    /// Get the PackingDensity of the Protein
    #[pyo3(name = "GetPackingDensity")]
    fn get_packing_density(&mut self, py: Python<'_>) -> PyResult<f64> {
        let inner = &mut self.inner;
        py.detach(|| inner.packing_density()).py()
    }
}

pub fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyDoubleCubicLatticeVolume>()
}
