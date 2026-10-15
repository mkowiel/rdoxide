//! Torsion candidates and the Z-matrix used to place atoms from internal
//! coordinates (upstream `DistGeom/ZMatrix{,Utils}.h` and
//! `computeZMatrixCoords` of `DistGeomUtils.cpp`).
//!
//! Upstream stores fixed torsion values in an `unordered_set`, so the order
//! in which values are drawn depends on the standard library. Here the values
//! are kept sorted, which makes the draws reproducible.

use std::f64::consts::PI;

use crate::error::{RDKitError, Result};
use crate::geometry::Point3D;
use crate::rdgeneral::random::MinStdRand;

const TORSION_PRECISION: f64 = 1.0e9;

fn dequantize(val: i64) -> f64 {
    val as f64 / TORSION_PRECISION
}

/// Maps a torsion to [-pi, pi) and truncates it to a fixed-point integer.
pub fn quantize(val: f64) -> i64 {
    let mut val = (val + PI) % (2.0 * PI);
    if val < 0.0 {
        val += 2.0 * PI;
    }
    val -= PI;
    (val * TORSION_PRECISION) as i64
}

/// A range of torsions that may wrap around +-pi.
#[derive(Debug, Clone, Copy)]
pub struct TorsionRange {
    pub lower: f64,
    pub upper: f64,
    q_lower: i64,
    q_upper: i64,
}

impl TorsionRange {
    pub fn new(lower: f64, upper: f64) -> Result<Self> {
        if lower > upper {
            return Err(RDKitError::precondition("Invalid range"));
        }
        let (mut q_lower, mut q_upper) = (quantize(lower), quantize(upper));
        if q_lower == q_upper {
            // the full range
            q_upper = (PI * TORSION_PRECISION) as i64;
            q_lower = -((PI * TORSION_PRECISION) as i64);
        }
        Ok(TorsionRange {
            lower,
            upper,
            q_lower,
            q_upper,
        })
    }

    pub fn from_quantized(q_lower: i64, q_upper: i64) -> Self {
        TorsionRange {
            lower: dequantize(q_lower),
            upper: dequantize(q_upper),
            q_lower,
            q_upper,
        }
    }

    pub fn q_contains(&self, value: i64) -> bool {
        if self.q_lower <= self.q_upper {
            value >= self.q_lower && value <= self.q_upper
        } else {
            // a range across pi and -pi
            (value < self.q_lower) != (value > self.q_upper)
        }
    }

    pub fn contains(&self, value: f64) -> bool {
        self.q_contains(quantize(value))
    }

    pub fn width(&self) -> f64 {
        let res = self.upper - self.lower;
        if res < 0.0 { res + 2.0 * PI } else { res }
    }

    pub fn sample(&self, rng: &mut MinStdRand) -> f64 {
        self.lower + (self.upper - self.lower) * rng.next_f64()
    }
}

impl PartialEq for TorsionRange {
    fn eq(&self, other: &Self) -> bool {
        self.q_lower == other.q_lower && self.q_upper == other.q_upper
    }
}

/// A set of distinct torsion values (quantized, sorted).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TorsionValues {
    content: Vec<i64>,
}

impl TorsionValues {
    pub fn new(values: &[f64]) -> Self {
        let mut v = TorsionValues::default();
        for &x in values {
            v.insert(x);
        }
        v
    }

    pub fn is_empty(&self) -> bool {
        self.content.is_empty()
    }

    pub fn len(&self) -> usize {
        self.content.len()
    }

    pub fn iter(&self) -> impl Iterator<Item = i64> + '_ {
        self.content.iter().copied()
    }

    pub fn insert(&mut self, val: f64) -> bool {
        self.insert_q(quantize(val))
    }

    pub fn insert_q(&mut self, val: i64) -> bool {
        match self.content.binary_search(&val) {
            Ok(_) => false,
            Err(pos) => {
                self.content.insert(pos, val);
                true
            }
        }
    }

    /// The values the two sets have in common (upstream's "merge" keeps the
    /// intersection in the left-hand set).
    pub fn merge(&self, other: &TorsionValues) -> TorsionValues {
        TorsionValues {
            content: self
                .content
                .iter()
                .copied()
                .filter(|v| other.content.binary_search(v).is_ok())
                .collect(),
        }
    }

    pub fn sample(&self, rng: &mut MinStdRand) -> f64 {
        let idx = (rng.next_f64() * self.len() as f64) as usize;
        dequantize(self.content[idx])
    }

    pub fn first_value(&self) -> f64 {
        dequantize(self.content[0])
    }

    pub fn contains(&self, value: f64) -> bool {
        self.content.binary_search(&quantize(value)).is_ok()
    }
}

/// A range or a set of values.
#[derive(Debug, Clone, PartialEq)]
pub enum TorsionCandidates {
    Range(TorsionRange),
    Values(TorsionValues),
}

impl From<TorsionRange> for TorsionCandidates {
    fn from(r: TorsionRange) -> Self {
        TorsionCandidates::Range(r)
    }
}

impl From<TorsionValues> for TorsionCandidates {
    fn from(v: TorsionValues) -> Self {
        TorsionCandidates::Values(v)
    }
}

fn merge_ranges(lhs: &TorsionRange, rhs: &TorsionRange) -> Result<TorsionCandidates> {
    if lhs.upper.min(rhs.upper) >= lhs.lower.max(rhs.lower) {
        Ok(TorsionRange::new(lhs.lower.max(rhs.lower), lhs.upper.min(rhs.upper))?.into())
    } else {
        Ok(TorsionRange::new(lhs.lower.min(rhs.lower), lhs.upper.max(rhs.upper))?.into())
    }
}

fn merge_range_values(range: &TorsionRange, values: &TorsionValues) -> Result<TorsionCandidates> {
    // the fixed values within the range; without any, the span of the values
    let mut merged = TorsionValues::default();
    for v in values.iter() {
        if range.q_contains(v) {
            merged.insert_q(v);
        }
    }
    if !merged.is_empty() {
        return Ok(merged.into());
    }
    let min = values
        .iter()
        .min()
        .ok_or_else(|| RDKitError::precondition("empty torsion values"))?;
    let max = values.iter().max().unwrap_or(min);
    merge_ranges(range, &TorsionRange::from_quantized(min, max))
}

impl TorsionCandidates {
    /// Combines two sets of candidates into one that satisfies both as far as
    /// possible.
    pub fn merge(&self, other: &TorsionCandidates) -> Result<TorsionCandidates> {
        use TorsionCandidates::*;
        match (self, other) {
            (Range(a), Range(b)) => merge_ranges(a, b),
            (Range(r), Values(v)) | (Values(v), Range(r)) => merge_range_values(r, v),
            (Values(a), Values(b)) => Ok(a.merge(b).into()),
        }
    }

    pub fn sample(&self, rng: &mut MinStdRand) -> f64 {
        match self {
            TorsionCandidates::Range(r) => r.sample(rng),
            TorsionCandidates::Values(v) => v.sample(rng),
        }
    }

    /// Whether `self` is more constrained than `other`.
    pub fn less(&self, other: &TorsionCandidates) -> bool {
        use TorsionCandidates::*;
        match (self, other) {
            (Range(a), Range(b)) => a.width() < b.width(),
            (Range(_), Values(_)) => false,
            (Values(_), Range(_)) => true,
            (Values(a), Values(b)) => a.len() < b.len(),
        }
    }

    pub fn contains(&self, value: f64) -> bool {
        match self {
            TorsionCandidates::Range(r) => r.contains(value),
            TorsionCandidates::Values(v) => v.contains(value),
        }
    }
}

/// The internal coordinates that place one atom.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct InternalCoordinateReference {
    pub length: Option<f64>,
    pub bond_ref: Option<usize>,
    pub angle: Option<f64>,
    pub angle_ref: Option<usize>,
    pub torsion: Option<TorsionCandidates>,
    pub torsion_ref: Option<usize>,
}

/// A torsion that is the offset of another one.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TorsionDependence {
    pub offset: f64,
    pub reference: usize,
}

/// A Z-matrix whose torsions are ranges or sets of values, and which also
/// records torsions that depend on others (direct neighbors that are not the
/// bond reference and come earlier in the matrix).
#[derive(Debug, Clone)]
pub struct ZMatrix {
    order: Vec<usize>,
    refs: Vec<Option<InternalCoordinateReference>>,
    dependences: Vec<Option<TorsionDependence>>,
}

/// One row of a [`ZMatrix`].
#[derive(Debug, Clone, Copy)]
pub struct ZMatrixRow<'a> {
    pub atom_idx: usize,
    pub internal: &'a InternalCoordinateReference,
    pub torsion_dependence: Option<&'a TorsionDependence>,
}

impl ZMatrix {
    pub fn new(n_elements: usize) -> Self {
        ZMatrix {
            order: Vec::with_capacity(n_elements),
            refs: vec![None; n_elements],
            dependences: vec![None; n_elements],
        }
    }

    pub fn len(&self) -> usize {
        self.order.len()
    }

    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }

    pub fn add_element(
        &mut self,
        idx: usize,
        internal: InternalCoordinateReference,
        torsion_dep: Option<TorsionDependence>,
    ) -> Result<()> {
        if idx >= self.refs.len() {
            return Err(RDKitError::precondition("Index out of range"));
        }
        if self.refs[idx].is_some() {
            return Err(RDKitError::precondition("element already set"));
        }
        self.order.push(idx);
        self.refs[idx] = Some(internal);
        self.dependences[idx] = torsion_dep;
        Ok(())
    }

    pub fn references(&self, atom_idx: usize) -> Result<&InternalCoordinateReference> {
        self.refs
            .get(atom_idx)
            .ok_or_else(|| RDKitError::precondition("Invalid atom index"))?
            .as_ref()
            .ok_or_else(|| {
                RDKitError::precondition("Accessing an ZMatrix Element that was not set before")
            })
    }

    pub fn torsion_reference(&self, atom_idx: usize) -> Result<Option<&TorsionDependence>> {
        self.dependences
            .get(atom_idx)
            .map(|d| d.as_ref())
            .ok_or_else(|| RDKitError::precondition("Invalid atom index"))
    }

    pub fn update_torsion(&mut self, atom_idx: usize, torsion: TorsionRange) -> Result<()> {
        match self.refs.get_mut(atom_idx) {
            Some(Some(r)) => {
                r.torsion = Some(torsion.into());
                Ok(())
            }
            _ => Err(RDKitError::precondition("Invalid atom index")),
        }
    }

    pub fn invert_torsion_dependence(&mut self, atom_idx: usize) -> Result<()> {
        match self.dependences.get_mut(atom_idx) {
            Some(Some(d)) => {
                d.offset = -d.offset;
                Ok(())
            }
            Some(None) => Err(RDKitError::precondition("Atom has no torsion dependence")),
            None => Err(RDKitError::precondition("Invalid atom index")),
        }
    }

    /// The row at position `idx` of the order in which atoms are placed.
    pub fn row(&self, idx: usize) -> Result<ZMatrixRow<'_>> {
        let atom_idx = *self
            .order
            .get(idx)
            .ok_or_else(|| RDKitError::precondition("Invalid row index"))?;
        Ok(ZMatrixRow {
            atom_idx,
            internal: self.references(atom_idx)?,
            torsion_dependence: self.dependences[atom_idx].as_ref(),
        })
    }

    pub fn rows(&self) -> impl DoubleEndedIterator<Item = ZMatrixRow<'_>> {
        (0..self.order.len()).map(|i| self.row(i).expect("rows are set"))
    }
}

fn position_from_references(
    ref1: &Point3D,
    ref2: &Point3D,
    ref3: &Point3D,
    length: f64,
    angle: f64,
    torsion: f64,
) -> Result<Point3D> {
    let mut v_1_2 = *ref1 - *ref2;
    let mut v_2_3 = *ref2 - *ref3;
    v_1_2.normalize()?;
    v_2_3.normalize()?;
    let c_prod = v_2_3.cross_product(&v_1_2);
    let dot_p = v_1_2.dot_product(&v_2_3);
    let denom = (1.0 - dot_p.powi(2)).sqrt();
    let n_cp = c_prod / denom;
    let cp2 = n_cp.cross_product(&v_1_2);
    let new_pos = (-v_1_2 * angle.cos()
        + cp2 * (angle.sin() * torsion.cos())
        + n_cp * (angle.sin() * torsion.sin()))
        * length;
    Ok(*ref1 + new_pos)
}

/// Places the atoms of a Z-matrix, drawing torsions from their candidates.
/// `coordinates` has one point per atom.
pub fn compute_z_matrix_coords(
    zmat: &ZMatrix,
    coordinates: &mut [Point3D],
    rng: &mut MinStdRand,
) -> Result<()> {
    let num_atoms = coordinates.len();
    let missing = |what: &str| RDKitError::invariant(format!("Z-matrix row without {what}"));
    if num_atoms == 0 {
        return Ok(());
    }
    // the first atom at the origin
    coordinates[zmat.row(0)?.atom_idx] = Point3D::new(0.0, 0.0, 0.0);
    if num_atoms == 1 {
        return Ok(());
    }
    // the second atom on the z axis
    let r1 = zmat.row(1)?;
    coordinates[r1.atom_idx] = Point3D::new(
        0.0,
        0.0,
        r1.internal.length.ok_or_else(|| missing("a length"))?,
    );
    if num_atoms == 2 {
        return Ok(());
    }
    // the third atom in the xz plane
    let r2 = zmat.row(2)?;
    {
        let bl = r2.internal.length.ok_or_else(|| missing("a length"))?;
        let ba = r2.internal.angle.ok_or_else(|| missing("an angle"))?;
        let bond_ref = r2
            .internal
            .bond_ref
            .ok_or_else(|| missing("a bond reference"))?;
        let z_pos = coordinates[bond_ref].z - bl * ba.cos();
        coordinates[r2.atom_idx] = Point3D::new(bl * ba.sin(), 0.0, z_pos);
    }
    for row in zmat.rows().skip(3) {
        let internals = row.internal;
        let torsion_ref = match row.torsion_dependence {
            Some(dep) => dep.reference,
            None => internals
                .torsion_ref
                .ok_or_else(|| missing("a torsion reference"))?,
        };
        let ref1 = coordinates[internals
            .bond_ref
            .ok_or_else(|| missing("a bond reference"))?];
        let ref2 = coordinates[internals
            .angle_ref
            .ok_or_else(|| missing("an angle reference"))?];
        let ref3 = coordinates[torsion_ref];
        let torsion = match row.torsion_dependence {
            Some(dep) => dep.offset,
            None => internals
                .torsion
                .as_ref()
                .ok_or_else(|| missing("a torsion"))?
                .sample(rng),
        };
        coordinates[row.atom_idx] = position_from_references(
            &ref1,
            &ref2,
            &ref3,
            internals.length.ok_or_else(|| missing("a length"))?,
            internals.angle.ok_or_else(|| missing("an angle"))?,
            torsion,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantization_wraps_around() {
        assert_eq!(quantize(0.0), 0);
        assert_eq!(quantize(2.0 * PI), quantize(0.0));
        assert!((dequantize(quantize(PI - 0.1)) - (PI - 0.1)).abs() < 1e-8);
        assert!((dequantize(quantize(PI + 0.1)) - (-PI + 0.1)).abs() < 1e-8);
    }

    #[test]
    fn ranges_across_pi() {
        let r = TorsionRange::new(PI - 0.2, PI + 0.2).unwrap();
        assert!(r.contains(PI));
        assert!(r.contains(-PI + 0.1));
        assert!(!r.contains(0.0));
        assert!((r.width() - 0.4).abs() < 1e-12 || r.width() > 0.0);
        // a degenerate range means the full circle
        let full = TorsionRange::new(1.0, 1.0).unwrap();
        assert!(full.contains(-3.0) && full.contains(3.0));
        assert!(TorsionRange::new(1.0, 0.0).is_err());
    }

    #[test]
    fn merging_candidates() {
        let a: TorsionCandidates = TorsionValues::new(&[0.0, PI]).into();
        let b: TorsionCandidates = TorsionValues::new(&[PI]).into();
        let m = a.merge(&b).unwrap();
        assert_eq!(m, TorsionCandidates::Values(TorsionValues::new(&[PI])));
        // fixed values inside a range win
        let r: TorsionCandidates = TorsionRange::new(-0.5, 0.5).unwrap().into();
        let m = r.merge(&a).unwrap();
        assert_eq!(m, TorsionCandidates::Values(TorsionValues::new(&[0.0])));
        // two ranges: the intersection, or the hull when they are disjoint
        let r1 = TorsionRange::new(0.0, 1.0).unwrap();
        let r2 = TorsionRange::new(0.5, 2.0).unwrap();
        match TorsionCandidates::from(r1).merge(&r2.into()).unwrap() {
            TorsionCandidates::Range(m) => {
                assert!((m.lower - 0.5).abs() < 1e-12 && (m.upper - 1.0).abs() < 1e-12)
            }
            _ => panic!("expected a range"),
        }
        assert!(a.less(&r) && !r.less(&a));
    }

    #[test]
    fn places_a_chain() {
        let mut z = ZMatrix::new(4);
        let plain = |length,
                     bond_ref,
                     angle,
                     angle_ref,
                     torsion: Option<TorsionCandidates>,
                     torsion_ref| {
            InternalCoordinateReference {
                length,
                bond_ref,
                angle,
                angle_ref,
                torsion,
                torsion_ref,
            }
        };
        z.add_element(0, InternalCoordinateReference::default(), None)
            .unwrap();
        z.add_element(1, plain(Some(1.5), Some(0), None, None, None, None), None)
            .unwrap();
        z.add_element(
            2,
            plain(Some(1.5), Some(1), Some(1.9), Some(0), None, None),
            None,
        )
        .unwrap();
        z.add_element(
            3,
            plain(
                Some(1.5),
                Some(2),
                Some(1.9),
                Some(1),
                Some(TorsionValues::new(&[PI]).into()),
                Some(0),
            ),
            None,
        )
        .unwrap();
        let mut pts = vec![Point3D::new(0.0, 0.0, 0.0); 4];
        let mut rng = MinStdRand::new(1);
        compute_z_matrix_coords(&z, &mut pts, &mut rng).unwrap();
        let d = |a: usize, b: usize| (pts[a] - pts[b]).length();
        assert!((d(0, 1) - 1.5).abs() < 1e-9);
        assert!((d(1, 2) - 1.5).abs() < 1e-9);
        assert!((d(2, 3) - 1.5).abs() < 1e-9);
        // the trans torsion puts atoms 0 and 3 far apart
        assert!(d(0, 3) > 3.5);
        assert!(
            z.add_element(3, InternalCoordinateReference::default(), None)
                .is_err()
        );
    }
}
