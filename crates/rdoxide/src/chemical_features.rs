//! Ported from RDKit `Code/ChemicalFeatures/FreeChemicalFeature.{h,cpp}`:
//! chemical features that are not tied to a molecule (typically derived from
//! pharmacophores and site maps).

use crate::geometry::Point3D;
use crate::{RDKitError, Result};

/// Version tag of the current pickle format (`ci_FEAT_VERSION`).
const FEAT_VERSION: u32 = 0x0020;
/// Older pickles lack the id field.
const FEAT_VERSION_NO_ID: u32 = 0x0010;

/// `ChemicalFeatures::FreeChemicalFeature`
#[derive(Debug, Clone, PartialEq)]
pub struct FreeChemicalFeature {
    pub id: i32,
    pub family: String,
    pub feature_type: String,
    pub position: Point3D,
}

impl Default for FreeChemicalFeature {
    fn default() -> Self {
        FreeChemicalFeature {
            id: -1,
            family: String::new(),
            feature_type: String::new(),
            position: Point3D::new(0.0, 0.0, 0.0),
        }
    }
}

impl FreeChemicalFeature {
    pub fn new(family: &str, feature_type: &str, position: Point3D, id: i32) -> Self {
        FreeChemicalFeature {
            id,
            family: family.to_string(),
            feature_type: feature_type.to_string(),
            position,
        }
    }

    /// `toString`: the binary pickle (version, id, NUL-terminated family and
    /// type, position; little endian like the platforms RDKit pickles on).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend(FEAT_VERSION.to_le_bytes());
        out.extend(self.id.to_le_bytes());
        for text in [&self.family, &self.feature_type] {
            out.extend((text.len() as u32 + 1).to_le_bytes());
            out.extend(text.as_bytes());
            out.push(0);
        }
        for v in [self.position.x, self.position.y, self.position.z] {
            out.extend(v.to_le_bytes());
        }
        out
    }

    /// `initFromString`
    pub fn from_bytes(pickle: &[u8]) -> Result<Self> {
        let mut r = Reader {
            data: pickle,
            pos: 0,
        };
        let version = r.u32()?;
        if version != FEAT_VERSION && version != FEAT_VERSION_NO_ID {
            return Err(RDKitError::Runtime(
                "Unknown version type for FreeChemicalFeature".into(),
            ));
        }
        let id = if version == FEAT_VERSION {
            r.u32()? as i32
        } else {
            -1
        };
        let family = r.c_string()?;
        let feature_type = r.c_string()?;
        let position = Point3D::new(r.f64()?, r.f64()?, r.f64()?);
        Ok(FreeChemicalFeature {
            id,
            family,
            feature_type,
            position,
        })
    }
}

struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl Reader<'_> {
    fn take(&mut self, n: usize) -> Result<&[u8]> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&e| e <= self.data.len())
            .ok_or_else(|| RDKitError::Runtime("Truncated FreeChemicalFeature pickle".into()))?;
        let bytes = &self.data[self.pos..end];
        self.pos = end;
        Ok(bytes)
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(
            self.take(4)?.try_into().unwrap_or_default(),
        ))
    }
    fn f64(&mut self) -> Result<f64> {
        Ok(f64::from_le_bytes(
            self.take(8)?.try_into().unwrap_or_default(),
        ))
    }
    /// A length-prefixed string; like the C string upstream reads, it ends
    /// at the first NUL.
    fn c_string(&mut self) -> Result<String> {
        let n = self.u32()? as usize;
        let bytes = self.take(n)?;
        let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
        Ok(String::from_utf8_lossy(&bytes[..end]).into_owned())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pickle_round_trip() {
        let f =
            FreeChemicalFeature::new("Donor", "SingleAtomDonor", Point3D::new(1.0, -2.5, 3.25), 7);
        assert_eq!(FreeChemicalFeature::from_bytes(&f.to_bytes()).unwrap(), f);
        let d = FreeChemicalFeature::default();
        assert_eq!(FreeChemicalFeature::from_bytes(&d.to_bytes()).unwrap(), d);
    }

    #[test]
    fn bad_pickles_are_errors() {
        assert!(FreeChemicalFeature::from_bytes(&[]).is_err());
        assert!(FreeChemicalFeature::from_bytes(&[1, 0, 0, 0]).is_err());
        let mut bytes = FreeChemicalFeature::default().to_bytes();
        bytes.truncate(bytes.len() - 3);
        assert!(FreeChemicalFeature::from_bytes(&bytes).is_err());
    }
}
