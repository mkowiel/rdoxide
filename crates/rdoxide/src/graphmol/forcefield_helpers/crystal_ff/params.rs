//! The collections of experimental torsion-angle preferences (upstream
//! `ExpTorsionAngleCollection` of `TorsionPreferences.cpp`): SMARTS patterns
//! with cosine-series constants (ETKDG versions 1 and 2) or Gaussian fits
//! (version 4).

use std::sync::OnceLock;

use crate::error::{RDKitError, Result};
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::rdgeneral::common_properties as cp;

/// A torsion preference as a sum of six cosines.
#[derive(Debug)]
pub struct ExpTorsionAngle {
    pub torsion_idx: usize,
    pub smarts: String,
    pub v: Vec<f64>,
    pub signs: Vec<i32>,
    pub pattern: ROMol,
    /// the pattern atoms that carry the atom map numbers 1 to 4
    pub idx: [u32; 4],
}

/// A torsion preference as a sum of Gaussians.
#[derive(Debug)]
pub struct GaussianExpTorsionAngle {
    pub torsion_idx: usize,
    pub smarts: String,
    pub heights: Vec<f64>,
    pub positions: Vec<f64>,
    pub widths: Vec<f64>,
    pub pattern: ROMol,
    pub idx: [u32; 4],
}

/// Parses the SMARTS and finds the atoms that carry the map numbers 1 to 4.
fn parse_pattern(smarts: &str) -> Result<(ROMol, [u32; 4])> {
    let pattern = mol_from_smarts(smarts, &SmartsParserParams::default())?
        .ok_or_else(|| RDKitError::Value(format!("bad torsion SMARTS: {smarts}")))?;
    let mut idx = [0u32; 4];
    for i in 0..pattern.num_atoms() {
        let num = pattern
            .atom(i)
            .props
            .get(cp::molAtomMapNumber)
            .and_then(|v| v.as_int());
        if let Some(num) = num {
            if (1..5).contains(&num) {
                idx[num as usize - 1] = i;
            }
        }
    }
    Ok((pattern, idx))
}

fn parse_cosine(data: &[&str]) -> Result<Vec<ExpTorsionAngle>> {
    let mut params = Vec::new();
    for line in data.iter().flat_map(|d| d.lines()) {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut tokens = line.split(' ').filter(|t| !t.is_empty());
        let smarts = tokens.next().expect("a SMARTS").to_string();
        let (mut signs, mut v) = (Vec::new(), Vec::new());
        for _ in 0..6 {
            let bad = || RDKitError::Value(format!("bad torsion parameters: {line}"));
            signs.push(tokens.next().ok_or_else(bad)?.parse().map_err(|_| bad())?);
            v.push(tokens.next().ok_or_else(bad)?.parse().map_err(|_| bad())?);
        }
        let (pattern, idx) = parse_pattern(&smarts)?;
        params.push(ExpTorsionAngle {
            torsion_idx: params.len(),
            smarts,
            v,
            signs,
            pattern,
            idx,
        });
    }
    Ok(params)
}

fn parse_gaussian(data: &[&str]) -> Result<Vec<GaussianExpTorsionAngle>> {
    let mut params = Vec::new();
    for line in data.iter().flat_map(|d| d.lines()) {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        let bad = || RDKitError::Value(format!("bad Gaussian torsion parameters: {line}"));
        if cols.len() != 4 {
            return Err(bad());
        }
        let numbers = |col: &str| -> Result<Vec<f64>> {
            col.split(',')
                .map(|v| v.parse().map_err(|_| bad()))
                .collect()
        };
        let (pattern, idx) = parse_pattern(cols[0])?;
        params.push(GaussianExpTorsionAngle {
            torsion_idx: params.len(),
            smarts: cols[0].to_string(),
            heights: numbers(cols[1])?,
            positions: numbers(cols[2])?,
            widths: numbers(cols[3])?,
            pattern,
            idx,
        });
    }
    Ok(params)
}

const V1: &str = include_str!("../../../../data/crystalff/torsion_v1.txt");
const V2: &str = include_str!("../../../../data/crystalff/torsion_v2.txt");
const SMALL_RINGS: &str = include_str!("../../../../data/crystalff/torsion_smallrings.txt");
const MACROCYCLES: &str = include_str!("../../../../data/crystalff/torsion_macrocycles.txt");
const GAUSS_V4: &str = include_str!("../../../../data/crystalff/gauss_v4.tsv");
const GAUSS_SMALL_RINGS: &str = include_str!("../../../../data/crystalff/gauss_smallrings.tsv");
const GAUSS_MACROCYCLES: &str = include_str!("../../../../data/crystalff/gauss_macrocycles.tsv");

/// The cosine-series preferences of ETKDG `version` 1 or 2, optionally with
/// the small-ring and macrocycle torsions.
pub fn cosine_torsion_params(
    version: u32,
    use_small_ring_torsions: bool,
    use_macrocycle_torsions: bool,
) -> Result<&'static [ExpTorsionAngle]> {
    let base = match version {
        1 => V1,
        2 => V2,
        _ => return Err(RDKitError::Value("ETversion must be 1 or 2.".into())),
    };
    static CACHE: [OnceLock<Vec<ExpTorsionAngle>>; 8] = [const { OnceLock::new() }; 8];
    let slot = (version as usize - 1) * 4
        + usize::from(use_small_ring_torsions) * 2
        + usize::from(use_macrocycle_torsions);
    let mut data = vec![base];
    if use_small_ring_torsions {
        data.push(SMALL_RINGS);
    }
    if use_macrocycle_torsions {
        data.push(MACROCYCLES);
    }
    if let Some(params) = CACHE[slot].get() {
        return Ok(params);
    }
    let parsed = parse_cosine(&data)?;
    Ok(CACHE[slot].get_or_init(|| parsed))
}

/// The Gaussian-fit preferences of ETKDG version 4.
pub fn gaussian_torsion_params(
    version: u32,
    use_small_ring_torsions: bool,
    use_macrocycle_torsions: bool,
) -> Result<&'static [GaussianExpTorsionAngle]> {
    if version != 4 {
        return Err(RDKitError::Value(
            "Currently, only ETVersions 4 is supported together with the Gaussian fit functional form."
                .into(),
        ));
    }
    static CACHE: [OnceLock<Vec<GaussianExpTorsionAngle>>; 4] = [const { OnceLock::new() }; 4];
    let slot = usize::from(use_small_ring_torsions) * 2 + usize::from(use_macrocycle_torsions);
    let mut data = vec![GAUSS_V4];
    if use_small_ring_torsions {
        data.push(GAUSS_SMALL_RINGS);
    }
    if use_macrocycle_torsions {
        data.push(GAUSS_MACROCYCLES);
    }
    if let Some(params) = CACHE[slot].get() {
        return Ok(params);
    }
    let parsed = parse_gaussian(&data)?;
    Ok(CACHE[slot].get_or_init(|| parsed))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn collections_parse() {
        let v2 = cosine_torsion_params(2, false, false).unwrap();
        assert!(v2.len() > 300);
        assert!(v2.iter().all(|p| p.signs.len() == 6 && p.v.len() == 6));
        // the first pattern is [O:1]=[C:2]!@;-[O:3]~[CH0:4]
        assert_eq!(v2[0].idx, [0, 1, 2, 3]);
        let all = cosine_torsion_params(2, true, true).unwrap();
        assert!(all.len() > v2.len() + 400);
        assert!(cosine_torsion_params(3, false, false).is_err());

        let g = gaussian_torsion_params(4, false, false).unwrap();
        assert!(g.len() > 300);
        assert!(g.iter().all(|p| {
            p.heights.len() == p.positions.len() && p.heights.len() == p.widths.len()
        }));
        assert!(gaussian_torsion_params(2, false, false).is_err());
        // macrocycle entries use PI expressions
        let all = gaussian_torsion_params(4, true, true).unwrap();
        assert!(all.iter().any(|p| {
            p.positions
                .iter()
                .any(|&x| (x - 2.0 * std::f64::consts::PI).abs() < 1e-12)
        }));
    }
}
