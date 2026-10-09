//! Ultrafast shape recognition descriptors (upstream
//! `GraphMol/Descriptors/USRDescriptor.cpp`): USR (Ballester and Richards,
//! Proc. R. Soc. A 2007) and the pharmacophore-aware USRCAT (Schreyer and
//! Blundell, J. Cheminf. 2012).
//!
//! The distances of all atoms to four reference points (centroid, closest and
//! farthest atom to it, farthest atom to that one) are summarized by their
//! first three moments: mean, standard deviation and the cubic root of the
//! skewness. USRCAT repeats this for the atoms of four feature classes
//! (hydrophobic, aromatic, acceptor and donor) with the same reference points.
//!
//! Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see
//! LICENSES/RDKit-ChemReactions-Novartis-BSD-3-Clause.txt).

use std::sync::OnceLock;

use crate::geometry::Point3D;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::fp::fmuladd;
use crate::{RDKitError, Result};

/// Number of descriptor values per atom selection (4 distributions x 3 moments).
pub const USR_BLOCK: usize = 12;

/// Feature definitions from <https://bitbucket.org/aschreyer/usrcat>.
const FEATURE_SMARTS: [&str; 4] = [
    // hydrophobic
    "[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]",
    // aromatic
    "[a]",
    // acceptor
    "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),\
$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]",
    // donor
    "[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]",
];

fn distances(coords: &[Point3D], point: Point3D) -> Vec<f64> {
    coords.iter().map(|c| (*c - point).length()).collect()
}

fn centroid(coords: &[Point3D]) -> Result<Point3D> {
    if coords.is_empty() {
        return Err(RDKitError::precondition("no coordinates"));
    }
    let mut pt = Point3D::new(0.0, 0.0, 0.0);
    for c in coords {
        pt += *c;
    }
    Ok(pt / coords.len() as f64)
}

/// Index of the first extreme value according to `better(candidate, best)`.
fn extreme_id(v: &[f64], better: impl Fn(f64, f64) -> bool) -> usize {
    let mut best = 0;
    for (i, &x) in v.iter().enumerate().skip(1) {
        if better(x, v[best]) {
            best = i;
        }
    }
    best
}

/// Mean, standard deviation and cubic root of the skewness of `dist`.
fn moments(dist: &[f64]) -> [f64; 3] {
    let mut m = [0.0; 3];
    let n = dist.len() as f64;
    if dist.is_empty() {
        return m;
    }
    for d in dist {
        m[0] += d;
    }
    m[0] /= n;
    for d in dist {
        let diff = d - m[0];
        m[1] = fmuladd(diff, diff, m[1]);
        m[2] = fmuladd(diff * diff, diff, m[2]);
    }
    m[1] = (m[1] / n).sqrt();
    m[2] /= n;
    m[2] = if m[1] == 0.0 {
        0.0
    } else {
        (m[2] / (m[1] * m[1] * m[1])).cbrt()
    };
    m
}

/// `calcUSRDistributions`: the four distance distributions and the reference
/// points they are measured from.
pub fn calc_usr_distributions(coords: &[Point3D]) -> Result<([Vec<f64>; 4], [Point3D; 4])> {
    let ctd = centroid(coords)?;
    let d_ctd = distances(coords, ctd);
    let catc = coords[extreme_id(&d_ctd, |x, best| x < best)];
    let fatc = coords[extreme_id(&d_ctd, |x, best| x > best)];
    let d_catc = distances(coords, catc);
    let d_fatc = distances(coords, fatc);
    let fatf = coords[extreme_id(&d_fatc, |x, best| x > best)];
    let d_fatf = distances(coords, fatf);
    Ok(([d_ctd, d_catc, d_fatc, d_fatf], [ctd, catc, fatc, fatf]))
}

/// `calcUSRDistributionsFromPoints`: the distances of `coords` to every point.
pub fn calc_usr_distributions_from_points(coords: &[Point3D], points: &[Point3D]) -> Vec<Vec<f64>> {
    points.iter().map(|p| distances(coords, *p)).collect()
}

/// `calcUSRFromDistributions`: three moments of every distribution.
pub fn calc_usr_from_distributions(dist: &[Vec<f64>]) -> Vec<f64> {
    dist.iter().flat_map(|d| moments(d)).collect()
}

/// `calcUSRScore`: similarity in (0, 1] of two descriptors made of blocks of
/// 12 values, one weight per block.
pub fn calc_usr_score(d1: &[f64], d2: &[f64], weights: &[f64]) -> Result<f64> {
    if d1.len() != d2.len() {
        return Err(RDKitError::precondition(
            "descriptors must have the same size",
        ));
    }
    if weights.len() != d1.len() / USR_BLOCK {
        return Err(RDKitError::precondition("size of weights not correct"));
    }
    let mut score = 1.0;
    for (w, (b1, b2)) in weights
        .iter()
        .zip(d1.chunks_exact(USR_BLOCK).zip(d2.chunks_exact(USR_BLOCK)))
    {
        let mut tmp = 0.0;
        for (a, b) in b1.iter().zip(b2) {
            tmp += (a - b).abs();
        }
        tmp /= USR_BLOCK as f64;
        score = fmuladd(*w, tmp, score);
    }
    Ok(1.0 / score)
}

/// Positions of the atoms of conformer `conf_id`, after the checks shared by
/// USR and USRCAT.
fn checked_coords(mol: &ROMol, conf_id: i32) -> Result<Vec<Point3D>> {
    if mol.num_atoms() < 3 {
        return Err(RDKitError::Value(
            "Number of atoms must be greater than 3".into(),
        ));
    }
    if mol.num_conformers() == 0 {
        return Err(RDKitError::Conformer(
            "No conformations available on this molecule".into(),
        ));
    }
    Ok(mol.conformer(conf_id)?.positions().to_vec())
}

/// `USR(mol, descriptor, confId)`: the 12 values of the whole molecule.
pub fn usr(mol: &ROMol, conf_id: i32) -> Result<Vec<f64>> {
    let coords = checked_coords(mol, conf_id)?;
    let (dist, _) = calc_usr_distributions(&coords)?;
    Ok(calc_usr_from_distributions(&dist))
}

/// The atoms of every feature class, found with the USRCAT SMARTS.
fn feature_atom_ids(mol: &ROMol) -> Result<Vec<Vec<u32>>> {
    static PATTERNS: OnceLock<Vec<ROMol>> = OnceLock::new();
    if PATTERNS.get().is_none() {
        let params = SmartsParserParams::default();
        let parsed = FEATURE_SMARTS
            .iter()
            .map(|sma| {
                mol_from_smarts(sma, &params)?
                    .ok_or_else(|| RDKitError::invariant("bad USRCAT feature SMARTS"))
            })
            .collect::<Result<Vec<_>>>()?;
        let _ = PATTERNS.set(parsed);
    }
    let params = SubstructMatchParameters::default();
    PATTERNS
        .get()
        .expect("patterns were just initialized")
        .iter()
        .map(|query| {
            let matches = substruct_match(mol, query, &params)?;
            Ok(matches
                .iter()
                .flat_map(|m| m.iter().map(|&(_, atom)| atom as u32))
                .collect())
        })
        .collect()
}

/// `USRCAT(mol, descriptor, atomIds, confId)`: the whole-molecule block
/// followed by one block per atom selection. Without `atom_ids` the four
/// default feature classes are used.
pub fn usrcat(mol: &ROMol, atom_ids: Option<&[Vec<u32>]>, conf_id: i32) -> Result<Vec<f64>> {
    let coords = checked_coords(mol, conf_id)?;
    let selections = match atom_ids {
        Some(ids) => ids.to_vec(),
        None => feature_atom_ids(mol)?,
    };
    let (dist, points) = calc_usr_distributions(&coords)?;
    let mut descriptor = calc_usr_from_distributions(&dist);
    for selection in &selections {
        let reduced = selection
            .iter()
            .map(|&idx| {
                coords.get(idx as usize).copied().ok_or_else(|| {
                    RDKitError::Index(format!("atom index {} out of range", idx as u64 + 1))
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let dist = calc_usr_distributions_from_points(&reduced, &points);
        descriptor.extend(calc_usr_from_distributions(&dist));
    }
    Ok(descriptor)
}
