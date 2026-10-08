//! Atropisomer perception (`RDKit::Atropisomers`).
//!
//! Ported from RDKit `Code/GraphMol/Atropisomers.cpp` (2026.09.1). Only the
//! detection half is ported; wedging from atropisomers comes with the
//! MolBlock writer.
//!
//! Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::BTreeSet;

use crate::Result;
use crate::geometry::Point3D;
use crate::graphmol::atom::HybridizationType;
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::conjug_hybrid::{set_conjugation, set_hybridization};
use crate::graphmol::romol::ROMol;
use crate::graphmol::stereo_group::StereoGroup;

const REALLY_SMALL_BOND_LEN: f64 = 0.0000001;

/// `(atom, neighbor bonds)` for one end of a candidate bond.
type AtropAtomAndBondVec = (u32, Vec<u32>);

fn is_atrop(s: BondStereo) -> bool {
    matches!(s, BondStereo::StereoAtropCW | BondStereo::StereoAtropCCW)
}

/// `canHaveDirection`
pub(crate) fn can_have_direction(t: BondType) -> bool {
    matches!(t, BondType::Single | BondType::Aromatic)
}

fn other(mol: &ROMol, bond: u32, atom: u32) -> u32 {
    let b = mol.bond(bond);
    if b.begin_atom_idx() == atom {
        b.end_atom_idx()
    } else {
        b.begin_atom_idx()
    }
}

/// `getAtropisomerAtomsAndBonds`
fn get_atropisomer_atoms_and_bonds(mol: &ROMol, bond: u32) -> Option<[AtropAtomAndBondVec; 2]> {
    let b = mol.bond(bond);
    let mut res = [
        (b.begin_atom_idx(), Vec::new()),
        (b.end_atom_idx(), Vec::new()),
    ];
    for end in res.iter_mut() {
        for nb in mol.atom_bond_indices(end.0) {
            if nb != bond {
                end.1.push(nb);
            }
        }
        if end.1.is_empty() {
            return None;
        }
        if end.1.len() == 2 && other(mol, end.1[1], end.0) < other(mol, end.1[0], end.0) {
            end.1.swap(0, 1);
        }
    }
    Some(res)
}

/// `getBondFrameOfReference`
fn get_bond_frame_of_reference(
    mol: &ROMol,
    bond: u32,
    conf: &Conformer,
) -> Result<Option<(Point3D, Point3D, Point3D)>> {
    let b = mol.bond(bond);
    let mut x = conf.atom_pos(b.end_atom_idx())? - conf.atom_pos(b.begin_atom_idx())?;
    if x.length() < REALLY_SMALL_BOND_LEN {
        return Ok(None);
    }
    x.normalize()?;
    if !conf.is_3d() {
        let mut y = Point3D::new(-x.y, x.x, 0.0);
        y.normalize()?;
        return Ok(Some((x, y, Point3D::new(0.0, 0.0, 1.0))));
    }
    let mut z = if x.x.abs() > REALLY_SMALL_BOND_LEN || x.y.abs() > REALLY_SMALL_BOND_LEN {
        Point3D::new(0.0, 0.0, 1.0)
    } else {
        Point3D::new(1.0, 0.0, 0.0)
    };
    let mut y = z.cross_product(&x);
    z = x.cross_product(&y);
    y.normalize()?;
    z.normalize()?;
    Ok(Some((x, y, z)))
}

/// `getAtropIsomerEndVect`
fn get_atrop_isomer_end_vect(
    mol: &ROMol,
    end: &AtropAtomAndBondVec,
    y: &Point3D,
    z: &Point3D,
    conf: &Conformer,
) -> Result<Option<Point3D>> {
    let center = conf.atom_pos(end.0)?;
    let v = conf.atom_pos(other(mol, end.1[0], end.0))? - center;
    let mut bond_vec = Point3D::new(0.0, v.dot_product(y), v.dot_product(z));
    if end.1.len() == 2 {
        let o = conf.atom_pos(other(mol, end.1[1], end.0))? - center;
        let other_vec = Point3D::new(0.0, o.dot_product(y), o.dot_product(z));
        if bond_vec.length() < REALLY_SMALL_BOND_LEN {
            bond_vec = -other_vec;
        } else if bond_vec.dot_product(&other_vec) > REALLY_SMALL_BOND_LEN {
            crate::rd_log_warning!(
                "Both bonds on one end of an atropisomer are on the same side - atoms is : {}",
                end.0
            );
            return Ok(None);
        }
    }
    if bond_vec.length() < REALLY_SMALL_BOND_LEN {
        crate::rd_log_warning!(
            "Could not find a bond on one end of an atropisomer that is not co-linear - atoms are : {}",
            end.0
        );
        return Ok(None);
    }
    bond_vec.normalize()?;
    Ok(Some(bond_vec))
}

/// `getBondDir`; `None` when the wedging on this end is inconsistent.
fn get_bond_dir(mol: &ROMol, bond: u32, end: &AtropAtomAndBondVec) -> Option<BondDir> {
    let wedge = |d: BondDir| {
        if matches!(d, BondDir::BeginWedge | BondDir::BeginDash) {
            d
        } else {
            BondDir::None
        }
    };
    let d1 = wedge(mol.bond(end.1[0]).bond_dir());
    let d2 = if end.1.len() == 2 {
        wedge(mol.bond(end.1[1]).bond_dir())
    } else {
        BondDir::None
    };
    if d1 != BondDir::None && d2 != BondDir::None && d1 == d2 {
        let b = mol.bond(bond);
        crate::rd_log_warning!(
            "The bonds on one end of an atropisomer are both UP or both DOWN - atoms are: {} {}",
            b.begin_atom_idx(),
            b.end_atom_idx()
        );
        return None;
    }
    if d1 == BondDir::BeginWedge || d2 == BondDir::BeginDash {
        return Some(BondDir::BeginWedge);
    }
    if d1 == BondDir::BeginDash || d2 == BondDir::BeginWedge {
        return Some(BondDir::BeginDash);
    }
    Some(BondDir::None)
}

/// `DetectAtropisomerChiralityOneBond`
fn detect_one_bond(mol: &mut ROMol, bond: u32, conf: Option<&Conformer>) -> Result<()> {
    let Some(ends) = get_atropisomer_atoms_and_bonds(mol, bond) else {
        return Ok(());
    };
    for end in &ends {
        for &eb in &end.1 {
            if mol.bond(eb).bond_dir() == BondDir::Unknown {
                return Ok(());
            }
        }
    }
    let (bb, be) = (
        mol.bond(bond).begin_atom_idx(),
        mol.bond(bond).end_atom_idx(),
    );

    let Some(conf) = conf else {
        let Some(d1) = get_bond_dir(mol, bond, &ends[0]) else {
            return Ok(());
        };
        let Some(d2) = get_bond_dir(mol, bond, &ends[1]) else {
            return Ok(());
        };
        if d1 == d2 {
            crate::rd_log_warning!(
                "inconsistent bond wedging for an atropisomer.  Atoms are: {bb} {be}"
            );
            return Ok(());
        }
        if d1 == BondDir::BeginWedge || d2 == BondDir::BeginDash {
            mol.bond_mut(bond).set_stereo(BondStereo::StereoAtropCCW)?;
        } else if d1 == BondDir::BeginDash || d2 == BondDir::BeginWedge {
            mol.bond_mut(bond).set_stereo(BondStereo::StereoAtropCW)?;
        }
        return Ok(());
    };

    let Some((_x, y, z)) = get_bond_frame_of_reference(mol, bond, conf)? else {
        crate::rd_log_warning!(
            "Failed to get a frame of reference along an atropisomer bond - atoms are: {bb} {be}"
        );
        return Ok(());
    };
    let mut bond_vecs = [Point3D::default(); 2];
    for (i, end) in ends.iter().enumerate() {
        if !conf.is_3d() {
            let Some(dir) = get_bond_dir(mol, bond, end) else {
                return Ok(());
            };
            let Some(v) = get_atrop_isomer_end_vect(mol, end, &y, &z, conf)? else {
                return Ok(());
            };
            bond_vecs[i] = v;
            if dir == BondDir::BeginWedge {
                bond_vecs[i].y *= 0.707;
                bond_vecs[i].z = bond_vecs[i].y.abs();
            } else if dir == BondDir::BeginDash {
                bond_vecs[i].y *= 0.707;
                bond_vecs[i].z = -bond_vecs[i].y.abs();
            }
        } else {
            let center = conf.atom_pos(end.0)?;
            let t = conf.atom_pos(other(mol, end.1[0], end.0))? - center;
            bond_vecs[i] = Point3D::new(0.0, t.dot_product(&y), t.dot_product(&z));
            if end.1.len() == 2 {
                let t = conf.atom_pos(other(mol, end.1[1], end.0))? - center;
                let other_vec = Point3D::new(0.0, t.dot_product(&y), t.dot_product(&z));
                if bond_vecs[i].length() < REALLY_SMALL_BOND_LEN {
                    bond_vecs[i] = -other_vec;
                } else if bond_vecs[i].dot_product(&other_vec) > REALLY_SMALL_BOND_LEN {
                    crate::rd_log_warning!(
                        "Both bonds on one end of an atropisomer are on the same side - atoms are: {bb} {be}"
                    );
                    return Ok(());
                }
            }
            if bond_vecs[i].length() < REALLY_SMALL_BOND_LEN {
                crate::rd_log_warning!(
                    "Failed to find a bond on one end of an atropisomer that is NOT co-linear - atoms are: {bb} {be}"
                );
                return Ok(());
            }
        }
    }
    let cp = bond_vecs[1].cross_product(&bond_vecs[0]);
    if cp.x > REALLY_SMALL_BOND_LEN {
        mol.bond_mut(bond).set_stereo(BondStereo::StereoAtropCCW)?;
    } else if cp.x < -REALLY_SMALL_BOND_LEN {
        mol.bond_mut(bond).set_stereo(BondStereo::StereoAtropCW)?;
    } else {
        crate::rd_log_warning!(
            "The 2 defining bonds for an atropisomer are co-planar - atoms are: {bb} {be}"
        );
    }
    Ok(())
}

/// `Atropisomers::cleanupAtropisomerStereoGroups`
pub fn cleanup_atropisomer_stereo_groups(mol: &mut ROMol) -> Result<()> {
    let mut new = Vec::new();
    for sg in mol.stereo_groups() {
        let mut ok_atoms = Vec::new();
        let mut ok_bonds: Vec<u32> = Vec::new();
        for &atom in sg.atoms() {
            let mut found = false;
            for b in mol.atom_bond_indices(atom) {
                if is_atrop(mol.bond(b).stereo()) {
                    found = true;
                    if !ok_bonds.contains(&b) {
                        ok_bonds.push(b);
                    }
                }
            }
            if !found {
                ok_atoms.push(atom);
            }
        }
        if ok_bonds.is_empty() {
            new.push(sg.clone());
        } else {
            new.push(StereoGroup::new(sg.group_type(), ok_atoms, ok_bonds, 0)?);
        }
    }
    mol.set_stereo_groups(new)
}

/// `Atropisomers::detectAtropisomerChirality(mol, conf, replaceExistingTags)`;
/// `conf` is the position of the conformer in `mol.conformers()`.
pub fn detect_atropisomer_chirality(
    mol: &mut ROMol,
    conf: Option<usize>,
    replace_existing_tags: bool,
) -> Result<()> {
    // upstream keeps a std::set<Bond*>; bonds are allocated in index order,
    // so index order is the closest deterministic equivalent.
    let mut bonds_to_try = BTreeSet::new();
    let conf_is_3d = conf.is_some_and(|c| mol.conformers()[c].is_3d());
    if replace_existing_tags && conf_is_3d {
        for b in mol.bonds() {
            if is_atrop(b.stereo()) {
                bonds_to_try.insert(b.idx());
            }
        }
    }
    for b in mol.bonds() {
        if can_have_direction(b.bond_type())
            && matches!(b.bond_dir(), BondDir::BeginDash | BondDir::BeginWedge)
        {
            for nb in mol.atom_bond_indices(b.begin_atom_idx()) {
                if nb == b.idx() {
                    continue;
                }
                if !replace_existing_tags && is_atrop(mol.bond(nb).stereo()) {
                    continue;
                }
                bonds_to_try.insert(nb);
            }
        }
    }
    if replace_existing_tags {
        for &b in &bonds_to_try {
            if is_atrop(mol.bond(b).stereo()) {
                mol.bond_mut(b).set_stereo(BondStereo::StereoNone)?;
            }
        }
    }
    if bonds_to_try.is_empty() {
        return Ok(());
    }
    let mut any = false;
    for &b in &bonds_to_try {
        let (bg, en) = (mol.bond(b).begin_atom_idx(), mol.bond(b).end_atom_idx());
        for a in [bg, en] {
            if mol.atom(a).needs_update_property_cache() {
                mol.update_atom_property_cache(a, false)?;
            }
        }
        let bond = mol.bond(b);
        let (bd, ed) = (mol.total_degree(bg)?, mol.total_degree(en)?);
        if bond.bond_type() == BondType::Single
            && bond.stereo() != BondStereo::StereoAny
            && (2..=3).contains(&bd)
            && (2..=3).contains(&ed)
        {
            any = true;
            break;
        }
    }
    if !any {
        return Ok(());
    }
    let needs_update = mol.needs_update_property_cache()
        || mol
            .atoms()
            .iter()
            .any(|a| a.atomic_num() != 0 && a.hybridization() == HybridizationType::Unspecified);
    if needs_update {
        mol.update_property_cache(false)?;
        set_conjugation(mol);
        set_hybridization(mol);
    }
    let conf_copy = conf.map(|c| mol.conformers()[c].clone());
    for &b in &bonds_to_try {
        let bond = mol.bond(b);
        if bond.bond_type() != BondType::Single
            || bond.stereo() == BondStereo::StereoAny
            || mol.atom(bond.begin_atom_idx()).hybridization() != HybridizationType::SP2
            || mol.atom(bond.end_atom_idx()).hybridization() != HybridizationType::SP2
        {
            continue;
        }
        detect_one_bond(mol, b, conf_copy.as_ref())?;
    }
    Ok(())
}

// ------------------------------------------------------------------ wedging

use crate::graphmol::file_parsers::mol_file_writer::{WedgeBonds, WedgeInfo};

/// `getBondDirForAtropisomerNoConf`
fn bond_dir_no_conf(stereo: BondStereo, which_end: usize, which_bond: usize) -> BondDir {
    let mut flips = 0;
    if stereo == BondStereo::StereoAtropCW {
        flips += 1;
    }
    if which_bond == 1 {
        flips += 1;
    }
    if which_end == 1 {
        flips += 1;
    }
    if flips % 2 == 1 {
        BondDir::BeginDash
    } else {
        BondDir::BeginWedge
    }
}

/// `getBondDirForAtropisomer2d`
fn bond_dir_2d(
    vecs: &[Point3D; 2],
    stereo: BondStereo,
    which_end: usize,
    which_bond: usize,
) -> BondDir {
    let mut flips = 0;
    if stereo == BondStereo::StereoAtropCCW {
        flips += 1;
    }
    if which_bond == 1 {
        flips += 1;
    }
    if which_end == 1 {
        flips += 1;
    }
    if vecs[1 - which_end].y < 0.0 {
        flips += 1;
    }
    if flips % 2 == 1 {
        BondDir::BeginWedge
    } else {
        BondDir::BeginDash
    }
}

/// `getBondDirForAtropisomer3d`
fn bond_dir_3d(mol: &ROMol, bond: u32, conf: &Conformer) -> Result<BondDir> {
    let b = mol.bond(bond);
    let dz = conf.atom_pos(b.end_atom_idx())?.z - conf.atom_pos(b.begin_atom_idx())?.z;
    Ok(if dz > REALLY_SMALL_BOND_LEN {
        BondDir::BeginWedge
    } else {
        BondDir::BeginDash
    })
}

/// `bond->setEndAtom(begin); bond->setBeginAtom(atom)`: only the bond's
/// own begin/end indices change (the graph is untouched).
fn make_begin(mol: &mut ROMol, bond: u32, atom: u32) {
    let b = mol.bond_mut(bond);
    if b.begin_atom_idx() != atom {
        let old = b.begin_atom_idx();
        b.set_end_atom_idx_unchecked(old);
        b.set_begin_atom_idx_unchecked(atom);
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Geom {
    NoConf,
    TwoD,
    ThreeD,
}

/// `WedgeBondFromAtropisomerOneBond{NoConf,2d,3d}`
fn wedge_one_bond(
    mol: &mut ROMol,
    bond: u32,
    conf: Option<&Conformer>,
    geom: Geom,
    wedge: &mut WedgeBonds,
) -> Result<bool> {
    let Some(ends) = get_atropisomer_atoms_and_bonds(mol, bond) else {
        return Ok(false);
    };
    for end in &ends {
        for &eb in &end.1 {
            if mol.bond(eb).bond_dir() == BondDir::Unknown {
                return Ok(false);
            }
        }
    }
    let (bb, be) = (
        mol.bond(bond).begin_atom_idx(),
        mol.bond(bond).end_atom_idx(),
    );
    let stereo = mol.bond(bond).stereo();
    let mut vecs = [Point3D::default(); 2];
    if geom == Geom::TwoD {
        let conf = conf.expect("2D conformer");
        let Some((_x, y, z)) = get_bond_frame_of_reference(mol, bond, conf)? else {
            crate::rd_log_warning!(
                "Cound not get a frame of reference for an atropisomer bond - atoms are: {bb} {be}"
            );
            return Ok(false);
        };
        for (i, end) in ends.iter().enumerate() {
            let Some(v) = get_atrop_isomer_end_vect(mol, end, &y, &z, conf)? else {
                return Ok(false);
            };
            if v.length() < REALLY_SMALL_BOND_LEN {
                crate::rd_log_warning!(
                    "Failed to get a representative vector for the defining bond of an atropisomer - atoms are: {bb} {be}"
                );
                return Ok(false);
            }
            vecs[i] = v;
        }
    }
    let atrop_type = mol.bond(bond).bond_type();

    // bonds that already carry a wedge from the atropisomer atoms
    let mut use_bonds: Vec<(usize, usize)> = Vec::new();
    for (we, end) in ends.iter().enumerate() {
        for (wb, &eb) in end.1.iter().enumerate() {
            let b = mol.bond(eb);
            // upstream checks canHaveDirection of the atropisomer bond for
            // NoConf/2d and of the neighbor bond for 3d
            let dir_ok = if geom == Geom::ThreeD {
                can_have_direction(b.bond_type())
            } else {
                can_have_direction(atrop_type)
            };
            if matches!(b.bond_dir(), BondDir::BeginWedge | BondDir::BeginDash)
                && b.begin_atom_idx() == end.0
                && dir_ok
            {
                use_bonds.push((we, wb));
            }
        }
    }
    if !use_bonds.is_empty() {
        for (we, wb) in use_bonds {
            let eb = ends[we].1[wb];
            let d = match geom {
                Geom::NoConf => bond_dir_no_conf(stereo, we, wb),
                Geom::TwoD => bond_dir_2d(&vecs, stereo, we, wb),
                Geom::ThreeD => bond_dir_3d(mol, eb, conf.expect("3D conformer"))?,
            };
            mol.bond_mut(eb).set_bond_dir(d);
        }
        return Ok(true);
    }

    let mut best: Option<(usize, u32)> = None;
    let mut best_single = false;
    let mut best_ring_count = if geom == Geom::ThreeD {
        u32::MAX
    } else {
        i32::MAX as u32
    };
    let mut largest_ring = 0u32;
    let mut best_dir = BondDir::None;
    for (we, end) in ends.iter().enumerate() {
        for (wb, &bt) in end.1.iter().enumerate() {
            let skip_key = if geom == Geom::ThreeD { bond } else { bt };
            if !can_have_direction(mol.bond(bt).bond_type()) || wedge.contains_key(&skip_key) {
                continue;
            }
            if geom == Geom::ThreeD {
                make_begin(mol, bt, end.0);
            }
            let b = mol.bond(bt);
            if b.bond_dir() != BondDir::None {
                if b.begin_atom_idx() == end.0 {
                    if geom == Geom::TwoD
                        && !matches!(b.bond_dir(), BondDir::BeginWedge | BondDir::BeginDash)
                    {
                        continue;
                    }
                    crate::rd_log_warning!(
                        "Wedge or hash bond found on atropisomer where not expected - atoms are: {bb} {be}"
                    );
                    return Ok(false);
                }
                continue;
            }
            let is_single = b.bond_type() == BondType::Single;
            let mut ring_count = mol.ring_info().num_bond_rings(bt)?;
            let mut ring_size = 0;
            if geom != Geom::NoConf {
                if ring_count == 0 {
                    ring_count = 10;
                } else {
                    ring_size = mol.ring_info().min_bond_ring_size(bt)?;
                    if ring_size > 8 {
                        ring_size = 0;
                    }
                }
            }
            let dir_for = |mol: &ROMol| -> Result<BondDir> {
                Ok(match geom {
                    Geom::NoConf => bond_dir_no_conf(stereo, we, wb),
                    Geom::TwoD => bond_dir_2d(&vecs, stereo, we, wb),
                    Geom::ThreeD => bond_dir_3d(mol, bt, conf.expect("3D conformer"))?,
                })
            };
            if ring_count > best_ring_count {
                continue;
            } else if ring_count < best_ring_count
                || (geom != Geom::NoConf && ring_size > largest_ring)
            {
                best = Some((we, bt));
                best_ring_count = ring_count;
                largest_ring = ring_size;
                best_single = is_single;
                best_dir = dir_for(mol)?;
            } else if best_single && !is_single {
                continue;
            } else if !best_single && is_single {
                best = Some((we, bt));
                best_ring_count = ring_count;
                best_single = true;
                best_dir = dir_for(mol)?;
            } else {
                let d = dir_for(mol)?;
                if best_dir == BondDir::None
                    || (best_dir == BondDir::BeginDash && d == BondDir::BeginWedge)
                {
                    best = Some((we, bt));
                    best_ring_count = ring_count;
                    best_single = is_single;
                    best_dir = d;
                }
            }
        }
    }
    match best {
        Some((we, bt)) => {
            make_begin(mol, bt, ends[we].0);
            mol.bond_mut(bt).set_bond_dir(best_dir);
            wedge.insert(bt, WedgeInfo::Atropisomer(best_dir));
            Ok(true)
        }
        None => {
            crate::rd_log_warning!(
                "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are: {bb} {be}"
            );
            Ok(false)
        }
    }
}

/// `Atropisomers::wedgeBondsFromAtropisomers`
pub fn wedge_bonds_from_atropisomers(
    mol: &mut ROMol,
    conf: Option<&Conformer>,
    wedge: &mut WedgeBonds,
) -> Result<()> {
    if !mol.ring_info().is_sssr_or_better() {
        crate::graphmol::find_rings::find_sssr(mol, false, false)?;
    }
    for bi in 0..mol.num_bonds() {
        let b = mol.bond(bi);
        if b.bond_type() != BondType::Single || !is_atrop(b.stereo()) {
            continue;
        }
        let (bd, ed) = (
            mol.total_degree(b.begin_atom_idx())?,
            mol.total_degree(b.end_atom_idx())?,
        );
        if !(2..=3).contains(&bd) || !(2..=3).contains(&ed) {
            continue;
        }
        let geom = match conf {
            Some(c) if c.is_3d() => Geom::ThreeD,
            Some(_) => Geom::TwoD,
            None => Geom::NoConf,
        };
        wedge_one_bond(mol, bi, conf, geom, wedge)?;
    }
    Ok(())
}
