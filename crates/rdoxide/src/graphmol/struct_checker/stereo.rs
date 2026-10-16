//! Stereochemistry checks, atom clashes and cis/trans perception.
//!
//! Ported from RDKit `Code/GraphMol/StructChecker/Stereo.cpp` (2026.09.1).
//! Upstream reports problems on `std::cerr`, which the `eprintln!`s keep.
//!
//! Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use super::utilities::{Neighbourhood, get_mol_atom_points, setup_neighbourhood};
use crate::Result;
use crate::geometry::Point3D;
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondDir, BondType};
use crate::graphmol::mol_ops::remove_stereochemistry;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::{RDProps, RDValue, common_properties as cp};

pub const ODD: i32 = 1;
pub const EVEN: i32 = 2;
pub const UNMARKED: i32 = 3;

pub const ALLENE_PARITY: i32 = -2;
pub const ILLEGAL_REPRESENTATION: i32 = -1;
pub const UNDEFINED_PARITY: i32 = 0;
pub const ODD_PARITY: i32 = 1;
pub const EVEN_PARITY: i32 = 2;

pub fn invert_parity(p: i32) -> i32 {
    if p == 0 { 0 } else { 3 - p }
}

pub const CIS: u32 = 1;
pub const TRANS: u32 = 2;

pub const EITHER_BOND_FOUND: i32 = 1;
pub const STEREO_BOND_AT_NON_STEREO_ATOM: i32 = 2;
pub const ZEROED_Z_COORDINATES: i32 = 4;
pub const CONVERTED_TO_2D: i32 = 8;

// (upstream's own approximation)
#[allow(clippy::approx_constant)]
const PI: f64 = 3.14159265359;
const ANGLE_EPSILON: f64 = 5.0 * PI / 180.0;
const EPS: f64 = 0.0000001;

#[derive(Debug, Clone, Copy, Default)]
struct NPoint {
    x: f64,
    y: f64,
    z: f64,
    number: i32,
}

/// The angle between the vectors (x1,y1) and (x2,y2) at (0,0).
fn angle(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    let l1 = (x1 * x1 + y1 * y1).sqrt();
    let l2 = (x2 * x2 + y2 * y2).sqrt();
    if l1 < 0.00001 || l2 < 0.00001 {
        return 0.0;
    }
    // safeguard against round off errors
    let cos_alpha = ((x1 * x2 + y1 * y2) / (l1 * l2)).clamp(-1.0, 1.0);
    let sin_alpha = (x1 * y2 - x2 * y1) / (l1 * l2);
    let mut result = cos_alpha.acos();
    if sin_alpha < 0.0 {
        result = 2.0 * PI - result;
    }
    result
}

/// The signed volume of the tetrahedron defined by four points.
fn volume(t: &[NPoint; 4]) -> f64 {
    let (ax, ay, az) = (t[1].x - t[0].x, t[1].y - t[0].y, t[1].z - t[0].z);
    let (bx, by, bz) = (t[2].x - t[0].x, t[2].y - t[0].y, t[2].z - t[0].z);
    let (cx, cy, cz) = (t[3].x - t[0].x, t[3].y - t[0].y, t[3].z - t[0].z);
    ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx)
}

fn mol_name(mol: &ROMol) -> String {
    match mol.get_prop(cp::_Name) {
        Some(RDValue::String(s)) => s.clone(),
        Some(v) => v.to_rd_string(),
        None => String::new(),
    }
}

fn is_wedge_or_dash(d: BondDir) -> bool {
    d == BondDir::BeginWedge || d == BondDir::BeginDash
}

/// Is the atom a potential stereocentre (for the dubious stereo bonds)?
fn stereogenic(element: u32, n_ligands: usize, nmulti: u32, allene: Option<bool>) -> bool {
    (element == 6 && n_ligands > 2 && n_ligands <= 4 && nmulti == 0)
        || allene.is_some_and(|is_allene| {
            element == 6 && (2..=3).contains(&n_ligands) && nmulti == 1 && is_allene
        })
        || (element == 16 && n_ligands > 2 && n_ligands <= 4)
        || (element == 7 && n_ligands > 3 && n_ligands <= 4 && nmulti == 0)
        || (element == 15 && n_ligands > 2 && n_ligands <= 4)
        || (element == 14 && n_ligands > 2 && n_ligands <= 4 && nmulti == 0)
}

/// `DubiousStereochemistry`: checks for ill-defined stereochemistry; the
/// result is a set of the `*_FOUND`/`*_AT_NON_STEREO_ATOM` bits.
pub fn dubious_stereochemistry(mol: &ROMol) -> i32 {
    let mut result = 0;
    let neighbour_array = setup_neighbourhood(mol);

    // look for EITHER bonds
    for i in 0..mol.num_bonds() {
        if mol.bond(i).bond_dir() == BondDir::Unknown {
            eprintln!("bond {i}: Either bond found");
            result |= EITHER_BOND_FOUND;
        }
    }
    // look for stereo bonds to non-stereogenic atoms
    for (i, nbp) in neighbour_array.iter().enumerate() {
        let mut nmulti = 0;
        let mut is_allene = false;
        for j in 0..nbp.bonds.len() {
            let bond = mol.bond(nbp.bonds[j]);
            if bond.bond_type() != BondType::Single {
                let mut ndb = 0;
                nmulti += 1;
                if bond.bond_type() == BondType::Double {
                    let jatom = nbp.atoms[j] as usize;
                    for &jj in &neighbour_array[jatom].bonds {
                        if mol.bond(jj).bond_type() == BondType::Double {
                            ndb += 1;
                        }
                    }
                }
                if ndb == 2 {
                    is_allene = true;
                }
            }
        }
        let element = mol.atom(i as u32).atomic_num();
        let n_ligands = nbp.bonds.len();
        if !stereogenic(element, n_ligands, nmulti, Some(is_allene)) {
            for j in 0..n_ligands {
                let bj = mol.bond(nbp.bonds[j]);
                if bj.begin_atom_idx() as usize == i && is_wedge_or_dash(bj.bond_dir()) {
                    crate::rd_log_warning!(
                        "{:>10}    atom {:>3} : {}\n",
                        mol_name(mol),
                        i,
                        "stereobond to non-stereogenic atom"
                    );
                    result |= STEREO_BOND_AT_NON_STEREO_ATOM;
                }
            }
        }
    }
    result
}

/// `FixDubious3DMolecule`: checks if the structure has 3D coordinates and/or
/// flat sp3 centres with stereo bonds and converts the designation to 2D.
/// Real 3D structures without stereo designations go through untouched.
pub fn fix_dubious_3d_molecule(mol: &mut ROMol) -> Result<i32> {
    let mut result = 0;
    let (atom_point, non_zero_z) = get_mol_atom_points(mol, false)?;
    // At first check if this is a trivial case i.e. designated '2D'
    if !non_zero_z {
        return Ok(0);
    }
    let nstereo = (0..mol.num_bonds())
        .filter(|&i| is_wedge_or_dash(mol.bond(i).bond_dir()))
        .count();
    if nstereo == 0 {
        return Ok(0);
    }

    // average bond length to use in the significance testing of the volume
    let mut length = 0.0;
    for i in 0..mol.num_bonds() {
        let bond = mol.bond(i);
        let a0 = atom_point[bond.begin_atom_idx() as usize];
        let a1 = atom_point[bond.end_atom_idx() as usize];
        length += (a0.x - a1.x) * (a0.x - a1.x)
            + (a0.y - a1.y) * (a0.y - a1.y)
            + (a0.z - a1.z) * (a0.z - a1.z);
    }
    length = (length / f64::from(mol.num_bonds())).sqrt();

    // check if there is a flat sp3 carbon
    let neighbour_array = setup_neighbourhood(mol);
    let mut nflat_sp3 = 0;
    for (i, nbp) in neighbour_array.iter().enumerate() {
        if nbp.atoms.len() < 3 {
            continue;
        }
        let element = mol.atom(i as u32).atomic_num();
        // (upstream tests the element with a condition that is never true)

        let any_designation = (0..mol.num_bonds()).any(|j| {
            let bond = mol.bond(j);
            bond.bond_dir() == BondDir::BeginWedge
                || (bond.bond_dir() == BondDir::BeginDash && i as u32 == bond.begin_atom_idx())
        });
        if any_designation {
            continue; // no stereo designation
        }

        let n_ligands = nbp.bonds.len();
        let mut tetra = [NPoint::default(); 4];
        tetra[0] = NPoint {
            x: atom_point[i].x,
            y: atom_point[i].y,
            z: atom_point[i].z,
            number: 0,
        };
        let mut i1 = 0;
        while i1 < n_ligands {
            if mol.bond(nbp.bonds[i1]).bond_type() != BondType::Single
                && element != 16
                && element != 15
            {
                break;
            }
            i1 += 1;
        }
        if i1 >= n_ligands {
            continue; // upstream: "multiple bond found => no sp3 carbon"
        }
        let point = |k: usize| NPoint {
            x: atom_point[k].x,
            y: atom_point[k].y,
            z: atom_point[k].z,
            number: 0,
        };
        let stereo = |k: usize| is_wedge_or_dash(mol.bond(nbp.bonds[k]).bond_dir());
        let mut stereo_triple = 0;
        let mut i1 = 0;
        while i1 < n_ligands {
            tetra[1] = point(i1);
            if stereo(i1) {
                stereo_triple |= 1;
            }
            let mut i2 = i1 + 1;
            while i2 < n_ligands {
                tetra[2] = point(i2);
                if stereo(i2) {
                    stereo_triple |= 2;
                }
                let mut i3 = i2 + 1;
                while i3 < n_ligands {
                    tetra[3] = point(i3);
                    if stereo(i3) {
                        stereo_triple |= 4;
                    }
                    let vol = volume(&tetra).abs();
                    if stereo_triple == 0 {
                        i3 += 1;
                        continue;
                    }
                    if vol < 0.01 * length * length * length {
                        nflat_sp3 += 1;
                        break;
                    }
                    stereo_triple &= !4;
                    i3 += 1;
                }
                stereo_triple &= !2;
                if i3 >= n_ligands {
                    break;
                }
                i2 += 1;
            }
            stereo_triple &= !1;
            if i2 >= n_ligands {
                break;
            }
            i1 += 1;
        }
    }

    // (upstream zeroes the z of its local copy of the coordinates only)
    result |= ZEROED_Z_COORDINATES;
    if non_zero_z && nstereo > 0 && nflat_sp3 > 0 {
        remove_stereochemistry(mol)?;
        result |= CONVERTED_TO_2D;
    }
    Ok(result)
}

/// `RemoveDubiousStereochemistry`: removes ill-defined stereo descriptors.
pub fn remove_dubious_stereochemistry(mol: &mut ROMol) {
    let neighbour_array = setup_neighbourhood(mol);
    for i in 0..mol.num_bonds() {
        if mol.bond(i).bond_dir() == BondDir::Unknown {
            mol.bond_mut(i).set_bond_dir(BondDir::None);
        }
    }
    // remove the stereo marks to non-stereogenic atoms
    for (i, nbp) in neighbour_array.iter().enumerate() {
        let nmulti = nbp
            .bonds
            .iter()
            .filter(|&&b| mol.bond(b).bond_type() != BondType::Single)
            .count() as u32;
        let element = mol.atom(i as u32).atomic_num();
        let n_ligands = nbp.bonds.len();
        if !stereogenic(element, n_ligands, nmulti, None) {
            for &b in &nbp.bonds {
                let bj = mol.bond(b);
                if bj.begin_atom_idx() as usize == i && is_wedge_or_dash(bj.bond_dir()) {
                    mol.bond_mut(b).set_bond_dir(BondDir::None);
                }
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct StereoBond {
    x: f64,
    y: f64,
    direction: BondDir,
    number: i32,
    angle: f64,
}

impl Default for StereoBond {
    fn default() -> Self {
        StereoBond {
            x: 0.0,
            y: 0.0,
            direction: BondDir::None,
            number: 0,
            angle: 0.0,
        }
    }
}

/// Insertion sort of a tetrahedron by the atom numbers.
fn sort_tetrahedron(t: &mut [NPoint; 4]) {
    for i in 1..4 {
        for j in (1..=i).rev() {
            if t[j].number < t[j - 1].number {
                t.swap(j, j - 1);
            } else {
                break;
            }
        }
    }
}

/// The stereo parity defined by three ligands.
fn atom3_parity(ligands: &[StereoBond; 4]) -> i32 {
    let maxnum = ligands[..3].iter().map(|l| l.number).max().unwrap_or(0);
    let mut reference = None;
    for (i, l) in ligands[..3].iter().enumerate() {
        if l.direction != BondDir::None {
            if reference.is_none() {
                reference = Some(i);
            } else {
                eprintln!("three attachments with more than 2 stereobonds");
                return ILLEGAL_REPRESENTATION;
            }
        }
    }
    let Some(reference) = reference else {
        return UNDEFINED_PARITY;
    };
    let (a, b) = match reference {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    };
    let ang = angle(ligands[a].x, ligands[a].y, ligands[b].x, ligands[b].y);
    if ang < ANGLE_EPSILON || (PI - ang).abs() < ANGLE_EPSILON {
        eprintln!("three attachments colinearity violation");
        return ILLEGAL_REPRESENTATION;
    }
    let mut tetrahedron = [NPoint::default(); 4];
    tetrahedron[0].number = maxnum + 1;
    for i in 0..3 {
        tetrahedron[i + 1].x = ligands[i].x;
        tetrahedron[i + 1].y = ligands[i].y;
        tetrahedron[i + 1].z = match ligands[i].direction {
            BondDir::BeginWedge => 1.0,
            BondDir::BeginDash => -1.0,
            BondDir::None => 0.0,
            _ => {
                eprintln!("three attachments illegal bond symbol");
                return ILLEGAL_REPRESENTATION;
            }
        };
        tetrahedron[i + 1].number = ligands[i].number;
    }
    sort_tetrahedron(&mut tetrahedron);
    if volume(&tetrahedron) > 0.0 {
        EVEN_PARITY
    } else {
        ODD_PARITY
    }
}

/// The stereo parity defined by four ligands; the central atom is at 0/0/0.
fn atom4_parity(ligands: &[StereoBond; 4]) -> i32 {
    let mut tetrahedron = [NPoint::default(); 4];
    let (mut nup, mut ndown) = (0, 0);
    for i in 0..4 {
        tetrahedron[i].x = ligands[i].x;
        tetrahedron[i].y = ligands[i].y;
        tetrahedron[i].z = 0.0;
        tetrahedron[i].number = ligands[i].number;
        match ligands[i].direction {
            BondDir::BeginWedge => {
                nup += 1;
                tetrahedron[i].z = 1.0;
            }
            BondDir::BeginDash => {
                ndown += 1;
                tetrahedron[i].z = -1.0;
            }
            BondDir::None => {}
            _ => {
                eprintln!("illegal bond symbol");
                return ILLEGAL_REPRESENTATION;
            }
        }
    }
    if nup == 0 && ndown == 0 {
        return UNDEFINED_PARITY;
    }
    if nup > 2 || ndown > 2 {
        eprintln!("too many stereobonds");
        return ILLEGAL_REPRESENTATION;
    }
    if nup + ndown == 1 {
        // check for 'umbrellas'
        let ij = ligands
            .iter()
            .position(|l| is_wedge_or_dash(l.direction))
            .unwrap_or(4);
        let mut nopposite = 0;
        for j in 0..4 {
            if ij == j {
                continue;
            } else if ligands[ij].x * ligands[j].x + ligands[ij].y * ligands[j].y < 0.0 {
                nopposite += 1;
            }
        }
        if nopposite > 2 {
            eprintln!("umbrella");
            return ILLEGAL_REPRESENTATION;
        }
    }
    for i in 0..2 {
        let (d1, d2) = (ligands[i].direction, ligands[i + 2].direction);
        if (d1 == BondDir::BeginWedge && d2 == BondDir::BeginDash)
            || (d1 == BondDir::BeginDash && d2 == BondDir::BeginWedge)
        {
            eprintln!("up/down");
            return ILLEGAL_REPRESENTATION;
        }
    }
    for i in 0..4 {
        let (d1, d2) = (ligands[i].direction, ligands[(i + 1) % 4].direction);
        if (d1 == BondDir::BeginWedge && d2 == BondDir::BeginWedge)
            || (d1 == BondDir::BeginDash && d2 == BondDir::BeginDash)
        {
            eprintln!("adjacent like");
            return ILLEGAL_REPRESENTATION;
        }
    }
    for i in 0..4 {
        if ligands[i].direction == BondDir::None
            && ligands[(i + 1) % 4].direction == BondDir::None
            && ligands[(i + 2) % 4].direction == BondDir::None
        {
            let ang = angle(
                ligands[i].x - ligands[(i + 1) % 4].x,
                ligands[i].y - ligands[(i + 1) % 4].y,
                ligands[(i + 2) % 4].x - ligands[(i + 1) % 4].x,
                ligands[(i + 2) % 4].y - ligands[(i + 1) % 4].y,
            );
            if ang < 185.0 * PI / 180.0 {
                eprintln!("colinearity or triangle rule");
                return ILLEGAL_REPRESENTATION;
            }
        }
    }
    sort_tetrahedron(&mut tetrahedron);
    if volume(&tetrahedron) > 0.0 {
        EVEN_PARITY
    } else {
        ODD_PARITY
    }
}

/// `AtomParity`: the stereo parity of the atom relative to its numbering;
/// `nbp` are its immediate neighbours.
///
/// Upstream mixes the indices of the neighbours with those of the molecule
/// in several places (the bond and the coordinates of ligand i are those of
/// bond/atom i of the molecule); that is kept.
pub fn atom_parity(mol: &ROMol, iatom: u32, nbp: &Neighbourhood) -> Result<i32> {
    let mut ligands = [StereoBond::default(); 4];
    let (mut multiple, mut allene, mut stereo) = (false, false, false);
    if nbp.atoms.len() < 3 || nbp.atoms.len() > 4 {
        return Ok(ILLEGAL_REPRESENTATION);
    }
    let (atom_point, _) = get_mol_atom_points(mol, false)?;
    for i in 0..nbp.bonds.len() {
        let bi = mol.bond(i as u32);
        if bi.bond_type() != BondType::Single {
            multiple = true;
            // check if the multiple bond is part of an allene
            let jatom = nbp.atoms[i] + 1;
            let mut ndb = 0;
            for j in 0..mol.num_bonds() {
                let bond = mol.bond(j);
                if (bond.begin_atom_idx() == jatom || bond.end_atom_idx() == jatom)
                    && bond.bond_type() == BondType::Double
                {
                    ndb += 1;
                }
            }
            if ndb == 2 {
                allene = true;
            }
        }
        ligands[i].x = atom_point[i].x - atom_point[iatom as usize].x;
        ligands[i].y = atom_point[i].y - atom_point[iatom as usize].y;
        ligands[i].number = nbp.atoms[i] as i32 + 1;
        if bi.begin_atom_idx() == iatom {
            ligands[i].direction = bi.bond_dir();
            if is_wedge_or_dash(ligands[i].direction) {
                stereo = true;
            }
        } else {
            ligands[i].direction = BondDir::None;
        }
    }
    let element = mol.atom(iatom).atomic_num();
    if multiple && stereo && element != 15 && element != 16 {
        if allene {
            return Ok(ALLENE_PARITY);
        }
        eprintln!("stereobond at unsaturated atom");
        return Ok(ILLEGAL_REPRESENTATION);
    } else if multiple && element != 16 {
        return Ok(UNDEFINED_PARITY);
    }

    ligands[0].angle = 0.0; // the angle relative to the first ligand
    for i in 1..nbp.atoms.len() {
        ligands[i].angle = angle(ligands[0].x, ligands[0].y, ligands[i].x, ligands[i].y);
    }
    for i in 2..nbp.atoms.len() {
        // sort the ligands
        for j in (2..=i).rev() {
            if ligands[j].angle < ligands[j - 1].angle {
                ligands.swap(j, j - 1);
            } else {
                break;
            }
        }
    }
    Ok(if nbp.atoms.len() == 3 {
        atom3_parity(&ligands)
    } else {
        atom4_parity(&ligands)
    })
}

/// `CheckStereo`: checks if all potential stereocentres are either
/// completely undefined or attributed with hashes and wedges according to
/// the MDL rules.
pub fn check_stereo(mol: &ROMol) -> Result<bool> {
    let mut result = true;
    let mut center_defined = false;
    let neighbour_array = setup_neighbourhood(mol);
    for (i, nbp) in neighbour_array.iter().enumerate() {
        let element = mol.atom(i as u32).atomic_num();
        let n_ligands = nbp.bonds.len();
        if n_ligands > 2 && n_ligands <= 4 && matches!(element, 6 | 16 | 7 | 8 | 15 | 14) {
            let parity = atom_parity(mol, i as u32, nbp)?;
            if matches!(
                mol.atom(i as u32).chiral_tag(),
                ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW
            ) {
                center_defined = true;
            }
            if parity == ILLEGAL_REPRESENTATION {
                result = false;
            } else if parity == EVEN_PARITY || parity == ODD_PARITY || parity == ALLENE_PARITY {
                center_defined = true;
            } else {
                for j in 0..nbp.bonds.len() {
                    let bond = mol.bond(j as u32);
                    if bond.begin_atom_idx() as usize == i && is_wedge_or_dash(bond.bond_dir()) {
                        eprintln!("stereobond to nonstereogenic");
                        result = false;
                    }
                }
            }
        }
    }
    let chiral_flag = mol
        .get_prop(cp::_MolFileChiralFlag)
        .and_then(RDValue::as_uint)
        .unwrap_or(0);
    if chiral_flag != 0 && !center_defined {
        eprintln!("chiral flag, no stereocenter");
        result = false;
    }
    Ok(result)
}

/// `AtomClash`: checks if any two atoms come closer than `clash_limit` of
/// the median bond length or if an atom is too close to the line between
/// two bonded atoms.
pub fn atom_clash(mol: &ROMol, clash_limit: f64) -> Result<bool> {
    let (atom_point, _) = get_mol_atom_points(mol, true)?;
    if mol.num_bonds() == 0 {
        return Ok(false);
    }
    let sq = |p: &Point3D, q: &Point3D| (p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y);
    // the median of the square of the bond length (quick/dirty)
    let mut blengths: Vec<f64> = (0..mol.num_bonds())
        .map(|i| {
            let bond = mol.bond(i);
            sq(
                &atom_point[bond.begin_atom_idx() as usize],
                &atom_point[bond.end_atom_idx() as usize],
            )
        })
        .collect();
    for i in 1..blengths.len() {
        for j in (0..i).rev() {
            if blengths[j] > blengths[j + 1] {
                blengths.swap(j, j + 1);
            } else {
                break;
            }
        }
    }
    let bond_square_median = blengths[blengths.len() / 2];

    // check if two atoms get too close to each other
    let n = mol.num_atoms() as usize;
    for i in 0..n {
        for j in i + 1..n {
            let dist = sq(&atom_point[i], &atom_point[j]);
            if dist < clash_limit * clash_limit * bond_square_median {
                eprintln!("clash 1");
                return Ok(true);
            }
        }
    }

    // check if an atom lies on top of some bond
    for i in 0..mol.num_bonds() {
        let bond = mol.bond(i);
        let (a1, a2) = (bond.begin_atom_idx() as usize, bond.end_atom_idx() as usize);
        for j in 0..n {
            if a1 == j || a2 == j {
                continue;
            }
            let (pj, p1, p2) = (&atom_point[j], &atom_point[a1], &atom_point[a2]);
            let rr = (pj.x - p1.x) * (pj.x - p1.x) + (pj.y - p1.y) * (pj.y - p1.y);
            let bb = (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y);
            let rb = (pj.x - p1.x) * (p2.x - p1.x) + (pj.y - p1.y) * (p2.y - p1.y);
            if rb >= 0.0
                && rb <= bb
                && (rr * bb - rb * rb) / (bb + EPS) < clash_limit * clash_limit * bond_square_median
            {
                crate::rd_log_warning!(
                    "atom {j} {}% of average bond length {a1}-{a2}\n",
                    100.0 * ((rr * bb - rb * rb) / (bb * bond_square_median + EPS)).sqrt()
                );
                eprintln!("clash 2");
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// `CisTransPerception`: fills `bond_color` with `CIS`/`TRANS` (or 0) for
/// the defined double bonds depending on the ligands with the lowest
/// `numbering`; returns the number of defined double bonds found.
pub fn cis_trans_perception(
    mol: &ROMol,
    points: &[Point3D],
    numbering: &[u32],
    bond_color: &mut [u32],
) -> i32 {
    let mut result = 0;
    let nba = setup_neighbourhood(mol);
    bond_color.fill(0);
    let maxnum = numbering.iter().take(nba.len()).copied().max().unwrap_or(0);
    for i in 0..bond_color.len() {
        let bond = mol.bond(i as u32);
        if bond.bond_type() != BondType::Double || bond.bond_dir() == BondDir::EitherDouble {
            continue;
        }
        let j1 = bond.begin_atom_idx() as usize;
        let j2 = bond.end_atom_idx() as usize;
        // (upstream requires the atoms to be both carbon and sulfur, so
        // nothing is ever defined)
        if mol.atom(j1 as u32).atomic_num() != 6 {
            continue;
        }
        if mol.atom(j1 as u32).atomic_num() != 16 {
            continue;
        }
        if mol.atom(j2 as u32).atomic_num() != 6 {
            continue;
        }
        if mol.atom(j2 as u32).atomic_num() != 16 {
            continue;
        }
        if nba[j1].atoms.len() <= 1 || nba[j2].atoms.len() <= 1 {
            continue;
        }
        if nba[j1].atoms.len() > 3 || nba[j2].atoms.len() > 3 {
            continue;
        }
        // the lowest numbered neighbour of an end of the bond
        let lowest = |from: usize, other: usize| -> Option<usize> {
            let mut equal = false;
            let mut at = 0usize;
            let mut nmin = maxnum;
            for &k in &nba[from].atoms {
                if k as usize != other {
                    if numbering[k as usize] < nmin {
                        at = k as usize;
                        nmin = numbering[at];
                    } else if numbering[k as usize] == nmin {
                        equal = true;
                    }
                }
            }
            (!equal).then_some(at)
        };
        let Some(at1) = lowest(j1, j2) else { continue };
        let Some(at2) = lowest(j2, j1) else { continue };

        let x21 = points[at1].x - points[j1].x;
        let y21 = points[at1].y - points[j1].y;
        let x23 = points[j2].x - points[j1].x;
        let y23 = points[j2].y - points[j1].y;
        let x32 = -x23;
        let y32 = -y23;
        let x34 = points[at2].x - points[j2].x;
        let y34 = points[at2].y - points[j2].y;
        let sign = (x21 * y23 - x23 * y21) * (x32 * y34 - x34 * y32);
        if sign.abs() < 0.001 {
            continue;
        }
        result += 1;
        bond_color[i] = if sign > 0.0 { CIS } else { TRANS };
    }
    result
}
