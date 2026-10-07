//! Ported from RDKit `Code/GraphMol/NontetrahedralStereo.cpp`.
//!
//! Copyright (C) 2022 Greg Landrum and NextMove Software. Licensed under the
//! BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! Permutation tables were transcribed mechanically from upstream.

use crate::graphmol::atom::ChiralType;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::RDValue;
use crate::rdgeneral::common_properties as cp;

const SWAP_SQUAREPLANAR_TABLE: [[u8; 6]; 4] = [
    [0, 0, 0, 0, 0, 0],
    [3, 1, 2, 2, 1, 3],
    [2, 3, 1, 1, 3, 2],
    [1, 2, 3, 3, 2, 1],
];

const SWAP_TRIGONALBIPYRAMIDAL_TABLE: [[u8; 10]; 21] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [9, 20, 17, 2, 2, 2, 7, 2, 6, 3],
    [11, 15, 18, 1, 1, 1, 8, 1, 5, 4],
    [10, 19, 4, 18, 4, 8, 4, 5, 4, 1],
    [12, 16, 3, 17, 3, 7, 3, 6, 3, 2],
    [13, 6, 16, 20, 7, 6, 6, 3, 2, 6],
    [14, 5, 19, 15, 8, 5, 5, 4, 1, 5],
    [8, 14, 10, 11, 5, 4, 1, 8, 8, 8],
    [7, 13, 12, 9, 6, 3, 2, 7, 7, 7],
    [1, 11, 11, 8, 15, 18, 11, 11, 14, 10],
    [3, 12, 7, 12, 16, 12, 17, 13, 12, 9],
    [2, 9, 9, 7, 20, 17, 9, 9, 13, 12],
    [4, 10, 8, 10, 19, 10, 18, 14, 10, 11],
    [5, 8, 14, 14, 14, 19, 15, 10, 11, 14],
    [6, 7, 13, 13, 13, 16, 20, 12, 9, 13],
    [20, 2, 20, 6, 9, 20, 13, 17, 20, 16],
    [19, 4, 5, 19, 10, 14, 19, 19, 18, 15],
    [18, 18, 1, 4, 18, 11, 10, 15, 19, 18],
    [17, 17, 2, 3, 17, 9, 12, 20, 16, 17],
    [16, 3, 6, 16, 12, 13, 16, 16, 17, 20],
    [15, 1, 15, 5, 11, 15, 14, 18, 15, 19],
];

const SWAP_OCTAHEDRAL_TABLE: [[u8; 15]; 31] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [17, 16, 30, 21, 2, 14, 2, 10, 25, 8, 2, 22, 4, 7, 3],
    [7, 3, 25, 22, 1, 4, 1, 8, 30, 10, 1, 21, 14, 17, 16],
    [18, 2, 29, 16, 22, 15, 16, 26, 11, 9, 21, 16, 6, 5, 1],
    [15, 18, 19, 28, 14, 2, 8, 14, 27, 14, 10, 24, 1, 6, 5],
    [14, 17, 20, 15, 27, 16, 9, 28, 15, 15, 23, 11, 7, 3, 4],
    [16, 14, 18, 26, 24, 17, 29, 18, 13, 19, 12, 18, 3, 4, 7],
    [2, 15, 17, 23, 25, 18, 30, 12, 17, 20, 17, 13, 5, 1, 6],
    [23, 26, 11, 12, 10, 10, 4, 2, 29, 1, 14, 20, 10, 13, 9],
    [24, 25, 10, 11, 13, 11, 5, 30, 16, 3, 19, 15, 12, 11, 8],
    [20, 29, 9, 13, 8, 8, 14, 1, 26, 2, 4, 23, 8, 12, 11],
    [19, 30, 8, 9, 12, 9, 15, 25, 3, 16, 24, 5, 13, 9, 10],
    [22, 27, 13, 8, 11, 13, 28, 7, 18, 21, 6, 17, 9, 10, 13],
    [21, 28, 12, 10, 9, 12, 27, 17, 6, 22, 18, 7, 11, 8, 12],
    [5, 6, 24, 27, 4, 1, 10, 4, 28, 4, 8, 19, 2, 18, 15],
    [4, 7, 23, 5, 28, 3, 11, 27, 5, 5, 20, 9, 17, 16, 14],
    [6, 1, 26, 3, 21, 5, 3, 29, 9, 11, 22, 3, 18, 15, 2],
    [1, 5, 7, 20, 30, 6, 25, 13, 7, 23, 7, 12, 15, 2, 18],
    [3, 4, 6, 29, 19, 7, 26, 6, 12, 24, 13, 6, 16, 14, 17],
    [11, 24, 4, 30, 18, 25, 23, 24, 22, 6, 9, 14, 21, 24, 20],
    [10, 23, 5, 17, 29, 26, 24, 21, 23, 7, 15, 8, 23, 22, 19],
    [13, 22, 28, 1, 16, 27, 22, 20, 24, 12, 3, 2, 19, 23, 22],
    [12, 21, 27, 2, 3, 28, 21, 23, 19, 13, 16, 1, 24, 20, 21],
    [8, 20, 15, 7, 26, 29, 19, 22, 20, 17, 5, 10, 20, 21, 24],
    [9, 19, 14, 25, 6, 30, 20, 19, 21, 18, 11, 4, 22, 19, 23],
    [30, 9, 2, 24, 7, 19, 17, 11, 1, 29, 30, 28, 27, 30, 26],
    [29, 8, 16, 6, 23, 20, 18, 3, 10, 30, 27, 29, 29, 28, 25],
    [28, 12, 22, 14, 5, 21, 13, 15, 4, 28, 26, 30, 25, 29, 28],
    [27, 13, 21, 4, 15, 22, 12, 5, 14, 27, 29, 25, 30, 26, 27],
    [26, 10, 3, 18, 20, 23, 6, 16, 8, 25, 28, 26, 26, 27, 30],
    [25, 11, 1, 19, 17, 24, 7, 9, 2, 26, 25, 27, 28, 25, 29],
];

const SQUAREPLANAR_ACROSS: [[u8; 4]; 4] = [[4, 4, 4, 4], [2, 3, 0, 1], [1, 0, 3, 2], [3, 2, 1, 0]];

const TRIGONALBIPYRAMIDAL_ACROSS: [[u8; 5]; 21] = [
    [5, 5, 5, 5, 5],
    [4, 5, 5, 5, 0],
    [4, 5, 5, 5, 0],
    [3, 5, 5, 0, 5],
    [3, 5, 5, 0, 5],
    [2, 5, 0, 5, 5],
    [2, 5, 0, 5, 5],
    [1, 0, 5, 5, 5],
    [1, 0, 5, 5, 5],
    [5, 4, 5, 5, 1],
    [5, 3, 5, 1, 5],
    [5, 4, 5, 5, 1],
    [5, 3, 5, 1, 5],
    [5, 2, 1, 5, 5],
    [5, 2, 1, 5, 5],
    [5, 5, 4, 5, 2],
    [5, 5, 3, 2, 5],
    [5, 5, 5, 4, 3],
    [5, 5, 5, 4, 3],
    [5, 5, 3, 2, 5],
    [5, 5, 4, 5, 2],
];

const OCTAHEDRAL_ACROSS: [[u8; 6]; 31] = [
    [6, 6, 6, 6, 6, 6],
    [5, 3, 4, 1, 2, 0],
    [5, 3, 4, 1, 2, 0],
    [4, 3, 5, 1, 0, 2],
    [5, 4, 3, 2, 1, 0],
    [4, 5, 3, 2, 0, 1],
    [3, 4, 5, 0, 1, 2],
    [3, 5, 4, 0, 2, 1],
    [5, 2, 1, 4, 3, 0],
    [4, 2, 1, 5, 0, 3],
    [5, 2, 1, 4, 3, 0],
    [4, 2, 1, 5, 0, 3],
    [3, 2, 1, 0, 5, 4],
    [3, 2, 1, 0, 5, 4],
    [5, 4, 3, 2, 1, 0],
    [4, 5, 3, 2, 0, 1],
    [4, 3, 5, 1, 0, 2],
    [3, 5, 4, 0, 2, 1],
    [3, 4, 5, 0, 1, 2],
    [2, 4, 0, 5, 1, 3],
    [2, 5, 0, 4, 3, 1],
    [2, 3, 0, 1, 5, 4],
    [2, 3, 0, 1, 5, 4],
    [2, 5, 0, 4, 3, 1],
    [2, 4, 0, 5, 1, 3],
    [1, 0, 4, 5, 2, 3],
    [1, 0, 5, 4, 3, 2],
    [1, 0, 3, 2, 5, 4],
    [1, 0, 3, 2, 5, 4],
    [1, 0, 5, 4, 3, 2],
    [1, 0, 4, 5, 2, 3],
];

const TRIGONALBIPYRAMIDAL_AXIAL: [[u8; 2]; 21] = [
    [5, 5],
    [0, 4],
    [0, 4],
    [0, 3],
    [0, 3],
    [0, 2],
    [0, 2],
    [0, 1],
    [0, 1],
    [1, 4],
    [1, 4],
    [1, 3],
    [1, 3],
    [1, 2],
    [1, 2],
    [2, 4],
    [2, 3],
    [3, 4],
    [3, 4],
    [2, 3],
    [2, 4],
];

fn swap_generic(perm: u32, x: u32, y: u32, offsets: &[u32], max: u32, table: &[&[u8]]) -> u32 {
    if x == y {
        return perm;
    }
    let (lo, hi) = if x < y { (x, y) } else { (y, x) };
    if hi > max {
        return 0;
    }
    let swapidx = offsets[lo as usize] + (hi - 1);
    table
        .get(perm as usize)
        .map_or(0, |row| u32::from(row[swapidx as usize]))
}

fn swap_squareplanar(perm: u32, x: u32, y: u32) -> u32 {
    let t: Vec<&[u8]> = SWAP_SQUAREPLANAR_TABLE.iter().map(|r| &r[..]).collect();
    swap_generic(perm, x, y, &[0, 2, 3], 3, &t)
}

fn swap_trigonalbipyramidal(perm: u32, x: u32, y: u32) -> u32 {
    let t: Vec<&[u8]> = SWAP_TRIGONALBIPYRAMIDAL_TABLE
        .iter()
        .map(|r| &r[..])
        .collect();
    swap_generic(perm, x, y, &[0, 3, 5, 6], 4, &t)
}

fn swap_octahedral(perm: u32, x: u32, y: u32) -> u32 {
    let t: Vec<&[u8]> = SWAP_OCTAHEDRAL_TABLE.iter().map(|r| &r[..]).collect();
    swap_generic(perm, x, y, &[0, 4, 7, 9, 10], 5, &t)
}

/// `getMaxNbors`
pub fn get_max_nbors(tag: ChiralType) -> u32 {
    match tag {
        ChiralType::TetrahedralCW | ChiralType::TetrahedralCCW | ChiralType::Tetrahedral => 4,
        ChiralType::Allene => 2,
        ChiralType::SquarePlanar => 4,
        ChiralType::TrigonalBipyramidal => 5,
        ChiralType::Octahedral => 6,
        _ => {
            crate::rd_log_warning!(
                "Warning: unexpected chiral tag getMaxNbors(): {}",
                tag as u8
            );
            0
        }
    }
}

fn permutation(mol: &ROMol, cen: u32) -> u32 {
    mol.atom(cen)
        .props
        .get(cp::_chiralPermutation)
        .and_then(RDValue::as_uint)
        .unwrap_or(0)
}

/// `hasNonTetrahedralStereo`
pub fn has_non_tetrahedral_stereo(mol: &ROMol, cen: u32) -> bool {
    matches!(
        mol.atom(cen).chiral_tag(),
        ChiralType::SquarePlanar | ChiralType::TrigonalBipyramidal | ChiralType::Octahedral
    )
}

/// `getChiralAcrossBond(cen, qry)`: index of the bond trans to `qry`.
pub fn get_chiral_across_bond(mol: &ROMol, cen: u32, qry: u32) -> Option<u32> {
    let tag = mol.atom(cen).chiral_tag();
    let perm = permutation(mol, cen) as usize;
    if perm == 0 {
        return None;
    }
    let ref_max = get_max_nbors(tag) as usize;
    let mut refs = Vec::with_capacity(6);
    let mut found: Option<usize> = None;
    for b in mol.atom_bond_indices(cen) {
        if refs.len() == ref_max {
            return None;
        }
        if b == qry {
            found = Some(refs.len());
        }
        refs.push(b);
    }
    let f = found?;
    let across = match tag {
        ChiralType::SquarePlanar => SQUAREPLANAR_ACROSS.get(perm).map_or(4, |r| r[f] as usize),
        ChiralType::TrigonalBipyramidal => TRIGONALBIPYRAMIDAL_ACROSS
            .get(perm)
            .map_or(5, |r| r[f] as usize),
        ChiralType::Octahedral => OCTAHEDRAL_ACROSS.get(perm).map_or(6, |r| r[f] as usize),
        _ => return None,
    };
    refs.get(across).copied()
}

/// `getChiralAcrossAtom(cen, qryAtom)`
pub fn get_chiral_across_atom(mol: &ROMol, cen: u32, qry_atom: u32) -> Option<u32> {
    let b = mol.bond_idx_between(cen, qry_atom)?;
    let across = get_chiral_across_bond(mol, cen, b)?;
    mol.bond(across).other_atom_idx(cen).ok()
}

/// `isTrigonalBipyramidalAxialBond`: 1 / -1 for the two axial bonds, else 0.
pub fn is_trigonal_bipyramidal_axial_bond(mol: &ROMol, cen: u32, qry: u32) -> i32 {
    let atom = mol.atom(cen);
    if mol.atom_degree(cen) > 5 || atom.chiral_tag() != ChiralType::TrigonalBipyramidal {
        return 0;
    }
    let perm = permutation(mol, cen) as usize;
    if perm == 0 || perm > 20 {
        return 0;
    }
    for (count, b) in mol.atom_bond_indices(cen).enumerate() {
        if b == qry {
            if count == TRIGONALBIPYRAMIDAL_AXIAL[perm][0] as usize {
                return 1;
            }
            if count == TRIGONALBIPYRAMIDAL_AXIAL[perm][1] as usize {
                return -1;
            }
            return 0;
        }
    }
    0
}

/// `getTrigonalBipyramidalAxialBond(cen, axial)`
pub fn get_trigonal_bipyramidal_axial_bond(mol: &ROMol, cen: u32, axial: i32) -> Option<u32> {
    if mol.atom(cen).chiral_tag() != ChiralType::TrigonalBipyramidal || mol.atom_degree(cen) > 5 {
        return None;
    }
    let perm = permutation(mol, cen) as usize;
    if perm == 0 || perm > 20 {
        return None;
    }
    let idx = if axial != -1 {
        TRIGONALBIPYRAMIDAL_AXIAL[perm][0]
    } else {
        TRIGONALBIPYRAMIDAL_AXIAL[perm][1]
    };
    mol.atom_bond_indices(cen).nth(idx as usize)
}

/// `getIdealAngleBetweenLigands`
pub fn get_ideal_angle_between_ligands(mol: &ROMol, cen: u32, lig1: u32, lig2: u32) -> f64 {
    match mol.atom(cen).chiral_tag() {
        ChiralType::SquarePlanar | ChiralType::Octahedral => {
            if get_chiral_across_atom(mol, cen, lig1) == Some(lig2) {
                180.0
            } else {
                90.0
            }
        }
        ChiralType::TrigonalBipyramidal => {
            let axial = |l: u32| {
                mol.bond_idx_between(cen, l)
                    .map_or(0, |b| is_trigonal_bipyramidal_axial_bond(mol, cen, b))
            };
            if get_chiral_across_atom(mol, cen, lig1) == Some(lig2) {
                180.0
            } else if axial(lig1) != 0 || axial(lig2) != 0 {
                90.0
            } else {
                120.0
            }
        }
        _ => 0.0,
    }
}

/// `getChiralPermutation(cen, probe, inverse)`; a probe entry of -1 marks an
/// implicit/missing ligand.
pub fn get_chiral_permutation(mol: &ROMol, cen: u32, probe: &[i32], inverse: bool) -> u32 {
    let perm0 = mol
        .atom(cen)
        .props
        .get(cp::_chiralPermutation)
        .and_then(RDValue::as_int)
        .unwrap_or(-1);
    if perm0 <= 0 {
        return 0;
    }
    let swap_func: fn(u32, u32, u32) -> u32 = match mol.atom(cen).chiral_tag() {
        ChiralType::Octahedral if probe.len() <= 6 => swap_octahedral,
        ChiralType::TrigonalBipyramidal if probe.len() <= 5 => swap_trigonalbipyramidal,
        ChiralType::SquarePlanar if probe.len() <= 4 => swap_squareplanar,
        _ => return 0,
    };
    let mut perm = perm0 as u32;
    let mut order = vec![u32::MAX; mol.num_bonds() as usize];
    let mut nbr_idx = 0u32;
    for b in mol.atom_bond_indices(cen) {
        order[b as usize] = nbr_idx;
        nbr_idx += 1;
    }
    let mut nbr_perm: Vec<u32> = (0..nbr_idx).collect();
    let mut probe_perm: Vec<u32> = probe
        .iter()
        .map(|&v| if v < 0 { u32::MAX } else { order[v as usize] })
        .collect();
    if nbr_perm.len() < probe_perm.len() {
        nbr_perm.resize(probe_perm.len(), u32::MAX);
    }
    if nbr_perm.len() != probe_perm.len() {
        return 0;
    }
    if inverse {
        std::mem::swap(&mut nbr_perm, &mut probe_perm);
    }
    for i in 0..probe_perm.len().saturating_sub(1) {
        let pval = probe_perm[i];
        if nbr_perm[i] == pval {
            continue;
        }
        let Some(tgt) = (i..nbr_perm.len()).find(|&j| nbr_perm[j] == pval) else {
            return 0;
        };
        perm = swap_func(perm, i as u32, tgt as u32);
        nbr_perm.swap(tgt, i);
    }
    perm
}
