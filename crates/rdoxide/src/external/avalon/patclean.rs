//! Port of `patclean.c`: cleaning, rotating and stereo enforcing of molecules with substructure templates
//! (`ForceStereoTemplate`, `TemplateClean`, `TemplateRotate`). Copyright 2001-2011 Novartis Pharma AG.
//!
//! The templates are modified by the substructure matcher (aromaticity and topography flags) as in C.

use std::io::Write;

use super::geometry::{USE_INWARDS, next_substituent_point, point_set_match_transformation};
use super::graph::{combine_rings, ring_list};
use super::layout::{RUBBER_BOND, flood_color};
use super::perceive::{
    CHAIN, DOUBLE, Neighbourhood, RING, SINGLE, perceive_ring_bonds, setup_neighbourhood,
};
use super::reaccs::{Atom, Bond, CIS_TRANS_EITHER, Molecule, NONE};
use super::ssmatch::ss_match;
use super::stereo::{UNDEFINED_PARITY, add_msg, atom_parity, cis_trans_perception};
use super::utilities::{
    DOWN, UP, compute_implicit_h, flip_stereo_symbols, flood_fill_fragments, reset_colors,
};

/// `MATCH_COLOR` and `BORDER_COLOR`
const MATCH_COLOR: i32 = NONE + 1;
const BORDER_COLOR: i32 = NONE + 2;
/// `ZERO_COUNT`
const ZERO_COUNT: i32 = 1;

/// `MoleculeRings`: the basis rings, rubber bonds are ignored
fn molecule_rings(mp: &Molecule) -> Vec<super::graph::BondSet> {
    if mp.bonds.is_empty() {
        return Vec::new();
    }
    let mut idummy = mp.atoms.len() as u32;
    let bonds: Vec<[u32; 2]> = mp
        .bonds
        .iter()
        .map(|b| {
            if b.bond_type & RUBBER_BOND == 0 {
                [b.atoms[0] as u32, b.atoms[1] as u32]
            } else {
                idummy += 2;
                [idummy - 1, idummy]
            }
        })
        .collect();
    let mut list = ring_list(&bonds);
    combine_rings(&mut list);
    list
}

/// `OverlapStrain`
fn overlap_strain(kernel: &Molecule, new_coords: &[[f64; 2]]) -> f64 {
    let mut result = 0.0;
    for ap in &kernel.atoms {
        for c in new_coords {
            let (dx, dy) = (f64::from(ap.x) - c[0], f64::from(ap.y) - c[1]);
            result += 1.0 / (0.01 + dx * dx + dy * dy);
        }
    }
    result
}

/// `ColorSubstituents`: colours the substituents of the coloured part, one colour per connected
/// substituent, starting after `border_color`; returns a number >= the largest colour used
fn color_substituents(mp: &mut Molecule, border_color: i32) -> i32 {
    let mut max_color = border_color + 1;
    for a in &mut mp.atoms {
        if a.color == NONE {
            a.color = max_color;
            max_color += 1;
        }
    }
    let mut real: Vec<i32> = (0..max_color).collect();
    let merge = |real: &mut Vec<i32>, c1: i32, c2: i32| {
        if real[c1 as usize] < real[c2 as usize] {
            real[c2 as usize] = real[c1 as usize];
        } else {
            real[c1 as usize] = real[c2 as usize];
        }
    };
    // bonds joining at one kernel atom are connected
    for i in 0..mp.atoms.len() {
        let c = mp.atoms[i].color;
        if c == NONE || c > border_color {
            continue;
        }
        let mut col1 = NONE;
        for bp in &mp.bonds {
            if bp.color != border_color
                || (bp.atoms[0] != i as i32 + 1 && bp.atoms[1] != i as i32 + 1)
            {
                continue;
            }
            if bp.atoms[0] == i as i32 + 1 {
                if col1 == NONE {
                    col1 = mp.atoms[(bp.atoms[1] - 1) as usize].color;
                } else {
                    let col2 = mp.atoms[(bp.atoms[1] - 1) as usize].color;
                    merge(&mut real, col1, col2);
                }
            }
            if bp.atoms[1] == i as i32 + 1 {
                if col1 == NONE {
                    col1 = mp.atoms[(bp.atoms[0] - 1) as usize].color;
                } else {
                    let col2 = mp.atoms[(bp.atoms[0] - 1) as usize].color;
                    merge(&mut real, col1, col2);
                }
            }
        }
    }
    loop {
        let mut changed = false;
        for bp in &mp.bonds {
            if bp.color != NONE {
                continue;
            }
            let col1 = mp.atoms[(bp.atoms[0] - 1) as usize].color;
            let col2 = mp.atoms[(bp.atoms[1] - 1) as usize].color;
            if real[col1 as usize] == real[col2 as usize] {
                continue;
            }
            merge(&mut real, col1, col2);
            changed = true;
        }
        if !changed {
            break;
        }
    }
    for i in 0..mp.bonds.len() {
        if mp.bonds[i].color == NONE {
            let c = mp.atoms[(mp.bonds[i].atoms[0] - 1) as usize].color;
            mp.bonds[i].color = real[c as usize];
        }
    }
    for a in &mut mp.atoms {
        if a.color >= border_color {
            a.color = real[a.color as usize];
        }
    }
    max_color
}

fn new_points_of(mp: &Molecule, icol: i32) -> Vec<[f64; 2]> {
    mp.atoms
        .iter()
        .filter(|a| a.color == icol)
        .map(|a| [f64::from(a.x), f64::from(a.y)])
        .collect()
}

/// `TransformColoredSubstituent`: moves the substituent coloured `icol` such that it fits the template
fn transform_colored_substituent(
    mp: &mut Molecule,
    ssp: &Molecule,
    matched: &[usize],
    icol: i32,
) -> bool {
    let (mut tpoints, mut fpoints): (Vec<[f64; 2]>, Vec<[f64; 2]>) = (Vec::new(), Vec::new());
    let mut snumber: Vec<i32> = Vec::new();
    let mut sub_atoms: Vec<usize> = Vec::new();
    let mut sat = 0;
    for bp in &mp.bonds {
        if bp.color != BORDER_COLOR {
            continue;
        }
        // at1: substituent atom, at2: kernel atom, sat: template atom of at2
        let (at1, at2) = if mp.atoms[(bp.atoms[0] - 1) as usize].color == icol {
            ((bp.atoms[0] - 1) as usize, (bp.atoms[1] - 1) as usize)
        } else if mp.atoms[(bp.atoms[1] - 1) as usize].color == icol {
            ((bp.atoms[1] - 1) as usize, (bp.atoms[0] - 1) as usize)
        } else {
            continue;
        };
        sat = match matched.iter().position(|&m| m == at2) {
            Some(s) => s,
            None => return false, // "internal error" in C
        };
        sub_atoms.push(at1);
        fpoints.push([f64::from(mp.atoms[at2].x), f64::from(mp.atoms[at2].y)]);
        tpoints.push([f64::from(ssp.atoms[sat].x), f64::from(ssp.atoms[sat].y)]);
        snumber.push(sat as i32);
    }
    let npoints = fpoints.len();
    if npoints == 0 {
        return false;
    } else if npoints == 1 {
        let at1 = sub_atoms[0];
        let edges: Vec<[u32; 2]> = ssp
            .bonds
            .iter()
            .map(|b| [(b.atoms[0] - 1) as u32, (b.atoms[1] - 1) as u32])
            .collect();
        let coordinates: Vec<[f64; 2]> = ssp
            .atoms
            .iter()
            .map(|a| [f64::from(a.x), f64::from(a.y)])
            .collect();
        let mut point = [0.0, 0.0];
        next_substituent_point(&mut point, &coordinates, &edges, sat, USE_INWARDS);
        tpoints.push(point);
        fpoints.push([f64::from(mp.atoms[at1].x), f64::from(mp.atoms[at1].y)]);
        snumber.push(-1);
    } else if snumber.windows(2).all(|w| w[0] == w[1]) {
        // only one attachment atom: map its environment, too (the spiral case)
        let sat = snumber[0] as usize;
        for bp in &ssp.bonds {
            let at1 = if (bp.atoms[0] - 1) as usize == sat {
                (bp.atoms[1] - 1) as usize
            } else if (bp.atoms[1] - 1) as usize == sat {
                (bp.atoms[0] - 1) as usize
            } else {
                continue;
            };
            let m = &mp.atoms[matched[at1]];
            fpoints.push([f64::from(m.x), f64::from(m.y)]);
            tpoints.push([f64::from(ssp.atoms[at1].x), f64::from(ssp.atoms[at1].y)]);
            snumber.push(at1 as i32);
        }
    }
    let np = fpoints.len();
    let mut new_points = new_points_of(mp, icol);
    point_set_match_transformation(&mut new_points, &fpoints, &tpoints, np, true);
    let strain = overlap_strain(ssp, &new_points);
    let mut new_points = new_points_of(mp, icol);
    point_set_match_transformation(&mut new_points, &fpoints, &tpoints, np, false);
    if strain < overlap_strain(ssp, &new_points) {
        new_points = new_points_of(mp, icol);
        point_set_match_transformation(&mut new_points, &fpoints, &tpoints, np, true);
    }
    let mut k = 0;
    for a in &mut mp.atoms {
        if a.color == icol {
            a.x = new_points[k][0] as f32;
            a.y = new_points[k][1] as f32;
            k += 1;
        }
    }
    true
}

/// `AdjustToSubstructureTemplate`: puts the match on the template coordinates and moves the substituents
fn adjust_to_substructure_template(
    mp: &mut Molecule,
    ssp: &Molecule,
    matched: &[usize],
    log: &mut Option<std::fs::File>,
) -> bool {
    let (na, nb) = (mp.atoms.len(), mp.bonds.len());
    if na == 0 || ssp.atoms.is_empty() {
        return false;
    }
    let saved: Vec<(f32, f32, f32)> = mp.atoms.iter().map(|a| (a.x, a.y, a.z)).collect();
    let Some(nba) = setup_neighbourhood(mp) else {
        return false;
    };
    let mut ring_size = vec![0usize; nb];
    for ring in molecule_rings(mp) {
        for i in 0..nb {
            if ring.bond_set.is_member(i) && (ring_size[i] == 0 || ring_size[i] > ring.cardinality)
            {
                ring_size[i] = ring.cardinality;
            }
        }
    }
    let mut err = "";
    let parity: Vec<i32> = (0..na)
        .map(|i| atom_parity(mp, i + 1, &nba[i], &mut err))
        .collect();
    let mut numbering: Vec<i32> = (0..na as i32).collect();
    let _ = cis_trans_perception(mp, &numbering);
    let cistrans: Vec<i32> = mp.bonds.iter().map(|b| b.color).collect();

    let n_match = matched.len() as f32;
    let (mut oldx, mut oldy) = (0.0f32, 0.0f32);
    for &m in matched {
        oldx += mp.atoms[m].x;
        oldy += mp.atoms[m].y;
    }
    oldx /= n_match;
    oldy /= n_match;

    reset_colors(mp);
    for &at1 in matched {
        mp.atoms[at1].color = MATCH_COLOR;
        for bp in &mut mp.bonds {
            let at2 = if (bp.atoms[0] - 1) as usize == at1 {
                (bp.atoms[1] - 1) as usize
            } else if (bp.atoms[1] - 1) as usize == at1 {
                (bp.atoms[0] - 1) as usize
            } else {
                continue;
            };
            bp.color = if matched.contains(&at2) {
                MATCH_COLOR
            } else {
                BORDER_COLOR
            };
        }
    }
    let max_color = color_substituents(mp, BORDER_COLOR);
    for icol in BORDER_COLOR + 1..=max_color {
        transform_colored_substituent(mp, ssp, matched, icol);
    }
    reset_colors(mp);
    for i in 0..ssp.atoms.len() {
        mp.atoms[matched[i]].x = ssp.atoms[i].x;
        mp.atoms[matched[i]].y = ssp.atoms[i].y;
    }

    let mut stereo_ok = true;
    for i in 0..na {
        let ph = atom_parity(mp, i + 1, &nba[i], &mut err);
        if ph != parity[i] {
            // try all possible combinations of swapping
            for swap_flags in 0..(1 << 4) {
                for j in 0..nba[i].atoms.len() {
                    let bp = &mut mp.bonds[nba[i].bonds[j]];
                    if bp.stereo_symbol == UP || bp.stereo_symbol == DOWN {
                        bp.stereo_symbol = if swap_flags & (1 << j) != 0 { UP } else { DOWN };
                    }
                }
                if atom_parity(mp, i + 1, &nba[i], &mut err) == parity[i] {
                    break;
                }
            }
            if atom_parity(mp, i + 1, &nba[i], &mut err) != parity[i] {
                stereo_ok = false;
            }
        }
    }
    for (i, n) in numbering.iter_mut().enumerate() {
        *n = i as i32;
    }
    let _ = cis_trans_perception(mp, &numbering);
    // double bonds that need to flip are flipped if possible
    for i in 0..nb {
        let (ai1, ai2) = (
            (mp.bonds[i].atoms[0] - 1) as usize,
            (mp.bonds[i].atoms[1] - 1) as usize,
        );
        let cth = mp.bonds[i].color;
        if ring_size[i] == 0 && cistrans[i] != cth {
            mp.atoms[ai1].color = -1;
            mp.atoms[ai2].color = 1;
            flood_color(mp, &nba, ai2, 1);
            mp.atoms[ai1].color = 0;
            let p1 = [f64::from(mp.atoms[ai1].x), f64::from(mp.atoms[ai1].y)];
            let p2 = [f64::from(mp.atoms[ai2].x), f64::from(mp.atoms[ai2].y)];
            let r12 = [p2[0] - p1[0], p2[1] - p1[1]];
            for ap in &mut mp.atoms {
                if ap.color == 1 {
                    let p = [f64::from(ap.x), f64::from(ap.y)];
                    let q = ((p[0] - p1[0]) * r12[0] + (p[1] - p1[1]) * r12[1])
                        / (r12[0] * r12[0] + r12[1] * r12[1]);
                    ap.x = (p1[0] - (p[0] - p1[0]) + 2.0 * r12[0] * q) as f32;
                    ap.y = (p1[1] - (p[1] - p1[1]) + 2.0 * r12[1] * q) as f32;
                }
            }
            flip_stereo_symbols(mp, 1);
            for ap in &mut mp.atoms {
                ap.color = NONE;
            }
        }
    }
    let _ = cis_trans_perception(mp, &numbering);
    for i in 0..nb {
        if cistrans[i] != mp.bonds[i].color {
            stereo_ok = false;
        }
    }
    if !stereo_ok {
        for (a, s) in mp.atoms.iter_mut().zip(&saved) {
            (a.x, a.y, a.z) = *s;
        }
        let msg = format!("template clean failed for molecule '{}'", mp.name);
        match log {
            Some(f) => {
                let _ = writeln!(f, "{msg}");
            }
            None => eprintln!("{msg}"),
        }
        return false;
    }
    let (mut newx, mut newy) = (0.0f32, 0.0f32);
    for &m in matched {
        newx += mp.atoms[m].x;
        newy += mp.atoms[m].y;
    }
    newx /= n_match;
    newy /= n_match;
    for &m in matched {
        mp.atoms[m].color = MATCH_COLOR;
    }
    flood_fill_fragments(mp);
    for ap in &mut mp.atoms {
        if ap.color != NONE {
            ap.x += oldx - newx;
            ap.y += oldy - newy;
            ap.color = NONE;
        }
    }
    reset_colors(mp);
    // C reads `ring_size` after freeing it: the first 16 bytes of the freed block hold the allocator's
    // bookkeeping, so the first four entries are non-zero
    for i in 0..nb {
        let bp = &mp.bonds[i];
        if bp.stereo_symbol != CIS_TRANS_EITHER || bp.bond_type != DOUBLE {
            continue;
        }
        if ring_size[i] == 0 && i >= 4 {
            continue;
        }
        let (at1, at2) = ((bp.atoms[0] - 1) as usize, (bp.atoms[1] - 1) as usize);
        if matched.contains(&at1) && matched.contains(&at2) {
            mp.bonds[i].stereo_symbol = NONE;
        }
    }
    true
}

/// `AlignMoleculeWithSubstructure`: the RMS deviation of the aligned match; `mp` is moved if `transform`
fn align_molecule_with_substructure(
    mp: &mut Molecule,
    ssp: &Molecule,
    matched: &[usize],
    transform: bool,
) -> f64 {
    let mut mpoints: Vec<[f64; 2]> = mp
        .atoms
        .iter()
        .map(|a| [f64::from(a.x), f64::from(a.y)])
        .collect();
    let spoints: Vec<[f64; 2]> = ssp
        .atoms
        .iter()
        .map(|a| [f64::from(a.x), f64::from(a.y)])
        .collect();
    let rpoints: Vec<[f64; 2]> = (0..ssp.atoms.len())
        .map(|i| {
            let a = &mp.atoms[matched[i]];
            [f64::from(a.x), f64::from(a.y)]
        })
        .collect();
    point_set_match_transformation(&mut mpoints, &rpoints, &spoints, ssp.atoms.len(), false);
    let mut fit = 0.0;
    for i in 0..ssp.atoms.len() {
        let diff = mpoints[matched[i]][0] - spoints[i][0];
        fit += diff * diff;
        let diff = mpoints[matched[i]][1] - spoints[i][1];
        fit += diff * diff;
    }
    if transform {
        for (a, p) in mp.atoms.iter_mut().zip(&mpoints) {
            a.x = p[0] as f32;
            a.y = p[1] as f32;
        }
    }
    f64::sqrt(fit)
}

/// `FlipMolecule(mp, NONE)`: mirrors the molecule at the vertical axis through its center (and swaps the wedges)
fn flip_molecule(mp: &mut Molecule) {
    if mp.atoms.is_empty() {
        return;
    }
    let mut xcenter = 0.0f64;
    for a in &mp.atoms {
        xcenter += f64::from(a.x);
    }
    xcenter /= mp.atoms.len() as f64;
    for a in &mut mp.atoms {
        a.x = (xcenter - f64::from(a.x)) as f32;
    }
    flip_stereo_symbols(mp, NONE);
}

/// `ClosestMatch`: the geometrically closest match of `ssp`; `mp` is aligned (and flipped) to it
fn closest_match(mp: &mut Molecule, ssp: &mut Molecule) -> Option<Vec<usize>> {
    let matches = ss_match(mp, ssp, false);
    if matches.is_empty() {
        return None;
    }
    let (mut best_fit, mut best) = (1.0e17, 0);
    let mut best_flip = false;
    for (k, m) in matches.iter().enumerate() {
        let fit = align_molecule_with_substructure(mp, ssp, m, false);
        if fit < best_fit {
            best_fit = fit;
            best = k;
        }
    }
    flip_molecule(mp);
    for (k, m) in matches.iter().enumerate() {
        let fit = align_molecule_with_substructure(mp, ssp, m, false);
        if fit < best_fit {
            best_fit = fit;
            best_flip = true;
            best = k;
        }
    }
    if !best_flip {
        flip_molecule(mp);
    }
    let best_match = matches[best].clone();
    align_molecule_with_substructure(mp, ssp, &best_match, true);
    Some(best_match)
}

/// `AddStereoHydrogen`: adds a hydrogen at atom `iatom` (1 based) with a single bond of `stereo_symbol`
fn add_stereo_hydrogen(mp: &mut Molecule, iatom: usize, stereo_symbol: i32) {
    let coordinates: Vec<[f64; 2]> = mp
        .atoms
        .iter()
        .map(|a| [f64::from(a.x), f64::from(a.y)])
        .collect();
    let edges: Vec<[u32; 2]> = mp
        .bonds
        .iter()
        .map(|b| [(b.atoms[0] - 1) as u32, (b.atoms[1] - 1) as u32])
        .collect();
    let mut point = [0.0, 0.0];
    next_substituent_point(&mut point, &coordinates, &edges, iatom - 1, USE_INWARDS);
    mp.atoms.push(Atom {
        x: point[0] as f32,
        y: point[1] as f32,
        z: 0.0,
        atom_symbol: "H".into(),
        ..Atom::default()
    });
    mp.bonds.push(Bond {
        atoms: [iatom as i32, mp.atoms.len() as i32],
        bond_type: SINGLE,
        stereo_symbol,
        ..Bond::default()
    });
}

/// `ForceStereoTemplate`: 1 if the stereo symbols of `ssp` were applied, 0 if there was no match, -1 if a
/// possible match could not be applied
pub fn force_stereo_template(mp: &mut Molecule, ssp: &mut Molecule, msgs: &mut Vec<String>) -> i32 {
    if mp.atoms.is_empty() || ssp.atoms.is_empty() {
        return 0;
    }
    let Some(best_match) = closest_match(mp, ssp) else {
        return 0;
    };
    let Some(nba) = setup_neighbourhood(mp) else {
        return 0;
    };
    perceive_ring_bonds(mp);
    reset_colors(ssp);
    for bp in &ssp.bonds {
        if bp.stereo_symbol == UP || bp.stereo_symbol == DOWN {
            ssp.atoms[(bp.atoms[0] - 1) as usize].color += 1;
        }
    }
    let (mut problem_found, mut match_applied) = (false, false);
    let mut err = "";
    for (i, &iatom) in best_match.iter().enumerate() {
        let nbp: &Neighbourhood = &nba[iatom];
        if ssp.atoms[i].color == NONE
            || atom_parity(mp, iatom + 1, nbp, &mut err) != UNDEFINED_PARITY
        {
            continue;
        }
        for j in 0..ssp.bonds.len() {
            let bp = &ssp.bonds[j];
            if !(bp.atoms[0] == i as i32 + 1
                && (bp.stereo_symbol == UP || bp.stereo_symbol == DOWN))
            {
                continue;
            }
            let symbol = bp.stereo_symbol;
            let jatom = best_match[(bp.atoms[1] - 1) as usize];
            for k in 0..mp.bonds.len() {
                let mut bm = mp.bonds[k].clone();
                // the C `else` belongs to the inner `if`
                if bm.atoms[1] == iatom as i32 + 1 && bm.atoms[0] == jatom as i32 + 1 {
                    if bm.stereo_symbol == NONE {
                        add_msg(
                            msgs,
                            format!(
                                "ordering of bond bond {}-{} was fixed",
                                iatom + 1,
                                jatom + 1
                            ),
                        );
                        mp.bonds[k].atoms[0] = bm.atoms[1];
                        mp.bonds[k].atoms[1] = bm.atoms[0];
                        bm = mp.bonds[k].clone();
                    } else {
                        problem_found = true;
                        add_msg(
                            msgs,
                            format!("bond {}-{} could not be fixed", iatom + 1, jatom + 1),
                        );
                        break;
                    }
                }
                if bm.atoms[0] == iatom as i32 + 1 && bm.atoms[1] == jatom as i32 + 1 {
                    if bm.stereo_symbol != NONE {
                        problem_found = true;
                        add_msg(
                            msgs,
                            format!("bond {}-{} could not be fixed", iatom + 1, jatom + 1),
                        );
                    } else if bm.topography != RING {
                        add_msg(msgs, format!("bond {}-{} was fixed", iatom + 1, jatom + 1));
                        mp.bonds[k].stereo_symbol = symbol;
                        match_applied = true;
                    } else {
                        add_stereo_hydrogen(mp, iatom + 1, if symbol == UP { DOWN } else { UP });
                        add_msg(
                            msgs,
                            format!(
                                "hydrogen was added for ring bond {}-{}",
                                iatom + 1,
                                jatom + 1
                            ),
                        );
                        match_applied = true;
                    }
                    break;
                }
            }
            break;
        }
    }
    for b in &mut mp.bonds {
        b.dummy = 0;
        b.topography = 0;
    }
    if problem_found {
        -1
    } else if match_applied {
        1
    } else {
        0
    }
}

/// `ScaleByTemplate`: scales `mp` to the average bond length of `ssp` unless the factor exceeds `slimit`
fn scale_by_template(mp: &mut Molecule, ssp: &Molecule, mut slimit: f64) {
    if mp.atoms.is_empty() || ssp.atoms.is_empty() {
        return;
    }
    let avg = |m: &Molecule| {
        let mut bl = 0.0;
        for bp in &m.bonds {
            let (a, b) = (
                &m.atoms[(bp.atoms[0] - 1) as usize],
                &m.atoms[(bp.atoms[1] - 1) as usize],
            );
            // SQR() of a float difference is a float product, and so is their sum
            let (dx, dy) = (a.x - b.x, a.y - b.y);
            bl += f64::from(dx * dx + dy * dy);
        }
        if m.bonds.is_empty() {
            bl
        } else {
            f64::sqrt(bl / m.bonds.len() as f64)
        }
    };
    let mbl = avg(mp);
    let sbl = avg(ssp);
    if slimit > 1.0 {
        slimit = 1.0 / slimit;
    }
    if sbl < mbl && sbl < mbl * slimit {
        return;
    }
    if sbl > mbl && sbl * slimit > mbl {
        return;
    }
    for a in &mut mp.atoms {
        a.x = (f64::from(a.x) * (sbl / mbl)) as f32;
        a.y = (f64::from(a.y) * (sbl / mbl)) as f32;
        a.z = (f64::from(a.z) * (sbl / mbl)) as f32;
    }
}

/// `TemplateClean`: cleans the part of `mp` matching `ssp` to the template coordinates; true if modified
pub fn template_clean(
    mp: &mut Molecule,
    ssp: &mut Molecule,
    log: &mut Option<std::fs::File>,
) -> bool {
    if mp.atoms.is_empty() || ssp.atoms.is_empty() {
        return false;
    }
    let mut h_count = vec![0; mp.atoms.len() + 1];
    compute_implicit_h(mp, &mut h_count);
    for bp in &mp.bonds {
        let (a0, a1) = ((bp.atoms[0] - 1) as usize, (bp.atoms[1] - 1) as usize);
        if mp.atoms[a0].atom_symbol == "H" {
            h_count[bp.atoms[1] as usize] += 1;
        } else if mp.atoms[a1].atom_symbol == "H" {
            h_count[bp.atoms[0] as usize] += 1;
        }
    }
    for i in 0..mp.atoms.len() {
        if h_count[i + 1] >= 0 {
            mp.atoms[i].query_h_count = ZERO_COUNT + h_count[i + 1];
        }
    }
    let best_match = closest_match(mp, ssp);
    for a in &mut mp.atoms {
        a.query_h_count = NONE;
    }
    for b in &mut mp.bonds {
        b.topography = NONE;
    }
    let Some(best_match) = best_match else {
        return false;
    };
    scale_by_template(mp, ssp, 0.1);
    adjust_to_substructure_template(mp, ssp, &best_match, log)
}

/// `TemplateRotate`: rotates `mp` to resemble the template (a clean if all template chain bonds match chain bonds)
pub fn template_rotate(
    mp: &mut Molecule,
    ssp: &mut Molecule,
    log: &mut Option<std::fs::File>,
) -> bool {
    if mp.atoms.is_empty() || ssp.atoms.is_empty() {
        return false;
    }
    let mut chained = ssp.clone();
    perceive_ring_bonds(&mut chained);
    for bp in &mut chained.bonds {
        if bp.topography == RING {
            bp.dummy = 0;
        } else {
            bp.topography = CHAIN;
        }
    }
    if !ss_match(mp, &mut chained, true).is_empty() {
        template_clean(mp, &mut chained, log);
        for b in &mut mp.bonds {
            b.dummy = 0;
            b.topography = 0;
        }
        return true;
    }
    let best_match = closest_match(mp, ssp);
    for b in &mut mp.bonds {
        b.dummy = 0;
        b.topography = 0;
    }
    best_match.is_some()
}
