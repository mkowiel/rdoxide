//! `utilities.c`, `casutils.c`: atom symbol matching, implicit hydrogens, CAS valences.
//! Copyright 2001-2011 Novartis Pharma AG.

use super::perceive::{
    AROMATIC, DOUBLE, Neighbourhood, SINGLE, TRIPLE, ring_state, setup_neighbourhood,
};
use super::reaccs::{DOUBLET, Molecule, NONE};
use super::scan::Scan;
use super::tables::*;

fn in_table(sym: &str, table: &[&str]) -> bool {
    table.contains(&sym)
}

/// `AtomSymbolMatch`: `atsym` is in the comma separated list `pattern`, whose entries may be
/// the pseudo symbols `alk`, `gr2`, `hal`, ... (lower case) or a plain symbol; the whole pattern
/// may also be `A`, `Qs`, `G`, `ONS`, `X`, `Q`, `M` or `Ami`
pub fn atom_symbol_match(atsym: &str, pattern: &str) -> bool {
    for tok in pattern.split(',').filter(|t| !t.is_empty()) {
        if tok.starts_with(|c: char| c.is_ascii_lowercase()) {
            let table: Option<&[&str]> = match tok {
                "alk" => Some(ALKALI_METALS),
                "gr2" => Some(GR2),
                "gr3" => Some(GR3),
                "gr4" => Some(GR4),
                "ons" => Some(ONS_TABLE),
                "on2" => Some(ON2),
                "hal" => Some(HALOGENES),
                "ha2" => Some(HA2),
                "trn" => Some(TRANSITION_METALS),
                "tra" => Some(TRA),
                "trb" => Some(TRB),
                "tm1" => Some(TM1),
                "tm2" => Some(TM2),
                "tm3" => Some(TM3),
                "tm4" => Some(TM4),
                "tm5" => Some(TM5),
                "tm6" => Some(TM6),
                "tm7" => Some(TM7),
                "tm8" => Some(TM8),
                "lan" => Some(LANTHANOIDS),
                "ami" => Some(AMINO_ACIDS),
                _ => None,
            };
            if table.is_some_and(|t| in_table(atsym, t)) {
                return true;
            }
        }
        if atsym == tok {
            return true;
        }
    }
    match pattern {
        "A" => atsym != "H",
        "Qs" => in_table(atsym, NON_METAL_SMALL_SOLUTION),
        "G" => in_table(atsym, HC_TABLE),
        "ONS" => in_table(atsym, ONS_TABLE),
        "X" => in_table(atsym, HALOGENES),
        "Q" => in_table(atsym, NON_METAL_HETERO_ELEMENTS),
        "M" => in_table(atsym, METALS),
        "Ami" => in_table(atsym, AMINO_ACIDS),
        _ => false,
    }
}

/// `CASValence`: 0 for unknown elements
pub fn cas_valence(symbol: &str) -> i32 {
    CAS_VALENCE
        .iter()
        .find(|(s, _)| *s == symbol)
        .map_or(0, |e| e.1)
}

/// `valence_table_entry`: symbol, from, to, step, lone pairs, pair deficit
const VALENCE_TABLE: &[(&str, i32, i32, i32, i32, i32)] = &[
    ("C", 4, 4, 1, 0, 0),
    ("H", -1, 1, 2, 0, 0),
    ("N", 3, 5, 2, 1, 0),
    ("O", 2, 2, 1, 2, 0),
    ("Cl", 1, 7, 2, 0, 0),
    ("P", 3, 5, 2, 1, 0),
    ("S", 2, 6, 2, 2, 0),
    ("F", 1, 1, 1, 1, 0),
    ("H", -1, 1, 2, 0, 0),
    ("Li", -1, 1, 2, 0, 0),
    ("Na", -1, 1, 2, 0, 0),
    ("K", -1, 1, 2, 0, 0),
    ("Rb", -1, 1, 2, 0, 0),
    ("Cs", -1, 1, 2, 0, 0),
    ("Be", -2, 2, 2, 0, 0),
    ("Mg", -2, 2, 2, 0, 0),
    ("B", 3, 3, 1, 0, 1),
    ("Al", -3, 3, 2, 0, 1),
    ("Ga", -3, 3, 2, 0, 1),
    ("In", -3, 3, 2, 0, 1),
    ("Tl", -3, 3, 2, 0, 0),
    ("Si", 4, 4, 1, 0, 0),
    ("As", 3, 5, 2, 0, 0),
    ("Sb", 3, 5, 2, 0, 0),
    ("Bi", 3, 5, 2, 0, 0),
    ("Se", 2, 6, 2, 0, 0),
    ("Te", 2, 6, 2, 0, 0),
    ("La", -3, 3, 2, 0, 1),
    ("Ce", -3, 3, 2, 0, 1),
    ("Pr", -3, 3, 2, 0, 1),
    ("Nd", -3, 3, 2, 0, 1),
    ("Pm", -3, 3, 2, 0, 1),
    ("Sm", -3, 3, 2, 0, 1),
    ("Eu", -3, 3, 2, 0, 1),
    ("Gd", -3, 3, 2, 0, 1),
    ("Tb", -3, 3, 2, 0, 1),
    ("Dy", -3, 3, 2, 0, 1),
    ("Ho", -3, 3, 2, 0, 1),
    ("Er", -3, 3, 2, 0, 1),
    ("Tm", -3, 3, 2, 0, 1),
    ("Yb", -3, 3, 2, 0, 1),
    ("Lu", -3, 3, 2, 0, 1),
    ("Br", 1, 7, 2, 0, 0),
    ("I", 1, 7, 2, 0, 0),
];

/// `ImplicitHydrogens`
pub fn implicit_hydrogens(
    symbol: &str,
    nsingle: i32,
    naromatic: i32,
    ndouble: i32,
    ntriple: i32,
    radical: bool,
    charge: i32,
) -> i32 {
    let mut bond_electrons = nsingle + 2 * ndouble + 3 * ntriple;
    if radical {
        bond_electrons += 1;
    }
    bond_electrons += match naromatic {
        0 => 0,
        1 => 2,
        2 => 3,
        3 => 4,
        n => n + 1,
    };
    for &(name, from, to, step, lone_pairs, pair_deficit) in VALENCE_TABLE {
        if name != symbol {
            continue;
        }
        let mut val = from;
        while val <= to {
            if charge == 0 {
                let h = val - bond_electrons;
                if h >= 0 {
                    return h;
                }
            } else if charge > 0 {
                let h = val - bond_electrons + charge;
                if h >= 0 {
                    return if lone_pairs > 0 { h } else { 0 };
                }
            } else {
                let mut h = val - bond_electrons - charge;
                if h >= 0 {
                    if pair_deficit < h {
                        h = pair_deficit;
                    }
                    return if lone_pairs > 0 { 0 } else { h };
                }
            }
            val += step;
        }
    }
    0
}

/// bond stereo symbols of `mdl.h`
pub const UP: i32 = 0x01;
pub const DOWN: i32 = 0x06;
pub const EITHER: i32 = 0x04;

/// `which_hydrogens` flags of `utilities.h`
pub const NON_STEREO_HYDROGENS: i32 = 1;
pub const ANCILLARY_STEREO: i32 = 2;
pub const ESSENTIAL_STEREO: i32 = 4;
pub const ISOTOPE_HYDROGENS: i32 = 8;
pub const NO_QUERY_H_COUNT: i32 = 16;
pub const ALL_HYDROGENS: i32 =
    NON_STEREO_HYDROGENS | ANCILLARY_STEREO | ESSENTIAL_STEREO | ISOTOPE_HYDROGENS;

/// `ComputeImplicitH`: fills only the zero entries of `h_count` (index origin 1)
pub fn compute_implicit_h(mp: &Molecule, h_count: &mut [i32]) {
    let n = mp.atoms.len();
    let (mut single, mut aromatic, mut double, mut triple) = (
        vec![0; n + 1],
        vec![0; n + 1],
        vec![0; n + 1],
        vec![0; n + 1],
    );
    for b in &mp.bonds {
        let tgt = match b.bond_type {
            DOUBLE => &mut double,
            TRIPLE => &mut triple,
            AROMATIC => &mut aromatic,
            _ => &mut single,
        };
        tgt[b.atoms[0] as usize] += 1;
        tgt[b.atoms[1] as usize] += 1;
    }
    for (i, a) in mp.atoms.iter().enumerate() {
        if h_count[i + 1] == 0 {
            let h = implicit_hydrogens(
                &a.atom_symbol,
                single[i + 1],
                aromatic[i + 1],
                double[i + 1],
                triple[i + 1],
                a.radical == DOUBLET,
                a.charge,
            );
            h_count[i + 1] = h.max(0);
        }
    }
}

/// `ResetColors` (atoms and bonds)
pub fn reset_colors(mp: &mut Molecule) {
    for a in &mut mp.atoms {
        a.color = NONE;
    }
    for b in &mut mp.bonds {
        b.color = NONE;
    }
}

/// `SortBondsWithNumbering`
pub fn sort_bonds_with_numbering(mp: &mut Molecule, numbering: &[i32]) {
    for b in &mut mp.bonds {
        if b.bond_type == SINGLE && b.stereo_symbol != NONE {
            continue;
        }
        if numbering[(b.atoms[0] - 1) as usize] > numbering[(b.atoms[1] - 1) as usize] {
            b.atoms.swap(0, 1);
        }
    }
    for i in 1..mp.bonds.len() {
        for j in (0..i).rev() {
            let key = |k: usize, p: usize| numbering[(mp.bonds[k].atoms[p] - 1) as usize];
            if key(j, 0) > key(j + 1, 0)
                || (key(j, 0) == key(j + 1, 0) && key(j, 1) > key(j + 1, 1))
            {
                mp.bonds.swap(j, j + 1);
            } else {
                break;
            }
        }
    }
}

/// `SortNeighboursByNumbering`
pub fn sort_neighbours_by_numbering(nbp: &mut [Neighbourhood], numbering: &[i32]) {
    for nb in nbp.iter_mut() {
        for j in 1..nb.atoms.len() {
            for jj in (0..j).rev() {
                if numbering[nb.atoms[jj]] > numbering[nb.atoms[jj + 1]] {
                    nb.atoms.swap(jj, jj + 1);
                    nb.bonds.swap(jj, jj + 1);
                } else {
                    break;
                }
            }
        }
    }
}

/// `GetNumProperty`: the value stored for atom `atno` in the `tag` property lines
pub fn get_num_property(prop_lines: &[String], tag: &str, atno: i32) -> Option<i32> {
    for line in prop_lines {
        if !line.starts_with(tag) {
            continue;
        }
        let rest = line.get(tag.len()..).unwrap_or("");
        let nentries = Scan::new(rest.as_bytes()).int(3).unwrap_or(0);
        let mut sc = Scan::new(rest.get(3..).unwrap_or("").as_bytes());
        let mut vals = Vec::new();
        while vals.len() < 16 {
            let Some(v) = sc.int(0) else { break };
            vals.push(v);
        }
        for n in 0..nentries.max(0) as usize {
            if 2 * n + 1 < vals.len() && vals[2 * n] == atno {
                return Some(vals[2 * n + 1]);
            }
        }
    }
    None
}

const TRIVIAL: i32 = 1;
const ANCILLARY: i32 = 2;
const ESSENTIAL: i32 = 3;

/// `PrepareCenterForHDrop`: arranges wedge bonds around `iatom` before its hydrogens are dropped
pub fn prepare_center_for_h_drop(
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    iatom: usize,
    atom_status: &[i32],
) {
    if mp.atoms[iatom].atom_symbol == "H" || nbp[iatom].atoms.len() < 3 {
        return;
    }
    let nb = &nbp[iatom];
    let (mut nup, mut ndown, mut nhydrogens) = (0, 0, 0);
    let (mut last_up, mut last_down, mut last_h) = (usize::MAX, usize::MAX, usize::MAX);
    let mut stereo_symbol = NONE;
    let mut bph: Option<usize> = None;
    for i in 0..nb.atoms.len() {
        let (ai, bi) = (nb.atoms[i], nb.bonds[i]);
        let bp = &mp.bonds[bi];
        if bp.stereo_symbol == UP && bp.atoms[0] == iatom as i32 + 1 {
            nup += 1;
            last_up = ai;
        }
        if bp.stereo_symbol == DOWN && bp.atoms[0] == iatom as i32 + 1 {
            ndown += 1;
            last_down = ai;
        }
        let ap = &mp.atoms[ai];
        if ap.atom_symbol == "H"
            && ap.mass_difference == NONE
            && ap.charge == NONE
            && ap.radical == NONE
        {
            nhydrogens += 1;
            last_h = ai;
            stereo_symbol = bp.stereo_symbol;
            bph = Some(bi);
        }
    }
    let (cx, cy) = (f64::from(mp.atoms[iatom].x), f64::from(mp.atoms[iatom].y));
    let rel =
        |mp: &Molecule, k: usize| (f64::from(mp.atoms[k].x) - cx, f64::from(mp.atoms[k].y) - cy);
    if nhydrogens == 0 {
        return;
    } else if nhydrogens > 0 && nup + ndown == 0 {
        for &ai in &nb.atoms {
            let ap = &mut mp.atoms[ai];
            if ap.atom_symbol == "H" {
                ap.color = if ap.mass_difference == NONE && ap.charge == NONE && ap.radical == NONE
                {
                    TRIVIAL
                } else {
                    ESSENTIAL
                };
            }
        }
        return;
    } else if nhydrogens >= 2 && nup + ndown != 0 {
        for i in 0..nb.atoms.len() {
            let bp = &mut mp.bonds[nb.bonds[i]];
            if bp.atoms[0] == iatom as i32 + 1 {
                bp.stereo_symbol = NONE;
            }
            mp.atoms[nb.atoms[i]].color = TRIVIAL;
        }
        return;
    }
    // exactly one hydrogen from here on
    let min_cosa_others = |mp: &Molecule| -> f64 {
        let mut best = 2.0;
        for i in 0..nb.atoms.len() {
            if nb.atoms[i] == last_h {
                continue;
            }
            let (x1, y1) = rel(mp, nb.atoms[i]);
            for j in i + 1..nb.atoms.len() {
                if nb.atoms[j] == last_h {
                    continue;
                }
                let (xh, yh) = rel(mp, nb.atoms[j]);
                let mut cosa = xh * x1 + yh * y1;
                let norm = (xh * xh + yh * yh).sqrt() * (x1 * x1 + y1 * y1).sqrt();
                if norm > 0.0001 {
                    cosa /= norm;
                } else {
                    continue;
                }
                if cosa < best {
                    best = cosa;
                }
            }
        }
        best
    };
    if stereo_symbol != UP && stereo_symbol != DOWN {
        let best_cosa = min_cosa_others(mp);
        if nup + ndown == 0 || best_cosa > -0.8 {
            mp.atoms[last_h].color = TRIVIAL;
        } else {
            mp.atoms[last_h].color = ANCILLARY;
        }
    } else if nup + ndown == 1 {
        let (xh, yh) = rel(mp, last_h);
        let mut best_index = 0;
        let mut best_cosa = -2.0;
        for i in 0..nb.atoms.len() {
            if nb.atoms[i] == last_h {
                continue;
            }
            let bp = &mp.bonds[nb.bonds[i]];
            if bp.stereo_symbol != NONE || bp.bond_type != SINGLE {
                continue;
            }
            let (x1, y1) = rel(mp, nb.atoms[i]);
            let mut cosa = xh * x1 + yh * y1;
            let norm = (xh * xh + yh * yh).sqrt() * (x1 * x1 + y1 * y1).sqrt();
            if norm > 0.001 {
                cosa /= norm;
            } else {
                continue;
            }
            if mp.atoms[nb.atoms[i]].atom_symbol != "C" {
                cosa *= 1.0001;
            }
            if atom_status[nb.atoms[i]] == 0 {
                cosa *= 1.0002;
            }
            if cosa > best_cosa {
                best_cosa = cosa;
                best_index = i;
            }
        }
        if best_cosa < 0.05 {
            mp.atoms[last_h].color = ESSENTIAL;
        } else {
            let best_cosa = min_cosa_others(mp);
            let bp = &mut mp.bonds[nb.bonds[best_index]];
            bp.stereo_symbol = if stereo_symbol == UP { DOWN } else { UP };
            if bp.atoms[1] == iatom as i32 + 1 {
                bp.atoms.swap(0, 1);
            }
            if let Some(b) = bph {
                mp.bonds[b].stereo_symbol = NONE;
            }
            mp.atoms[last_h].color = if best_cosa > -0.8 { TRIVIAL } else { ANCILLARY };
        }
    } else if nup + ndown == 2 && nup * ndown == 1 {
        let ipartner = if stereo_symbol == UP {
            last_down
        } else {
            last_up
        };
        let (xh, yh) = rel(mp, last_h);
        let (x1, y1) = rel(mp, ipartner);
        let (xi, yi) = ((xh + x1) / 2.0, (yh + y1) / 2.0);
        let mut ref_cosa = xi * x1 + yi * y1;
        let norm = (xi * xi + yi * yi).sqrt() * (x1 * x1 + y1 * y1).sqrt();
        if norm > 0.001 {
            ref_cosa /= norm;
        }
        let mut problem = false;
        for i in 0..nb.atoms.len() {
            if nb.atoms[i] == last_h || nb.atoms[i] == ipartner {
                continue;
            }
            if mp.bonds[nb.bonds[i]].bond_type != SINGLE {
                continue;
            }
            let (x2, y2) = rel(mp, nb.atoms[i]);
            let mut cosa = xi * x2 + yi * y2;
            let norm = (xi * xi + yi * yi).sqrt() * (x2 * x2 + y2 * y2).sqrt();
            if norm > 0.001 {
                cosa /= norm;
            }
            if cosa > ref_cosa {
                problem = true;
            }
        }
        if problem {
            mp.atoms[last_h].color = ESSENTIAL;
        } else {
            mp.atoms[last_h].color = ANCILLARY;
            if let Some(b) = bph {
                mp.bonds[b].stereo_symbol = NONE;
            }
        }
    } else {
        mp.atoms[last_h].color = ESSENTIAL;
    }
}

/// `MakeSomeHydrogensImplicit`
pub fn make_some_hydrogens_implicit(mp: &mut Molecule, which_hydrogens: i32) {
    for a in &mut mp.atoms {
        a.color = NONE;
    }
    let nbp = setup_neighbourhood(mp).unwrap_or_default();
    let (atom_status, _) = ring_state(mp);
    if nbp.len() == mp.atoms.len() {
        for i in 0..mp.atoms.len() {
            prepare_center_for_h_drop(mp, &nbp, i, &atom_status);
        }
    }
    let mut hnr = 1;
    while hnr <= mp.atoms.len() {
        if mp.atoms[hnr - 1].atom_symbol != "H" {
            hnr += 1;
            continue;
        }
        if which_hydrogens & ISOTOPE_HYDROGENS == 0 && mp.atoms[hnr - 1].mass_difference != 0 {
            hnr += 1;
            continue;
        }
        let (mut nlig, mut partner, mut bond_index) = (0, 0, 0);
        for (i, b) in mp.bonds.iter().enumerate() {
            if b.atoms[0] == hnr as i32 {
                partner = b.atoms[1];
                bond_index = i;
                nlig += 1;
            }
            if b.atoms[1] == hnr as i32 {
                partner = b.atoms[0];
                bond_index = i;
                nlig += 1;
            }
        }
        if nlig != 1 {
            hnr += 1;
            continue;
        }
        if mp.atoms[(partner - 1) as usize]
            .atom_symbol
            .starts_with('R')
        {
            hnr += 1;
            continue;
        }
        let ss = mp.bonds[bond_index].stereo_symbol;
        if which_hydrogens & (ANCILLARY_STEREO | ESSENTIAL_STEREO) == 0
            && ss != EITHER
            && ss != NONE
        {
            hnr += 1;
            continue;
        }
        let color = mp.atoms[hnr - 1].color;
        if which_hydrogens & ESSENTIAL_STEREO == 0 && color == ESSENTIAL {
            hnr += 1;
            continue;
        }
        if which_hydrogens & ANCILLARY_STEREO == 0 && color == ANCILLARY {
            hnr += 1;
            continue;
        }
        if which_hydrogens & NO_QUERY_H_COUNT == 0 {
            let p = &mut mp.atoms[(partner - 1) as usize];
            if p.query_h_count == 0 {
                p.query_h_count = 2;
            } else {
                p.query_h_count += 1;
            }
        }
        mp.atoms.remove(hnr - 1);
        mp.bonds.remove(bond_index);
        for b in &mut mp.bonds {
            for k in 0..2 {
                if b.atoms[k] > hnr as i32 {
                    b.atoms[k] -= 1;
                }
            }
        }
        for sl in &mut mp.symbol_lists {
            if sl.atom > hnr as i32 {
                sl.atom -= 1;
            }
        }
        for line in &mut mp.prop_lines {
            if line.starts_with("M  RGP") {
                let nrgp = Scan::new(line.get(6..).unwrap_or("").as_bytes())
                    .int(0)
                    .unwrap_or(0);
                let mut pos = 9;
                for _ in 0..nrgp.max(0) {
                    let mut sc = Scan::new(line.get(pos..).unwrap_or("").as_bytes());
                    if let (Some(mut ratom), Some(rval)) = (sc.int(0), sc.int(0)) {
                        if ratom > hnr as i32 {
                            ratom -= 1;
                            let end = (pos + 8).min(line.len());
                            line.replace_range(pos..end, &format!(" {ratom:3} {rval:3}"));
                        }
                    }
                    pos += 8;
                }
            }
        }
    }
    for a in &mut mp.atoms {
        a.color = NONE;
    }
}

/// `MakeHydrogensImplicit`
pub fn make_hydrogens_implicit(mp: &mut Molecule) {
    make_some_hydrogens_implicit(mp, ALL_HYDROGENS);
}

/// `SortNeighbourhood`: sorts the ligands by descending bond type, then descending atom
/// symbol, then descending query H count
pub fn sort_neighbourhood(np: &mut Neighbourhood, mp: &Molecule) {
    for i in 1..np.atoms.len() {
        for j in (0..i).rev() {
            let (btj, btj1) = (
                mp.bonds[np.bonds[j]].bond_type,
                mp.bonds[np.bonds[j + 1]].bond_type,
            );
            let (aj, aj1) = (&mp.atoms[np.atoms[j]], &mp.atoms[np.atoms[j + 1]]);
            let ord = aj.atom_symbol.as_bytes().cmp(aj1.atom_symbol.as_bytes());
            if btj < btj1
                || (btj == btj1 && ord == std::cmp::Ordering::Greater)
                || (btj == btj1
                    && ord == std::cmp::Ordering::Equal
                    && aj.query_h_count < aj1.query_h_count)
            {
                np.bonds.swap(j, j + 1);
                np.atoms.swap(j, j + 1);
            } else {
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symbol_groups() {
        assert!(atom_symbol_match("Cl", "hal"));
        assert!(atom_symbol_match("N", "O,S,P,N,L"));
        assert!(atom_symbol_match("C", "A"));
        assert!(!atom_symbol_match("H", "A"));
        assert!(atom_symbol_match("O", "Q"));
        assert_eq!(cas_valence("C"), 4);
        assert_eq!(cas_valence("??x"), 0);
        assert_eq!(implicit_hydrogens("N", 3, 0, 0, 0, false, 0), 0);
        assert_eq!(implicit_hydrogens("O", 1, 0, 0, 0, false, 0), 1);
    }
}
