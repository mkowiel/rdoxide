//! Port of Avalon Toolkit's `layout.c` (2D coordinate generation), in slices. This slice
//! holds the colouring helpers and the ring system layout. Where C mixes `float` fields with
//! `double` arithmetic the same single precision steps are taken here.

#![allow(clippy::approx_constant)]

use super::geometry::transform_points;
use super::graph::{BondSet, CRand, combine_rings, ring_list};
use super::perceive::{Neighbourhood, ring_state};
use super::reaccs::{Molecule, NONE};
use super::utilities::{atom_symbol_match, flip_stereo_symbols};

pub const RUBBER_BOND: i32 = 0x40;
pub const DONT_FLIP_BOND: i32 = 0x20;
pub const KEEP_POSITION: i32 = 0x4000;

const STDBOND: f64 = 1.514;
const PI: f64 = 3.14159265359;

/// `RandomCoordinates`: random coordinates for the atoms without `KEEP_POSITION`
pub fn random_coordinates(mp: &mut Molecule, rng: &mut CRand) {
    let n = mp.atoms.len() as f64;
    for a in &mut mp.atoms {
        if a.color & KEEP_POSITION == 0 {
            a.x = (STDBOND * n.sqrt() * f64::from(rng.next_value()) / 2147483647.0) as f32;
            a.y = (STDBOND * n.sqrt() * f64::from(rng.next_value()) / 2147483647.0) as f32;
            a.z = 0.0;
        }
    }
}

/// `RecolorMolecule`: every atom and bond gets its own colour
pub fn recolor_molecule(mp: &mut Molecule) {
    let mut ncolor = 1;
    for a in &mut mp.atoms {
        a.color = ncolor;
        ncolor += 1;
    }
    for b in &mut mp.bonds {
        b.color = ncolor;
        ncolor += 1;
    }
}

/// `PutColorIntoValue`
pub fn put_color_into_value(mp: &mut Molecule) {
    for a in &mut mp.atoms {
        a.value = a.color as f32;
    }
}

/// `FloodColor`: colours the fragment containing atom `aindex` (not crossing rubber bonds) and
/// returns the number of recoloured atoms
pub fn flood_color(mp: &mut Molecule, nbp: &[Neighbourhood], aindex: usize, color: i32) -> i32 {
    mp.atoms[aindex].color = color;
    let mut result = 1;
    for i in 0..nbp[aindex].atoms.len() {
        let (at, bo) = (nbp[aindex].atoms[i], nbp[aindex].bonds[i]);
        if mp.atoms[at].color == NONE && mp.bonds[bo].bond_type & RUBBER_BOND == 0 {
            result += flood_color(mp, nbp, at, color);
        }
    }
    result
}

/// `FloodClearColor`: clears the colour of the connected atoms that have `color`
pub fn flood_clear_color(
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    aindex: usize,
    color: i32,
) -> i32 {
    if mp.atoms[aindex].color != color {
        return 0;
    }
    mp.atoms[aindex].color = NONE;
    let mut result = 1;
    for i in 0..nbp[aindex].atoms.len() {
        let (at, bo) = (nbp[aindex].atoms[i], nbp[aindex].bonds[i]);
        if mp.atoms[at].color == color && mp.bonds[bo].bond_type & RUBBER_BOND == 0 {
            result += flood_clear_color(mp, nbp, at, color);
        }
    }
    result
}

/// mirrors the atom `i` on the line through `p1` with direction `r12`
fn mirror_on_bond(mp: &mut Molecule, i: usize, p1: [f64; 2], r12: [f64; 2]) {
    let p = [f64::from(mp.atoms[i].x), f64::from(mp.atoms[i].y)];
    let q =
        ((p[0] - p1[0]) * r12[0] + (p[1] - p1[1]) * r12[1]) / (r12[0] * r12[0] + r12[1] * r12[1]);
    mp.atoms[i].x = (p1[0] - (p[0] - p1[0]) + 2.0 * r12[0] * q) as f32;
    mp.atoms[i].y = (p1[1] - (p[1] - p1[1]) + 2.0 * r12[1] * q) as f32;
}

/// the single precision strain term of `MakeBondTrans`
fn strain(mp: &Molecule, ia: &[usize], sa: &[i32], ib: &[usize], sb: &[i32]) -> f64 {
    let mut s = 0.0f64;
    for i in 0..ia.len() {
        for j in 0..ib.len() {
            let (a, b) = (&mp.atoms[ia[i]], &mp.atoms[ib[j]]);
            let d: f32 = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
            let den: f32 = 1.0 + d;
            let term: f32 = (sa[i] * sb[j]) as f32 / den;
            s += f64::from(term);
        }
    }
    s
}

/// `MakeBondTrans`: flips the ligands at one end of bond `bi` if that makes it trans
pub fn make_bond_trans(mp: &mut Molecule, nbp: &[Neighbourhood], bi: usize) {
    let atom_colors: Vec<i32> = mp.atoms.iter().map(|a| a.color).collect();
    for a in &mut mp.atoms {
        a.color = NONE;
    }
    let (aia, aib) = (
        (mp.bonds[bi].atoms[0] - 1) as usize,
        (mp.bonds[bi].atoms[1] - 1) as usize,
    );
    mp.atoms[aia].color = 1;
    mp.atoms[aib].color = 2;
    let (mut indexa, mut sizea, mut indexb, mut sizeb) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for &l in &nbp[aia].atoms {
        if mp.atoms[l].color == NONE {
            indexa.push(l);
            sizea.push(flood_color(mp, nbp, l, 3));
            flood_clear_color(mp, nbp, l, 3);
        }
    }
    for &l in &nbp[aib].atoms {
        if mp.atoms[l].color == NONE {
            indexb.push(l);
            sizeb.push(flood_color(mp, nbp, l, 3));
            flood_clear_color(mp, nbp, l, 3);
        }
    }
    let strain1 = strain(mp, &indexa, &sizea, &indexb, &sizeb);
    let p1 = [f64::from(mp.atoms[aia].x), f64::from(mp.atoms[aia].y)];
    let p2 = [f64::from(mp.atoms[aib].x), f64::from(mp.atoms[aib].y)];
    let r12 = [p2[0] - p1[0], p2[1] - p1[1]];
    for &i in &indexa {
        mirror_on_bond(mp, i, p1, r12);
    }
    let strain2 = strain(mp, &indexa, &sizea, &indexb, &sizeb);
    for &i in &indexa {
        mirror_on_bond(mp, i, p1, r12);
    }
    if strain1 > strain2 {
        for &i in &indexa {
            flood_color(mp, nbp, i, 3);
        }
        for i in 0..mp.atoms.len() {
            if mp.atoms[i].color == 3 {
                mirror_on_bond(mp, i, p1, r12);
            }
        }
    }
    for (a, c) in mp.atoms.iter_mut().zip(atom_colors) {
        a.color = c;
    }
}

/// `InvertFragmentColor`: multiplies the colours of the atoms with `color` that are reachable from
/// `aindex` by -1
pub fn invert_fragment_color(
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    aindex: usize,
    color: i32,
) -> i32 {
    let mut result = 0;
    for i in 0..nbp[aindex].atoms.len() {
        let (at, bo) = (nbp[aindex].atoms[i], nbp[aindex].bonds[i]);
        if mp.atoms[at].color == color && mp.bonds[bo].bond_type & RUBBER_BOND == 0 {
            result += 1;
            mp.atoms[at].color *= -1;
            result += invert_fragment_color(mp, nbp, at, color);
        }
    }
    result
}

/// `MoleculeRings`: the combined basis rings, ignoring rubber bonds
fn molecule_rings(mp: &Molecule) -> Vec<BondSet> {
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
                idummy += 1;
                let a = idummy;
                idummy += 1;
                [a, idummy]
            }
        })
        .collect();
    let mut rings = ring_list(&bonds);
    combine_rings(&mut rings);
    rings
}

/// `struct ring_node_t`
#[derive(Debug, Clone, Copy, Default)]
pub struct RingNode {
    pub xfirst: f32,
    pub yfirst: f32,
    pub xlast: f32,
    pub ylast: f32,
    pub aifirst: i32,
    pub ailast: i32,
    pub color: i32,
    pub d: f32,
    pub angle: f32,
}

/// `FillRingNodeTable`: the table of colour segments on the perimeter of `ring`
pub fn fill_ring_node_table(ring: &BondSet, mp: &Molecule) -> Vec<RingNode> {
    if mp.bonds.is_empty() {
        return Vec::new();
    }
    let mut bonds: Vec<[i32; 2]> = Vec::new();
    for (i, b) in mp.bonds.iter().enumerate() {
        if ring.bond_set.is_member(i) {
            bonds.push(b.atoms);
        }
    }
    let nrnt = bonds.len();
    for i in 1..nrnt {
        for j in i..nrnt {
            if bonds[j][0] == bonds[i - 1][1] {
                bonds.swap(i, j);
                break;
            } else if bonds[j][1] == bonds[i - 1][1] {
                bonds[j].swap(0, 1);
                bonds.swap(i, j);
                break;
            }
        }
    }
    let acolor = |a: i32| mp.atoms[(a - 1) as usize].color;
    let mut i = 1;
    while i < nrnt {
        if acolor(bonds[i - 1][0]) != acolor(bonds[i][0]) {
            break;
        }
        i += 1;
    }
    if i == nrnt {
        return Vec::new();
    }
    let from = i - 1;
    let mut i = nrnt - 1;
    while i > from {
        if acolor(bonds[i][0]) != acolor(bonds[from][0]) {
            break;
        }
        i -= 1;
    }
    let to = i;
    let mut table = vec![RingNode {
        aifirst: bonds[to][1],
        ailast: bonds[from][0],
        color: acolor(bonds[from][0]),
        ..Default::default()
    }];
    for i in from + 1..=to {
        if acolor(bonds[i][0]) != acolor(bonds[i - 1][0]) {
            let mut j = i;
            while j < to {
                if acolor(bonds[j][0]) != acolor(bonds[j][1]) {
                    break;
                }
                j += 1;
            }
            table.push(RingNode {
                aifirst: bonds[i][0],
                ailast: bonds[j][0],
                color: acolor(bonds[i][0]),
                ..Default::default()
            });
        }
    }
    for r in &mut table {
        let (f, l) = (
            &mp.atoms[(r.aifirst - 1) as usize],
            &mp.atoms[(r.ailast - 1) as usize],
        );
        r.xfirst = f.x;
        r.yfirst = f.y;
        r.xlast = l.x;
        r.ylast = l.y;
    }
    table
}

/// `LayoutRingSegment`: places the colour segments of a ring on a circle
pub fn layout_ring_segment(mp: &mut Molecule, rnt: &mut [RingNode]) {
    let nrnt = rnt.len();
    let (mut xc, mut yc) = (0.0f64, 0.0f64);
    let mut n = 0;
    for r in rnt.iter().skip(1) {
        for a in &mp.atoms {
            if a.color == r.color {
                xc += f64::from(a.x);
                yc += f64::from(a.y);
                n += 1;
            }
        }
    }
    xc /= f64::from(n);
    yc /= f64::from(n);

    let mut len = 0.0f64;
    for r in rnt.iter_mut() {
        len += STDBOND;
        let dx: f32 = r.xfirst - r.xlast;
        let dy: f32 = r.yfirst - r.ylast;
        let s: f32 = dx * dx + dy * dy;
        r.d = f64::from(s).sqrt() as f32;
        len += f64::from(r.d);
    }
    let mut peri = 0.0f64;
    for r in rnt.iter_mut() {
        peri += STDBOND * (PI * STDBOND / len) / (PI * STDBOND / len).sin();
        r.angle = (2.0 * PI * f64::from(r.d) / len) as f32;
        if f64::from(r.angle) > 0.0001 {
            let half: f32 = r.angle / 2.0;
            let prod: f32 = r.d * half;
            peri += f64::from(prod) / f64::from(half).sin();
        }
    }
    let r = peri / (2.0 * PI);

    let mut points = vec![[0.0f64; 2]; mp.atoms.len()];
    let mut angle = 0.0f64;
    for i in 0..nrnt {
        let x1 = r * angle.sin() + xc;
        let y1 = r * angle.cos() + yc;
        angle += f64::from(rnt[i].angle);
        let x2 = r * angle.sin() + xc;
        let y2 = r * angle.cos() + yc;
        angle += 2.0 * PI * STDBOND / len;

        let mut p1 = [f64::from(rnt[i].xfirst), f64::from(rnt[i].yfirst)];
        let mut p2 = [f64::from(rnt[i].xlast), f64::from(rnt[i].ylast)];
        let p1p = [x1, y1];
        let mut p2p = [x2, y2];
        let mut npts = 0;
        for a in &mp.atoms {
            if a.color == rnt[i].color {
                points[npts] = [f64::from(a.x), f64::from(a.y)];
                npts += 1;
            }
        }
        let sq = |v: f64| v * v;
        if sq(p1[0] - p2[0]) + sq(p1[1] - p2[1]) > 0.01 * sq(STDBOND) {
            // ring fusion
            let mut ps = [0.0f64; 2];
            let mut npts = 0;
            for a in &mp.atoms {
                if a.color == rnt[i].color {
                    ps[0] += f64::from(a.x);
                    ps[1] += f64::from(a.y);
                    npts += 1;
                }
            }
            ps[0] /= f64::from(npts);
            ps[1] /= f64::from(npts);
            if (ps[0] - p1[0]) * (p2[1] - p1[1]) - (ps[1] - p1[1]) * (p2[0] - p1[0]) > 0.0 {
                p1[1] *= -1.0;
                p2[1] *= -1.0;
                for pt in points.iter_mut().take(npts as usize) {
                    pt[1] *= -1.0;
                }
                flip_stereo_symbols(mp, rnt[i].color);
            }
        } else {
            // single atom or spiro
            let mut ps = [0.0f64; 2];
            let mut npts = 0;
            for b in &mp.bonds {
                let c0 = mp.atoms[(b.atoms[0] - 1) as usize].color;
                let c1 = mp.atoms[(b.atoms[1] - 1) as usize].color;
                let col = rnt[i].color;
                if (b.atoms[0] == rnt[i].aifirst && c1 == col)
                    || (b.atoms[1] == rnt[i].aifirst && c0 == col)
                    || (b.atoms[0] == rnt[i].ailast && c1 == col)
                    || (b.atoms[1] == rnt[i].ailast && c0 == col)
                {
                    ps[0] += f64::from(mp.atoms[(b.atoms[0] - 1) as usize].x);
                    ps[1] += f64::from(mp.atoms[(b.atoms[0] - 1) as usize].y);
                    npts += 1;
                    ps[0] += f64::from(mp.atoms[(b.atoms[1] - 1) as usize].x);
                    ps[1] += f64::from(mp.atoms[(b.atoms[1] - 1) as usize].y);
                    npts += 1;
                }
            }
            if npts > 0 {
                ps[0] /= f64::from(npts);
                ps[1] /= f64::from(npts);
                p2 = ps;
                p2p[0] = p1p[0]
                    + (p1p[0] - xc) * (sq(ps[0] - p1[0]) + sq(ps[1] - p1[1])).sqrt()
                        / (sq(p1p[0] - xc) + sq(p1p[1] - yc)).sqrt();
                p2p[1] = p1p[1]
                    + (p1p[1] - yc) * (sq(ps[0] - p1[0]) + sq(ps[1] - p1[1])).sqrt()
                        / (sq(p1p[0] - xc) + sq(p1p[1] - yc)).sqrt();
            }
        }
        transform_points(&mut points, p1, p2, p1p, p2p);
        let mut npts = 0;
        for a in &mut mp.atoms {
            if a.color == rnt[i].color {
                a.x = points[npts][0] as f32;
                a.y = points[npts][1] as f32;
                npts += 1;
            }
        }
    }
}

/// `MergeConnectedFragmentColors` including the assignment typo of the C code
/// (`if (color = col2)`), which paints every atom with the smaller colour
pub fn merge_connected_fragment_colors(mp: &mut Molecule) {
    loop {
        let mut found = None;
        for b in &mp.bonds {
            let c1 = mp.atoms[(b.atoms[0] - 1) as usize].color;
            let c2 = mp.atoms[(b.atoms[1] - 1) as usize].color;
            if c1 != c2 {
                found = Some((c1, c2));
                break;
            }
        }
        let Some((mut col1, mut col2)) = found else {
            break;
        };
        if col1 > col2 {
            std::mem::swap(&mut col1, &mut col2);
        }
        for a in &mut mp.atoms {
            a.color = col2;
            if col2 != 0 {
                a.color = col1;
            }
        }
    }
}

/// `MergeColors`: merges the colours of the table entries into the first one
pub fn merge_colors(mp: &mut Molecule, rnt: &[RingNode]) {
    for r in rnt.iter().skip(1) {
        for a in &mut mp.atoms {
            if a.color == r.color {
                a.color = rnt[0].color;
            }
        }
    }
}

/// `HideMetalComplexProblems`: bonds to metals in several rings and bonds to hydrogens become
/// rubber bonds
pub fn hide_metal_complex_problems(mp: &mut Molecule) {
    if mp.bonds.is_empty() {
        return;
    }
    if !mp
        .atoms
        .iter()
        .any(|a| atom_symbol_match(&a.atom_symbol, "trn,H,D,T"))
    {
        return;
    }
    let (_, bond_status) = ring_state(mp);
    for i in 0..mp.bonds.len() {
        let s0 = mp.atoms[(mp.bonds[i].atoms[0] - 1) as usize]
            .atom_symbol
            .clone();
        let s1 = mp.atoms[(mp.bonds[i].atoms[1] - 1) as usize]
            .atom_symbol
            .clone();
        if bond_status[i] > 1 {
            if atom_symbol_match(&s0, "trn") || atom_symbol_match(&s1, "trn") {
                mp.bonds[i].bond_type |= RUBBER_BOND;
            }
        } else if ["H", "D", "T"]
            .iter()
            .any(|p| atom_symbol_match(&s0, p) || atom_symbol_match(&s1, p))
        {
            mp.bonds[i].bond_type |= RUBBER_BOND;
        }
    }
}

/// `LayoutRingSort`: six-membered rings first, then by increasing size
pub fn layout_ring_sort(rings: &mut [BondSet]) {
    for _ in 0..rings.len() {
        for i in 0..rings.len().saturating_sub(1) {
            let (c, n) = (rings[i].cardinality, rings[i + 1].cardinality);
            if (c != 6 && n == 6) || (c != 6 && n != 6 && c > n) {
                rings.swap(i, i + 1);
            }
        }
    }
}

/// `LayoutRings`: coordinates for the atoms of ring systems; atoms placed relative to each
/// other end up with the same colour
pub fn layout_rings(mp: &mut Molecule) {
    hide_metal_complex_problems(mp);
    let mut rings = molecule_rings(mp);
    layout_ring_sort(&mut rings);
    let mut i = 0;
    while i < rings.len() {
        let mut rnt = fill_ring_node_table(&rings[i], mp);
        if rnt.is_empty() {
            i += 1;
        } else {
            layout_ring_segment(mp, &mut rnt);
            merge_colors(mp, &rnt);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::avalon::perceive::{DOUBLE, setup_neighbourhood};
    use crate::external::avalon::reaccsio::mol_str_to_mol;
    use crate::external::avalon::scan::fmt_gp;

    fn seeded(seed: u32, mp: &Molecule) -> Molecule {
        let mut c = mp.clone();
        random_coordinates(&mut c, &mut CRand::new(seed));
        c
    }

    fn dump(tag: &str, c: &Molecule) -> String {
        let g = |v: f32| fmt_gp(f64::from(v), 9);
        let mut out = format!("{tag}\n");
        for (i, a) in c.atoms.iter().enumerate() {
            out += &format!(
                "a {} {} {} {} {} {}\n",
                i + 1,
                g(a.x),
                g(a.y),
                g(a.z),
                a.color,
                a.value as i32
            );
        }
        for (i, b) in c.bonds.iter().enumerate() {
            out += &format!(
                "b {} {} {} {} {}\n",
                i + 1,
                b.atoms[0],
                b.atoms[1],
                b.bond_type,
                b.stereo_symbol
            );
        }
        out
    }

    /// the first part of the layout must leave the same coordinates and colours as the C code
    #[test]
    fn layout_part1_matches_c() {
        let text = include_str!("../../../../../tests/data/avalon_layout1.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let mp = mol_str_to_mol(input).expect("molecule");
            let mut out = String::new();

            let mut c = seeded(7, &mp);
            out += &dump("random", &c);
            recolor_molecule(&mut c);
            put_color_into_value(&mut c);
            out += &dump("recolor", &c);
            layout_rings(&mut c);
            out += &dump("rings", &c);

            let mut c = seeded(3, &mp);
            recolor_molecule(&mut c);
            merge_connected_fragment_colors(&mut c);
            out += &dump("merge", &c);

            let mut c = mp.clone();
            hide_metal_complex_problems(&mut c);
            out += &dump("hide", &c);

            let mut c = seeded(5, &mp);
            let nbp = setup_neighbourhood(&c).unwrap();
            for i in 0..c.bonds.len() {
                if c.bonds[i].bond_type == DOUBLE {
                    make_bond_trans(&mut c, &nbp, i);
                }
            }
            out += &dump("trans", &c);

            let mut c = mp.clone();
            let nbp = setup_neighbourhood(&c).unwrap();
            for a in &mut c.atoms {
                a.color = 0;
            }
            if !c.atoms.is_empty() {
                out += &format!("flood={}\n", flood_color(&mut c, &nbp, 0, 5));
                out += &format!("invert={}\n", invert_fragment_color(&mut c, &nbp, 0, 5));
                out += &dump("flooded", &c);
                let last = c.atoms.len() - 1;
                out += &format!("clear={}\n", flood_clear_color(&mut c, &nbp, last, 5));
                out += &dump("cleared", &c);
            }
            if out != want {
                bad.push(format!("{input}\nwant:\n{want}\ngot:\n{out}"));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }
}
