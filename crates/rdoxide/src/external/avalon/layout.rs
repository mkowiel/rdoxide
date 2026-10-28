//! Port of Avalon Toolkit's `layout.c` (2D coordinate generation), in slices. This slice
//! holds the colouring helpers and the ring system layout. Where C mixes `float` fields with
//! `double` arithmetic the same single precision steps are taken here.

#![allow(clippy::approx_constant)]

use super::geometry::{
    H_GEOMETRY, NO_RAY_CROSSING, USE_INWARDS, next_substituent_point, transform_points,
};
use super::graph::{BondSet, CRand, combine_rings, ring_list};
use super::perceive::{DOUBLE, Neighbourhood, SINGLE, TRIPLE, ring_state};
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

/// `GetColoredGraph`: the edges (renumbered) and the coordinates of the atoms with `color`;
/// `numbers[i]` is the node index of atom `i`
pub fn get_colored_graph(mp: &Molecule, color: i32) -> (Vec<[u32; 2]>, Vec<[f64; 2]>, Vec<usize>) {
    let mut edges = Vec::new();
    for b in &mp.bonds {
        if b.bond_type & RUBBER_BOND == 0
            && mp.atoms[(b.atoms[0] - 1) as usize].color == color
            && mp.atoms[(b.atoms[1] - 1) as usize].color == color
        {
            edges.push([b.atoms[0] as u32, b.atoms[1] as u32]);
        }
    }
    let mut coords = Vec::new();
    let mut numbers = vec![0usize; mp.atoms.len()];
    for (i, a) in mp.atoms.iter().enumerate() {
        if a.color == color {
            coords.push([f64::from(a.x), f64::from(a.y)]);
            numbers[i] = coords.len() - 1;
        }
    }
    for e in &mut edges {
        e[0] = numbers[(e[0] - 1) as usize] as u32;
        e[1] = numbers[(e[1] - 1) as usize] as u32;
    }
    (edges, coords, numbers)
}

/// `ColorStrain`: the repulsion measure between the parts coloured `col1` and `col2`
pub fn color_strain(mp: &Molecule, col1: i32, col2: i32) -> f64 {
    let mut result = 0.0f64;
    for (i1, a1) in mp.atoms.iter().enumerate() {
        if a1.color != col1 {
            continue;
        }
        for (i2, a2) in mp.atoms.iter().enumerate() {
            if a2.color == col2 && i1 != i2 {
                let dx: f32 = (a1.x - a2.x) * (a1.x - a2.x);
                let dy: f32 = (a1.y - a2.y) * (a1.y - a2.y);
                result += 1.0 / (0.1 + f64::from(dx) + f64::from(dy));
            }
        }
    }
    result
}

/// `ImproveBondByFlip`: flips one side of bond `bi` if the strain improves enough
pub fn improve_bond_by_flip(mp: &mut Molecule, bi: usize, threshold: f64) -> bool {
    if mp.bonds[bi].bond_type & DONT_FLIP_BOND != 0 {
        return false;
    }
    let (a0, a1) = (
        (mp.bonds[bi].atoms[0] - 1) as usize,
        (mp.bonds[bi].atoms[1] - 1) as usize,
    );
    let strain1 = color_strain(mp, mp.atoms[a0].color, mp.atoms[a1].color);
    let p1 = [f64::from(mp.atoms[a0].x), f64::from(mp.atoms[a0].y)];
    let p2 = [f64::from(mp.atoms[a1].x), f64::from(mp.atoms[a1].y)];
    let r12 = [p2[0] - p1[0], p2[1] - p1[1]];
    let col = mp.atoms[a1].color;
    for i in 0..mp.atoms.len() {
        if mp.atoms[i].color == col {
            mirror_on_bond(mp, i, p1, r12);
        }
    }
    let strain2 = color_strain(mp, mp.atoms[a0].color, mp.atoms[a1].color);
    if strain1 <= strain2 * (1.0 + threshold) + threshold {
        let col = mp.atoms[a1].color;
        for i in 0..mp.atoms.len() {
            if mp.atoms[i].color == col {
                mirror_on_bond(mp, i, p1, r12);
            }
        }
        return false;
    }
    true
}

/// `ImproveBondByStretch`: stretches bond `bi` by a third if the two sides collide
pub fn improve_bond_by_stretch(mp: &mut Molecule, bi: usize) {
    let (a0, a1) = (
        (mp.bonds[bi].atoms[0] - 1) as usize,
        (mp.bonds[bi].atoms[1] - 1) as usize,
    );
    let (col1, col2) = (mp.atoms[a0].color, mp.atoms[a1].color);
    let mut collision = false;
    for ap1 in &mp.atoms {
        if ap1.color != col1 {
            continue;
        }
        for ap2 in &mp.atoms {
            let d: f32 = (ap1.x - ap2.x) * (ap1.x - ap2.x) + (ap1.y - ap2.y) * (ap1.y - ap2.y);
            if ap2.color == col2 && f64::from(d) < 0.05 * STDBOND * STDBOND {
                collision = true;
                break;
            }
        }
    }
    if !collision {
        return;
    }
    let r12 = [
        f64::from(mp.atoms[a1].x - mp.atoms[a0].x),
        f64::from(mp.atoms[a1].y - mp.atoms[a0].y),
    ];
    let col = mp.atoms[a1].color;
    for a in &mut mp.atoms {
        if a.color == col {
            a.x = (f64::from(a.x) + 0.3 * r12[0]) as f32;
            a.y = (f64::from(a.y) + 0.3 * r12[1]) as f32;
        }
    }
}

/// one flip trial of `ImproveMoleculeByBondFlip` around the bond `bi`; `swap_ring_end` temporarily
/// exchanges the bond's atoms. Returns whether the flip was kept (None: nothing to flip).
fn flip_trial(
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    bi: usize,
    threshold: f64,
    swap_ring_end: bool,
) -> Option<bool> {
    let (ai0, ai1) = (
        (mp.bonds[bi].atoms[0] - 1) as usize,
        (mp.bonds[bi].atoms[1] - 1) as usize,
    );
    let oldcolor = mp.atoms[ai0].color;
    mp.atoms[ai0].color = 0;
    mp.atoms[ai1].color *= -1;
    let ninvert = invert_fragment_color(mp, nbp, ai1, oldcolor);
    if ninvert > 0 {
        mp.atoms[ai0].color = oldcolor;
        if swap_ring_end {
            mp.bonds[bi].atoms = [ai1 as i32 + 1, ai0 as i32 + 1];
        }
        let improved = improve_bond_by_flip(mp, bi, threshold);
        for a in &mut mp.atoms {
            if a.color == -oldcolor {
                a.color = oldcolor;
            }
        }
        if swap_ring_end {
            mp.bonds[bi].atoms = [ai0 as i32 + 1, ai1 as i32 + 1];
        }
        Some(improved)
    } else {
        mp.atoms[ai0].color = oldcolor;
        mp.atoms[ai1].color *= -1;
        None
    }
}

/// `ImproveMoleculeByBondFlip`: flips around bonds until no improvement is found; returns the
/// number of improvement steps
pub fn improve_molecule_by_bond_flip(
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    ring_size: &[i32],
    is_ring_atom: &[i32],
    frag_color: i32,
    only_single_bonds: bool,
) -> i32 {
    let mut ret = 0;
    for _ in 0..2 {
        let mut changed = false;
        for i in 0..mp.bonds.len() {
            let (ai0, ai1) = (
                (mp.bonds[i].atoms[0] - 1) as usize,
                (mp.bonds[i].atoms[1] - 1) as usize,
            );
            if mp.atoms[ai0].color & KEEP_POSITION != 0 || mp.atoms[ai1].color & KEEP_POSITION != 0
            {
                continue;
            }
            if only_single_bonds && mp.bonds[i].bond_type != SINGLE {
                continue;
            }
            if mp.bonds[i].bond_type & RUBBER_BOND == 0
                && ring_size[i] == 0
                && nbp[ai0].atoms.len() > 1
                && nbp[ai1].atoms.len() > 1
                && (frag_color == 0
                    || (mp.atoms[ai0].color == frag_color && mp.atoms[ai1].color == frag_color))
            {
                if let Some(true) = flip_trial(mp, nbp, i, 0.1, is_ring_atom[ai1] != 0) {
                    if frag_color == 0 {
                        changed = true;
                    }
                    ret += 1;
                }
            }
        }
        if mp.bonds.len() > 60 {
            changed = false;
        }
        if !changed {
            break;
        }
    }
    if mp.atoms.len() > 60 {
        return ret;
    }
    for pass in 0..2 {
        // pass 0: chain bonds at atoms with more than two neighbours; pass 1: all chain bonds
        if pass == 1 && mp.atoms.len() > 30 {
            return ret;
        }
        for _ in 0..2 {
            let mut changed = false;
            for i in 0..mp.bonds.len() {
                let (ai0, ai1) = (
                    (mp.bonds[i].atoms[0] - 1) as usize,
                    (mp.bonds[i].atoms[1] - 1) as usize,
                );
                let cond = mp.bonds[i].bond_type & RUBBER_BOND == 0
                    && is_ring_atom[ai0] == 0
                    && is_ring_atom[ai1] == 0
                    && nbp[ai0].atoms.len() > 1
                    && nbp[ai1].atoms.len() > 1
                    && (pass == 1 || nbp[ai1].atoms.len() > 2 || nbp[ai0].atoms.len() > 2)
                    && (frag_color == 0 || mp.atoms[ai0].color == frag_color/* C tests atom 0 twice here */);
                if cond {
                    if only_single_bonds && mp.bonds[i].bond_type != SINGLE {
                        continue;
                    }
                    if let Some(true) = flip_trial(mp, nbp, i, 1.0, false) {
                        if frag_color == 0 {
                            changed = true;
                        }
                        ret += 1;
                    }
                }
            }
            if !changed {
                break;
            }
        }
    }
    ret
}

/// `ImproveSPAtoms`: makes sp carbon and nitrogen atoms linear (only atom number `atom` if != 0)
pub fn improve_sp_atoms(mp: &mut Molecule, nbp: &[Neighbourhood], atom: usize) {
    let colors: Vec<i32> = mp.atoms.iter().map(|a| a.color).collect();
    for i in 0..mp.atoms.len() {
        if nbp[i].atoms.len() != 2 {
            continue;
        }
        let sym = mp.atoms[i].atom_symbol.as_str();
        if sym != "C" && sym != "N" {
            continue;
        }
        let (b0, b1) = (
            mp.bonds[nbp[i].bonds[0]].bond_type,
            mp.bonds[nbp[i].bonds[1]].bond_type,
        );
        if !((b0 == SINGLE && b1 == TRIPLE)
            || (b1 == SINGLE && b0 == TRIPLE)
            || (b0 == DOUBLE && b1 == DOUBLE))
        {
            continue;
        }
        if !(atom == 0 || atom == i + 1) {
            continue;
        }
        for a in &mut mp.atoms {
            a.color = NONE;
        }
        mp.atoms[i].color = 1;
        mp.atoms[nbp[i].atoms[0]].color = 2;
        mp.atoms[nbp[i].atoms[1]].color = 3;
        loop {
            let mut changed = false;
            for j in 0..mp.bonds.len() {
                let (x0, x1) = (
                    (mp.bonds[j].atoms[0] - 1) as usize,
                    (mp.bonds[j].atoms[1] - 1) as usize,
                );
                if mp.atoms[x0].color > 0 && mp.atoms[x1].color == 0 {
                    changed = true;
                    mp.atoms[x1].color = mp.atoms[x0].color;
                } else if mp.atoms[x1].color > 0 && mp.atoms[x0].color == 0 {
                    changed = true;
                    mp.atoms[x0].color = mp.atoms[x1].color;
                }
            }
            if !changed {
                break;
            }
        }
        if mp.atoms[nbp[i].atoms[0]].color != mp.atoms[nbp[i].atoms[1]].color {
            let ap = &mp.atoms[i];
            let p1 = [f64::from(ap.x), f64::from(ap.y)];
            let n0 = &mp.atoms[nbp[i].atoms[0]];
            let n1 = &mp.atoms[nbp[i].atoms[1]];
            let p2 = [f64::from(n0.x), f64::from(n0.y)];
            let p2p = [2.0 * p1[0] - f64::from(n1.x), 2.0 * p1[1] - f64::from(n1.y)];
            let (_, mut coords, numbers) = get_colored_graph(mp, 2);
            transform_points(&mut coords, p1, p2, p1, p2p);
            for j in 0..mp.atoms.len() {
                if mp.atoms[j].color == 2 {
                    mp.atoms[j].x = coords[numbers[j]][0] as f32;
                    mp.atoms[j].y = coords[numbers[j]][1] as f32;
                }
            }
        }
    }
    for (a, c) in mp.atoms.iter_mut().zip(colors) {
        a.color = c;
    }
}

/// `CountColor`: the number of atoms with `color`
pub fn count_color(mp: &Molecule, color: i32) -> usize {
    mp.atoms.iter().filter(|a| a.color == color).count()
}

/// the indices of the carbonyl O, amide N and second C ligand of the large-ring lactam
/// carbonyl carbon `i`, if it is one (shared by `LayoutLargeLactams` and `FlipRingCarbonyls`)
fn lactam_atoms(
    mp: &Molecule,
    i: usize,
    is_ring_atom: &[i32],
    ring_size: &[i32],
    nbp: &[Neighbourhood],
) -> Option<(usize, usize, usize)> {
    let nb = &nbp[i];
    if nb.atoms.len() != 3 || is_ring_atom[i] == 0 || mp.atoms[i].atom_symbol != "C" {
        return None;
    }
    let jo = (0..3).find(|&j| {
        mp.bonds[nb.bonds[j]].bond_type == DOUBLE
            && mp.atoms[nb.atoms[j]].atom_symbol == "O"
            && nbp[nb.atoms[j]].atoms.len() == 1
    })?;
    let nrings = |a: usize| {
        nbp[a]
            .atoms
            .iter()
            .filter(|&&x| is_ring_atom[x] != 0)
            .count()
    };
    let jn = (0..3).find(|&j| {
        mp.bonds[nb.bonds[j]].bond_type == SINGLE
            && mp.atoms[nb.atoms[j]].atom_symbol == "N"
            && ring_size[nb.bonds[j]] > 12
    })?;
    if nrings(nb.atoms[jn]) > 2 {
        return None;
    }
    let jx = (0..3).find(|&j| {
        mp.bonds[nb.bonds[j]].bond_type == SINGLE && mp.atoms[nb.atoms[j]].atom_symbol == "C"
    })?;
    if nrings(nb.atoms[jx]) > 2 {
        return None;
    }
    Some((nb.atoms[jo], nb.atoms[jn], nb.atoms[jx]))
}

/// `LayoutLargeLactams`: makes O=C-N of large lactam amides trigonal
pub fn layout_large_lactams(
    mp: &mut Molecule,
    is_ring_atom: &[i32],
    ring_size: &[i32],
    nbp: &[Neighbourhood],
) {
    for i in 0..mp.atoms.len() {
        let Some((o, n, x)) = lactam_atoms(mp, i, is_ring_atom, ring_size, nbp) else {
            continue;
        };
        let color = mp.atoms[i].color;
        mp.atoms[i].x = 0.0;
        mp.atoms[i].y = 0.0;
        mp.atoms[o].x = 1.514;
        mp.atoms[o].y = 0.0;
        mp.atoms[o].color = color;
        mp.atoms[n].x = (-1.514 * (3.1415926f64 / 3.0).cos()) as f32;
        mp.atoms[n].y = (1.514 * (3.1415926f64 / 3.0).sin()) as f32;
        mp.atoms[n].color = color;
        mp.atoms[x].x = (-1.514 * (3.1415926f64 / 3.0).cos()) as f32;
        mp.atoms[x].y = (-1.514 * (3.1415926f64 / 3.0).sin()) as f32;
        mp.atoms[x].color = color;
    }
}

/// `FlipRingCarbonyls`: lets laid-out large lactam carbonyls point into the ring
pub fn flip_ring_carbonyls(
    mp: &mut Molecule,
    is_ring_atom: &[i32],
    ring_size: &[i32],
    nbp: &[Neighbourhood],
) {
    for i in 0..mp.atoms.len() {
        let Some((o, n, x)) = lactam_atoms(mp, i, is_ring_atom, ring_size, nbp) else {
            continue;
        };
        let c = mp.atoms[i].color;
        if c != mp.atoms[o].color || c != mp.atoms[n].color || c != mp.atoms[x].color {
            continue;
        }
        let sx: f32 = mp.atoms[x].x + mp.atoms[n].x;
        let sy: f32 = mp.atoms[x].y + mp.atoms[n].y;
        let (xref, yref) = (f64::from(sx) / 2.0, f64::from(sy) / 2.0);
        for k in [i, o] {
            mp.atoms[k].x = (xref - (f64::from(mp.atoms[k].x) - xref)) as f32;
            mp.atoms[k].y = (yref - (f64::from(mp.atoms[k].y) - yref)) as f32;
        }
    }
}

/// `SproutRingSubstituents`: attaches the first atom of non-ring substituents to laid-out
/// rings in a symmetrical way (only at ring atoms with exclusively non-ring attachments)
pub fn sprout_ring_substituents(
    mp: &mut Molecule,
    is_ring_atom: &[i32],
    ring_size: &[i32],
    nbp: &[Neighbourhood],
) {
    if mp.bonds.is_empty() {
        return;
    }
    for j in 0..mp.atoms.len() {
        if is_ring_atom[j] == 0 {
            continue;
        }
        let mut subst: Vec<usize> = Vec::new();
        let mut n_rbonds = 0;
        let mut problem = NONE;
        let mut min_rsize = 100;
        for i in 0..nbp[j].atoms.len() {
            let bo = nbp[j].bonds[i];
            if mp.bonds[bo].bond_type & RUBBER_BOND != 0 {
                continue;
            }
            if ring_size[bo] != 0 {
                if ring_size[bo] < min_rsize {
                    min_rsize = ring_size[bo];
                }
                n_rbonds += 1;
            } else {
                let iatom = nbp[j].atoms[i];
                if is_ring_atom[iatom] != 0 {
                    problem = problem.max(NONE + 1);
                } else {
                    if count_color(mp, mp.atoms[iatom].color) != 1 {
                        problem = problem.max(NONE + 2);
                    }
                    subst.push(iatom);
                }
            }
        }
        let n_subst = subst.len();
        if n_subst == 0 || n_rbonds != 2 || problem != NONE {
            continue;
        }
        let (edges, coords, numbers) = get_colored_graph(mp, mp.atoms[j].color);
        let seed = numbers[j];
        let p1 = coords[seed];
        let mut p2 = [0.0f64; 2];
        let flags = if min_rsize <= 8 { 0 } else { USE_INWARDS };
        next_substituent_point(&mut p2, &coords, &edges, seed, flags);
        let d = [p2[0] - p1[0], p2[1] - p1[1]];
        let ns = n_subst as f64;
        if n_subst > 1 {
            for i in 1..n_subst {
                for k in (0..i).rev() {
                    if nbp[subst[k]].atoms.len() < nbp[subst[k + 1]].atoms.len() {
                        subst.swap(k, k + 1);
                    } else {
                        break;
                    }
                }
            }
            let ang = (0.5 * PI / ns - PI / 2.0) * 80.0 / 90.0;
            let (cosa, sina) = (ang.cos(), ang.sin());
            let dp = [
                cosa * d[0] - sina * d[1] + p1[0],
                sina * d[0] + cosa * d[1] + p1[1],
            ];
            let mut c = [0.0f64; 2];
            for co in &coords {
                c[0] += co[0];
                c[1] += co[1];
            }
            c[0] /= coords.len() as f64;
            c[1] /= coords.len() as f64;
            if (dp[0] - c[0]) * (dp[0] - c[0]) + (dp[1] - c[1]) * (dp[1] - c[1])
                < (p2[0] - c[0]) * (p2[0] - c[0]) + (p2[1] - c[1]) * (p2[1] - c[1])
            {
                subst.swap(0, n_subst - 1);
            }
        }
        let color = mp.atoms[j].color;
        for (i, &sa) in subst.iter().enumerate() {
            let ang = ((i as f64 + 0.5) * PI / ns - PI / 2.0) * 80.0 / 90.0;
            let (cosa, sina) = (ang.cos(), ang.sin());
            mp.atoms[sa].x = (cosa * d[0] - sina * d[1] + p1[0]) as f32;
            mp.atoms[sa].y = (sina * d[0] + cosa * d[1] + p1[1]) as f32;
            mp.atoms[sa].color = color;
        }
    }
}

/// the two-fragment move shared by `LinkRingFragments` and `LinkRemainingFragments`: the
/// smaller fragment (colour `col1`, bond atom 0) is moved next to the reference fragment
fn move_fragment(
    mp: &mut Molecule,
    bi: usize,
    col1: i32,
    col2: i32,
    flags_ref: i32,
    flags_move: i32,
) {
    let (a0, a1) = (
        (mp.bonds[bi].atoms[0] - 1) as usize,
        (mp.bonds[bi].atoms[1] - 1) as usize,
    );
    let (edges, coords, numbers) = get_colored_graph(mp, col2);
    let seed = numbers[a1];
    let p2p = coords[seed];
    let mut p1p = [0.0f64; 2];
    next_substituent_point(&mut p1p, &coords, &edges, seed, flags_ref);
    let (edges, mut coords, numbers) = get_colored_graph(mp, col1);
    let seed = numbers[a0];
    let p1 = coords[seed];
    let mut p2 = [0.0f64; 2];
    next_substituent_point(&mut p2, &coords, &edges, seed, flags_move);
    transform_points(&mut coords, p1, p2, p1p, p2p);
    for i in 0..mp.atoms.len() {
        if mp.atoms[i].color == col1 {
            mp.atoms[i].x = coords[numbers[i]][0] as f32;
            mp.atoms[i].y = coords[numbers[i]][1] as f32;
        }
    }
}

/// orders the bond so that the larger fragment comes second; returns (col1, col2)
fn order_bond_by_fragment_size(mp: &mut Molecule, bi: usize) -> (i32, i32) {
    let mut col1 = mp.atoms[(mp.bonds[bi].atoms[0] - 1) as usize].color;
    let mut col2 = mp.atoms[(mp.bonds[bi].atoms[1] - 1) as usize].color;
    if count_color(mp, col1) > count_color(mp, col2) {
        mp.bonds[bi].atoms.swap(0, 1);
        std::mem::swap(&mut col1, &mut col2);
    }
    (col1, col2)
}

/// `LinkRingFragments`: links directly connected ring fragments
pub fn link_ring_fragments(mp: &mut Molecule, is_ring_atom: &[i32]) {
    if mp.bonds.is_empty() {
        return;
    }
    for j in 0..mp.bonds.len() {
        let (a0, a1) = (
            (mp.bonds[j].atoms[0] - 1) as usize,
            (mp.bonds[j].atoms[1] - 1) as usize,
        );
        if mp.atoms[a0].color == mp.atoms[a1].color
            || mp.bonds[j].bond_type & RUBBER_BOND != 0
            || is_ring_atom[a0] == 0
            || is_ring_atom[a1] == 0
        {
            continue;
        }
        let (col1, col2) = order_bond_by_fragment_size(mp, j);
        move_fragment(mp, j, col1, col2, NO_RAY_CROSSING | USE_INWARDS, 0);
        improve_bond_by_flip(mp, j, 0.01);
        improve_bond_by_stretch(mp, j);
        for a in &mut mp.atoms {
            if a.color == col2 {
                a.color = col1;
            }
        }
    }
}

/// `LinkRemainingFragments`: links the fragments that are not laid out yet
pub fn link_remaining_fragments(mp: &mut Molecule) {
    if mp.bonds.is_empty() {
        return;
    }
    for j in 0..mp.bonds.len() {
        let (a0, a1) = (
            (mp.bonds[j].atoms[0] - 1) as usize,
            (mp.bonds[j].atoms[1] - 1) as usize,
        );
        if mp.atoms[a0].color != mp.atoms[a1].color
            && mp.bonds[j].bond_type & RUBBER_BOND != 0
            && (atom_symbol_match(&mp.atoms[a0].atom_symbol, "H,D,T")
                || atom_symbol_match(&mp.atoms[a1].atom_symbol, "H,D,T"))
        {
            mp.bonds[j].bond_type &= !RUBBER_BOND;
        }
    }
    for j in 0..mp.bonds.len() {
        let (a0, a1) = (
            (mp.bonds[j].atoms[0] - 1) as usize,
            (mp.bonds[j].atoms[1] - 1) as usize,
        );
        if mp.atoms[a0].color == mp.atoms[a1].color || mp.bonds[j].bond_type & RUBBER_BOND != 0 {
            continue;
        }
        let (col1, col2) = order_bond_by_fragment_size(mp, j);
        let (b0, b1) = (
            (mp.bonds[j].atoms[0] - 1) as usize,
            (mp.bonds[j].atoms[1] - 1) as usize,
        );
        let h_flags = |sym: &str| {
            if sym == "H" {
                USE_INWARDS | H_GEOMETRY
            } else {
                USE_INWARDS
            }
        };
        let flags_ref = h_flags(&mp.atoms[b1].atom_symbol);
        let flags_move = h_flags(&mp.atoms[b0].atom_symbol);
        move_fragment(mp, j, col1, col2, flags_ref, flags_move);
        for a in &mut mp.atoms {
            if a.color == col1 {
                a.color = col2;
            }
        }
    }
}

/// `BranchQuality`: how well atom `ai` is suited to be laid out next (-1: not at all)
pub fn branch_quality(mp: &Molecule, ai: usize, nbp: &[Neighbourhood], used_atoms: &[i32]) -> i32 {
    let mut nbranch = vec![0i32; nbp[ai].atoms.len()];
    let mut already_layouted = 0;
    let mut can_be_candidate = false;
    let mut placed_neighbour = false;
    for j in 0..nbp[ai].atoms.len() {
        if mp.bonds[nbp[ai].bonds[j]].bond_type & RUBBER_BOND != 0 {
            continue;
        }
        if mp.atoms[ai].color != mp.atoms[nbp[ai].atoms[j]].color {
            can_be_candidate = true;
        } else {
            already_layouted += 1;
        }
        if used_atoms[nbp[ai].atoms[j]] != 0 {
            placed_neighbour = true;
        }
        let nbph = &nbp[nbp[ai].atoms[j]];
        nbranch[j] = nbph
            .bonds
            .iter()
            .filter(|&&b| mp.bonds[b].bond_type & RUBBER_BOND == 0)
            .count() as i32;
    }
    if !can_be_candidate || already_layouted > 1 {
        return -1;
    }
    let (mut nsingle, mut nmore, mut ntotal) = (0, 0, 0);
    for &nb in &nbranch {
        if nb == 1 {
            nsingle += 1;
        } else if nb > 0 {
            nmore += 1;
        }
        ntotal += nb;
    }
    if placed_neighbour {
        1000 + nsingle * 100 + nmore * 10 + ntotal
    } else {
        nsingle * 100 + nmore * 10 + ntotal
    }
}

/// `ColorSize`: the number of atoms with colour `col`
pub fn color_size(mp: &Molecule, col: i32) -> usize {
    count_color(mp, col)
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

    fn graph_dump(mp: &Molecule) -> String {
        let col = mp.atoms[0].color;
        let (edges, coords, _) = get_colored_graph(mp, col);
        let mut out = format!("graph {} {}\n", edges.len(), coords.len());
        for e in &edges {
            out += &format!("e {} {}\n", e[0], e[1]);
        }
        for p in &coords {
            out += &format!("p {} {}\n", fmt_gp(p[0], 17), fmt_gp(p[1], 17));
        }
        out
    }

    /// the second part of the layout must leave the same coordinates as the C code
    #[test]
    fn layout_part2_matches_c() {
        let text = include_str!("../../../../../tests/data/avalon_layout2.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let mp = mol_str_to_mol(input).expect("molecule");
            let mut out = String::new();

            let mut base = seeded(7, &mp);
            recolor_molecule(&mut base);
            layout_rings(&mut base);
            out += &dump("base", &base);
            for i in 0..base.atoms.len().min(3) {
                out += &format!("count{i}={}\n", count_color(&base, base.atoms[i].color));
            }
            let nbp = setup_neighbourhood(&base).unwrap();
            let (is_ring_atom, ring_size) = ring_state(&base);
            for k in 0..4 {
                let mut c = base.clone();
                let frag = if k & 1 != 0 { base.atoms[0].color } else { 0 };
                let ret = improve_molecule_by_bond_flip(
                    &mut c,
                    &nbp,
                    &ring_size,
                    &is_ring_atom,
                    frag,
                    k & 2 != 0,
                );
                out += &format!("improve{k}={ret}\n");
                out += &dump("improved", &c);
            }
            let mut c = base.clone();
            improve_sp_atoms(&mut c, &nbp, 0);
            out += &dump("sp", &c);
            for i in 0..base.bonds.len().min(8) {
                let mut c = base.clone();
                let a0 = (base.bonds[i].atoms[0] - 1) as usize;
                let a1 = (base.bonds[i].atoms[1] - 1) as usize;
                out += &format!(
                    "strain{i}={}\n",
                    fmt_gp(color_strain(&c, c.atoms[a0].color, c.atoms[a1].color), 17)
                );
                out += &format!(
                    "flip{i}={}\n",
                    i32::from(improve_bond_by_flip(&mut c, i, 0.1))
                );
                out += &dump("flipped", &c);
                let mut c = base.clone();
                improve_bond_by_stretch(&mut c, i);
                out += &dump("stretched", &c);
            }
            if !base.atoms.is_empty() {
                out += &graph_dump(&base);
            }
            if out != want {
                bad.push(format!("{input}\nwant:\n{want}\ngot:\n{out}"));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }

    /// the third part of the layout must leave the same coordinates as the C code
    #[test]
    fn layout_part3_matches_c() {
        let text = include_str!("../../../../../tests/data/avalon_layout3.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let mp = mol_str_to_mol(input).expect("molecule");
            let mut base = seeded(7, &mp);
            recolor_molecule(&mut base);
            layout_rings(&mut base);
            let nbp = setup_neighbourhood(&base).unwrap();
            let (as_, bs) = ring_state(&base);
            let mut out = dump("base", &base);
            let rs: Vec<i32> = bs.iter().map(|&b| if b != 0 { 14 } else { 0 }).collect();

            let mut c = base.clone();
            layout_large_lactams(&mut c, &as_, &rs, &nbp);
            out += &dump("lactam", &c);
            flip_ring_carbonyls(&mut c, &as_, &rs, &nbp);
            out += &dump("carbonyl", &c);
            let mut c = base.clone();
            flip_ring_carbonyls(&mut c, &as_, &rs, &nbp);
            out += &dump("carbonyl0", &c);
            let mut c = base.clone();
            sprout_ring_substituents(&mut c, &as_, &bs, &nbp);
            out += &dump("sprout", &c);
            let mut c = base.clone();
            link_ring_fragments(&mut c, &as_);
            out += &dump("linkring", &c);
            let mut c = base.clone();
            link_remaining_fragments(&mut c);
            out += &dump("linkrest", &c);
            let mut c = base.clone();
            sprout_ring_substituents(&mut c, &as_, &bs, &nbp);
            link_ring_fragments(&mut c, &as_);
            link_remaining_fragments(&mut c);
            out += &dump("pipeline", &c);
            let used: Vec<i32> = (0..base.atoms.len() as i32).map(|i| i % 2).collect();
            let mut c = base.clone();
            sprout_ring_substituents(&mut c, &as_, &bs, &nbp);
            for i in 0..c.atoms.len() {
                out += &format!(
                    "bq{i}={} size={}\n",
                    branch_quality(&c, i, &nbp, &used),
                    color_size(&c, c.atoms[i].color)
                );
            }
            if out != want {
                bad.push(format!("{input}\nwant:\n{want}\ngot:\n{out}"));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }
}
