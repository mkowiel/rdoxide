//! Port of Avalon Toolkit's `layout.c` (2D coordinate generation), in slices. This slice
//! holds the colouring helpers and the ring system layout. Where C mixes `float` fields with
//! `double` arithmetic the same single precision steps are taken here.

#![allow(clippy::approx_constant)]

use super::geometry::{
    H_GEOMETRY, NO_RAY_CROSSING, USE_INWARDS, next_substituent_point, transform_points,
};
use super::graph::{BondSet, CRand, combine_rings_rng, ring_list};
use super::perceive::{DOUBLE, Neighbourhood, SINGLE, TRIPLE, ring_state, setup_neighbourhood};
use super::reaccs::{CIS_TRANS_EITHER, CIS_TRANS_SWAPPED, Molecule, NONE};
use super::utilities::{DOWN, UP, atom_symbol_match, flip_stereo_symbols, reset_colors};

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
fn molecule_rings(mp: &Molecule, rng: &mut CRand) -> Vec<BondSet> {
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
    combine_rings_rng(&mut rings, rng);
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
    layout_rings_rng(mp, &mut CRand::new(1));
}

/// `LayoutRings` that leaves the C `rand()` state (reseeded by the ring perception) in `rng`
pub fn layout_rings_rng(mp: &mut Molecule, rng: &mut CRand) {
    hide_metal_complex_problems(mp);
    let mut rings = molecule_rings(mp, rng);
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

/// `LayoutChainAtoms`: arranges the environments of atoms that have only chain substituents
pub fn layout_chain_atoms(mp: &mut Molecule, is_ring_atom: &[i32], nbp: &[Neighbourhood]) {
    if mp.bonds.is_empty() {
        return;
    }
    let n = mp.atoms.len();
    let mut used_atoms = vec![0i32; n];
    loop {
        let mut quality = -1;
        let mut ibest = 0;
        for i in 0..n {
            if is_ring_atom[i] == 0 {
                let q = branch_quality(mp, i, nbp, &used_atoms);
                if quality < q {
                    ibest = i;
                    quality = q;
                }
            }
        }
        if quality < 0 {
            break;
        }
        used_atoms[ibest] = 1;

        let (mut index, mut bindex, mut col, mut size): (
            Vec<usize>,
            Vec<usize>,
            Vec<i32>,
            Vec<i32>,
        ) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let oldcolor = mp.atoms[ibest].color;
        mp.atoms[ibest].color = -1;
        for i in 0..nbp[ibest].atoms.len() {
            if mp.bonds[nbp[ibest].bonds[i]].bond_type & RUBBER_BOND == 0 {
                let a = nbp[ibest].atoms[i];
                index.push(a);
                bindex.push(nbp[ibest].bonds[i]);
                let c = mp.atoms[a].color;
                col.push(c);
                size.push(100 * color_size(mp, c) as i32 + nbp[a].atoms.len() as i32);
            }
        }
        let nneigh = index.len();
        for i in 1..nneigh {
            if col[i] == oldcolor {
                index.swap(i, 0);
                bindex.swap(i, 0);
                col.swap(i, 0);
                size.swap(i, 0);
            }
        }
        loop {
            let mut i = 1;
            while i + 1 < nneigh {
                if (size[i] - size[i - 1]).abs() + (size[i + 1] - size[(i + 2) % nneigh]).abs()
                    < (size[i] - size[(i + 2) % nneigh]).abs() + (size[i + 1] - size[i - 1]).abs()
                {
                    index.swap(i, i + 1);
                    bindex.swap(i, i + 1);
                    col.swap(i, i + 1);
                    size.swap(i, i + 1);
                    break;
                }
                i += 1;
            }
            if i + 1 >= nneigh {
                break;
            }
        }

        if nneigh == 1 {
            let (edges, coords, numbers) = get_colored_graph(mp, col[0]);
            let seed = numbers[index[0]];
            let mut p1 = [0.0f64; 2];
            next_substituent_point(&mut p1, &coords, &edges, seed, USE_INWARDS);
            mp.atoms[ibest].x = p1[0] as f32;
            mp.atoms[ibest].y = p1[1] as f32;
            mp.atoms[ibest].color = mp.atoms[index[0]].color;
        } else {
            let (edges, mut coords, numbers) = get_colored_graph(mp, col[0]);
            let nn = nneigh as f64;
            let p1p = [
                (2.0 * PI * 0.0 / nn).cos() * STDBOND,
                (2.0 * PI * 0.0 / nn).sin() * STDBOND,
            ];
            let p2p = [0.0, 0.0];
            let p1 = [
                f64::from(mp.atoms[index[0]].x),
                f64::from(mp.atoms[index[0]].y),
            ];
            let mut p2 = [0.0f64; 2];
            if col[0] == oldcolor {
                p2 = [f64::from(mp.atoms[ibest].x), f64::from(mp.atoms[ibest].y)];
            } else {
                let seed = numbers[index[0]];
                next_substituent_point(&mut p2, &coords, &edges, seed, USE_INWARDS);
            }
            transform_points(&mut coords, p1, p2, p1p, p2p);
            for j in 0..n {
                if mp.atoms[j].color == col[0] {
                    mp.atoms[j].x = coords[numbers[j]][0] as f32;
                    mp.atoms[j].y = coords[numbers[j]][1] as f32;
                }
            }
            mp.atoms[ibest].x = 0.0;
            mp.atoms[ibest].y = 0.0;

            for i in 1..nneigh {
                let (edges, mut coords, numbers) = get_colored_graph(mp, col[i]);
                let m = f64::from(nneigh.max(3) as i32);
                let p1p = [
                    (2.0 * PI * i as f64 / m).cos() * STDBOND,
                    (2.0 * PI * i as f64 / m).sin() * STDBOND,
                ];
                let p2p = [0.0, 0.0];
                let p1 = [
                    f64::from(mp.atoms[index[i]].x),
                    f64::from(mp.atoms[index[i]].y),
                ];
                let seed = numbers[index[i]];
                let mut p2 = [0.0f64; 2];
                next_substituent_point(&mut p2, &coords, &edges, seed, USE_INWARDS);
                transform_points(&mut coords, p1, p2, p1p, p2p);
                for j in 0..n {
                    if mp.atoms[j].color == col[i] {
                        mp.atoms[j].x = coords[numbers[j]][0] as f32;
                        mp.atoms[j].y = coords[numbers[j]][1] as f32;
                    }
                }
            }
            for &c in &col {
                for a in &mut mp.atoms {
                    if a.color == c {
                        a.color = oldcolor;
                    }
                }
            }
            mp.atoms[ibest].color = oldcolor;
        }

        if nneigh == 2
            && mp.bonds[bindex[0]].bond_type + mp.bonds[bindex[1]].bond_type > SINGLE + SINGLE
        {
            improve_sp_atoms(mp, nbp, ibest + 1);
        }
        for i in 0..nneigh {
            if mp.bonds[bindex[i]].bond_type == TRIPLE {
                improve_sp_atoms(mp, nbp, index[i] + 1);
            }
        }
        for i in 0..nneigh {
            if size[i] > 200 && is_ring_atom[index[i]] == 0 {
                make_bond_trans(mp, nbp, bindex[i]);
            }
        }
    }
}

/// `LayoutRubberFragments`: places the fragments connected by rubber bonds
pub fn layout_rubber_fragments(mp: &mut Molecule) {
    if mp.bonds.is_empty() {
        return;
    }
    loop {
        let mut found = None;
        for i in 0..mp.bonds.len() {
            let (a0, a1) = (
                (mp.bonds[i].atoms[0] - 1) as usize,
                (mp.bonds[i].atoms[1] - 1) as usize,
            );
            let not_shortcut = |a: &super::reaccs::Atom| a.atom_symbol != "R" || a.atext.is_empty();
            if mp.bonds[i].bond_type & RUBBER_BOND != 0
                && not_shortcut(&mp.atoms[a0])
                && not_shortcut(&mp.atoms[a1])
                && mp.atoms[a0].color > 0
                && mp.atoms[a1].color > 0
                && mp.atoms[a0].color != mp.atoms[a1].color
            {
                let (mut col1, mut col2) = (mp.atoms[a0].color, mp.atoms[a1].color);
                let (mut n1, mut n2) = (0, 0);
                for bh in &mp.bonds {
                    let (h0, h1) = ((bh.atoms[0] - 1) as usize, (bh.atoms[1] - 1) as usize);
                    if bh.atoms[0] == mp.bonds[i].atoms[0] && mp.atoms[h1].color == col1 {
                        n1 += 1;
                    }
                    if bh.atoms[1] == mp.bonds[i].atoms[0] && mp.atoms[h0].color == col1 {
                        n1 += 1;
                    }
                    if bh.atoms[0] == mp.bonds[i].atoms[1] && mp.atoms[h1].color == col2 {
                        n2 += 1;
                    }
                    if bh.atoms[1] == mp.bonds[i].atoms[1] && mp.atoms[h0].color == col2 {
                        n2 += 1;
                    }
                }
                if n2 == 0 && n1 > 0 {
                    found = Some((i, col1, col2));
                    break;
                } else if n1 == 0 && n2 > 0 {
                    std::mem::swap(&mut col1, &mut col2);
                    mp.bonds[i].atoms.swap(0, 1);
                    found = Some((i, col1, col2));
                    break;
                }
            }
        }
        let Some((bi, col1, col2)) = found else { break };
        let a0 = (mp.bonds[bi].atoms[0] - 1) as usize;
        let (edges, coords, numbers) = get_colored_graph(mp, col1);
        let mut p2p = [0.0f64; 2];
        next_substituent_point(&mut p2p, &coords, &edges, numbers[a0], USE_INWARDS);
        let p1p = coords[numbers[a0]];
        mp.bonds[bi].bond_type &= !RUBBER_BOND;

        let mut p2 = [0.0f64; 2];
        let mut n2 = 0;
        for a in &mp.atoms {
            if a.color == col2 {
                p2[0] += f64::from(a.x);
                p2[1] += f64::from(a.y);
                n2 += 1;
            }
        }
        p2[0] /= f64::from(n2);
        p2[1] /= f64::from(n2);

        let (edges, mut coords, numbers) = get_colored_graph(mp, col2);
        coords.push(p2);
        let nnodes = coords.len();
        let mut p1 = [0.0f64; 2];
        next_substituent_point(&mut p1, &coords, &edges, nnodes - 1, USE_INWARDS);
        let stretch_factor = (nnodes as f64 / 1.5).sqrt();
        p1[0] = p2[0] + stretch_factor * (p1[0] - p2[0]);
        p1[1] = p2[1] + stretch_factor * (p1[1] - p2[1]);
        p2p[0] = p1p[0] + stretch_factor * (p2p[0] - p1p[0]);
        p2p[1] = p1p[1] + stretch_factor * (p2p[1] - p1p[1]);
        transform_points(&mut coords, p1, p2, p1p, p2p);
        for j in 0..mp.atoms.len() {
            if mp.atoms[j].color == col2 {
                mp.atoms[j].x = coords[numbers[j]][0] as f32;
                mp.atoms[j].y = coords[numbers[j]][1] as f32;
                mp.atoms[j].color = -col1;
            }
        }
    }
    for a in &mut mp.atoms {
        if a.color < 0 {
            a.color *= -1;
        }
    }
}

/// `BestBondColor`: the highest bond colour and the number of bonds that have it
pub fn best_bond_color(mp: &Molecule) -> (i32, i32) {
    if mp.bonds.is_empty() {
        return (0, 0);
    }
    let mut best = mp.bonds[0].color;
    let mut nbest = 1;
    for b in &mp.bonds {
        if b.color > best {
            best = b.color;
            nbest = 1;
        } else if b.color == best {
            nbest += 1;
        }
    }
    (best, nbest)
}

/// `Diameter`: the extent of the fragment perpendicular to bond `bi` and the projection of the
/// direction to `center` (stored in `proj`); the bond is turned to point towards the center
pub fn diameter(mp: &mut Molecule, bi: usize, center: [f64; 2], proj: &mut f64) -> f64 {
    let pt = |mp: &Molecule, a: i32| {
        let at = &mp.atoms[(a - 1) as usize];
        [f64::from(at.x), f64::from(at.y)]
    };
    let mut p1 = pt(mp, mp.bonds[bi].atoms[0]);
    let mut p2 = pt(mp, mp.bonds[bi].atoms[1]);
    if (p2[0] - p1[0]) * (center[0] - (p1[0] + p2[0]) / 2.0)
        + (p2[1] - p1[1]) * (center[1] - (p1[1] + p2[1]) / 2.0)
        < 0.0
    {
        mp.bonds[bi].atoms.swap(0, 1);
        p1 = pt(mp, mp.bonds[bi].atoms[0]);
        p2 = pt(mp, mp.bonds[bi].atoms[1]);
    }
    let mut e12 = [-(p2[1] - p1[1]), p2[0] - p1[0]];
    let mut e11 = [
        (p2[0] - p1[0]) + (p2[1] - p1[1]),
        (p2[1] - p1[1]) - (p2[0] - p1[0]),
    ];
    let l = (e12[0] * e12[0] + e12[1] * e12[1]).sqrt();
    if l < 0.00001 {
        *proj = 0.0;
        return 0.0;
    }
    e12[0] /= l;
    e12[1] /= l;
    e11[0] /= l * 2.0f64.sqrt();
    e11[1] /= l * 2.0f64.sqrt();
    let color = mp.atoms[(mp.bonds[bi].atoms[0] - 1) as usize].color;
    let (mut minx, mut maxx) = (1.0e7f64, -1.0e7f64);
    for a in &mp.atoms {
        if a.color == color {
            let x = (f64::from(a.x) - p1[0]) * e12[0] + (f64::from(a.y) - p1[1]) * e12[1];
            if x < minx {
                minx = x;
            }
            if x > maxx {
                maxx = x;
            }
        }
    }
    *proj =
        (center[0] - (p1[0] + p2[0]) / 2.0) * e11[0] + (center[1] - (p1[1] + p2[1]) / 2.0) * e11[1];
    maxx - minx
}

/// `OrientFragment`: orients the fragment with `color` (landscape, preferred bond vertical)
pub fn orient_fragment(mp: &mut Molecule, color: i32, ring_size: &[i32], ring_count: &[i32]) {
    for i in 0..mp.bonds.len() {
        let (a0, a1) = (
            (mp.bonds[i].atoms[0] - 1) as usize,
            (mp.bonds[i].atoms[1] - 1) as usize,
        );
        mp.bonds[i].color = if mp.atoms[a0].color != color || mp.atoms[a1].color != color {
            0
        } else {
            1 + ring_size[i] + 10 * ring_count[i]
        };
    }
    let (best_color, _) = best_bond_color(mp);
    if best_color == 0 {
        return;
    }
    let mut center = [0.0f64; 2];
    if best_color == 1 {
        let mut n = 0.0f64;
        for a in &mp.atoms {
            if a.color == color {
                n += 1.0;
                center[0] += f64::from(a.x);
                center[1] += f64::from(a.y);
            }
        }
        center[0] /= n;
        center[1] /= n;
    } else {
        for i in 0..mp.bonds.len() {
            if ring_size[i] > 0 {
                for end in 0..2 {
                    let a = (mp.bonds[i].atoms[end] - 1) as usize;
                    if mp.atoms[a].color == color {
                        mp.atoms[a].color *= -1;
                    }
                }
            }
        }
        let mut n = 0.0f64;
        for (i, a) in mp.atoms.iter().enumerate() {
            if a.color == -color {
                // `200*ap->x` is int * float, i.e. a single precision product
                let w: f32 = if a.atom_symbol == "C" { 200.0 } else { 1000.0 };
                n += f64::from(w);
                center[0] += f64::from(w * a.x);
                center[1] += f64::from(w * a.y);
            } else if a.color == color {
                let w = 1.0 + 0.01 * i as f64;
                n += w;
                center[0] += w * f64::from(a.x);
                center[1] += w * f64::from(a.y);
            }
        }
        center[0] /= n;
        center[1] /= n;
    }

    let mut bestscore = -1.0f64;
    let mut bpbest = 0;
    let mut projection = 0.0;
    for i in 0..mp.bonds.len() {
        if mp.bonds[i].color == best_color {
            let d = diameter(mp, i, center, &mut projection);
            if d * d + 0.1 * projection * projection > bestscore {
                bestscore = d * d + 0.1 * projection * projection;
                bpbest = i;
            }
        }
    }
    diameter(mp, bpbest, center, &mut projection);
    for a in &mut mp.atoms {
        if a.color < 0 {
            a.color *= -1;
        }
    }
    let first = &mp.atoms[(mp.bonds[bpbest].atoms[0] - 1) as usize];
    let p = [f64::from(first.x), f64::from(first.y)];
    for a in &mut mp.atoms {
        if a.color == color {
            a.x = (f64::from(a.x) - p[0]) as f32;
            a.y = (f64::from(a.y) - p[1]) as f32;
        }
    }
    let second = &mp.atoms[(mp.bonds[bpbest].atoms[1] - 1) as usize];
    let p = [f64::from(second.x), f64::from(second.y)];
    let cosa = p[1] / (p[0] * p[0] + p[1] * p[1]).sqrt();
    let sina = -p[0] / (p[0] * p[0] + p[1] * p[1]).sqrt();
    for a in &mut mp.atoms {
        if a.color == color {
            let q = [f64::from(a.x), f64::from(a.y)];
            a.x = (cosa * q[0] + sina * q[1]) as f32;
            a.y = (-sina * q[0] + cosa * q[1]) as f32;
        }
    }
    if projection < 0.0 {
        for a in &mut mp.atoms {
            if a.color == color {
                a.x = -a.x;
            }
        }
    }
}

/// `LayoutFragments`: places the fragments (atoms with the same colour) next to each other
pub fn layout_fragments(mp: &mut Molecule, ring_size: &[i32], ring_count: &[i32]) {
    struct Frag {
        xll: f64,
        yll: f64,
        xur: f64,
        yur: f64,
        xoffset: f64,
        yoffset: f64,
        color: i32,
        natoms: usize,
    }
    // the C list is built by prepending
    let mut flist: Vec<Frag> = Vec::new();
    for a in &mp.atoms {
        if let Some(f) = flist.iter_mut().find(|f| f.color == a.color) {
            f.natoms += 1;
        } else {
            flist.insert(
                0,
                Frag {
                    xll: 1.0e7,
                    yll: 1.0e7,
                    xur: -1.0e7,
                    yur: -1.0e7,
                    xoffset: 0.0,
                    yoffset: 0.0,
                    color: a.color,
                    natoms: 1,
                },
            );
        }
    }
    for _ in 0..flist.len().saturating_sub(1) {
        for k in 0..flist.len() - 1 {
            if flist[k].natoms < flist[k + 1].natoms {
                let (n, c) = (flist[k + 1].natoms, flist[k + 1].color);
                flist[k + 1].natoms = flist[k].natoms;
                flist[k + 1].color = flist[k].color;
                flist[k].natoms = n;
                flist[k].color = c;
            }
        }
    }
    for k in 0..flist.len() {
        if flist[k].natoms > 1 && flist[k].color & KEEP_POSITION == 0 {
            orient_fragment(mp, flist[k].color, ring_size, ring_count);
        }
    }
    for f in &mut flist {
        for a in &mp.atoms {
            if a.color == f.color {
                f.xll = f.xll.min(f64::from(a.x));
                f.yll = f.yll.min(f64::from(a.y));
                f.xur = f.xur.max(f64::from(a.x));
                f.yur = f.yur.max(f64::from(a.y));
            }
        }
    }
    for f in &mut flist {
        f.xll -= STDBOND / 2.0;
        f.yll -= STDBOND / 2.0;
        f.xur += STDBOND / 2.0;
        f.yur += STDBOND / 2.0;
        f.xoffset = -f.xll;
        f.yoffset = -f.yll;
    }
    let mut xoffset = 0.0;
    for f in &mut flist {
        f.xoffset += xoffset;
        xoffset += f.xur - f.xll;
    }
    for a in &mut mp.atoms {
        if let Some(f) = flist.iter().find(|f| f.color == a.color) {
            a.x = (f64::from(a.x) + f.xoffset) as f32;
            a.y = (f64::from(a.y) + f.yoffset) as f32;
        }
    }
}

const TRANS_MASK: i32 = 0x10;
const CIS_MASK: i32 = 0x20;
const ALL_BOND_TYPES: i32 = 0xF;

/// `Volume` of four points (x, y, z)
fn volume4(t: &[[f64; 3]; 4]) -> f64 {
    let (ax, ay, az) = (t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]);
    let (bx, by, bz) = (t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]);
    let (cx, cy, cz) = (t[3][0] - t[0][0], t[3][1] - t[0][1], t[3][2] - t[0][2]);
    ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx)
}

/// `LayoutAtomStereo`: turns the stereo parities of the atoms into wedged and hashed bonds
pub fn layout_atom_stereo(
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    is_ring_atom: &[i32],
    ring_size: &[i32],
) {
    const CHAIN_ONLY: i32 = 0;
    const ALLOW_RING_BONDS: i32 = 1;
    const USE_ALL: i32 = 2;
    reset_colors(mp);
    for i in 0..mp.atoms.len() {
        if mp.atoms[i].stereo_parity == NONE {
            continue;
        }
        if nbp[i].atoms.len() != 4 {
            continue;
        }
        let (x, y) = (f64::from(mp.atoms[i].x), f64::from(mp.atoms[i].y));
        let mut tetra = [[0.0f64; 3]; 4];
        let mut number = [0usize; 4];
        let mut ligands = [0usize; 4];
        for j in 0..4 {
            number[j] = nbp[i].bonds[j];
            ligands[j] = nbp[i].atoms[j];
            tetra[j][0] = f64::from(mp.atoms[nbp[i].atoms[j]].x);
            tetra[j][1] = f64::from(mp.atoms[nbp[i].atoms[j]].y);
        }
        let mut parity = 1;
        for j in 1..4 {
            for k in (0..j).rev() {
                if ligands[k] > ligands[k + 1] {
                    parity *= -1;
                    ligands.swap(k, k + 1);
                } else {
                    break;
                }
            }
        }
        let mut maxvol = 0.0f64;
        for irun in CHAIN_ONLY..=USE_ALL {
            maxvol = 0.0;
            let mut d = [0i32; 4];
            // whether stereo mark `v` (-1, 0, 1) is allowed on ligand `j`
            let allowed = |mp: &Molecule, j: usize, v: i32| -> bool {
                if v == 0 {
                    return true;
                }
                let bp = &mp.bonds[number[j]];
                if (bp.bond_type != 1) || bp.stereo_symbol != NONE {
                    return false;
                }
                let h = if bp.atoms[0] == i as i32 + 1 {
                    bp.atoms[1]
                } else {
                    bp.atoms[0]
                };
                let h = (h - 1) as usize;
                if irun <= CHAIN_ONLY && is_ring_atom[h] != 0 {
                    return false;
                }
                if irun < ALLOW_RING_BONDS && ring_size[number[j]] != 0 {
                    return false;
                }
                if irun < USE_ALL
                    && (mp.atoms[h].stereo_parity != NONE || mp.atoms[h].color != NONE)
                {
                    return false;
                }
                true
            };
            let zval =
                |j: usize, v: i32| f64::from(v) * 0.1 / nbp[nbp[i].atoms[j]].atoms.len() as f64;
            for i0 in -1..=1 {
                if !allowed(mp, 0, i0) {
                    continue;
                }
                tetra[0][2] = zval(0, i0);
                for i1 in -1..=1 {
                    if !allowed(mp, 1, i1) {
                        continue;
                    }
                    tetra[1][2] = zval(1, i1);
                    for i2 in -1..=1 {
                        if !allowed(mp, 2, i2) {
                            continue;
                        }
                        tetra[2][2] = zval(2, i2);
                        for i3 in -1..=1 {
                            if !allowed(mp, 3, i3) {
                                continue;
                            }
                            tetra[3][2] = zval(3, i3);
                            let cnt = i0 * i0 + i1 * i1 + i2 * i2 + i3 * i3;
                            if !(1..=2).contains(&cnt) {
                                continue;
                            }
                            let v = [i0, i1, i2, i3];
                            if cnt == 1 {
                                // umbrella check
                                let k = (0..4).find(|&q| v[q] != 0).unwrap();
                                let mut nclose = 0;
                                for j in 0..4 {
                                    if j != k
                                        && (tetra[k][0] - x) * (tetra[j][0] - x)
                                            + (tetra[k][1] - y) * (tetra[j][1] - y)
                                            > ((tetra[k][0] - x) * (tetra[k][0] - x)
                                                + (tetra[k][1] - y) * (tetra[k][1] - y))
                                                * 0.1
                                    {
                                        nclose += 1;
                                    }
                                }
                                if nclose < 1 || (nclose < 2 && irun < USE_ALL) {
                                    continue;
                                }
                            }
                            if cnt == 2 {
                                // opposing stereo bonds
                                let nz: Vec<usize> = (0..4).filter(|&q| v[q] != 0).collect();
                                let (j, k) = (nz[0], nz[1]);
                                let dot = (tetra[j][0] - x) * (tetra[k][0] - x)
                                    + (tetra[j][1] - y) * (tetra[k][1] - y);
                                if dot
                                    < -0.3
                                        * ((tetra[j][0] - x) * (tetra[j][0] - x)
                                            + (tetra[j][1] - y) * (tetra[j][1] - y))
                                {
                                    continue;
                                }
                                if dot
                                    < -0.3
                                        * ((tetra[k][0] - x) * (tetra[k][0] - x)
                                            + (tetra[k][1] - y) * (tetra[k][1] - y))
                                {
                                    continue;
                                }
                            }
                            if irun <= ALLOW_RING_BONDS && cnt > 1 {
                                continue;
                            }
                            if maxvol < volume4(&tetra) / f64::from(cnt) {
                                maxvol = volume4(&tetra) / f64::from(cnt);
                                d = v;
                            }
                        }
                    }
                }
            }
            if maxvol > 0.00001 {
                for j in 0..4 {
                    if d[j] == 0 {
                        continue;
                    }
                    if mp.atoms[i].stereo_parity == 1 {
                        d[j] *= -1;
                    }
                    if parity < 0 {
                        d[j] *= -1;
                    }
                    let bp = &mut mp.bonds[number[j]];
                    if bp.atoms[0] != i as i32 + 1 {
                        bp.atoms.swap(0, 1);
                    }
                    bp.stereo_symbol = if d[j] > 0 { UP } else { DOWN };
                }
                mp.atoms[i].stereo_parity = NONE;
                mp.atoms[i].color = NONE + 1;
                break;
            }
        }
        let _ = maxvol;
    }
    reset_colors(mp);
}

/// `ClearFlipFlags`
pub fn clear_flip_flags(mp: &mut Molecule) {
    for b in &mut mp.bonds {
        b.bond_type &= !DONT_FLIP_BOND;
    }
}

/// `ClearDBStereoInSmallRings`
pub fn clear_db_stereo_in_small_rings(mp: &mut Molecule, ring_size: &[i32]) {
    for (i, b) in mp.bonds.iter_mut().enumerate() {
        if ALL_BOND_TYPES & b.bond_type == DOUBLE && ring_size[i] < 8 && ring_size[i] > 0 {
            b.stereo_symbol = NONE;
        }
    }
}

/// the cis/trans product of `LayoutBondStereo` for the reference ligands `l1`, `l2`
fn db_q(mp: &Molecule, ai1: usize, ai2: usize, l1: usize, l2: usize) -> f64 {
    let (a1, a2) = (&mp.atoms[ai1], &mp.atoms[ai2]);
    let p1 = [
        f64::from(mp.atoms[l1].x - a1.x),
        f64::from(mp.atoms[l1].y - a1.y),
    ];
    let p2 = [
        f64::from(mp.atoms[l2].x - a2.x),
        f64::from(mp.atoms[l2].y - a2.y),
    ];
    let r12 = [f64::from(a2.x - a1.x), f64::from(a2.y - a1.y)];
    (p1[0] * r12[1] - p1[1] * r12[0]) * (p2[0] * r12[1] - p2[1] * r12[0])
}

fn db_ok(symbol: i32, q: f64) -> bool {
    ((symbol & TRANS_MASK) != 0 && q < 0.0) || ((symbol & CIS_MASK) != 0 && q > 0.0)
}

/// rotates the double bond atoms around their center (`first`: +30 degrees, else -30 degrees)
fn rotate_bond(
    mp: &mut Molecule,
    ai1: usize,
    ai2: usize,
    pdb1: [f64; 2],
    pdb2: [f64; 2],
    pcenter: [f64; 2],
    first: bool,
) {
    let (m00, m01, m10, m11) = if first {
        (
            (PI / 6.0).cos(),
            (PI / 6.0).sin(),
            -(PI / 6.0).sin(),
            (PI / 6.0).cos(),
        )
    } else {
        (
            (PI / 6.0).cos(),
            -(PI / 6.0).sin(),
            (PI / 6.0).sin(),
            (PI / 6.0).cos(),
        )
    };
    for (ai, p) in [(ai1, pdb1), (ai2, pdb2)] {
        mp.atoms[ai].x =
            (pcenter[0] + m00 * (p[0] - pcenter[0]) + m01 * (p[1] - pcenter[1])) as f32;
        mp.atoms[ai].y =
            (pcenter[1] + m10 * (p[0] - pcenter[0]) + m11 * (p[1] - pcenter[1])) as f32;
    }
}

/// `LayoutBondStereo`: flips double bonds (or marks them swapped) to match their stereo marks
pub fn layout_bond_stereo(mp: &mut Molecule, nbp: &[Neighbourhood], ring_size: &[i32]) {
    let atom_colors: Vec<i32> = mp.atoms.iter().map(|a| a.color).collect();
    for a in &mut mp.atoms {
        a.color = NONE;
    }
    let natoms = mp.atoms.len();
    for i in 0..mp.bonds.len() {
        if ALL_BOND_TYPES & mp.bonds[i].bond_type != DOUBLE {
            continue;
        }
        if ring_size[i] < 8 && ring_size[i] > 0 {
            mp.bonds[i].stereo_symbol = NONE;
            continue;
        }
        let (ai1, ai2) = (
            (mp.bonds[i].atoms[0] - 1) as usize,
            (mp.bonds[i].atoms[1] - 1) as usize,
        );
        let sym = mp.bonds[i].stereo_symbol;
        if sym == NONE || sym == CIS_TRANS_SWAPPED || sym == CIS_TRANS_EITHER {
            continue;
        }
        if nbp[ai1].atoms.len() < 2
            || nbp[ai1].atoms.len() > 3
            || nbp[ai2].atoms.len() < 2
            || nbp[ai2].atoms.len() > 3
        {
            mp.bonds[i].stereo_symbol = NONE;
            continue;
        }
        let mut l1 = natoms;
        for &a in &nbp[ai1].atoms {
            if a != ai2 && a < l1 {
                l1 = a;
            }
        }
        let mut l2 = natoms;
        for &a in &nbp[ai2].atoms {
            if a != ai1 && a < l2 {
                l2 = a;
            }
        }
        let mut q = db_q(mp, ai1, ai2, l1, l2);
        if db_ok(mp.bonds[i].stereo_symbol, q) {
            mp.bonds[i].stereo_symbol = NONE;
            continue;
        }
        if ring_size[i] != 0 && ring_size[i] < 10 {
            mp.bonds[i].stereo_symbol = CIS_TRANS_SWAPPED;
            continue;
        } else if ring_size[i] != 0 {
            let pdb1 = [f64::from(mp.atoms[ai1].x), f64::from(mp.atoms[ai1].y)];
            let pdb2 = [f64::from(mp.atoms[ai2].x), f64::from(mp.atoms[ai2].y)];
            let pcenter = [0.5 * (pdb1[0] + pdb2[0]), 0.5 * (pdb1[1] + pdb2[1])];
            let (mut pdb1_new, mut pdb2_new) = ([0.0f64; 2], [0.0f64; 2]);
            rotate_bond(mp, ai1, ai2, pdb1, pdb2, pcenter, true);
            let mut lastq = 0.0f64;
            q = db_q(mp, ai1, ai2, l1, l2);
            if db_ok(mp.bonds[i].stereo_symbol, q) {
                mp.bonds[i].stereo_symbol = NONE;
                pdb1_new = [f64::from(mp.atoms[ai1].x), f64::from(mp.atoms[ai1].y)];
                pdb2_new = [f64::from(mp.atoms[ai2].x), f64::from(mp.atoms[ai2].y)];
                lastq = q.abs();
            }
            rotate_bond(mp, ai1, ai2, pdb1, pdb2, pcenter, false);
            q = db_q(mp, ai1, ai2, l1, l2);
            if db_ok(mp.bonds[i].stereo_symbol, q) {
                if mp.bonds[i].stereo_symbol != NONE && lastq < q.abs() {
                    pdb1_new = [f64::from(mp.atoms[ai1].x), f64::from(mp.atoms[ai1].y)];
                    pdb2_new = [f64::from(mp.atoms[ai2].x), f64::from(mp.atoms[ai2].y)];
                }
                mp.bonds[i].stereo_symbol = NONE;
            }
            if lastq > 0.0 {
                mp.atoms[ai1].x = pdb1_new[0] as f32;
                mp.atoms[ai1].y = pdb1_new[1] as f32;
                mp.atoms[ai2].x = pdb2_new[0] as f32;
                mp.atoms[ai2].y = pdb2_new[1] as f32;
                continue;
            }
            mp.atoms[ai1].x = pdb1[0] as f32;
            mp.atoms[ai1].y = pdb1[1] as f32;
            mp.atoms[ai2].x = pdb2[0] as f32;
            mp.atoms[ai2].y = pdb2[1] as f32;
            mp.bonds[i].stereo_symbol = CIS_TRANS_SWAPPED;
            continue;
        }

        // a chain bond that has to be flipped
        mp.atoms[ai1].color = -1;
        mp.atoms[ai2].color = 1;
        flood_color(mp, nbp, ai2, 1);
        mp.atoms[ai1].color = 0;
        let p1 = [f64::from(mp.atoms[ai1].x), f64::from(mp.atoms[ai1].y)];
        let p2 = [f64::from(mp.atoms[ai2].x), f64::from(mp.atoms[ai2].y)];
        let r12 = [p2[0] - p1[0], p2[1] - p1[1]];
        for j in 0..natoms {
            if mp.atoms[j].color == 1 {
                mirror_on_bond(mp, j, p1, r12);
            }
        }
        flip_stereo_symbols(mp, 1);
        reset_colors(mp);
        mp.bonds[i].stereo_symbol = NONE;
    }
    for (a, c) in mp.atoms.iter_mut().zip(atom_colors) {
        a.color = c;
    }
}

/// `MakeLandscape`: rotates the molecule by 90 degrees if it is higher than wide
pub fn make_landscape(mp: &mut Molecule) {
    if mp.atoms.is_empty() {
        return;
    }
    let (mut xmin, mut xmax, mut ymin, mut ymax) =
        (mp.atoms[0].x, mp.atoms[0].x, mp.atoms[0].y, mp.atoms[0].y);
    for a in mp.atoms.iter().skip(1) {
        if xmin > a.x {
            xmin = a.x;
        }
        if xmax < a.x {
            xmax = a.x;
        }
        if ymin > a.y {
            ymin = a.y;
        }
        if ymax < a.y {
            ymax = a.y;
        }
    }
    if xmax - xmin < ymax - ymin {
        for a in &mut mp.atoms {
            let h = a.x;
            a.x = a.y;
            a.y = -h;
        }
    }
}

const MAX_COLLISIONS: usize = 4;
const SYM_EPS: f64 = 0.02;
const NO_SYM_EPS: f64 = 0.004;
const STRETCH: f64 = 1.4;

/// difference of two float coordinates; C subtracts in single precision
fn fd(a: f32, b: f32) -> f64 {
    f64::from(a - b)
}

/// `BondLengthSQR`
pub fn bond_length_sqr(mp: &Molecule, bi: usize) -> f64 {
    let b = &mp.bonds[bi];
    let dx = fd(
        mp.atoms[(b.atoms[0] - 1) as usize].x,
        mp.atoms[(b.atoms[1] - 1) as usize].x,
    );
    let dy = fd(
        mp.atoms[(b.atoms[0] - 1) as usize].y,
        mp.atoms[(b.atoms[1] - 1) as usize].y,
    );
    dx * dx + dy * dy
}

/// `AtomBondClash`: clash of atom `ai` with bond `bi` within `eps`; returns the clash and its distance
pub fn atom_bond_clash(mp: &Molecule, ai: usize, bi: usize, eps: f64) -> (bool, f64) {
    let b = &mp.bonds[bi];
    let a0 = &mp.atoms[(b.atoms[0] - 1) as usize];
    let a1 = &mp.atoms[(b.atoms[1] - 1) as usize];
    let ap = &mp.atoms[ai];
    let bx = fd(a1.x, a0.x);
    let by = fd(a1.y, a0.y);
    let ax = fd(ap.x, a0.x);
    let ay = fd(ap.y, a0.y);
    let a2 = ax * ax + ay * ay;
    let b2 = bx * bx + by * by;
    let ab = ax * bx + ay * by;
    let amb2 = (ax - bx) * (ax - bx) + (ay - by) * (ay - by);
    if b2 <= 0.0 {
        return (false, STDBOND);
    }
    let y2 = a2 - ab * ab / b2;
    if y2 >= eps * b2 {
        (false, STDBOND)
    } else if 0.0 <= ab && ab <= b2 {
        (true, y2 / eps)
    } else if a2 <= eps * b2 {
        (true, a2 / eps)
    } else if amb2 <= eps * b2 {
        (true, amb2 / eps)
    } else {
        (false, STDBOND)
    }
}

/// `CountCollisions`: number of collisions and the atom and bond of the closest one
pub fn count_collisions(mp: &Molecule) -> (usize, Option<usize>, Option<usize>) {
    let mut ncol = 0;
    let mut dist = 100.0 * STDBOND;
    let mut apcol = None;
    let mut bpcol = None;
    for (i, ap) in mp.atoms.iter().enumerate() {
        let eps = if ap.charge != 0
            || ap.radical != 0
            || ap.mass_difference != 0
            || ap.atom_symbol != "C"
        {
            SYM_EPS
        } else {
            NO_SYM_EPS
        };
        for (j, bp) in mp.bonds.iter().enumerate() {
            if bp.bond_type & RUBBER_BOND != 0 {
                continue;
            }
            if bp.atoms[0] as usize == i + 1 || bp.atoms[1] as usize == i + 1 {
                continue;
            }
            if ap.color != mp.atoms[(bp.atoms[0] - 1) as usize].color {
                continue;
            }
            if ap.color != mp.atoms[(bp.atoms[1] - 1) as usize].color {
                continue;
            }
            let (clash, disth) = atom_bond_clash(mp, i, j, eps);
            if clash {
                if disth < 100.0 * STDBOND {
                    ncol += 1;
                }
                if disth < dist {
                    dist = disth;
                    apcol = Some(i);
                    bpcol = Some(j);
                }
            }
        }
    }
    (ncol, apcol, bpcol)
}

fn flood_side(mp: &mut Molecule, nbp: &[Neighbourhood], bi: usize, root: usize) {
    let col = mp.atoms[root].color;
    loop {
        let mut changed = false;
        for (i, nb) in nbp.iter().enumerate().take(mp.atoms.len()) {
            if mp.atoms[i].color != col {
                continue;
            }
            for j in 0..nb.atoms.len() {
                if mp.atoms[nb.atoms[j]].color == 0
                    && bi != nb.bonds[j]
                    && mp.bonds[nb.bonds[j]].bond_type & RUBBER_BOND == 0
                {
                    mp.atoms[nb.atoms[j]].color = col;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
}

/// `BondLinksAtoms`: true if bond `bi` is the only link between atoms `a1` and `a2`;
/// the two sides are left colored with `col1` and `col2`
pub fn bond_links_atoms(
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    bi: usize,
    a1: usize,
    col1: i32,
    a2: usize,
    col2: i32,
) -> bool {
    for a in &mut mp.atoms {
        a.color = 0;
    }
    mp.atoms[a1].color = col1;
    flood_side(mp, nbp, bi, a1);
    if mp.atoms[a2].color == col1 {
        return false;
    }
    let b = &mp.bonds[bi];
    let (b0, b1) = ((b.atoms[0] - 1) as usize, (b.atoms[1] - 1) as usize);
    if mp.atoms[b0].color == mp.atoms[b1].color {
        return false;
    }
    mp.atoms[a2].color = col2;
    flood_side(mp, nbp, bi, a2);
    !(mp.atoms[b0].color == 0 || mp.atoms[b1].color == 0)
}

/// `StretchBond`: moves the fragments colored `frag_col1` and `frag_col2` apart along bond `bi`
pub fn stretch_bond(
    mp: &mut Molecule,
    bi: usize,
    mut frag_col1: i32,
    mut frag_col2: i32,
    stretch_factor: f64,
) {
    let b = &mp.bonds[bi];
    let (b0, b1) = ((b.atoms[0] - 1) as usize, (b.atoms[1] - 1) as usize);
    if mp.atoms[b0].color == frag_col2 && mp.atoms[b1].color == frag_col1 {
        std::mem::swap(&mut frag_col1, &mut frag_col2);
    }
    let r12 = [
        fd(mp.atoms[b1].x, mp.atoms[b0].x),
        fd(mp.atoms[b1].y, mp.atoms[b0].y),
    ];
    if r12[0] * r12[0] < 0.0000001 && r12[1] * r12[1] < 0.0000001 {
        return;
    }
    let sf = stretch_factor.sqrt();
    for ap in &mut mp.atoms {
        if ap.color == frag_col1 {
            ap.x = (f64::from(ap.x) - (sf - 1.0) * r12[0]) as f32;
            ap.y = (f64::from(ap.y) - (sf - 1.0) * r12[1]) as f32;
        } else if ap.color == frag_col2 {
            ap.x = (f64::from(ap.x) + (sf - 1.0) * r12[0]) as f32;
            ap.y = (f64::from(ap.y) + (sf - 1.0) * r12[1]) as f32;
        }
    }
}

/// `ImproveCollisions`: removes collisions by stretching chain bonds; `rng` is C's `rand()`
pub fn improve_collisions(
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    is_ring_atom: &[i32],
    ring_count: &[i32],
    rng: &mut CRand,
) {
    let mut atom_colors: Option<Vec<i32>> = None;
    'outer: for _ in 0..MAX_COLLISIONS {
        let (ncol, apcol, bpcol) = count_collisions(mp);
        if ncol == 0 {
            break;
        }
        if atom_colors.is_none() {
            atom_colors = Some(mp.atoms.iter().map(|a| a.color).collect());
        }
        let apcol = apcol.expect("collision atom");
        let bpcol = bpcol.expect("collision bond");
        let mut moved = false;
        'classes: for k in 1..=4 {
            for j in 0..mp.bonds.len() {
                if ring_count[j] > 0 {
                    continue;
                }
                let (i0, i1) = (
                    (mp.bonds[j].atoms[0] - 1) as usize,
                    (mp.bonds[j].atoms[1] - 1) as usize,
                );
                match k {
                    1 => {
                        if is_ring_atom[i0] == 0 || is_ring_atom[i1] == 0 {
                            continue;
                        }
                    }
                    2 => {
                        if (is_ring_atom[i0] != 0) == (is_ring_atom[i1] != 0) {
                            continue;
                        }
                    }
                    3 => {
                        if is_ring_atom[i0] != 0 || is_ring_atom[i1] != 0 {
                            continue;
                        }
                        if nbp[i0].atoms.len() < 3 || nbp[i1].atoms.len() < 3 {
                            continue;
                        }
                    }
                    _ => {
                        if is_ring_atom[i0] != 0 || is_ring_atom[i1] != 0 {
                            continue;
                        }
                        if nbp[i0].atoms.len() >= 3 && nbp[i1].atoms.len() >= 3 {
                            continue;
                        }
                    }
                }
                let a2 = (mp.bonds[bpcol].atoms[0] - 1) as usize;
                if bond_links_atoms(mp, nbp, j, apcol, 1, a2, 2) {
                    stretch_bond(mp, j, 1, 2, STRETCH);
                    let (ncolh, _, _) = count_collisions(mp);
                    if ncol <= ncolh {
                        stretch_bond(mp, j, 1, 2, 1.0 / STRETCH);
                    } else if ncolh > 0 {
                        moved = true;
                        break 'classes;
                    } else {
                        break 'outer;
                    }
                }
            }
        }
        if !moved {
            // a hard case: nudge the atom randomly by up to a tenth of a bond
            let r1 = rng.next_value();
            let r2 = rng.next_value();
            let dx = STDBOND * f64::from(r1 % 100 - 50) / 500.0;
            let dy = STDBOND * f64::from(r2 % 100 - 50) / 500.0;
            let a = &mut mp.atoms[apcol];
            a.x = (f64::from(a.x) + dx) as f32;
            a.y = (f64::from(a.y) + dy) as f32;
        }
    }
    if let Some(cols) = atom_colors {
        for (a, c) in mp.atoms.iter_mut().zip(cols) {
            a.color = c;
        }
    }
}

/// `ScaleByFixedFragments`: scales such that already laid-out bonds average `STDBOND`
pub fn scale_by_fixed_fragments(mp: &mut Molecule) {
    let mut scale = 0.0f64;
    let mut n = 0;
    for b in &mp.bonds {
        let a1 = &mp.atoms[(b.atoms[0] - 1) as usize];
        let a2 = &mp.atoms[(b.atoms[1] - 1) as usize];
        if a1.color == a2.color {
            // SQR of a float difference is evaluated in single precision
            let dx = a1.x - a2.x;
            let dy = a1.y - a2.y;
            scale += f64::from(dx * dx + dy * dy);
            n += 1;
        }
    }
    if scale == 0.0 {
        return;
    }
    scale /= f64::from(n);
    scale = STDBOND / scale.sqrt();
    for a in &mut mp.atoms {
        a.x = (f64::from(a.x) * scale) as f32;
        a.y = (f64::from(a.y) * scale) as f32;
    }
}

fn is_chain_h(
    mp: &Molecule,
    nbp: &[Neighbourhood],
    ring_count: &[i32],
    i: usize,
    j: usize,
) -> bool {
    mp.atoms[nbp[i].atoms[j]].atom_symbol == "H"
        && ring_count[nbp[i].bonds[j]] == 0
        && nbp[nbp[i].atoms[j]].atoms.len() == 1
}

/// `FixLinearBridgeHeads`: deflects bridgeheads with a linear pair of ring bonds
pub fn fix_linear_bridge_heads(mp: &mut Molecule, ring_count: &[i32], nbp: &[Neighbourhood]) {
    for i in 0..mp.atoms.len() {
        let nl = nbp[i].atoms.len();
        if nl != 4 {
            continue;
        }
        let mut nring = 0;
        let mut nh = 0;
        let mut jh = 0;
        for j in 0..nl {
            if is_chain_h(mp, nbp, ring_count, i, j) {
                nh += 1;
                jh = j;
            } else if ring_count[nbp[i].bonds[j]] > 0 {
                nring += 1;
            }
        }
        if nring != 3 || nh != 1 {
            continue;
        }
        let p = [f64::from(mp.atoms[i].x), f64::from(mp.atoms[i].y)];
        let mut i3 = 0;
        let mut found = false;
        'search: for j in 0..nl {
            if j == jh {
                continue;
            }
            let q1 = &mp.atoms[nbp[i].atoms[j]];
            let p1 = [f64::from(q1.x), f64::from(q1.y)];
            let l1 = ((p1[0] - p[0]) * (p1[0] - p[0]) + (p1[1] - p[1]) * (p1[1] - p[1])).sqrt();
            for jj in j + 1..nl {
                if jj == jh {
                    continue;
                }
                let q2 = &mp.atoms[nbp[i].atoms[jj]];
                let p2 = [f64::from(q2.x), f64::from(q2.y)];
                let l2 = ((p2[0] - p[0]) * (p2[0] - p[0]) + (p2[1] - p[1]) * (p2[1] - p[1])).sqrt();
                let mut cosa = (p1[0] - p[0]) * (p2[0] - p[0]) + (p1[1] - p[1]) * (p2[1] - p[1]);
                cosa /= l1 * l2;
                if cosa < -0.99 {
                    found = true;
                    for k in 0..4 {
                        if j != k && jj != k && jh != k {
                            i3 = nbp[i].atoms[k];
                        }
                    }
                    break 'search;
                }
            }
        }
        if !found {
            continue;
        }
        let p3 = [f64::from(mp.atoms[i3].x), f64::from(mp.atoms[i3].y)];
        let ap = &mut mp.atoms[i];
        ap.x = (f64::from(ap.x) + 0.25 * (p3[0] - p[0])) as f32;
        ap.y = (f64::from(ap.y) + 0.25 * (p3[1] - p[1])) as f32;
    }
}

/// `FixInwardHydrogens`: moves hydrogens that point into a ring to the outside
pub fn fix_inward_hydrogens(mp: &mut Molecule, ring_count: &[i32], nbp: &[Neighbourhood]) {
    for i in 0..mp.atoms.len() {
        let nl = nbp[i].atoms.len();
        if nl != 4 {
            continue;
        }
        let mut nring = 0;
        let mut nh = 0;
        let mut nbulk = 100;
        let (mut ih, mut ichain, mut ibest) = (0, 0, 0);
        for j in 0..nl {
            if is_chain_h(mp, nbp, ring_count, i, j) {
                nh += 1;
                ih = nbp[i].atoms[j];
            } else if ring_count[nbp[i].bonds[j]] > 0 {
                nring += 1;
                if nbulk > nbp[nbp[i].atoms[j]].atoms.len() {
                    nbulk = nbp[nbp[i].atoms[j]].atoms.len();
                    ibest = nbp[i].atoms[j];
                }
            } else {
                ichain = nbp[i].atoms[j];
            }
        }
        if nring != 2 || nh != 1 {
            continue;
        }
        let p = [f64::from(mp.atoms[i].x), f64::from(mp.atoms[i].y)];
        let ph = [f64::from(mp.atoms[ih].x), f64::from(mp.atoms[ih].y)];
        let pchain = [f64::from(mp.atoms[ichain].x), f64::from(mp.atoms[ichain].y)];
        let sprod = (ph[0] - p[0]) * (pchain[0] - p[0]) + (ph[1] - p[1]) * (pchain[1] - p[1]);
        if sprod >= 0.0 {
            continue;
        }
        let psum = [
            f64::from(mp.atoms[ibest].x) - p[0] + pchain[0] - p[0],
            f64::from(mp.atoms[ibest].y) - p[1] + pchain[1] - p[1],
        ];
        let lenh = ((ph[0] - p[0]) * (ph[0] - p[0]) + (ph[1] - p[1]) * (ph[1] - p[1])).sqrt();
        let lensum = (psum[0] * psum[0] + psum[1] * psum[1]).sqrt();
        if lensum < 0.25 * lenh {
            continue;
        }
        mp.atoms[ih].x = (p[0] + psum[0] * lenh / lensum) as f32;
        mp.atoms[ih].y = (p[1] + psum[1] * lenh / lensum) as f32;
    }
}

pub const CENTER_DISTANCE: i32 = 1;
pub const TOP_BOTTOM_SEPARATION: i32 = 2;

fn is_plain_h(a: &crate::external::avalon::reaccs::Atom) -> bool {
    a.atom_symbol.as_bytes().first() == Some(&b'H') && a.atom_symbol.len() == 1
}

/// `TopBottomSeparation`: separation of the top and bottom attachments perpendicular to the fragment diameter
pub fn top_bottom_separation(
    mp: &Molecule,
    fragment_color: i32,
    top_color: i32,
    bottom_color: i32,
) -> f64 {
    let mut dmax = 0.0f64;
    let (mut dia1, mut dia2) = (usize::MAX, usize::MAX);
    for i in 0..mp.atoms.len() {
        let a1 = &mp.atoms[i];
        if a1.color != fragment_color || is_plain_h(a1) {
            continue;
        }
        for j in i + 1..mp.atoms.len() {
            let a2 = &mp.atoms[j];
            if a2.color != fragment_color || is_plain_h(a2) {
                continue;
            }
            // SQR of float differences is evaluated in single precision
            let dx = a1.x - a2.x;
            let dy = a1.y - a2.y;
            let d = f64::from(dx * dx + dy * dy);
            if d > dmax {
                dia1 = i;
                dia2 = j;
                dmax = d;
            }
        }
    }
    if dia1 == usize::MAX || dia2 == usize::MAX {
        return 0.0;
    }
    let p = [f64::from(mp.atoms[dia1].x), f64::from(mp.atoms[dia1].y)];
    let mut edia = [
        fd(mp.atoms[dia2].x, mp.atoms[dia1].x),
        fd(mp.atoms[dia2].y, mp.atoms[dia1].y),
    ];
    let d = (edia[0] * edia[0] + edia[1] * edia[1]).sqrt();
    edia[0] /= d;
    edia[1] /= d;
    let eper = [edia[1], -edia[0]];
    let mut total = 0.0;
    let (edges, coords, numbers) = get_colored_graph(mp, fragment_color);
    for b in &mp.bonds {
        let mut i1 = (b.atoms[0] - 1) as usize;
        let mut i2 = (b.atoms[1] - 1) as usize;
        let (c1, c2) = (mp.atoms[i1].color, mp.atoms[i2].color);
        if c1 != fragment_color && c2 != fragment_color {
            continue;
        }
        if c1 == fragment_color && c2 == fragment_color {
            continue;
        }
        if c2 == fragment_color {
            std::mem::swap(&mut i1, &mut i2);
        }
        let c2 = mp.atoms[i2].color;
        if c2 != top_color && c2 != bottom_color {
            continue;
        }
        let seed = numbers[i1];
        let mut p1 = [0.0, 0.0];
        next_substituent_point(&mut p1, &coords, &edges, seed, 0);
        let t = eper[0] * (p1[0] - p[0]) + eper[1] * (p1[1] - p[1]);
        if c2 == top_color {
            total += t;
        } else {
            total -= t;
        }
    }
    total.abs() * dmax
}

/// `AttachmentCenterDistance`: distance of the attachment atoms of fragment `color` from its center
pub fn attachment_center_distance(mp: &Molecule, color: i32, attachment_color: i32) -> f64 {
    let col = |b: &crate::external::avalon::reaccs::Bond, k: usize| {
        mp.atoms[(b.atoms[k] - 1) as usize].color
    };
    // the atom of the bond inside the fragment if the bond leaves it
    let leaving = |b: &crate::external::avalon::reaccs::Bond| -> Option<usize> {
        if col(b, 0) == color && col(b, 1) != color {
            Some((b.atoms[0] - 1) as usize)
        } else if col(b, 1) == color && col(b, 0) != color {
            Some((b.atoms[1] - 1) as usize)
        } else {
            None
        }
    };
    let mut att_dist = 0.0f64;
    for i in 0..mp.bonds.len() {
        let Some(ai1) = leaving(&mp.bonds[i]) else {
            continue;
        };
        for j in i + 1..mp.bonds.len() {
            let Some(ai2) = leaving(&mp.bonds[j]) else {
                continue;
            };
            let dx = mp.atoms[ai1].x - mp.atoms[ai2].x;
            let dy = mp.atoms[ai1].y - mp.atoms[ai2].y;
            let dist = f64::from(dx * dx + dy * dy).sqrt();
            let h = 4.0 * 1.54 - dist;
            let m = if 0.0 < h { 0.0 } else { h };
            att_dist += m * m;
        }
    }
    let (mut xc, mut yc) = (0.0f64, 0.0f64);
    let mut n = 0;
    for a in &mp.atoms {
        if a.color == color {
            xc += f64::from(a.x);
            yc += f64::from(a.y);
            n += 1;
        }
    }
    if n == 0 {
        return 0.0;
    }
    xc /= f64::from(n);
    yc /= f64::from(n);
    let mut result = 0.0;
    for b in &mp.bonds {
        let (a0, a1) = (
            &mp.atoms[(b.atoms[0] - 1) as usize],
            &mp.atoms[(b.atoms[1] - 1) as usize],
        );
        let at = if a0.color == color && a1.color == attachment_color {
            a0
        } else if a1.color == color && a0.color == attachment_color {
            a1
        } else {
            continue;
        };
        result += (f64::from(at.x) - xc) * (f64::from(at.x) - xc);
        result += (f64::from(at.y) - yc) * (f64::from(at.y) - yc);
    }
    result - 10.0 * att_dist
}

/// `FloodInvertColor`: negates the colour of the part of the fragment reachable from `start` without crossing `flip_bond`
pub fn flood_invert_color(
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    start: usize,
    fragment_color: i32,
    flip_bond: usize,
) -> usize {
    if mp.atoms[start].color != fragment_color {
        return 0;
    }
    let mut nflip = 1;
    mp.atoms[start].color *= -1;
    for i in 0..nbp[start].atoms.len() {
        if flip_bond == nbp[start].bonds[i] {
            continue;
        }
        let ai = nbp[start].atoms[i];
        if mp.atoms[ai].color != fragment_color {
            continue;
        }
        nflip += flood_invert_color(mp, nbp, ai, fragment_color, flip_bond);
    }
    nflip
}

/// `ChangeAtomColors`: returns the number of atoms recolored
pub fn change_atom_colors(mp: &mut Molecule, from: i32, to: i32) -> usize {
    let mut n = 0;
    for a in &mut mp.atoms {
        if a.color == from {
            a.color = to;
            n += 1;
        }
    }
    n
}

fn mirror_negative(mp: &mut Molecule, color: i32, p1: [f64; 2], r12: [f64; 2]) {
    for ap in &mut mp.atoms {
        if ap.color == -color {
            let p = [f64::from(ap.x), f64::from(ap.y)];
            let q = ((p[0] - p1[0]) * r12[0] + (p[1] - p1[1]) * r12[1])
                / (r12[0] * r12[0] + r12[1] * r12[1]);
            ap.x = (p1[0] - (p[0] - p1[0]) + 2.0 * r12[0] * q) as f32;
            ap.y = (p1[1] - (p[1] - p1[1]) + 2.0 * r12[1] * q) as f32;
        }
    }
}

fn fragment_score(
    mp: &Molecule,
    improvement_type: i32,
    color: i32,
    col1: i32,
    col2: i32,
) -> Option<f64> {
    let color_strain = -color_strain(mp, color, color);
    match improvement_type {
        CENTER_DISTANCE => Some(attachment_center_distance(mp, color, col1) + color_strain),
        TOP_BOTTOM_SEPARATION => Some(top_bottom_separation(mp, color, col1, col2) + color_strain),
        _ => None,
    }
}

/// `ImproveFragmentByBondFlip`: flips one side of bond `bi` and keeps it if the score improves; returns the score change
pub fn improve_fragment_by_bond_flip(
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    bi: usize,
    improvement_type: i32,
    col1: i32,
    col2: i32,
    forced_undo: bool,
) -> f64 {
    if mp.bonds[bi].bond_type & DONT_FLIP_BOND != 0 {
        return 0.0;
    }
    let (b0, b1) = (
        (mp.bonds[bi].atoms[0] - 1) as usize,
        (mp.bonds[bi].atoms[1] - 1) as usize,
    );
    if mp.atoms[b0].color != mp.atoms[b1].color {
        return 0.0;
    }
    let color = mp.atoms[b0].color;
    let Some(strain1) = fragment_score(mp, improvement_type, color, col1, col2) else {
        return 0.0;
    };
    flood_invert_color(mp, nbp, b0, color, bi);
    if mp.atoms[b1].color == -color {
        change_atom_colors(mp, -color, color);
        return 0.0;
    }
    let p1 = [f64::from(mp.atoms[b0].x), f64::from(mp.atoms[b0].y)];
    let p2 = [f64::from(mp.atoms[b1].x), f64::from(mp.atoms[b1].y)];
    let r12 = [p2[0] - p1[0], p2[1] - p1[1]];
    mirror_negative(mp, color, p1, r12);
    change_atom_colors(mp, -color, color);
    let strain2 = fragment_score(mp, improvement_type, color, col1, col2).unwrap_or(0.0);
    if strain1 > strain2 || forced_undo {
        flood_invert_color(mp, nbp, b0, color, bi);
        mirror_negative(mp, color, p1, r12);
        change_atom_colors(mp, -color, color);
    }
    strain2 - strain1
}

/// `ImproveFragmentByAllBondFlips`: repeatedly applies the best flip of fragment `color`
pub fn improve_fragment_by_all_bond_flips(
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    color: i32,
    improvement_type: i32,
    col1: i32,
    col2: i32,
) -> bool {
    let mut result = false;
    loop {
        let mut best_improvement = 0.0;
        let mut best_bond = None;
        for i in 0..mp.bonds.len() {
            let b = &mp.bonds[i];
            if b.bond_type & DONT_FLIP_BOND != 0
                || mp.atoms[(b.atoms[0] - 1) as usize].color != color
                || mp.atoms[(b.atoms[1] - 1) as usize].color != color
            {
                continue;
            }
            let improvement =
                improve_fragment_by_bond_flip(mp, nbp, i, improvement_type, col1, col2, true);
            if improvement > best_improvement + 0.001 {
                best_improvement = improvement;
                best_bond = Some(i);
            }
        }
        match best_bond {
            Some(bi) => {
                improve_fragment_by_bond_flip(mp, nbp, bi, improvement_type, col1, col2, false);
                result = true;
            }
            None => return result,
        }
    }
}

pub const SHORTCUT_NONE: i32 = 0;
pub const SHORTCUT_START: i32 = 1;
pub const SHORTCUT_END: i32 = 2;
pub const SHORTCUT_INNER: i32 = 4;
pub const SHORTCUT_BRANCH: i32 = 8;
pub const SHORTCUT_RING: i32 = 16;
pub const SHORTCUT_SINGLE: i32 = 32;

fn is_shortcut_atom(a: &crate::external::avalon::reaccs::Atom) -> bool {
    a.atom_symbol == "R" && !a.atext.is_empty()
}

/// `IsShortcutLikeAtom`: R atoms with text, and terminal single-symbol atoms attached to such shortcuts
pub fn is_shortcut_like_atom(mp: &Molecule, i: usize, nbp: &[Neighbourhood]) -> bool {
    let ap = &mp.atoms[i];
    if is_shortcut_atom(ap) {
        // repeat atoms are not rendered as strings
        if ap.atext.bytes().any(|c| c.is_ascii_digit()) {
            return false;
        }
        return nbp[i].atoms.iter().any(|&j| is_shortcut_atom(&mp.atoms[j]));
    }
    if ap.atom_symbol.len() > 1 || nbp[i].atoms.len() > 1 {
        return false;
    }
    // C reads the (zeroed) first neighbour slot of an atom without neighbours
    let n0 = nbp[i].atoms.first().copied().unwrap_or(0);
    let nap = &mp.atoms[n0];
    if is_shortcut_atom(nap) {
        if nap.atext.bytes().any(|c| c.is_ascii_digit()) {
            return false;
        }
        return nbp[n0]
            .atoms
            .iter()
            .any(|&j| is_shortcut_atom(&mp.atoms[j]));
    }
    false
}

const SEQ: i32 = SHORTCUT_INNER | SHORTCUT_START | SHORTCUT_END;

/// `PerceiveShortcutClasses`: fills `shortcut_class` and tells whether a chain of shortcuts was found
pub fn perceive_shortcut_classes(
    mp: &Molecule,
    is_ring_atom: &[i32],
    shortcut_class: &mut [i32],
    nbp: &[Neighbourhood],
    max_line_size: usize,
) -> bool {
    let n = mp.atoms.len();
    let mut used = vec![false; n];
    for i in 0..n {
        shortcut_class[i] = SHORTCUT_NONE;
        if !is_shortcut_like_atom(mp, i, nbp) {
            used[i] = true;
            continue;
        }
        let nshortcuts = nbp[i]
            .atoms
            .iter()
            .filter(|&&j| is_shortcut_like_atom(mp, j, nbp))
            .count();
        if nshortcuts > 2 {
            shortcut_class[i] |= SHORTCUT_BRANCH;
        } else if nshortcuts == 0 {
            shortcut_class[i] |= SHORTCUT_SINGLE;
        }
        if is_ring_atom[i] != 0 && shortcut_class[i] & SHORTCUT_BRANCH != 0 {
            shortcut_class[i] |= SHORTCUT_RING;
        }
        if shortcut_class[i] != SHORTCUT_NONE {
            used[i] = true;
        }
    }
    let mut chain_found = false;
    for i in 0..n {
        if used[i] {
            continue;
        }
        let cands: Vec<usize> = (0..nbp[i].atoms.len())
            .filter(|&j| {
                mp.bonds[nbp[i].bonds[j]].bond_type != NONE
                    && !used[nbp[i].atoms[j]]
                    && is_shortcut_like_atom(mp, nbp[i].atoms[j], nbp)
            })
            .collect();
        if cands.len() == 2 {
            chain_found = true;
            shortcut_class[i] |= SHORTCUT_INNER;
        } else if cands.len() == 1 {
            if nbp[i].atoms[cands[0]] > i {
                shortcut_class[i] |= SHORTCUT_START;
            } else {
                shortcut_class[i] |= SHORTCUT_END;
            }
        } else {
            used[i] = true;
        }
    }
    // split long sequences
    let tlen = |i: usize| mp.atoms[i].atext.len();
    loop {
        let mut next_start: Option<usize> = None;
        for i in 0..n {
            if used[i] {
                continue;
            }
            if shortcut_class[i] & SHORTCUT_START != 0 {
                next_start = Some(i);
                break;
            } else if next_start.is_none() && shortcut_class[i] & SEQ != 0 {
                next_start = Some(i);
            }
        }
        let Some(mut cur) = next_start else {
            break;
        };
        shortcut_class[cur] = SHORTCUT_START;
        let mut line_size = tlen(cur);
        loop {
            used[cur] = true;
            let nb = &nbp[cur];
            let found = (0..nb.atoms.len())
                .find(|&j| !used[nb.atoms[j]] && shortcut_class[nb.atoms[j]] & SEQ != 0);
            let Some(j) = found else {
                if shortcut_class[cur] == SHORTCUT_START {
                    shortcut_class[cur] = SHORTCUT_START;
                } else {
                    shortcut_class[cur] = SHORTCUT_END;
                }
                break;
            };
            let nxt = nb.atoms[j];
            line_size += tlen(cur);
            if max_line_size > 1
                && line_size > max_line_size
                && nb.atoms.len() == 2
                && nbp[nxt].atoms.len() == 2
                && shortcut_class[nxt] & SHORTCUT_END == 0
            {
                line_size = tlen(cur);
                let singleton = !nbp[nxt]
                    .atoms
                    .iter()
                    .any(|&a| !used[a] && shortcut_class[a] & SEQ != 0);
                if !singleton {
                    shortcut_class[cur] = SHORTCUT_END;
                    shortcut_class[nxt] = SHORTCUT_START;
                }
            }
            cur = nxt;
        }
    }
    chain_found
}

/// `LayoutShortcutChains`: arranges chains of acyclic shortcut atoms linearly and rubber-bands their links
pub fn layout_shortcut_chains(mp: &mut Molecule, shortcut_class: &[i32], nbp: &[Neighbourhood]) {
    let n = mp.atoms.len();
    let mut used = vec![false; n];
    let half = |a: &crate::external::avalon::reaccs::Atom| {
        let len = if a.atext.is_empty() {
            a.atom_symbol.len()
        } else {
            a.atext.len()
        };
        0.5 * 1.54 * f64::from(3 + len as i32) / 4.0
    };
    loop {
        let mut start = None;
        let (mut col, mut x, mut y) = (0, 0f32, 0f32);
        for i in 0..n {
            if shortcut_class[i] & (SHORTCUT_START | SHORTCUT_INNER | SHORTCUT_END) == 0 {
                used[i] = true;
            }
            if used[i] {
                continue;
            }
            if shortcut_class[i] & SHORTCUT_START != 0 {
                let ap = &mp.atoms[i];
                col = ap.color;
                x = (f64::from(ap.x) + half(ap)) as f32;
                y = ap.y;
                start = Some(i);
                break;
            }
        }
        let Some(mut i) = start else {
            break;
        };
        used[i] = true;
        loop {
            let nb = &nbp[i];
            let found = (0..nb.atoms.len())
                .find(|&j| !used[nb.atoms[j]] && shortcut_class[nb.atoms[j]] & SEQ != 0);
            let Some(j) = found else {
                break;
            };
            i = nb.atoms[j];
            used[i] = true;
            let h = half(&mp.atoms[i]);
            let ap = &mut mp.atoms[i];
            ap.color = col;
            x = (f64::from(x) + h) as f32;
            ap.x = x;
            ap.y = y;
            x = (f64::from(x) + h) as f32;
            if shortcut_class[i] & SHORTCUT_END != 0 {
                break;
            }
        }
    }
    for b in &mut mp.bonds {
        let (c0, c1) = (
            shortcut_class[(b.atoms[0] - 1) as usize],
            shortcut_class[(b.atoms[1] - 1) as usize],
        );
        if c0 & SEQ != 0 && c1 == SHORTCUT_NONE {
            b.bond_type |= RUBBER_BOND;
        }
        if c1 & SEQ != 0 && c0 == SHORTCUT_NONE {
            b.bond_type |= RUBBER_BOND;
        }
        if c1 & (SHORTCUT_START | SHORTCUT_END) != 0 && c0 & (SHORTCUT_START | SHORTCUT_END) != 0 {
            b.bond_type |= RUBBER_BOND;
        }
    }
}

const FRAME_HEIGHT: f64 = 1.54 * 4.0 / 3.0;
pub const FT_SHORTCUTS: i32 = 1;
pub const FT_ATOMS: i32 = 2;
const MAX_FRAGMENT_NEIGHBOURS: usize = 3;

/// `struct fragment_t`: a fragment (one color) to be placed
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Fragment {
    pub color: i32,
    pub first_index: usize,
    pub xcenter: f64,
    pub ycenter: f64,
    pub radius: f64,
    pub xoffset: f64,
    pub yoffset: f64,
    pub fragment_type: i32,
    pub nsymbols: i32,
    pub nneighbours: usize,
    pub neighbour_cols: [i32; MAX_FRAGMENT_NEIGHBOURS],
    pub nsneighbours: i32,
    pub used: bool,
}

/// the neighbour bookkeeping of `CollectFragmentInfo` for fragment `col` seeing `other` through a bond
fn connect_fragments(frags: &mut [Fragment], col: i32, other: i32) {
    let Some(fi) = frags.iter().position(|f| f.color == col) else {
        return;
    };
    let nn = frags[fi].nneighbours;
    let mut j = 0;
    while j < nn && j < MAX_FRAGMENT_NEIGHBOURS {
        if frags[fi].neighbour_cols[j] == other {
            break;
        }
        j += 1;
    }
    if j < nn {
        return;
    }
    if let Some(h) = frags.iter().find(|f| f.color == other) {
        if h.fragment_type == FT_SHORTCUTS {
            frags[fi].nsneighbours += 1;
        }
    }
    if j == MAX_FRAGMENT_NEIGHBOURS {
        return;
    }
    frags[fi].neighbour_cols[nn] = other;
    frags[fi].nneighbours += 1;
}

/// `CollectFragmentInfo`: shortcut fragments, each followed by the atom fragments linked to it
pub fn collect_fragment_info(
    mp: &Molecule,
    shortcut_class: &[i32],
    nbp: &[Neighbourhood],
) -> Vec<Vec<Fragment>> {
    if mp.atoms.is_empty() {
        return Vec::new();
    }
    // list in creation order; C prepends, so the list proper is the reverse
    let mut frags: Vec<Fragment> = Vec::new();
    for (i, ap) in mp.atoms.iter().enumerate() {
        let shortcut = is_shortcut_like_atom(mp, i, nbp) && shortcut_class[i] != SHORTCUT_SINGLE;
        if let Some(fp) = frags.iter_mut().find(|f| f.color == ap.color) {
            if shortcut {
                fp.fragment_type = FT_SHORTCUTS;
            }
            fp.xcenter += f64::from(ap.x);
            fp.ycenter += f64::from(ap.y);
            fp.nsymbols += 1;
        } else {
            frags.push(Fragment {
                color: ap.color,
                first_index: i,
                xcenter: f64::from(ap.x),
                ycenter: f64::from(ap.y),
                nsymbols: 1,
                fragment_type: if shortcut { FT_SHORTCUTS } else { FT_ATOMS },
                radius: FRAME_HEIGHT * 1.0 / 3.0,
                ..Fragment::default()
            });
        }
    }
    frags.reverse();
    for fp in &mut frags {
        fp.xcenter /= f64::from(fp.nsymbols);
        fp.ycenter /= f64::from(fp.nsymbols);
    }
    for fp in &mut frags {
        fp.radius = FRAME_HEIGHT * 1.0 / 3.0;
        if fp.fragment_type == FT_SHORTCUTS {
            continue;
        }
        for ap in &mp.atoms {
            if fp.color != ap.color {
                continue;
            }
            let dx = fp.xcenter - f64::from(ap.x);
            let dy = fp.ycenter - f64::from(ap.y);
            let d = (dx * dx + dy * dy).sqrt();
            if d > fp.radius {
                fp.radius = d;
            }
        }
    }
    for b in &mp.bonds {
        let col1 = mp.atoms[(b.atoms[0] - 1) as usize].color;
        let col2 = mp.atoms[(b.atoms[1] - 1) as usize].color;
        if col1 == col2 {
            continue;
        }
        connect_fragments(&mut frags, col1, col2);
        connect_fragments(&mut frags, col2, col1);
    }
    for fp in &mut frags {
        for i in 1..fp.nneighbours {
            for j in (0..i).rev() {
                if fp.neighbour_cols[j] > fp.neighbour_cols[j + 1] {
                    fp.neighbour_cols.swap(j, j + 1);
                } else {
                    break;
                }
            }
        }
    }
    // shortcut fragments go to the array, the others are collected in reversed order
    let mut farray: Vec<Vec<Fragment>> = Vec::new();
    let mut atoms_list: Vec<Fragment> = Vec::new();
    for f in frags {
        if f.fragment_type == FT_SHORTCUTS {
            farray.push(vec![f]);
        } else {
            atoms_list.insert(0, f);
        }
    }
    for i in 1..farray.len() {
        for j in (0..i).rev() {
            if farray[j + 1][0].first_index < farray[j][0].first_index {
                farray.swap(j, j + 1);
            } else {
                break;
            }
        }
    }
    for f in atoms_list {
        let target = farray.iter().position(|chain| {
            let c = chain[0].color;
            f.neighbour_cols[..f.nneighbours].contains(&c)
        });
        match target {
            Some(i) => farray[i].insert(1, f),
            None => farray.push(vec![f]),
        }
    }
    farray
}

/// `ColorScore`: how well the connections of the `moving_color` atoms to the `fixed_color` atoms have standard length
pub fn color_score(mp: &Molecule, moving_color: i32, fixed_color: i32) -> f64 {
    let mut result = 0.0;
    for b in &mp.bonds {
        let mut ap1 = &mp.atoms[(b.atoms[0] - 1) as usize];
        let mut ap2 = &mp.atoms[(b.atoms[1] - 1) as usize];
        if ap1.color == ap2.color {
            continue;
        }
        if ap1.color == moving_color && ap2.color == fixed_color {
            std::mem::swap(&mut ap1, &mut ap2);
        } else if !(ap2.color == moving_color && ap1.color == fixed_color) {
            continue;
        }
        let dx = ap1.x - ap2.x;
        let dy = ap1.y - ap2.y;
        let d = f64::from(dx * dx + dy * dy).sqrt();
        result += (d - 1.54) * (d - 1.54);
    }
    result
}

const MAX_ATTACHMENTS: usize = 10;

/// `struct attachment_t`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Attachment {
    pub iattachment: usize,
    pub iligand: usize,
    pub cligand: i32,
    pub nneighbours: i32,
    pub has_symbol: bool,
    pub x: f64,
    pub y: f64,
    pub xlig: f64,
    pub ylig: f64,
    pub xdir: f64,
    pub ydir: f64,
}

/// `CollectAttachments`: the off-fragment attachments of the fragment `color` (at most 10)
pub fn collect_attachments(mp: &Molecule, color: i32, nbp: &[Neighbourhood]) -> Vec<Attachment> {
    let mut atts: Vec<Attachment> = Vec::new();
    for b in &mp.bonds {
        let mut i1 = (b.atoms[0] - 1) as usize;
        let mut i2 = (b.atoms[1] - 1) as usize;
        if mp.atoms[i1].color != color && mp.atoms[i2].color != color {
            continue;
        }
        if mp.atoms[i1].color == mp.atoms[i2].color {
            continue;
        }
        if mp.atoms[i2].color == color {
            std::mem::swap(&mut i1, &mut i2);
        }
        let (a1, a2) = (&mp.atoms[i1], &mp.atoms[i2]);
        atts.push(Attachment {
            iattachment: i1,
            iligand: i2,
            cligand: a2.color,
            has_symbol: a1.atom_symbol != "C",
            x: f64::from(a1.x),
            y: f64::from(a1.y),
            xlig: f64::from(a2.x),
            ylig: f64::from(a2.y),
            ..Attachment::default()
        });
        if atts.len() >= MAX_ATTACHMENTS {
            break;
        }
    }
    for att in &mut atts {
        let nb = &nbp[att.iattachment];
        for &a in &nb.atoms {
            if mp.atoms[a].color == color {
                att.nneighbours += 1;
            } else {
                att.xdir = f64::from(mp.atoms[a].x);
                att.ydir = f64::from(mp.atoms[a].y);
            }
        }
        if att.nneighbours > 0 {
            let (edges, coords, numbers) = get_colored_graph(mp, color);
            let seed = numbers[att.iattachment];
            let mut p = [0.0, 0.0];
            next_substituent_point(&mut p, &coords, &edges, seed, 0);
            att.xdir = p[0] - att.x;
            att.ydir = p[1] - att.y;
        }
    }
    atts
}

/// `struct tfm_t`: rotation, optionally preceded by a flip along the Y axis
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tfm {
    pub cos_alpha: f64,
    pub sin_alpha: f64,
    pub do_flip: bool,
}

/// `TfmPoint`
pub fn tfm_point(x: f64, y: f64, from: [f64; 2], to: [f64; 2], t: &Tfm) -> [f64; 2] {
    let mut x = x - from[0];
    let y = y - from[1];
    if t.do_flip {
        x = -x;
    }
    [
        x * t.cos_alpha - y * t.sin_alpha + to[0],
        x * t.sin_alpha + y * t.cos_alpha + to[1],
    ]
}

/// `ScoreTransformation`
#[allow(clippy::too_many_arguments)]
pub fn score_transformation(
    t: &Tfm,
    atts: &[Attachment],
    from: [f64; 2],
    to: [f64; 2],
    mp: &Molecule,
    color: i32,
    bottom_color: i32,
    top_color: i32,
    target_height: f64,
) -> f64 {
    let na = atts.len() as i32;
    let mut score = 0.0f64;
    let (mut yamin, mut yamax) = (1.0e10f64, -1.0e10f64);
    let (mut xamin, mut xamax) = (1.0e10f64, -1.0e10f64);
    let mut natoms = 0;
    for ap in &mp.atoms {
        if ap.color != color {
            continue;
        }
        natoms += 1;
        let [xnew, ynew] = tfm_point(f64::from(ap.x), f64::from(ap.y), from, to, t);
        yamin = yamin.min_c(ynew);
        yamax = yamax.max_c(ynew);
        xamin = xamin.min_c(xnew);
        xamax = xamax.max_c(xnew);
        for aph in &mp.atoms {
            if aph.color != color && (aph.color == top_color || aph.color == bottom_color) {
                // C multiplies the x difference by (ap->x - ap->x), which is zero; the y term is a plain square.
                // Products of float differences are formed in single precision.
                #[allow(clippy::eq_op)]
                let fx: f32 = (ap.x - aph.x) * (ap.x - ap.x);
                let sq: f32 = (ap.y - aph.y) * (ap.y - aph.y);
                score += 0.1 / (0.1 + f64::from(fx) + f64::from(sq));
                score += 0.2 / (0.1 + f64::from(sq));
            }
        }
    }
    let (mut ymin, mut ymax) = (1.0e10f64, -1.0e10f64);
    for att in atts {
        let l = &mp.atoms[att.iligand];
        let [_, ynew] = tfm_point(f64::from(l.x), f64::from(l.y), from, to, t);
        ymin = ymin.min_c(ynew);
        ymax = ymax.max_c(ynew);
    }
    for att in atts {
        let [xnew, ynew] = tfm_point(att.x, att.y, from, to, t);
        if att.nneighbours <= 1 && att.has_symbol && na > 1 {
            score += 1.0 * ((xnew - att.xlig) * (xnew - att.xlig));
            score += 1.0 * ((ynew - att.ylig) * (ynew - att.ylig));
        } else {
            let [xdir, ydir] = tfm_point(att.xdir, att.ydir, [0.0, 0.0], [0.0, 0.0], t);
            score += (xnew + xdir - att.xlig) * (xnew + xdir - att.xlig);
            score += (ynew + ydir - att.ylig) * (ynew + ydir - att.ylig);
        }
    }
    let nat = f64::from(natoms);
    let naf = f64::from(na);
    if natoms <= 2 {
        0.01 * score + 1.0 * nat * ((yamax - yamin) - 0.00001 * (xamax - xamin))
    } else if na <= 1 {
        1.00 * score
            + 0.0 * naf * (ymax - ymin)
            + 0.0 * nat * ((yamax - yamin) - 0.1 * (xamax - xamin))
    } else {
        let h = target_height - (yamax - yamin);
        1.00 * score
            + 0.0100 * ((ymax - ymin) * (ymax - ymin))
            + 0.0100 * nat * ((yamax - yamin) - (xamax - xamin) + 1.0 * (h * h))
    }
}

/// C's `if (a < m) m = a` / `if (a > m) m = a` as methods that keep `m` for NaN
trait MinMaxC {
    fn min_c(self, a: f64) -> f64;
    fn max_c(self, a: f64) -> f64;
}

impl MinMaxC for f64 {
    fn min_c(self, a: f64) -> f64 {
        if a < self { a } else { self }
    }
    fn max_c(self, a: f64) -> f64 {
        if a > self { a } else { self }
    }
}

const DEGREE_STEP: i32 = 15;

/// `FindBestTransformation`: the best rotation (with or without flip) in 15 degree steps and its score
#[allow(clippy::too_many_arguments)]
pub fn find_best_transformation(
    atts: &[Attachment],
    from: [f64; 2],
    to: [f64; 2],
    mp: &Molecule,
    color: i32,
    bottom_color: i32,
    top_color: i32,
    target_height: f64,
) -> (Tfm, f64) {
    let mut t = Tfm {
        cos_alpha: 1.0,
        sin_alpha: 0.0,
        do_flip: false,
    };
    let mut best = t;
    let mut score = score_transformation(
        &t,
        atts,
        from,
        to,
        mp,
        color,
        bottom_color,
        top_color,
        target_height,
    );
    let angle = |i: i32| f64::from(i * DEGREE_STEP) * 3.14159265358979 / 180.0;
    for i in 1..=360 / DEGREE_STEP {
        t.cos_alpha = angle(i).cos();
        t.sin_alpha = -angle(i).sin();
        let ns = score_transformation(
            &t,
            atts,
            from,
            to,
            mp,
            color,
            bottom_color,
            top_color,
            target_height,
        );
        if ns < score {
            best = t;
            score = ns;
        }
    }
    t.do_flip = true;
    for i in 0..=360 / DEGREE_STEP {
        t.cos_alpha = angle(i).cos();
        t.sin_alpha = -angle(i).sin();
        let ns = score_transformation(
            &t,
            atts,
            from,
            to,
            mp,
            color,
            bottom_color,
            top_color,
            target_height,
        );
        if ns < score {
            best = t;
            score = ns;
        }
    }
    (best, score)
}

/// `AlignFragment`: places the `color` fragment between the shortcut rows `top_y` and `bottom_y`;
/// returns the spacing the rows need
#[allow(clippy::too_many_arguments)]
pub fn align_fragment(
    mp: &mut Molecule,
    nbp: &[Neighbourhood],
    color: i32,
    top_y: f64,
    bottom_y: f64,
    top_color: i32,
    bottom_color: i32,
    mut min_spacing: f64,
    nsneighbours: i32,
) -> f64 {
    let atts = collect_attachments(mp, color, nbp);
    let na = atts.len() as f64;
    let n = mp.atoms.iter().filter(|a| a.color == color).count();
    let (mut xcenter, mut ycenter, mut xnew) = (0.0f64, 0.0f64, 0.0f64);
    for att in &atts {
        xcenter += att.x;
        ycenter += att.y;
        xnew += att.xlig;
    }
    xcenter /= na;
    ycenter /= na;
    xnew /= na;
    // the y coordinate of the sequence atoms is not used: both branches set it from the row positions
    let ynew;
    if nsneighbours == 1 {
        if atts.len() > 1 && n > 2 {
            ynew = bottom_y + (na + 0.1 * (n as f64 - 2.0)).sqrt() * 1.54;
        } else {
            ynew = bottom_y + 1.54;
        }
    } else {
        ynew = (top_y + bottom_y) / 2.0;
    }
    let (tfm, _) = find_best_transformation(
        &atts,
        [xcenter, ycenter],
        [xnew, ynew],
        mp,
        color,
        bottom_color,
        top_color,
        top_y - bottom_y,
    );
    let (mut ymin, mut ymax) = (ynew, ynew);
    for ap in &mut mp.atoms {
        if ap.color != color {
            continue;
        }
        let [x, y] = tfm_point(
            f64::from(ap.x),
            f64::from(ap.y),
            [xcenter, ycenter],
            [xnew, ynew],
            &tfm,
        );
        ap.x = x as f32;
        ap.y = y as f32;
        if ymin > f64::from(ap.y) {
            ymin = f64::from(ap.y);
        }
        if ymax < f64::from(ap.y) {
            ymax = f64::from(ap.y);
        }
    }
    for att in &atts {
        let [_, y] = tfm_point(
            att.x + att.xdir,
            att.y + att.ydir,
            [xcenter, ycenter],
            [xnew, ynew],
            &tfm,
        );
        if ymin > y {
            ymin = y;
        }
        if ymax < y {
            ymax = y;
        }
    }
    let h = FRAME_HEIGHT + ymax - ymin;
    min_spacing = if min_spacing > h { min_spacing } else { h };
    min_spacing
}

/// `FragmentRenderedHeight`
pub fn fragment_rendered_height(mp: &Molecule, color: i32) -> f64 {
    let (mut ymin, mut ymax) = (1e10f64, -1e10f64);
    let mut nmatch = 0;
    for ap in &mp.atoms {
        if ap.color != color {
            continue;
        }
        nmatch += 1;
        if ymin > f64::from(ap.y) {
            ymin = f64::from(ap.y);
        }
        if ymax < f64::from(ap.y) {
            ymax = f64::from(ap.y);
        }
    }
    if nmatch == 0 {
        return 1.5 * FRAME_HEIGHT;
    }
    (ymax - ymin) + 1.5 * FRAME_HEIGHT
}

/// `FragmentHeight`: the extent of the fragment perpendicular to its largest dimension
pub fn fragment_height(mp: &Molecule, color: i32) -> f64 {
    let mut dmax = -1.0f64;
    let mut ex = [0.0f64; 2];
    let mut r1 = [0.0f64; 2];
    for i1 in 0..mp.atoms.len() {
        let a1 = &mp.atoms[i1];
        if a1.color != color {
            continue;
        }
        for a2 in &mp.atoms[i1 + 1..] {
            if a2.color != color {
                continue;
            }
            let dx = a1.x - a2.x;
            let dy = a1.y - a2.y;
            let d = f64::from(dx * dx + dy * dy);
            if d > dmax {
                dmax = d;
                ex = [fd(a2.x, a1.x), fd(a2.y, a1.y)];
                r1 = [f64::from(a1.x), f64::from(a1.y)];
            }
        }
    }
    if dmax <= 0.0 {
        return 1.5 * FRAME_HEIGHT;
    }
    dmax = dmax.sqrt();
    ex[0] /= dmax;
    ex[1] /= dmax;
    let ey = [ex[1], -ex[0]];
    let (mut ymin, mut ymax) = (1.0e10f64, -1.0e10f64);
    for a in &mp.atoms {
        if a.color != color {
            continue;
        }
        let lambda = ey[0] * (r1[0] - f64::from(a.x)) + ey[1] * (r1[1] - f64::from(a.y));
        if lambda < ymin {
            ymin = lambda;
        }
        if lambda > ymax {
            ymax = lambda;
        }
    }
    (ymax - ymin) + 1.5 * FRAME_HEIGHT
}

fn shift_color_y(mp: &mut Molecule, color: i32, dy: f64) {
    for ap in &mut mp.atoms {
        if ap.color == color {
            ap.y = (f64::from(ap.y) - dy) as f32;
        }
    }
}

/// `AlignShortcutsAndFragments`: puts shortcut strings horizontally in rows and places the connected fragments between them
pub fn align_shortcuts_and_fragments(
    mp: &mut Molecule,
    shortcut_class: &[i32],
    nbp: &[Neighbourhood],
) {
    let mut farray = collect_fragment_info(mp, shortcut_class, nbp);
    let nfrag = farray.len();
    let mut first_color = -1;
    for chain in &mut farray {
        for f in chain {
            f.used = false;
        }
    }
    // shortcut fragments go into separate rows
    let mut yoffset = 0.0f64;
    for ifrag in 0..nfrag {
        if ifrag != 0 {
            yoffset += FRAME_HEIGHT;
        }
        if farray[ifrag][0].fragment_type != FT_SHORTCUTS {
            continue;
        }
        let fcolor = farray[ifrag][0].color;
        if first_color == -1 {
            first_color = fcolor;
        }
        farray[ifrag][0].used = true;
        let mut placed = false;
        for i in 0..mp.atoms.len() {
            if mp.atoms[i].color != fcolor {
                continue;
            }
            if shortcut_class[i] & SHORTCUT_START != 0 {
                let col = mp.atoms[i].color;
                let Some(&lig) = nbp[i].atoms.iter().find(|&&a| {
                    col == mp.atoms[a].color
                        && shortcut_class[a] & (SHORTCUT_INNER | SHORTCUT_END) != 0
                }) else {
                    continue;
                };
                let (_, mut coords, numbers) = get_colored_graph(mp, col);
                let p1 = [f64::from(mp.atoms[i].x), f64::from(mp.atoms[i].y)];
                let p2 = [f64::from(mp.atoms[lig].x), f64::from(mp.atoms[lig].y)];
                let p1p = [0.0, -yoffset];
                let p2p = [1.54, -yoffset];
                farray[ifrag][0].ycenter = -yoffset;
                transform_points(&mut coords, p1, p2, p1p, p2p);
                for j in 0..mp.atoms.len() {
                    if mp.atoms[j].color == col {
                        mp.atoms[j].x = coords[numbers[j]][0] as f32;
                        mp.atoms[j].y = coords[numbers[j]][1] as f32;
                    }
                }
                placed = true;
                break;
            }
        }
        if !placed {
            break;
        }
    }

    // move adjacent shortcut rows apart to fit the fragments between them
    for ifrag in 0..nfrag.saturating_sub(1) {
        if farray[ifrag][0].fragment_type != FT_SHORTCUTS {
            continue;
        }
        let gap = farray[ifrag][0].ycenter - farray[ifrag + 1][0].ycenter;
        let mut deltah = gap;
        for k in 1..farray[ifrag].len() {
            let (fcolor, nsn) = (farray[ifrag][k].color, farray[ifrag][k].nsneighbours);
            if nsn == 1 {
                improve_fragment_by_all_bond_flips(
                    mp,
                    nbp,
                    fcolor,
                    CENTER_DISTANCE,
                    farray[ifrag][0].color,
                    -1,
                );
                continue;
            }
            improve_fragment_by_all_bond_flips(
                mp,
                nbp,
                fcolor,
                TOP_BOTTOM_SEPARATION,
                farray[ifrag][0].color,
                farray[ifrag + 1][0].color,
            );
            let new_deltah = fragment_height(mp, fcolor);
            if new_deltah > deltah {
                deltah = new_deltah;
            }
        }
        if deltah > gap {
            deltah -= gap;
            for jfrag in ifrag + 1..nfrag {
                farray[jfrag][0].ycenter -= deltah;
                shift_color_y(mp, farray[jfrag][0].color, deltah);
            }
        }
    }

    // make space before the rows
    for ifrag in 1..nfrag {
        if farray[ifrag][0].fragment_type != FT_SHORTCUTS {
            continue;
        }
        let gap = farray[ifrag - 1][0].ycenter - farray[ifrag][0].ycenter;
        let mut deltah = gap;
        for k in 1..farray[ifrag].len() {
            if farray[ifrag][k].nsneighbours > 1 {
                continue;
            }
            let new_deltah = fragment_height(mp, farray[ifrag][k].color);
            if new_deltah > deltah {
                deltah = new_deltah;
            }
        }
        if deltah > gap {
            deltah -= gap;
            for jfrag in ifrag..nfrag {
                farray[jfrag][0].ycenter -= deltah;
                shift_color_y(mp, farray[jfrag][0].color, deltah);
            }
        }
    }

    // fragments that connect two neighbouring rows
    for ifrag in 0..nfrag {
        if farray[ifrag][0].fragment_type != FT_SHORTCUTS {
            continue;
        }
        for k in 1..farray[ifrag].len() {
            let f = farray[ifrag][k].clone();
            if f.nsneighbours == 1 {
                // placed in the next round
            } else if f.nsneighbours >= 2
                && ifrag + 1 < nfrag
                && f.neighbour_cols[0] == farray[ifrag][0].color
            {
                let deltah = farray[ifrag][0].ycenter - farray[ifrag + 1][0].ycenter;
                align_fragment(
                    mp,
                    nbp,
                    f.color,
                    farray[ifrag][0].ycenter,
                    farray[ifrag + 1][0].ycenter,
                    f.neighbour_cols[0],
                    f.neighbour_cols[1],
                    deltah,
                    f.nsneighbours,
                );
                let mut new_deltah = fragment_rendered_height(mp, f.color);
                if new_deltah > deltah {
                    new_deltah -= deltah;
                    shift_color_y(mp, f.color, new_deltah / 2.0);
                    for jfrag in ifrag + 1..nfrag {
                        farray[jfrag][0].ycenter -= new_deltah;
                        shift_color_y(mp, farray[jfrag][0].color, new_deltah);
                    }
                }
                farray[ifrag][k].used = true;
            }
        }
    }

    // fragments with a single attachment are placed above their row
    for ifrag in 0..nfrag {
        if farray[ifrag][0].fragment_type != FT_SHORTCUTS {
            continue;
        }
        for k in 1..farray[ifrag].len() {
            let f = farray[ifrag][k].clone();
            if f.nsneighbours != 1 {
                continue;
            }
            let deltah = if ifrag > 0 {
                farray[ifrag - 1][0].ycenter - farray[ifrag][0].ycenter
            } else {
                FRAME_HEIGHT
            };
            let top_y = if ifrag > 0 {
                farray[ifrag - 1][0].ycenter
            } else {
                fragment_height(mp, f.color) - FRAME_HEIGHT / 2.0
            };
            align_fragment(
                mp,
                nbp,
                f.color,
                top_y,
                farray[ifrag][0].ycenter,
                0,
                f.neighbour_cols[0],
                deltah,
                f.nsneighbours,
            );
            let mut new_deltah = fragment_rendered_height(mp, f.color);
            if ifrag > 0 && new_deltah > deltah {
                new_deltah -= deltah;
                for kk in 1..farray[ifrag].len() {
                    shift_color_y(mp, farray[ifrag][kk].color, new_deltah);
                }
                for jfrag in ifrag..nfrag {
                    farray[jfrag][0].ycenter -= new_deltah;
                    shift_color_y(mp, farray[jfrag][0].color, new_deltah);
                }
            }
            farray[ifrag][k].used = true;
        }
    }

    // everything laid out gets the color of the first row and stays in place
    for chain in &farray {
        for f in chain {
            if f.used {
                for ap in &mut mp.atoms {
                    if f.color == ap.color {
                        ap.color = first_color;
                    }
                }
            }
        }
    }
    for ap in &mut mp.atoms {
        if ap.color == first_color {
            ap.color |= KEEP_POSITION;
        }
    }
}

/// `LayoutMolecule`: the complete layout of `mp`; fragments pre-laid out by the caller carry a common color
pub fn layout_molecule(mp: &Molecule) -> Molecule {
    let mut result = mp.clone();
    let mut rng = CRand::new(1);
    random_coordinates(&mut result, &mut rng);
    scale_by_fixed_fragments(&mut result);

    let mut nbp = setup_neighbourhood(&result).expect("neighbourhood");
    for nb in &mut nbp {
        let mut pairs: Vec<(usize, usize)> = nb
            .atoms
            .iter()
            .copied()
            .zip(nb.bonds.iter().copied())
            .collect();
        pairs.sort_by_key(|p| p.0);
        nb.atoms = pairs.iter().map(|p| p.0).collect();
        nb.bonds = pairs.iter().map(|p| p.1).collect();
    }
    for nb in &nbp {
        if nb.atoms.len() == 1 {
            let b = &mut result.bonds[nb.bonds[0]];
            if ALL_BOND_TYPES & b.bond_type == DOUBLE {
                b.stereo_symbol = NONE;
            }
        }
    }

    let mut is_ring_atom = vec![0i32; result.atoms.len()];
    let mut ring_count = vec![0i32; result.bonds.len()];
    let mut ring_size = vec![0i32; result.bonds.len()];
    for ring in molecule_rings(&result, &mut rng) {
        for i in 0..result.bonds.len() {
            if ring.bond_set.is_member(i) {
                let (a0, a1) = (
                    result.bonds[i].atoms[0] as usize - 1,
                    result.bonds[i].atoms[1] as usize - 1,
                );
                is_ring_atom[a0] = 1;
                is_ring_atom[a1] = 1;
                if ring_size[i] == 0 || ring_size[i] > ring.cardinality as i32 {
                    ring_size[i] = ring.cardinality as i32;
                }
                result.bonds[i].rsize_flags |= 9.max(ring.cardinality) as u32;
                ring_count[i] += 1;
            }
        }
    }

    let mut shortcut_class = vec![0i32; result.atoms.len()];
    let process_shortcuts =
        perceive_shortcut_classes(&result, &is_ring_atom, &mut shortcut_class, &nbp, 30);
    if process_shortcuts {
        layout_shortcut_chains(&mut result, &shortcut_class, &nbp);
    }
    layout_large_lactams(&mut result, &is_ring_atom, &ring_size, &nbp);
    layout_rings_rng(&mut result, &mut rng);
    flip_ring_carbonyls(&mut result, &is_ring_atom, &ring_size, &nbp);
    link_ring_fragments(&mut result, &is_ring_atom);
    sprout_ring_substituents(&mut result, &is_ring_atom, &ring_size, &nbp);
    layout_chain_atoms(&mut result, &is_ring_atom, &nbp);
    fix_linear_bridge_heads(&mut result, &ring_count, &nbp);
    link_remaining_fragments(&mut result);
    fix_inward_hydrogens(&mut result, &ring_count, &nbp);
    improve_molecule_by_bond_flip(&mut result, &nbp, &ring_size, &is_ring_atom, 0, false);
    clear_flip_flags(&mut result);
    improve_collisions(&mut result, &nbp, &is_ring_atom, &ring_count, &mut rng);
    layout_rubber_fragments(&mut result);
    if process_shortcuts {
        align_shortcuts_and_fragments(&mut result, &shortcut_class, &nbp);
    }
    merge_connected_fragment_colors(&mut result);
    layout_fragments(&mut result, &ring_size, &ring_count);
    layout_bond_stereo(&mut result, &nbp, &ring_size);
    improve_molecule_by_bond_flip(&mut result, &nbp, &ring_size, &is_ring_atom, 0, true);
    clear_flip_flags(&mut result);
    layout_atom_stereo(&mut result, &nbp, &is_ring_atom, &ring_size);
    clear_db_stereo_in_small_rings(&mut result, &ring_size);
    for b in &mut result.bonds {
        b.bond_type &= ALL_BOND_TYPES;
    }
    result
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

    fn dump2(tag: &str, c: &Molecule) -> String {
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
                "b {} {} {} {} {} {}\n",
                i + 1,
                b.atoms[0],
                b.atoms[1],
                b.bond_type,
                b.stereo_symbol,
                b.color
            );
        }
        out
    }

    /// the fourth part of the layout must leave the same coordinates as the C code
    #[test]
    fn layout_part4_matches_c() {
        let text = include_str!("../../../../../tests/data/avalon_layout4.jsonl");
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
            let rs: Vec<i32> = bs.iter().map(|&b| if b != 0 { 6 } else { 0 }).collect();
            sprout_ring_substituents(&mut base, &as_, &rs, &nbp);
            link_ring_fragments(&mut base, &as_);
            let mut out = dump2("base", &base);

            let mut c = base.clone();
            layout_chain_atoms(&mut c, &as_, &nbp);
            out += &dump2("chain", &c);
            let mut c = base.clone();
            layout_rubber_fragments(&mut c);
            out += &dump2("rubber", &c);
            let mut c = base.clone();
            layout_chain_atoms(&mut c, &as_, &nbp);
            layout_rubber_fragments(&mut c);
            link_remaining_fragments(&mut c);
            out += &dump2("linked", &c);
            layout_fragments(&mut c, &rs, &bs);
            out += &dump2("fragments", &c);
            let mut c = base.clone();
            layout_fragments(&mut c, &rs, &bs);
            out += &dump2("fragments0", &c);
            if !base.atoms.is_empty() {
                let mut c = base.clone();
                orient_fragment(&mut c, base.atoms[0].color, &rs, &bs);
                out += &dump2("orient", &c);
            }
            let mut c = base.clone();
            for (i, b) in c.bonds.iter_mut().enumerate() {
                b.color = ((i * 7) % 5) as i32;
            }
            let (best, nb) = best_bond_color(&c);
            out += &format!("best={best}\nnbest={nb}\n");
            for i in 0..base.bonds.len().min(5) {
                let mut c = base.clone();
                let mut proj = -99.0;
                let d = diameter(&mut c, i, [1.0, 2.0], &mut proj);
                out += &format!("diam{i}={} proj={}\n", fmt_gp(d, 17), fmt_gp(proj, 17));
                out += &dump2("diam", &c);
            }
            if out != want {
                bad.push(format!("{input}\nwant:\n{want}\ngot:\n{out}"));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }

    fn dump3(tag: &str, c: &Molecule) -> String {
        let g = |v: f32| fmt_gp(f64::from(v), 9);
        let mut out = format!("{tag}\n");
        for (i, a) in c.atoms.iter().enumerate() {
            out += &format!(
                "a {} {} {} {} {} {} {}\n",
                i + 1,
                g(a.x),
                g(a.y),
                g(a.z),
                a.color,
                a.value as i32,
                a.stereo_parity
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

    /// the stereo part of the layout must leave the same wedges and coordinates as the C code
    #[test]
    fn layout_stereo_matches_c() {
        use crate::external::avalon::smi2mol::{TRUE_DB_STEREO, smi_to_mol};
        let text = include_str!("../../../../../tests/data/avalon_layout5.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let Some(mp) = smi_to_mol(input, TRUE_DB_STEREO).unwrap() else {
                assert_eq!(want, "NULL\n");
                continue;
            };
            let mut base = seeded(7, &mp);
            recolor_molecule(&mut base);
            layout_rings(&mut base);
            let nbp = setup_neighbourhood(&base).unwrap();
            let (as_, bs) = ring_state(&base);
            let mut rs = vec![0i32; base.bonds.len()];
            for i in 0..base.bonds.len() {
                if bs[i] != 0 {
                    let col = base.atoms[(base.bonds[i].atoms[0] - 1) as usize].color;
                    rs[i] = (0..base.atoms.len())
                        .filter(|&j| as_[j] != 0 && base.atoms[j].color == col)
                        .count() as i32;
                }
            }
            sprout_ring_substituents(&mut base, &as_, &rs, &nbp);
            link_ring_fragments(&mut base, &as_);
            layout_chain_atoms(&mut base, &as_, &nbp);
            layout_rubber_fragments(&mut base);
            link_remaining_fragments(&mut base);
            let mut out = dump3("base", &base);

            let mut c = base.clone();
            layout_bond_stereo(&mut c, &nbp, &rs);
            out += &dump3("bondstereo", &c);
            let mut c = base.clone();
            layout_atom_stereo(&mut c, &nbp, &as_, &rs);
            out += &dump3("atomstereo", &c);
            layout_bond_stereo(&mut c, &nbp, &rs);
            clear_db_stereo_in_small_rings(&mut c, &rs);
            clear_flip_flags(&mut c);
            make_landscape(&mut c);
            out += &dump3("all", &c);
            let mut c = base.clone();
            clear_db_stereo_in_small_rings(&mut c, &rs);
            out += &dump3("cleardb", &c);
            if out != want {
                bad.push(format!("{input}\nwant:\n{want}\ngot:\n{out}"));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }

    fn collisions_line(tag: &str, c: &Molecule) -> String {
        let (n, a, b) = count_collisions(c);
        let ix = |o: Option<usize>| o.map_or(-1, |v| v as i64);
        format!("{tag} {n} {} {}\n", ix(a), ix(b))
    }

    /// collisions, stretching, scaling and bridgehead fixes must agree with the C code
    #[test]
    fn layout_collisions_match_c() {
        use crate::external::avalon::smi2mol::{TRUE_DB_STEREO, smi_to_mol};
        let text = include_str!("../../../../../tests/data/avalon_layout6.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let Some(mp) = smi_to_mol(input, TRUE_DB_STEREO).unwrap() else {
                assert_eq!(want, "NULL\n");
                continue;
            };
            let mut base = seeded(7, &mp);
            recolor_molecule(&mut base);
            layout_rings(&mut base);
            let nbp = setup_neighbourhood(&base).unwrap();
            let (as_, bs) = ring_state(&base);
            let mut rs = vec![0i32; base.bonds.len()];
            for i in 0..base.bonds.len() {
                if bs[i] != 0 {
                    let col = base.atoms[(base.bonds[i].atoms[0] - 1) as usize].color;
                    rs[i] = (0..base.atoms.len())
                        .filter(|&j| as_[j] != 0 && base.atoms[j].color == col)
                        .count() as i32;
                }
            }
            sprout_ring_substituents(&mut base, &as_, &rs, &nbp);
            link_ring_fragments(&mut base, &as_);
            layout_chain_atoms(&mut base, &as_, &nbp);
            layout_rubber_fragments(&mut base);
            link_remaining_fragments(&mut base);
            let mut out = dump3("base", &base);
            out += &collisions_line("collisions", &base);
            if !base.bonds.is_empty() {
                out += &format!("bondlen {}\n", fmt_gp(bond_length_sqr(&base, 0), 9));
            } else {
                out += "bondlen 0\n";
            }
            let mut c = base.clone();
            improve_collisions(&mut c, &nbp, &as_, &bs, &mut CRand::new(11));
            out += &dump3("improve", &c);
            let mut c = base.clone();
            scale_by_fixed_fragments(&mut c);
            out += &dump3("scale", &c);
            let mut c = base.clone();
            fix_linear_bridge_heads(&mut c, &bs, &nbp);
            out += &dump3("linear", &c);
            let mut c = base.clone();
            fix_inward_hydrogens(&mut c, &bs, &nbp);
            out += &dump3("inward", &c);
            let mut c = seeded(5, &mp);
            recolor_molecule(&mut c);
            for a in &mut c.atoms {
                a.color = 1;
            }
            out += &collisions_line("rcollisions", &c);
            improve_collisions(&mut c, &nbp, &as_, &bs, &mut CRand::new(11));
            out += &dump3("rimprove", &c);
            scale_by_fixed_fragments(&mut c);
            out += &dump3("rscale", &c);
            if out != want {
                bad.push(format!("{input}\nwant:\n{want}\ngot:\n{out}"));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }

    /// fragment bond flipping and its scores must agree with the C code
    #[test]
    fn layout_fragment_flips_match_c() {
        use crate::external::avalon::smi2mol::{TRUE_DB_STEREO, smi_to_mol};
        let text = include_str!("../../../../../tests/data/avalon_layout7.jsonl");
        let g = |v: f64| fmt_gp(v, 9);
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let Some(mp) = smi_to_mol(input, TRUE_DB_STEREO).unwrap() else {
                assert_eq!(want, "NULL\n");
                continue;
            };
            let mut base = seeded(7, &mp);
            recolor_molecule(&mut base);
            layout_rings(&mut base);
            let nbp = setup_neighbourhood(&base).unwrap();
            let (as_, bs) = ring_state(&base);
            let mut rs = vec![0i32; base.bonds.len()];
            for i in 0..base.bonds.len() {
                if bs[i] != 0 {
                    let col = base.atoms[(base.bonds[i].atoms[0] - 1) as usize].color;
                    rs[i] = (0..base.atoms.len())
                        .filter(|&j| as_[j] != 0 && base.atoms[j].color == col)
                        .count() as i32;
                }
            }
            sprout_ring_substituents(&mut base, &as_, &rs, &nbp);
            link_ring_fragments(&mut base, &as_);
            layout_chain_atoms(&mut base, &as_, &nbp);
            layout_rubber_fragments(&mut base);
            link_remaining_fragments(&mut base);
            let mut out = dump3("base", &base);
            let mut c = base.clone();
            let nring = as_.iter().filter(|&&x| x != 0).count();
            let n = c.atoms.len();
            for (i, a) in c.atoms.iter_mut().enumerate() {
                a.color = if nring > 0 {
                    if as_[i] != 0 {
                        2
                    } else if i % 2 == 1 {
                        4
                    } else {
                        3
                    }
                } else if i < n.div_ceil(2) {
                    2
                } else if i % 2 == 1 {
                    4
                } else {
                    3
                };
            }
            out += &format!("acd {}\n", g(attachment_center_distance(&c, 2, 3)));
            out += &format!("acd2 {}\n", g(attachment_center_distance(&c, 2, 4)));
            out += &format!("tbs {}\n", g(top_bottom_separation(&c, 2, 3, 4)));
            out += &format!("tbs2 {}\n", g(top_bottom_separation(&c, 2, 4, 3)));
            for t in 1..=2 {
                for k in 0..c.bonds.len().min(14) {
                    for fu in 0..2 {
                        let mut d = c.clone();
                        let r = improve_fragment_by_bond_flip(&mut d, &nbp, k, t, 3, 4, fu == 1);
                        out += &format!("flip t{t} b{k} u{fu} {}\n", g(r));
                        if r != 0.0 {
                            out += &dump3("flipped", &d);
                        }
                    }
                }
                let mut d = c.clone();
                let r = improve_fragment_by_all_bond_flips(&mut d, &nbp, 2, t, 3, 4);
                out += &format!("all t{t} {}\n", i32::from(r));
                out += &dump3("allflip", &d);
            }
            if out != want {
                bad.push(format!("{input}\nwant:\n{want}\ngot:\n{out}"));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }

    /// shortcut perception and linear shortcut layout must agree with the C code
    #[test]
    fn layout_shortcuts_match_c() {
        let text = include_str!("../../../../../tests/data/avalon_layout8.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let Some(mp) = mol_str_to_mol(input) else {
                assert_eq!(want, "NULL\n");
                continue;
            };
            let mut c = mp.clone();
            let nbp = setup_neighbourhood(&c).unwrap();
            let (as_, _bs) = ring_state(&c);
            recolor_molecule(&mut c);
            let mut out = String::from("like");
            for i in 0..c.atoms.len() {
                out += &format!(" {}", i32::from(is_shortcut_like_atom(&c, i, &nbp)));
            }
            out += "\n";
            for mls in [0usize, 8, 14] {
                let mut d = c.clone();
                let mut sc = vec![0i32; d.atoms.len()];
                let found = perceive_shortcut_classes(&d, &as_, &mut sc, &nbp, mls);
                out += &format!("classes {mls} {}:", i32::from(found));
                for x in &sc {
                    out += &format!(" {x}");
                }
                out += "\n";
                layout_shortcut_chains(&mut d, &sc, &nbp);
                out += &dump3("chains", &d);
            }
            if out != want {
                bad.push(format!("{}\nwant:\n{want}\ngot:\n{out}", v["name"]));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }

    /// fragment bookkeeping, attachments and alignment must agree with the C code
    #[test]
    fn layout_fragment_alignment_matches_c() {
        let text = include_str!("../../../../../tests/data/avalon_layout9.jsonl");
        let g = |v: f64| fmt_gp(v, 9);
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let Some(mp) = mol_str_to_mol(input) else {
                assert_eq!(want, "NULL\n");
                continue;
            };
            let mut c = mp.clone();
            let nbp = setup_neighbourhood(&c).unwrap();
            let (as_, _bs) = ring_state(&c);
            recolor_molecule(&mut c);
            let mut sc = vec![0i32; c.atoms.len()];
            perceive_shortcut_classes(&c, &as_, &mut sc, &nbp, 14);
            let fa = collect_fragment_info(&c, &sc, &nbp);
            let mut out = format!("nf {}\n", fa.len());
            for (i, chain) in fa.iter().enumerate() {
                out += &format!("chain {i}\n");
                for f in chain {
                    out += &format!(
                        "f {} {} {} {} {} {} {} {} {} {} {},{},{} {} {}\n",
                        f.color,
                        f.first_index,
                        g(f.xcenter),
                        g(f.ycenter),
                        g(f.radius),
                        g(f.xoffset),
                        g(f.yoffset),
                        f.fragment_type,
                        f.nsymbols,
                        f.nneighbours,
                        f.neighbour_cols[0],
                        f.neighbour_cols[1],
                        f.neighbour_cols[2],
                        f.nsneighbours,
                        i32::from(f.used)
                    );
                }
            }
            for k in 0..c.atoms.len() {
                let col = c.atoms[k].color;
                if c.atoms[..k].iter().any(|a| a.color == col) {
                    continue;
                }
                let at = collect_attachments(&c, col, &nbp);
                out += &format!(
                    "color {col} fh {} frh {} na {}\n",
                    g(fragment_height(&c, col)),
                    g(fragment_rendered_height(&c, col)),
                    at.len()
                );
                for a in &at {
                    out += &format!(
                        "att {} {} {} {} {} {} {} {} {} {} {}\n",
                        a.iattachment,
                        a.iligand,
                        a.cligand,
                        a.nneighbours,
                        i32::from(a.has_symbol),
                        g(a.x),
                        g(a.y),
                        g(a.xlig),
                        g(a.ylig),
                        g(a.xdir),
                        g(a.ydir)
                    );
                }
                if let Some(a) = c.atoms.iter().find(|a| a.color != col) {
                    out += &format!("score {} {}\n", a.color, g(color_score(&c, col, a.color)));
                }
                if !at.is_empty() {
                    let topc = at[0].cligand;
                    let botc = at[at.len() - 1].cligand;
                    for ns in 1..=2 {
                        let mut d = c.clone();
                        let sp = align_fragment(&mut d, &nbp, col, 2.0, -2.0, topc, botc, 1.0, ns);
                        out += &format!("align {ns} {}\n", g(sp));
                        out += &dump3("aligned", &d);
                    }
                }
            }
            if out != want {
                bad.push(format!("{}\nwant:\n{want}\ngot:\n{out}", v["name"]));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }

    /// the complete `LayoutMolecule` must agree with the C code
    #[test]
    fn layout_molecule_matches_c() {
        use crate::external::avalon::smi2mol::{TRUE_DB_STEREO, smi_to_mol};
        let text = include_str!("../../../../../tests/data/avalon_layout10.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let mol = if let Some(smi) = input.strip_prefix("SMILES:") {
                smi_to_mol(smi, TRUE_DB_STEREO).unwrap()
            } else {
                mol_str_to_mol(input)
            };
            let Some(mut mp) = mol else {
                assert_eq!(want, "NULL\n");
                continue;
            };
            recolor_molecule(&mut mp);
            let mut out = dump3("input", &mp);
            let r = layout_molecule(&mp);
            out += &dump3("layout", &r);
            if out != want {
                bad.push(format!("{}\nwant:\n{want}\ngot:\n{out}", v["name"]));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }

    /// `SMIToMOL` with `DO_LAYOUT` and coordinate strings must agree with the C code
    #[test]
    fn smi_to_mol_layout_matches_c() {
        use crate::external::avalon::smi2mol::smi_to_mol;
        let text = include_str!("../../../../../tests/data/avalon_layout11.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let flags = v["flags"].as_i64().unwrap() as i32;
            let out = match smi_to_mol(input, flags).unwrap() {
                Some(mp) => dump3("mol", &mp),
                None => "NULL\n".to_string(),
            };
            if out != want {
                bad.push(format!("{flags} {input}\nwant:\n{want}\ngot:\n{out}"));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }
}
