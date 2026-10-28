//! Port of Avalon Toolkit's `canonizer.c`: canonical atom numbering used for canonical SMILES.

use super::graph::CRand;
use super::pattern::{AugmentedAtom, string_to_augmented_atom, transform_augmented_atoms};
use super::perceive::{DOUBLE, Neighbourhood, SINGLE, TRIPLE, ring_state, set_ring_size_flags};
use super::reaccs::{CIS_TRANS_EITHER, CIS_TRANS_SWAPPED, Molecule, NONE, number_of};
use super::utilities::{
    DOWN, UP, atom_symbol_match, compute_implicit_h, flood_fill_fragments, get_num_property,
    reset_colors, strip_molecule,
};

/// `TO_MONOMER` of `smi2mol.h`
pub const TO_MONOMER: i32 = 0x1000;

/// `PropertySplit`: refines the class numbers `pre_num` with `split_prop`; true if they changed
fn property_split(atom_index: &mut [usize], pre_num: &mut [i32], split_prop: &[i64]) -> bool {
    let n = atom_index.len();
    for i in 1..n {
        for j in (0..i).rev() {
            if pre_num[atom_index[j]] == pre_num[atom_index[j + 1]]
                && split_prop[atom_index[j]] > split_prop[atom_index[j + 1]]
            {
                atom_index.swap(j, j + 1);
            } else {
                break;
            }
        }
    }
    let mut new_num = vec![0i32; n];
    let mut last_class = 0;
    new_num[atom_index[0]] = last_class;
    for i in 1..n {
        if pre_num[atom_index[i - 1]] != pre_num[atom_index[i]]
            || split_prop[atom_index[i - 1]] != split_prop[atom_index[i]]
        {
            new_num[atom_index[i]] = i as i32;
            last_class = i as i32;
        } else {
            new_num[atom_index[i]] = last_class;
        }
    }
    let mut result = false;
    for i in 0..n {
        if pre_num[i] != new_num[i] {
            result = true;
        }
        pre_num[i] = new_num[i];
    }
    result
}

/// `ComputeNewGraphProperties`
fn compute_new_graph_properties(
    mp: &Molecule,
    nbp: &[Neighbourhood],
    split_prop: &mut [i64],
    pre_num: &[i32],
) {
    for i in 0..mp.atoms.len() {
        let (mut sum, mut prod) = (0u32, 1u32);
        for j in 0..nbp[i].atoms.len() {
            let bt = mp.bonds[nbp[i].bonds[j]].bond_type;
            let pn = pre_num[nbp[i].atoms[j]];
            sum = sum.wrapping_add(pn.wrapping_mul(bt) as u32);
            prod = prod.wrapping_mul(pn.wrapping_add(bt) as u32);
            prod &= 0xFFFF;
        }
        split_prop[i] = i64::from(0xFFFF & sum.wrapping_mul(0x100).wrapping_add(prod));
    }
}

#[derive(Clone, Copy, Default)]
struct NPoint {
    x: f64,
    y: f64,
    z: f64,
}

fn volume(t: &[NPoint; 4]) -> f64 {
    let (ax, ay, az) = (t[1].x - t[0].x, t[1].y - t[0].y, t[1].z - t[0].z);
    let (bx, by, bz) = (t[2].x - t[0].x, t[2].y - t[0].y, t[2].z - t[0].z);
    let (cx, cy, cz) = (t[3].x - t[0].x, t[3].y - t[0].y, t[3].z - t[0].z);
    ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx)
}

fn sort_neighbourhood_by_rank(nbp: &mut [Neighbourhood], result: &[i32]) {
    for nb in nbp.iter_mut() {
        for j in 1..nb.atoms.len() {
            for jj in (0..j).rev() {
                if result[nb.atoms[jj]] > result[nb.atoms[jj + 1]] {
                    nb.atoms.swap(jj, jj + 1);
                    nb.bonds.swap(jj, jj + 1);
                } else {
                    break;
                }
            }
        }
    }
}

/// the tetrahedron of the ligands of atom `i` with the wedge heights of the bonds that start
/// at `i`; three ligands use the center as fourth point
fn stereo_tetra(mp: &Molecule, nb: &Neighbourhood, i: usize) -> [NPoint; 4] {
    let mut tetra = [NPoint::default(); 4];
    for j in 0..nb.atoms.len() {
        let a = &mp.atoms[nb.atoms[j]];
        tetra[j].x = f64::from(a.x);
        tetra[j].y = f64::from(a.y);
        tetra[j].z = 0.0;
        let bp = &mp.bonds[nb.bonds[j]];
        if (bp.atoms[0] - 1) as usize != i {
            continue;
        }
        if bp.stereo_symbol == UP {
            tetra[j].z = 1.0;
        }
        if bp.stereo_symbol == DOWN {
            tetra[j].z = -1.0;
        }
    }
    if nb.atoms.len() == 3 {
        tetra[3].x = f64::from(mp.atoms[i].x);
        tetra[3].y = f64::from(mp.atoms[i].y);
        tetra[3].z = 0.0;
    }
    tetra
}

/// cross products of the two ligand vectors of a double bond with the bond vector
fn db_q(mp: &Molecule, ai1: usize, l1: usize, ai2: usize, l2: usize, swapped: bool) -> f64 {
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
    let mut q = (p1[0] * r12[1] - p1[1] * r12[0]) * (p2[0] * r12[1] - p2[1] * r12[0]);
    if !swapped {
        q *= -1.0;
    }
    q
}

/// `CanonicalNumbering`: canonical numbers `0..n` for the atoms of `mp`, refined on atom and
/// bond types and on stereochemistry. As in C it orders `nbp` by canonical rank and may clear
/// redundant stereo marks in `mp`.
pub fn canonical_numbering(mp: &mut Molecule, nbp: &mut [Neighbourhood], flags: i32) -> Vec<i32> {
    let n = mp.atoms.len();
    if n == 0 {
        return Vec::new();
    }
    let mut result = vec![0i32; n];
    let mut atom_index: Vec<usize> = (0..n).collect();
    let mut pre_num = vec![0i32; n];
    let (atom_status, bond_status) = ring_state(mp);
    let mut split_prop = vec![0i64; n];

    for i in 0..n {
        let ap = &mp.atoms[i];
        split_prop[i] = nbp[i].atoms.len() as i64 * 0x100;
        let tmp = i64::from(number_of(&ap.atom_symbol));
        if tmp <= 114 {
            split_prop[i] += tmp;
        } else {
            split_prop[i] -= tmp;
        }
        if ap.atom_symbol.starts_with("R#") && (flags & TO_MONOMER) != 0 {
            split_prop[i] += 0x4000;
        }
    }
    property_split(&mut atom_index, &mut pre_num, &split_prop);

    set_ring_size_flags(mp, 12, nbp);
    for i in 0..n {
        split_prop[i] = i64::from(mp.atoms[i].rsize_flags & 0x3FFF);
    }
    property_split(&mut atom_index, &mut pre_num, &split_prop);

    for i in 0..n {
        let ap = &mp.atoms[i];
        let mut v = i64::from(ap.mass_difference);
        if ap.atom_symbol == "D" {
            v = 1;
        }
        if ap.atom_symbol == "T" {
            v = 2;
        }
        v *= 4;
        v += i64::from(ap.charge);
        v *= 4;
        v += i64::from(ap.radical);
        v *= 4;
        v += i64::from(ap.query_h_count);
        if ap.atom_symbol.starts_with("R#") {
            let rgnum = get_num_property(&mp.prop_lines, "M  RGP", i as i32 + 1).unwrap_or(
                if ap.mass_difference > 0 {
                    ap.mass_difference
                } else {
                    0
                },
            );
            v += 0x1000 * i64::from(rgnum);
        }
        split_prop[i] = v;
    }
    property_split(&mut atom_index, &mut pre_num, &split_prop);

    let mut changed;
    'next_cycle: loop {
        loop {
            compute_new_graph_properties(mp, nbp, &mut split_prop, &pre_num);
            changed = property_split(&mut atom_index, &mut pre_num, &split_prop);
            if !changed {
                break;
            }
        }
        for i in 0..n {
            result[atom_index[i]] = i as i32;
        }
        sort_neighbourhood_by_rank(nbp, &result);

        // refine by stereo parity
        for i in 0..n {
            split_prop[i] = 0;
            let nl = nbp[i].atoms.len();
            if !(3..=4).contains(&nl) {
                continue;
            }
            let nequal = (1..nl)
                .filter(|&j| pre_num[nbp[i].atoms[j - 1]] == pre_num[nbp[i].atoms[j]])
                .count();
            if nequal > 0 {
                continue;
            }
            let vol = volume(&stereo_tetra(mp, &nbp[i], i));
            if vol > 0.0001 {
                split_prop[i] = 1;
            }
            if vol < -0.0001 {
                split_prop[i] = -1;
            }
        }
        if property_split(&mut atom_index, &mut pre_num, &split_prop) {
            continue 'next_cycle;
        }

        sort_neighbourhood_by_rank(nbp, &result);

        // split DB neighbours on CIS/TRANS properties
        split_prop.iter_mut().for_each(|p| *p = 0);
        for bp in &mp.bonds {
            if bp.bond_type != DOUBLE || bp.stereo_symbol == CIS_TRANS_EITHER {
                continue;
            }
            let (ai1, ai2) = ((bp.atoms[0] - 1) as usize, (bp.atoms[1] - 1) as usize);
            if nbp[ai1].atoms.len() <= 1 || nbp[ai2].atoms.len() <= 1 {
                continue;
            }
            for &l1 in &nbp[ai1].atoms {
                if l1 == ai2 {
                    continue;
                }
                for &l2 in &nbp[ai2].atoms {
                    if l2 == ai1 {
                        continue;
                    }
                    let q = db_q(mp, ai1, l1, ai2, l2, bp.stereo_symbol == CIS_TRANS_SWAPPED);
                    let s = i64::from(pre_num[l1] + pre_num[l2]);
                    if q < 0.00001 {
                        split_prop[l1] += s;
                        split_prop[l2] += s;
                    } else if q > 0.00001 {
                        split_prop[l1] -= s;
                        split_prop[l2] -= s;
                    }
                }
            }
        }
        if property_split(&mut atom_index, &mut pre_num, &split_prop) {
            continue 'next_cycle;
        }

        // split on DB atoms by CIS/TRANS properties
        split_prop.iter_mut().for_each(|p| *p = 0);
        for bp in &mp.bonds {
            if bp.bond_type != DOUBLE || bp.stereo_symbol == CIS_TRANS_EITHER {
                continue;
            }
            let (ai1, ai2) = ((bp.atoms[0] - 1) as usize, (bp.atoms[1] - 1) as usize);
            let Some(&l1) = nbp[ai1].atoms.iter().find(|&&a| a != ai2) else {
                continue;
            };
            let Some(&l2) = nbp[ai2].atoms.iter().find(|&&a| a != ai1) else {
                continue;
            };
            let q = db_q(mp, ai1, l1, ai2, l2, bp.stereo_symbol == CIS_TRANS_SWAPPED);
            if q < 0.00001 {
                split_prop[ai1] -= 1;
                split_prop[ai2] -= 1;
            } else if q > 0.00001 {
                split_prop[ai1] += 1;
                split_prop[ai2] += 1;
            }
        }
        if property_split(&mut atom_index, &mut pre_num, &split_prop) {
            continue 'next_cycle;
        }

        // refine ring centers by stereo ligand classes
        split_prop.iter_mut().for_each(|p| *p = 0);
        for i in 0..n {
            let nl = nbp[i].atoms.len();
            if !(3..=4).contains(&nl) {
                continue;
            }
            let nequal = (1..nl)
                .filter(|&j| pre_num[nbp[i].atoms[j - 1]] == pre_num[nbp[i].atoms[j]])
                .count();
            if nequal != 1 {
                continue;
            }
            let mut vol = volume(&stereo_tetra(mp, &nbp[i], i));
            if vol > 0.0001 {
                vol = 1.0;
            } else if vol < -0.0001 {
                vol = -1.0;
            } else {
                continue;
            }
            let mut split_in_ring = false;
            for j in 1..nl {
                if pre_num[nbp[i].atoms[j - 1]] == pre_num[nbp[i].atoms[j]] {
                    if atom_status[nbp[i].atoms[j]] != 0 {
                        split_in_ring = true;
                    }
                    if atom_status[nbp[i].atoms[j - 1]] != 0 {
                        split_in_ring = true;
                    }
                    if !split_in_ring {
                        break;
                    }
                    split_prop[nbp[i].atoms[j - 1]] += vol as i64;
                    split_prop[nbp[i].atoms[j]] += -vol as i64;
                    break;
                }
            }
        }
        if property_split(&mut atom_index, &mut pre_num, &split_prop) {
            continue 'next_cycle;
        }

        // clear stereo parity for remaining unresolved redundant centers
        for i in 0..n {
            let nl = nbp[i].atoms.len();
            if !(3..=4).contains(&nl) {
                continue;
            }
            let nrb = nbp[i].bonds.iter().filter(|&&b| bond_status[b] > 0).count();
            let nequal = (1..nl)
                .filter(|&j| pre_num[nbp[i].atoms[j - 1]] == pre_num[nbp[i].atoms[j]])
                .count();
            if nequal == 0 {
                continue;
            }
            let (nrb, nequal) = (nrb as i32, nequal as i32);
            if nrb > 2 && nrb - nequal > 1 {
                continue;
            }
            if nl == 3 && nrb == 2 && nrb - nequal == 1 {
                continue;
            }
            if nl >= 3 && nrb == 3 && nrb - nequal == 1 {
                continue;
            }
            for j in 0..nl {
                let bp = &mut mp.bonds[nbp[i].bonds[j]];
                if (bp.atoms[0] - 1) as usize != i {
                    continue;
                }
                if bp.stereo_symbol == UP {
                    bp.stereo_symbol = NONE;
                }
                if bp.stereo_symbol == DOWN {
                    bp.stereo_symbol = NONE;
                }
            }
        }

        // clear DB stereo with symmetrical substitution on at least one end
        for bi in 0..mp.bonds.len() {
            if mp.bonds[bi].bond_type != DOUBLE {
                continue;
            }
            let mut nequal = 0;
            for end in 0..2 {
                let ai = (mp.bonds[bi].atoms[end] - 1) as usize;
                for j in 1..nbp[ai].atoms.len() {
                    if pre_num[nbp[ai].atoms[j - 1]] == pre_num[nbp[ai].atoms[j]] {
                        nequal += 1;
                    }
                }
            }
            if nequal > 0 && mp.bonds[bi].stereo_symbol == NONE {
                changed = true;
                mp.bonds[bi].stereo_symbol = CIS_TRANS_EITHER;
            }
        }
        if changed {
            continue 'next_cycle;
        }

        // find a remaining tie if any
        let mut min_tie = n as i32;
        let mut min_tie_index: Option<usize> = None;
        for i in 0..n {
            if min_tie <= pre_num[i] {
                continue;
            }
            if (i + 1..n).any(|j| pre_num[i] == pre_num[j]) {
                min_tie = pre_num[i];
                min_tie_index = Some(i);
            }
        }
        if let Some(mti) = min_tie_index {
            for i in 0..n {
                split_prop[i] = i64::from(pre_num[i] == min_tie && i != mti);
            }
            if property_split(&mut atom_index, &mut pre_num, &split_prop) {
                continue 'next_cycle;
            }
        }
        break;
    }
    result
}

/// `ScrambleMolecule`: renumbers `mp` pseudo-randomly (for tests of order independence);
/// `seed` is the C file-level `seed`, which is advanced on each call
pub fn scramble_molecule(mp: &mut Molecule, seed: &mut u32) {
    let mut rng = CRand::new(*seed);
    *seed = seed.wrapping_mul(4).wrapping_add(3);
    let n = mp.atoms.len();
    let mut numbering: Vec<usize> = (0..n).collect();
    for _ in 0..2 * n {
        let i1 = (n as f64 * (f64::from(rng.next_value()) / 2147483647.0)) as usize;
        let i2 = (n as f64 * (f64::from(rng.next_value()) / 2147483647.0)) as usize;
        numbering.swap(i1, i2);
        mp.atoms.swap(numbering[i1], numbering[i2]);
    }
    for b in &mut mp.bonds {
        b.atoms[0] = numbering[(b.atoms[0] - 1) as usize] as i32 + 1;
        b.atoms[1] = numbering[(b.atoms[1] - 1) as usize] as i32 + 1;
    }
}

/// `SplitIonicBonds`: converts bonds between main group metals and hetero atoms to charge
/// separated atoms; true if a bond was deleted
pub fn split_ionic_bonds(mp: &mut Molecule) -> bool {
    let mut deleted = false;
    for i in 0..mp.bonds.len() {
        if mp.bonds[i].bond_type != SINGLE {
            continue;
        }
        let (ai1, ai2) = (
            (mp.bonds[i].atoms[0] - 1) as usize,
            (mp.bonds[i].atoms[1] - 1) as usize,
        );
        let (s1, s2) = (&mp.atoms[ai1].atom_symbol, &mp.atoms[ai2].atom_symbol);
        if atom_symbol_match(s1, "alk,gr2") && atom_symbol_match(s2, "Q") {
            mp.bonds[i].bond_type = -1;
            mp.atoms[ai1].charge += 1;
            mp.atoms[ai2].charge -= 1;
            deleted = true;
        } else if atom_symbol_match(s2, "alk,gr2") && atom_symbol_match(s1, "Q") {
            mp.bonds[i].bond_type = -1;
            mp.atoms[ai2].charge += 1;
            mp.atoms[ai1].charge -= 1;
            deleted = true;
        }
    }
    if !deleted {
        return false;
    }
    let good_atoms = vec![true; mp.atoms.len() + 1];
    let good_bonds: Vec<bool> = mp.bonds.iter().map(|b| b.bond_type >= 0).collect();
    strip_molecule(mp, &good_atoms, &good_bonds);
    true
}

const TFM_STRINGS: [[&str; 2]; 3] = [
    ["N(-C,N,O)(=C,N,O)(=O:1)", "N+1(-C,N,O)(=C,N,O)(-O-1)"],
    ["N(=C,N,O)(#C,N:1)", "N+1(=C,N,O)(=C,N-1)"],
    ["B(-A)(-A)(-A)(-A)", "B-1(-A)(-A)(-A)(-A)"],
];

/// `ForceOctet`: converts bipolar bonds to double bonds, neutralizes charged non-metals and
/// applies the octet-rule transformations; `msgs` receives the transformation messages
pub fn force_octet(mp: &mut Molecule, msgs: &mut Vec<String>) {
    const ZERO_COUNT: i32 = 1;
    let mut h_count = vec![0i32; mp.atoms.len() + 1];
    compute_implicit_h(mp, &mut h_count);
    for (i, ap) in mp.atoms.iter().enumerate() {
        if ap.query_h_count != NONE {
            h_count[i + 1] = ap.query_h_count - ZERO_COUNT;
        }
    }
    for bi in 0..mp.bonds.len() {
        let (mut ai1, mut ai2) = (
            (mp.bonds[bi].atoms[0] - 1) as usize,
            (mp.bonds[bi].atoms[1] - 1) as usize,
        );
        if mp.atoms[ai1].charge <= 0 && mp.atoms[ai2].charge <= 0 {
            continue;
        }
        if mp.atoms[ai2].charge > 0 {
            std::mem::swap(&mut ai1, &mut ai2);
        }
        if mp.atoms[ai2].charge < 0 {
            match mp.bonds[bi].bond_type {
                SINGLE => mp.bonds[bi].bond_type = DOUBLE,
                DOUBLE => mp.bonds[bi].bond_type = TRIPLE,
                _ => continue,
            }
            mp.atoms[ai1].charge -= 1;
            mp.atoms[ai2].charge += 1;
        } else if h_count[ai2 + 1] > 0 && atom_symbol_match(&mp.atoms[ai2].atom_symbol, "Q") {
            mp.atoms[ai1].charge -= 1;
            h_count[ai2 + 1] -= 1;
            mp.bonds[bi].bond_type = DOUBLE;
        }
    }
    for (i, ap) in mp.atoms.iter_mut().enumerate() {
        if ap.charge == 0 || !atom_symbol_match(&ap.atom_symbol, "Q") {
            continue;
        }
        if ap.charge < 0 {
            if ap.query_h_count == NONE {
                h_count[i + 1] -= ap.charge;
                ap.charge = 0;
            }
        } else {
            let mut tmp = ap.charge;
            if h_count[i + 1] < ap.charge {
                tmp = h_count[i + 1];
            }
            h_count[i + 1] -= tmp;
            ap.charge -= tmp;
        }
        if ap.query_h_count != NONE {
            ap.query_h_count = h_count[i + 1] + ZERO_COUNT;
        }
    }
    let table: Vec<[AugmentedAtom; 2]> = TFM_STRINGS
        .iter()
        .map(|p| {
            [
                string_to_augmented_atom(p[0]).expect("rule"),
                string_to_augmented_atom(p[1]).expect("rule"),
            ]
        })
        .collect();
    let _ = transform_augmented_atoms(mp, &table, msgs);
}

/// `RemoveInorganicFragments`: removes fully inorganic fragments if the molecule has carbon
pub fn remove_inorganic_fragments(mp: &mut Molecule) {
    reset_colors(mp);
    let mut ncarbon = 0;
    for a in &mut mp.atoms {
        if a.atom_symbol == "C" {
            ncarbon += 1;
            a.color = 1;
        }
    }
    if ncarbon > 0 {
        flood_fill_fragments(mp);
        let mut good_atoms = vec![false; mp.atoms.len() + 1];
        for (i, a) in mp.atoms.iter().enumerate() {
            good_atoms[i + 1] = a.color == 1;
        }
        let good_bonds: Vec<bool> = mp
            .bonds
            .iter()
            .map(|b| {
                mp.atoms[(b.atoms[0] - 1) as usize].color != NONE
                    || mp.atoms[(b.atoms[1] - 1) as usize].color != NONE
            })
            .collect();
        strip_molecule(mp, &good_atoms, &good_bonds);
    }
    reset_colors(mp);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::avalon::perceive::setup_neighbourhood;
    use crate::external::avalon::reaccsio::{mol_str_to_mol, mol_to_mol_str};

    /// numbering, neighbourhood order and stereo clean-up must be what the C code produces
    #[test]
    fn canonical_numbering_matches_c() {
        let text = include_str!("../../../../../tests/data/avalon_canon.jsonl");
        let mut bad = Vec::new();
        let mut n = 0;
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let mp = mol_str_to_mol(input).expect("molecule");
            let mut out = String::new();
            for flags in [0, TO_MONOMER] {
                let mut c = mp.clone();
                let mut nbp = setup_neighbourhood(&c).unwrap();
                let num = canonical_numbering(&mut c, &mut nbp, flags);
                out += &format!("flags={flags}\nnum:");
                for x in &num {
                    out += &format!(" {x}");
                }
                out += "\nnbp:";
                for nb in &nbp {
                    out += " |";
                    for (a, b) in nb.atoms.iter().zip(&nb.bonds) {
                        out += &format!(" {a}/{b}");
                    }
                }
                out += "\n";
                out += &mol_to_mol_str(&c);
            }
            n += 1;
            if out != want {
                bad.push(format!("{input}\nwant:\n{want}\ngot:\n{out}"));
            }
        }
        assert!(
            bad.is_empty(),
            "{} of {n} differ; first:\n{}",
            bad.len(),
            bad[0]
        );
    }

    fn show(tag: &str, c: &Molecule, msgs: &mut Vec<String>) -> String {
        let mut out = format!("{tag} msgs:\n");
        for m in msgs.drain(..) {
            out += &m;
            out += "\n";
        }
        out += &format!(
            "{tag} name=[{}] init=[{}] prog=[{}] comment=[{}] nlists={} nprops={} chiral={}\n",
            c.name,
            c.user_initials,
            c.program_name,
            c.comment,
            c.n_atom_lists,
            c.n_props,
            c.chiral_flag
        );
        for s in &c.symbol_lists {
            out += &format!("SL {} {} {}\n", s.atom, s.string, s.logic);
        }
        for p in &c.prop_lines {
            out += &format!("PL [{p}]\n");
        }
        out + &mol_to_mol_str(c)
    }

    /// the clean-up functions must leave the same molecule as the C code
    #[test]
    fn strip_functions_match_c() {
        let text = include_str!("../../../../../tests/data/avalon_strip.jsonl");
        let mut bad = Vec::new();
        let mut n = 0;
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let mp = mol_str_to_mol(input).expect("molecule");
            let mut msgs = Vec::new();
            let mut out = String::new();
            let mut c = mp.clone();
            out += &format!("split={}\n", i32::from(split_ionic_bonds(&mut c)));
            out += &show("split", &c, &mut msgs);
            let mut c = mp.clone();
            force_octet(&mut c, &mut msgs);
            out += &show("octet", &c, &mut msgs);
            let mut c = mp.clone();
            remove_inorganic_fragments(&mut c);
            out += &show("inorg", &c, &mut msgs);
            for (tag, k) in [("strip", 0), ("strip2", 1)] {
                let mut c = mp.clone();
                let mut ga = vec![false; c.atoms.len() + 1];
                for i in 0..c.atoms.len() {
                    ga[i + 1] = if k == 0 { i % 4 != 3 } else { i % 3 != 1 };
                }
                let gb: Vec<bool> = c
                    .bonds
                    .iter()
                    .enumerate()
                    .map(|(i, b)| {
                        ga[b.atoms[0] as usize] && ga[b.atoms[1] as usize] && (k == 0 || i % 5 != 4)
                    })
                    .collect();
                strip_molecule(&mut c, &ga, &gb);
                out += &show(tag, &c, &mut msgs);
            }
            let mut seed = 5;
            let mut c = mp.clone();
            scramble_molecule(&mut c, &mut seed);
            out += &show("scr1", &c, &mut msgs);
            let mut c = mp.clone();
            scramble_molecule(&mut c, &mut seed);
            out += &show("scr2", &c, &mut msgs);
            n += 1;
            if out != want {
                bad.push(format!("{input}\nwant:\n{want}\ngot:\n{out}"));
            }
        }
        assert!(
            bad.is_empty(),
            "{} of {n} differ; first:\n{}",
            bad.len(),
            bad[0]
        );
    }
}
