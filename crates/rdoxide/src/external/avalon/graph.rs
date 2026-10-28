//! `graph.c`: basis rings of a graph given as a bond list. Copyright 2010 Novartis Institutes for
//! BioMedical Research Inc.

use super::set::BitSet;

/// `bond_set_node`
#[derive(Debug, Clone)]
pub struct BondSet {
    pub cardinality: usize,
    pub bond_set: BitSet,
}

/// `srand(1)` / `rand()` of glibc (the generator behind the reference results)
#[derive(Debug)]
pub struct CRand {
    f: usize,
    b: usize,
    state: [u32; 31],
}

impl CRand {
    pub fn new(seed: u32) -> Self {
        let mut st = [0u32; 31];
        let mut word = if seed == 0 { 1 } else { seed } as i32 as i64;
        st[0] = word as u32;
        for s in st.iter_mut().skip(1) {
            let hi = word / 127773;
            let lo = word % 127773;
            word = 16807 * lo - 2836 * hi;
            if word < 0 {
                word += 2147483647;
            }
            *s = word as u32;
        }
        let mut c = CRand {
            f: 3,
            b: 0,
            state: st,
        };
        for _ in 0..310 {
            c.next_value();
        }
        c
    }

    pub fn next_value(&mut self) -> i32 {
        self.state[self.f] = self.state[self.f].wrapping_add(self.state[self.b]);
        let result = (self.state[self.f] >> 1) & 0x7fff_ffff;
        self.f = (self.f + 1) % 31;
        self.b = (self.b + 1) % 31;
        result as i32
    }
}

fn other_atom(bond: [u32; 2], atom: usize) -> usize {
    (bond[0] + bond[1]) as usize - atom
}

/// `RingList`: bonds with atom 0 are ignored; the newest ring comes first
pub fn ring_list(bonds: &[[u32; 2]]) -> Vec<BondSet> {
    let nbonds = bonds.len();
    let natoms = bonds.iter().flatten().copied().max().unwrap_or(0) as usize;
    const UNLINKED: i64 = -1;
    let mut color = vec![0i32; natoms + 1];
    let mut link = vec![UNLINKED; natoms + 1];
    let mut result: Vec<BondSet> = Vec::new();
    let mut next_color = 0;
    for b in 0..nbonds {
        let (mut at1, mut at2) = (bonds[b][0] as usize, bonds[b][1] as usize);
        if at1 == 0 || at2 == 0 {
            continue;
        }
        if color[at1] == 0 && color[at2] == 0 {
            next_color += 1;
            link[at2] = b as i64;
            color[at1] = next_color;
            color[at2] = next_color;
        } else if color[at1] == 0 {
            color[at1] = color[at2];
            link[at1] = b as i64;
        } else if color[at2] == 0 {
            color[at2] = color[at1];
            link[at2] = b as i64;
        } else if color[at1] != color[at2] {
            let (new_color, old_color) = (color[at1], color[at2]);
            for c in color.iter_mut() {
                if *c == old_color {
                    *c = new_color;
                }
            }
            let mut tmp_link = link[at2];
            link[at2] = b as i64;
            while tmp_link != UNLINKED {
                at2 = other_atom(bonds[tmp_link as usize], at2);
                std::mem::swap(&mut tmp_link, &mut link[at2]);
            }
        } else {
            let mut trace = at1;
            let mut level1 = 0;
            while link[trace] != UNLINKED {
                trace = other_atom(bonds[link[trace] as usize], trace);
                level1 += 1;
            }
            trace = at2;
            let mut level2 = 0;
            while link[trace] != UNLINKED {
                trace = other_atom(bonds[link[trace] as usize], trace);
                level2 += 1;
            }
            if level1 > level2 {
                std::mem::swap(&mut level1, &mut level2);
                std::mem::swap(&mut at1, &mut at2);
            }
            let mut p = BondSet {
                cardinality: 1,
                bond_set: BitSet::new(nbonds),
            };
            p.bond_set.put(b);
            for _ in 0..level2 - level1 {
                p.bond_set.put(link[at2] as usize);
                p.cardinality += 1;
                at2 = other_atom(bonds[link[at2] as usize], at2);
            }
            while at1 != at2 {
                p.bond_set.put(link[at1] as usize);
                p.cardinality += 1;
                at1 = other_atom(bonds[link[at1] as usize], at1);
                p.bond_set.put(link[at2] as usize);
                p.cardinality += 1;
                at2 = other_atom(bonds[link[at2] as usize], at2);
            }
            result.insert(0, p);
        }
    }
    result
}

/// `SortRings`: the C code's exchange sort (not stable), descending by... ascending cardinality
pub fn sort_rings(list: &mut [BondSet]) {
    for i in 0..list.len() {
        for j in i + 1..list.len() {
            if list[i].cardinality > list[j].cardinality {
                let (a, b) = list.split_at_mut(j);
                std::mem::swap(&mut a[i], &mut b[0]);
            }
        }
    }
}

/// `CombineRings`
pub fn combine_rings(list: &mut [BondSet]) {
    combine_rings_rng(list, &mut CRand::new(1));
}

/// `CombineRings` that leaves the C `rand()` state in `rng`: the function does `srand(1)` unless the ring list is empty
pub fn combine_rings_rng(list: &mut [BondSet], rng: &mut CRand) {
    if list.is_empty() {
        return;
    }
    *rng = CRand::new(1);
    let mut ntoggle = 0;
    let mut set = BitSet::new(list[0].bond_set.max_member());
    loop {
        let mut changed = false;
        sort_rings(list);
        for i in 0..list.len() {
            for j in i + 1..list.len() {
                if list[i].bond_set.intersection_is_empty(&list[j].bond_set) {
                    continue;
                }
                set.copy_from(&list[i].bond_set);
                set.xor_with(&list[j].bond_set);
                let size = set.cardinality();
                let (c1, c2) = (list[i].cardinality, list[j].cardinality);
                if size > 0 && (size <= c1 || size <= c2) {
                    let k = if c1 > c2 { i } else { j };
                    let ck = list[k].cardinality;
                    if ck > size || (rng.next_value() / 10) % 2 != 0 {
                        if ck == size {
                            ntoggle += 1;
                        } else {
                            ntoggle = 0;
                        }
                        changed = true;
                        list[k].cardinality = size;
                        std::mem::swap(&mut list[k].bond_set, &mut set);
                    }
                }
            }
        }
        if ntoggle > 4 {
            changed = false;
        }
        if !changed {
            break;
        }
    }
}

/// `ProperRingPairs`
pub fn proper_ring_pairs(base: &[BondSet], maxnode: usize, bonds: &[[u32; 2]]) -> Vec<BondSet> {
    let mut result = Vec::new();
    if base.is_empty() {
        return result;
    }
    let mut set = BitSet::new(base[0].bond_set.max_member());
    for i in 0..base.len() {
        for j in i + 1..base.len() {
            if base[i].bond_set.intersection_is_empty(&base[j].bond_set) {
                continue;
            }
            set.copy_from(&base[i].bond_set);
            set.intersect_with(&base[j].bond_set);
            let mut touched = vec![0; maxnode + 1];
            for b in 0..set.max_member() {
                if set.is_member(b) {
                    touched[bonds[b][0] as usize] += 1;
                    touched[bonds[b][1] as usize] += 1;
                }
            }
            if touched.iter().filter(|&&t| t == 1).count() == 2 {
                let mut ph = BondSet {
                    cardinality: 0,
                    bond_set: BitSet::new(base[i].bond_set.max_member()),
                };
                ph.bond_set.copy_from(&base[i].bond_set);
                ph.bond_set.xor_with(&base[j].bond_set);
                ph.cardinality = ph.bond_set.cardinality();
                result.insert(0, ph);
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn glibc_rand_sequence() {
        let mut r = CRand::new(1);
        assert_eq!(
            [r.next_value(), r.next_value(), r.next_value()],
            [1804289383, 846930886, 1681692777]
        );
    }

    #[test]
    fn two_fused_rings() {
        // naphthalene-like: bonds of two six rings sharing one bond
        let bonds: Vec<[u32; 2]> = vec![
            [1, 2],
            [2, 3],
            [3, 4],
            [4, 5],
            [5, 6],
            [6, 1],
            [5, 7],
            [7, 8],
            [8, 9],
            [9, 10],
            [10, 6],
        ];
        let mut rings = ring_list(&bonds);
        assert_eq!(rings.len(), 2);
        combine_rings(&mut rings);
        let mut sizes: Vec<_> = rings.iter().map(|r| r.cardinality).collect();
        sizes.sort();
        assert_eq!(sizes, [6, 6]);
    }
}
