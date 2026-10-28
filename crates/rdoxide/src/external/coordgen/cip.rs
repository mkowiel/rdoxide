//! The stereochemistry code of `sketcherMinimizerAtom` and `sketcherMinimizerBond`: the CIP ranking used to decide the
//! R/S and Z/E descriptors of the sketcher (`CIPPriority` and its helpers), reading and writing the wedge based
//! stereochemistry, and the double bond stereo calls that depend on it.

use std::collections::{HashMap, VecDeque};

use super::maths::signed_angle;
use super::model::{AtomId, BondId, BondStereo, Chirality, Model, MolId};

/// `sketcherMinimizerAtomPriority`
#[derive(Debug, Clone, Copy)]
pub struct AtomPriority {
    pub a: AtomId,
    pub priority: u32,
}

/// the maps that the `CIPAtom`s of one side of a comparison share
#[derive(Debug, Default, Clone)]
struct CipState {
    scores: HashMap<AtomId, i32>,
    medals: HashMap<AtomId, Vec<i32>>,
    visited: HashMap<AtomId, i32>,
}

impl CipState {
    fn score(&self, a: AtomId) -> i32 {
        self.scores.get(&a).copied().unwrap_or(0)
    }
    fn medals_of(&self, a: AtomId) -> Vec<i32> {
        self.medals.get(&a).cloned().unwrap_or_default()
    }
}

/// `CIPAtom`: the atoms of one level of the hierarchical digraph below one parent (`None`: dummy atom)
#[derive(Debug, Clone)]
struct CipAtom {
    these_atoms: Vec<(i32, Option<AtomId>)>,
    parent: AtomId,
    all_parents: Vec<AtomId>,
    /// which of the two states the maps belong to
    side: usize,
}

/// `CIPAtom::operator<` (and `isBetter` when `friends` is given): has `lhs` priority over `rhs`
fn cip_less(
    model: &Model,
    st: &[CipState; 2],
    lhs: &CipAtom,
    rhs: &CipAtom,
    friends: Option<&HashMap<AtomId, u32>>,
) -> bool {
    debug_assert_eq!(lhs.all_parents.len(), rhs.all_parents.len());
    let (sl, sr) = (&st[lhs.side], &st[rhs.side]);
    for i in 0..lhs.all_parents.len() {
        let (pl, pr) = (lhs.all_parents[i], rhs.all_parents[i]);
        if let Some(m) = friends {
            let (ml, mr) = (
                m.get(&pl).copied().unwrap_or(0),
                m.get(&pr).copied().unwrap_or(0),
            );
            if ml > mr {
                return true;
            }
            if ml < mr {
                return false;
            }
        }
        let (zl, zr) = (model.atoms[pl].atomic_number, model.atoms[pr].atomic_number);
        if zl > zr {
            return true;
        }
        if zl < zr {
            return false;
        }
        if sl.score(pl) < sr.score(pr) {
            return true;
        }
        if sl.score(pl) > sr.score(pr) {
            return false;
        }
        let (meds, meds2) = (sl.medals_of(pl), sr.medals_of(pr));
        let s = meds.len().min(meds2.len());
        for mm in 0..s {
            if meds[mm] > meds2[mm] {
                return true;
            }
            if meds[mm] < meds2[mm] {
                return false;
            }
        }
        if meds.len() > meds2.len() {
            return true;
        }
        if meds2.len() > meds.len() {
            return false;
        }
    }
    let siz = lhs.these_atoms.len().min(rhs.these_atoms.len());
    for i in 0..siz {
        if lhs.these_atoms[i].0 > rhs.these_atoms[i].0 {
            return true;
        }
        if lhs.these_atoms[i].0 < rhs.these_atoms[i].0 {
            return false;
        }
    }
    lhs.these_atoms.len() > rhs.these_atoms.len()
}

/// `CIPAtom::operator==`
fn cip_equal(model: &Model, st: &[CipState; 2], lhs: &CipAtom, rhs: &CipAtom) -> bool {
    debug_assert_eq!(lhs.all_parents.len(), rhs.all_parents.len());
    for i in 0..lhs.all_parents.len() {
        let (pl, pr) = (lhs.all_parents[i], rhs.all_parents[i]);
        if model.atoms[pl].atomic_number != model.atoms[pr].atomic_number {
            return false;
        }
        if st[lhs.side].score(pl) != st[rhs.side].score(pr) {
            return false;
        }
    }
    lhs.these_atoms.len() == rhs.these_atoms.len()
        && lhs
            .these_atoms
            .iter()
            .zip(&rhs.these_atoms)
            .all(|(a, b)| a.0 == b.0)
}

/// `std::stable_sort` with `operator<`
fn stable_sort_cip(model: &Model, st: &[CipState; 2], v: &mut [CipAtom]) {
    v.sort_by(|a, b| {
        if cip_less(model, st, a, b, None) {
            std::cmp::Ordering::Less
        } else if cip_less(model, st, b, a, None) {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    });
}

impl Model {
    /// a temporary atom that is not part of any molecule list (`dummyH`: hydrogen, `dummyLP`: lone pair)
    fn push_dummy(&mut self, atomic_number: i32, mol: Option<MolId>) -> AtomId {
        self.atoms.push(super::model::Atom {
            atomic_number,
            molecule: mol,
            ..super::model::Atom::default()
        });
        self.atoms.len() - 1
    }

    /// `matchCIPSequence`: whether the two sequences are an even permutation of each other (pads the shorter one with 3)
    pub fn match_cip_sequence(v1: &mut Vec<i32>, v2: &mut Vec<i32>) -> bool {
        if v1.len() < v2.len() {
            v1.push(3);
        } else if v2.len() < v1.len() {
            v2.push(3);
        }
        let out_of_place = v1.iter().zip(v2.iter()).filter(|(a, b)| a != b).count();
        if out_of_place == 2 {
            return false;
        } else if out_of_place == 4 {
            let n1 = v1[0];
            let index2 = v2.iter().position(|&x| x == n1).unwrap_or(0);
            if v1[index2] != v2[0] {
                return false;
            }
        }
        true
    }

    /// `sketcherMinimizer::sameRing(at1, at2, at3)`: the smallest ring shared by the three atoms
    pub fn same_ring3(&self, at1: AtomId, at2: AtomId, at3: AtomId) -> Option<usize> {
        let (r1, r2, r3) = (
            &self.atoms[at1].rings,
            &self.atoms[at2].rings,
            &self.atoms[at3].rings,
        );
        if r1.is_empty() || r2.is_empty() || r3.is_empty() {
            return None;
        }
        let mut r: Option<usize> = None;
        for pass in 0..2 {
            for &ring in r1 {
                if pass == 0 && self.rings[ring].is_macrocycle() {
                    continue;
                }
                for &ring2 in r2 {
                    if ring != ring2 {
                        continue;
                    }
                    for &ring3 in r3 {
                        if ring3 == ring2 {
                            match r {
                                None => r = Some(ring2),
                                Some(cur)
                                    if self.rings[ring2].atoms.len()
                                        < self.rings[cur].atoms.len() =>
                                {
                                    r = Some(ring2)
                                }
                                _ => {}
                            }
                        }
                    }
                }
            }
        }
        r
    }

    /// `CIPPriority`: the atom of higher CIP priority as seen from `center`, `None` if they are equivalent
    pub fn cip_priority(&self, at1: AtomId, at2: AtomId, center: AtomId) -> Option<AtomId> {
        let (z1, z2) = (self.atoms[at1].atomic_number, self.atoms[at2].atomic_number);
        if z1 > z2 {
            return Some(at1);
        } else if z2 > z1 {
            return Some(at2);
        }
        let c = &self.atoms[center];
        if c.has_ring_chirality
            && !c.ignore_ring_chirality
            && self.same_ring3(center, at1, at2).is_some()
        {
            return if self.atoms[at1].general_use_n > self.atoms[at2].general_use_n {
                Some(at1)
            } else {
                Some(at2)
            };
        }
        let mut st: [CipState; 2] = [CipState::default(), CipState::default()];
        for (side, at) in [at1, at2].into_iter().enumerate() {
            st[side].visited.insert(center, 1);
            st[side].visited.insert(at, 2);
        }
        let mk = |at: AtomId, side: usize| CipAtom {
            these_atoms: vec![(self.atoms[at].atomic_number, Some(at))],
            parent: center,
            all_parents: vec![center, at],
            side,
        };
        let mut an1 = vec![mk(at1, 0)];
        let mut an2 = vec![mk(at2, 1)];
        while !an1.is_empty() || !an2.is_empty() {
            stable_sort_cip(self, &st, &mut an1);
            stable_sort_cip(self, &st, &mut an2);
            self.assign_medals(&mut st, &an1);
            self.assign_medals(&mut st, &an2);
            self.choose_first_and_sort_accordingly(&st, &mut an1);
            self.choose_first_and_sort_accordingly(&st, &mut an2);
            self.finalize_scores(&mut st, &an1);
            self.finalize_scores(&mut st, &an2);
            let nn = an1.len().min(an2.len());
            for i in 0..nn {
                if cip_less(self, &st, &an1[i], &an2[i], None) {
                    return Some(at1);
                }
                if cip_less(self, &st, &an2[i], &an1[i], None) {
                    return Some(at2);
                }
            }
            if an1.len() > an2.len() {
                return Some(at1);
            }
            if an2.len() > an1.len() {
                return Some(at2);
            }
            an1 = self.expand_one_level(&mut st, &an1);
            an2 = self.expand_one_level(&mut st, &an2);
        }
        None
    }

    fn choose_first_and_sort_accordingly(&self, st: &[CipState; 2], v: &mut Vec<CipAtom>) {
        if v.len() < 2 {
            return;
        }
        let mut copy = std::mem::take(v);
        let mut friends_mask: HashMap<AtomId, u32> = HashMap::new();
        while !copy.is_empty() {
            let mut best_i = 0;
            for i in 1..copy.len() {
                if cip_less(self, st, &copy[i], &copy[best_i], Some(&friends_mask)) {
                    best_i = i;
                }
            }
            let new_best = copy.remove(best_i);
            for &p in &new_best.all_parents {
                *friends_mask.entry(p).or_insert(0) |= 1u32.wrapping_shl(copy.len() as u32);
            }
            v.push(new_best);
        }
    }

    fn expand_one_level(&self, st: &mut [CipState; 2], old: &[CipAtom]) -> Vec<CipAtom> {
        let mut new_v = Vec::new();
        for an in old {
            for aa in 0..an.these_atoms.len() {
                let Some(a) = an.these_atoms[aa].1 else {
                    continue; // dummy atom
                };
                let state = &mut st[an.side];
                let mut all_parents = an.all_parents.clone();
                all_parents.push(a);
                let mut these_ats: Vec<(i32, Option<AtomId>)> = Vec::new();
                let at = &self.atoms[a];
                for n in 0..at.bonds.len() {
                    let neigh = at.neighbors[n];
                    let ne = &self.atoms[neigh];
                    if neigh != an.parent {
                        let vn = state.visited.get(&neigh).copied().unwrap_or(0);
                        let va = state.visited.get(&a).copied().unwrap_or(0);
                        // hydrogens and atoms closing a ring become ghosts
                        let ghost = ne.atomic_number == 1 || (vn != 0 && vn != va + 1);
                        these_ats.push((ne.atomic_number, if ghost { None } else { Some(neigh) }));
                        if !ghost {
                            state.visited.insert(neigh, va + 1);
                        }
                    }
                    // ghosts for multiple order bonds, even to the parent
                    match self.bonds[at.bonds[n]].bond_order {
                        2 => these_ats.push((ne.atomic_number, None)),
                        3 => {
                            these_ats.push((ne.atomic_number, None));
                            these_ats.push((ne.atomic_number, None));
                        }
                        _ => {}
                    }
                }
                for _ in 0..at.implicit_hs {
                    these_ats.push((1, None));
                }
                these_ats.sort();
                these_ats.reverse();
                new_v.push(CipAtom {
                    these_atoms: these_ats,
                    parent: a,
                    all_parents,
                    side: an.side,
                });
            }
        }
        new_v
    }

    fn assign_medals(&self, st: &mut [CipState; 2], v: &[CipAtom]) {
        if v.is_empty() {
            return;
        }
        let side = v[0].side;
        let mut equal_to_previous = vec![false; v.len()];
        for i in 1..v.len() {
            equal_to_previous[i] = cip_equal(self, st, &v[i], &v[i - 1]);
        }
        let mut medal_lvl = 0usize;
        for i in 0..v.len() {
            if i > 0 && equal_to_previous[i] {
                medal_lvl -= 1;
            }
            for &p in &v[i].all_parents {
                let mut medals_v = st[side].medals_of(p);
                while medals_v.len() < medal_lvl {
                    medals_v.push(0);
                }
                if medals_v.len() > medal_lvl {
                    medals_v[medal_lvl] += 1;
                } else {
                    medals_v.push(1);
                }
                st[side].medals.insert(p, medals_v);
            }
            medal_lvl += 1;
        }
    }

    fn finalize_scores(&self, st: &mut [CipState; 2], v: &[CipAtom]) {
        if v.is_empty() {
            return;
        }
        let side = v[0].side;
        let mut equal_to_previous = vec![false; v.len()];
        for i in 1..v.len() {
            equal_to_previous[i] = cip_equal(self, st, &v[i], &v[i - 1]);
        }
        st[side].scores.clear();
        let mut score = 1;
        for i in 0..v.len() {
            if i > 0 && equal_to_previous[i] {
                score -= 1;
            }
            for &p in &v[i].all_parents {
                if st[side].score(p) == 0 {
                    st[side].scores.insert(p, score);
                }
            }
            score += 1;
        }
        st[side].medals.clear();
    }

    /// `setCIPPriorities`: ranks four substituents (priority 0 is the highest); false if two of them are equivalent
    pub fn set_cip_priorities(&self, atom_priorities: &mut [AtomPriority], center: AtomId) -> bool {
        for p in atom_priorities.iter_mut() {
            p.priority = 3;
        }
        if atom_priorities.len() != 4 {
            return false;
        }
        for i in 0..3 {
            for j in i + 1..4 {
                let res = self.cip_priority(atom_priorities[i].a, atom_priorities[j].a, center);
                if res == Some(atom_priorities[i].a) {
                    atom_priorities[i].priority -= 1;
                } else if res == Some(atom_priorities[j].a) {
                    atom_priorities[j].priority -= 1;
                }
            }
        }
        let mut found = [false; 4];
        for p in atom_priorities.iter() {
            if found[p.priority as usize] {
                return false; // two atoms have the same priority
            }
            found[p.priority as usize] = true;
        }
        true
    }

    /// `orderAtomPriorities`: orders the substituents trying to keep long chains in position 2 and 3 and side
    /// substituents in 1 and 4
    pub fn order_atom_priorities(
        &mut self,
        atom_priorities: &mut Vec<AtomPriority>,
        center: AtomId,
    ) {
        debug_assert_eq!(atom_priorities.len(), 4);
        let mut weights = [0.0f32; 4];
        let mol = self.atoms[center].molecule.expect("center in a molecule");
        for i in 0..4 {
            for a in self.molecules[mol].atoms.clone() {
                self.atoms[a].general_use_visited = false;
            }
            let mut q = VecDeque::new();
            q.push_back(atom_priorities[i].a);
            self.atoms[center].general_use_visited = true;
            self.atoms[atom_priorities[i].a].general_use_visited = true;
            let mut counter = 0;
            while let Some(at) = q.pop_front() {
                counter += 1;
                for n in self.atoms[at].neighbors.clone() {
                    if !self.atoms[n].general_use_visited {
                        q.push_back(n);
                        self.atoms[n].general_use_visited = true;
                    }
                }
            }
            weights[i] = counter as f32;
            let ap = atom_priorities[i].a;
            if let Some(b) = self.bond_to(center, ap) {
                let bd = &self.bonds[b];
                if bd.bond_order == 2 {
                    weights[i] = (f64::from(weights[i]) - 0.25) as f32;
                }
                if self.atoms[center].atomic_number == 16 && bd.bond_order == 2 {
                    weights[i] += 2000.0;
                }
                if self.share_a_ring(bd.start_atom, bd.end_atom).is_some() {
                    weights[i] += 500.0;
                }
            }
            let a = &self.atoms[ap];
            if a.atomic_number == 6 {
                weights[i] = (f64::from(weights[i]) + 0.5) as f32;
            }
            if a.atomic_number == 1 {
                weights[i] = (f64::from(weights[i]) - 0.5) as f32;
            }
            if a.is_shared_and_inner && !self.atoms[center].is_shared_and_inner {
                weights[i] -= 2000.0;
            }
            if self.atoms[center].cross_layout && a.neighbors.len() > 1 {
                weights[i] += 200.0;
            }
            if a.has_stereochemistry_set {
                weights[i] += 10000.0;
            }
            if a.bonds.iter().any(|&b| self.bonds[b].bond_order == 2) {
                weights[i] += 100.0;
            }
        }
        let mut lowest_weight = weights[0];
        let mut index = 0;
        for i in 1..4 {
            if weights[i] < lowest_weight {
                lowest_weight = weights[i];
                index = i;
            }
        }
        let first_atom = atom_priorities.remove(index);
        let mut w: Vec<f32> = weights.to_vec();
        w.remove(index);
        index = 0;
        lowest_weight = w[0];
        for i in 1..3 {
            if w[i] < lowest_weight {
                lowest_weight = w[i];
                index = i;
            }
        }
        let second_atom = atom_priorities.remove(index);
        let c = &self.atoms[center];
        if (c.atomic_number != 16 && c.atomic_number != 15) || c.neighbors.len() != 4 {
            atom_priorities.push(second_atom);
            atom_priorities.insert(0, first_atom);
        } else {
            atom_priorities.insert(1, second_atom);
            atom_priorities.insert(0, first_atom);
        }
    }

    /// neighbors of the atom ordered by the angle to the first one; the angles between consecutive ones
    fn angle_ordered_neighbors(
        &self,
        a: AtomId,
        dummy_first_angle: bool,
    ) -> (Vec<AtomId>, Vec<BondId>, Vec<f32>) {
        let _ = dummy_first_angle;
        let at = &self.atoms[a];
        let mut neighs = at.neighbors.clone();
        let mut bbonds = at.bonds.clone();
        let mut ordered_neighs = Vec::new();
        let mut ordered_bonds = Vec::new();
        let mut angles = Vec::new();
        let mut last_popped_index = 0;
        let mut last_popped_atom = neighs[last_popped_index];
        ordered_neighs.push(last_popped_atom);
        neighs.remove(last_popped_index);
        ordered_bonds.push(bbonds.remove(last_popped_index));
        while !neighs.is_empty() {
            let mut smallest_angle = 361.0f32;
            for (i, &n) in neighs.iter().enumerate() {
                let mut new_angle = signed_angle(
                    self.atoms[last_popped_atom].coordinates,
                    at.coordinates,
                    self.atoms[n].coordinates,
                );
                if new_angle < 0.0 {
                    new_angle += 360.0;
                }
                if new_angle < smallest_angle {
                    smallest_angle = new_angle;
                    last_popped_index = i;
                }
            }
            angles.push(smallest_angle);
            last_popped_atom = neighs[last_popped_index];
            ordered_neighs.push(last_popped_atom);
            neighs.remove(last_popped_index);
            ordered_bonds.push(bbonds.remove(last_popped_index));
        }
        (ordered_neighs, ordered_bonds, angles)
    }

    /// `readStereochemistry`: 0 if not assigned, 1 if R, -1 if S (in coordgen coordinates)
    pub fn read_stereochemistry(&mut self, a: AtomId, read_only: bool) -> i32 {
        if !read_only {
            self.atoms[a].has_ring_chirality = false;
            self.atoms[a].is_stereogenic = false;
        }
        if !self.can_be_chiral(a) {
            return 0;
        }
        let mol = self.atoms[a].molecule;
        let base = self.atoms.len();
        let dummy_h = self.push_dummy(1, mol);
        let dummy_lp = self.push_dummy(0, mol);
        let result = self.read_stereochemistry_inner(a, read_only, dummy_h, dummy_lp);
        self.atoms.truncate(base);
        result
    }

    fn read_stereochemistry_inner(
        &mut self,
        a: AtomId,
        read_only: bool,
        dummy_h: AtomId,
        dummy_lp: AtomId,
    ) -> i32 {
        let (mut ordered_neighs, ordered_bonds, mut angles) =
            self.angle_ordered_neighbors(a, false);
        let mut total_angle = 0.0f32;
        for &angle in &angles {
            total_angle += angle;
        }
        angles.push(360.0 - total_angle);
        let mut semiplane = false;
        let mut central_atom = None;
        if angles.len() == 3 {
            for i in 0..angles.len() {
                if angles[i] > 180.0 {
                    semiplane = true;
                    let prec_i = if i == 0 { angles.len() - 1 } else { i - 1 };
                    central_atom = Some(ordered_neighs[prec_i]);
                }
            }
        }
        let mut wedge_n: i32 = -1;
        let mut dashed_n: i32 = -1;
        let mut give_up = false;
        for (i, &ob) in ordered_bonds.iter().enumerate() {
            let b = &self.bonds[ob];
            if b.has_stereochemistry_display {
                if !b.is_reversed && b.start_atom != a {
                    continue;
                }
                if b.is_reversed && b.end_atom != a {
                    continue;
                }
                let mut wedge = b.is_wedge;
                if b.is_reversed {
                    wedge = !wedge;
                }
                if b.start_atom != a {
                    wedge = !wedge;
                }
                if wedge {
                    if wedge_n == -1 {
                        wedge_n = i as i32;
                    } else {
                        give_up = true;
                    }
                } else if dashed_n == -1 {
                    dashed_n = i as i32;
                } else {
                    give_up = true;
                }
            }
        }
        let mut start_index: i32 = 0;
        let mut invert = false;
        if dashed_n == -1 && wedge_n == -1 {
            give_up = true;
        } else if dashed_n == -1 {
            start_index = wedge_n;
        } else if wedge_n == -1 {
            start_index = dashed_n;
            invert = true;
        } else {
            if ordered_bonds.len() == 3 {
                return 0;
            }
            if wedge_n - dashed_n == 1 || wedge_n - dashed_n == -3 {
                start_index = wedge_n;
            } else if wedge_n - dashed_n == -1 || wedge_n - dashed_n == 3 {
                start_index = dashed_n;
                invert = true;
            } else {
                give_up = true;
            }
        }
        let n_of_hs = self.atoms[a].implicit_hs;
        let total_substituents_n = self.atoms[a].neighbors.len() as i32 + n_of_hs;
        let z = self.atoms[a].atomic_number;
        let si = start_index as usize;
        if ordered_neighs.len() == 3 && n_of_hs == 1 {
            if semiplane && central_atom == Some(ordered_neighs[si]) {
                invert = !invert;
            }
            ordered_neighs.insert(si, dummy_h);
            invert = !invert;
        }
        if ordered_neighs.len() == 3 && (z == 7 || z == 16) {
            ordered_neighs.insert(si, dummy_lp);
            invert = !invert;
        }
        let mut atom_priorities: Vec<AtomPriority> = ordered_neighs
            .iter()
            .map(|&n| AtomPriority { a: n, priority: 0 })
            .collect();
        if atom_priorities.len() != 4 {
            return 0;
        }
        let mut canonical = vec![0, 1, 2, 3];
        self.atoms[a].rs_priorities.clear();
        let mut is_stereocenter = self.set_cip_priorities(&mut atom_priorities, a);
        if !is_stereocenter && !read_only && !self.atoms[a].ignore_ring_chirality {
            self.atoms[a].has_ring_chirality = true;
            is_stereocenter = self.set_cip_priorities(&mut atom_priorities, a);
        }
        if !is_stereocenter && !read_only {
            self.atoms[a].has_ring_chirality = false;
        }
        if !is_stereocenter {
            give_up = true;
        } else if !read_only {
            self.atoms[a].is_stereogenic = true;
        }
        if total_substituents_n < 4 && z != 7 && z != 16 && !read_only {
            self.atoms[a].is_stereogenic = false;
        }
        if !self.atoms[a].is_stereogenic && !read_only {
            self.atoms[a].has_ring_chirality = false;
        }
        let mut rs: Vec<i32> = Vec::new();
        for n in self.atoms[a].neighbors.clone() {
            if let Some(p) = atom_priorities.iter().find(|p| p.a == n) {
                rs.push(p.priority as i32);
            }
        }
        // may pad `rs` with a trailing 3 as the C++ code does
        let matched = Self::match_cip_sequence(&mut canonical, &mut rs);
        self.atoms[a].rs_priorities = rs;
        self.atoms[a].clockwise_invert = !matched;
        if !give_up {
            let mut out_of_place = 0;
            for i in 0..atom_priorities.len() {
                let mut n = si + i;
                if n > 3 {
                    n -= 4;
                }
                if atom_priorities[n].priority != i as u32 {
                    out_of_place += 1;
                }
            }
            if out_of_place == 2 {
                invert = !invert;
            } else if out_of_place == 4 {
                let mut n2 = atom_priorities[si].priority as usize + si;
                if n2 > 3 {
                    n2 -= 4;
                }
                if atom_priorities[n2].priority != 0 {
                    invert = !invert;
                }
            }
            return if invert { -1 } else { 1 };
        }
        0
    }

    /// `getRelativeStereo`
    pub fn get_relative_stereo(
        &mut self,
        a: AtomId,
        looking_from: AtomId,
        atom1: AtomId,
        atom2: AtomId,
    ) -> Chirality {
        self.read_stereochemistry(a, false);
        let rs = self.atoms[a].rs_priorities.clone();
        if rs.len() < 3 {
            return Chirality::Unspecified;
        }
        let mut priorities = vec![3, 3, 3, 3];
        for (nn, &n) in self.atoms[a].neighbors.iter().enumerate() {
            if n == atom1 {
                priorities[0] = rs[nn];
            } else if n == atom2 {
                priorities[1] = rs[nn];
            } else if n == looking_from {
                priorities[3] = rs[nn];
            } else {
                priorities[2] = rs[nn];
            }
        }
        let mut can = vec![0, 1, 2, 3];
        let matched = Self::match_cip_sequence(&mut priorities, &mut can);
        let is_r = self.atoms[a].is_r;
        let is_clockwise = if matched { !is_r } else { is_r };
        if is_clockwise {
            Chirality::Clockwise
        } else {
            Chirality::CounterClockwise
        }
    }

    /// `setAbsoluteStereoFromChiralityInfo`: false on an inconsistent description
    pub fn set_absolute_stereo_from_chirality_info(&mut self, a: AtomId) -> bool {
        let info = self.atoms[a].chirality_info;
        if info.direction == Chirality::Unspecified {
            return true;
        }
        self.read_stereochemistry(a, false);
        let rs = self.atoms[a].rs_priorities.clone();
        if rs.len() < 3 {
            return false;
        }
        let mut priorities = vec![5, 5, 5, 5];
        let mut at3 = false;
        for (nn, &n) in self.atoms[a].neighbors.iter().enumerate() {
            if Some(n) == info.atom1 {
                priorities[0] = rs[nn];
            } else if Some(n) == info.atom2 {
                priorities[1] = rs[nn];
            } else if Some(n) == info.looking_from {
                priorities[3] = rs[nn];
            } else if at3 {
                return false; // more than 1 atom not matching
            } else {
                at3 = true;
                priorities[2] = rs[nn];
            }
        }
        let mut adding_hn = 0;
        for p in priorities.iter_mut() {
            if *p == 5 {
                *p = 3;
                adding_hn += 1;
            }
        }
        if adding_hn > 1 {
            return false;
        }
        let mut invert = false;
        let mut can = vec![0, 1, 2, 3];
        if !Self::match_cip_sequence(&mut priorities, &mut can) {
            invert = !invert;
        }
        let mut is_r = info.direction != Chirality::Clockwise;
        if invert {
            is_r = !is_r;
        }
        self.atoms[a].is_r = is_r;
        self.atoms[a].has_stereochemistry_set = true;
        true
    }

    /// `writeStereoChemistry`: sets the wedges of the bonds from `has_stereochemistry_set` and `is_r`
    pub fn write_stereo_chemistry(&mut self, a: AtomId) {
        if !self.atoms[a].has_stereochemistry_set {
            return;
        }
        if !self.can_be_chiral(a) {
            self.atoms[a].has_stereochemistry_set = false;
            return;
        }
        let mol = self.atoms[a].molecule;
        let base = self.atoms.len();
        let dummy_h = self.push_dummy(1, mol);
        let dummy_lp = self.push_dummy(0, mol);
        self.write_stereo_chemistry_inner(a, dummy_h, dummy_lp);
        self.atoms.truncate(base);
    }

    fn write_stereo_chemistry_inner(&mut self, a: AtomId, dummy_h: AtomId, dummy_lp: AtomId) {
        let mut dummy = dummy_h;
        let (ordered_neighs, ordered_bonds, mut angles) = self.angle_ordered_neighbors(a, true);
        let z = self.atoms[a].atomic_number;
        if (z == 7 || z == 16) && self.atoms[a].implicit_hs == 0 && ordered_bonds.len() == 3 {
            dummy = dummy_lp;
        }
        let mut four = true;
        let total_angle: f32 = angles.iter().fold(0.0f32, |s, &x| s + x);
        angles.push(360.0 - total_angle);
        let mut atom_priorities: Vec<AtomPriority> = ordered_neighs
            .iter()
            .map(|&n| AtomPriority { a: n, priority: 0 })
            .collect();
        if atom_priorities.len() == 3 {
            four = false;
            atom_priorities.push(AtomPriority {
                a: dummy,
                priority: 0,
            });
        }
        let mut is_stereocenter = self.set_cip_priorities(&mut atom_priorities, a);
        if !is_stereocenter && !self.atoms[a].ignore_ring_chirality {
            self.atoms[a].has_ring_chirality = true;
            is_stereocenter = self.set_cip_priorities(&mut atom_priorities, a);
        }
        if !is_stereocenter {
            self.atoms[a].has_ring_chirality = false;
        }
        let mut ordered_atom_priorities = atom_priorities.clone();
        self.order_atom_priorities(&mut ordered_atom_priorities, a);
        if !is_stereocenter {
            for b in self.atoms[a].bonds.clone() {
                self.bonds[b].has_stereochemistry_display = false;
            }
            return;
        }
        if !four && ordered_atom_priorities[0].a == dummy {
            let first = ordered_atom_priorities[0];
            ordered_atom_priorities.push(first);
            let fourth = ordered_atom_priorities[3];
            ordered_atom_priorities.insert(0, fourth);
            ordered_atom_priorities.remove(1);
            ordered_atom_priorities.remove(3);
        }
        let mut start_index = 0;
        for (i, p) in atom_priorities.iter().enumerate() {
            if p.a == ordered_atom_priorities[0].a {
                start_index = i;
                break;
            }
        }
        for _ in 0..start_index {
            let first = atom_priorities.remove(0);
            atom_priorities.push(first);
        }
        if four && atom_priorities[1].a == ordered_atom_priorities[3].a {
            let first = atom_priorities.remove(0);
            atom_priorities.push(first);
        }
        let mut main_atom = None;
        if four {
            if atom_priorities[0].a == ordered_atom_priorities[0].a {
                main_atom = Some(atom_priorities[0].a);
            }
            if atom_priorities[3].a == ordered_atom_priorities[0].a {
                main_atom = Some(atom_priorities[3].a);
            }
        }
        let mut invert = false;
        let b1 = self.bond_to(a, atom_priorities[0].a);
        self.set_wedge_display(a, b1, atom_priorities[0].a, four, main_atom, !invert, true);
        if four {
            let b2 = self.bond_to(a, atom_priorities[3].a);
            self.set_wedge_display(a, b2, atom_priorities[3].a, four, main_atom, invert, false);
        }
        if let Some(b3) = self.bond_to(a, atom_priorities[1].a) {
            self.bonds[b3].has_stereochemistry_display = false;
        }
        if let Some(b4) = self.bond_to(a, atom_priorities[2].a) {
            self.bonds[b4].has_stereochemistry_display = false;
        }
        let mut read_s = self.read_stereochemistry(a, true);
        read_s = -read_s; // isR is in sketcher coordinates, readStereochemistry in coordgen coordinates (flipped y)
        let is_r = self.atoms[a].is_r;
        if (read_s == -1 && is_r) || (read_s == 1 && !is_r) {
            invert = true;
            self.set_wedge_display(a, b1, atom_priorities[0].a, four, main_atom, !invert, true);
            if four {
                let b2 = self.bond_to(a, atom_priorities[3].a);
                self.set_wedge_display(a, b2, atom_priorities[3].a, four, main_atom, invert, false);
            }
        }
    }

    /// the wedge assignment of `writeStereoChemistry` for the bond to one substituent
    #[allow(clippy::too_many_arguments)]
    fn set_wedge_display(
        &mut self,
        a: AtomId,
        bond: Option<BondId>,
        other: AtomId,
        four: bool,
        main_atom: Option<AtomId>,
        wedge: bool,
        first: bool,
    ) {
        let Some(b) = bond else {
            return;
        };
        let cond = if first {
            !four
                || main_atom == Some(other)
                || (self.share_a_ring(a, other).is_none()
                    && !self.atoms[other].has_stereochemistry_set)
        } else {
            main_atom == Some(other)
                || (self.share_a_ring(a, other).is_none()
                    && !self.atoms[other].has_stereochemistry_set)
        };
        if cond {
            let reverse = self.bonds[b].start_atom != a;
            let bd = &mut self.bonds[b];
            bd.is_wedge = wedge;
            bd.has_stereochemistry_display = true;
            bd.is_reversed = reverse;
        } else {
            self.bonds[b].has_stereochemistry_display = false;
        }
    }

    // ---- bond ----

    /// `startAtomCIPFirstNeighbor`
    pub fn bond_start_atom_cip_first_neighbor(&self, b: BondId) -> Option<AtomId> {
        let bd = &self.bonds[b];
        if bd.bond_order != 2 {
            return None;
        }
        self.cip_first_neighbor(bd.start_atom, bd.end_atom)
    }

    /// `endAtomCIPFirstNeighbor`
    pub fn bond_end_atom_cip_first_neighbor(&self, b: BondId) -> Option<AtomId> {
        let bd = &self.bonds[b];
        if bd.bond_order != 2 {
            return None;
        }
        self.cip_first_neighbor(bd.end_atom, bd.start_atom)
    }

    fn cip_first_neighbor(&self, atom: AtomId, other: AtomId) -> Option<AtomId> {
        let nb = &self.atoms[atom].neighbors;
        if nb.len() == 2 {
            Some(if nb[0] == other { nb[1] } else { nb[0] })
        } else if nb.len() == 3 {
            let ats: Vec<AtomId> = nb.iter().copied().filter(|&n| n != other).collect();
            if ats.len() == 2 {
                self.cip_priority(ats[0], ats[1], atom)
            } else {
                None
            }
        } else {
            None
        }
    }

    /// `setAbsoluteStereoFromStereoInfo`
    pub fn bond_set_absolute_stereo_from_stereo_info(&mut self, b: BondId) {
        let bd = &self.bonds[b];
        if self.bond_is_stereo(b) && bd.stereo.atom1.is_some() && bd.stereo.atom2.is_some() {
            let first_start = self.bond_start_atom_cip_first_neighbor(b);
            let first_end = self.bond_end_atom_cip_first_neighbor(b);
            if first_start.is_some() && first_end.is_some() {
                let st = self.bonds[b].stereo;
                let mut invert = false;
                if st.atom1 != first_start && st.atom1 != first_end {
                    invert = !invert;
                }
                if st.atom2 != first_start && st.atom2 != first_end {
                    invert = !invert;
                }
                let mut setting_is_z = st.stereo == BondStereo::Cis;
                if invert {
                    setting_is_z = !setting_is_z;
                }
                self.bonds[b].is_z = setting_is_z;
            }
        }
        if self.bonds[b].stereo.stereo == BondStereo::Unspecified {
            self.bonds[b].ignore_ze = true;
        }
    }

    /// `checkStereoChemistry`
    pub fn bond_check_stereo_chemistry(&self, b: BondId) -> bool {
        if !self.bond_is_stereo(b) || self.bond_is_in_small_ring(b) {
            return true;
        }
        let Some(first_start) = self.bond_start_atom_cip_first_neighbor(b) else {
            return true;
        };
        let Some(first_end) = self.bond_end_atom_cip_first_neighbor(b) else {
            return true;
        };
        let bd = &self.bonds[b];
        super::maths::same_side(
            self.atoms[first_start].coordinates,
            self.atoms[first_end].coordinates,
            self.atoms[bd.start_atom].coordinates,
            self.atoms[bd.end_atom].coordinates,
        ) == bd.is_z
    }

    /// `markedAsCis`
    pub fn bond_marked_as_cis(&self, b: BondId, atom1: AtomId, atom2: AtomId) -> bool {
        let f1 = self.bond_start_atom_cip_first_neighbor(b);
        let f2 = self.bond_end_atom_cip_first_neighbor(b);
        let mut cis = self.bonds[b].is_z;
        if Some(atom1) != f1 && Some(atom1) != f2 {
            cis = !cis;
        }
        if Some(atom2) != f1 && Some(atom2) != f2 {
            cis = !cis;
        }
        cis
    }
}

#[cfg(test)]
mod tests {
    use super::super::testutil::build;
    use super::*;

    fn run(input: &str) -> String {
        let (mut model, mol) = build(input);
        // the wedge column of the bond lines
        let bond_lines = input.lines().filter(|l| l.starts_with("B "));
        for (b, l) in bond_lines.enumerate() {
            let w: i32 = l.split(' ').nth(4).map_or(0, |x| x.parse().unwrap());
            if w != 0 {
                model.bonds[b].has_stereochemistry_display = true;
                model.bonds[b].is_wedge = w == 1;
            }
        }
        model.force_update_struct(mol);
        let mut out = String::new();
        let n = model.molecules[mol].atoms.len();
        let atoms = model.molecules[mol].atoms.clone();
        let idx = |a: Option<AtomId>| {
            a.map_or(-1, |a| atoms.iter().position(|&x| x == a).unwrap() as i64)
        };
        for (i, &a) in atoms.iter().enumerate() {
            let nb = model.atoms[a].neighbors.clone();
            for j in 0..nb.len().saturating_sub(1) {
                for k in j + 1..nb.len() {
                    let r = model.cip_priority(nb[j], nb[k], a);
                    out += &format!("p {i} {j} {k} {}\n", idx(r));
                }
            }
        }
        for (i, &a) in atoms.iter().enumerate() {
            let r = model.read_stereochemistry(a, false);
            let at = &model.atoms[a];
            out += &format!(
                "r {i} {r} stereogenic {} ring {} inv {} rs",
                at.is_stereogenic as i32, at.has_ring_chirality as i32, at.clockwise_invert as i32
            );
            for p in &at.rs_priorities {
                out += &format!(" {p}");
            }
            out += "\n";
        }
        for (i, &a) in atoms.iter().enumerate() {
            model.atoms[a].has_stereochemistry_set = true;
            model.atoms[a].is_r = i % 2 == 0;
            model.write_stereo_chemistry(a);
            let at = &model.atoms[a];
            out += &format!(
                "w {i} set {} ring {}\n",
                at.has_stereochemistry_set as i32, at.has_ring_chirality as i32
            );
        }
        let bonds = model.molecules[mol].bonds.clone();
        for (i, &b) in bonds.iter().enumerate() {
            let f = model.bond_start_atom_cip_first_neighbor(b);
            let e = model.bond_end_atom_cip_first_neighbor(b);
            let bd = &model.bonds[b];
            out += &format!(
                "b {i} disp {} wedge {} rev {} first {} {}\n",
                bd.has_stereochemistry_display as i32,
                bd.is_wedge as i32,
                bd.is_reversed as i32,
                idx(f),
                idx(e)
            );
        }
        let _ = n;
        out
    }

    /// CIP ranking and stereo reading/writing must be the ones of coordgenlibs (`tests/data/coordgen_cip.txt`)
    #[test]
    fn cip_matches_cpp() {
        let text = include_str!("../../../../../tests/data/coordgen_cip.txt");
        let mut bad = Vec::new();
        let mut n = 0;
        for rec in text.split("$$\n").filter(|r| !r.is_empty()) {
            let end = rec.find("END\n").unwrap() + 4;
            let (input, want) = rec.split_at(end);
            let got = run(input);
            n += 1;
            if got != want {
                let pos = got
                    .lines()
                    .zip(want.lines())
                    .position(|(a, b)| a != b)
                    .unwrap_or(0);
                bad.push(format!(
                    "{}\nline {pos}\nwant {:?}\ngot  {:?}",
                    input.lines().next().unwrap(),
                    want.lines().nth(pos),
                    got.lines().nth(pos)
                ));
            }
        }
        assert!(n > 50, "{n}");
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }
}
