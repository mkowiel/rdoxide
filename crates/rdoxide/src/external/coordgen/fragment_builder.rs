//! `CoordgenFragmentBuilder`: initial coordinates of a fragment: ring systems (templates, central rings, side rings),
//! then the non-ring atoms by breadth first search. Macrocycles (`CoordgenMacrocycleBuilder`) and the clash search
//! of newly placed atoms (`CoordgenMinimizer::checkForClashes`) are not ported yet and stop with `unimplemented!`.

use std::collections::{BTreeSet, VecDeque};

use super::cip::AtomPriority;
use super::fragment::DofKind;
use super::maths::{BONDLENGTH, PointF, SKETCHER_EPSILON};
use super::minimizer::Minimizer;
use super::model::{AtomId, BondId, FragmentId, Model, RingId};
use super::templates::{alignment_matrix, templates};

const MACROCYCLE: usize = super::maths::MACROCYCLE;
const PERFECTLY_PLANAR_SYSTEM_SCORE: f32 = 50.0;
const NON_PLANAR_SYSTEM_SCORE: f32 = 1000.0;
const UNTREATABLE_SYSTEM_PLANARITY_SCORE: f32 = 200000.0;

const MACROCYCLE_CENTRAL_RING_SCORE: usize = 1000;
const NUMBER_OF_FUSED_RINGS_CENTRAL_RING_SCORE: usize = 40;
const NUMBER_OF_FUSION_ATOMS_CENTRAL_RING_SCORE: usize = 15;
const NEIGHBOR_ALREADY_BUILT_RING_SCORE: usize = 100000;

const PI: f64 = std::f64::consts::PI;

impl Model {
    /// `initializeCoordinates`
    pub fn initialize_coordinates(&mut self, fragment: FragmentId) {
        assert!(!self.fragments[fragment].atoms().is_empty());
        self.build_fragment(fragment);
        self.fragment_store_coordinate_information(fragment);
    }

    /// `sketcherMinimizerFragment::storeCoordinateInformation`
    pub fn fragment_store_coordinate_information(&mut self, f: FragmentId) {
        let mut coords = std::collections::BTreeMap::new();
        let mut origin = PointF::new(0.0, 0.0);
        let mut angle = 0.0f32;
        let frag = &self.fragments[f];
        if let Some(b) = frag.bond_to_parent {
            let bd = &self.bonds[b];
            origin = self.atoms[bd.end_atom].coordinates;
            let s = self.atoms[bd.start_atom].coordinates;
            angle = (s.y - origin.y).atan2(-s.x + origin.x);
        } else if !frag.constrained && !frag.fixed {
            origin = self.atoms[frag.atoms()[0]].coordinates;
        }
        let cosine = (-angle).cos();
        let sine = (-angle).sin();
        for &at in frag.atoms() {
            let mut c = self.atoms[at].coordinates - origin;
            c.rotate(sine, cosine);
            coords.insert(at, c);
        }
        for &child in &frag.children {
            let at = self.bonds[self.fragments[child].bond_to_parent.expect("child bond")].end_atom;
            let mut c = self.atoms[at].coordinates - origin;
            c.rotate(sine, cosine);
            coords.insert(at, c);
        }
        self.fragments[f].coordinates = coords;
    }

    /// `rotateMainFragment`
    fn rotate_main_fragment(&mut self, f: FragmentId) {
        let frag = &self.fragments[f];
        if frag.fixed || frag.is_templated || !frag.constrained {
            return;
        }
        let mut old_center = PointF::new(0.0, 0.0);
        let mut new_center = PointF::new(0.0, 0.0);
        let mut constrained_atoms: Vec<AtomId> = Vec::new();
        for &a in frag.atoms() {
            if self.atoms[a].constrained {
                constrained_atoms.push(a);
            }
        }
        let children = frag.children.clone();
        for &child in &children {
            let atom = self.bonds[self.fragments[child].bond_to_parent.unwrap()].end_atom;
            if self.atoms[atom].constrained {
                constrained_atoms.push(atom);
            }
        }
        for &a in &constrained_atoms {
            old_center += self.atoms[a].template_coordinates;
            new_center += self.atoms[a].coordinates;
        }
        if !constrained_atoms.is_empty() {
            old_center = old_center / constrained_atoms.len() as f32;
            new_center = new_center / constrained_atoms.len() as f32;
        }
        let mut v1 = Vec::new();
        let mut v2 = Vec::new();
        for &a in &constrained_atoms {
            v2.push(self.atoms[a].coordinates - new_center);
            v1.push(self.atoms[a].template_coordinates - old_center);
        }
        let rot = alignment_matrix(&v1, &v2);
        let rotate =
            |p: PointF| PointF::new(p.x * rot[0] + p.y * rot[1], p.x * rot[2] + p.y * rot[3]);
        for a in self.fragments[f].atoms().to_vec() {
            let v = rotate(self.atoms[a].coordinates - new_center);
            self.set_coordinates(a, v + old_center);
        }
        for child in children {
            let at = self.bonds[self.fragments[child].bond_to_parent.unwrap()].end_atom;
            let v = rotate(self.atoms[at].coordinates - new_center);
            self.set_coordinates(at, v + old_center);
            self.atoms[at].coordinates_set = false;
        }
    }

    /// `findTemplate`
    fn find_template(&mut self, rings: &[RingId]) -> bool {
        let mut atoms: Vec<AtomId> = Vec::new();
        let mut bonds: Vec<BondId> = Vec::new();
        let mut visited: BTreeSet<AtomId> = BTreeSet::new();
        let mut old_indices = Vec::new();
        for &r in rings {
            for a in self.rings[r].atoms.clone() {
                if visited.insert(a) {
                    old_indices.push(self.atoms[a].general_use_n);
                    self.atoms[a].general_use_n = atoms.len() as i32;
                    atoms.push(a);
                }
            }
        }
        for &r in rings {
            for &b in &self.rings[r].bonds {
                if !bonds.contains(&b) {
                    bonds.push(b);
                }
            }
        }
        let mut found_template = false;
        for temp in templates() {
            if let Some(mapping) = self.compare_with_template(&atoms, &bonds, temp) {
                found_template = true;
                if let Some(&a0) = atoms.first() {
                    let f = self.atoms[a0].fragment.expect("fragment");
                    self.fragments[f].is_templated = true;
                }
                for (i, &a) in atoms.iter().enumerate() {
                    self.set_coordinates(a, temp.coordinates[mapping[i]] * BONDLENGTH);
                    self.atoms[a].rigid = true;
                }
                for &r in rings {
                    self.rings[r].coordinates_generated = true;
                }
                break;
            }
        }
        for (i, &old) in old_indices.iter().enumerate() {
            self.atoms[atoms[i]].general_use_n = old;
        }
        found_template
    }

    /// `findCentralRingOfSystem`
    fn find_central_ring_of_system(&self, rings: &[RingId]) -> RingId {
        let mut highest: Option<RingId> = None;
        let mut high_score = 0usize;
        for &r in rings {
            let ring = &self.rings[r];
            let mut priority = 0usize;
            if ring
                .fused_with
                .iter()
                .any(|&n| self.rings[n].coordinates_generated)
            {
                priority += NEIGHBOR_ALREADY_BUILT_RING_SCORE;
            }
            if ring.is_macrocycle() {
                priority += MACROCYCLE_CENTRAL_RING_SCORE;
            }
            if ring.atoms.len() == 6 {
                priority += 10;
            }
            priority += ring.atoms.len();
            priority += NUMBER_OF_FUSED_RINGS_CENTRAL_RING_SCORE * ring.fused_with.len();
            for fa in &ring.fusion_atoms {
                priority += NUMBER_OF_FUSION_ATOMS_CENTRAL_RING_SCORE * fa.len();
            }
            if priority > high_score || highest.is_none() {
                highest = Some(r);
                high_score = priority;
            }
        }
        highest.unwrap_or(rings[0])
    }

    /// `generateCoordinatesCentralRings`
    fn generate_coordinates_central_rings(&mut self, mut rings: Vec<RingId>) {
        if rings.len() == 1 {
            self.build_ring(rings[0]);
            return;
        }
        if self.find_template(&rings) {
            return;
        }
        let planarity_score = self.new_score_planarity(&rings);
        if planarity_score < NON_PLANAR_SYSTEM_SCORE {
            if planarity_score > PERFECTLY_PLANAR_SYSTEM_SCORE {
                self.find_template(&rings);
            }
            while !rings.is_empty() {
                let central = self.find_central_ring_of_system(&rings);
                self.build_ring(central);
                rings.retain(|&r| r != central);
            }
            // `maybeMinimizeRings` of the (now empty) list does nothing
        } else if planarity_score > UNTREATABLE_SYSTEM_PLANARITY_SCORE {
        } else {
            let first_ring = self.find_central_ring_of_system(&rings);
            self.open_cycle_and_generate_coords(first_ring);
            let a0 = self.rings[first_ring].atoms[0];
            let mol = self.atoms[a0].molecule.expect("molecule");
            self.molecules[mol].require_minimization = true;
        }
    }

    /// `newScorePlanarity`: above 1000 the system is not planar
    pub fn new_score_planarity(&self, rings: &[RingId]) -> f32 {
        let mut score = 0.0f32;
        for &r in rings {
            let ring = &self.rings[r];
            if ring.is_macrocycle() && self.find_bond_to_open(r).is_none() {
                continue;
            }
            if ring.is_macrocycle() {
                for &other in &ring.fused_with {
                    if self.rings[other].is_macrocycle() {
                        score += NON_PLANAR_SYSTEM_SCORE;
                    }
                }
            }
            for &b in &ring.bonds {
                let n = self.bonds[b].rings.len();
                if n > 2 {
                    score += (NON_PLANAR_SYSTEM_SCORE as usize * (n - 2)) as f32;
                }
            }
            for &a in &ring.atoms {
                let at = &self.atoms[a];
                if at.neighbors.len() > 3 {
                    let mut angle = 0.0f32;
                    for &rr in &at.rings {
                        angle += (PI - (2.0 * PI / self.rings[rr].atoms.len() as f64)) as f32;
                    }
                    if f64::from(angle) >= 1.99 * PI {
                        score += NON_PLANAR_SYSTEM_SCORE;
                    }
                }
            }
        }
        score
    }

    /// `getSharedAtomsWithAlreadyDrawnRing`
    fn shared_atoms_with_drawn_ring(
        &self,
        ring: RingId,
    ) -> (Option<RingId>, Vec<AtomId>, Option<BondId>) {
        let mut parent: Option<RingId> = None;
        for &i in &self.rings[ring].fused_with {
            if self.rings[i].coordinates_generated {
                if let Some(p) = parent {
                    if self.ring_fusion_atoms_with(i, ring).len()
                        < self.ring_fusion_atoms_with(p, ring).len()
                        || self.rings[i].size() < self.rings[p].size()
                    {
                        continue;
                    }
                }
                parent = Some(i);
            }
        }
        let mut fusion_atoms = Vec::new();
        let mut fusion_bond = None;
        if let Some(p) = parent {
            for (i, &f) in self.rings[p].fused_with.iter().enumerate() {
                if f == ring {
                    fusion_atoms = self.rings[p].fusion_atoms[i].clone();
                }
            }
            for &b in &self.rings[p].fusion_bonds {
                let bd = &self.bonds[b];
                if self.rings[ring].contains_atom(bd.start_atom)
                    || self.rings[ring].contains_atom(bd.end_atom)
                {
                    fusion_bond = Some(b);
                    break;
                }
            }
        }
        (parent, fusion_atoms, fusion_bond)
    }

    /// `orderChainOfAtoms`
    pub(super) fn order_chain_of_atoms(&self, atoms: &[AtomId], start: AtomId) -> Vec<AtomId> {
        let mut ordered = Vec::new();
        let mut still: BTreeSet<AtomId> = atoms.iter().copied().collect();
        let mut atom_to_add = start;
        while ordered.len() < atoms.len() {
            ordered.push(atom_to_add);
            still.remove(&atom_to_add);
            if ordered.len() >= atoms.len() {
                break;
            }
            for &n in &self.atoms[atom_to_add].neighbors {
                if still.contains(&n) {
                    atom_to_add = n;
                    break;
                }
            }
        }
        ordered
    }

    /// `listOfCoordinatesFromListofRingAtoms`: regular polygon
    pub(super) fn polygon_coordinates(n_atoms: usize) -> Vec<PointF> {
        let a = (2.0 * PI / n_atoms as f64) as f32;
        let mut out = Vec::new();
        let mut coords = PointF::new(0.0, 0.0);
        for n in 0..n_atoms {
            out.push(coords);
            let angle = a * n as f32;
            coords += PointF::new(angle.cos() * BONDLENGTH, -angle.sin() * BONDLENGTH);
        }
        out
    }

    /// `buildRing`
    fn build_ring(&mut self, ring: RingId) {
        if self.rings[ring].coordinates_generated {
            return;
        }
        let (parent, fusion_atoms, fusion_bond) = self.shared_atoms_with_drawn_ring(ring);
        let ring_atoms = self.rings[ring].atoms.clone();
        let atoms = self.order_chain_of_atoms(&ring_atoms, ring_atoms[0]);
        let mut coords = if atoms.len() >= MACROCYCLE {
            self.new_macrocycle(ring, &atoms)
        } else {
            Self::polygon_coordinates(atoms.len())
        };
        let Some(parent) = parent else {
            for (i, &a) in atoms.iter().enumerate() {
                self.set_coordinates(a, coords[i]);
            }
            self.rings[ring].coordinates_generated = true;
            return;
        };
        if fusion_atoms.len() < 2 {
            let (pivot, pivot_on_parent) = if let Some(&f0) = fusion_atoms.first() {
                (f0, f0)
            } else {
                let fb = &self.bonds[fusion_bond.expect("fusion bond")];
                if atoms.contains(&fb.start_atom) {
                    (fb.start_atom, fb.end_atom)
                } else {
                    (fb.end_atom, fb.start_atom)
                }
            };
            debug_assert!(self.rings[parent].contains_atom(pivot_on_parent));
            let pivot_coordinates = atoms
                .iter()
                .position(|&a| a == pivot)
                .map(|i| coords[i])
                .unwrap_or_default();
            for c in coords.iter_mut() {
                *c -= pivot_coordinates;
            }
            let mut sum = PointF::new(0.0, 0.0);
            for &c in &coords {
                sum += c;
            }
            let mut center = sum / coords.len() as f32;
            center.normalize();
            let mut neighbors_mean = PointF::new(0.0, 0.0);
            for &n in &self.atoms[pivot_on_parent].neighbors {
                if self.rings[parent].contains_atom(n) {
                    neighbors_mean += self.atoms[n].coordinates;
                }
            }
            neighbors_mean = neighbors_mean * 0.5;
            let pivot_atom_coordinates = self.atoms[pivot_on_parent].coordinates;
            let mut direction = pivot_atom_coordinates - neighbors_mean;
            direction.normalize();
            let mut bond_vector = PointF::new(0.0, 0.0);
            if fusion_bond.is_some() {
                bond_vector = direction * BONDLENGTH;
            }
            for c in coords.iter_mut() {
                let mut coordinates = *c;
                coordinates.rotate(center.y, center.x);
                coordinates.rotate(-direction.y, direction.x);
                *c = coordinates + pivot_atom_coordinates + bond_vector;
            }
            for (i, &a) in atoms.iter().enumerate() {
                self.set_coordinates(a, coords[i]);
            }
            if let Some(fb) = fusion_bond {
                if !self.bond_check_stereo_chemistry(fb) {
                    for &a in &atoms {
                        self.mirror_coordinates(a, fb);
                    }
                }
            }
        } else {
            // two lists of coordinates, mirror of each other
            let mut coords2: Vec<PointF> = coords.iter().map(|c| PointF::new(c.x, -c.y)).collect();
            let idx = |a: AtomId| atoms.iter().position(|&x| x == a).unwrap();
            let last_i = fusion_atoms.len() - 1;
            let (f0, fl) = (fusion_atoms[0], fusion_atoms[last_i]);
            let mean = (self.atoms[f0].coordinates + self.atoms[fl].coordinates) * 0.5;
            let mut new_mean = (coords[idx(f0)] + coords[idx(fl)]) * 0.5;
            let v = self.atoms[fl].coordinates - self.atoms[f0].coordinates;
            let mut distance = v.length();
            if distance < SKETCHER_EPSILON {
                distance = SKETCHER_EPSILON;
            }
            let sine = v.y / distance;
            let cosine = v.x / distance;
            let mut new_v = coords[idx(fl)] - coords[idx(f0)];
            let mut new_distance = new_v.length();
            let mut new_sine = new_v.y / new_distance;
            let mut new_cosine = new_v.x / new_distance;
            for c in coords.iter_mut() {
                let mut d = *c - new_mean;
                d.rotate(new_sine, new_cosine);
                d.rotate(-sine, cosine);
                *c = mean + d;
            }
            // the second list is derived from the mirrored (not the rotated) coordinates
            new_v = coords2[idx(fl)] - coords2[idx(f0)];
            new_mean = (coords2[idx(f0)] + coords2[idx(fl)]) * 0.5;
            new_distance = new_v.length();
            new_sine = new_v.y / new_distance;
            new_cosine = new_v.x / new_distance;
            for c in coords2.iter_mut() {
                let mut d = *c - new_mean;
                d.rotate(new_sine, new_cosine);
                d.rotate(-sine, cosine);
                *c = mean + d;
            }
            let center = self.ring_find_center(parent);
            let (mut distance1, mut distance2) = (0.0f32, 0.0f32);
            for (i, &a) in atoms.iter().enumerate() {
                if fusion_atoms.contains(&a) {
                    continue;
                }
                distance1 += (coords[i] - center).length();
                distance2 += (coords2[i] - center).length();
            }
            let target = if distance1 < distance2 {
                &coords2
            } else {
                &coords
            };
            for (i, &a) in atoms.iter().enumerate() {
                self.set_coordinates(a, target[i]);
            }
        }
        self.rings[ring].coordinates_generated = true;
    }

    /// `buildNonRingAtoms`
    fn build_non_ring_atoms(&mut self, fragment: FragmentId) {
        let mut visited: BTreeSet<AtomId> = BTreeSet::new();
        let mut q: VecDeque<AtomId> = VecDeque::new();
        for &a in self.fragments[fragment].atoms() {
            if !self.atoms[a].rings.is_empty() {
                q.push_back(a);
                visited.insert(a);
            }
        }
        if q.is_empty() {
            let frag = &self.fragments[fragment];
            let btp = frag.bond_to_parent;
            if frag.parent.is_some() {
                let s = self.bonds[btp.unwrap()].start_atom;
                self.set_coordinates(s, PointF::new(-BONDLENGTH, 0.0));
            }
            let start = if self.fragments[fragment].parent.is_some() {
                self.bonds[btp.unwrap()].end_atom
            } else {
                self.fragments[fragment].atoms()[0]
            };
            self.set_coordinates(start, PointF::new(0.0, 0.0));
            q.push_back(start);
            visited.insert(start);
        }
        while !q.is_empty() {
            self.generate_coordinates_neighbors_of_first_atom_in_queue(
                &mut q,
                &mut visited,
                fragment,
            );
        }
    }

    /// `neighborsAnglesAtCenter` (`m_evenAngles` is false)
    fn neighbors_angles_at_center(&self, atom: AtomId) -> Vec<f32> {
        let at = &self.atoms[atom];
        let mut angle_division = at.neighbors.len();
        let mut angles = Vec::new();
        if at.neighbors.len() == 2 {
            if at.atomic_number == 6
                || !self.atoms[at.neighbors[0]].cross_layout
                || !self.atoms[at.neighbors[1]].cross_layout
            {
                angle_division = 3;
            }
            if self.bonds[at.bonds[0]].bond_order + self.bonds[at.bonds[1]].bond_order >= 4 {
                angle_division = 2;
            }
        } else if at.neighbors.len() == 4 && !at.cross_layout {
            angles.push((PI / 3.0) as f32);
            angles.push((PI * 0.5) as f32);
            angles.push((2.0 * PI / 3.0) as f32);
            angles.push((PI * 0.5) as f32);
        }
        if angles.is_empty() {
            for _ in 0..at.neighbors.len() {
                angles.push((2.0 * PI / angle_division as f64) as f32);
            }
        }
        angles
    }

    /// `initializeVariablesForNeighboursCoordinates`
    fn init_variables_for_neighbors(
        &mut self,
        atom: AtomId,
        visited: &BTreeSet<AtomId>,
        start_coordinates: &mut PointF,
        ordered_neighbors: &mut Vec<AtomId>,
        angles: &mut Vec<f32>,
    ) {
        if !self.atoms[atom].rings.is_empty() {
            self.init_variables_for_ring_atom(atom, start_coordinates, ordered_neighbors, angles);
            return;
        }
        ordered_neighbors.clear();
        if self.atoms[atom].neighbors.len() != 4 {
            *ordered_neighbors = self.atoms[atom].neighbors.clone();
        } else {
            let mut priorities: Vec<AtomPriority> = self.atoms[atom]
                .neighbors
                .iter()
                .map(|&n| AtomPriority { a: n, priority: 0 })
                .collect();
            self.order_atom_priorities(&mut priorities, atom);
            *ordered_neighbors = priorities.iter().map(|p| p.a).collect();
        }
        let mut start_n = 0;
        for (i, n) in ordered_neighbors.iter().enumerate() {
            if visited.contains(n) {
                *start_coordinates = self.atoms[*n].coordinates - self.atoms[atom].coordinates;
                start_n = i;
                break;
            }
        }
        ordered_neighbors.rotate_left(start_n);
        *angles = self.neighbors_angles_at_center(atom);
    }

    /// `initializeVariablesForNeighboursCoordinatesRingAtom`
    fn init_variables_for_ring_atom(
        &self,
        atom: AtomId,
        start_coordinates: &mut PointF,
        ordered_neighbors: &mut Vec<AtomId>,
        angles: &mut Vec<f32>,
    ) {
        let at = &self.atoms[atom];
        let mut ring_neighbors: Vec<(f32, AtomId)> = Vec::new();
        ordered_neighbors.clear();
        for &neigh in &at.neighbors {
            if self.share_a_ring(neigh, atom).is_some() {
                let nc = self.atoms[neigh].coordinates;
                let mut ang = (nc.y - at.coordinates.y).atan2(nc.x - at.coordinates.x);
                if ang < 0.0 {
                    ang += (2.0 * PI) as f32;
                }
                ring_neighbors.push((ang, neigh));
            } else {
                ordered_neighbors.push(neigh);
            }
        }
        ring_neighbors.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut gaps = Vec::new();
        let mut scaled_gaps = Vec::new();
        let fragment_rings = self.fragments[at.fragment.expect("fragment")]
            .rings()
            .to_vec();
        for i in 0..ring_neighbors.len() {
            let next = (i + 1) % ring_neighbors.len();
            let mut gap = ring_neighbors[next].0 - ring_neighbors[i].0;
            if gap < 0.0 {
                gap += (2.0 * PI) as f32;
            }
            gaps.push(gap);
            let mut rin = false;
            let mid = ring_neighbors[i].0 + gap * 0.5;
            let sine = (-mid).sin();
            let cosine = (-mid).cos();
            let mut p = PointF::new(BONDLENGTH * 0.1, 0.0);
            p.rotate(sine, cosine);
            let point = at.coordinates + p;
            for &r in &fragment_rings {
                if self.rings[r].is_macrocycle() {
                    continue;
                }
                if self.ring_contains_point(r, point) {
                    rin = true;
                    break;
                }
            }
            let mut scaled_gap = gap;
            if f64::from(gap) > PI {
                scaled_gap *= 10.0;
            } else if rin {
                scaled_gap *= 0.2;
            }
            scaled_gaps.push(scaled_gap);
        }
        let mut best_i = 0;
        for i in 0..gaps.len() {
            if scaled_gaps[i] > scaled_gaps[best_i] {
                best_i = i;
            }
        }
        angles.clear();
        let biggest_gap = gaps[best_i];
        *start_coordinates = self.atoms[ring_neighbors[best_i].1].coordinates - at.coordinates;
        for _ in 0..ordered_neighbors.len() {
            angles.push(-biggest_gap / (ordered_neighbors.len() + 1) as f32);
        }
    }

    /// `generateCoordinatesNeighborsOfFirstAtomInQueue`
    fn generate_coordinates_neighbors_of_first_atom_in_queue(
        &mut self,
        queue: &mut VecDeque<AtomId>,
        visited: &mut BTreeSet<AtomId>,
        fragment: FragmentId,
    ) {
        let at = queue.pop_front().expect("queue");
        let mut initial = PointF::new(-BONDLENGTH, 0.0);
        let mut ordered_neighbors = Vec::new();
        let mut angles = Vec::new();
        self.init_variables_for_neighbors(
            at,
            visited,
            &mut initial,
            &mut ordered_neighbors,
            &mut angles,
        );
        for (n, &neigh) in ordered_neighbors.iter().enumerate() {
            if visited.contains(&neigh) {
                continue;
            }
            let s = angles[n].sin();
            let c = angles[n].cos();
            initial.rotate(s, c);
            let target = self.atoms[at].coordinates + initial;
            self.set_coordinates(neigh, target);
            if self.atoms[neigh].fragment == Some(fragment) {
                visited.insert(neigh);
                queue.push_back(neigh);
            } else {
                self.atoms[neigh].coordinates_set = false;
            }
            if self.atoms[at].needs_check_for_clashes {
                self.atoms[neigh].needs_check_for_clashes = true;
                Minimizer::check_for_clashes(self, neigh);
            }
            let at_fragment = self.atoms[at].fragment.expect("fragment");
            for dof in self.fragment_dofs_of_atom(at_fragment, at) {
                if self.dofs[dof].fragment == fragment {
                    self.dof_add_atom(dof, neigh);
                }
            }
        }
        self.avoid_ze_inversions(at, visited);
        self.maybe_add_macrocycle_dof(at);
        for neighbor in self.atoms[at].neighbors.clone() {
            if self.share_a_ring(at, neighbor).is_none()
                && self.atoms[at].fragment == self.atoms[neighbor].fragment
            {
                let dof = self.new_dof(DofKind::ScaleAtoms, at);
                self.dof_add_atom(dof, neighbor);
                let f = self.atoms[at].fragment.unwrap();
                self.fragments[f].dofs.push(dof);
            }
        }
    }

    /// `maybeAddMacrocycleDOF`
    fn maybe_add_macrocycle_dof(&mut self, atom: AtomId) {
        let at = &self.atoms[atom];
        if at.rings.len() == 1 && self.rings[at.rings[0]].is_macrocycle() && at.neighbors.len() == 3
        {
            for &b in &at.bonds {
                if self.bond_is_stereo(b) && !self.bond_is_terminal(b) {
                    return;
                }
            }
            for neighbor in at.neighbors.clone() {
                if self.share_a_ring(atom, neighbor).is_none() {
                    let dof = self.new_dof(DofKind::InvertBond, atom);
                    self.dof_add_atom(dof, neighbor);
                    let f = self.atoms[atom].fragment.unwrap();
                    self.fragments[f].dofs.push(dof);
                }
            }
        }
    }

    /// `avoidZEInversions`
    fn avoid_ze_inversions(&mut self, at: AtomId, visited: &BTreeSet<AtomId>) {
        if !self.atoms[at].rings.is_empty() {
            return;
        }
        let mut double_bond = None;
        let mut to_mirror = Vec::new();
        for i in 0..self.atoms[at].bonds.len() {
            let (b, n) = (self.atoms[at].bonds[i], self.atoms[at].neighbors[i]);
            if self.bond_is_stereo(b) && visited.contains(&n) {
                double_bond = Some(b);
            } else {
                to_mirror.push(n);
            }
        }
        let Some(db) = double_bond else {
            return;
        };
        if !to_mirror.is_empty() {
            if self.bond_start_atom_cip_first_neighbor(db).is_none() {
                return;
            }
            if self.bond_end_atom_cip_first_neighbor(db).is_none() {
                return;
            }
            if !self.bond_check_stereo_chemistry(db) {
                for a in to_mirror {
                    self.mirror_coordinates(a, db);
                }
            }
        }
    }

    /// `CoordgenMinimizer::avoidInternalClashes`
    fn avoid_internal_clashes(&mut self, fragment: FragmentId) {
        let fragment_atoms = self.fragments[fragment].atoms().to_vec();
        for &a in &fragment_atoms {
            if self.atoms[a].neighbors.len() != 1
                || self.atoms[a].needs_check_for_clashes
                || self.atoms[a].fixed
            {
                continue;
            }
            if !self.fragment_dofs_of_atom(fragment, a).is_empty() {
                continue;
            }
            for &a2 in &fragment_atoms {
                if a == a2
                    || !self.fragment_dofs_of_atom(fragment, a2).is_empty()
                    || self.bond_to(a, a2).is_some()
                {
                    continue;
                }
                let dx = self.atoms[a2].coordinates.x - self.atoms[a].coordinates.x;
                if dx > BONDLENGTH * 0.5 || dx < -BONDLENGTH * 0.5 {
                    continue;
                }
                let dy = self.atoms[a2].coordinates.y - self.atoms[a].coordinates.y;
                if dy > BONDLENGTH * 0.5 || dy < -BONDLENGTH * 0.5 {
                    continue;
                }
                let square_d = dx * dx + dy * dy;
                if square_d > BONDLENGTH * 0.5 * BONDLENGTH * 0.5 {
                    continue;
                }
                let n0 = self.atoms[a].neighbors[0];
                let mut vec = self.atoms[a].coordinates - self.atoms[n0].coordinates;
                vec = vec * 0.3;
                self.atoms[a].coordinates -= vec;
                if self.atoms[a2].neighbors.len() == 1 {
                    self.atoms[a2].coordinates += vec;
                    self.atoms[a2].coordinates.round(2);
                }
            }
        }
    }

    /// `buildFragment`
    fn build_fragment(&mut self, fragment: FragmentId) {
        self.build_rings(fragment);
        self.build_non_ring_atoms(fragment);
        self.avoid_internal_clashes(fragment);
        // fallbackIfNanCoordinates
        let atoms = self.fragments[fragment].atoms().to_vec();
        let has_nan = atoms
            .iter()
            .any(|&a| self.atoms[a].coordinates.x.is_nan() || self.atoms[a].coordinates.y.is_nan());
        if has_nan && atoms.iter().all(|&a| self.has_valid_3d_coordinates(a)) {
            for &a in &atoms {
                let (x, y) = (self.atoms[a].x3d, self.atoms[a].y3d);
                self.set_coordinates(a, PointF::new(x * 35.0, -y * 35.0));
            }
        }
        if self.fragments[fragment].parent.is_none() && self.fragments[fragment].constrained {
            self.rotate_main_fragment(fragment);
        }
        if self.fragments[fragment].fixed {
            self.fragment_set_all_coordinates_to_template(fragment);
        }
    }

    /// `initializeFusedRingInformation`
    fn initialize_fused_ring_information(&mut self, fragment: FragmentId) {
        if self.fragments[fragment].rings().len() < 2 {
            return;
        }
        let fragment_atoms = self.fragments[fragment].atoms().to_vec();
        for &atom in &fragment_atoms {
            let ar = self.atoms[atom].rings.clone();
            if ar.len() > 1 {
                for i in 0..ar.len() {
                    for j in i + 1..ar.len() {
                        let (ri, rj) = (ar[i], ar[j]);
                        let mut found = false;
                        for r in 0..self.rings[ri].fused_with.len() {
                            if self.rings[ri].fused_with[r] == rj {
                                self.rings[ri].fusion_atoms[r].push(atom);
                                for rr in 0..self.rings[rj].fused_with.len() {
                                    if self.rings[rj].fused_with[rr] == ri {
                                        self.rings[rj].fusion_atoms[rr].push(atom);
                                        break;
                                    }
                                }
                                found = true;
                                break;
                            }
                        }
                        if !found {
                            self.rings[ri].fused_with.push(rj);
                            self.rings[rj].fused_with.push(ri);
                            self.rings[ri].fusion_atoms.push(vec![atom]);
                            self.rings[rj].fusion_atoms.push(vec![atom]);
                        }
                    }
                }
            }
        }
        for b in self.fragments[fragment].bonds().to_vec() {
            let bd = &self.bonds[b];
            let (s, e) = (bd.start_atom, bd.end_atom);
            if bd.bond_order != 1
                && !self.atoms[s].rings.is_empty()
                && !self.atoms[e].rings.is_empty()
                && self.share_a_ring(s, e).is_none()
            {
                for r in self.atoms[s].rings.clone() {
                    for r2 in self.atoms[e].rings.clone() {
                        self.rings[r].fused_with.push(r2);
                        self.rings[r2].fused_with.push(r);
                        self.rings[r].fusion_atoms.push(Vec::new());
                        self.rings[r2].fusion_atoms.push(Vec::new());
                        self.rings[r].fusion_bonds.push(b);
                        self.rings[r2].fusion_bonds.push(b);
                    }
                }
            }
        }
        let fragment_rings = self.fragments[fragment].rings().to_vec();
        for &r in &fragment_rings {
            for i in 0..self.rings[r].fused_with.len() {
                let fusion_atoms = self.rings[r].fusion_atoms[i].clone();
                if fusion_atoms.len() <= 2 {
                    continue;
                }
                let mut start_atom = None;
                for &a in &fusion_atoms {
                    let counter = self.atoms[a]
                        .neighbors
                        .iter()
                        .filter(|n| fusion_atoms.contains(n))
                        .count();
                    if counter == 1 {
                        start_atom = Some(a);
                        break;
                    }
                }
                let start = start_atom.expect("fusion chain start");
                self.rings[r].fusion_atoms[i] = self.order_chain_of_atoms(&fusion_atoms, start);
            }
        }
        for &r in &fragment_rings {
            for ii in 0..self.rings[r].fused_with.len() {
                if self.rings[r].fusion_atoms[ii].len() > 2 {
                    for a in self.rings[r].fusion_atoms[ii].clone() {
                        let mut found = false;
                        for &n in &self.atoms[a].neighbors {
                            let mut nn = 0;
                            for &nr in &self.atoms[n].rings {
                                if nr == r || nr == self.rings[r].fused_with[ii] {
                                    nn += 1;
                                }
                            }
                            if nn == 1 {
                                found = true;
                                break;
                            }
                        }
                        if !found {
                            self.atoms[a].is_shared_and_inner = true;
                        }
                    }
                }
            }
        }
    }

    /// `simplifyRingSystem`: (side rings in the order they are pushed on the stack, central rings)
    fn simplify_ring_system(&mut self, all_rings: &[RingId]) -> (Vec<RingId>, Vec<RingId>) {
        let mut chosen: BTreeSet<RingId> = BTreeSet::new();
        let mut side_rings = Vec::new();
        let mut found = true;
        while found {
            found = false;
            for &r in all_rings {
                if chosen.contains(&r) || self.rings[r].is_macrocycle() {
                    continue;
                }
                let mut n = 0;
                for (counter, &fused) in self.rings[r].fused_with.iter().enumerate() {
                    if !chosen.contains(&fused) {
                        n += 1;
                        if self.rings[fused].is_macrocycle() {
                            n += 1;
                        }
                        let fa = self.rings[r].fusion_atoms[counter].len();
                        if fa > 3 {
                            n += 1;
                        }
                        if fa == 3 && self.rings[r].size() == 4 && self.rings[fused].size() == 4 {
                            n += 1;
                        }
                    }
                }
                if n == 1 {
                    side_rings.push(r);
                    found = true;
                    chosen.insert(r);
                    self.rings[r].side = true;
                }
            }
        }
        let central = all_rings
            .iter()
            .copied()
            .filter(|r| !chosen.contains(r))
            .collect();
        (side_rings, central)
    }

    /// `buildRings`
    fn build_rings(&mut self, fragment: FragmentId) {
        if self.fragments[fragment].rings().is_empty() {
            return;
        }
        self.initialize_fused_ring_information(fragment);
        let fragment_rings = self.fragments[fragment].rings().to_vec();
        if fragment_rings.len() == 1 {
            self.build_ring(fragment_rings[0]);
        } else {
            let (side_rings, central) = self.simplify_ring_system(&fragment_rings);
            self.generate_coordinates_central_rings(central);
            // the stack pops the last pushed ring first
            for r in side_rings.into_iter().rev() {
                self.build_ring(r);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{build, list};
    use super::*;
    use std::fmt::Write;

    fn run(input: &str) -> String {
        let (mut model, mol) = build(input);
        model.force_update_struct(mol);
        model.split_into_fragments(mol);
        let m = &model.molecules[mol];
        let (atoms, rings, frags) = (m.atoms.clone(), m.rings.clone(), m.fragments.clone());
        if rings.iter().any(|&r| model.rings[r].is_macrocycle())
            || frags
                .iter()
                .any(|&f| model.new_score_planarity(model.fragments[f].rings()) >= 1000.0)
        {
            return "skip\n".to_string();
        }
        for &f in &frags {
            model.initialize_coordinates(f);
        }
        let ai = |a: AtomId| atoms.iter().position(|&x| x == a).unwrap();
        let mut out = String::new();
        for (i, &a) in atoms.iter().enumerate() {
            let at = &model.atoms[a];
            writeln!(
                out,
                "a {i} {:08x} {:08x} rigid {} inner {} set {}",
                at.coordinates.x.to_bits(),
                at.coordinates.y.to_bits(),
                at.rigid as i32,
                at.is_shared_and_inner as i32,
                at.coordinates_set as i32
            )
            .unwrap();
        }
        for (i, &r) in rings.iter().enumerate() {
            let ring = &model.rings[r];
            write!(
                out,
                "r {i} gen {} side {} fused",
                ring.coordinates_generated as i32, ring.side as i32
            )
            .unwrap();
            for (k, &f) in ring.fused_with.iter().enumerate() {
                write!(out, " {}:", rings.iter().position(|&x| x == f).unwrap()).unwrap();
                for &a in &ring.fusion_atoms[k] {
                    write!(out, "{},", ai(a)).unwrap();
                }
            }
            out.push('\n');
        }
        for (i, &f) in frags.iter().enumerate() {
            let fr = &model.fragments[f];
            write!(out, "f {i} templated {} dofs", fr.is_templated as i32).unwrap();
            for &d in &fr.dofs {
                let dof = &model.dofs[d];
                let tier = dof.kind.tier();
                write!(out, " [{tier}").unwrap();
                for &a in &dof.atoms {
                    write!(out, " {}", ai(a)).unwrap();
                }
                out.push(']');
            }
            out.push_str(" coords");
            // C++ orders the map by pointer; compare by atom index
            let mut c: Vec<(usize, PointF)> =
                fr.coordinates.iter().map(|(&a, &p)| (ai(a), p)).collect();
            c.sort_by_key(|x| x.0);
            for (a, p) in c {
                write!(out, " {a}:{:08x}:{:08x}", p.x.to_bits(), p.y.to_bits()).unwrap();
            }
            out.push('\n');
        }
        let _ = list;
        out
    }

    /// equal up to the last bits of the platform's `sin`/`cos` (the 8 digit hex words are float coordinates)
    fn same(got: &str, want: &str) -> bool {
        let tokens = |s: &str| -> Vec<String> {
            s.split(|c: char| c.is_whitespace() || c == ':')
                .map(String::from)
                .collect()
        };
        let (g, w) = (tokens(got), tokens(want));
        g.len() == w.len()
            && g.iter().zip(&w).all(|(a, b)| {
                let hex = |t: &str| t.len() == 8 && t.chars().all(|c| c.is_ascii_hexdigit());
                if a == b {
                    return true;
                }
                if hex(a) && hex(b) {
                    let f = |t: &str| f32::from_bits(u32::from_str_radix(t, 16).unwrap());
                    return (f(a) - f(b)).abs() < 1e-3;
                }
                false
            })
    }

    /// initial fragment coordinates must be the ones CoordgenFragmentBuilder makes (`tests/data/coordgen_fbuild.txt`)
    #[test]
    fn fragment_builder_matches_cpp() {
        let text = include_str!("../../../../../tests/data/coordgen_fbuild.txt");
        let mut bad = Vec::new();
        let mut n = 0;
        for rec in text.split("$$\n").filter(|r| !r.is_empty()) {
            let end = rec.find("END\n").unwrap() + 4;
            let (input, want) = rec.split_at(end);
            let got = run(input);
            n += 1;
            if !same(&got, want) {
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
        assert!(n > 60, "{n}");
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }
}
