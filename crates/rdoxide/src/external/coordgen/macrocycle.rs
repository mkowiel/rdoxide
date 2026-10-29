//! `CoordgenMacrocycleBuilder`: the atoms of a macrocycle are placed on the contour of a polyomino (a shape made of
//! hexagons) that satisfies the double bond and fused ring constraints; when none fits, the ring is opened at one bond,
//! laid out as a chain and pulled closed again.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use super::driver::Sketcher;
use super::maths::{BONDLENGTH, MACROCYCLE, PointF, same_side};
use super::minimizer::{Interaction, Stretch};
use super::model::{AtomId, BondId, Model, RingId};
use super::polyomino::{HexCoords, Polyomino, VertexCoords};

const MAX_MACROCYCLES: i32 = 40;
const PATH_FAILED: i32 = -1000;
const HETEROATOM_RESTRAINT: i32 = 1;
const SUBSTITUTED_ATOM_RESTRAINT: i32 = 10;
const SUBSTITUENTS_TO_SAME_VERTEX_RESTRAINT: i32 = 200;
const SUBSTITUENT_CLASH_WITH_PATH_RESTRAINT: i32 = 400;
const SQRT3HALF: f64 = 0.866_025_403_784_438_6;

/// `doubleBondConstraint`
struct DoubleBondConstraint {
    trans: bool,
    atom1: usize,
    atom2: usize,
    previous_atom: usize,
    following_atom: usize,
}

/// `ringConstraint`
struct RingConstraint {
    atom: usize,
    ring: RingId,
    force_outside: bool,
}

/// `pathRestraints`
struct PathRestraints {
    hetero_atoms: Vec<usize>,
    substituted_atoms: Vec<(usize, i32)>,
}

/// `pathConstraints`
struct PathConstraints {
    double_bonds: Vec<DoubleBondConstraint>,
    ring_constraints: Vec<RingConstraint>,
}

/// `coordsOfVertex`
fn coords_of_vertex(v: VertexCoords) -> PointF {
    let bl = f64::from(BONDLENGTH);
    PointF::new(
        (bl * SQRT3HALF * f64::from(v.x) - bl * SQRT3HALF * f64::from(v.z)) as f32,
        (-bl * 0.5 * f64::from(v.x) + bl * f64::from(v.y) + -bl * 0.5 * f64::from(v.z)) as f32,
    )
}

/// `acceptableShapeScore`
fn acceptable_shape_score(number_of_atoms: i32) -> i32 {
    if number_of_atoms < 10 {
        0
    } else {
        number_of_atoms * SUBSTITUTED_ATOM_RESTRAINT / 2
    }
}

/// `removeDuplicates`
fn remove_duplicates(pols: &[Polyomino]) -> Vec<Polyomino> {
    let mut out: Vec<Polyomino> = Vec::new();
    for pol in pols {
        if !out.iter().any(|j| pol.is_the_same_as(j)) {
            out.push(pol.clone());
        }
    }
    out
}

/// `listOfEquivalent`: polyominoes with the same number of vertices from removing hexagons with 3 neighbors
fn list_of_equivalent(p: &Polyomino) -> Vec<Polyomino> {
    let mut out = Vec::new();
    let pentagon_vs = p.pentagon_vertices.len();
    for &c in &p.list {
        if p.is_equivalent_without(c) {
            let mut new_p = p.clone();
            new_p.pentagon_vertices.clear();
            new_p.remove_hex(c);
            for _ in 0..pentagon_vs {
                new_p.mark_one_vertex_as_pentagon();
            }
            out.push(new_p);
        }
    }
    out
}

/// `buildSquaredShapes`
fn build_squared_shapes(mut tot_vertices: i32) -> Vec<Polyomino> {
    // new shapes go to the front: rounder shapes (more space inside) are evaluated first
    let mut out: Vec<Polyomino> = Vec::new();
    let mut pentagon = false;
    if tot_vertices % 2 > 0 {
        pentagon = true;
        tot_vertices += 1;
    }
    if tot_vertices % 4 == 0 {
        if tot_vertices >= 12 {
            let xandy = tot_vertices / 4;
            let mid = xandy / 2;
            for i in 1..mid {
                let i2 = xandy - i;
                if i2 % 2 == 0 && i > 1 {
                    let mut p = Polyomino::new();
                    p.build_ragged_bigger_box_shape(i, i2, pentagon);
                    out.insert(0, p);
                }
                if i % 2 == 0 && i2 > 1 {
                    let mut p = Polyomino::new();
                    p.build_ragged_bigger_box_shape(i2, i, pentagon);
                    out.insert(0, p);
                }
            }
        }
    } else {
        let xandy = (tot_vertices + 2) / 4;
        let mid = xandy / 2;
        for i in 1..mid + 1 {
            let mut p = Polyomino::new();
            p.build_skewed_box_shape(i, xandy - i, pentagon);
            out.insert(0, p);
            let i2 = xandy - i;
            if i < 2 || i2 < 2 {
                continue;
            }
            let mut p = Polyomino::new();
            p.build_ragged_box_shape(i, i2, pentagon);
            out.insert(0, p);
            let mut p = Polyomino::new();
            p.build_ragged_box_shape(i2, i, pentagon);
            out.insert(0, p);
            if i2 % 2 > 0 {
                let mut p = Polyomino::new();
                p.build_ragged_bigger_box_shape(i, i2, pentagon);
                out.insert(0, p);
            }
            if i % 2 > 0 {
                let mut p = Polyomino::new();
                p.build_ragged_bigger_box_shape(i2, i, pentagon);
                out.insert(0, p);
            }
            if i > 2 && i2 % 2 > 0 {
                let mut p = Polyomino::new();
                p.build_ragged_smaller_box_shape(i, i2, pentagon);
                out.insert(0, p);
            }
            if i2 > 2 && i % 2 > 0 {
                let mut p = Polyomino::new();
                p.build_ragged_smaller_box_shape(i2, i, pentagon);
                out.insert(0, p);
            }
        }
    }
    out
}

/// `getLowestPeriod`
fn get_lowest_period(neighbors: &[i32]) -> usize {
    for period in 1..neighbors.len() {
        let mut has_difference = false;
        for i in 0..neighbors.len() {
            let i2 = (i + period) % neighbors.len();
            if neighbors[i] != neighbors[i2] {
                has_difference = true;
                break;
            }
        }
        if !has_difference {
            return period;
        }
    }
    neighbors.len()
}

/// `checkDoubleBoundConstraints`
fn check_double_bond_constraints(
    constraints: &[DoubleBondConstraint],
    vertices: &[VertexCoords],
    start_i: usize,
) -> bool {
    let n = vertices.len();
    for c in constraints {
        let p1 = coords_of_vertex(vertices[(start_i + c.previous_atom) % n]);
        let p2 = coords_of_vertex(vertices[(start_i + c.atom1) % n]);
        let mut counter = start_i + c.atom2;
        if counter >= n {
            counter -= n;
        }
        let p3 = coords_of_vertex(vertices[counter]);
        let mut counter = start_i + c.following_atom;
        if counter >= n {
            counter -= n;
        }
        let p4 = coords_of_vertex(vertices[counter]);
        if same_side(p1, p4, p2, p3) == c.trans {
            return false;
        }
    }
    true
}

/// `checkRingConstraints`: all shared atoms with a fused small ring must map to the same hexagon
fn check_ring_constraints(
    constraints: &[RingConstraint],
    p: &Polyomino,
    path: &[VertexCoords],
    neighbor_ns: &[i32],
    start_i: usize,
) -> bool {
    let mut allowed_hexs: BTreeMap<RingId, Vec<HexCoords>> = BTreeMap::new();
    for rc in constraints {
        let counter = (rc.atom + start_i) % path.len();
        if rc.force_outside && neighbor_ns[counter] != 1 {
            return false;
        }
        let new_pos = p.free_vertex_neighbor_positions(path[counter]);
        let old_pos = allowed_hexs.entry(rc.ring).or_default().clone();
        let next_pos: Vec<HexCoords> = if old_pos.is_empty() {
            new_pos
        } else {
            new_pos
                .into_iter()
                .filter(|c| old_pos.contains(c))
                .collect()
        };
        if next_pos.is_empty() {
            return false;
        }
        allowed_hexs.insert(rc.ring, next_pos);
    }
    true
}

/// `scorePathRestraints`
fn score_path_restraints(
    pr: &PathRestraints,
    p: &Polyomino,
    path: &[VertexCoords],
    neighbor_ns: &[i32],
    start_i: usize,
) -> i32 {
    let mut score = 0;
    let mut used: BTreeSet<(i32, i32, i32)> = BTreeSet::new();
    for &hetero in &pr.hetero_atoms {
        let counter = (hetero + start_i) % neighbor_ns.len();
        if neighbor_ns[counter] == 1 {
            score -= HETEROATOM_RESTRAINT; // heteroatom placed towards outside
        }
    }
    for (i, &(atom, children)) in pr.substituted_atoms.iter().enumerate() {
        let counter = (atom + start_i) % neighbor_ns.len();
        if neighbor_ns[counter] == 2 {
            score -= SUBSTITUTED_ATOM_RESTRAINT * children; // substituted atoms placed towards inside
            // upstream indexes the path with the position in the list of substituted atoms
            let s = p.coordinates_of_substituent(path[i]);
            let key = (s.x, s.y, s.z);
            if used.contains(&key) {
                score -= SUBSTITUENTS_TO_SAME_VERTEX_RESTRAINT;
            }
            if path.contains(&s) {
                score -= SUBSTITUENT_CLASH_WITH_PATH_RESTRAINT;
            }
            used.insert(key);
        }
    }
    score
}

/// `scorePath`
fn score_path(
    p: &Polyomino,
    path: &[VertexCoords],
    neighbor_ns: &[i32],
    start_i: usize,
    pc: &PathConstraints,
    pr: &PathRestraints,
) -> i32 {
    if !check_ring_constraints(&pc.ring_constraints, p, path, neighbor_ns, start_i)
        || !check_double_bond_constraints(&pc.double_bonds, path, start_i)
    {
        return PATH_FAILED;
    }
    score_path_restraints(pr, p, path, neighbor_ns, start_i)
}

/// `matchPolyomino`: (matched, best start, best score)
fn match_polyomino(p: &Polyomino, pc: &PathConstraints, pr: &PathRestraints) -> (bool, usize, i32) {
    let path = p.get_path();
    let neighbor_ns: Vec<i32> = path
        .iter()
        .map(|&v| p.hexagons_at_vertex(v) as i32)
        .collect();
    let mut best_start = 0;
    let mut best_score = PATH_FAILED;
    let mut start_i = 0;
    while start_i < get_lowest_period(&neighbor_ns) {
        let score = score_path(p, &path, &neighbor_ns, start_i, pc, pr);
        if score > best_score {
            best_score = score;
            best_start = start_i;
            if best_score == 0 {
                break;
            }
        }
        start_i += 1;
    }
    (best_score > PATH_FAILED, best_start, best_score)
}

/// `matchPolyominoes`: (matched, best polyomino, best score, best start)
fn match_polyominoes(
    pols: &[Polyomino],
    pc: &PathConstraints,
    pr: &PathRestraints,
    checked_macrocycles: &mut i32,
) -> (bool, usize, i32, usize) {
    let mut best_start = 0;
    let mut best_p = 0;
    let mut best_score = PATH_FAILED;
    let mut matched = false;
    for (i, pol) in pols.iter().enumerate() {
        let (m, start_of_this, score_of_this) = match_polyomino(pol, pc, pr);
        if m {
            matched = true;
            if score_of_this > best_score {
                best_score = score_of_this;
                best_start = start_of_this;
                best_p = i;
                if best_score == 0 {
                    return (true, best_p, best_score, best_start);
                }
            }
        }
        let c = *checked_macrocycles;
        *checked_macrocycles += 1;
        if c > MAX_MACROCYCLES {
            break;
        }
    }
    (matched, best_p, best_score, best_start)
}

impl Model {
    /// `sketcherMinimizer::getBond`
    fn bond_between(&self, a1: AtomId, a2: AtomId) -> Option<BondId> {
        self.atoms[a1].bonds.iter().copied().find(|&b| {
            let bd = &self.bonds[b];
            (bd.start_atom == a1 && bd.end_atom == a2) || (bd.start_atom == a2 && bd.end_atom == a1)
        })
    }

    /// `getRingConstraints`
    fn get_ring_constraints(&self, atoms: &[AtomId]) -> Vec<RingConstraint> {
        let mut out = Vec::new();
        for (i, &a) in atoms.iter().enumerate() {
            if self.atoms[a].rings.len() > 1 {
                for &r in &self.atoms[a].rings {
                    // this excludes the current cycle and all fused macrocycles
                    if self.rings[r].atoms.len() < MACROCYCLE {
                        let mut force_outside = false;
                        for &n in &self.atoms[a].neighbors {
                            if !atoms.contains(&n) {
                                if self.rings[r].contains_atom(n) {
                                    force_outside = true;
                                }
                                break;
                            }
                        }
                        out.push(RingConstraint {
                            atom: i,
                            ring: r,
                            force_outside,
                        });
                    }
                }
            }
        }
        out
    }

    /// `getDoubleBondConstraints`
    fn get_double_bond_constraints(&self, atoms: &[AtomId]) -> Vec<DoubleBondConstraint> {
        let mut out = Vec::new();
        if atoms.len() >= MACROCYCLE {
            let n = atoms.len();
            for i in 0..n {
                let mut index = i;
                let mut index2 = (i + 1) % n;
                let Some(b) = self.bond_between(atoms[i], atoms[index2]) else {
                    break;
                };
                let bd = &self.bonds[b];
                if bd.bond_order != 2 {
                    continue;
                }
                let mut small_ring_bond = false;
                if bd.rings.len() > 1 {
                    for &ring in &bd.rings {
                        if self.rings[ring].atoms.len() < MACROCYCLE {
                            small_ring_bond = true;
                            break;
                        }
                    }
                }
                if small_ring_bond {
                    continue;
                }
                let mut previous_i = (i + n - 1) % n;
                let mut following_i = (i + 2) % n;
                // is the atom trans in the ring? (isZ stores the absolute chirality)
                let mut is_trans = !bd.is_z;
                if bd.start_atom != atoms[i] {
                    // atoms[i] could be the end atom of b
                    std::mem::swap(&mut previous_i, &mut following_i);
                    index = index2;
                    index2 = i;
                }
                if self.bond_start_atom_cip_first_neighbor(b) != Some(atoms[previous_i]) {
                    is_trans = !is_trans;
                }
                if self.bond_end_atom_cip_first_neighbor(b) != Some(atoms[following_i]) {
                    is_trans = !is_trans;
                }
                out.push(DoubleBondConstraint {
                    trans: is_trans,
                    atom1: index,
                    atom2: index2,
                    previous_atom: previous_i,
                    following_atom: following_i,
                });
            }
        }
        out
    }

    /// `getNumberOfChildren`
    fn get_number_of_children(&self, a: AtomId, parent: AtomId) -> i32 {
        let mut n = 0;
        let mut visited: BTreeSet<AtomId> = BTreeSet::new();
        let mut q = VecDeque::new();
        q.push_back(a);
        visited.insert(parent);
        while let Some(this_a) = q.pop_front() {
            visited.insert(this_a);
            n += 1;
            for &nb in &self.atoms[this_a].neighbors {
                if visited.contains(&nb) {
                    continue;
                }
                q.push_back(nb);
            }
        }
        n
    }

    /// `getPathRestraints`
    fn get_path_restraints(&self, atoms: &[AtomId]) -> PathRestraints {
        let mut hetero_atoms = Vec::new();
        let mut substituted_atoms = Vec::new();
        let len = atoms.len();
        for i in 0..len {
            if self.atoms[atoms[i]].atomic_number != 6 {
                hetero_atoms.push(i);
            }
            if self.atoms[atoms[i]].neighbors.len() != 2 {
                let mut tot_n = 0;
                let prev_i = (i + len - 1) % len;
                let post_i = (i + 1) % len;
                for &n in &self.atoms[atoms[i]].neighbors {
                    if n == atoms[prev_i] || n == atoms[post_i] {
                        continue;
                    }
                    tot_n += self.get_number_of_children(n, atoms[i]);
                }
                substituted_atoms.push((i, tot_n));
            }
        }
        PathRestraints {
            hetero_atoms,
            substituted_atoms,
        }
    }

    /// `findBondToOpen`
    pub fn find_bond_to_open(&self, ring: RingId) -> Option<BondId> {
        let mut best_bond = None;
        let mut best_score = 0.0f32;
        for &bond in &self.rings[ring].bonds {
            let bd = &self.bonds[bond];
            let mut score = 0.0f32;
            if self.rings[ring].is_macrocycle() {
                if bd.bond_order != 1 {
                    continue;
                }
                let next_to_multiple_bond = [bd.start_atom, bd.end_atom].iter().any(|&a| {
                    self.atoms[a]
                        .bonds
                        .iter()
                        .any(|&ob| self.bond_is_stereo(ob))
                });
                if next_to_multiple_bond {
                    continue;
                }
            }
            score += (bd.rings.len() * 10) as f32;
            score += self.atoms[bd.start_atom].neighbors.len() as f32;
            score += self.atoms[bd.end_atom].neighbors.len() as f32;
            score /= bd.crossing_bond_penalty_multiplier;
            if best_bond.is_none() || score < best_score {
                best_score = score;
                best_bond = Some(bond);
            }
        }
        best_bond
    }

    /// `openCycleAndGenerateCoords`: the coordinates of the whole molecule from opening one bond of the macrocycle
    pub fn open_cycle_and_generate_coords(&mut self, ring: RingId) -> bool {
        let Some(bond_to_break) = self.find_bond_to_open(ring) else {
            return false;
        };
        let mut nm = Model::new();
        nm.force_open_macrocycles = true;
        let mut sk = Sketcher::new(self.precision);
        sk.minimizer.skip_minimization = true;
        let min_mol = nm.add_molecule();
        let first_atom = self.rings[ring].atoms[0];
        let mol = self.atoms[first_atom].molecule.expect("molecule");
        let atoms = self.molecules[mol].atoms.clone();
        let mut atom_map: BTreeMap<AtomId, AtomId> = BTreeMap::new();
        for &at in &atoms {
            let na = nm.add_atom(min_mol);
            atom_map.insert(at, na);
            let src = &self.atoms[at];
            let dst = &mut nm.atoms[na];
            dst.template_coordinates = src.coordinates;
            dst.coordinates = src.coordinates;
            dst.atomic_number = src.atomic_number;
            dst.general_use_n = src.general_use_n;
            dst.implicit_hs = src.implicit_hs;
        }
        for &bo in &self.molecules[mol].bonds.clone() {
            if bo == bond_to_break {
                continue;
            }
            let bd = &self.bonds[bo];
            let nb = nm.add_bond(min_mol, atom_map[&bd.start_atom], atom_map[&bd.end_atom]);
            nm.bonds[nb].bond_order = bd.bond_order;
            nm.bonds[nb].is_z = bd.is_z;
            nm.bonds[nb].ignore_ze = bd.ignore_ze;
        }
        sk.initialize(&mut nm, min_mol);
        sk.find_fragments(&mut nm);
        for &mm in &sk.molecules.clone() {
            nm.build_molecule_from_fragments(mm, true);
        }
        let bb = &self.bonds[bond_to_break];
        let (bs, be) = (atom_map[&bb.start_atom], atom_map[&bb.end_atom]);
        let broken = nm.add_bond(min_mol, bs, be);
        nm.bonds[broken].bond_order = bb.bond_order;
        nm.bonds[broken].is_z = bb.is_z;
        nm.force_update_struct(min_mol);
        let extra = vec![Interaction::Stretch(Stretch::new(bs, be))];
        sk.minimizer
            .avoid_clashes_of_molecule_with(&mut nm, min_mol, extra);

        for &atom in &atoms {
            let other = atom_map[&atom];
            if nm.atoms[other].rigid {
                self.atoms[atom].rigid = true;
            }
            self.set_coordinates(atom, nm.atoms[other].coordinates);
        }
        let fragment = self.atoms[first_atom].fragment.expect("fragment");
        for r in self.fragments[fragment].rings().to_vec() {
            self.rings[r].coordinates_generated = true;
        }
        for br in self.bonds[bond_to_break].rings.clone() {
            if !self.rings[br].is_macrocycle() && br != ring {
                self.rings[br].coordinates_generated = false;
            }
        }
        true
    }

    /// `newMacrocycle`: coordinates for the atoms of a ring of at least `MACROCYCLE` atoms
    pub fn new_macrocycle(&mut self, ring: RingId, atoms: &[AtomId]) -> Vec<PointF> {
        let natoms = atoms.len() as i32;
        let mut p = Polyomino::new();
        p.build_with_vertices_n(natoms);
        let mut pols = vec![p];
        let square_pols = build_squared_shapes(natoms);
        for sp in &square_pols {
            assert_eq!(sp.get_path().len(), atoms.len());
        }
        pols.extend(square_pols);
        pols = remove_duplicates(&pols);
        let mut found = false;
        let pr = self.get_path_restraints(atoms);
        let pc = PathConstraints {
            double_bonds: self.get_double_bond_constraints(atoms),
            ring_constraints: self.get_ring_constraints(atoms),
        };
        let mut score_of_chosen = PATH_FAILED;
        let mut start_of_chosen = 0;
        let mut chosen_p = pols[0].clone();
        let acceptable_score = acceptable_shape_score(natoms);
        let mut checked_macrocycles = 0;
        if !self.force_open_macrocycles {
            loop {
                let (f, best_p, best_score, best_start) =
                    match_polyominoes(&pols, &pc, &pr, &mut checked_macrocycles);
                found = f;
                if best_score > score_of_chosen {
                    start_of_chosen = best_start;
                    score_of_chosen = best_score;
                    chosen_p = pols[best_p].clone();
                    if best_score > acceptable_score {
                        break;
                    }
                }
                if checked_macrocycles > MAX_MACROCYCLES {
                    break;
                }
                pols = pols.iter().flat_map(list_of_equivalent).collect();
                pols = remove_duplicates(&pols);
                if pols.is_empty() {
                    break;
                }
            }
        }
        let mol = self.atoms[atoms[0]].molecule.expect("molecule");
        if found {
            let path = chosen_p.get_path();
            for (n, &a) in atoms.iter().enumerate() {
                let v = path[(n + start_of_chosen) % path.len()];
                if !self.atoms[a].coordinates_set {
                    self.set_coordinates(a, coords_of_vertex(v));
                }
            }
            if !chosen_p.pentagon_vertices.is_empty() {
                self.molecules[mol].require_minimization = true;
            }
        } else {
            // could not find a shape: fallback methods
            if !self.open_cycle_and_generate_coords(ring) {
                let coords = Self::polygon_coordinates(atoms.len());
                for (i, &a) in atoms.iter().enumerate() {
                    self.set_coordinates(a, coords[i]);
                }
            }
            self.molecules[mol].require_minimization = true;
        }
        atoms.iter().map(|&a| self.atoms[a].coordinates).collect()
    }
}
