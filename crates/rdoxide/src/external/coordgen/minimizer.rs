//! `CoordgenMinimizer` and the force field interactions (`sketcherMinimizer*Interaction`): stretch, bend, clash, constraint
//! and E/Z constraint terms, the interactions of a molecule, and the steepest descent `run`.

#![allow(clippy::manual_clamp)]

use super::maths::{
    BONDLENGTH, PointF, SKETCHER_EPSILON, project_point_on_line, same_side, squared_distance,
    squared_distance_point_segment, unsigned_angle,
};
use super::model::{AtomId, Model, MolId};

const FORCE_MULTIPLIER: f32 = 0.3;
const ITERATION_HISTORY_SIZE: usize = 100;
const MAX_NET_ENERGY_CHANGE: f32 = 20.0;
const CONSTRAINT_SCALE: f32 = 0.5;

/// `sketcherMinimizerStretchInteraction`
#[derive(Debug, Clone)]
pub struct Stretch {
    pub atom1: AtomId,
    pub atom2: AtomId,
    pub k: f32,
    pub rest_v: f32,
    pub tolerance: f32,
}

/// `sketcherMinimizerBendInteraction`
#[derive(Debug, Clone)]
pub struct Bend {
    pub atom1: AtomId,
    pub atom2: AtomId,
    pub atom3: AtomId,
    pub k: f32,
    pub k2: f32,
    pub rest_v: f32,
    pub is_ring: bool,
}

/// `sketcherMinimizerClashInteraction`
#[derive(Debug, Clone)]
pub struct Clash {
    pub atom1: AtomId,
    pub atom2: AtomId,
    pub atom3: AtomId,
    pub k: f32,
    pub k2: f32,
    pub rest_v: f32,
    squared_distance: f32,
}

/// `sketcherMinimizerConstraintInteraction`
#[derive(Debug, Clone)]
pub struct Constraint {
    pub atom1: AtomId,
    pub k: f32,
    pub origin: PointF,
}

/// `sketcherMinimizerEZConstrainInteraction`
#[derive(Debug, Clone)]
pub struct EzConstrain {
    pub atom1: AtomId,
    pub atom2: AtomId,
    pub atom3: AtomId,
    pub atom4: AtomId,
    pub is_z: bool,
    pub force_movement: bool,
}

#[derive(Debug, Clone)]
pub enum Interaction {
    Stretch(Stretch),
    Bend(Bend),
    Clash(Clash),
    Constraint(Constraint),
    EzConstrain(EzConstrain),
}

impl Stretch {
    pub fn new(atom1: AtomId, atom2: AtomId) -> Self {
        Stretch {
            atom1,
            atom2,
            k: 1.0,
            rest_v: 50.0,
            tolerance: 0.0,
        }
    }

    fn energy(&self, m: &Model, e: &mut f32) {
        let l = m.atoms[self.atom1].coordinates - m.atoms[self.atom2].coordinates;
        let dr = (l.x * l.x + l.y * l.y).sqrt() - self.rest_v;
        *e += 0.5 * self.k * dr * dr;
    }

    fn score(&self, m: &mut Model, total: &mut f32) {
        self.energy(m, total);
        let mut l = m.atoms[self.atom1].coordinates - m.atoms[self.atom2].coordinates;
        let len = l.length();
        let dr = if len < self.rest_v - self.tolerance {
            self.rest_v - self.tolerance - len
        } else if len > self.rest_v + self.tolerance {
            self.rest_v + self.tolerance - len
        } else {
            return;
        };
        let short_bond_threshold = self.rest_v * 0.4;
        let mut penalty_for_very_short_bonds = short_bond_threshold - len;
        if penalty_for_very_short_bonds < 0.0 {
            penalty_for_very_short_bonds = 0.0;
        }
        if len > SKETCHER_EPSILON {
            l = l / len;
        }
        l = l * (self.k * dr + penalty_for_very_short_bonds * 10.0);
        m.atoms[self.atom1].force += l;
        m.atoms[self.atom2].force -= l;
    }
}

impl Bend {
    pub fn new(atom1: AtomId, atom2: AtomId, atom3: AtomId) -> Self {
        Bend {
            atom1,
            atom2,
            atom3,
            k: 1.0,
            k2: 0.05,
            rest_v: 120.0,
            is_ring: false,
        }
    }

    fn energy(&self, m: &Model, e: &mut f32) {
        let d_a = self.angle(m) - self.rest_v;
        *e += 0.5 * self.k * self.k2 * d_a * d_a * 10.0;
    }

    /// the angle in degrees
    fn angle(&self, m: &Model) -> f32 {
        let (p1, p2, p3) = (
            m.atoms[self.atom1].coordinates,
            m.atoms[self.atom2].coordinates,
            m.atoms[self.atom3].coordinates,
        );
        let (v1x, v1y) = (p1.x - p2.x, p1.y - p2.y);
        let (v2x, v2y) = (p3.x - p2.x, p3.y - p2.y);
        let mut d = (v1x * v1x + v1y * v1y).sqrt() * (v2x * v2x + v2y * v2y).sqrt();
        if d < SKETCHER_EPSILON {
            d = SKETCHER_EPSILON;
        }
        let mut cosine = (v1x * v2x + v1y * v2y) / d;
        if cosine < -1.0 {
            cosine = -1.0;
        } else if cosine > 1.0 {
            cosine = 1.0;
        }
        (f64::from(cosine.acos() * 180.0) / std::f64::consts::PI) as f32
    }

    fn score(&self, m: &mut Model, total: &mut f32) {
        let mut a = self.angle(m);
        if a < 0.0 {
            a = -a;
        }
        let mut target = self.rest_v;
        if target > 180.0 {
            target = 360.0 - target;
        }
        let d_a = target - a;
        self.energy(m, total);
        let (p1, p2, p3) = (
            m.atoms[self.atom1].coordinates,
            m.atoms[self.atom2].coordinates,
            m.atoms[self.atom3].coordinates,
        );
        let (v1x, v1y) = (p1.x - p2.x, p1.y - p2.y);
        let (v2x, v2y) = (p3.x - p2.x, p3.y - p2.y);
        let (v3x, v3y) = (p3.x - p1.x, p3.y - p1.y);
        let newk2 = self.k2;
        let mut n1 = PointF::new(v1y, -v1x);
        let mut n2 = PointF::new(v2y, -v2x);
        if n1.x * v3x + n1.y * v3y > 0.0 {
            n1 = n1 * -1.0;
        }
        if n2.x * v3x + n2.y * v3y < 0.0 {
            n2 = n2 * -1.0;
        }
        let mut q1 = (n1.x * n1.x + n1.y * n1.y).sqrt();
        if q1 < SKETCHER_EPSILON {
            q1 = SKETCHER_EPSILON;
        }
        let mut q2 = (n2.x * n2.x + n2.y * n2.y).sqrt();
        if q2 < SKETCHER_EPSILON {
            q2 = SKETCHER_EPSILON;
        }
        n1 = n1 / q1;
        n2 = n2 / q2;
        n1 = n1 * (self.k * newk2 * d_a);
        n2 = n2 * (self.k * newk2 * d_a);
        m.atoms[self.atom1].force += n1;
        m.atoms[self.atom3].force += n2;
        m.atoms[self.atom2].force -= n1 + n2;
    }
}

impl Clash {
    pub fn new(atom1: AtomId, atom2: AtomId, atom3: AtomId) -> Self {
        Clash {
            atom1,
            atom2,
            atom3,
            k: 1.0,
            k2: 0.1,
            rest_v: 900.0,
            squared_distance: 0.0,
        }
    }

    fn energy(&mut self, m: &Model, e: &mut f32) {
        self.squared_distance = squared_distance_point_segment(
            m.atoms[self.atom2].coordinates,
            m.atoms[self.atom1].coordinates,
            m.atoms[self.atom3].coordinates,
        )
        .0;
        if self.squared_distance > self.rest_v {
            return;
        }
        let dr = self.rest_v - self.squared_distance;
        if dr > 0.0 {
            *e += 0.5 * self.k * self.k2 * dr;
        }
    }

    pub fn score(&mut self, m: &mut Model, total: &mut f32, skip_force: bool) {
        self.energy(m, total);
        if skip_force || self.squared_distance > self.rest_v {
            return;
        }
        let atom_p = m.atoms[self.atom2].coordinates;
        let projection = project_point_on_line(
            atom_p,
            m.atoms[self.atom1].coordinates,
            m.atoms[self.atom3].coordinates,
        );
        let mut f = atom_p - projection;
        f.normalize();
        f = f * ((self.rest_v - self.squared_distance) * self.k * self.k2);
        m.atoms[self.atom2].force += f;
        m.atoms[self.atom1].force -= f * 0.5;
        m.atoms[self.atom3].force -= f * 0.5;
    }
}

impl EzConstrain {
    fn inversion(&self, m: &Model) -> bool {
        same_side(
            m.atoms[self.atom1].coordinates,
            m.atoms[self.atom4].coordinates,
            m.atoms[self.atom2].coordinates,
            m.atoms[self.atom3].coordinates,
        ) != self.is_z
    }

    fn score(&self, m: &mut Model, total: &mut f32) {
        if !self.inversion(m) {
            return;
        }
        *total += 5000.0;
        let c = |a: AtomId| m.atoms[a].coordinates;
        let projection1 = project_point_on_line(c(self.atom1), c(self.atom2), c(self.atom3));
        let projection2 = project_point_on_line(c(self.atom4), c(self.atom2), c(self.atom3));
        let mut side_atom = self.atom1;
        let mut double_bond_atom = self.atom2;
        let mut projection = projection1;
        if squared_distance(c(self.atom1), projection1)
            > squared_distance(c(self.atom4), projection2)
        {
            side_atom = self.atom4;
            double_bond_atom = self.atom3;
            projection = projection2;
        }
        let mut force = projection - m.atoms[side_atom].coordinates;
        if self.force_movement {
            m.atoms[side_atom].coordinates += force;
            m.atoms[double_bond_atom].coordinates -= force;
            m.atoms[side_atom].force = PointF::new(0.0, 0.0);
            m.atoms[double_bond_atom].force = PointF::new(0.0, 0.0);
        } else {
            force.normalize();
            force = force * 10.0;
            m.atoms[side_atom].force += force;
            m.atoms[double_bond_atom].force -= force;
        }
    }
}

impl Interaction {
    /// `score`: adds the energy and applies the forces
    pub fn score(&mut self, m: &mut Model, total: &mut f32, skip_force: bool) {
        match self {
            Interaction::Stretch(i) => i.score(m, total),
            Interaction::Bend(i) => i.score(m, total),
            Interaction::Clash(i) => i.score(m, total, skip_force),
            Interaction::Constraint(i) => {
                *total += i.k * squared_distance(m.atoms[i.atom1].coordinates, i.origin);
            }
            Interaction::EzConstrain(i) => i.score(m, total),
        }
    }
}

/// `CoordgenMinimizer`
#[derive(Debug, Clone)]
pub struct Minimizer {
    pub atoms: Vec<AtomId>,
    pub molecules: Vec<MolId>,
    pub max_iterations: usize,
    pub interactions: Vec<Interaction>,
    /// indices in `interactions` (`_intramolecularClashInteractions`)
    pub intramolecular_clash_interactions: Vec<usize>,
    pub even_angles: bool,
}

impl Default for Minimizer {
    fn default() -> Self {
        Minimizer {
            atoms: Vec::new(),
            molecules: Vec::new(),
            max_iterations: 1000,
            interactions: Vec::new(),
            intramolecular_clash_interactions: Vec::new(),
            even_angles: false,
        }
    }
}

impl Minimizer {
    pub fn clear_interactions(&mut self) {
        self.interactions.clear();
        self.intramolecular_clash_interactions.clear();
    }

    /// `scoreInteractions`
    pub fn score_interactions(&mut self, m: &mut Model) -> f32 {
        let mut total = 0.0f32;
        for i in &mut self.interactions {
            i.score(m, &mut total, false);
        }
        total
    }

    /// `applyForces`: false when the atoms hardly moved
    pub fn apply_forces(&self, m: &mut Model, maxd: f32) -> bool {
        let delta = 0.001f32;
        let mut distance = 0.0f32;
        for &a in &self.atoms {
            if m.atoms[a].fixed {
                continue;
            }
            let mut displacement = m.atoms[a].force * FORCE_MULTIPLIER;
            if displacement.x.is_nan() || displacement.y.is_nan() {
                displacement = PointF::new(0.0, 0.0);
            }
            let mut dsquare = displacement.x * displacement.x + displacement.y * displacement.y;
            if dsquare < SKETCHER_EPSILON {
                dsquare = SKETCHER_EPSILON;
            }
            if dsquare > maxd * maxd {
                displacement = displacement * (maxd / dsquare.sqrt());
            }
            m.atoms[a].coordinates += displacement;
            distance += displacement.square_length();
            m.atoms[a].force = PointF::new(0.0, 0.0);
        }
        distance >= delta
    }

    /// `run`. The coordinates of the lowest energy are never restored: upstream never updates its running minimum.
    pub fn run(&mut self, m: &mut Model) {
        if self.interactions.is_empty() {
            self.setup_interactions(m, false);
        }
        let mut energies = vec![0.0f32; self.max_iterations];
        for iterations in 0..self.max_iterations {
            energies[iterations] = self.score_interactions(m);
            if !self.apply_forces(m, 0.1) {
                break;
            }
            if iterations < 2 * ITERATION_HISTORY_SIZE {
                continue;
            }
            if energies[iterations - ITERATION_HISTORY_SIZE] - energies[iterations]
                < MAX_NET_ENERGY_CHANGE
            {
                break;
            }
        }
    }

    /// `setupInteractions`
    pub fn setup_interactions(&mut self, m: &Model, intrafragment_clashes: bool) {
        self.clear_interactions();
        for mol in self.molecules.clone() {
            self.add_interactions_of_molecule(m, mol, intrafragment_clashes);
        }
    }

    /// `addInteractionsOfMolecule`
    pub fn add_interactions_of_molecule(
        &mut self,
        m: &Model,
        mol: MolId,
        intrafragment_clashes: bool,
    ) {
        self.add_clash_interactions_of_molecule(m, mol, intrafragment_clashes);
        self.add_stretch_interactions_of_molecule(m, mol);
        self.add_bend_interactions_of_molecule(m, mol);
        self.add_chiral_inversion_constraints_of_molecule(m, mol);
    }

    fn fragment_has_no_dofs_at(m: &Model, a: AtomId) -> bool {
        let f = m.atoms[a].fragment.expect("atom in a fragment");
        m.fragment_dofs_of_atom(f, a).is_empty()
    }

    fn add_clash_interactions_of_molecule(
        &mut self,
        m: &Model,
        mol: MolId,
        intrafragment_clashes: bool,
    ) {
        let atoms = &m.molecules[mol].atoms;
        let bonds = &m.molecules[mol].bonds;
        if atoms.len() <= 1 {
            return;
        }
        for &atom in atoms {
            for &bond in bonds {
                let at2 = atom;
                let at1 = m.bonds[bond].start_atom;
                let at3 = m.bonds[bond].end_atom;
                if at1 == at2 || at1 == at3 || at2 == at3 {
                    continue;
                }
                if Self::fragment_has_no_dofs_at(m, at1)
                    && Self::fragment_has_no_dofs_at(m, at2)
                    && Self::fragment_has_no_dofs_at(m, at3)
                    && !intrafragment_clashes
                {
                    if m.atoms[at1].fragment == m.atoms[at2].fragment {
                        continue;
                    }
                    if m.atoms[at3].fragment == m.atoms[at2].fragment {
                        continue;
                    }
                }
                if m.atoms[at2].fixed && m.atoms[at1].fixed && m.atoms[at3].fixed {
                    continue;
                }
                if m.atoms[at1].neighbors.contains(&at2) {
                    continue;
                }
                if m.atoms[at3].neighbors.contains(&at2) {
                    continue;
                }
                if !(m.atoms[at1].rigid && m.atoms[at2].rigid && m.atoms[at3].rigid) {
                    let mut interaction = Clash::new(at1, at2, at3);
                    let mut rest_vk = 0.8f32;
                    if m.atoms[at2].atomic_number == 6 && m.atoms[at2].charge == 0 {
                        rest_vk -= 0.1;
                    }
                    if m.atoms[at1].atomic_number == 6
                        && m.atoms[at1].charge == 0
                        && m.atoms[at3].atomic_number == 6
                        && m.atoms[at3].charge == 0
                    {
                        rest_vk -= 0.1;
                    }
                    interaction.rest_v = (BONDLENGTH * rest_vk) * (BONDLENGTH * rest_vk);
                    self.intramolecular_clash_interactions
                        .push(self.interactions.len());
                    self.interactions.push(Interaction::Clash(interaction));
                }
            }
        }
    }

    fn add_stretch_interactions_of_molecule(&mut self, m: &Model, mol: MolId) {
        for &bo in &m.molecules[mol].bonds {
            let (at1, at2) = (m.bonds[bo].start_atom, m.bonds[bo].end_atom);
            let mut interaction = Stretch::new(at1, at2);
            interaction.k *= 0.1;
            interaction.rest_v = BONDLENGTH;
            if m.atoms[at1].rigid && m.atoms[at2].rigid {
                let v = m.atoms[at2].coordinates - m.atoms[at1].coordinates;
                interaction.rest_v = v.length();
            }
            if let Some(r) = m.share_a_ring(at1, at2) {
                if !m.rings[r].is_macrocycle() {
                    interaction.k *= 50.0;
                }
            }
            self.interactions.push(Interaction::Stretch(interaction));
        }
    }

    fn add_bend_interactions_of_molecule(&mut self, m: &Model, mol: MolId) {
        for &at in &m.molecules[mol].atoms {
            let mut interactions: Vec<Bend> = Vec::new();
            let mut ring_interactions: Vec<usize> = Vec::new();
            let mut non_ring_interactions: Vec<usize> = Vec::new();
            let nbonds = m.atoms[at].neighbors.len();
            let mut inverted_macrocycle_bond = false;
            if nbonds > 1 {
                let ordered_neighs = m.clockwise_ordered_neighbors(at);
                let mut angle = 120.0f32;
                if nbonds == 2
                    && m.bonds[m.atoms[at].bonds[0]].bond_order
                        + m.bonds[m.atoms[at].bonds[1]].bond_order
                        > 3
                {
                    angle = 180.0;
                }
                if nbonds > 2 {
                    angle = 360.0 / nbonds as f32;
                }
                for i in 0..nbonds {
                    let j = (i + nbonds - 1) % nbonds;
                    if nbonds == 2 && i == 1 {
                        continue;
                    }
                    let (at1, at2, at3) = (ordered_neighs[i], at, ordered_neighs[j]);
                    let mut interaction = Bend::new(at1, at2, at3);
                    interaction.rest_v = angle;
                    let idx = interactions.len();
                    if let Some(r) = m.same_ring3(at, ordered_neighs[i], ordered_neighs[j]) {
                        let ring = &m.rings[r];
                        if !ring.is_macrocycle() {
                            let mut extra_atoms: i32 = 0;
                            for k in 0..ring.fused_with.len() {
                                if m.rings[ring.fused_with[k]].is_macrocycle() {
                                    continue;
                                }
                                if ring.fusion_atoms[k].len() > 2 {
                                    extra_atoms += (m.rings[ring.fused_with[k]].atoms.len() as i64
                                        - ring.fusion_atoms[k].len() as i64)
                                        as i32;
                                }
                            }
                            interaction.is_ring = true;
                            interaction.k *= 100.0;
                            interaction.rest_v = (180.0
                                - (360.0 / (ring.size() as i32 + extra_atoms) as f64))
                                as f32;
                            ring_interactions.push(idx);
                        } else {
                            if nbonds == 3 {
                                let other = ordered_neighs
                                    .iter()
                                    .copied()
                                    .find(|&a| a != at1 && a != at3);
                                if let Some(other) = other {
                                    if same_side(
                                        m.atoms[at3].coordinates,
                                        m.atoms[other].coordinates,
                                        m.atoms[at1].coordinates,
                                        m.atoms[at2].coordinates,
                                    ) {
                                        inverted_macrocycle_bond = true;
                                    }
                                }
                            }
                            let mut fused_to_ring = ordered_neighs.len() > 2;
                            for &a in &ordered_neighs {
                                if m.share_a_ring(at, a).is_none() {
                                    fused_to_ring = false;
                                    break;
                                }
                            }
                            if fused_to_ring || inverted_macrocycle_bond {
                                ring_interactions.push(idx);
                            } else {
                                non_ring_interactions.push(idx);
                            }
                        }
                    } else {
                        non_ring_interactions.push(idx);
                    }
                    if m.atoms[interaction.atom1].rigid
                        && m.atoms[interaction.atom2].rigid
                        && m.atoms[interaction.atom3].rigid
                    {
                        interaction.rest_v = unsigned_angle(
                            m.atoms[interaction.atom1].coordinates,
                            m.atoms[interaction.atom2].coordinates,
                            m.atoms[interaction.atom3].coordinates,
                        );
                    }
                    interactions.push(interaction);
                }
            }
            if ring_interactions.len() != 1 || non_ring_interactions.len() != 2 {
                inverted_macrocycle_bond = false;
            }
            if !ring_interactions.is_empty() {
                let mut total_angle_in_rings = 0.0f32;
                for &i in &ring_interactions {
                    total_angle_in_rings += interactions[i].rest_v;
                }
                if inverted_macrocycle_bond {
                    total_angle_in_rings = 360.0 - total_angle_in_rings;
                }
                for &i in &non_ring_interactions {
                    interactions[i].rest_v =
                        (360.0 - total_angle_in_rings) / non_ring_interactions.len() as f32;
                }
            } else if non_ring_interactions.len() == 4 {
                if m.atoms[at].cross_layout || self.even_angles {
                    for &i in &non_ring_interactions {
                        interactions[i].rest_v = 90.0;
                    }
                } else {
                    let mut index_of_biggest_angle = 0;
                    let mut biggest_angle = 0.0f32;
                    for (counter, &i) in non_ring_interactions.iter().enumerate() {
                        let b = &interactions[i];
                        let angle = unsigned_angle(
                            m.atoms[b.atom1].coordinates,
                            m.atoms[b.atom2].coordinates,
                            m.atoms[b.atom3].coordinates,
                        );
                        if angle > biggest_angle {
                            biggest_angle = angle;
                            index_of_biggest_angle = counter;
                        }
                    }
                    interactions[non_ring_interactions[index_of_biggest_angle]].rest_v = 120.0;
                    interactions[non_ring_interactions[(index_of_biggest_angle + 1) % 4]].rest_v =
                        90.0;
                    interactions[non_ring_interactions[(index_of_biggest_angle + 2) % 4]].rest_v =
                        60.0;
                    interactions[non_ring_interactions[(index_of_biggest_angle + 3) % 4]].rest_v =
                        90.0;
                }
            } else if non_ring_interactions.len() > 4 {
                for &i in &non_ring_interactions {
                    interactions[i].rest_v = (360.0f64 / non_ring_interactions.len() as f64) as f32;
                }
            }
            for i in interactions {
                if !(m.atoms[i.atom1].fixed && m.atoms[i.atom2].fixed && m.atoms[i.atom3].fixed) {
                    self.interactions.push(Interaction::Bend(i));
                }
            }
        }
    }

    fn add_chiral_inversion_constraints_of_molecule(&mut self, m: &Model, mol: MolId) {
        for &r in &m.molecules[mol].rings {
            if !m.rings[r].is_macrocycle() {
                continue;
            }
            let ring_atoms = m.rings[r].atoms.clone();
            let atoms = m.order_chain_of_atoms(&ring_atoms, ring_atoms[0]);
            let size = atoms.len();
            for i in 0..size {
                let a1 = (i + size - 1) % size;
                let a11 = (i + size - 2) % size;
                let a2 = (i + 1) % size;
                let bond = m.bond_to(atoms[a1], atoms[i]).expect("ring bond");
                if m.bond_is_stereo(bond) {
                    let cis = m.bond_marked_as_cis(bond, atoms[a11], atoms[a2]);
                    self.interactions
                        .push(Interaction::EzConstrain(EzConstrain {
                            atom1: atoms[a11],
                            atom2: atoms[a1],
                            atom3: atoms[i],
                            atom4: atoms[a2],
                            is_z: cis,
                            force_movement: false,
                        }));
                }
            }
        }
    }

    /// `addConstrainedInteractionsOfMolecule`
    pub fn add_constrained_interactions_of_molecule(&mut self, m: &Model, mol: MolId) {
        for &a in &m.molecules[mol].atoms {
            if m.atoms[a].constrained {
                self.intramolecular_clash_interactions
                    .push(self.interactions.len());
                self.interactions.push(Interaction::Constraint(Constraint {
                    atom1: a,
                    k: CONSTRAINT_SCALE,
                    origin: m.atoms[a].template_coordinates,
                }));
            }
        }
    }

    /// `minimizeMolecule`: minimizes and returns to the previous coordinates if a stereo bond got inverted
    pub fn minimize_molecule(&mut self, m: &mut Model, mol: MolId) {
        let atoms = m.molecules[mol].atoms.clone();
        let previous: Vec<PointF> = atoms.iter().map(|&a| m.atoms[a].coordinates).collect();
        self.clear_interactions();
        self.add_interactions_of_molecule(m, mol, true);
        self.run(m);
        for b in m.molecules[mol].bonds.clone() {
            if !m.bond_check_stereo_chemistry(b) {
                for (i, &a) in atoms.iter().enumerate() {
                    m.set_coordinates(a, previous[i]);
                }
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::testutil::build;
    use super::*;
    use std::fmt::Write;

    fn run(input: &str) -> String {
        let (mut model, mol) = build(input);
        model.force_update_struct(mol);
        model.split_into_fragments(mol);
        let mut min = Minimizer {
            atoms: model.molecules[mol].atoms.clone(),
            molecules: vec![mol],
            ..Minimizer::default()
        };
        min.minimize_molecule(&mut model, mol);
        let e = min.score_interactions(&mut model);
        let count = |f: fn(&Interaction) -> bool| min.interactions.iter().filter(|i| f(i)).count();
        let mut out = String::new();
        writeln!(
            out,
            "interactions {} clash {} stretch {} bend {} energy {:08x}",
            min.interactions.len(),
            min.intramolecular_clash_interactions.len(),
            count(|i| matches!(i, Interaction::Stretch(_))),
            count(|i| matches!(i, Interaction::Bend(_))),
            e.to_bits()
        )
        .unwrap();
        for (i, &a) in model.molecules[mol].atoms.iter().enumerate() {
            let c = model.atoms[a].coordinates;
            writeln!(out, "a {i} {:08x} {:08x}", c.x.to_bits(), c.y.to_bits()).unwrap();
        }
        out
    }

    /// Equal up to the platform's `acos`. The minimization stops on thresholds, so other libms can end in a slightly
    /// different state: off macOS (where the corpus was made) coordinates must agree within 2 units (a bond is 50) and
    /// the energy is not compared.
    fn same(got: &str, want: &str) -> bool {
        let strict = cfg!(target_os = "macos");
        let (gl, wl): (Vec<&str>, Vec<&str>) = (got.lines().collect(), want.lines().collect());
        gl.len() == wl.len()
            && gl.iter().zip(&wl).all(|(g, w)| {
                let (gt, wt): (Vec<&str>, Vec<&str>) = (
                    g.split_whitespace().collect(),
                    w.split_whitespace().collect(),
                );
                let energy_line = w.starts_with("interactions");
                gt.len() == wt.len()
                    && gt.iter().zip(&wt).all(|(a, b)| {
                        let hex =
                            |t: &str| t.len() == 8 && t.chars().all(|c| c.is_ascii_hexdigit());
                        if a == b {
                            return true;
                        }
                        if hex(a) && hex(b) {
                            if energy_line && !strict {
                                return true;
                            }
                            let f = |t: &str| f32::from_bits(u32::from_str_radix(t, 16).unwrap());
                            let (x, y) = (f(a), f(b));
                            return if strict {
                                (x - y).abs() <= 1e-3 * (1.0 + x.abs().max(y.abs()))
                            } else {
                                (x - y).abs() <= 2.0
                            };
                        }
                        false
                    })
            })
    }

    /// the minimized coordinates must be the ones CoordgenMinimizer makes (`tests/data/coordgen_minim.txt`)
    #[test]
    fn minimizer_matches_cpp() {
        let text = include_str!("../../../../../tests/data/coordgen_minim.txt");
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
