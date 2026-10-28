//! The degrees of freedom of a fragment (`CoordgenFragmentDOF` and subclasses: states, `apply`, penalties),
//! `sketcherMinimizerFragment::setCoordinates`, the placement of the fragments of a molecule
//! (`buildMoleculeFromFragments`) and the choice of the side of a fragment (`alignWithParentDirection`).

use super::fragment::{Dof, DofKind};
use super::maths::{
    BONDLENGTH, PointF, SKETCHER_EPSILON, dot_product, mirror_point, round_to_two_decimal_digits,
};
use super::model::{AtomId, BondId, FragmentId, Model, MolId};

const INVERTED_MACROCYCLE_BOND_PENALTY: f32 = 100.0;
const SCALE_ATOMS_PENALTY: f32 = 50.0;
const ROTATE_FRAGMENT_PENALTY: f32 = 400.0;
const BREAK_CHAIN_PENALTY: f32 = 10.0;
const CHANGE_PARENT_BOND_PENALTY: f32 = 200.0;
const FLIP_CONSTRAINED_FRAGMENT_PENALTY: f32 = 1000.0;

const SCORE_MULTIPLIER_FOR_DOUBLE_BONDS: f32 = 0.82;
const SCORE_MULTIPLIER_FOR_SINGLE_BONDED_HETEROATOMS: f32 = 0.9;
const SCORE_MULTIPLIER_FOR_FRAGMENTS: f32 = 0.1;

impl Model {
    fn dof(&self, d: usize) -> &Dof {
        &self.dofs[d]
    }

    /// `numberOfStates`
    pub fn dof_number_of_states(&self, d: usize) -> u16 {
        let dof = self.dof(d);
        match dof.kind {
            DofKind::FlipFragment => {
                if self.fragments[dof.fragment].parent.is_none() {
                    1
                } else {
                    2
                }
            }
            DofKind::ChangeParentBondLength => 7,
            DofKind::RotateFragment => {
                if self.fragments[dof.fragment].parent.is_none() {
                    1
                } else {
                    5
                }
            }
            DofKind::ScaleAtoms | DofKind::InvertBond => 2,
        }
    }

    /// `changeState`
    pub fn dof_change_state(&mut self, d: usize) {
        let n = self.dof_number_of_states(d);
        let dof = &mut self.dofs[d];
        dof.current_state = (dof.current_state + 1) % n;
    }

    /// `getCurrentPenalty`
    pub fn dof_current_penalty(&self, d: usize) -> f32 {
        let dof = self.dof(d);
        let state = dof.current_state;
        let half = f32::from(state.div_ceil(2));
        match dof.kind {
            DofKind::FlipFragment => {
                let frag = &self.fragments[dof.fragment];
                let mut penalty = 0.0f32;
                if state != 0 && frag.constrained_flip {
                    penalty += FLIP_CONSTRAINED_FRAGMENT_PENALTY;
                }
                if frag.is_chain && frag.parent.is_some_and(|p| self.fragments[p].is_chain) {
                    penalty += BREAK_CHAIN_PENALTY;
                }
                penalty
            }
            DofKind::ChangeParentBondLength if state != 0 => CHANGE_PARENT_BOND_PENALTY * half,
            DofKind::RotateFragment if state != 0 => ROTATE_FRAGMENT_PENALTY * half,
            DofKind::ScaleAtoms if state != 0 => SCALE_ATOMS_PENALTY * dof.atoms.len() as f32,
            DofKind::InvertBond if state != 0 => INVERTED_MACROCYCLE_BOND_PENALTY,
            _ => 0.0,
        }
    }

    /// `apply`
    pub fn dof_apply(&mut self, d: usize) {
        let dof = self.dofs[d].clone();
        let state = dof.current_state;
        if state == 0 {
            return;
        }
        let half = i32::from(state.div_ceil(2));
        match dof.kind {
            DofKind::FlipFragment => {
                for &a in self.fragments[dof.fragment].coordinates.keys() {
                    self.atoms[a].coordinates.y = -self.atoms[a].coordinates.y;
                }
            }
            DofKind::ChangeParentBondLength => {
                let mut scale = 1.6f64.powf(f64::from(half)) as f32;
                if state % 2 == 0 {
                    scale = 1.0 / scale;
                }
                let move_by = BONDLENGTH * (scale - 1.0);
                for &a in self.fragments[dof.fragment].coordinates.keys() {
                    self.atoms[a].coordinates.x += move_by;
                }
            }
            DofKind::RotateFragment => {
                let mut angle = (std::f64::consts::PI / 180.0 * 15.0 * f64::from(half)) as f32;
                if state % 2 == 0 {
                    angle = -angle;
                }
                let (sine, cosine) = (angle.sin(), angle.cos());
                let origin = PointF::new(-BONDLENGTH, 0.0);
                let keys: Vec<AtomId> = self.fragments[dof.fragment]
                    .coordinates
                    .keys()
                    .copied()
                    .collect();
                for a in keys {
                    let mut coords = self.atoms[a].coordinates - origin;
                    coords.rotate(sine, cosine);
                    self.set_coordinates(a, coords + origin);
                }
            }
            DofKind::ScaleAtoms => {
                let pivot = self.atoms[dof.pivot.expect("pivot")].coordinates;
                let scale = 0.4f32;
                for &a in &dof.atoms {
                    let distance = self.atoms[a].coordinates - pivot;
                    self.set_coordinates(a, distance * scale + pivot);
                }
            }
            DofKind::InvertBond => {
                let pivot_atom = dof.pivot.expect("pivot");
                let bound = dof.atoms[0];
                let pivot = self.atoms[pivot_atom].coordinates;
                let bond_dir = self.atoms[bound].coordinates - pivot;
                let normal = PointF::new(bond_dir.y, -bond_dir.x);
                let point1 = pivot + normal;
                let point2 = pivot - normal;
                for &a in &dof.atoms {
                    let m = mirror_point(self.atoms[a].coordinates, point1, point2);
                    self.set_coordinates(a, m);
                }
            }
        }
    }

    /// `sketcherMinimizerFragment::setCoordinates`: places the fragment at `position` rotated by `angle`
    pub fn fragment_set_coordinates(&mut self, f: FragmentId, position: PointF, angle: f32) {
        let (sine, cosine) = (angle.sin(), angle.cos());
        let coords: Vec<(AtomId, PointF)> = self.fragments[f]
            .coordinates
            .iter()
            .map(|(&a, &p)| (a, p))
            .collect();
        debug_assert_eq!(
            coords.len(),
            self.fragments[f].atoms().len() + self.fragments[f].children.len()
        );
        for &(a, p) in &coords {
            self.set_coordinates(a, p);
        }
        for d in self.fragments[f].dofs.clone() {
            self.dof_apply(d);
        }
        for &(a, _) in &coords {
            let mut c = self.atoms[a].coordinates;
            c.rotate(sine, cosine);
            self.set_coordinates(a, c + position);
        }
    }

    /// `getAllTerminalBonds`
    fn all_terminal_bonds(&self, f: FragmentId) -> Vec<BondId> {
        let frag = &self.fragments[f];
        let mut bonds = Vec::new();
        for &b in frag.bonds() {
            let bd = &self.bonds[b];
            if self.atoms[bd.start_atom].neighbors.len() == 1
                || self.atoms[bd.end_atom].neighbors.len() == 1
            {
                bonds.push(b);
            }
        }
        for &child in &frag.children {
            bonds.push(self.fragments[child].bond_to_parent.expect("child bond"));
        }
        if frag.parent.is_some() {
            bonds.push(frag.bond_to_parent.expect("parent bond"));
        }
        bonds
    }

    fn single_bonded_heteroatom(&self, b: BondId) -> bool {
        let bd = &self.bonds[b];
        (self.atoms[bd.start_atom].neighbors.len() == 1
            && self.atoms[bd.start_atom].atomic_number != 6)
            || (self.atoms[bd.end_atom].neighbors.len() == 1
                && self.atoms[bd.end_atom].atomic_number != 6)
    }

    /// `findDirectionsToAlignWith`
    fn find_directions_to_align_with(&self, f: FragmentId) -> Vec<(PointF, f32)> {
        let mut chain_dirs = Vec::new();
        let frag = &self.fragments[f];
        let btp = &self.bonds[frag.bond_to_parent.expect("parent bond")];
        let origin =
            (self.atoms[btp.start_atom].coordinates + self.atoms[btp.end_atom].coordinates) * 0.5;
        let parent = frag.parent.expect("parent");
        for b in self.all_terminal_bonds(parent) {
            let bond = &self.bonds[b];
            if self.atoms[bond.end_atom].fragment == Some(f) {
                continue;
            }
            let mut direction = origin
                - (self.atoms[bond.start_atom].coordinates + self.atoms[bond.end_atom].coordinates)
                    * 0.5;
            direction.normalize();
            let mut score = 1.0f32;
            if bond.bond_order == 2 {
                score *= SCORE_MULTIPLIER_FOR_DOUBLE_BONDS;
            }
            if self.single_bonded_heteroatom(b) {
                score *= SCORE_MULTIPLIER_FOR_SINGLE_BONDED_HETEROATOMS;
            }
            if self.atoms[bond.end_atom].fragment != Some(parent)
                || self.atoms[bond.start_atom].fragment != Some(parent)
            {
                let end_frag = self.atoms[bond.end_atom].fragment.expect("fragment");
                score = self.fragments[end_frag].longest_chain_from_here
                    * SCORE_MULTIPLIER_FOR_FRAGMENTS;
                if let Some(pp) = self.fragments[parent].parent {
                    if self.atoms[bond.start_atom].fragment == Some(pp) {
                        score *= 100.0;
                    }
                }
            }
            chain_dirs.push((direction, score));
        }
        chain_dirs
    }

    /// `testAlignment`
    fn test_alignment(direction: PointF, templat: (PointF, f32)) -> f32 {
        let mut dot = dot_product(direction, templat.0);
        if dot < 0.0 {
            dot = 0.0;
        }
        let mut score = dot * dot;
        if dot > 1.0 - SKETCHER_EPSILON {
            score += 1000.0;
        }
        score * templat.1
    }

    /// `scoreDirections`: whether the fragment has to be inverted
    fn score_directions(&self, f: FragmentId, angle: f32, directions: &[(PointF, f32)]) -> bool {
        let (sine, cosine) = (angle.sin(), angle.cos());
        let mut best_score = 0.0f32;
        let mut invert = false;
        let frag = &self.fragments[f];
        for b in self.all_terminal_bonds(f) {
            let bond = &self.bonds[b];
            if self.atoms[bond.start_atom].fragment != Some(f) {
                continue;
            }
            let local = |a: AtomId| frag.coordinates.get(&a).copied().unwrap_or_default();
            let mut plain = (local(bond.start_atom) + local(bond.end_atom)) * 0.5
                - PointF::new(-BONDLENGTH * 0.5, 0.0);
            plain.normalize();
            let mut inverted = PointF::new(plain.x, -plain.y);
            plain.rotate(sine, cosine);
            inverted.rotate(sine, cosine);
            let mut modifier = 1.0f32;
            if bond.bond_order == 2 {
                modifier *= SCORE_MULTIPLIER_FOR_DOUBLE_BONDS;
            }
            if self.single_bonded_heteroatom(b) {
                modifier *= SCORE_MULTIPLIER_FOR_SINGLE_BONDED_HETEROATOMS;
            }
            if self.atoms[bond.end_atom].fragment != Some(f) {
                let end_frag = self.atoms[bond.end_atom].fragment.expect("fragment");
                modifier = self.fragments[end_frag].longest_chain_from_here
                    * SCORE_MULTIPLIER_FOR_FRAGMENTS;
            }
            for &direction in directions {
                let score_plain = Self::test_alignment(plain, direction) * modifier;
                if score_plain > best_score {
                    best_score = score_plain;
                    invert = false;
                }
                let score_inverted = Self::test_alignment(inverted, direction) * modifier;
                if score_inverted > best_score {
                    best_score = score_inverted;
                    invert = true;
                }
            }
        }
        invert
    }

    /// `RMSD`
    fn rmsd(templates: &[PointF], points: &[PointF]) -> f32 {
        let counter = templates.len();
        let mut total = 0.0f32;
        for i in 0..templates.len() {
            let diff = templates[i] - points[i];
            total += diff.x * diff.x + diff.y * diff.y;
        }
        if counter > 0 {
            total /= counter as f32;
        }
        total.sqrt()
    }

    /// `alignWithParentDirectionConstrained`
    fn align_with_parent_direction_constrained(
        &self,
        f: FragmentId,
        position: PointF,
        angle: f32,
    ) -> bool {
        let (sine, cosine) = (angle.sin(), angle.cos());
        let (mut templates, mut plain, mut flipped) = (Vec::new(), Vec::new(), Vec::new());
        // upstream iterates a map ordered by pointer; atom index order is used
        for (&a, &p) in &self.fragments[f].coordinates {
            if self.atoms[a].constrained {
                let mut plain_atom = p;
                let mut flipped_atom = PointF::new(p.x, -p.y);
                plain_atom.rotate(sine, cosine);
                flipped_atom.rotate(sine, cosine);
                templates.push(self.atoms[a].template_coordinates);
                plain.push(plain_atom + position);
                flipped.push(flipped_atom + position);
            }
        }
        let score_plain = round_to_two_decimal_digits(Self::rmsd(&templates, &plain));
        let score_flipped = round_to_two_decimal_digits(Self::rmsd(&templates, &flipped));
        score_flipped < score_plain
    }

    /// `alignWithParentDirection`: decides on which side of its parent bond the fragment is drawn
    pub fn align_with_parent_direction(&mut self, f: FragmentId, position: PointF, angle: f32) {
        if self.fragments[f].fixed {
            return;
        }
        let invert = if self.fragments[f].constrained {
            self.align_with_parent_direction_constrained(f, position, angle)
        } else {
            let directions = self.find_directions_to_align_with(f);
            self.score_directions(f, angle, &directions)
        };
        if invert {
            for p in self.fragments[f].coordinates.values_mut() {
                p.y = -p.y;
            }
            for a in self.fragments[f].atoms().to_vec() {
                if self.atoms[a].has_stereochemistry_set {
                    for b in self.atoms[a].bonds.clone() {
                        self.bonds[b].is_wedge = !self.bonds[b].is_wedge;
                    }
                }
            }
        }
    }

    /// `CoordgenMinimizer::buildMoleculeFromFragments`
    pub fn build_molecule_from_fragments(&mut self, mol: MolId, first_time: bool) {
        for f in self.molecules[mol].fragments.clone() {
            let mut angle = 0.0f32;
            let mut position = PointF::new(0.0, 0.0);
            if self.fragments[f].parent.is_some() {
                let b = self.fragments[f].bond_to_parent.expect("parent bond");
                let p1 = self.atoms[self.bonds[b].start_atom].coordinates;
                let p2 = self.atoms[self.bonds[b].end_atom].coordinates;
                let p = p2 - p1;
                angle = (-p.y).atan2(p.x);
                position = p2;
                if first_time {
                    self.align_with_parent_direction(f, position, angle);
                }
            }
            self.fragment_set_coordinates(f, position, angle);
        }
    }
}
