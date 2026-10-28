//! `sketcherMinimizerFragment`: a rigid fragment of a molecule (ring systems and the atoms bound to them by
//! non-rotatable bonds). The degrees of freedom (`CoordgenFragmentDOF`) follow with the minimizer.

use std::collections::BTreeMap;

use super::maths::PointF;
use super::model::{AtomId, BondId, FragmentId, Model, RingId};

#[derive(Debug, Clone, Default)]
pub struct Fragment {
    atoms: Vec<AtomId>,
    bonds: Vec<BondId>,
    rings: Vec<RingId>,
    pub parent: Option<FragmentId>,
    pub inter_fragment_bonds: Vec<BondId>,
    pub children: Vec<FragmentId>,
    /// local coordinates of the atoms (and of the first atom of every child), keyed by atom
    pub coordinates: BTreeMap<AtomId, PointF>,
    pub bond_to_parent_coordinates_start: PointF,
    pub bond_to_parent_coordinates_end: PointF,
    pub fixed: bool,
    pub is_templated: bool,
    pub constrained: bool,
    pub constrained_flip: bool,
    pub is_chain: bool,
    pub bond_to_parent: Option<BondId>,
    pub longest_chain_from_here: f32,
    pub number_of_children_atoms: usize,
    pub number_of_children_atoms_rank: f32,
}

impl Fragment {
    pub fn atoms(&self) -> &[AtomId] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[BondId] {
        &self.bonds
    }

    pub fn rings(&self) -> &[RingId] {
        &self.rings
    }
}

impl Model {
    pub fn new_fragment(&mut self) -> FragmentId {
        self.fragments.push(Fragment::default());
        self.fragments.len() - 1
    }

    /// `addAtom`: the atom belongs to the fragment from now on
    pub fn fragment_add_atom(&mut self, f: FragmentId, atom: AtomId) {
        self.fragments[f].atoms.push(atom);
        self.atoms[atom].fragment = Some(f);
    }

    pub fn fragment_add_bond(&mut self, f: FragmentId, bond: BondId) {
        self.fragments[f].bonds.push(bond);
    }

    pub fn fragment_add_ring(&mut self, f: FragmentId, ring: RingId) {
        self.fragments[f].rings.push(ring);
    }

    pub fn fragment_add_inter_fragment_bond(&mut self, f: FragmentId, bond: BondId) {
        self.fragments[f].inter_fragment_bonds.push(bond);
    }

    /// `totalWeight`
    pub fn fragment_total_weight(&self, f: FragmentId) -> u32 {
        self.fragments[f]
            .atoms
            .iter()
            .map(|&a| self.atoms[a].atomic_number + self.atoms[a].implicit_hs)
            .sum::<i32>() as u32
    }

    pub fn fragment_count_double_bonds(&self, f: FragmentId) -> u32 {
        self.fragments[f]
            .bonds
            .iter()
            .filter(|&&b| self.bonds[b].bond_order == 2)
            .count() as u32
    }

    pub fn fragment_count_heavy_atoms(&self, f: FragmentId) -> u32 {
        self.fragments[f]
            .atoms
            .iter()
            .filter(|&&a| self.atoms[a].atomic_number != 6)
            .count() as u32
    }

    pub fn fragment_count_constrained_atoms(&self, f: FragmentId) -> u32 {
        self.fragments[f]
            .atoms
            .iter()
            .filter(|&&a| self.atoms[a].constrained)
            .count() as u32
    }

    pub fn fragment_count_fixed_atoms(&self, f: FragmentId) -> u32 {
        self.fragments[f]
            .atoms
            .iter()
            .filter(|&&a| self.atoms[a].fixed)
            .count() as u32
    }

    /// `setAllCoordinatesToTemplate`
    pub fn fragment_set_all_coordinates_to_template(&mut self, f: FragmentId) {
        for a in self.fragments[f].atoms.clone() {
            self.atom_set_coordinates_to_template(a);
        }
        if let Some(b) = self.fragments[f].bond_to_parent {
            let (s, e) = (self.bonds[b].start_atom, self.bonds[b].end_atom);
            self.atom_set_coordinates_to_template(s);
            self.atom_set_coordinates_to_template(e);
        }
        for child in self.fragments[f].children.clone() {
            let b = self.fragments[child]
                .bond_to_parent
                .expect("child with a parent bond");
            let (s, e) = (self.bonds[b].start_atom, self.bonds[b].end_atom);
            self.atom_set_coordinates_to_template(s);
            self.atom_set_coordinates_to_template(e);
        }
    }

    /// `setCoordinatesToTemplate` of an atom
    pub fn atom_set_coordinates_to_template(&mut self, a: AtomId) {
        let c = self.atoms[a].template_coordinates;
        self.set_coordinates(a, c);
    }
}
