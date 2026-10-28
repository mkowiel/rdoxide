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
    /// ids (in `Model::dofs`) of the degrees of freedom of this fragment (`m_dofs`)
    pub dofs: Vec<usize>,
    /// the degrees of freedom that move an atom of this fragment (`m_dofsForAtom`)
    pub atom_dofs: BTreeMap<AtomId, Vec<usize>>,
}

/// the kinds of `CoordgenFragmentDOF` the fragment builder creates
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DofKind {
    /// `CoordgenFlipFragmentDOF`
    FlipFragment,
    /// `CoordgenChangeParentBondLengthFragmentDOF`
    ChangeParentBondLength,
    /// `CoordgenRotateFragmentDOF`
    RotateFragment,
    /// `CoordgenScaleAtomsDOF`
    ScaleAtoms,
    /// `CoordgenInvertBondDOF`
    InvertBond,
}

impl DofKind {
    /// `tier()`
    pub fn tier(self) -> i32 {
        match self {
            DofKind::FlipFragment => 0,
            DofKind::InvertBond => 1,
            DofKind::ChangeParentBondLength => 2,
            DofKind::RotateFragment => 3,
            DofKind::ScaleAtoms => 4,
        }
    }
}

/// `CoordgenFragmentDOF`
#[derive(Debug, Clone)]
pub struct Dof {
    pub kind: DofKind,
    pub fragment: FragmentId,
    pub pivot: Option<AtomId>,
    pub atoms: Vec<AtomId>,
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
    /// `CoordgenFragmentDOF::addAtom`
    pub fn dof_add_atom(&mut self, dof: usize, atom: AtomId) {
        self.dofs[dof].atoms.push(atom);
        let f = self.atoms[atom].fragment.expect("atom in a fragment");
        self.fragments[f]
            .atom_dofs
            .entry(atom)
            .or_default()
            .push(dof);
    }

    /// `getDofsOfAtom`
    pub fn fragment_dofs_of_atom(&self, f: FragmentId, atom: AtomId) -> Vec<usize> {
        self.fragments[f]
            .atom_dofs
            .get(&atom)
            .cloned()
            .unwrap_or_default()
    }

    /// creates a DOF of the fragment of `pivot` (not yet added to the fragment)
    pub fn new_dof(&mut self, kind: DofKind, pivot: AtomId) -> usize {
        let fragment = self.atoms[pivot].fragment.expect("pivot in a fragment");
        self.dofs.push(Dof {
            kind,
            fragment,
            pivot: Some(pivot),
            atoms: Vec::new(),
        });
        self.dofs.len() - 1
    }

    pub fn new_fragment(&mut self) -> FragmentId {
        let f = self.fragments.len();
        self.fragments.push(Fragment::default());
        // the constructor adds the fragment-wide degrees of freedom
        for kind in [
            DofKind::FlipFragment,
            DofKind::ChangeParentBondLength,
            DofKind::RotateFragment,
        ] {
            self.dofs.push(Dof {
                kind,
                fragment: f,
                pivot: None,
                atoms: Vec::new(),
            });
            let id = self.dofs.len() - 1;
            self.fragments[f].dofs.push(id);
        }
        f
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
