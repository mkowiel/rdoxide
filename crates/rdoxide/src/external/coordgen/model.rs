//! The model classes of coordgenlibs: `sketcherMinimizerAtom`, `Bond`, `Ring` and `Molecule` (the parts that do not
//! depend on the CIP ranking or on residues). The C++ object graph with raw pointers becomes an arena: `Model` owns all
//! atoms, bonds, rings and molecules and the objects refer to each other by index.

use std::collections::{HashSet, VecDeque};

use super::maths::{MACROCYCLE, PointF, SKETCHER_EPSILON, mirror_point, signed_angle};

pub type AtomId = usize;
pub type BondId = usize;
pub type RingId = usize;
pub type MolId = usize;
/// index of a fragment in the fragment list of the sketcher
pub type FragmentId = usize;

/// `COORDINATES_LIMIT` and `INVALID_COORDINATES`
pub const COORDINATES_LIMIT: f32 = 10000000.0;
pub const INVALID_COORDINATES: f32 = COORDINATES_LIMIT + 1.0;

/// `sketcherMinimizerAtomChiralityInfo::sketcherMinimizerChirality`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Chirality {
    Clockwise,
    CounterClockwise,
    #[default]
    Unspecified,
}

/// `sketcherMinimizerAtomChiralityInfo`
#[derive(Debug, Clone, Copy, Default)]
pub struct AtomChiralityInfo {
    pub looking_from: Option<AtomId>,
    pub atom1: Option<AtomId>,
    pub atom2: Option<AtomId>,
    pub direction: Chirality,
}

/// `sketcherMinimizerAtom`
#[derive(Debug, Clone)]
pub struct Atom {
    pub cross_layout: bool,
    pub fixed: bool,
    pub constrained: bool,
    pub rigid: bool,
    pub is_shared_and_inner: bool,
    pub hidden: bool,
    pub atomic_number: i32,
    pub charge: i32,
    pub valence: i32,
    pub general_use_n: i32,
    pub general_use_n2: i32,
    pub chm_n: i32,
    pub general_use_visited: bool,
    pub general_use_visited2: bool,
    pub clockwise_invert: bool,
    pub ignore_ring_chirality: bool,
    pub rs_priorities: Vec<i32>,
    pub implicit_hs: i32,
    pub molecule: Option<MolId>,
    pub fragment: Option<FragmentId>,
    pub neighbors: Vec<AtomId>,
    pub bonds: Vec<BondId>,
    pub rings: Vec<RingId>,
    pub pseudo_z: f32,
    pub x3d: f32,
    pub y3d: f32,
    pub z3d: f32,
    pub is_clashing: bool,
    pub needs_check_for_clashes: bool,
    pub visited: bool,
    pub coordinates_set: bool,
    pub is_r: bool,
    pub has_stereochemistry_set: bool,
    pub is_stereogenic: bool,
    pub has_ring_chirality: bool,
    pub coordinates: PointF,
    pub template_coordinates: PointF,
    pub force: PointF,
    pub chirality_info: AtomChiralityInfo,
}

impl Default for Atom {
    fn default() -> Self {
        Atom {
            cross_layout: false,
            fixed: false,
            constrained: false,
            rigid: false,
            is_shared_and_inner: false,
            hidden: false,
            atomic_number: 6,
            charge: 0,
            valence: -10,
            general_use_n: -1,
            general_use_n2: -1,
            chm_n: -1,
            general_use_visited: false,
            general_use_visited2: false,
            clockwise_invert: false,
            ignore_ring_chirality: false,
            rs_priorities: Vec::new(),
            implicit_hs: -1,
            molecule: None,
            fragment: None,
            neighbors: Vec::new(),
            bonds: Vec::new(),
            rings: Vec::new(),
            pseudo_z: 0.0,
            x3d: INVALID_COORDINATES,
            y3d: INVALID_COORDINATES,
            z3d: INVALID_COORDINATES,
            is_clashing: false,
            needs_check_for_clashes: false,
            visited: false,
            coordinates_set: false,
            is_r: true,
            has_stereochemistry_set: false,
            is_stereogenic: false,
            has_ring_chirality: false,
            coordinates: PointF::default(),
            template_coordinates: PointF::default(),
            force: PointF::default(),
            chirality_info: AtomChiralityInfo::default(),
        }
    }
}

/// `sketcherMinimizerBondStereoInfo::sketcherMinimizerBondStereo`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BondStereo {
    Cis,
    Trans,
    #[default]
    Unspecified,
}

/// `sketcherMinimizerBondStereoInfo`
#[derive(Debug, Clone, Copy, Default)]
pub struct BondStereoInfo {
    pub atom1: Option<AtomId>,
    pub atom2: Option<AtomId>,
    pub stereo: BondStereo,
}

/// `sketcherMinimizerBond`
#[derive(Debug, Clone)]
pub struct Bond {
    pub start_atom: AtomId,
    pub end_atom: AtomId,
    pub bond_order: i32,
    pub skip: bool,
    pub is_ze_active: bool,
    pub is_z: bool,
    pub chm_n: i32,
    pub stereo: BondStereoInfo,
    pub is_wedge: bool,
    pub is_reversed: bool,
    pub has_stereochemistry_display: bool,
    pub sssr_visited: bool,
    pub sssr_parent_at_start: bool,
    pub ignore_ze: bool,
    pub sssr_parent: Option<BondId>,
    pub crossing_bond_penalty_multiplier: f32,
    pub rings: Vec<RingId>,
}

impl Bond {
    pub fn new(start_atom: AtomId, end_atom: AtomId) -> Self {
        Bond {
            start_atom,
            end_atom,
            bond_order: 1,
            skip: false,
            is_ze_active: false,
            is_z: false,
            chm_n: -1,
            stereo: BondStereoInfo::default(),
            is_wedge: false,
            is_reversed: false,
            has_stereochemistry_display: false,
            sssr_visited: false,
            sssr_parent_at_start: true,
            ignore_ze: false,
            sssr_parent: None,
            crossing_bond_penalty_multiplier: 1.0,
            rings: Vec::new(),
        }
    }
}

/// `sketcherMinimizerRing`
#[derive(Debug, Clone, Default)]
pub struct Ring {
    pub fused_with: Vec<RingId>,
    pub fusion_atoms: Vec<Vec<AtomId>>,
    pub fusion_bonds: Vec<BondId>,
    pub visited: bool,
    pub coordinates_generated: bool,
    /// not central
    pub side: bool,
    pub atoms: Vec<AtomId>,
    pub bonds: Vec<BondId>,
}

impl Ring {
    pub fn size(&self) -> usize {
        self.atoms.len()
    }

    pub fn is_macrocycle(&self) -> bool {
        self.size() >= MACROCYCLE
    }

    pub fn contains_atom(&self, a: AtomId) -> bool {
        self.atoms.contains(&a)
    }

    pub fn contains_bond(&self, b: BondId) -> bool {
        self.bonds.contains(&b)
    }

    pub fn same_as(&self, ring: &Ring) -> bool {
        self.bonds.len() == ring.bonds.len() && self.bonds.iter().all(|&b| ring.contains_bond(b))
    }
}

/// `sketcherMinimizerMolecule`
#[derive(Debug, Clone, Default)]
pub struct Molecule {
    pub atoms: Vec<AtomId>,
    pub bonds: Vec<BondId>,
    pub rings: Vec<RingId>,
    pub proximity_relations: Vec<BondId>,
    pub fragments: Vec<FragmentId>,
    pub main_fragment: Option<FragmentId>,
    pub fixed: bool,
    pub has_fixed_fragments: bool,
    pub has_constrained_fragments: bool,
    pub need_to_align_non_ring_atoms: bool,
    pub need_to_align_whole_molecule: bool,
    pub is_placed: bool,
    pub require_minimization: bool,
}

/// the arena that owns the objects of one sketcher
#[derive(Debug, Clone, Default)]
pub struct Model {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub rings: Vec<Ring>,
    pub molecules: Vec<Molecule>,
    pub fragments: Vec<super::fragment::Fragment>,
}

impl Model {
    pub fn new() -> Self {
        Model::default()
    }

    pub fn add_molecule(&mut self) -> MolId {
        self.molecules.push(Molecule::default());
        self.molecules.len() - 1
    }

    /// `addNewAtom`
    pub fn add_atom(&mut self, mol: MolId) -> AtomId {
        self.atoms.push(Atom {
            molecule: Some(mol),
            ..Atom::default()
        });
        let id = self.atoms.len() - 1;
        self.molecules[mol].atoms.push(id);
        id
    }

    /// `addNewBond`
    pub fn add_bond(&mut self, mol: MolId, at1: AtomId, at2: AtomId) -> BondId {
        self.bonds.push(Bond::new(at1, at2));
        let id = self.bonds.len() - 1;
        self.molecules[mol].bonds.push(id);
        id
    }

    /// `setCoordinates`: the position is rounded to two decimals
    pub fn set_coordinates(&mut self, a: AtomId, mut coords: PointF) {
        coords.round(2);
        let at = &mut self.atoms[a];
        at.coordinates = coords;
        at.coordinates_set = true;
    }

    pub fn coordinates(&self, a: AtomId) -> PointF {
        self.atoms[a].coordinates
    }

    // ---- atom ----

    /// `expectedValence`
    pub fn expected_valence(atomic_number: i32) -> i32 {
        match atomic_number {
            1 => 1,
            5 => 3,
            6 => 4,
            7 => 3,
            8 => 2,
            9 => 1,
            14 => 4,
            15 => 3,
            16 => 2,
            17 => 1,
            34 => 2,
            35 => 1,
            53 => 1,
            _ => 4,
        }
    }

    /// `findHsNumber`
    pub fn find_hs_number(&self, a: AtomId) -> i32 {
        let at = &self.atoms[a];
        let mut valence = at.valence;
        if valence == -10 {
            valence = Self::expected_valence(at.atomic_number);
        }
        let n_bond_orders: i32 = at.bonds.iter().map(|&b| self.bonds[b].bond_order).sum();
        // sulfite & sulfate, phosphorus
        if at.atomic_number == 16 || at.atomic_number == 15 {
            let mut n_os = 0;
            for i in 0..at.neighbors.len() {
                if self.atoms[at.neighbors[i]].atomic_number == 8
                    && self.bonds[at.bonds[i]].bond_order == 2
                {
                    n_os += 1;
                }
            }
            let limit = if at.atomic_number == 16 { 3 } else { 2 };
            if n_os < limit {
                valence += n_os * 2;
            }
        }
        (valence - n_bond_orders + at.charge).clamp(0, 4)
    }

    /// `shareARing`: a ring shared by the two atoms, a non-macrocycle if possible
    pub fn share_a_ring(&self, atom1: AtomId, atom2: AtomId) -> Option<RingId> {
        let (r1, r2) = (&self.atoms[atom1].rings, &self.atoms[atom2].rings);
        if r1.is_empty() || r2.is_empty() {
            return None;
        }
        for &ring in r1 {
            if self.rings[ring].is_macrocycle() {
                continue;
            }
            if r2.contains(&ring) {
                return Some(ring);
            }
        }
        r1.iter().copied().find(|ring| r2.contains(ring))
    }

    /// `bondTo`
    pub fn bond_to(&self, a: AtomId, other: AtomId) -> Option<BondId> {
        let at = &self.atoms[a];
        at.neighbors
            .iter()
            .position(|&n| n == other)
            .map(|i| at.bonds[i])
    }

    /// `hasNoStereoActiveBonds`
    pub fn has_no_stereo_active_bonds(&self, a: AtomId) -> bool {
        !self.atoms[a].bonds.iter().any(|&b| self.bond_is_stereo(b))
    }

    /// `canBeChiral`
    pub fn can_be_chiral(&self, a: AtomId) -> bool {
        let at = &self.atoms[a];
        let n = at.neighbors.len();
        if at.atomic_number == 16 && n == 3 {
            return true;
        }
        if at.atomic_number == 7 && (n == 3 || n == 4) {
            return true;
        }
        if n != 3 && n != 4 {
            return false;
        }
        n as i32 + at.implicit_hs == 4
    }

    /// `clockwiseOrderedNeighbors`: ties are broken by the atom index (the C++ code compares the pointers)
    pub fn clockwise_ordered_neighbors(&self, a: AtomId) -> Vec<AtomId> {
        let at = &self.atoms[a];
        let mut ranked: Vec<(f32, AtomId)> = Vec::with_capacity(at.neighbors.len());
        for &neighbor in &at.neighbors {
            let mut new_angle = signed_angle(
                self.atoms[at.neighbors[0]].coordinates,
                at.coordinates,
                self.atoms[neighbor].coordinates,
            );
            if new_angle.is_nan() {
                new_angle = 361.0;
            } else if new_angle < 0.0 {
                new_angle += 360.0;
            }
            ranked.push((new_angle, neighbor));
        }
        ranked.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap().then(x.1.cmp(&y.1)));
        ranked.into_iter().map(|r| r.1).collect()
    }

    /// `getSingleAdditionVector`: points to the most free region around the atom
    pub fn single_addition_vector(&self, a: AtomId) -> PointF {
        let at = &self.atoms[a];
        let mut out = PointF::new(0.0, 0.0);
        let mut totalf = 0.0f32;
        if !at.neighbors.is_empty() {
            for &neighbor in &at.neighbors {
                let mut f = 1.0f32;
                if self.share_a_ring(a, neighbor).is_some() {
                    f = 4.0;
                }
                out += f * (self.atoms[neighbor].coordinates - at.coordinates);
                totalf += f;
            }
            out = out / totalf;
        }
        out * -1.0
    }

    /// static `getSingleAdditionVector(ats)`
    pub fn single_addition_vector_of(&self, ats: &[AtomId]) -> PointF {
        let mut out = PointF::new(0.0, 0.0);
        let mut n = 0;
        for &a in ats {
            for &ne in &self.atoms[a].neighbors {
                if self.atoms[ne].neighbors.len() > 1 && !ats.contains(&ne) {
                    out += self.atoms[ne].coordinates - self.atoms[a].coordinates;
                    n += 1;
                }
            }
        }
        if n > 0 {
            out = out / n as f32;
        } else {
            return PointF::new(50.0, 0.0);
        }
        out * -1.0
    }

    /// `isMetal`
    pub fn is_metal(z: u32) -> bool {
        matches!(z,
            3 | 4 | 11 | 12 | 19 | 20 | 37 | 38 | 55 | 56 | 87 | 88
            | 21..=30 | 39..=48 | 72..=80 | 104..=112 | 57..=71 | 89..=103
            | 13 | 31 | 49 | 81 | 32 | 50 | 82 | 51 | 83 | 84)
    }

    /// `mirrorCoordinates`: mirrors the atom at the line through the bond
    pub fn mirror_coordinates(&mut self, at: AtomId, bond: BondId) {
        let b = &self.bonds[bond];
        let p = mirror_point(
            self.atoms[at].coordinates,
            self.atoms[b.start_atom].coordinates,
            self.atoms[b.end_atom].coordinates,
        );
        self.set_coordinates(at, p);
    }

    /// `hasValid3DCoordinates`
    pub fn has_valid_3d_coordinates(&self, a: AtomId) -> bool {
        let at = &self.atoms[a];
        at.x3d < INVALID_COORDINATES && at.y3d < INVALID_COORDINATES && at.z3d < INVALID_COORDINATES
    }

    /// `getSubmolecule`: the atoms reachable from `a` without passing `excluded` (breadth first)
    pub fn submolecule(&self, a: AtomId, excluded: AtomId) -> Vec<AtomId> {
        let mut sub = Vec::new();
        let mut visited: HashSet<AtomId> = HashSet::new();
        visited.insert(excluded);
        let mut q = VecDeque::new();
        q.push_back(a);
        visited.insert(a);
        while let Some(atom) = q.pop_front() {
            sub.push(atom);
            for &n in &self.atoms[atom].neighbors {
                if visited.insert(n) {
                    q.push_back(n);
                }
            }
        }
        sub
    }

    // ---- bond ----

    pub fn bond_is_in_small_ring(&self, b: BondId) -> bool {
        self.bonds[b]
            .rings
            .iter()
            .any(|&r| !self.rings[r].is_macrocycle())
    }

    pub fn bond_is_in_macrocycle(&self, b: BondId) -> bool {
        self.bonds[b]
            .rings
            .iter()
            .any(|&r| self.rings[r].is_macrocycle())
    }

    pub fn bond_is_terminal(&self, b: BondId) -> bool {
        let bd = &self.bonds[b];
        self.atoms[bd.start_atom].bonds.len() == 1 || self.atoms[bd.end_atom].bonds.len() == 1
    }

    /// `isStereo`: a double bond with a Z/E form
    pub fn bond_is_stereo(&self, b: BondId) -> bool {
        let bd = &self.bonds[b];
        if bd.bond_order != 2 || bd.ignore_ze {
            return false;
        }
        match self.share_a_ring(bd.start_atom, bd.end_atom) {
            Some(r) => self.rings[r].is_macrocycle(),
            None => true,
        }
    }

    /// `isInterFragment`: is the bond a rotatable single bond between two non-terminal atoms
    pub fn bond_is_inter_fragment(&self, b: BondId) -> bool {
        let bd = &self.bonds[b];
        if self.atoms[bd.start_atom].bonds.len() == 1 || self.atoms[bd.end_atom].bonds.len() == 1 {
            return false;
        }
        if self.share_a_ring(bd.start_atom, bd.end_atom).is_some() {
            return false;
        }
        !self.bond_is_stereo(b)
    }

    /// `flip`: mirrors the smaller side of the bond at the bond and toggles the wedges inside it
    pub fn bond_flip(&mut self, b: BondId) {
        let (start, end) = (self.bonds[b].start_atom, self.bonds[b].end_atom);
        let mol = self.atoms[start].molecule.expect("atom in a molecule");
        let total = self.molecules[mol].atoms.len();
        let mut atoms = self.submolecule(start, end);
        if atoms.len() > total / 2 {
            atoms = self.submolecule(end, start);
        }
        let all_bonds = self.molecules[mol].bonds.clone();
        for &atom in &atoms {
            self.mirror_coordinates(atom, b);
        }
        for bond in all_bonds {
            let (s, e) = (self.bonds[bond].start_atom, self.bonds[bond].end_atom);
            if atoms.contains(&s) && atoms.contains(&e) {
                self.bonds[bond].is_wedge = !self.bonds[bond].is_wedge;
            }
        }
    }

    // ---- ring ----

    pub fn ring_find_center(&self, r: RingId) -> PointF {
        let ring = &self.rings[r];
        let mut o = PointF::new(0.0, 0.0);
        for &a in &ring.atoms {
            o += self.atoms[a].coordinates;
        }
        o / ring.atoms.len() as f32
    }

    pub fn ring_is_benzene(&self, r: RingId) -> bool {
        let ring = &self.rings[r];
        if ring.atoms.len() != 6 || ring.atoms.iter().any(|&a| self.atoms[a].atomic_number != 6) {
            return false;
        }
        ring.atoms.iter().all(|&a| {
            self.atoms[a]
                .bonds
                .iter()
                .any(|&b| self.bonds[b].bond_order == 2)
        })
    }

    /// `isAromatic`: not chemically accurate, but good enough for the minimizer
    pub fn ring_is_aromatic(&self, r: RingId) -> bool {
        let ring = &self.rings[r];
        let bonds = ring.bonds.len();
        let double_bonds = ring
            .bonds
            .iter()
            .filter(|&&b| self.bonds[b].bond_order == 2)
            .count();
        let mut nso_count = 0;
        for &a in &ring.atoms {
            let an = self.atoms[a].atomic_number;
            let double_bound = self.atoms[a]
                .bonds
                .iter()
                .any(|&b| self.bonds[b].bond_order == 2);
            if !double_bound && (an == 8 || an == 7 || an == 16) {
                nso_count += 1;
            }
        }
        (bonds == 6 && double_bonds == 3) || (bonds == 5 && double_bonds == 2 && nso_count == 1)
    }

    pub fn ring_is_fused_with(&self, r: RingId, other: RingId) -> bool {
        self.rings[r].fused_with.contains(&other)
    }

    pub fn ring_fusion_atoms_with(&self, r: RingId, other: RingId) -> Vec<AtomId> {
        let ring = &self.rings[r];
        ring.fused_with
            .iter()
            .position(|&f| f == other)
            .map(|i| ring.fusion_atoms[i].clone())
            .unwrap_or_default()
    }

    /// `contains`: is the point inside the ring polygon (crossing number)
    pub fn ring_contains_point(&self, r: RingId, p: PointF) -> bool {
        let mut n = 0;
        for &b in &self.rings[r].bonds {
            let s = self.atoms[self.bonds[b].start_atom].coordinates;
            let e = self.atoms[self.bonds[b].end_atom].coordinates;
            if (p.y < s.y && p.y > e.y) || (p.y > s.y && p.y < e.y) {
                let mut v = e - s;
                if v.y > SKETCHER_EPSILON || v.y < -SKETCHER_EPSILON {
                    v = v * ((p.y - s.y) / v.y);
                    v += s;
                    if p.x > v.x {
                        n += 1;
                    }
                }
            }
        }
        n % 2 != 0
    }

    // ---- molecule ----

    pub fn total_charge(&self, mol: MolId) -> i32 {
        self.molecules[mol]
            .atoms
            .iter()
            .map(|&a| self.atoms[a].charge)
            .sum()
    }

    /// `boundingBox`: (min, max) of the atom coordinates
    pub fn bounding_box(&self, mol: MolId) -> (PointF, PointF) {
        let atoms = &self.molecules[mol].atoms;
        let (mut min, mut max) = (PointF::new(0.0, 0.0), PointF::new(0.0, 0.0));
        if let Some(&first) = atoms.first() {
            min = self.atoms[first].coordinates;
            max = min;
            for &a in atoms {
                let c = self.atoms[a].coordinates;
                if c.x < min.x {
                    min.x = c.x;
                }
                if c.y < min.y {
                    min.y = c.y;
                }
                if c.x > max.x {
                    max.x = c.x;
                }
                if c.y > max.y {
                    max.y = c.y;
                }
            }
        }
        (min, max)
    }

    pub fn molecule_center(&self, mol: MolId) -> PointF {
        let atoms = &self.molecules[mol].atoms;
        if atoms.is_empty() {
            return PointF::new(0.0, 0.0);
        }
        let mut c = PointF::new(0.0, 0.0);
        for &a in atoms {
            c += self.atoms[a].coordinates;
        }
        c / atoms.len() as f32
    }

    /// `assignBondsAndNeighbors`
    pub fn assign_bonds_and_neighbors(&mut self, mol: MolId) {
        let atoms = self.molecules[mol].atoms.clone();
        let bonds = self.molecules[mol].bonds.clone();
        for &a in &atoms {
            let at = &mut self.atoms[a];
            at.bonds.clear();
            at.neighbors.clear();
            at.rings.clear();
        }
        for &b in &bonds {
            self.bonds[b].rings.clear();
            let (s, e) = (self.bonds[b].start_atom, self.bonds[b].end_atom);
            self.atoms[s].bonds.push(b);
            self.atoms[e].neighbors.push(s);
            self.atoms[e].bonds.push(b);
            self.atoms[s].neighbors.push(e);
        }
        for &a in &atoms {
            if self.atoms[a].implicit_hs == -1 {
                self.atoms[a].implicit_hs = self.find_hs_number(a);
            }
        }
    }

    /// `forceUpdateStruct`: neighbors, bonds and the ring system of the molecule
    pub fn force_update_struct(&mut self, mol: MolId) {
        self.assign_bonds_and_neighbors(mol);
        self.find_rings(mol);
        let atoms = self.molecules[mol].atoms.clone();
        for b in self.molecules[mol].bonds.clone() {
            for ring in self.bonds[b].rings.clone() {
                for end in [self.bonds[b].start_atom, self.bonds[b].end_atom] {
                    if !self.atoms[end].rings.contains(&ring) {
                        self.atoms[end].rings.push(ring);
                    }
                }
            }
        }
        for &a in &atoms {
            for ring in self.atoms[a].rings.clone() {
                self.rings[ring].atoms.push(a);
            }
        }
    }

    /// `findRings`: the rings are closed from breadth first searches starting at every bond
    pub fn find_rings(&mut self, mol: MolId) {
        self.molecules[mol].rings.clear();
        let bonds = self.molecules[mol].bonds.clone();
        for &bond in &bonds {
            for &b in &bonds {
                let bd = &mut self.bonds[b];
                bd.sssr_visited = false;
                bd.sssr_parent = None;
                bd.sssr_parent_at_start = true;
            }
            let mut q = VecDeque::new();
            self.bonds[bond].sssr_visited = true;
            q.push_back(bond);
            let mut closed_ring = false;
            while !q.is_empty() && !closed_ring {
                let last_bond = q.pop_front().unwrap();
                let pivot = if self.bonds[last_bond].sssr_parent_at_start {
                    self.bonds[last_bond].end_atom
                } else {
                    self.bonds[last_bond].start_atom
                };
                for j in 0..self.atoms[pivot].bonds.len() {
                    let next = self.atoms[pivot].bonds[j];
                    if next == last_bond {
                        continue;
                    }
                    if self.bonds[next].sssr_visited {
                        if next == bond {
                            let ring = self.close_ring(last_bond);
                            self.add_ring(mol, ring);
                            closed_ring = true;
                        }
                    } else {
                        if self.bonds[next].end_atom == pivot {
                            self.bonds[next].sssr_parent_at_start = false;
                        }
                        self.bonds[next].sssr_parent = Some(last_bond);
                        self.bonds[next].sssr_visited = true;
                        q.push_back(next);
                    }
                }
            }
        }
        for ring in self.molecules[mol].rings.clone() {
            for b in self.rings[ring].bonds.clone() {
                self.bonds[b].rings.push(ring);
            }
        }
    }

    /// `closeRing`: a new ring from the parent chain of the bond
    fn close_ring(&mut self, bond: BondId) -> Ring {
        let mut ring = Ring::default();
        let mut last = Some(bond);
        while let Some(b) = last {
            ring.bonds.push(b);
            last = self.bonds[b].sssr_parent;
        }
        ring
    }

    /// `addRing`: only rings that are not in the list yet are kept
    fn add_ring(&mut self, mol: MolId, ring: Ring) {
        if self.molecules[mol]
            .rings
            .iter()
            .any(|&r| self.rings[r].same_as(&ring))
        {
            return;
        }
        self.rings.push(ring);
        let id = self.rings.len() - 1;
        self.molecules[mol].rings.push(id);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fmt::Write;

    fn bits(f: f32) -> u32 {
        f.to_bits()
    }

    fn list(ids: &[usize], all: &[usize]) -> String {
        ids.iter()
            .map(|i| format!(" {}", all.iter().position(|x| x == i).unwrap()))
            .collect()
    }

    /// builds the molecule of an input block and prints what the driver `tools/coordgen_drv/model.cpp` prints
    fn run(input: &str) -> String {
        let mut model = Model::new();
        let mol = model.add_molecule();
        let mut lines = input.lines();
        let _ = lines.next();
        for line in lines {
            let f: Vec<&str> = line.split(' ').collect();
            match f[0] {
                "A" => {
                    let a = model.add_atom(mol);
                    model.atoms[a].atomic_number = f[1].parse().unwrap();
                    model.atoms[a].charge = f[2].parse().unwrap();
                    let x = f32::from_bits(u32::from_str_radix(f[3], 16).unwrap());
                    let y = f32::from_bits(u32::from_str_radix(f[4], 16).unwrap());
                    model.set_coordinates(a, PointF::new(x, y));
                }
                "B" => {
                    let b = model.add_bond(mol, f[1].parse().unwrap(), f[2].parse().unwrap());
                    model.bonds[b].bond_order = f[3].parse().unwrap();
                }
                _ => {}
            }
        }
        model.force_update_struct(mol);
        let atoms = model.molecules[mol].atoms.clone();
        let bonds = model.molecules[mol].bonds.clone();
        let rings = model.molecules[mol].rings.clone();
        let mut out = String::new();
        for (i, &a) in atoms.iter().enumerate() {
            let at = &model.atoms[a];
            write!(
                out,
                "a {i} hs {} nb{} bd{} rg{}",
                at.implicit_hs,
                list(&at.neighbors, &atoms),
                list(&at.bonds, &bonds),
                list(&at.rings, &rings)
            )
            .unwrap();
            write!(
                out,
                " chiral {} stereoActive {}",
                i32::from(model.can_be_chiral(a)),
                i32::from(!model.has_no_stereo_active_bonds(a))
            )
            .unwrap();
            write!(
                out,
                " cw{}",
                list(&model.clockwise_ordered_neighbors(a), &atoms)
            )
            .unwrap();
            if let Some(&n0) = at.neighbors.first() {
                write!(out, " sub{}", list(&model.submolecule(a, n0), &atoms)).unwrap();
            }
            let v = model.single_addition_vector(a);
            writeln!(out, " sav {:08x} {:08x}", bits(v.x), bits(v.y)).unwrap();
        }
        for (i, &b) in bonds.iter().enumerate() {
            writeln!(
                out,
                "b {i} small {} macro {} term {} inter {} stereo {} rings{}",
                i32::from(model.bond_is_in_small_ring(b)),
                i32::from(model.bond_is_in_macrocycle(b)),
                i32::from(model.bond_is_terminal(b)),
                i32::from(model.bond_is_inter_fragment(b)),
                i32::from(model.bond_is_stereo(b)),
                list(&model.bonds[b].rings, &rings)
            )
            .unwrap();
        }
        for (i, &r) in rings.iter().enumerate() {
            let c = model.ring_find_center(r);
            writeln!(
                out,
                "r {i} size {} benz {} arom {} mac {} atoms{} bonds{} center {:08x} {:08x}",
                model.rings[r].size(),
                i32::from(model.ring_is_benzene(r)),
                i32::from(model.ring_is_aromatic(r)),
                i32::from(model.rings[r].is_macrocycle()),
                list(&model.rings[r].atoms, &atoms),
                list(&model.rings[r].bonds, &bonds),
                bits(c.x),
                bits(c.y)
            )
            .unwrap();
        }
        let (mn, mx) = model.bounding_box(mol);
        let ce = model.molecule_center(mol);
        writeln!(
            out,
            "m charge {} bbox {:08x} {:08x} {:08x} {:08x} center {:08x} {:08x}",
            model.total_charge(mol),
            bits(mn.x),
            bits(mn.y),
            bits(mx.x),
            bits(mx.y),
            bits(ce.x),
            bits(ce.y)
        )
        .unwrap();
        for (i, &b) in bonds.iter().enumerate() {
            if !model.bonds[b].rings.is_empty() || model.bond_is_terminal(b) {
                continue;
            }
            let saved: Vec<PointF> = atoms.iter().map(|&a| model.atoms[a].coordinates).collect();
            model.bond_flip(b);
            write!(out, "flip {i}").unwrap();
            for &a in &atoms {
                let c = model.atoms[a].coordinates;
                write!(out, " {:08x} {:08x}", bits(c.x), bits(c.y)).unwrap();
            }
            write!(out, " wedge").unwrap();
            for &bb in &bonds {
                write!(out, " {}", i32::from(model.bonds[bb].is_wedge)).unwrap();
            }
            out.push('\n');
            for (&a, c) in atoms.iter().zip(saved) {
                model.atoms[a].coordinates = c;
            }
            for &bb in &bonds {
                model.bonds[bb].is_wedge = false;
            }
        }
        out
    }

    /// the model must answer what the C++ classes answer (`tests/data/coordgen_model.txt`)
    #[test]
    fn model_matches_cpp() {
        let text = include_str!("../../../../../tests/data/coordgen_model.txt");
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
        assert!(n > 60, "{n}");
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }
}
