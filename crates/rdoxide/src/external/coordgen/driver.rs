//! The driver of coordgenlibs (`sketcherMinimizer`): canonical ordering, splitting into molecules, fragment
//! initialization, the whole `runGenerateCoordinates` pipeline (fragments, clash avoidance, best rotation, flip,
//! arrangement of several molecules). Residue (protein) layouts, proximity relations and the Maestro reader are not
//! part of the RDKit use of the library and are not ported.

use std::collections::VecDeque;

use super::clashes::{alpha_cs, amino_ns, cheto_cs};
use super::maths::{BONDLENGTH, PointF, SKETCHER_EPSILON, round_to_two_decimal_digits};
use super::minimizer::Minimizer;
use super::model::{AtomId, BondId, FragmentId, Model, MolId};
use super::templates::morgan_scores;

const MAX_NUMBER_OF_RINGS: usize = 40;
const PI: f64 = std::f64::consts::PI;

/// `sketcherMinimizer`
#[derive(Debug, Clone)]
pub struct Sketcher {
    pub minimizer: Minimizer,
    pub atoms: Vec<AtomId>,
    pub bonds: Vec<BondId>,
    pub molecules: Vec<MolId>,
    pub fragments: Vec<FragmentId>,
    pub independent_fragments: Vec<FragmentId>,
}

impl Sketcher {
    /// `sketcherMinimizer(precision)`
    pub fn new(precision: f32) -> Self {
        Sketcher {
            minimizer: Minimizer {
                precision,
                ..Minimizer::default()
            },
            atoms: Vec::new(),
            bonds: Vec::new(),
            molecules: Vec::new(),
            fragments: Vec::new(),
            independent_fragments: Vec::new(),
        }
    }

    /// `canonicalOrdering`: atoms and bonds in an order that does not depend on the input order
    fn canonical_ordering(m: &mut Model, mol: MolId) {
        let atoms = m.molecules[mol].atoms.clone();
        let bonds = m.molecules[mol].bonds.clone();
        for (i, &a) in atoms.iter().enumerate() {
            m.atoms[a].general_use_n = i as i32;
        }
        let pairs: Vec<(usize, usize)> = bonds
            .iter()
            .map(|&b| {
                (
                    m.atoms[m.bonds[b].start_atom].general_use_n as usize,
                    m.atoms[m.bonds[b].end_atom].general_use_n as usize,
                )
            })
            .collect();
        let (_, mut scores) = morgan_scores(atoms.len(), &pairs);
        if scores.len() != atoms.len() {
            return;
        }
        for (i, s) in scores.iter_mut().enumerate() {
            *s *= 100;
            *s += m.atoms[atoms[i]].atomic_number;
        }
        for &a in &atoms {
            m.atoms[a].neighbors.clear();
            m.atoms[a].bonds.clear();
        }
        for &b in &bonds {
            let (s, e) = (m.bonds[b].start_atom, m.bonds[b].end_atom);
            m.atoms[s].neighbors.push(e);
            m.atoms[e].neighbors.push(s);
            m.atoms[s].bonds.push(b);
            m.atoms[e].bonds.push(b);
        }
        let mut new_atoms = Vec::new();
        let mut new_bonds = Vec::new();
        for (i, &a) in atoms.iter().enumerate() {
            m.atoms[a].general_use_n = i as i32;
            m.atoms[a].general_use_visited = false;
        }
        for &b in &bonds {
            m.bonds[b].sssr_visited = false;
        }
        loop {
            let mut score_max_i: Option<usize> = None;
            for i in 0..scores.len() {
                if m.atoms[atoms[i]].general_use_visited {
                    continue;
                }
                match score_max_i {
                    None => score_max_i = Some(i),
                    Some(j) if scores[i] > scores[j] => score_max_i = Some(i),
                    _ => {}
                }
            }
            let Some(start) = score_max_i else {
                break;
            };
            let mut q = VecDeque::new();
            q.push_back(atoms[start]);
            m.atoms[atoms[start]].general_use_visited = true;
            while let Some(at) = q.pop_front() {
                new_atoms.push(at);
                loop {
                    let mut neigh_i: Option<usize> = None;
                    for i in 0..m.atoms[at].neighbors.len() {
                        if m.bonds[m.atoms[at].bonds[i]].sssr_visited {
                            continue;
                        }
                        match neigh_i {
                            None => neigh_i = Some(i),
                            Some(j) => {
                                let nj = m.atoms[m.atoms[at].neighbors[j]].general_use_n as usize;
                                let ni = m.atoms[m.atoms[at].neighbors[i]].general_use_n as usize;
                                if scores[nj] < scores[ni] {
                                    neigh_i = Some(i);
                                }
                            }
                        }
                    }
                    let Some(ni) = neigh_i else {
                        break;
                    };
                    let n = m.atoms[at].neighbors[ni];
                    if !m.atoms[n].general_use_visited {
                        m.atoms[n].general_use_visited = true;
                        q.push_back(n);
                    }
                    let b = m.atoms[at].bonds[ni];
                    m.bonds[b].sssr_visited = true;
                    new_bonds.push(b);
                }
            }
        }
        m.molecules[mol].atoms = new_atoms;
        m.molecules[mol].bonds = new_bonds;
    }

    /// `splitIntoMolecules`: one molecule per connected component
    fn split_into_molecules(m: &mut Model, mol: MolId, mols: &mut Vec<MolId>) {
        if m.molecules[mol].atoms.is_empty() {
            mols.push(mol);
            return;
        }
        for a in m.molecules[mol].atoms.clone() {
            m.atoms[a].general_use_visited = false;
        }
        let mut q = VecDeque::new();
        let first = m.molecules[mol].atoms[0];
        q.push_back(first);
        if let Some(&a) = m.molecules[mol].atoms.iter().find(|&&a| !m.atoms[a].hidden) {
            q.push_back(a);
        }
        while let Some(a) = q.pop_front() {
            m.atoms[a].general_use_visited = true;
            for n in m.atoms[a].neighbors.clone() {
                if !m.atoms[n].general_use_visited && !m.atoms[n].hidden {
                    q.push_back(n);
                }
            }
        }
        let any_unvisited = m.molecules[mol]
            .atoms
            .iter()
            .any(|&a| !m.atoms[a].general_use_visited && !m.atoms[a].hidden);
        if !any_unvisited {
            mols.push(mol);
            for &mm in mols.iter() {
                for a in m.molecules[mm].atoms.clone() {
                    m.atoms[a].general_use_visited = false;
                }
            }
            return;
        }
        let new_mol = m.add_molecule();
        let rings = std::mem::take(&mut m.molecules[mol].rings);
        let (moved, kept): (Vec<_>, Vec<_>) = rings
            .into_iter()
            .partition(|&r| !m.atoms[m.rings[r].atoms[0]].general_use_visited);
        m.molecules[mol].rings = kept;
        m.molecules[new_mol].rings = moved;
        let bonds = std::mem::take(&mut m.molecules[mol].bonds);
        let (moved, kept): (Vec<_>, Vec<_>) = bonds
            .into_iter()
            .partition(|&b| !m.atoms[m.bonds[b].start_atom].general_use_visited);
        m.molecules[mol].bonds = kept;
        m.molecules[new_mol].bonds = moved;
        let atoms = std::mem::take(&mut m.molecules[mol].atoms);
        let (moved, kept): (Vec<_>, Vec<_>) = atoms
            .into_iter()
            .partition(|&a| !m.atoms[a].general_use_visited);
        for &a in &moved {
            m.atoms[a].molecule = Some(new_mol);
        }
        m.molecules[mol].atoms = kept;
        m.molecules[new_mol].atoms = moved;
        mols.push(mol);
        Self::split_into_molecules(m, new_mol, mols);
    }

    /// `flagCrossAtoms`
    fn flag_cross_atoms(&self, m: &mut Model) {
        for &a in &self.atoms {
            if m.atoms[a].atomic_number == 16 || m.atoms[a].atomic_number == 15 {
                m.atoms[a].cross_layout = true;
            }
        }
        for &a in &self.atoms {
            if m.atoms[a].cross_layout {
                continue;
            }
            let cross = m.atoms[a]
                .neighbors
                .iter()
                .filter(|&&n| m.atoms[n].neighbors.len() > 3)
                .count();
            if cross > 2 {
                m.atoms[a].cross_layout = true;
            }
        }
    }

    /// `initialize` (zero order bonds, hidden atoms and residues are not supported)
    pub fn initialize(&mut self, m: &mut Model, min_mol: MolId) {
        for &b in &m.molecules[min_mol].bonds {
            assert!(
                m.bonds[b].bond_order != 0,
                "proximity relations (zero order bonds) are not ported"
            );
        }
        let bonds: Vec<BondId> = m.molecules[min_mol]
            .bonds
            .iter()
            .copied()
            .filter(|&b| !m.bonds[b].skip)
            .collect();
        m.molecules[min_mol].bonds = bonds;
        Self::canonical_ordering(m, min_mol);
        self.atoms = m.molecules[min_mol].atoms.clone();
        self.bonds = m.molecules[min_mol].bonds.clone();
        m.force_update_struct(min_mol);
        let mut mols = Vec::new();
        Self::split_into_molecules(m, min_mol, &mut mols);
        self.molecules = mols;
        self.flag_cross_atoms(m);
        self.minimizer.atoms = self.atoms.clone();
        self.minimizer.molecules = self.molecules.clone();
    }

    /// `structurePassSanityCheck`
    fn structure_pass_sanity_check(&self, m: &Model) -> bool {
        !self.atoms.is_empty()
            && self
                .molecules
                .iter()
                .all(|&mol| m.molecules[mol].rings.len() <= MAX_NUMBER_OF_RINGS)
    }

    fn assign_number_of_children_atoms_from_here(m: &mut Model, f: FragmentId) {
        let mut cumulated = 0usize;
        let mut cumulated_ranks = 0.0f32;
        let mut children_atoms = 0usize;
        for child in m.fragments[f].children.clone() {
            Self::assign_number_of_children_atoms_from_here(m, child);
            cumulated += m.fragments[child].number_of_children_atoms;
            cumulated_ranks += m.fragments[child].number_of_children_atoms_rank;
            children_atoms += m.fragments[child].atoms().len();
        }
        m.fragments[f].number_of_children_atoms = cumulated + children_atoms;
        m.fragments[f].number_of_children_atoms_rank =
            0.01f32 * cumulated_ranks + children_atoms as f32;
    }

    fn assign_longest_chain_from_here(m: &mut Model, f: FragmentId) {
        let mut longest = 0.0f32;
        for child in m.fragments[f].children.clone() {
            Self::assign_longest_chain_from_here(m, child);
            if m.fragments[child].longest_chain_from_here > longest {
                longest = m.fragments[child].longest_chain_from_here;
            }
        }
        let mut position_from_parent = PointF::new(0.0, 0.0);
        if let Some(p) = m.fragments[f].parent {
            let end = m.bonds[m.fragments[f].bond_to_parent.expect("parent bond")].end_atom;
            position_from_parent = m.fragments[p]
                .coordinates
                .get(&end)
                .copied()
                .unwrap_or_default();
        }
        m.fragments[f].longest_chain_from_here = longest + position_from_parent.length();
    }

    /// `findFragments` (with `initializeFragments`)
    fn find_fragments(&mut self, m: &mut Model) {
        assert!(!self.molecules.is_empty());
        for mol in self.molecules.clone() {
            m.split_into_fragments(mol);
            if m.molecules[mol].fragments.is_empty() {
                continue;
            }
            self.fragments
                .extend(m.molecules[mol].fragments.iter().copied());
            if let Some(f) = m.molecules[mol].main_fragment {
                self.independent_fragments.push(f);
            }
        }
        if self.fragments.is_empty() {
            return;
        }
        for &indf in &self.independent_fragments {
            Self::assign_number_of_children_atoms_from_here(m, indf);
        }
        for &f in &self.fragments {
            m.initialize_coordinates(f);
        }
        for &indf in &self.independent_fragments {
            Self::assign_longest_chain_from_here(m, indf);
        }
    }

    /// `writeStereoChemistry` and `assignPseudoZ`
    fn write_stereo_chemistry(&self, m: &mut Model) {
        for &a in &self.atoms {
            if m.atoms[a].has_stereochemistry_set {
                m.write_stereo_chemistry(a);
            }
        }
        for &mol in &self.molecules {
            for a in m.molecules[mol].atoms.clone() {
                m.atoms[a].general_use_visited = false;
            }
            while let Some(start) = m.molecules[mol]
                .atoms
                .iter()
                .copied()
                .find(|&a| !m.atoms[a].general_use_visited)
            {
                let mut q = VecDeque::new();
                q.push_back(start);
                while let Some(last) = q.pop_front() {
                    m.atoms[last].general_use_visited = true;
                    for i in 0..m.atoms[last].neighbors.len() {
                        let n = m.atoms[last].neighbors[i];
                        if m.atoms[n].general_use_visited {
                            continue;
                        }
                        let mut z = m.atoms[last].pseudo_z;
                        let b = &m.bonds[m.atoms[last].bonds[i]];
                        if b.has_stereochemistry_display {
                            let forward = (b.start_atom == last && !b.is_reversed)
                                || (b.end_atom == last && b.is_reversed);
                            let backward = (b.start_atom == last && b.is_reversed)
                                || (b.end_atom == last && !b.is_reversed);
                            if b.is_wedge {
                                if forward {
                                    z += 1.0;
                                } else if backward {
                                    z -= 1.0;
                                }
                            } else if forward {
                                z -= 1.0;
                            } else if backward {
                                z += 1.0;
                            }
                        }
                        m.atoms[n].pseudo_z = z;
                        q.push_back(n);
                    }
                }
            }
        }
    }

    /// `maybeFlipPeptides`
    fn maybe_flip_peptides(m: &Model, atoms: &[AtomId], score_x: &mut f32) {
        let cheto = cheto_cs(m, atoms);
        let amino = amino_ns(m, atoms);
        for alpha in alpha_cs(m, atoms, &cheto, &amino) {
            let (mut amino_n, mut cheto_c) = (None, None);
            for &n in &m.atoms[alpha].neighbors {
                if amino.contains(&n) {
                    amino_n = Some(n);
                } else if cheto.contains(&n) {
                    cheto_c = Some(n);
                }
            }
            if let (Some(an), Some(cc)) = (amino_n, cheto_c) {
                let direction = m.atoms[an].coordinates - m.atoms[cc].coordinates;
                const PEPTIDE_SCORE: f32 = 100.0;
                if direction.x > 0.0 {
                    *score_x -= PEPTIDE_SCORE;
                } else {
                    *score_x += PEPTIDE_SCORE;
                }
            }
        }
    }

    /// `maybeFlip`
    fn maybe_flip(&self, m: &mut Model) {
        for &mol in &self.molecules {
            if m.molecules[mol].has_fixed_fragments || m.molecules[mol].has_constrained_fragments {
                continue;
            }
            let atoms = m.molecules[mol].atoms.clone();
            if atoms.len() < 2 {
                continue;
            }
            let (mut score_y, mut score_x) = (0.0f32, 0.0f32);
            Self::maybe_flip_peptides(m, &atoms, &mut score_x);
            let mut cent = PointF::new(0.0, 0.0);
            for &a in &atoms {
                cent += m.atoms[a].coordinates;
            }
            cent = cent / atoms.len() as f32;
            for f in m.molecules[mol].fragments.clone() {
                let rings = m.fragments[f].rings().to_vec();
                if rings.len() == 2 {
                    let (mut r1, mut r2) = (rings[0], rings[1]);
                    let (mut c1, mut c2) = (m.ring_find_center(r1), m.ring_find_center(r2));
                    if c1.x - c2.x > SKETCHER_EPSILON {
                        std::mem::swap(&mut c1, &mut c2);
                        std::mem::swap(&mut r1, &mut r2);
                    }
                    if m.ring_is_benzene(r2) && !m.ring_is_benzene(r1) {
                        score_x -= 20.0;
                    } else if m.ring_is_benzene(r1) && !m.ring_is_benzene(r2) {
                        score_x += 20.0;
                    }
                }
                if rings.len() > 3 {
                    let mut center = PointF::new(0.0, 0.0);
                    let mut weighted_center = PointF::new(0.0, 0.0);
                    let mut total_n = 0usize;
                    let mut total_rings = 0i32;
                    for &r in &rings {
                        let size = m.rings[r].atoms.len();
                        if size < 4 {
                            continue;
                        }
                        let c = m.ring_find_center(r);
                        center += c;
                        weighted_center += c * size as f32;
                        total_n += size;
                        total_rings += 1;
                    }
                    if total_rings != 0 && total_n != 0 {
                        center = center / total_rings as f32;
                        weighted_center = weighted_center / total_n as f32;
                    }
                    if weighted_center.y - center.y < -SKETCHER_EPSILON {
                        score_y += 50.0;
                    } else if weighted_center.y - center.y > SKETCHER_EPSILON {
                        score_y -= 50.0;
                    }
                    if weighted_center.x - center.x < -SKETCHER_EPSILON {
                        score_x += 50.0;
                    } else if weighted_center.x - center.x > SKETCHER_EPSILON {
                        score_x -= 50.0;
                    }
                }
            }
            let (mut minx, mut miny, mut maxx, mut maxy) =
                (9999.0f32, 9999.0f32, -9999.0f32, -9999.0f32);
            for &a in &atoms {
                let (x, y) = (m.atoms[a].coordinates.x, m.atoms[a].coordinates.y);
                if x < minx {
                    minx = x;
                } else if x > maxx {
                    maxx = x;
                }
                if y < miny {
                    miny = y;
                } else if y > maxy {
                    maxy = y;
                }
            }
            let meanx = (maxx + minx) * 0.5;
            let meany = (maxy + miny) * 0.5;
            if meanx - cent.x > SKETCHER_EPSILON {
                score_x -= 0.5;
            } else if meanx - cent.x < -SKETCHER_EPSILON {
                score_x += 0.5;
            }
            if meany - cent.y > SKETCHER_EPSILON {
                score_y += 0.5;
            } else if meany - cent.y < -SKETCHER_EPSILON {
                score_y -= 0.5;
            }
            for b in m.molecules[mol].bonds.clone() {
                let bd = &m.bonds[b];
                if bd.bond_order != 2 {
                    continue;
                }
                let (s, e) = (bd.start_atom, bd.end_atom);
                let (sn, en) = (m.atoms[s].neighbors.len(), m.atoms[e].neighbors.len());
                let diff = if sn == 1 && en > 1 {
                    Some(m.atoms[s].coordinates.y - m.atoms[e].coordinates.y)
                } else if en == 1 && sn > 1 {
                    Some(m.atoms[e].coordinates.y - m.atoms[s].coordinates.y)
                } else {
                    None
                };
                if let Some(diff) = diff {
                    if diff > SKETCHER_EPSILON {
                        score_y += 1.0;
                    } else if diff < -SKETCHER_EPSILON {
                        score_y -= 1.0;
                    }
                }
            }
            if 0.0 > score_y {
                for &a in &atoms {
                    m.atoms[a].coordinates.y = -m.atoms[a].coordinates.y;
                }
            }
            if 0.0 > score_x {
                for &a in &atoms {
                    m.atoms[a].coordinates.x = -m.atoms[a].coordinates.x;
                }
            }
        }
    }

    /// `addToVector`: (weight, angle) pairs ordered by angle
    fn add_to_vector(weight: f32, angle: f32, angles: &mut Vec<(f32, f32)>) {
        let mut angle = round_to_two_decimal_digits(angle);
        while angle <= 0.0 {
            angle += PI as f32;
        }
        for i in 0..angles.len() {
            if angles[i].1 < angle - SKETCHER_EPSILON {
                if i == angles.len() - 1 {
                    angles.push((weight, angle));
                    break;
                }
            } else if angles[i].1 - angle < SKETCHER_EPSILON
                && angles[i].1 - angle > -SKETCHER_EPSILON
            {
                angles[i].0 += weight;
                break;
            } else {
                angles.insert(i, (weight, angle));
                break;
            }
        }
        if angles.is_empty() {
            angles.push((weight, angle));
        }
    }

    /// `bestRotation`: rotates every molecule so that most bonds point in preferred directions
    fn best_rotation(&self, m: &mut Model) {
        for &mol in &self.molecules {
            let mut angles: Vec<(f32, f32)> = Vec::new();
            if m.molecules[mol].has_fixed_fragments || m.molecules[mol].has_constrained_fragments {
                continue;
            }
            let atoms = m.molecules[mol].atoms.clone();
            // peptide backbones are made horizontal
            let cheto = cheto_cs(m, &atoms);
            let amino = amino_ns(m, &atoms);
            for alpha in alpha_cs(m, &atoms, &cheto, &amino) {
                let (mut amino_n, mut cheto_c) = (None, None);
                for &n in &m.atoms[alpha].neighbors {
                    if amino.contains(&n) {
                        amino_n = Some(n);
                    } else if cheto.contains(&n) {
                        cheto_c = Some(n);
                    }
                }
                if let (Some(an), Some(cc)) = (amino_n, cheto_c) {
                    let direction = m.atoms[an].coordinates - m.atoms[cc].coordinates;
                    Self::add_to_vector(1000.0, (-direction.y).atan2(direction.x), &mut angles);
                }
            }
            let increment = (PI / 6.0) as f32;
            for &a in &atoms {
                if !m.atoms[a].rings.is_empty() {
                    continue;
                }
                let nb = m.atoms[a].neighbors.clone();
                if nb.len() > 1 {
                    for i in 0..nb.len() - 1 {
                        for j in i + 1..nb.len() {
                            let mut weight = 6.0f32;
                            if m.atoms[nb[i]].neighbors.len() != 1 {
                                weight += 2.0;
                            }
                            if m.atoms[nb[j]].neighbors.len() != 1 {
                                weight += 2.0;
                            }
                            if m.atoms[nb[j]].atomic_number == 6 {
                                weight += 1.0;
                            }
                            if m.atoms[nb[j]].atomic_number == 6 {
                                weight += 1.0;
                            }
                            if m.atoms[nb[i]].charge == 0 {
                                weight += 1.0;
                            }
                            if m.atoms[nb[j]].charge == 0 {
                                weight += 1.0;
                            }
                            let p = m.atoms[nb[i]].coordinates - m.atoms[nb[j]].coordinates;
                            Self::add_to_vector(weight, (-p.y).atan2(p.x), &mut angles);
                        }
                    }
                }
            }
            for b in m.molecules[mol].bonds.clone() {
                let bd = &m.bonds[b];
                let (s, e) = (bd.start_atom, bd.end_atom);
                let p = m.atoms[e].coordinates - m.atoms[s].coordinates;
                let mut angle = round_to_two_decimal_digits((-p.y).atan2(p.x));
                while angle <= 0.0 {
                    angle += PI as f32;
                }
                let mut last_angle = angle;
                for i in 0..6 {
                    let mut weight = if i == 1 || i == 5 {
                        5.0f32
                    } else if i == 0 || i == 3 {
                        1.5
                    } else {
                        1.0
                    };
                    let (sn, en) = (m.atoms[s].neighbors.len(), m.atoms[e].neighbors.len());
                    if bd.bond_order == 2 && i == 3 && (sn == 1 || en == 1) {
                        weight += 1.5;
                    }
                    if sn == 1 && en == 1 && i == 0 {
                        weight += 10.0;
                    }
                    Self::add_to_vector(weight, last_angle, &mut angles);
                    last_angle += increment;
                    if f64::from(last_angle) > PI {
                        last_angle -= PI as f32;
                    }
                }
            }
            for f in m.molecules[mol].fragments.clone() {
                let in_plane = m.fragments[f].rings().to_vec();
                match in_plane.len() {
                    2 => {
                        let mut p =
                            m.ring_find_center(in_plane[1]) - m.ring_find_center(in_plane[0]);
                        p.normalize();
                        Self::add_to_vector(25.0, (-p.y).atan2(p.x), &mut angles);
                    }
                    3 => {
                        let mut r1 = m.ring_find_center(in_plane[0]);
                        let mut r2 = m.ring_find_center(in_plane[1]);
                        for &r in &in_plane {
                            let ring = &m.rings[r];
                            if ring.fused_with.len() == 2
                                && ring.fusion_atoms[0].len() == 2
                                && ring.fusion_atoms[1].len() == 2
                            {
                                let c = |a: AtomId| m.atoms[a].coordinates;
                                r1 =
                                    (c(ring.fusion_atoms[0][0]) + c(ring.fusion_atoms[0][1])) * 0.5;
                                r2 =
                                    (c(ring.fusion_atoms[1][0]) + c(ring.fusion_atoms[1][1])) * 0.5;
                                break;
                            }
                        }
                        let p = r2 - r1;
                        Self::add_to_vector(50.0, (-p.y).atan2(p.x), &mut angles);
                    }
                    _ => {
                        for &r in in_plane.iter().filter(|&&r| m.rings[r].atoms.len() == 6) {
                            for fusion_ats in &m.rings[r].fusion_atoms {
                                if fusion_ats.len() == 2 {
                                    let p = m.atoms[fusion_ats[0]].coordinates
                                        - m.atoms[fusion_ats[1]].coordinates;
                                    let angle = (f64::from((-p.y).atan2(p.x)) - PI * 0.5) as f32;
                                    Self::add_to_vector(25.0, angle, &mut angles);
                                }
                            }
                        }
                    }
                }
            }
            if angles.len() > 1
                && f64::from(angles[angles.len() - 1].1 - angles[0].1)
                    >= PI - 2.0 * f64::from(SKETCHER_EPSILON)
            {
                angles[0].0 += angles[angles.len() - 1].0;
                angles.pop();
            }
            if !angles.is_empty() {
                let mut best_i = 0;
                for i in 0..angles.len() {
                    if angles[i].0 > angles[best_i].0 {
                        best_i = i;
                    }
                }
                let s = -angles[best_i].1.sin();
                let c = angles[best_i].1.cos();
                let mut center = PointF::new(0.0, 0.0);
                for &a in &atoms {
                    center += m.atoms[a].coordinates;
                }
                if !atoms.is_empty() {
                    center = center / atoms.len() as f32;
                }
                for &a in &atoms {
                    let mut v = m.atoms[a].coordinates - center;
                    v.rotate(s, c);
                    m.set_coordinates(a, center + v);
                }
            }
        }
    }

    /// `exploreGridAround` (without residues)
    #[allow(clippy::too_many_arguments)]
    fn explore_grid_around(
        &self,
        m: &Model,
        center_of_grid: PointF,
        levels: u32,
        grid_d: f32,
        dx: f32,
        dy: f32,
        distance_from_atoms: f32,
    ) -> PointF {
        let direction = PointF::new(0.0, 1.0);
        let mut place_next_to = center_of_grid;
        let mut distance_from_atoms = distance_from_atoms;
        for i in 0..levels {
            let step = (1 + i) as f32 * grid_d;
            let top = center_of_grid + PointF::new(0.0, step);
            let bottom = center_of_grid + PointF::new(0.0, -step);
            let right = center_of_grid + PointF::new(step, 0.0);
            let left = center_of_grid + PointF::new(-step, 0.0);
            let mut points = vec![center_of_grid, right, left, bottom, top];
            for j in 0..i {
                let o = grid_d * (j + 1) as f32;
                points.push(right + PointF::new(0.0, o));
                points.push(right - PointF::new(0.0, o));
                points.push(left + PointF::new(0.0, o));
                points.push(left - PointF::new(0.0, o));
                points.push(bottom + PointF::new(o, 0.0));
                points.push(bottom - PointF::new(o, 0.0));
                points.push(top + PointF::new(o, 0.0));
                points.push(top - PointF::new(o, 0.0));
            }
            points.push(center_of_grid + PointF::new(step, step));
            points.push(center_of_grid + PointF::new(step, -step));
            points.push(center_of_grid + PointF::new(-step, step));
            points.push(center_of_grid + PointF::new(-step, -step));
            let mut no_clash = true;
            if distance_from_atoms < 0.0 {
                distance_from_atoms = (f64::from(BONDLENGTH) * 1.8) as f32;
            }
            let mut direction_normal = PointF::new(-direction.y, direction.x);
            direction_normal.normalize();
            for &pc in &points {
                no_clash = true;
                let point = pc - center_of_grid;
                place_next_to = point.y * direction + point.x * direction_normal + center_of_grid;
                'mols: for &mm in &self.molecules {
                    if !m.molecules[mm].is_placed {
                        continue;
                    }
                    for &a in &m.molecules[mm].atoms {
                        let dist = distance_from_atoms;
                        let c = m.atoms[a].coordinates;
                        if c.x < place_next_to.x + dist + dx
                            && c.x > place_next_to.x - dist - dx
                            && c.y < place_next_to.y + dist + dy
                            && c.y > place_next_to.y - dist - dy
                        {
                            no_clash = false;
                            break 'mols;
                        }
                    }
                }
                if no_clash {
                    break;
                }
            }
            if no_clash {
                break;
            }
        }
        place_next_to
    }

    /// `arrangeMultipleMolecules`: places the molecules (non-counterions around the biggest one, then counterions next to
    /// an opposite charge)
    fn arrange_multiple_molecules(&self, m: &mut Model) {
        if self.molecules.len() <= 1 {
            return;
        }
        let mut max_i = 0;
        let mut max_size = m.molecules[self.molecules[0]].atoms.len();
        for (i, &mol) in self.molecules.iter().enumerate() {
            let size = m.molecules[mol].atoms.len();
            if size > max_size {
                max_i = i;
                max_size = size;
            }
        }
        let central = self.molecules[max_i];
        m.molecules[central].is_placed = true;
        for &a in &self.atoms {
            m.atoms[a].general_use_visited = false;
        }
        let center = m.molecule_center(central);
        let grid_d = BONDLENGTH;
        let mut found = true;
        while found {
            found = false;
            for &mol in &self.molecules {
                if m.molecules[mol].has_constrained_fragments {
                    m.molecules[mol].is_placed = true;
                }
                if m.molecules[mol].has_fixed_fragments {
                    m.molecules[mol].is_placed = true;
                }
                if m.molecules[mol].is_placed {
                    continue;
                }
                if m.total_charge(mol) == 0 {
                    let (cmin, cmax) = m.bounding_box(mol);
                    let (counterion_dx, counterion_dy) =
                        ((cmax.x - cmin.x) * 0.5, (cmax.y - cmin.y) * 0.5);
                    let counterion_center = (cmax + cmin) * 0.5;
                    found = true;
                    let place_next_to = self.explore_grid_around(
                        m,
                        center,
                        10,
                        grid_d,
                        counterion_dx,
                        counterion_dy,
                        -1.0,
                    );
                    for a in m.molecules[mol].atoms.clone() {
                        m.atoms[a].coordinates += place_next_to - counterion_center;
                    }
                    m.molecules[mol].is_placed = true;
                }
            }
        }
        found = true;
        while found {
            found = false;
            for &mol in &self.molecules {
                if m.molecules[mol].is_placed {
                    continue;
                }
                let charge = m.total_charge(mol);
                if charge != 0 {
                    let (cmin, cmax) = m.bounding_box(mol);
                    let mut place_next_to = center;
                    let (counterion_dx, counterion_dy) =
                        ((cmax.x - cmin.x) * 0.5, (cmax.y - cmin.y) * 0.5);
                    let counterion_center = (cmax + cmin) * 0.5;
                    found = true;
                    'search: for &mm in &self.molecules {
                        if !m.molecules[mm].is_placed {
                            continue;
                        }
                        for a in m.molecules[mm].atoms.clone() {
                            if m.atoms[a].charge == 0 || m.atoms[a].general_use_visited {
                                continue;
                            }
                            if m.atoms[a].charge * charge < 0 {
                                m.atoms[a].general_use_visited = true;
                                place_next_to = m.atoms[a].coordinates;
                                break 'search;
                            }
                        }
                    }
                    let place_next_to = self.explore_grid_around(
                        m,
                        place_next_to,
                        10,
                        grid_d,
                        counterion_dx,
                        counterion_dy,
                        (f64::from(BONDLENGTH) * 0.8) as f32,
                    );
                    for a in m.molecules[mol].atoms.clone() {
                        m.atoms[a].coordinates += place_next_to - counterion_center;
                    }
                    m.molecules[mol].is_placed = true;
                }
            }
        }
    }

    /// `runGenerateCoordinates`: true if the pose is clean
    pub fn run_generate_coordinates(&mut self, m: &mut Model) -> bool {
        let mut clean_pose = true;
        if self.structure_pass_sanity_check(m) {
            self.find_fragments(m);
            for &mol in &self.molecules {
                m.build_molecule_from_fragments(mol, true);
            }
            clean_pose = self.minimizer.avoid_clashes(m);
            self.best_rotation(m);
            self.maybe_flip(m);
            self.arrange_multiple_molecules(m);
            self.write_stereo_chemistry(m);
        }
        clean_pose
    }
}

#[cfg(test)]
mod tests {
    use super::super::model::{BondStereo, BondStereoInfo};
    use super::super::testutil::build;
    use super::*;
    use std::fmt::Write;
    use std::panic::{AssertUnwindSafe, catch_unwind};

    /// what RDKit's `CoordGen::addCoords` does: no input coordinates, double bond stereo, precision 0.01
    fn run(input: &str) -> String {
        let (mut model, mol) = build(input);
        let atoms = model.molecules[mol].atoms.clone();
        for &a in &atoms {
            model.set_coordinates(a, PointF::new(0.0, 0.0));
        }
        model.assign_bonds_and_neighbors(mol);
        for line in input.lines().filter(|l| l.starts_with("S ")) {
            let f: Vec<&str> = line.split(' ').collect();
            let b: usize = f[1].parse().unwrap();
            let (a1, a2): (usize, usize) = (f[2].parse().unwrap(), f[3].parse().unwrap());
            model.bonds[b].stereo = BondStereoInfo {
                atom1: Some(a1),
                atom2: Some(a2),
                stereo: if f[4] == "c" {
                    BondStereo::Cis
                } else {
                    BondStereo::Trans
                },
            };
            model.bond_set_absolute_stereo_from_stereo_info(b);
        }
        let mut sketcher = Sketcher::new(0.01);
        sketcher.initialize(&mut model, mol);
        let clean = sketcher.run_generate_coordinates(&mut model);
        let mut out = String::new();
        writeln!(out, "clean {}", i32::from(clean)).unwrap();
        for (i, &a) in atoms.iter().enumerate() {
            let c = model.atoms[a].coordinates;
            writeln!(out, "a {i} {:08x} {:08x}", c.x.to_bits(), c.y.to_bits()).unwrap();
        }
        out
    }

    /// Equal up to the platform's `sin`/`cos`/`acos`. The layout stops its minimization on thresholds, so other libms
    /// can end in a slightly different state: off macOS (where the corpus was made) coordinates must agree within 2 units
    /// (a bond is 50).
    fn same(got: &str, want: &str) -> bool {
        let strict = cfg!(target_os = "macos");
        let (gl, wl): (Vec<&str>, Vec<&str>) = (got.lines().collect(), want.lines().collect());
        gl.len() == wl.len()
            && gl.iter().zip(&wl).all(|(g, w)| {
                let (gt, wt): (Vec<&str>, Vec<&str>) = (
                    g.split_whitespace().collect(),
                    w.split_whitespace().collect(),
                );
                gt.len() == wt.len()
                    && gt.iter().zip(&wt).all(|(a, b)| {
                        let hex =
                            |t: &str| t.len() == 8 && t.chars().all(|c| c.is_ascii_hexdigit());
                        if a == b {
                            return true;
                        }
                        if hex(a) && hex(b) {
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

    /// the whole layout must be the one `runGenerateCoordinates` makes (`tests/data/coordgen_sketch.txt`)
    #[test]
    fn layout_matches_cpp() {
        let text = include_str!("../../../../../tests/data/coordgen_sketch.txt");
        let (mut bad, mut n, mut skipped) = (Vec::new(), 0, 0);
        for rec in text.split("$$\n").filter(|r| !r.is_empty()) {
            let end = rec.find("END\n").unwrap() + 4;
            let (input, want) = rec.split_at(end);
            let Ok(got) = catch_unwind(AssertUnwindSafe(|| run(input))) else {
                skipped += 1; // a macrocycle or a ring system for the macrocycle builder
                continue;
            };
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
        assert!(n > 50, "{n} compared, {skipped} skipped");
        assert!(
            bad.is_empty(),
            "{} of {n} differ; first:\n{}",
            bad.len(),
            bad[0]
        );
    }
}
