//! The clash search of `CoordgenMinimizer`: scoring of clashes, crossing bonds and degrees of freedom, the search over the
//! fragment degrees of freedom (`flipFragments`, `runSearch`, `growSolutions`), `avoidTerminalClashes`,
//! `checkForClashes` and `avoidClashesOfMolecule`.

use std::collections::{BTreeMap, BTreeSet};

use super::maths::{
    BONDLENGTH, MACROCYCLE, PointF, SKETCHER_EPSILON, intersection_of_segments, points_coincide,
};
use super::minimizer::{Clash, EzConstrain, Interaction, Minimizer};
use super::model::{AtomId, BondId, Model, MolId};

const CLASH_ENERGY_THRESHOLD: f32 = 10.0;
const SAME_SIDE_DPR_PENALTY: f32 = 100.0;
const STANDARD_CROSSING_BOND_PENALTY: f32 = 2500.0;
const TERMINAL_BOND_CROSSING_MULTIPLIER: f32 = 0.5;
const MACROCYCLE_BOND_CROSSING_MULTIPLIER: f32 = 8.0;
const RING_BOND_CROSSING_MULTIPLIER: f32 = 2.0;
const MAXIMUM_NUMBER_OF_SCORED_SOLUTIONS: f32 = 10000.0;
const REJECTED_SOLUTION_SCORE: f32 = 99999999.0;

/// `CoordgenDOFSolutions`: the scored states of a list of degrees of freedom
#[derive(Debug, Clone)]
pub struct DofSolutions {
    molecule: MolId,
    solutions: BTreeMap<Vec<u16>, f32>,
    pub all_dofs: Vec<usize>,
}

impl DofSolutions {
    pub fn new(molecule: MolId, all_dofs: Vec<usize>) -> Self {
        DofSolutions {
            molecule,
            solutions: BTreeMap::new(),
            all_dofs,
        }
    }

    pub fn current_solution(&self, m: &Model) -> Vec<u16> {
        self.all_dofs
            .iter()
            .map(|&d| m.dofs[d].current_state)
            .collect()
    }

    pub fn load_solution(&self, m: &mut Model, solution: &[u16]) {
        debug_assert_eq!(solution.len(), self.all_dofs.len());
        for (i, &d) in self.all_dofs.iter().enumerate() {
            m.dofs[d].current_state = solution[i];
        }
    }

    /// `findBestSolution`
    pub fn find_best_solution(&self) -> (Vec<u16>, f32) {
        let mut it = self.solutions.iter();
        let first = it.next().expect("a scored solution");
        let mut best = (first.0.clone(), *first.1);
        for (s, &score) in self.solutions.iter() {
            if score < best.1 {
                best = (s.clone(), score);
            }
        }
        best
    }
}

impl Minimizer {
    /// `scoreCurrentSolution`
    fn score_current_solution(&mut self, m: &mut Model, sol: &mut DofSolutions) -> f32 {
        let solution = sol.current_solution(m);
        if let Some(&score) = sol.solutions.get(&solution) {
            return score;
        }
        if sol.solutions.len() as f32 > MAXIMUM_NUMBER_OF_SCORED_SOLUTIONS * self.precision {
            return REJECTED_SOLUTION_SCORE;
        }
        m.build_molecule_from_fragments(sol.molecule, false);
        let result = self.score_clashes(m, sol.molecule, true, true);
        sol.solutions.insert(solution, result);
        result
    }

    /// `scoreClashes`
    pub fn score_clashes(
        &mut self,
        m: &mut Model,
        mol: MolId,
        residue_interactions: bool,
        score_proximity: bool,
    ) -> f32 {
        let mut e = 0.0f32;
        for &i in &self.intramolecular_clash_interactions {
            self.interactions[i].score(m, &mut e, true);
        }
        for &i in &self.extra_interactions {
            self.interactions[i].score(m, &mut e, true);
        }
        e += self.score_dofs(m, mol);
        e += self.score_cross_bonds(m, mol, residue_interactions);
        e += self.score_atoms_inside_rings(m);
        if score_proximity {
            e += self.score_proximity_relations_on_opposite_sides(m);
        }
        e
    }

    fn score_dofs(&self, m: &Model, mol: MolId) -> f32 {
        let mut e = 0.0f32;
        for &f in &m.molecules[mol].fragments {
            for &d in &m.fragments[f].dofs {
                e += m.dof_current_penalty(d);
            }
        }
        e
    }

    /// `bondsClash`
    fn bonds_clash(m: &Model, bond: BondId, bond2: BondId) -> bool {
        if bond == bond2 {
            return false;
        }
        let (b1, b2) = (&m.bonds[bond], &m.bonds[bond2]);
        if b1.start_atom == b2.start_atom
            || b1.start_atom == b2.end_atom
            || b1.end_atom == b2.start_atom
            || b1.end_atom == b2.end_atom
        {
            return false;
        }
        let c = |a: AtomId| m.atoms[a].coordinates;
        let (start1, start2, end1, end2) = (
            c(b1.start_atom),
            c(b2.start_atom),
            c(b1.end_atom),
            c(b2.end_atom),
        );
        if start1.x.max(end1.x) < start2.x.min(end2.x)
            || start1.y.max(end1.y) < start2.y.min(end2.y)
            || start1.x.min(end1.x) > start2.x.max(end2.x)
            || start1.y.min(end1.y) > start2.y.max(end2.y)
        {
            return false;
        }
        if points_coincide(start1, start2)
            || points_coincide(start1, end2)
            || points_coincide(end1, start2)
            || points_coincide(end1, end2)
        {
            return true;
        }
        intersection_of_segments(start1, end1, start2, end2).is_some()
    }

    /// `scoreCrossBonds` (without residue interactions)
    fn score_cross_bonds(&self, m: &Model, mol: MolId, _residue_interactions: bool) -> f32 {
        let mut out = 0.0f32;
        let bonds = &m.molecules[mol].bonds;
        if bonds.len() > 2 {
            for b in 0..bonds.len() - 1 {
                let b1 = bonds[b];
                for &b2 in &bonds[b + 1..] {
                    if m.atoms[m.bonds[b2].start_atom].molecule
                        != m.atoms[m.bonds[b1].start_atom].molecule
                    {
                        continue;
                    }
                    if Self::bonds_clash(m, b1, b2) {
                        let mut penalty = STANDARD_CROSSING_BOND_PENALTY
                            * m.bonds[b1].crossing_bond_penalty_multiplier
                            * m.bonds[b2].crossing_bond_penalty_multiplier;
                        if m.bond_is_terminal(b1) || m.bond_is_terminal(b2) {
                            penalty *= TERMINAL_BOND_CROSSING_MULTIPLIER;
                        }
                        if m.bond_is_in_macrocycle(b1) || m.bond_is_in_macrocycle(b2) {
                            penalty *= MACROCYCLE_BOND_CROSSING_MULTIPLIER;
                        }
                        if m.bond_is_in_small_ring(b1) || m.bond_is_in_small_ring(b2) {
                            penalty *= RING_BOND_CROSSING_MULTIPLIER;
                        }
                        out += penalty;
                    }
                }
            }
        }
        out
    }

    /// `scoreAtomsInsideRings`
    fn score_atoms_inside_rings(&self, m: &Model) -> f32 {
        let mut out = 0.0f32;
        let cut_off = BONDLENGTH;
        for &mol in &self.molecules {
            for &r in &m.molecules[mol].rings {
                let ring = &m.rings[r];
                if ring.atoms.len() > MACROCYCLE || ring.atoms.len() < 3 {
                    continue;
                }
                let c = m.ring_find_center(r);
                for &a in &m.molecules[mol].atoms {
                    if m.atoms[a].fragment == m.atoms[ring.atoms[0]].fragment {
                        continue;
                    }
                    let d = c - m.atoms[a].coordinates;
                    if d.x > cut_off || d.y > cut_off || d.x < -cut_off || d.y < -cut_off {
                        continue;
                    }
                    if d.square_length() > cut_off * cut_off {
                        continue;
                    }
                    let dist = d.length();
                    if dist < cut_off {
                        out += 50.0 + 100.0 * (1.0 - (dist / cut_off));
                    }
                }
            }
        }
        out
    }

    /// `scoreProximityRelationsOnOppositeSides`: molecules have no proximity relations without residues
    fn score_proximity_relations_on_opposite_sides(&self, m: &Model) -> f32 {
        let mut out = 0.0f32;
        for &mol in &self.molecules {
            if m.molecules[mol].atoms.len() < 2 {
                continue;
            }
            assert!(
                m.molecules[mol].proximity_relations.is_empty(),
                "proximity relations are not ported"
            );
        }
        let _ = SAME_SIDE_DPR_PENALTY;
        out += 0.0;
        out
    }

    /// `checkForClashes`: moves a newly placed atom a little if that relieves clashes with the bonds of its fragment
    pub fn check_for_clashes(m: &mut Model, a: AtomId) {
        if m.atoms[a].fixed {
            return;
        }
        let old = m.atoms[a].coordinates;
        let d = BONDLENGTH * 0.25;
        #[allow(clippy::approx_constant)] // upstream uses this literal
        let e = BONDLENGTH * 0.25 * 0.7071;
        let coords_vect = [
            old,
            old + PointF::new(d, 0.0),
            old + PointF::new(-d, 0.0),
            old + PointF::new(0.0, d),
            old + PointF::new(0.0, -d),
            old + PointF::new(e, -e),
            old + PointF::new(-e, -e),
            old + PointF::new(-e, e),
            old + PointF::new(e, e),
        ];
        let mut best_e = 999999.0f32;
        let mut best_i = 0;
        let fragment = m.atoms[a].fragment.expect("fragment");
        for (i, &cv) in coords_vect.iter().enumerate() {
            m.atoms[a].coordinates = cv;
            let mut clash = Clash::new(a, a, a);
            clash.rest_v = 300.0;
            let mut clash_e = 0.0f32;
            for &b in m.fragments[fragment].bonds() {
                let bd = &m.bonds[b];
                if !m.atoms[bd.start_atom].coordinates_set || !m.atoms[bd.end_atom].coordinates_set
                {
                    continue;
                }
                if bd.start_atom == a || bd.end_atom == a {
                    continue;
                }
                clash.atom1 = bd.start_atom;
                clash.atom2 = a;
                clash.atom3 = bd.end_atom;
                clash.energy(m, &mut clash_e);
            }
            for &b in &m.atoms[a].bonds {
                let bd = &m.bonds[b];
                for &atom in m.fragments[fragment].atoms() {
                    if atom == a {
                        continue;
                    }
                    if !m.atoms[bd.start_atom].coordinates_set
                        || !m.atoms[bd.end_atom].coordinates_set
                    {
                        continue;
                    }
                    if bd.start_atom == atom || bd.end_atom == atom {
                        continue;
                    }
                    clash.atom1 = bd.start_atom;
                    clash.atom2 = atom;
                    clash.atom3 = bd.end_atom;
                    clash.energy(m, &mut clash_e);
                }
            }
            if clash_e < SKETCHER_EPSILON {
                return;
            }
            if i == 0 {
                best_e = clash_e;
            }
            if clash_e < best_e {
                best_e = clash_e;
                best_i = i;
            }
        }
        m.set_coordinates(a, coords_vect[best_i]);
    }

    /// `avoidTerminalClashes`
    fn avoid_terminal_clashes(&mut self, m: &mut Model, mol: MolId, clash_e: &mut f32) {
        if *clash_e < 0.1 {
            return;
        }
        for bond in m.molecules[mol].bonds.clone() {
            if !m.bond_is_terminal(bond) {
                continue;
            }
            let (mut terminal, mut root) = (m.bonds[bond].end_atom, m.bonds[bond].start_atom);
            if m.atoms[terminal].bonds.len() != 1 {
                terminal = m.bonds[bond].start_atom;
                root = m.bonds[bond].end_atom;
            }
            if m.atoms[terminal].fixed {
                continue;
            }
            for bond2 in m.molecules[mol].bonds.clone() {
                if Self::bonds_clash(m, bond, bond2) {
                    let rc = m.atoms[root].coordinates;
                    let tc = m.atoms[terminal].coordinates;
                    m.set_coordinates(terminal, rc + (tc - rc) * 0.1);
                }
            }
        }
        *clash_e = self.score_clashes(m, mol, false, true);
    }

    /// `growSolutions`: false if a solution was rejected
    fn grow_solutions(
        &mut self,
        m: &mut Model,
        sol: &mut DofSolutions,
        all_scored: &mut BTreeSet<Vec<u16>>,
        current_tier: &mut i32,
        growing: &mut BTreeMap<Vec<u16>, f32>,
        best_score: &mut f32,
    ) -> bool {
        let old_growing = growing.clone();
        let best_score_for_run = *best_score;
        let mut best_solutions: Vec<(f32, Vec<u16>)> =
            growing.iter().map(|(s, &sc)| (sc, s.clone())).collect();
        best_solutions.sort_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.1.cmp(&b.1))
        });
        growing.clear();
        let mut max_n = (6.0 * self.precision) as i32;
        if max_n < 1 {
            max_n = 1;
        }
        let mut n = 0;
        for solution in &best_solutions {
            if n > max_n {
                break;
            }
            for dof in sol.all_dofs.clone() {
                if m.dofs[dof].kind.tier() > *current_tier {
                    continue;
                }
                sol.load_solution(m, &solution.1);
                for _ in 1..m.dof_number_of_states(dof) {
                    m.dof_change_state(dof);
                    let new_solution = sol.current_solution(m);
                    if !all_scored.contains(&new_solution) {
                        let score = self.score_current_solution(m, sol);
                        if score == REJECTED_SOLUTION_SCORE {
                            return false;
                        }
                        all_scored.insert(new_solution.clone());
                        if score < *best_score {
                            *best_score = score;
                        }
                        if score < best_score_for_run && score < REJECTED_SOLUTION_SCORE {
                            growing.insert(new_solution, score);
                        }
                    }
                }
            }
            n += 1;
        }
        if growing.is_empty() && *current_tier < 5 {
            *current_tier += 3;
            *growing = old_growing;
        }
        true
    }

    /// `runSearch`
    fn run_search(&mut self, m: &mut Model, mut tier: i32, sol: &mut DofSolutions) -> bool {
        let mut growing: BTreeMap<Vec<u16>, f32> = BTreeMap::new();
        let mut all_scored: BTreeSet<Vec<u16>> = BTreeSet::new();
        let mut best_score = self.score_current_solution(m, sol);
        growing.insert(sol.current_solution(m), best_score);
        let mut i = 0;
        loop {
            i += 1;
            let valid = self.grow_solutions(
                m,
                sol,
                &mut all_scored,
                &mut tier,
                &mut growing,
                &mut best_score,
            );
            if !(valid && !growing.is_empty() && i < 100) {
                break;
            }
        }
        let best = sol.find_best_solution();
        sol.load_solution(m, &best.0);
        best.1 < CLASH_ENERGY_THRESHOLD
    }

    /// `flipFragments`
    fn flip_fragments(&mut self, m: &mut Model, mol: MolId, clash_e: &mut f32) -> bool {
        if *clash_e < CLASH_ENERGY_THRESHOLD {
            return true;
        }
        let mut dofs = Vec::new();
        let mut fragments = m.molecules[mol].fragments.clone();
        fragments.reverse();
        for f in fragments {
            if !m.fragments[f].fixed {
                for &d in &m.fragments[f].dofs {
                    if m.dof_number_of_states(d) > 1 {
                        dofs.push(d);
                    }
                }
            }
        }
        let mut solutions = DofSolutions::new(mol, dofs);
        let clean_pose = self.run_search(m, 0, &mut solutions);
        m.build_molecule_from_fragments(mol, false);
        clean_pose
    }

    /// `addPeptideBondInversionConstraintsOfMolecule`: keeps the backbone of peptides in a straight trans line
    fn add_peptide_bond_inversion_constraints_of_molecule(&mut self, m: &Model, mol: MolId) {
        let atoms = &m.molecules[mol].atoms;
        let cheto_cs: BTreeSet<AtomId> = atoms
            .iter()
            .copied()
            .filter(|&a| {
                m.atoms[a].atomic_number == 6
                    && m.atoms[a].neighbors.iter().any(|&n| {
                        m.atoms[n].atomic_number == 8
                            && m.bond_to(a, n).is_some_and(|b| m.bonds[b].bond_order == 2)
                    })
            })
            .collect();
        if cheto_cs.len() < 2 {
            return;
        }
        let amino_ns: BTreeSet<AtomId> = atoms
            .iter()
            .copied()
            .filter(|&a| m.atoms[a].atomic_number == 7)
            .collect();
        if amino_ns.len() < 2 {
            return;
        }
        let alpha_cs: BTreeSet<AtomId> = atoms
            .iter()
            .copied()
            .filter(|&a| {
                m.atoms[a].atomic_number == 6
                    && !cheto_cs.contains(&a)
                    && m.atoms[a].neighbors.iter().any(|n| cheto_cs.contains(n))
                    && m.atoms[a].neighbors.iter().any(|n| amino_ns.contains(n))
            })
            .collect();
        if alpha_cs.len() < 2 {
            return;
        }
        let mut groups: Vec<[AtomId; 4]> = Vec::new();
        let mut four = |s1: &BTreeSet<AtomId>,
                        s2: &BTreeSet<AtomId>,
                        s3: &BTreeSet<AtomId>,
                        s4: &BTreeSet<AtomId>| {
            for &a1 in s1 {
                for &a2 in &m.atoms[a1].neighbors {
                    if !s2.contains(&a2) {
                        continue;
                    }
                    for &a3 in &m.atoms[a2].neighbors {
                        if !s3.contains(&a3) {
                            continue;
                        }
                        for &a4 in &m.atoms[a3].neighbors {
                            if s4.contains(&a4) {
                                groups.push([a1, a2, a3, a4]);
                            }
                        }
                    }
                }
            }
        };
        four(&cheto_cs, &amino_ns, &alpha_cs, &cheto_cs);
        four(&amino_ns, &alpha_cs, &cheto_cs, &amino_ns);
        four(&alpha_cs, &cheto_cs, &amino_ns, &alpha_cs);
        for g in groups {
            self.extra_interactions.push(self.interactions.len());
            self.interactions
                .push(Interaction::EzConstrain(EzConstrain {
                    atom1: g[0],
                    atom2: g[1],
                    atom3: g[2],
                    atom4: g[3],
                    is_z: false,
                    force_movement: false,
                }));
        }
    }

    /// `avoidClashesOfMolecule`: true if the pose is clean
    pub fn avoid_clashes_of_molecule(&mut self, m: &mut Model, mol: MolId) -> bool {
        self.clear_interactions();
        self.add_clash_interactions_of_molecule(m, mol, false);
        self.add_peptide_bond_inversion_constraints_of_molecule(m, mol);
        let mut clash_e = self.score_clashes(m, mol, true, true);
        let clean_pose = self.flip_fragments(m, mol, &mut clash_e);
        if !clean_pose {
            self.avoid_terminal_clashes(m, mol, &mut clash_e);
            m.molecules[mol].require_minimization = true;
        }
        if m.molecules[mol].require_minimization {
            self.minimize_molecule(m, mol);
        }
        clean_pose
    }

    /// `avoidClashes`
    pub fn avoid_clashes(&mut self, m: &mut Model) -> bool {
        let mut all_clean = true;
        for mol in self.molecules.clone() {
            let clean = self.avoid_clashes_of_molecule(m, mol);
            all_clean = all_clean && clean;
        }
        all_clean
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
        let mut min = Minimizer {
            atoms: atoms.clone(),
            molecules: vec![mol],
            ..Minimizer::default()
        };
        model.build_molecule_from_fragments(mol, true);
        let clean = min.avoid_clashes(&mut model);
        let e = min.score_clashes(&mut model, mol, false, true);
        let mut out = String::new();
        writeln!(out, "clean {} energy {:08x}", i32::from(clean), e.to_bits()).unwrap();
        for (i, &a) in atoms.iter().enumerate() {
            let c = model.atoms[a].coordinates;
            writeln!(out, "a {i} {:08x} {:08x}", c.x.to_bits(), c.y.to_bits()).unwrap();
        }
        out
    }

    /// Equal up to the platform's `sin`/`cos`/`acos`. The minimization stops on thresholds, so other libms can end in a
    /// slightly different state: off macOS (where the corpus was made) coordinates must agree within 2 units (a bond is
    /// 50) and the energies are not compared.
    fn same(got: &str, want: &str) -> bool {
        let strict = cfg!(target_os = "macos");
        let (gl, wl): (Vec<&str>, Vec<&str>) = (got.lines().collect(), want.lines().collect());
        gl.len() == wl.len()
            && gl.iter().zip(&wl).all(|(g, w)| {
                let (gt, wt): (Vec<&str>, Vec<&str>) = (
                    g.split_whitespace().collect(),
                    w.split_whitespace().collect(),
                );
                let energy_line = w.starts_with("interactions") || w.starts_with("clean");
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

    /// the clash avoidance must give the coordinates of `CoordgenMinimizer::avoidClashes` (`tests/data/coordgen_clash.txt`)
    #[test]
    fn avoid_clashes_matches_cpp() {
        let text = include_str!("../../../../../tests/data/coordgen_clash.txt");
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
