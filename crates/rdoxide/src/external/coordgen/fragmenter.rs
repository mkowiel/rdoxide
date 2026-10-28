//! `CoordgenFragmenter`: divides a molecule into rigid fragments and finds the main fragment.

use std::collections::{HashMap, VecDeque};

use super::model::{BondId, FragmentId, Model, MolId, RingId};

impl Model {
    /// `splitIntoFragments`
    pub fn split_into_fragments(&mut self, mol: MolId) {
        let mut fragments: Vec<FragmentId> = Vec::new();
        for a in self.molecules[mol].atoms.clone() {
            self.atoms[a].fragment = None;
        }
        // molecule with a single atom and no bonds
        if self.molecules[mol].atoms.len() == 1 {
            let f = self.new_fragment();
            let only = self.molecules[mol].atoms[0];
            self.fragment_add_atom(f, only);
            fragments.push(f);
        }
        for bond in self.molecules[mol].bonds.clone() {
            if self.bond_is_inter_fragment(bond) {
                self.process_inter_fragment_bond(bond, &mut fragments);
            } else {
                self.process_bond_internal_to_fragment(bond, &mut fragments);
            }
        }
        if !fragments.is_empty() {
            self.initialize_fragment_information(fragments, mol);
        }
    }

    fn process_inter_fragment_bond(&mut self, bond: BondId, fragments: &mut Vec<FragmentId>) {
        for atom in [self.bonds[bond].start_atom, self.bonds[bond].end_atom] {
            if self.atoms[atom].fragment.is_none() {
                let f = self.new_fragment();
                self.fragment_add_atom(f, atom);
                fragments.push(f);
            }
        }
    }

    fn process_bond_internal_to_fragment(&mut self, bond: BondId, fragments: &mut Vec<FragmentId>) {
        let (s, e) = (self.bonds[bond].start_atom, self.bonds[bond].end_atom);
        let (fs, fe) = (self.atoms[s].fragment, self.atoms[e].fragment);
        match (fs, fe) {
            (None, None) => {
                let f = self.new_fragment();
                self.fragment_add_atom(f, s);
                self.fragment_add_atom(f, e);
                fragments.push(f);
            }
            (Some(f), None) => self.fragment_add_atom(f, e),
            (None, Some(f)) => self.fragment_add_atom(f, s),
            (Some(f1), Some(f2)) if f1 != f2 => self.join_fragments(f1, f2, fragments),
            _ => {}
        }
    }

    /// `joinFragments`: all atoms of `fragment2` go to `fragment1`, `fragment2` is dropped
    fn join_fragments(
        &mut self,
        fragment1: FragmentId,
        fragment2: FragmentId,
        fragments: &mut Vec<FragmentId>,
    ) {
        for atom in self.fragments[fragment2].atoms().to_vec() {
            self.fragment_add_atom(fragment1, atom);
        }
        fragments.retain(|&f| f != fragment2);
    }

    fn add_bond_information(&mut self, bond: BondId) {
        let (s, e) = (self.bonds[bond].start_atom, self.bonds[bond].end_atom);
        let (fs, fe) = (
            self.atoms[s].fragment.unwrap(),
            self.atoms[e].fragment.unwrap(),
        );
        if fs == fe {
            self.fragment_add_bond(fs, bond);
        } else {
            self.fragment_add_inter_fragment_bond(fs, bond);
            self.fragment_add_inter_fragment_bond(fe, bond);
        }
    }

    fn add_ring_information(&mut self, ring: RingId) {
        let f = self.atoms[self.rings[ring].atoms[0]].fragment.unwrap();
        self.fragment_add_ring(f, ring);
    }

    /// `isChain`: can the fragment be part of a zig-zag chain
    fn fragment_is_chain(&self, f: FragmentId) -> bool {
        let fr = &self.fragments[f];
        if fr.atoms().len() > 3 {
            return false;
        }
        for &a in fr.atoms() {
            if self.atoms[a].bonds.len() > 3 || !self.atoms[a].rings.is_empty() {
                return false;
            }
        }
        fr.bonds().iter().all(|&b| self.bonds[b].bond_order <= 2)
    }

    fn initialize_fragment_information(&mut self, fragments: Vec<FragmentId>, mol: MolId) {
        for b in self.molecules[mol].bonds.clone() {
            self.add_bond_information(b);
        }
        for r in self.molecules[mol].rings.clone() {
            self.add_ring_information(r);
        }
        for &f in &fragments {
            self.fragments[f].is_chain = self.fragment_is_chain(f);
        }
        let mut n_constrained = 0;
        for &f in &fragments {
            let c = self.fragments[f]
                .atoms()
                .iter()
                .any(|&a| self.atoms[a].constrained);
            self.fragments[f].constrained = c;
            n_constrained += usize::from(c);
        }
        self.molecules[mol].has_constrained_fragments = n_constrained > 0;
        let mut n_fixed = 0;
        for &f in &fragments {
            let c = self.fragments[f]
                .atoms()
                .iter()
                .any(|&a| self.atoms[a].fixed);
            self.fragments[f].fixed = c;
            n_fixed += usize::from(c);
        }
        self.molecules[mol].has_fixed_fragments = n_fixed > 0;
        let main = self.find_main_fragment(&fragments);
        self.add_parent_relations_to_fragments(main, &fragments);
        let ordered = self.order_fragments(fragments, main);
        self.molecules[mol].main_fragment = Some(main);
        self.molecules[mol].fragments = ordered.clone();
        for &f in &ordered {
            self.set_flip_info(f);
        }
    }

    fn set_flip_info(&mut self, f: FragmentId) {
        self.fragments[f].constrained_flip = false;
        if self.fragments[f].atoms().len() == 1 && self.fragment_count_constrained_atoms(f) == 1 {
            for &child in &self.fragments[f].children.clone() {
                if self.fragments[child].constrained {
                    self.fragments[f].constrained_flip = true;
                }
            }
        } else {
            self.fragments[f].constrained_flip = self.fragment_count_constrained_atoms(f) > 1;
        }
    }

    /// `getValueOfCheck`: the value of one descriptor and whether this was the last one
    fn fragment_value_of_check(&self, f: FragmentId, check_n: i32) -> (usize, bool) {
        match check_n {
            0 => (self.fragment_count_fixed_atoms(f) as usize, false),
            1 => (self.fragment_count_constrained_atoms(f) as usize, false),
            2 => (self.fragments[f].rings().len(), false),
            3 => (self.fragments[f].atoms().len(), false),
            4 => (self.fragments[f].inter_fragment_bonds.len(), false),
            5 => (self.fragment_count_heavy_atoms(f) as usize, false),
            6 => (self.fragment_total_weight(f) as usize, false),
            7 => (self.fragment_count_double_bonds(f) as usize, false),
            _ => (0, true),
        }
    }

    /// `hasPriority`
    pub fn fragment_has_priority(&self, f1: FragmentId, f2: FragmentId) -> bool {
        let mut check_n = 0;
        loop {
            let (left, no_more1) = self.fragment_value_of_check(f1, check_n);
            let (right, no_more2) = self.fragment_value_of_check(f2, check_n);
            if left > right {
                return true;
            }
            if left < right {
                return false;
            }
            if no_more1 || no_more2 {
                return false;
            }
            check_n += 1;
        }
    }

    fn find_main_fragment(&self, fragments: &[FragmentId]) -> FragmentId {
        let mut main = fragments[0];
        for &f in &fragments[1..] {
            if self.fragment_has_priority(f, main) {
                main = f;
            }
        }
        self.consider_chains(fragments, main)
    }

    fn consider_chains(&self, fragments: &[FragmentId], main: FragmentId) -> FragmentId {
        if fragments
            .iter()
            .any(|&f| self.fragments[f].fixed || self.fragments[f].constrained)
        {
            return main;
        }
        let longest = self.find_longest_chain(fragments);
        if longest.len() >= self.acceptable_chain_length(main) {
            return longest[0];
        }
        main
    }

    fn acceptable_chain_length(&self, main: FragmentId) -> usize {
        match self.fragments[main].rings().len() {
            0 => 1,
            1 => 5,
            2 => 8,
            3 => 10,
            _ => 12,
        }
    }

    /// the fragment on the other side of an inter fragment bond
    fn other_fragment(&self, bond: BondId, fragment: FragmentId) -> FragmentId {
        let b = &self.bonds[bond];
        let fs = self.atoms[b.start_atom].fragment.unwrap();
        if fs != fragment {
            fs
        } else {
            self.atoms[b.end_atom].fragment.unwrap()
        }
    }

    fn find_longest_chain(&self, fragments: &[FragmentId]) -> Vec<FragmentId> {
        let mut longest: Vec<FragmentId> = Vec::new();
        for &fragment in fragments {
            if !self.fragments[fragment].is_chain {
                continue;
            }
            let mut chain_n = 0;
            for &b in &self.fragments[fragment].inter_fragment_bonds {
                if self.fragments[self.other_fragment(b, fragment)].is_chain {
                    chain_n += 1;
                }
            }
            if chain_n > 1 {
                continue; // in the middle of a chain
            }
            let mut q = VecDeque::new();
            let mut parent_map: HashMap<FragmentId, FragmentId> = HashMap::new();
            parent_map.insert(fragment, fragment);
            q.push_back(fragment);
            let mut last_fragment = None;
            while let Some(last) = q.pop_front() {
                last_fragment = Some(last);
                for &b in &self.fragments[last].inter_fragment_bonds {
                    let child = self.other_fragment(b, last);
                    if parent_map.contains_key(&child) || !self.fragments[child].is_chain {
                        continue;
                    }
                    parent_map.insert(child, last);
                    q.push_back(child);
                }
            }
            let mut chain: Vec<FragmentId> = Vec::new();
            if let Some(mut recursive) = last_fragment {
                while recursive != fragment {
                    chain.insert(0, recursive);
                    recursive = parent_map[&recursive];
                }
                chain.insert(0, recursive);
            }
            if chain.len() > longest.len() {
                longest = chain;
            }
        }
        longest
    }

    /// `addParentRelationsToFragments`: breadth first from the main fragment; the parent bond points from the parent
    fn add_parent_relations_to_fragments(&mut self, main: FragmentId, fragments: &[FragmentId]) {
        let mut q = VecDeque::new();
        q.push_back(main);
        while let Some(fragment) = q.pop_front() {
            for bond in self.fragments[fragment].inter_fragment_bonds.clone() {
                let b = &self.bonds[bond];
                let child = if self.atoms[b.start_atom].fragment == Some(fragment) {
                    self.atoms[b.end_atom].fragment.unwrap()
                } else {
                    self.atoms[b.start_atom].fragment.unwrap()
                };
                if Some(child) == self.fragments[fragment].parent {
                    continue;
                }
                self.fragments[fragment].children.push(child);
                self.fragments[child].parent = Some(fragment);
                self.fragments[child].bond_to_parent = Some(bond);
                q.push_back(child);
            }
        }
        for &fragment in fragments {
            if let Some(bond) = self.fragments[fragment].bond_to_parent {
                // start atom is always the parent's, end atom the child's
                if self.atoms[self.bonds[bond].end_atom].fragment != Some(fragment) {
                    let b = &mut self.bonds[bond];
                    std::mem::swap(&mut b.start_atom, &mut b.end_atom);
                    b.is_reversed = !b.is_reversed;
                }
            }
        }
    }

    fn order_fragments(&self, fragments: Vec<FragmentId>, main: FragmentId) -> Vec<FragmentId> {
        let mut q = VecDeque::new();
        let mut ordered = Vec::new();
        q.push_back(main);
        while let Some(f) = q.pop_front() {
            ordered.push(f);
            for &c in &self.fragments[f].children {
                q.push_back(c);
            }
        }
        debug_assert_eq!(fragments.len(), ordered.len());
        ordered
    }
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{build, list};
    use std::fmt::Write;

    /// prints what `tools/coordgen_drv/frag.cpp` prints
    fn run(input: &str) -> String {
        let (mut model, mol) = build(input);
        model.force_update_struct(mol);
        model.split_into_fragments(mol);
        let m = &model.molecules[mol];
        let (atoms, bonds, rings) = (m.atoms.clone(), m.bonds.clone(), m.rings.clone());
        let frags = m.fragments.clone();
        let fi =
            |f: Option<usize>| f.map_or(-1, |f| frags.iter().position(|&x| x == f).unwrap() as i64);
        let mut out = String::new();
        writeln!(
            out,
            "main {} hasFixed {} hasConstrained {} nfrag {}",
            fi(m.main_fragment),
            i32::from(m.has_fixed_fragments),
            i32::from(m.has_constrained_fragments),
            frags.len()
        )
        .unwrap();
        for (i, &f) in frags.iter().enumerate() {
            let fr = &model.fragments[f];
            write!(
                out,
                "f {i} atoms{} bonds{} rings{} inter{} children",
                list(fr.atoms(), &atoms),
                list(fr.bonds(), &bonds),
                list(fr.rings(), &rings),
                list(&fr.inter_fragment_bonds, &bonds)
            )
            .unwrap();
            for &c in &fr.children {
                write!(out, " {}", fi(Some(c))).unwrap();
            }
            writeln!(
                out,
                " parent {} btp {} chain {} fixed {} constrained {} cflip {} weight {} dbl {} heavy {} nconstr {} nfixed {}",
                fi(fr.parent),
                fr.bond_to_parent.map_or(-1, |b| bonds.iter().position(|&x| x == b).unwrap() as i64),
                i32::from(fr.is_chain),
                i32::from(fr.fixed),
                i32::from(fr.constrained),
                i32::from(fr.constrained_flip),
                model.fragment_total_weight(f),
                model.fragment_count_double_bonds(f),
                model.fragment_count_heavy_atoms(f),
                model.fragment_count_constrained_atoms(f),
                model.fragment_count_fixed_atoms(f)
            )
            .unwrap();
        }
        for (i, &b) in bonds.iter().enumerate() {
            let bd = &model.bonds[b];
            writeln!(
                out,
                "b {i} {} {} rev {}",
                atoms.iter().position(|&x| x == bd.start_atom).unwrap(),
                atoms.iter().position(|&x| x == bd.end_atom).unwrap(),
                i32::from(bd.is_reversed)
            )
            .unwrap();
        }
        out
    }

    /// the fragments must be the ones CoordgenFragmenter makes (`tests/data/coordgen_frag.txt`)
    #[test]
    fn fragments_match_cpp() {
        let text = include_str!("../../../../../tests/data/coordgen_frag.txt");
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
