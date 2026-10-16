//! The restructuring of the group tree that makes a Marvin molecule
//! palatable to RDKit (`prepSgroupsForRDKit` and what it triggers): multiple
//! groups are expanded, contracted superatoms are unfolded and multicenter
//! groups are dropped.

use super::defs::{AtomId, InAtomSet, MarvinMol, MarvinNode, NodeId, SgroupKind};
use crate::{RDKitError, Result};

fn file_parse(msg: impl Into<String>) -> RDKitError {
    RDKitError::FileParse(msg.into())
}

fn new_atom_missing() -> RDKitError {
    file_parse("unexpected error - new atom not found in parent")
}

fn new_bond_missing() -> RDKitError {
    file_parse("unexpected error - new bond not found in parent")
}

/// Inserts `item` after `after` in `list`, or at the end if `after` is not in
/// the list.
fn insert_after(list: &mut Vec<usize>, after: usize, item: usize) {
    match list.iter().position(|&x| x == after) {
        Some(p) => list.insert(p + 1, item),
        None => list.push(item),
    }
}

impl MarvinMol {
    /// `moveSgroup`: the group is looked up in its parent by `id` (so with
    /// duplicate ids another group may be the one that moves), appended to
    /// the new parent, and `moved` gets the new parent.
    pub fn move_sgroup(
        &mut self,
        moved: NodeId,
        new_parent: NodeId,
        remove_from_old_parent: bool,
    ) -> Result<()> {
        let parent = self.nodes[moved]
            .parent
            .ok_or_else(|| RDKitError::precondition("bad sgroupToMove parent pointer"))?;
        let id = &self.nodes[moved].id;
        let pos = self.nodes[parent]
            .children
            .iter()
            .position(|&c| self.nodes[c].id == *id)
            .ok_or_else(|| file_parse("Unexpected error - sgroup not found in its parent"))?;
        let found = self.nodes[parent].children[pos];
        self.nodes[new_parent].children.push(found);
        if remove_from_old_parent {
            self.nodes[parent].children.remove(pos);
        }
        self.nodes[moved].parent = Some(new_parent);
        Ok(())
    }

    /// `promoteChild`: moves the group to its grandparent.
    pub fn promote_child(&mut self, node: NodeId) -> Result<()> {
        let parent = self.nodes[node]
            .parent
            .ok_or_else(|| RDKitError::precondition("bad molToMove pointer"))?;
        let grandparent = self.nodes[parent]
            .parent
            .ok_or_else(|| RDKitError::precondition("bad molToMove parent pointer"))?;
        self.move_sgroup(node, grandparent, true)
    }

    /// `isSgroupInSetOfAtoms`
    pub fn is_sgroup_in_set_of_atoms(&self, node: NodeId, set: &[AtomId]) -> Result<InAtomSet> {
        let in_set = |a: AtomId| set.iter().any(|&s| self.atoms[s].id == self.atoms[a].id);
        if matches!(self.nodes[node].kind, SgroupKind::Superatom { .. }) {
            // in the set if the condensed atom in the parent is
            let parent = self.nodes[node].parent.expect("a group has a parent");
            let dummy = self.nodes[parent]
                .atoms
                .iter()
                .copied()
                .find(|&a| self.atoms[a].sgroup_ref == self.nodes[node].id)
                .ok_or_else(|| file_parse("No contracted atom found for a superatomSgroup"))?;
            return Ok(if in_set(dummy) {
                InAtomSet::In
            } else {
                InAtomSet::NotIn
            });
        }
        let mut is_in = false;
        let mut is_not_in = false;
        for &a in &self.nodes[node].atoms {
            if in_set(a) {
                is_in = true;
            } else {
                is_not_in = true;
            }
        }
        Ok(match (is_in, is_not_in) {
            (true, true) => InAtomSet::BothInAndNotIn,
            (true, false) => InAtomSet::In,
            _ => InAtomSet::NotIn,
        })
    }

    /// `prepSgroupsForRDKit`: fixes the hierarchy of the groups below `node`
    /// (a group that does not belong under its parent is promoted) and
    /// processes the groups that change the atoms.
    pub fn prep_sgroups_for_rdkit(&mut self, node: NodeId) -> Result<()> {
        // fixing one child can add groups at this level
        let mut done: Vec<String> = Vec::new();
        while let Some(child) = self.nodes[node]
            .children
            .iter()
            .copied()
            .find(|&c| !done.contains(&self.nodes[c].mol_id))
        {
            done.push(self.nodes[child].mol_id.clone());
            self.prep_sgroups_for_rdkit(child)?;
        }
        let Some(parent) = self.nodes[node].parent else {
            return Ok(());
        };
        // a child of the root cannot move up
        if self.nodes[parent].parent.is_some() {
            if self.nodes[parent].is_passive_role_for_expansion() {
                return self.promote_child(node);
            }
            let parent_atoms = self.nodes[parent].atoms.clone();
            match self.is_sgroup_in_set_of_atoms(node, &parent_atoms)? {
                InAtomSet::NotIn => return self.promote_child(node),
                InAtomSet::In => {}
                InAtomSet::BothInAndNotIn => {
                    return Err(file_parse(
                        "Child sGroup has atoms both in the parent and NOT in the parent",
                    ));
                }
            }
        }
        match self.nodes[node].kind {
            SgroupKind::Superatom { .. } => self.convert_from_one_superatom(node),
            SgroupKind::Multicenter { .. } => self.process_one_multicenter(node),
            SgroupKind::Multiple { .. } => self.expand_one_multiple(node),
            _ => Ok(()),
        }
    }

    /// The groups from `from` up to and including `to`.
    fn chain(&self, from: NodeId, to: NodeId) -> Result<Vec<NodeId>> {
        let mut out = Vec::new();
        let mut cur = from;
        loop {
            out.push(cur);
            if cur == to {
                return Ok(out);
            }
            cur = self.nodes[cur]
                .parent
                .ok_or_else(|| RDKitError::invariant("group chain ended early"))?;
        }
    }

    /// Removes `atom` from the atoms (and for multiple groups the parent
    /// atoms) of the groups on the chain.
    fn remove_atom_on_chain(&mut self, chain: &[NodeId], atom: AtomId) -> Result<()> {
        for &g in chain {
            let node = &mut self.nodes[g];
            let p = node
                .atoms
                .iter()
                .position(|&a| a == atom)
                .ok_or_else(|| RDKitError::invariant("atom not found in a parent group"))?;
            node.atoms.remove(p);
            if let SgroupKind::Multiple { parent_atoms, .. } = &mut node.kind {
                let p = parent_atoms
                    .iter()
                    .position(|&a| a == atom)
                    .ok_or_else(|| RDKitError::invariant("atom not found in parent atoms"))?;
                parent_atoms.remove(p);
            }
        }
        Ok(())
    }

    /// `processOneMulticenterSgroup`: the center atom and its bonds go.
    fn process_one_multicenter(&mut self, node: NodeId) -> Result<()> {
        let SgroupKind::Multicenter { center } = self.nodes[node].kind else {
            unreachable!()
        };
        let parent = self.nodes[node].parent.expect("a group has a parent");
        let actual = self.actual_parent(node);
        let center_id = self.atoms[center].id.clone();
        let bonds = &self.bonds;
        self.nodes[actual].bonds.retain(|&b| {
            bonds[b].atom_refs2[0] != center_id && bonds[b].atom_refs2[1] != center_id
        });
        // another multicenter group may have removed the atom already
        if self.nodes[actual].atoms.contains(&center) {
            let chain = self.chain(parent, actual)?;
            self.remove_atom_on_chain(&chain, center)?;
        }
        let id = self.nodes[node].id.clone();
        for g in self.chain(parent, actual)? {
            let p = self.nodes[g]
                .children
                .iter()
                .position(|&c| self.nodes[c].id == id)
                .ok_or_else(|| RDKitError::invariant("sgroup not found in its parent"))?;
            self.nodes[g].children.remove(p);
        }
        Ok(())
    }

    /// `copyMol`: copies a group below a multiple group that is being
    /// expanded; the atoms and bonds of the copy must already exist.
    fn copy_node(&mut self, node: NodeId, suffix: &str) -> Result<NodeId> {
        let parent = self.nodes[node].parent.expect("a group has a parent");
        let actual = self.actual_parent(node);
        // the group whose atoms and bonds the copy refers to
        let owner = match self.nodes[node].kind {
            SgroupKind::SruCoMod { .. } | SgroupKind::SuperatomExpanded { .. } => actual,
            _ => parent,
        };
        let find_atom = |m: &Self, id: &str| -> Result<AtomId> {
            m.nodes[owner]
                .atoms
                .iter()
                .copied()
                .find(|&a| m.atoms[a].id == id)
                .ok_or_else(new_atom_missing)
        };
        let find_bond = |m: &Self, id: &str| -> Result<usize> {
            m.nodes[owner]
                .bonds
                .iter()
                .copied()
                .find(|&b| m.bonds[b].id == id)
                .ok_or_else(new_bond_missing)
        };
        let src = self.nodes[node].clone();
        let copy_atoms = |m: &Self, list: &[AtomId]| -> Result<Vec<AtomId>> {
            list.iter()
                .map(|&a| find_atom(m, &format!("{}{suffix}", m.atoms[a].id)))
                .collect()
        };
        let copy_bonds = |m: &Self, list: &[usize]| -> Result<Vec<usize>> {
            list.iter()
                .map(|&b| find_bond(m, &format!("{}{suffix}", m.bonds[b].id)))
                .collect()
        };
        let mut out = MarvinNode {
            mol_id: format!("{}{suffix}", src.mol_id),
            id: format!("{}{suffix}", src.id),
            coordinate_precision: 6,
            atoms: Vec::new(),
            bonds: Vec::new(),
            children: Vec::new(),
            parent: Some(parent),
            kind: src.kind.clone(),
        };
        match &src.kind {
            SgroupKind::Mol | SgroupKind::Superatom { .. } => {
                return Err(RDKitError::precondition("cannot copy this kind of group"));
            }
            SgroupKind::SruCoMod { .. } => {
                // the copy has no id of its own
                out.id.clear();
                out.atoms = copy_atoms(self, &src.atoms)?;
                out.bonds = copy_bonds(self, &src.bonds)?;
            }
            SgroupKind::SuperatomExpanded { .. } => {
                out.atoms = copy_atoms(self, &src.atoms)?;
                out.bonds = copy_bonds(self, &src.bonds)?;
            }
            SgroupKind::Multiple {
                title,
                is_expanded,
                parent_atoms,
                ..
            } => {
                out.atoms = copy_atoms(self, &src.atoms)?;
                let new_parent_atoms = copy_atoms(self, parent_atoms)?;
                out.bonds = copy_bonds(self, &src.bonds)?;
                out.kind = SgroupKind::Multiple {
                    title: title.clone(),
                    is_expanded: *is_expanded,
                    parent_atoms: new_parent_atoms,
                    bonds_to_atoms_not_in_expanded_group: Vec::new(),
                };
            }
            SgroupKind::Multicenter { center } => {
                // the copy keeps the id of the original
                out.id = src.id.clone();
                out.atoms = copy_atoms(self, &src.atoms)?;
                let c = find_atom(self, &format!("{}{suffix}", self.atoms[*center].id))?;
                out.kind = SgroupKind::Multicenter { center: c };
            }
            SgroupKind::Data { .. } | SgroupKind::Generic { .. } | SgroupKind::Monomer { .. } => {
                out.atoms = copy_atoms(self, &src.atoms)?;
            }
        }
        self.nodes.push(out);
        let new = self.nodes.len() - 1;
        self.nodes[parent].children.push(new);
        Ok(new)
    }

    /// `expandOneMultipleSgroup`: the atoms of the group are replicated
    /// (MOL files list all of them) and the replicates are chained.
    fn expand_one_multiple(&mut self, node: NodeId) -> Result<()> {
        let SgroupKind::Multiple {
            title,
            is_expanded,
            parent_atoms,
            ..
        } = self.nodes[node].kind.clone()
        else {
            unreachable!()
        };
        if is_expanded {
            return Ok(());
        }
        let actual = self.actual_parent(node);
        let group_atoms = self.nodes[node].atoms.clone();
        let in_group = |m: &Self, id: &str| group_atoms.iter().any(|&a| m.atoms[a].id == id);

        let mut original_connectors: Vec<String> = Vec::new();
        let mut outside_connector = String::new();
        let mut bonds_in_group = Vec::new();
        let mut to_outside = Vec::new();
        for &b in &self.nodes[actual].bonds {
            let bond = &self.bonds[b];
            let first = in_group(self, &bond.atom_refs2[0]);
            let second = in_group(self, &bond.atom_refs2[1]);
            if first && second {
                bonds_in_group.push(b);
            } else if first {
                original_connectors.push(bond.atom_refs2[0].clone());
                to_outside.push(b);
                outside_connector = bond.atom_refs2[1].clone();
            } else if second {
                original_connectors.push(bond.atom_refs2[1].clone());
                to_outside.push(b);
                outside_connector = bond.atom_refs2[0].clone();
            }
        }
        if to_outside.len() != 2 && !to_outside.is_empty() {
            return Err(file_parse(
                "Error - there must be zero or two bonds from the group to atoms not in the group",
            ));
        }

        // the replicates go after the last atom of the group
        let highest = group_atoms
            .iter()
            .filter_map(|&a| {
                self.nodes[actual]
                    .atoms
                    .iter()
                    .position(|&x| self.atoms[x].id == self.atoms[a].id)
            })
            .max()
            .unwrap_or(0);
        let mut last_atom = *self.nodes[actual]
            .atoms
            .get(highest)
            .ok_or_else(|| file_parse("Multiple group without atoms"))?;

        let copy_count: i32 = title
            .trim()
            .parse()
            .map_err(|_| file_parse("invalid argument"))?;
        let parent = self.nodes[node].parent.expect("a group has a parent");
        let chain = self.chain(parent, actual)?;
        let mut connector1 = if to_outside.len() == 2 {
            original_connectors[1].clone()
        } else {
            String::new()
        };

        for copy_index in 2..=copy_count {
            let suffix = format!("_{}:{copy_index}", self.nodes[node].id);
            for &pa in &parent_atoms {
                let src = &self.atoms[pa];
                let mut copy = src.copy_with_id(format!("{}{suffix}", src.id));
                if !src.sgroup_ref.is_empty() && !copy.sgroup_ref_is_superatom {
                    copy.sgroup_ref = format!("{}{suffix}", src.sgroup_ref);
                }
                self.atoms.push(copy);
                let copy_id = self.atoms.len() - 1;
                for &g in &chain {
                    insert_after(&mut self.nodes[g].atoms, last_atom, copy_id);
                    if let SgroupKind::Multiple { parent_atoms, .. } = &mut self.nodes[g].kind {
                        insert_after(parent_atoms, last_atom, copy_id);
                    }
                }
                last_atom = copy_id;
                self.nodes[node].atoms.push(copy_id);
            }
            for &b in &bonds_in_group {
                let src = &self.bonds[b];
                let copy = src.copy_with(
                    format!("{}{suffix}", src.id),
                    format!("{}{suffix}", src.atom_refs2[0]),
                    format!("{}{suffix}", src.atom_refs2[1]),
                );
                self.bonds.push(copy);
                let id = self.bonds.len() - 1;
                self.nodes[actual].bonds.push(id);
            }
            let n_children = self.nodes[node].children.len();
            for i in (0..n_children).rev() {
                let Some(&child) = self.nodes[node].children.get(i) else {
                    break;
                };
                let copied = self.copy_node(child, &suffix)?;
                self.move_sgroup(copied, parent, true)?;
            }
            if to_outside.len() == 2 {
                let connector0 = format!("{}{suffix}", original_connectors[0]);
                let src = &self.bonds[to_outside[0]];
                let bond = src.copy_with(
                    format!("{}{suffix}", src.id),
                    connector0,
                    connector1.clone(),
                );
                self.bonds.push(bond);
                let id = self.bonds.len() - 1;
                self.nodes[actual].bonds.push(id);
                connector1 = format!("{}{suffix}", original_connectors[1]);
            }
        }

        // the bond to the second outside atom comes from the last replicate
        if to_outside.len() == 2 {
            self.bonds[to_outside[1]].atom_refs2 = [outside_connector, connector1];
        }
        for i in (0..self.nodes[node].children.len()).rev() {
            if let Some(&child) = self.nodes[node].children.get(i) {
                self.promote_child(child)?;
            }
        }
        if let SgroupKind::Multiple {
            is_expanded,
            bonds_to_atoms_not_in_expanded_group,
            ..
        } = &mut self.nodes[node].kind
        {
            *is_expanded = true;
            *bonds_to_atoms_not_in_expanded_group = to_outside;
        }
        Ok(())
    }

    /// `convertFromOneSuperAtom`: the atoms and bonds of a contracted
    /// superatom replace its dummy atom in the parent; a
    /// `SuperatomExpanded` group takes its place in the tree.
    fn convert_from_one_superatom(&mut self, node: NodeId) -> Result<()> {
        let SgroupKind::Superatom {
            title,
            attachment_points,
        } = self.nodes[node].kind.clone()
        else {
            unreachable!()
        };
        let parent = self.nodes[node].parent.expect("a group has a parent");
        let mut expanded = MarvinNode {
            mol_id: self.nodes[node].mol_id.clone(),
            id: self.nodes[node].id.clone(),
            coordinate_precision: 6,
            atoms: Vec::new(),
            bonds: Vec::new(),
            children: Vec::new(),
            parent: Some(parent),
            kind: SgroupKind::SuperatomExpanded { title },
        };
        self.nodes.push(std::mem::replace(
            &mut expanded,
            MarvinNode {
                mol_id: String::new(),
                id: String::new(),
                coordinate_precision: 6,
                atoms: Vec::new(),
                bonds: Vec::new(),
                children: Vec::new(),
                parent: None,
                kind: SgroupKind::Mol,
            },
        ));
        let expanded_id = self.nodes.len() - 1;
        self.nodes[parent].children.push(expanded_id);

        let actual = self.actual_parent(node);
        let sg_id = self.nodes[node].id.clone();
        let coord_exist = self.has_2d_coords(node);
        let dummy = self.nodes[actual]
            .atoms
            .iter()
            .copied()
            .find(|&a| self.atoms[a].sgroup_ref == sg_id)
            .ok_or_else(|| file_parse("No contracted atom found for a superatomSgroup"))?;

        // siblings that contain the dummy atom get the atoms of the superatom
        let node_atoms = self.nodes[node].atoms.clone();
        for sibling in self.nodes[parent].children.clone() {
            if sibling == node {
                continue;
            }
            let atoms = &mut self.nodes[sibling].atoms;
            if let Some(p) = atoms.iter().position(|&a| a == dummy) {
                atoms.remove(p);
                atoms.extend_from_slice(&node_atoms);
            }
        }

        // the bonds that lose the dummy atom
        let dummy_id = self.atoms[dummy].id.clone();
        let mut orphaned = 0usize;
        let (mut cx, mut cy) = (0.0f64, 0.0f64);
        for &b in &self.nodes[actual].bonds {
            let bond = &self.bonds[b];
            let attached = if bond.atom_refs2[0] == dummy_id {
                &bond.atom_refs2[1]
            } else if bond.atom_refs2[1] == dummy_id {
                &bond.atom_refs2[0]
            } else {
                continue;
            };
            let _ = attached;
            let point = attachment_points
                .iter()
                .find(|p| p.bond == bond.id)
                .ok_or_else(|| {
                    file_parse(
                        "No attachment point found for bond to the condensed atom in a superatomSgroup",
                    )
                })?;
            orphaned += 1;
            let idx = self.atom_index(node, &point.atom);
            if idx < 0 {
                return Err(file_parse(
                    "Attachment atom not found in the superatomsGroup",
                ));
            }
            if coord_exist {
                let a = &self.atoms[self.nodes[node].atoms[idx as usize]];
                cx += a.x2;
                cy += a.y2;
            }
        }
        if coord_exist {
            let (ox, oy) = if orphaned > 0 {
                (
                    self.atoms[dummy].x2 - cx / orphaned as f64,
                    self.atoms[dummy].y2 - cy / orphaned as f64,
                )
            } else {
                let n = node_atoms.len() as f64;
                let (sx, sy) = node_atoms.iter().fold((0.0, 0.0), |(sx, sy), &a| {
                    (sx + self.atoms[a].x2, sy + self.atoms[a].y2)
                });
                (self.atoms[dummy].x2 - sx / n, self.atoms[dummy].y2 - sy / n)
            };
            for &a in &node_atoms {
                self.atoms[a].x2 += ox;
                self.atoms[a].y2 += oy;
            }
        }

        // remove the dummy atom
        let chain = self.chain(parent, actual)?;
        self.remove_atom_on_chain(&chain, dummy)?;

        // add the atoms and bonds of the superatom to the parents
        for &a in &node_atoms {
            for &g in &chain {
                self.nodes[g].atoms.push(a);
                if let SgroupKind::Multiple { parent_atoms, .. } = &mut self.nodes[g].kind {
                    parent_atoms.push(a);
                }
            }
            self.nodes[expanded_id].atoms.push(a);
            self.atoms[a].sgroup_attachment_point.clear();
        }
        let node_bonds = self.nodes[node].bonds.clone();
        self.nodes[actual].bonds.extend_from_slice(&node_bonds);

        // repair the bonds the dummy atom left dangling
        for point in &attachment_points {
            let bond = self.nodes[actual]
                .bonds
                .iter()
                .copied()
                .find(|&b| self.bonds[b].id == point.bond)
                .ok_or_else(|| {
                    file_parse(
                        "Bond specification for an AttachmentPoint definition was not found in the bond array in MRV file",
                    )
                })?;
            self.nodes[expanded_id].bonds.push(bond);
            let mut fixed = false;
            for k in 0..2 {
                let present = self.nodes[actual]
                    .atoms
                    .iter()
                    .any(|&a| self.atoms[a].id == self.bonds[bond].atom_refs2[k]);
                if !present {
                    self.bonds[bond].atom_refs2[k] = point.atom.clone();
                    fixed = true;
                    break;
                }
            }
            if !fixed {
                return Err(file_parse(format!(
                    "Bond {} from attachment point did not have a missing atom in MRV file",
                    point.bond
                )));
            }
        }

        if let SgroupKind::Superatom {
            attachment_points, ..
        } = &mut self.nodes[node].kind
        {
            attachment_points.clear();
        }
        self.nodes[node].atoms.clear();
        self.nodes[node].bonds.clear();

        for child in self.nodes[node].children.clone() {
            self.move_sgroup(child, actual, false)?;
        }
        self.nodes[node].children.clear();
        let p = self.nodes[parent]
            .children
            .iter()
            .position(|&c| self.nodes[c].id == sg_id)
            .ok_or_else(|| RDKitError::invariant("superatom not found in its parent"))?;
        self.nodes[parent].children.remove(p);
        Ok(())
    }
}
