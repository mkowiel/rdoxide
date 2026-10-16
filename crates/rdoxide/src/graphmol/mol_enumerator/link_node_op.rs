//! Link node enumeration (upstream `MolEnumerator/LinkNode.cpp`): the group
//! of atoms between the two bonds of a link node is repeated.

use std::collections::BTreeMap;

use super::link_node::get_mol_link_nodes_mapped;
use super::{IDX_PROP_NAME, MolEnumeratorOp, preserve_orig_indices};
use crate::graphmol::chem_reactions::daylight::{ReactionSmartsParserParams, reaction_from_smarts};
use crate::graphmol::chem_reactions::runner::run_reactant;
use crate::graphmol::find_rings::get_mol_frags;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::log::LogStateSetter;
use crate::rdgeneral::{RDProps, RDValue, common_properties as cp};
use crate::{RDKitError, Result};

/// `LinkNodeOp`
#[derive(Debug, Clone, Default)]
pub struct LinkNodeOp {
    mol: Option<ROMol>,
    frame: Option<ROMol>,
    count_at_each_point: Vec<usize>,
    /// (isotope of the varying atom, of the first and of the second attachment)
    variations: Vec<(u32, u32, u32)>,
    point_ranges: Vec<(u32, u32)>,
    isotope_map: BTreeMap<u32, u32>,
}

impl LinkNodeOp {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_mol(mol: &ROMol) -> Result<Self> {
        let mut op = Self::new();
        op.init_from_mol(mol)?;
        Ok(op)
    }
}

impl MolEnumeratorOp for LinkNodeOp {
    fn init_from_mol(&mut self, mol: &ROMol) -> Result<()> {
        let mut mol = mol.clone();
        self.variations.clear();
        self.point_ranges.clear();
        self.isotope_map.clear();
        if !mol.has_prop(cp::molFileLinkNodes) {
            self.mol = Some(mol);
            return Ok(());
        }
        if !mol.has_prop(IDX_PROP_NAME) {
            preserve_orig_indices(&mut mol)?;
        }
        let mut atom_map: BTreeMap<u32, u32> = BTreeMap::new();
        for idx in 0..mol.num_atoms() {
            if let Some(RDValue::UnsignedInt(orig)) = mol.atom(idx).get_prop(IDX_PROP_NAME) {
                atom_map.insert(*orig, idx);
            }
        }
        if get_mol_frags(&mol).len() > 1 {
            return Err(RDKitError::Value(
                "LINKNODE enumeration only supported for molecules with a single fragment.".into(),
            ));
        }
        let mut frame = mol.clone();
        let nodes = get_mol_link_nodes_mapped(&frame, true, Some(&atom_map))?;
        for node in nodes {
            if node.n_bonds != 2 {
                return Err(RDKitError::Runtime(
                    "only link nodes with 2 bonds are currently supported".into(),
                ));
            }
            self.count_at_each_point
                .push((node.max_rep - node.min_rep + 1) as usize);
            self.point_ranges.push((node.min_rep, node.max_rep));
            let [(var, attach1), (_, attach2)] = [node.bond_atoms[0], node.bond_atoms[1]];
            for atom in [var, attach1, attach2] {
                let key = 1000 * (atom + 1);
                self.isotope_map
                    .entry(key)
                    .or_insert_with(|| frame.atom(atom).isotope());
            }
            for atom in [var, attach1, attach2] {
                frame.atom_mut(atom).set_isotope(1000 * (atom + 1));
            }
            self.variations
                .push((1000 * (var + 1), 1000 * (attach1 + 1), 1000 * (attach2 + 1)));
        }
        self.mol = Some(mol);
        self.frame = Some(frame);
        Ok(())
    }

    fn variation_counts(&self) -> Vec<usize> {
        self.count_at_each_point.clone()
    }

    fn apply(&self, which: &[usize]) -> Result<ROMol> {
        if self.mol.is_none() {
            return Err(RDKitError::precondition("no molecule"));
        }
        let Some(frame) = &self.frame else {
            return Err(RDKitError::precondition("not initialized"));
        };
        if which.len() != self.count_at_each_point.len() {
            return Err(RDKitError::Value(
                "bad element choice in enumeration".into(),
            ));
        }
        for (w, count) in which.iter().zip(&self.count_at_each_point) {
            if w >= count {
                return Err(RDKitError::Value("bad element value in enumeration".into()));
            }
        }
        let mut res = frame.clone();
        for (i, &w) in which.iter().enumerate() {
            let variation_idx = i + 1;
            let variation_count = self.point_ranges[i].0 as usize + w;
            let (var, attach1, attach2) = self.variations[i];
            let reacts = format!(
                "[{attach1}*:{}]-[{var}*:{}]-[{attach2}*:{}]",
                variation_idx * 100 + 1,
                variation_idx * 100,
                variation_idx * 100 + 2
            );
            let mut prods = format!("[*:{}]", variation_idx * 100 + 1);
            for _ in 0..variation_count - 1 {
                prods += &format!("-[*:{}]", variation_idx * 100);
            }
            prods += &format!(
                "-[{var}*:{}]-[*:{}]",
                variation_idx * 100,
                variation_idx * 100 + 2
            );
            let mut rxn = reaction_from_smarts(
                &format!("{reacts}>>{prods}"),
                &ReactionSmartsParserParams::new(),
            )?;
            rxn.init_reactant_matchers(true);
            let products = {
                let _blocker = LogStateSetter::new();
                run_reactant(&rxn, &res, 0, 1)?
            };
            let Some(first) = products.into_iter().next() else {
                return Err(RDKitError::invariant("no products from reaction"));
            };
            if first.len() != 1 {
                return Err(RDKitError::invariant("too many products from reaction"));
            }
            res = first.into_iter().next().unwrap();
        }
        for idx in 0..res.num_atoms() {
            let isotope = res.atom(idx).isotope();
            if isotope >= 1000
                && let Some(&orig) = self.isotope_map.get(&isotope)
            {
                res.atom_mut(idx).set_isotope(orig);
            }
        }
        Ok(res)
    }

    fn copy(&self) -> Box<dyn MolEnumeratorOp> {
        Box::new(self.clone())
    }
}
