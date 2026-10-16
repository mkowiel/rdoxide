//! Port of RDKit `Code/GraphMol/CIPLabeler`: stereo labels (`R`/`S`, `E`/`Z`,
//! `M`/`P`) from an implementation of the CIP rules that ranks the
//! substituents with hierarchical digraphs. It is a port of
//! <https://github.com/SiMolecule/centres> by John Mayfield, described in
//! Hanson, Musacchio, Mayfield, Vainio, Yerin, Redkin, *Algorithmic Analysis
//! of Cahn--Ingold--Prelog Rules of Stereochemistry: Proposals for Revised
//! Rules and a Guide for Machine Implementation*, J. Chem. Inf. Model. 2018,
//! 58, 1755-1765.
//!
//! Not ported: the old rule 5 (`Rule5`), which upstream replaces with
//! `Rule5New`, and the interruption with Ctrl-C.

mod budget;
mod cip_mol;
mod configs;
mod descriptor;
mod digraph;
mod mancude;
mod rules;

pub use descriptor::Descriptor;

use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondStereo, BondType};
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::{RDProps, common_properties as cp};
use crate::{RDKitError, Result};
use budget::{ScopedIterationBudget, ScopedPreliminaryBudget, is_max_iterations_exceeded};
use cip_mol::CipMol;
use configs::Configuration;
use digraph::{Digraph, NodeId};
use rules::Rules;
use std::collections::HashMap;

/// Whether an error is the one raised when the limit of recursive
/// comparisons is exceeded (`MaxIterationsExceeded`).
pub fn is_max_iterations_exceeded_error(e: &RDKitError) -> bool {
    is_max_iterations_exceeded(e)
}

/// A configuration with its digraph and whether it is to be labelled
struct Entry<'a> {
    config: Configuration,
    digraph: Digraph<'a>,
    selected: bool,
}

fn is_selected(selection: &[bool], index: u32) -> bool {
    selection.get(index as usize).copied().unwrap_or(false)
}

/// `findConfigs`
fn find_configs<'a>(mol: &'a CipMol, atoms: &[bool], bonds: &[bool]) -> Result<Vec<Entry<'a>>> {
    let mut entries = Vec::new();
    let mut push = |config: Configuration, selected: bool| -> Result<()> {
        let digraph = config.new_digraph(mol)?;
        entries.push(Entry {
            config,
            digraph,
            selected,
        });
        Ok(())
    };
    for index in 0..mol.num_atoms() {
        let tag = mol.atom(index).chiral_tag();
        if tag == ChiralType::TetrahedralCW || tag == ChiralType::TetrahedralCCW {
            let config = Configuration::tetrahedral(mol, index)?;
            if config.carriers().len() == 4 {
                push(config, is_selected(atoms, index))?;
            }
        }
    }
    for index in 0..mol.num_bonds() {
        let bond = mol.bond(index);
        let stereo = match bond.stereo() {
            BondStereo::StereoE => BondStereo::StereoTrans,
            BondStereo::StereoZ => BondStereo::StereoCis,
            s => s,
        };
        match stereo {
            BondStereo::StereoTrans | BondStereo::StereoCis => {
                if bond.bond_type() != BondType::Double {
                    continue;
                }
                let config = Configuration::sp2_bond(mol, index, stereo)?;
                if config.carriers().len() == 2 {
                    push(config, is_selected(bonds, index))?;
                }
            }
            BondStereo::StereoAtropCCW | BondStereo::StereoAtropCW => {
                if bond.bond_type() != BondType::Single {
                    continue;
                }
                let config = Configuration::atropisomer_bond(mol, index, stereo)?;
                if config.carriers().len() == 2 {
                    push(config, is_selected(bonds, index))?;
                }
            }
            _ => {}
        }
    }
    Ok(entries)
}

/// `labelAux`: labels the other configurations met in the digraph of the
/// center, which the rules 3 to 5 compare.
fn label_aux(entries: &mut [Entry<'_>], rules: &Rules, center: usize) -> Result<()> {
    let mut aux: Vec<(NodeId, usize)> = Vec::new();
    for other in 0..entries.len() {
        if other == center {
            continue;
        }
        let foci = entries[other].config.foci().to_vec();
        let digraph = &mut entries[center].digraph;
        if !foci.iter().any(|&f| digraph.seen_atom(f)) {
            continue;
        }
        for node in digraph.nodes_of_atom(foci[0])? {
            if digraph.node(node).is_duplicate() {
                continue;
            }
            let mut low = node;
            if foci.len() == 2 {
                for edge in digraph.edges_to_atom(node, foci[1])? {
                    let other_node = digraph.edge(edge).other(node)?;
                    if digraph.node(other_node).dist < digraph.node(node).dist {
                        low = other_node;
                    }
                }
            }
            if !digraph.node(low).is_duplicate() {
                aux.push((low, other));
            }
        }
    }
    // the farthest first
    let digraph_dist = |d: &Digraph<'_>, n: NodeId| d.node(n).dist;
    {
        let d = &entries[center].digraph;
        aux.sort_by_key(|&(n, _)| std::cmp::Reverse(digraph_dist(d, n)));
    }
    let mut queue: HashMap<NodeId, Descriptor> = HashMap::new();
    let mut prev = i32::MAX;
    for (node, other) in aux {
        let dist = digraph_dist(&entries[center].digraph, node);
        if dist < prev {
            for (n, desc) in queue.drain() {
                entries[center].digraph.set_node_aux(n, desc);
            }
            prev = dist;
        }
        // the digraph of the center and the configuration of another entry
        let (digraph, config) = if other < center {
            let (head, tail) = entries.split_at_mut(center);
            (&mut tail[0].digraph, &mut head[other].config)
        } else {
            let (head, tail) = entries.split_at_mut(other);
            (&mut head[center].digraph, &mut tail[0].config)
        };
        let label = config.label_in(digraph, node, rules)?;
        queue.insert(node, label);
    }
    for (n, desc) in queue {
        entries[center].digraph.set_node_aux(n, desc);
    }
    Ok(())
}

/// `label`
fn label(mol: &mut ROMol, entries: &mut [Entry<'_>], max_recursive_iterations: u32) -> Result<()> {
    let constitutional_rules = Rules::constitutional();
    let all_rules = Rules::all();
    let _budget = ScopedIterationBudget::new(max_recursive_iterations);
    for entry in entries.iter_mut().filter(|e| e.selected) {
        let _preliminary = ScopedPreliminaryBudget::new();
        match entry
            .config
            .label(&mut entry.digraph, &constitutional_rules)
        {
            Ok(desc) => {
                if desc != Descriptor::Unknown {
                    entry.config.set_primary_label(mol, desc)?;
                }
            }
            Err(e) if is_max_iterations_exceeded(&e) => {}
            Err(e) => return Err(e),
        }
    }
    for i in 0..entries.len() {
        if !entries[i].selected || entries[i].config.has_primary_label() {
            continue;
        }
        let entry = &mut entries[i];
        let mut desc = entry
            .config
            .label(&mut entry.digraph, &constitutional_rules)?;
        if desc != Descriptor::Unknown {
            entries[i].config.set_primary_label(mol, desc)?;
        } else {
            label_aux(entries, &all_rules, i)?;
            let entry = &mut entries[i];
            desc = entry.config.label(&mut entry.digraph, &all_rules)?;
            if desc != Descriptor::Unknown {
                entry.config.set_primary_label(mol, desc)?;
            }
        }
    }
    Ok(())
}

/// `clearSelectedCIPProperties`
fn clear_selected_cip_properties(
    mol: &mut ROMol,
    atoms: &[bool],
    bonds: &[bool],
    full_selection: bool,
) {
    mol.clear_prop(cp::_CIPComputed);
    for idx in 0..mol.num_atoms() {
        if full_selection || is_selected(atoms, idx) {
            let atom = mol.atom_mut(idx);
            atom.clear_prop(cp::_CIPCode);
            atom.clear_prop(cp::_CIPNeighborOrder);
        }
    }
    for idx in 0..mol.num_bonds() {
        if full_selection || is_selected(bonds, idx) {
            let bond = mol.bond_mut(idx);
            bond.clear_prop(cp::_CIPCode);
            bond.clear_prop(cp::_CIPNeighborOrder);
        }
    }
}

/// `assignCIPLabels` for the selected atoms and bonds.
///
/// Very symmetrical molecules can cause pseudo-infinite processing
/// (dodecahedrane, say); `max_recursive_iterations` limits the recursive
/// comparisons of the call (0: no limit), which fails with
/// "Max Iterations Exceeded in CIP label calculation" (see
/// [`is_max_iterations_exceeded_error`]) when exhausted.
pub fn assign_cip_labels_selected(
    mol: &mut ROMol,
    atoms: &[bool],
    bonds: &[bool],
    max_recursive_iterations: u32,
) -> Result<()> {
    if atoms.len() != mol.num_atoms() as usize {
        return Err(RDKitError::precondition(
            "Atoms bitset size does not match number of atoms",
        ));
    }
    if bonds.len() != mol.num_bonds() as usize {
        return Err(RDKitError::precondition(
            "Bonds bitset size does not match number of bonds",
        ));
    }
    let full_selection = atoms.iter().all(|&a| a) && bonds.iter().all(|&b| b);
    clear_selected_cip_properties(mol, atoms, bonds, full_selection);
    let cipmol = CipMol::new(mol)?;
    let mut entries = find_configs(&cipmol, atoms, bonds)?;
    label(mol, &mut entries, max_recursive_iterations)?;
    if full_selection {
        mol.set_prop(cp::_CIPComputed, true, true)?;
    }
    Ok(())
}

/// `assignCIPLabels`: labels every stereo center, double bond and
/// atropisomeric bond.
pub fn assign_cip_labels(mol: &mut ROMol, max_recursive_iterations: u32) -> Result<()> {
    let atoms = vec![true; mol.num_atoms() as usize];
    let bonds = vec![true; mol.num_bonds() as usize];
    assign_cip_labels_selected(mol, &atoms, &bonds, max_recursive_iterations)
}
