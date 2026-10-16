//! The stereo configurations that get labelled: tetrahedral centers, double
//! bonds and atropisomeric bonds (upstream
//! `CIPLabeler/configs/{Configuration,Tetrahedral,Sp2Bond,AtropisomerBond}.{h,cpp}`).
//!
//! Upstream a configuration owns its digraph. Here the digraph is kept
//! beside it, since labelling the auxiliary descriptors of one digraph needs
//! the other configurations.

use super::cip_mol::CipMol;
use super::descriptor::Descriptor;
use super::digraph::{BOND_DUPLICATE, Digraph, EdgeId, IMPL_HYDROGEN, NodeId};
use super::rules::Rules;
use crate::graphmol::atom::ChiralType;
use crate::graphmol::atropisomers::get_atropisomer_atoms_and_bonds;
use crate::graphmol::bond::{BondStereo, BondType};
use crate::graphmol::chirality::bond_dirs::find_stereo_atoms;
use crate::graphmol::romol::ROMol;
use crate::rdgeneral::{RDProps, common_properties as cp};
use crate::{RDKitError, Result};

/// `Configuration::parity4`: 2 if the targets are an even permutation of the
/// reference, 1 if odd, 0 if they are not a permutation of it.
pub fn parity4<T: PartialEq>(trg: &[T], reference: &[T]) -> Result<i32> {
    if reference.len() != 4 || trg.len() != reference.len() {
        return Err(RDKitError::Runtime(
            "Parity vectors must have size 4.".into(),
        ));
    }
    // positions of the reference elements within the targets
    let mut perm = [0usize; 4];
    for (i, r) in reference.iter().enumerate() {
        match trg.iter().position(|t| t == r) {
            Some(p) => perm[i] = p,
            None => return Ok(0),
        }
    }
    // upstream only accepts true permutations (equal elements make the
    // first match ambiguous)
    let mut seen = [false; 4];
    for &p in &perm {
        if seen[p] {
            return Ok(0);
        }
        seen[p] = true;
    }
    let inversions = (0..4)
        .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
        .filter(|&(i, j)| perm[i] > perm[j])
        .count();
    Ok(if inversions % 2 == 0 { 2 } else { 1 })
}

#[derive(Debug, Clone, Copy)]
enum Kind {
    Tetrahedral,
    Sp2Bond { bond: u32, cfg: BondStereo },
    AtropisomerBond { bond: u32, cfg: BondStereo },
}

#[derive(Debug)]
pub struct Configuration {
    kind: Kind,
    foci: Vec<u32>,
    carriers: Vec<Option<u32>>,
    ranked_anchors: Vec<u32>,
    has_primary_label: bool,
}

impl Configuration {
    /// `Tetrahedral::Tetrahedral`
    pub fn tetrahedral(mol: &CipMol, focus: u32) -> Result<Self> {
        let tag = mol.atom(focus).chiral_tag();
        if tag != ChiralType::TetrahedralCCW && tag != ChiralType::TetrahedralCW {
            return Err(RDKitError::invariant("bad config"));
        }
        let mut carriers: Vec<Option<u32>> = mol.neighbors_of(focus).map(Some).collect();
        if carriers.len() < 2 || carriers.len() > 4 {
            carriers.clear();
        } else {
            if carriers.len() < 4 {
                carriers.push(Some(focus));
            }
            if carriers.len() < 4 {
                carriers.push(None);
            }
        }
        Ok(Configuration {
            kind: Kind::Tetrahedral,
            foci: vec![focus],
            carriers,
            ranked_anchors: Vec::new(),
            has_primary_label: false,
        })
    }

    /// `Sp2Bond::Sp2Bond`
    pub fn sp2_bond(mol: &CipMol, bond: u32, cfg: BondStereo) -> Result<Self> {
        let b = mol.bond(bond);
        let (start, end) = (b.begin_atom_idx(), b.end_atom_idx());
        if cfg != BondStereo::StereoTrans && cfg != BondStereo::StereoCis {
            return Err(RDKitError::invariant("bad config"));
        }
        let mut config = Configuration {
            kind: Kind::Sp2Bond { bond, cfg },
            foci: vec![start, end],
            carriers: Vec::new(),
            ranked_anchors: Vec::new(),
            has_primary_label: false,
        };
        if b.bond_type() != BondType::Double {
            return Ok(config);
        }
        let stereo_atoms = find_stereo_atoms(mol.mol(), bond)?;
        if stereo_atoms.len() != 2 {
            return Err(RDKitError::invariant("incorrect number of stereo atoms"));
        }
        let is_valid_carrier = |focus: u32, other_focus: u32, carrier: i32| {
            carrier >= 0
                && (carrier as u32) < mol.num_atoms()
                && carrier as u32 != focus
                && carrier as u32 != other_focus
                && mol.neighbors_of(focus).any(|n| n == carrier as u32)
        };
        if stereo_atoms[0] == stereo_atoms[1]
            || !is_valid_carrier(start, end, stereo_atoms[0])
            || !is_valid_carrier(end, start, stereo_atoms[1])
        {
            return Ok(config);
        }
        config.carriers = vec![Some(stereo_atoms[0] as u32), Some(stereo_atoms[1] as u32)];
        Ok(config)
    }

    /// `AtropisomerBond::AtropisomerBond`
    pub fn atropisomer_bond(mol: &CipMol, bond: u32, cfg: BondStereo) -> Result<Self> {
        let b = mol.bond(bond);
        if cfg != BondStereo::StereoAtropCW && cfg != BondStereo::StereoAtropCCW {
            return Err(RDKitError::invariant("bad config"));
        }
        let mut config = Configuration {
            kind: Kind::AtropisomerBond { bond, cfg },
            foci: vec![b.begin_atom_idx(), b.end_atom_idx()],
            carriers: Vec::new(),
            ranked_anchors: Vec::new(),
            has_primary_label: false,
        };
        let Some(ends) = get_atropisomer_atoms_and_bonds(mol.mol(), bond) else {
            return Ok(config); // not an atropisomer
        };
        if ends[0].1.is_empty() || ends[1].1.is_empty() {
            return Ok(config);
        }
        let carrier = |end: &(u32, Vec<u32>)| mol.other_atom(end.1[0], end.0);
        config.carriers = vec![Some(carrier(&ends[0])?), Some(carrier(&ends[1])?)];
        Ok(config)
    }

    pub fn foci(&self) -> &[u32] {
        &self.foci
    }

    pub fn carriers(&self) -> &[Option<u32>] {
        &self.carriers
    }

    pub fn has_primary_label(&self) -> bool {
        self.has_primary_label
    }

    /// The digraph of the configuration
    pub fn new_digraph<'a>(&self, mol: &'a CipMol) -> Result<Digraph<'a>> {
        let atropisomer = matches!(self.kind, Kind::AtropisomerBond { .. });
        Digraph::new(mol, self.foci[0], atropisomer)
    }

    /// `setPrimaryLabel`: stores the label on the atom or bond.
    pub fn set_primary_label(&mut self, mol: &mut ROMol, desc: Descriptor) -> Result<()> {
        let anchors = self.ranked_anchors.clone();
        match self.kind {
            Kind::Tetrahedral => match desc {
                Descriptor::R | Descriptor::S | Descriptor::PseudoR | Descriptor::PseudoS => {
                    let atom = mol.atom_mut(self.foci[0]);
                    atom.set_prop(cp::_CIPCode, desc.as_str(), false)?;
                    atom.set_prop(cp::_CIPNeighborOrder, anchors, true)?;
                }
                Descriptor::SeqTrans
                | Descriptor::SeqCis
                | Descriptor::E
                | Descriptor::Z
                | Descriptor::M
                | Descriptor::P
                | Descriptor::PseudoM
                | Descriptor::PseudoP
                | Descriptor::Sp4
                | Descriptor::Tbpy5
                | Descriptor::Oc6 => {
                    return Err(RDKitError::Runtime(
                        "Received a Descriptor that is not supported for atoms".into(),
                    ));
                }
                _ => {
                    return Err(RDKitError::Runtime(
                        "Received an invalid Atom Descriptor".into(),
                    ));
                }
            },
            Kind::Sp2Bond { bond, cfg } => match desc {
                Descriptor::SeqTrans | Descriptor::E | Descriptor::SeqCis | Descriptor::Z => {
                    let carriers: Vec<i32> = self
                        .carriers
                        .iter()
                        .map(|c| c.map_or(-1, |c| c as i32))
                        .collect();
                    let b = mol.bond_mut(bond);
                    *b.stereo_atoms_mut() = carriers;
                    b.set_stereo(cfg)?;
                    b.set_prop(cp::_CIPCode, desc.as_str(), false)?;
                    b.set_prop(cp::_CIPNeighborOrder, anchors, true)?;
                }
                Descriptor::R
                | Descriptor::S
                | Descriptor::PseudoR
                | Descriptor::PseudoS
                | Descriptor::M
                | Descriptor::P
                | Descriptor::PseudoM
                | Descriptor::PseudoP
                | Descriptor::Sp4
                | Descriptor::Tbpy5
                | Descriptor::Oc6 => {
                    return Err(RDKitError::Runtime(
                        "Received a Descriptor that is not supported for double bonds".into(),
                    ));
                }
                _ => {
                    return Err(RDKitError::Runtime(
                        "Received an invalid Bond Descriptor".into(),
                    ));
                }
            },
            Kind::AtropisomerBond { bond, .. } => match desc {
                Descriptor::M | Descriptor::P | Descriptor::PseudoM | Descriptor::PseudoP => {
                    let b = mol.bond_mut(bond);
                    b.set_prop(cp::_CIPCode, desc.as_str(), false)?;
                    b.set_prop(cp::_CIPNeighborOrder, anchors, true)?;
                }
                Descriptor::R
                | Descriptor::S
                | Descriptor::PseudoR
                | Descriptor::PseudoS
                | Descriptor::Sp4
                | Descriptor::Tbpy5
                | Descriptor::Oc6
                | Descriptor::SeqTrans
                | Descriptor::E
                | Descriptor::SeqCis
                | Descriptor::Z => {
                    return Err(RDKitError::Runtime(
                        "Received a Descriptor that is not supported for atropisomer bonds".into(),
                    ));
                }
                _ => {
                    return Err(RDKitError::Runtime(
                        "Received an invalid Bond Descriptor".into(),
                    ));
                }
            },
        }
        self.has_primary_label = true;
        Ok(())
    }

    /// `label(comp)`: labels the configuration seen from its own digraph.
    pub fn label(&mut self, d: &mut Digraph<'_>, comp: &Rules) -> Result<Descriptor> {
        let root = d.original_root();
        if d.current_root() != root {
            d.change_root(root)?;
        }
        self.label_at(d, root, comp)
    }

    /// `label(node, digraph, comp)`: labels the configuration as it is met at
    /// a node of the digraph of another configuration.
    pub fn label_in(
        &mut self,
        d: &mut Digraph<'_>,
        node: NodeId,
        comp: &Rules,
    ) -> Result<Descriptor> {
        if matches!(self.kind, Kind::Tetrahedral) {
            d.change_root(node)?;
        }
        self.label_at(d, node, comp)
    }

    fn label_at(&mut self, d: &mut Digraph<'_>, node: NodeId, comp: &Rules) -> Result<Descriptor> {
        match self.kind {
            Kind::Tetrahedral => self.label_tetrahedral(d, node, comp),
            Kind::Sp2Bond { .. } | Kind::AtropisomerBond { .. } => self.label_bond(d, node, comp),
        }
    }

    fn label_tetrahedral(
        &mut self,
        d: &mut Digraph<'_>,
        node: NodeId,
        comp: &Rules,
    ) -> Result<Descriptor> {
        let focus = self.foci[0];
        let mut edges = d.edges(node)?;
        self.ranked_anchors.clear();
        if self.carriers.len() != 4 || edges.len() < 3 {
            return Ok(Descriptor::Ns);
        }
        let mut priority = comp.sort(d, node, &mut edges)?;
        let is_unique = priority.is_unique();
        if !is_unique && edges.len() == 4 {
            if comp.num_sub_rules() == 3 {
                return Ok(Descriptor::Unknown);
            }
            let partition = comp.groups(d, &edges)?;
            if partition.len() == 2 {
                d.set_rule6_ref(d.end_of(edges[1]).atom);
                priority = comp.sort(d, node, &mut edges)?;
                d.set_rule6_ref(None);
            } else if partition.len() == 1 {
                d.set_rule6_ref(d.end_of(edges[0]).atom);
                comp.sort(d, node, &mut edges)?;
                let nbrs1 = edges.clone();
                d.set_rule6_ref(d.end_of(edges[1]).atom);
                priority = comp.sort(d, node, &mut edges)?;
                d.set_rule6_ref(None);
                if parity4(&nbrs1, &edges)? == 1 {
                    return Ok(Descriptor::Unknown);
                }
            }
            if !priority.is_unique() {
                return Ok(Descriptor::Unknown);
            }
        } else if !is_unique {
            return Ok(Descriptor::Unknown);
        }

        let mut ordered: Vec<Option<u32>> = vec![None; 4];
        let mut idx = 0;
        for &edge in &edges {
            let end = d.end_of(edge);
            if end.is_set(BOND_DUPLICATE) || end.is_set(IMPL_HYDROGEN) {
                continue;
            }
            if idx >= ordered.len() {
                return Err(RDKitError::Runtime(
                    "Could not calculate parity! invalid atom index".into(),
                ));
            }
            let atom = end.atom.unwrap();
            ordered[idx] = Some(atom);
            self.ranked_anchors.push(atom);
            idx += 1;
        }
        if idx < 4 {
            ordered[idx] = Some(focus);
        }
        let parity = parity4(&ordered, &self.carriers)?;
        if parity == 0 {
            return Err(RDKitError::Runtime(
                "Could not calculate parity! Carrier mismatch".into(),
            ));
        }
        let mut config = d.mol().atom(focus).chiral_tag();
        if parity == 1 {
            config = if config == ChiralType::TetrahedralCCW {
                ChiralType::TetrahedralCW
            } else {
                ChiralType::TetrahedralCCW
            };
        }
        let pseudo = priority.is_pseudo_asymmetric();
        Ok(match config {
            ChiralType::TetrahedralCCW => {
                if pseudo {
                    Descriptor::PseudoS
                } else {
                    Descriptor::S
                }
            }
            ChiralType::TetrahedralCW => {
                if pseudo {
                    Descriptor::PseudoR
                } else {
                    Descriptor::R
                }
            }
            _ => Descriptor::Unknown,
        })
    }

    /// `Sp2Bond::label` and `AtropisomerBond::label`, which differ in how
    /// the carriers flip the configuration and in the descriptors.
    fn label_bond(
        &mut self,
        d: &mut Digraph<'_>,
        root1: NodeId,
        comp: &Rules,
    ) -> Result<Descriptor> {
        let (focus1, focus2) = (self.foci[0], self.foci[1]);
        let atropisomer = matches!(self.kind, Kind::AtropisomerBond { .. });
        let is_constitutional = comp.num_sub_rules() == 3;
        self.ranked_anchors.clear();

        let root1_edges = d.edges(root1)?;
        let Some(internal) = find_internal_edge(d, &root1_edges, focus1, focus2) else {
            return Ok(Descriptor::Unknown);
        };
        let root2 = d.edge(internal).other(root1)?;
        let mut edges1 = d.edges(root1)?;
        let mut edges2 = d.edges(root2)?;
        remove_internal_edges(d, &mut edges1, focus1, focus2);
        remove_internal_edges(d, &mut edges2, focus1, focus2);
        if atropisomer {
            remove_duplicates_and_hs(d, &mut edges1);
            remove_duplicates_and_hs(d, &mut edges2);
        }
        if self.carriers.len() != 2 || edges1.is_empty() || edges2.is_empty() {
            return Ok(Descriptor::Ns);
        }
        let mut carriers = self.carriers.clone();
        // whether the configuration is cis (or counter clockwise)
        let mut positive = match self.kind {
            Kind::Sp2Bond { cfg, .. } => cfg == BondStereo::StereoCis,
            Kind::AtropisomerBond { cfg, .. } => cfg == BondStereo::StereoAtropCCW,
            Kind::Tetrahedral => unreachable!(),
        };
        if d.node(root1).atom == Some(focus2) {
            carriers.swap(0, 1);
        }

        d.change_root(root1)?;
        let priority1 = comp.sort(d, root1, &mut edges1)?;
        if !priority1.is_unique() && !is_constitutional {
            return Ok(Descriptor::Unknown);
        }
        if edges1.len() > 1 && self.flips(&carriers, 0, d, &edges1) {
            positive = !positive;
        }

        d.change_root(root2)?;
        let priority2 = comp.sort(d, root2, &mut edges2)?;
        if !priority2.is_unique() || !priority1.is_unique() {
            return Ok(Descriptor::Unknown);
        }
        if edges2.len() > 1 && self.flips(&carriers, 1, d, &edges2) {
            positive = !positive;
        }

        let carrier1_idx = d.end_of(edges1[0]).atom_idx();
        let carrier2_idx = d.end_of(edges2[0]).atom_idx();
        if d.beg_of(edges1[0]).atom == Some(focus1) {
            self.ranked_anchors = vec![carrier1_idx, carrier2_idx];
        } else if d.beg_of(edges2[0]).atom == Some(focus1) {
            self.ranked_anchors = vec![carrier2_idx, carrier1_idx];
        }

        let pseudo = priority1.is_pseudo_asymmetric() != priority2.is_pseudo_asymmetric();
        Ok(match (atropisomer, positive, pseudo) {
            (false, true, true) => Descriptor::SeqCis,
            (false, true, false) => Descriptor::Z,
            (false, false, true) => Descriptor::SeqTrans,
            (false, false, false) => Descriptor::E,
            (true, true, true) => Descriptor::PseudoM,
            (true, true, false) => Descriptor::M,
            (true, false, true) => Descriptor::PseudoP,
            (true, false, false) => Descriptor::P,
        })
    }

    /// Whether the highest ranked edge at a side of the bond makes the
    /// configuration flip: double bonds flip when the carrier is not the
    /// highest ranked substituent, atropisomers when it is the second one.
    fn flips(
        &self,
        carriers: &[Option<u32>],
        side: usize,
        d: &Digraph<'_>,
        edges: &[EdgeId],
    ) -> bool {
        match self.kind {
            Kind::Sp2Bond { .. } => carriers[side] != d.end_of(edges[0]).atom,
            _ => carriers[side] == d.end_of(edges[1]).atom,
        }
    }
}

/// `findInternalEdge`: the edge between the two foci
fn find_internal_edge(d: &Digraph<'_>, edges: &[EdgeId], f1: u32, f2: u32) -> Option<EdgeId> {
    edges.iter().copied().find(|&e| {
        !d.beg_of(e).is_duplicate() && !d.end_of(e).is_duplicate() && is_internal_edge(d, e, f1, f2)
    })
}

fn is_internal_edge(d: &Digraph<'_>, edge: EdgeId, f1: u32, f2: u32) -> bool {
    let (beg, end) = (d.beg_of(edge).atom, d.end_of(edge).atom);
    (Some(f1) == beg && Some(f2) == end) || (Some(f1) == end && Some(f2) == beg)
}

fn remove_internal_edges(d: &Digraph<'_>, edges: &mut Vec<EdgeId>, f1: u32, f2: u32) {
    edges.retain(|&e| !is_internal_edge(d, e, f1, f2));
}

fn remove_duplicates_and_hs(d: &Digraph<'_>, edges: &mut Vec<EdgeId>) {
    edges.retain(|&e| !d.beg_of(e).is_duplicate_or_h() && !d.end_of(e).is_duplicate_or_h());
}
