//! Port of RDKit `Code/GraphMol/ScaffoldNetwork`: scaffold networks built by
//! breaking bonds of molecules (after Varin et al., J. Chem. Inf. Model.
//! 2011). Boost-serialization pickles are not supported.

use std::collections::{BTreeSet, VecDeque};
use std::sync::Arc;

use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::chem_reactions::ChemicalReaction;
use crate::graphmol::chem_reactions::daylight::{ReactionSmartsParserParams, reaction_from_smarts};
use crate::graphmol::chem_reactions::runner::run_reactant;
use crate::graphmol::chem_transforms::murcko_decompose;
use crate::graphmol::find_rings::fast_find_rings;
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::mol_standardize::fragment::LargestFragmentChooser;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::write::mol_to_smiles;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result};

/// `ScaffoldNetworkParams`
#[derive(Debug, Clone)]
pub struct ScaffoldNetworkParams {
    /// include scaffolds with all atoms replaced by dummies
    pub include_generic_scaffolds: bool,
    /// include scaffolds with all bonds replaced by single bonds
    pub include_generic_bond_scaffolds: bool,
    /// remove attachment points from scaffolds and include the result
    pub include_scaffolds_without_attachments: bool,
    /// include the version of the scaffold with attachment points
    pub include_scaffolds_with_attachments: bool,
    /// include the names of the input molecules
    pub include_names: bool,
    /// keep only the first fragment from the bond breaking rule
    pub keep_only_first_fragment: bool,
    /// do a pruning/flattening step before starting fragmenting
    pub prune_before_fragmenting: bool,
    /// remove isotopes when flattening
    pub flatten_isotopes: bool,
    /// remove chirality and bond stereo when flattening
    pub flatten_chirality: bool,
    /// keep only the largest fragment when doing flattening
    pub flatten_keep_largest: bool,
    /// keep track of the number of molecules each scaffold was found in
    pub collect_mol_counts: bool,
    /// the reactions used to fragment: one reactant, two products
    pub bond_breakers_rxns: Vec<Arc<ChemicalReaction>>,
}

const DEFAULT_BOND_BREAKER: &str = "[!#0;R:1]-!@[!#0:2]>>[*:1]-[#0].[#0]-[*:2]";

impl Default for ScaffoldNetworkParams {
    fn default() -> Self {
        Self::new(&[DEFAULT_BOND_BREAKER]).expect("default bond breaker")
    }
}

impl ScaffoldNetworkParams {
    /// The parameters with the given bond breaking reactions (SMARTS).
    pub fn new<S: AsRef<str>>(bond_breakers_smarts: &[S]) -> Result<Self> {
        let mut rxns = Vec::with_capacity(bond_breakers_smarts.len());
        for sma in bond_breakers_smarts {
            let mut rxn = reaction_from_smarts(sma.as_ref(), &ReactionSmartsParserParams::new())?;
            if rxn.num_reactant_templates() != 1 {
                return Err(RDKitError::Value(
                    "bond breaker reactions must have exactly one reactant".into(),
                ));
            }
            if rxn.num_product_templates() != 2 {
                return Err(RDKitError::Value(
                    "bond breaker reactions must have exactly two products".into(),
                ));
            }
            rxn.init_reactant_matchers(false);
            rxns.push(Arc::new(rxn));
        }
        Ok(ScaffoldNetworkParams {
            include_generic_scaffolds: true,
            include_generic_bond_scaffolds: false,
            include_scaffolds_without_attachments: true,
            include_scaffolds_with_attachments: true,
            include_names: false,
            keep_only_first_fragment: true,
            prune_before_fragmenting: true,
            flatten_isotopes: true,
            flatten_chirality: true,
            flatten_keep_largest: true,
            collect_mol_counts: true,
            bond_breakers_rxns: rxns,
        })
    }
}

/// `EdgeType`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeType {
    /// molecule -> fragment
    Fragment = 1,
    /// molecule -> generic molecule (all atoms are dummies)
    Generic = 2,
    /// molecule -> generic bond molecule (all bonds single)
    GenericBond = 3,
    /// molecule -> molecule with no attachment points
    RemoveAttachment = 4,
    /// molecule -> flattened molecule
    Initialize = 5,
}

impl EdgeType {
    pub fn name(self) -> &'static str {
        match self {
            EdgeType::Fragment => "Fragment",
            EdgeType::Generic => "Generic",
            EdgeType::GenericBond => "GenericBond",
            EdgeType::RemoveAttachment => "RemoveAttachment",
            EdgeType::Initialize => "Initialize",
        }
    }
}

impl EdgeType {
    fn from_code(code: u32) -> Option<Self> {
        Some(match code {
            1 => EdgeType::Fragment,
            2 => EdgeType::Generic,
            3 => EdgeType::GenericBond,
            4 => EdgeType::RemoveAttachment,
            5 => EdgeType::Initialize,
            _ => return None,
        })
    }
}

/// `NetworkEdge`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NetworkEdge {
    pub begin_idx: usize,
    pub end_idx: usize,
    pub edge_type: EdgeType,
}

impl std::fmt::Display for NetworkEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "NetworkEdge( {}->{}, type:{} )",
            self.begin_idx,
            self.end_idx,
            self.edge_type.name()
        )
    }
}

/// `ScaffoldNetwork`
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ScaffoldNetwork {
    /// SMILES of the scaffolds
    pub nodes: Vec<String>,
    /// number of times each scaffold was encountered
    pub counts: Vec<u32>,
    /// number of molecules each scaffold was found in
    pub mol_counts: Vec<u32>,
    pub edges: Vec<NetworkEdge>,
}

/// Header of a Boost text archive (library archive version 19).
const ARCHIVE_HEADER: &str = "22 serialization::archive 19";
/// Archive layout of `ScaffoldNetwork` (class version 1) and of the first
/// `NetworkEdge` element, as written by `boost::archive::text_oarchive`.
const NETWORK_PREFIX: &str = "0 1 0 0";
const EDGE_VECTOR_PREFIX: &str = "0 0";
const EDGE_ITEM_PREFIX: &str = "0 0";

impl ScaffoldNetwork {
    /// The network as the Boost text archive upstream pickles.
    pub fn to_archive(&self) -> String {
        use std::fmt::Write;
        let mut out = format!("{ARCHIVE_HEADER} {NETWORK_PREFIX} {} 0", self.nodes.len());
        for n in &self.nodes {
            let _ = write!(out, " {} {n}", n.len());
        }
        for counts in [&self.counts, &self.mol_counts] {
            let _ = write!(out, " {} 0", counts.len());
            for c in counts {
                let _ = write!(out, " {c}");
            }
        }
        let _ = write!(out, " {EDGE_VECTOR_PREFIX} {} 0", self.edges.len());
        if !self.edges.is_empty() {
            let _ = write!(out, " {EDGE_ITEM_PREFIX}");
        }
        for e in &self.edges {
            let _ = write!(out, " {} {} {}", e.begin_idx, e.end_idx, e.edge_type as u32);
        }
        out
    }

    /// Reads a Boost text archive written by [`Self::to_archive`] or upstream.
    pub fn from_archive(text: &str) -> Result<Self> {
        let bad = || RDKitError::Value("invalid ScaffoldNetwork archive".into());
        let rest = text.strip_prefix(ARCHIVE_HEADER).ok_or_else(bad)?;
        let mut cur = ArchiveReader { rest };
        // object id / class info, class version
        let _ = cur.number()?;
        let version = cur.number()?;
        cur.number()?;
        cur.number()?;
        let mut net = ScaffoldNetwork::default();
        let n = cur.number()?;
        cur.number()?; // item version
        for _ in 0..n {
            net.nodes.push(cur.string()?);
        }
        let n = cur.number()?;
        cur.number()?;
        for _ in 0..n {
            net.counts.push(cur.number()?);
        }
        if version > 0 {
            let n = cur.number()?;
            cur.number()?;
            for _ in 0..n {
                net.mol_counts.push(cur.number()?);
            }
        }
        cur.number()?;
        cur.number()?;
        let n = cur.number()?;
        cur.number()?;
        if n > 0 {
            cur.number()?;
            cur.number()?;
        }
        for _ in 0..n {
            let (b, e, t) = (cur.number()?, cur.number()?, cur.number()?);
            net.edges.push(NetworkEdge {
                begin_idx: b as usize,
                end_idx: e as usize,
                edge_type: EdgeType::from_code(t).ok_or_else(bad)?,
            });
        }
        Ok(net)
    }
}

/// Whitespace-separated token reader for [`ScaffoldNetwork::from_archive`].
struct ArchiveReader<'a> {
    rest: &'a str,
}

impl ArchiveReader<'_> {
    fn number(&mut self) -> Result<u32> {
        let t = self.rest.trim_start();
        let end = t.find(char::is_whitespace).unwrap_or(t.len());
        self.rest = &t[end..];
        t[..end]
            .parse()
            .map_err(|_| RDKitError::Value("invalid ScaffoldNetwork archive".into()))
    }

    /// A length-prefixed string.
    fn string(&mut self) -> Result<String> {
        let len = self.number()? as usize;
        let t = self.rest.strip_prefix(' ').unwrap_or(self.rest);
        let s = t
            .get(..len)
            .ok_or_else(|| RDKitError::Value("invalid ScaffoldNetwork archive".into()))?;
        self.rest = &t[len..];
        Ok(s.to_string())
    }
}

fn smiles(mol: &ROMol) -> Result<String> {
    mol_to_smiles(&mut mol.clone(), &Default::default())
}

/// `detail::makeScaffoldGeneric`
pub fn make_scaffold_generic(mol: &ROMol, do_atoms: bool, do_bonds: bool) -> ROMol {
    let mut res = mol.clone();
    if do_atoms {
        for idx in 0..res.num_atoms() {
            let atom = res.atom_mut(idx);
            atom.set_atomic_num(0);
            atom.set_num_explicit_hs(0);
            atom.set_no_implicit(false);
            atom.set_isotope(0);
        }
    }
    if do_bonds {
        for idx in 0..res.num_bonds() {
            let bond = res.bond_mut(idx);
            bond.set_bond_type(BondType::Single);
            let (b, e) = (bond.begin_atom_idx(), bond.end_atom_idx());
            res.atom_mut(b).set_is_aromatic(false);
            res.atom_mut(e).set_is_aromatic(false);
        }
    }
    res
}

/// `detail::removeAttachmentPoints`
pub fn remove_attachment_points(mol: &ROMol) -> Result<ROMol> {
    let mut res = mol.clone();
    res.begin_batch_edit()?;
    for idx in 0..res.num_atoms() {
        if res.atom(idx).atomic_num() == 0 && res.atom_degree(idx) == 1 {
            let nbr = res.atom_neighbors(idx).next().unwrap();
            let nbr_atom = res.atom_mut(nbr);
            if nbr_atom.is_aromatic() && nbr_atom.atomic_num() > 0 && nbr_atom.atomic_num() != 6 {
                nbr_atom.set_no_implicit(true);
                nbr_atom.set_num_explicit_hs(1);
            }
            res.remove_atom(idx)?;
        }
    }
    res.commit_batch_edit()?;
    Ok(res)
}

/// `detail::pruneMol`: the Murcko scaffold
pub fn prune_mol(mol: &ROMol) -> Result<ROMol> {
    let mut res = murcko_decompose(&mut mol.clone())?;
    res.update_property_cache(true)?;
    fast_find_rings(&mut res)?;
    Ok(res)
}

/// `detail::flattenMol`
pub fn flatten_mol(mol: &ROMol, params: &ScaffoldNetworkParams) -> Result<ROMol> {
    let mut res = if params.flatten_keep_largest {
        LargestFragmentChooser::default().choose(mol)?
    } else {
        mol.clone()
    };
    for idx in 0..res.num_atoms() {
        let atom = res.atom_mut(idx);
        if params.flatten_isotopes {
            atom.set_isotope(0);
        }
        if params.flatten_chirality && atom.chiral_tag() != ChiralType::Unspecified {
            atom.set_chiral_tag(ChiralType::Unspecified);
            if atom.no_implicit() && matches!(atom.atomic_num(), 6 | 7 | 15 | 16) {
                atom.set_no_implicit(false);
                atom.set_num_explicit_hs(0);
            }
        }
    }
    if params.flatten_chirality {
        for idx in 0..res.num_bonds() {
            let bond = res.bond_mut(idx);
            bond.set_bond_dir(BondDir::None);
            bond.set_stereo(BondStereo::StereoNone)?;
        }
    }
    Ok(res)
}

fn update_mol_props(mol: &mut RWMol) -> Result<()> {
    sanitize_mol(mol, SANITIZE_ALL).map_err(|(_, e)| e)
}

/// `detail::getMolFragments`: (SMILES of the parent, fragment) pairs
pub fn get_mol_fragments(
    mol: &ROMol,
    params: &ScaffoldNetworkParams,
) -> Result<Vec<(String, ROMol)>> {
    let mut res = Vec::new();
    let mut stack: VecDeque<ROMol> = VecDeque::new();
    stack.push_back(mol.clone());
    let mut seen: Vec<String> = Vec::new();
    while let Some(wmol) = stack.pop_front() {
        let parent_smi = smiles(&wmol)?;
        for rxn in &params.bond_breakers_rxns {
            for products in run_reactant(rxn, &wmol, 0, 1)? {
                let n_keep = if params.keep_only_first_fragment {
                    1
                } else {
                    2
                };
                for mut p in products.into_iter().take(n_keep) {
                    if !params.include_scaffolds_with_attachments {
                        p = remove_attachment_points(&p)?;
                    }
                    update_mol_props(&mut p)?;
                    let smi = smiles(&p)?;
                    if !seen.contains(&smi) {
                        stack.push_back(p.clone());
                        seen.push(smi);
                    }
                    res.push((parent_smi.clone(), p));
                }
            }
        }
    }
    Ok(res)
}

/// `addEntryIfMissing` for the nodes: the index of the entry, which is
/// counted when `counts` is given.
fn add_node(nodes: &mut Vec<String>, smi: &str, counts: Option<&mut Vec<u32>>) -> usize {
    let idx = match nodes.iter().position(|n| n == smi) {
        Some(i) => i,
        None => {
            nodes.push(smi.to_string());
            nodes.len() - 1
        }
    };
    if let Some(counts) = counts {
        if counts.len() < nodes.len() {
            counts.resize(nodes.len(), 0);
        }
        counts[idx] += 1;
    }
    idx
}

fn add_edge(edges: &mut Vec<NetworkEdge>, begin_idx: usize, end_idx: usize, edge_type: EdgeType) {
    let edge = NetworkEdge {
        begin_idx,
        end_idx,
        edge_type,
    };
    if !edges.contains(&edge) {
        edges.push(edge);
    }
}

/// The generic scaffolds of `mol` (atoms, and then bonds also), added to the
/// network below the node `from`.
fn add_generics(
    network: &mut ScaffoldNetwork,
    entries: &mut BTreeSet<usize>,
    params: &ScaffoldNetworkParams,
    mol: &ROMol,
    from: usize,
) -> Result<()> {
    let gsmi = smiles(&make_scaffold_generic(mol, true, false))?;
    let gidx = add_node(&mut network.nodes, &gsmi, Some(&mut network.counts));
    if params.collect_mol_counts {
        entries.insert(gidx);
    }
    add_edge(&mut network.edges, from, gidx, EdgeType::Generic);
    if params.include_generic_bond_scaffolds {
        let gbsmi = smiles(&make_scaffold_generic(mol, true, true))?;
        let gbidx = add_node(&mut network.nodes, &gbsmi, Some(&mut network.counts));
        if params.collect_mol_counts {
            entries.insert(gbidx);
        }
        if gidx != gbidx {
            add_edge(&mut network.edges, gidx, gbidx, EdgeType::GenericBond);
        }
    }
    Ok(())
}

/// `detail::addMolToNetwork`
pub fn add_mol_to_network(
    mol: &ROMol,
    network: &mut ScaffoldNetwork,
    params: &ScaffoldNetworkParams,
) -> Result<()> {
    let mut ismi = smiles(mol)?;
    if params.include_names
        && let Some(name) = mol.get_prop("_Name").map(RDValue::to_rd_string)
    {
        ismi += &format!(" {name}");
    }
    let mut fmol = flatten_mol(mol, params)?;
    if params.prune_before_fragmenting {
        fmol = prune_mol(&fmol)?;
    }
    let mut entries_from_mol: BTreeSet<usize> = BTreeSet::new();
    let fsmi = smiles(&fmol)?;
    if ismi != fsmi {
        let iidx = add_node(&mut network.nodes, &ismi, Some(&mut network.counts));
        let fidx = add_node(&mut network.nodes, &fsmi, Some(&mut network.counts));
        if params.collect_mol_counts {
            entries_from_mol.insert(iidx);
            entries_from_mol.insert(fidx);
        }
        add_edge(&mut network.edges, iidx, fidx, EdgeType::Initialize);
        if params.include_generic_scaffolds {
            add_generics(network, &mut entries_from_mol, params, &fmol, fidx)?;
        }
    } else {
        let midx = add_node(&mut network.nodes, &fsmi, Some(&mut network.counts));
        if params.collect_mol_counts {
            entries_from_mol.insert(midx);
        }
    }

    for (smi, frag_mol) in get_mol_fragments(&fmol, params)? {
        let iidx = add_node(&mut network.nodes, &smi, None);
        let lsmi = smiles(&frag_mol)?;
        let lidx = add_node(&mut network.nodes, &lsmi, Some(&mut network.counts));
        if params.collect_mol_counts {
            entries_from_mol.insert(lidx);
        }
        add_edge(&mut network.edges, iidx, lidx, EdgeType::Fragment);
        if params.include_generic_scaffolds {
            add_generics(network, &mut entries_from_mol, params, &frag_mol, lidx)?;
        }
        if params.include_scaffolds_with_attachments && params.include_scaffolds_without_attachments
        {
            let amol = remove_attachment_points(&frag_mol)?;
            let asmi = smiles(&amol)?;
            let aidx = add_node(&mut network.nodes, &asmi, Some(&mut network.counts));
            if params.collect_mol_counts {
                entries_from_mol.insert(aidx);
            }
            add_edge(&mut network.edges, lidx, aidx, EdgeType::RemoveAttachment);
            if params.include_generic_scaffolds {
                add_generics(network, &mut entries_from_mol, params, &amol, aidx)?;
            }
        }
    }
    if params.collect_mol_counts
        && let Some(&max_val) = entries_from_mol.last()
    {
        if network.mol_counts.len() <= max_val {
            network.mol_counts.resize(max_val + 1, 0);
        }
        for idx in entries_from_mol {
            network.mol_counts[idx] += 1;
        }
    }
    Ok(())
}

/// `updateScaffoldNetwork`: adds molecules to a network.
pub fn update_scaffold_network<M: std::borrow::Borrow<ROMol>>(
    mols: &[M],
    network: &mut ScaffoldNetwork,
    params: &ScaffoldNetworkParams,
) -> Result<()> {
    if !params.include_scaffolds_with_attachments && !params.include_scaffolds_without_attachments {
        return Err(RDKitError::Value(
            "must include at least one of scaffolds with attachments or scaffolds without attachments"
                .into(),
        ));
    }
    for mol in mols {
        add_mol_to_network(mol.borrow(), network, params)?;
    }
    Ok(())
}

/// `createScaffoldNetwork`
pub fn create_scaffold_network<M: std::borrow::Borrow<ROMol>>(
    mols: &[M],
    params: &ScaffoldNetworkParams,
) -> Result<ScaffoldNetwork> {
    let mut res = ScaffoldNetwork::default();
    update_scaffold_network(mols, &mut res, params)?;
    Ok(res)
}

/// `BRICSDefinitions`: the BRICS bond breaking rules
const BRICS_DEFINITIONS: &[&str] = &[
    "[$([C;D3]([#0,#6,#7,#8])(=O)):1]-;!@[$([O;D2]-;!@[#0,#6,#1]):2]>>[1*]-[*:1].[3*]-[*:2]",
    "[$([C;D3]([#0,#6,#7,#8])(=O)):1]-;!@[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):2]>>[1*]-[*:1].[5*]-[*:2]",
    "[$([C;D3]([#0,#6,#7,#8])(=O)):1]-;!@[$([N;R;$(N(@C(=O))@[C,N,O,S])]):2]>>[1*]-[*:1].[10*]-[*:2]",
    "[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([C;!D1;!$(C=*)]-;!@[#6]):2]>>[3*]-[*:1].[4*]-[*:2]",
    "[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[3*]-[*:1].[13*]-[*:2]",
    "[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[3*]-[*:1].[14*]-[*:2]",
    "[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[3*]-[*:1].[15*]-[*:2]",
    "[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([c;$(c(:c):c)]):2]>>[3*]-[*:1].[16*]-[*:2]",
    "[$([C;!D1;!$(C=*)]-;!@[#6]):1]-;!@[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):2]>>[4*]-[*:1].[5*]-[*:2]",
    "[$([C;!D1;!$(C=*)]-;!@[#6]):1]-;!@[$([S;D2](-;!@[#0,#6])):2]>>[4*]-[*:1].[11*]-[*:2]",
    "[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$([S;D4]([#6,#0])(=O)(=O)):2]>>[5*]-[*:1].[12*]-[*:2]",
    "[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[5*]-[*:1].[14*]-[*:2]",
    "[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$([c;$(c(:c):c)]):2]>>[5*]-[*:1].[16*]-[*:2]",
    "[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[5*]-[*:1].[13*]-[*:2]",
    "[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[5*]-[*:1].[15*]-[*:2]",
    "[$([C;D3;!R](=O)-;!@[#0,#6,#7,#8]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[6*]-[*:1].[13*]-[*:2]",
    "[$([C;D3;!R](=O)-;!@[#0,#6,#7,#8]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[6*]-[*:1].[14*]-[*:2]",
    "[$([C;D3;!R](=O)-;!@[#0,#6,#7,#8]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[6*]-[*:1].[15*]-[*:2]",
    "[$([C;D3;!R](=O)-;!@[#0,#6,#7,#8]):1]-;!@[$([c;$(c(:c):c)]):2]>>[6*]-[*:1].[16*]-[*:2]",
    "[$([C;D2,D3]-[#6]):1]=;!@[$([C;D2,D3]-[#6]):2]>>[7*]-[*:1].[7*]-[*:2]",
    "[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):2]>>[8*]-[*:1].[9*]-[*:2]",
    "[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([N;R;$(N(@C(=O))@[C,N,O,S])]):2]>>[8*]-[*:1].[10*]-[*:2]",
    "[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[8*]-[*:1].[13*]-[*:2]",
    "[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[8*]-[*:1].[14*]-[*:2]",
    "[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[8*]-[*:1].[15*]-[*:2]",
    "[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([c;$(c(:c):c)]):2]>>[8*]-[*:1].[16*]-[*:2]",
    "[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[9*]-[*:1].[13*]-[*:2]",
    "[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[9*]-[*:1].[14*]-[*:2]",
    "[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[9*]-[*:1].[15*]-[*:2]",
    "[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):1]-;!@[$([c;$(c(:c):c)]):2]>>[9*]-[*:1].[16*]-[*:2]",
    "[$([N;R;$(N(@C(=O))@[C,N,O,S])]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[10*]-[*:1].[13*]-[*:2]",
    "[$([N;R;$(N(@C(=O))@[C,N,O,S])]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[10*]-[*:1].[14*]-[*:2]",
    "[$([N;R;$(N(@C(=O))@[C,N,O,S])]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[10*]-[*:1].[15*]-[*:2]",
    "[$([N;R;$(N(@C(=O))@[C,N,O,S])]):1]-;!@[$([c;$(c(:c):c)]):2]>>[10*]-[*:1].[16*]-[*:2]",
    "[$([S;D2](-;!@[#0,#6])):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[11*]-[*:1].[13*]-[*:2]",
    "[$([S;D2](-;!@[#0,#6])):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[11*]-[*:1].[14*]-[*:2]",
    "[$([S;D2](-;!@[#0,#6])):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[11*]-[*:1].[15*]-[*:2]",
    "[$([S;D2](-;!@[#0,#6])):1]-;!@[$([c;$(c(:c):c)]):2]>>[11*]-[*:1].[16*]-[*:2]",
    "[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[13*]-[*:1].[14*]-[*:2]",
    "[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[13*]-[*:1].[15*]-[*:2]",
    "[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):1]-;!@[$([c;$(c(:c):c)]):2]>>[13*]-[*:1].[16*]-[*:2]",
    "[$([c;$(c(:[c,n,o,s]):[n,o,s])]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[14*]-[*:1].[14*]-[*:2]",
    "[$([c;$(c(:[c,n,o,s]):[n,o,s])]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[14*]-[*:1].[15*]-[*:2]",
    "[$([c;$(c(:[c,n,o,s]):[n,o,s])]):1]-;!@[$([c;$(c(:c):c)]):2]>>[14*]-[*:1].[16*]-[*:2]",
    "[$([C;$(C(-;@C)-;@C)]):1]-;!@[$([c;$(c(:c):c)]):2]>>[15*]-[*:1].[16*]-[*:2]",
    "[$([c;$(c(:c):c)]):1]-;!@[$([c;$(c(:c):c)]):2]>>[16*]-[*:1].[16*]-[*:2]",
];

/// `getBRICSNetworkParams`
pub fn get_brics_network_params() -> Result<ScaffoldNetworkParams> {
    let mut res = ScaffoldNetworkParams::new(BRICS_DEFINITIONS)?;
    res.keep_only_first_fragment = false;
    Ok(res)
}

#[cfg(test)]
mod archive_tests {
    use super::*;

    /// Written by upstream's Boost.Serialization (from `testPickleScaffoldNetwork.py`).
    const UPSTREAM: &str = "22 serialization::archive 19 0 1 0 0 7 0 21 O=C1CCCC(Cc2ccccc2)N1 9 *c1ccccc1 15 **1:*:*:*:*:*:1 13 *C1CCCC(=O)N1 13 **1****(=*)*1 21 O=C1CCCC(Cc2ccccn2)N1 9 *c1ccccn1 7 0 1 1 2 2 2 1 1 7 0 1 1 2 2 2 1 1 0 0 7 0 0 0 0 1 1 1 2 2 0 3 1 3 4 2 5 6 1 6 2 2 5 3 1";

    #[test]
    fn upstream_archive_round_trips() {
        let net = ScaffoldNetwork::from_archive(UPSTREAM).unwrap();
        assert_eq!(net.nodes.len(), 7);
        assert_eq!(net.edges.len(), 7);
        assert_eq!(net.nodes[2], "**1:*:*:*:*:*:1");
        assert_eq!(net.edges[3].edge_type, EdgeType::Generic);
        assert_eq!(net.to_archive(), UPSTREAM);
    }

    #[test]
    fn empty_network_round_trips() {
        let net = ScaffoldNetwork::default();
        assert_eq!(
            ScaffoldNetwork::from_archive(&net.to_archive()).unwrap(),
            net
        );
    }
}
