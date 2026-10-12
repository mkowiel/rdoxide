//! Tautomer enumeration, scoring and canonicalization.
//!
//! Ported from RDKit `Code/GraphMol/MolStandardize/Tautomer.{h,cpp}`. The
//! enumeration rules follow M. Sitzmann et al., "Tautomerism in Large
//! Databases.", JCAMD 24:521 (2010).

use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};
use std::ops::Bound;
use std::sync::{Arc, OnceLock};

use super::CleanupParameters;
use super::catalogs::TautomerCatalogParams;
use crate::error::{RDKitError, Result, SanitizeError};
use crate::graphmol::aromaticity::{AromaticityModel, set_aromaticity};
use crate::graphmol::atom::{ChiralType, HybridizationType};
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::conjug_hybrid::{set_conjugation, set_hybridization};
use crate::graphmol::find_rings::{SymmetrizeSSSRAlgorithm, fast_find_rings, symmetrize_sssr};
use crate::graphmol::kekulize::kekulize;
use crate::graphmol::mol_ops::adjust_hs;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::graphmol::substruct::{
    MatchVect, SubstructMatchParameters, substruct_match, substruct_match_count,
};
use crate::rdgeneral::common_properties;
use crate::rdgeneral::dict::RDProps;
use crate::{rd_log_error, rd_log_warning};

/// `tautomerScoringVersion`
pub const TAUTOMER_SCORING_VERSION: &str = "1.0.0";

fn smiles(mol: &mut ROMol) -> Result<String> {
    mol_to_smiles(mol, &SmilesWriteParams::default())
}

// ------------------------------------------------------------------- scoring

/// `TautomerScoringFunctions::SubstructTerm`: `score` per match of `smarts`.
#[derive(Debug, Clone)]
pub struct SubstructTerm {
    pub name: String,
    pub smarts: String,
    pub score: i32,
    /// empty when the SMARTS could not be parsed
    pub matcher: ROMol,
    /// atomic numbers that must be present for the term to match
    pub required_elements: Vec<u32>,
    /// bond-order agnostic pre-screen pattern
    pub connectivity_smarts: String,
    pub connectivity_matcher: ROMol,
}

impl PartialEq for SubstructTerm {
    fn eq(&self, rhs: &Self) -> bool {
        self.name == rhs.name && self.smarts == rhs.smarts && self.score == rhs.score
    }
}

fn smarts_or_empty(smarts: &str) -> ROMol {
    mol_from_smarts(smarts, &SmartsParserParams::default())
        .ok()
        .flatten()
        .unwrap_or_default()
}

impl SubstructTerm {
    pub fn new(
        name: &str,
        smarts: &str,
        score: i32,
        required_elements: Vec<u32>,
        connectivity_smarts: &str,
    ) -> Self {
        let connectivity_matcher = if connectivity_smarts.is_empty() {
            ROMol::default()
        } else {
            smarts_or_empty(connectivity_smarts)
        };
        SubstructTerm {
            name: name.to_string(),
            smarts: smarts.to_string(),
            score,
            matcher: smarts_or_empty(smarts),
            required_elements,
            connectivity_smarts: connectivity_smarts.to_string(),
            connectivity_matcher,
        }
    }
}

/// `getDefaultTautomerScoreSubstructs`
pub fn default_tautomer_score_substructs() -> &'static [SubstructTerm] {
    static TERMS: OnceLock<Vec<SubstructTerm>> = OnceLock::new();
    TERMS.get_or_init(|| {
        vec![
            SubstructTerm::new(
                "benzoquinone",
                "[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]",
                25,
                vec![6],
                "[#6]1(~[#6]~[#6]~[#6](~[#6]~[#6]~1)~[N,S,O])~[N,S,O]",
            ),
            SubstructTerm::new("oxim", "[#6]=[N][OH]", 4, vec![6, 7, 8], "[#6]~[#7]~[#8]"),
            SubstructTerm::new("C=O", "[#6]=,:[#8]", 2, vec![6, 8], "[#6]~[#8]"),
            SubstructTerm::new("N=O", "[#7]=,:[#8]", 2, vec![7, 8], "[#7]~[#8]"),
            SubstructTerm::new("P=O", "[#15]=,:[#8]", 2, vec![15, 8], "[#15]~[#8]"),
            SubstructTerm::new("C=hetero", "[C]=[!#1;!#6]", 1, vec![6], "[C]~[!#1;!#6]"),
            SubstructTerm::new(
                "C(=hetero)-hetero",
                "[C](=[!#1;!#6])[!#1;!#6]",
                2,
                vec![6],
                "[C](~[!#1;!#6])~[!#1;!#6]",
            ),
            SubstructTerm::new(
                "aromatic C = exocyclic N",
                "[c]=!@[N]",
                -1,
                vec![6, 7],
                "[c]~[N]",
            ),
            SubstructTerm::new("methyl", "[CX4H3]", 1, vec![6], ""),
            SubstructTerm::new(
                "guanidine terminal=N",
                "[#7]C(=[NR0])[#7H0]",
                1,
                vec![6, 7],
                "[#7]~[#6]~[#7]",
            ),
            SubstructTerm::new(
                "guanidine endocyclic=N",
                "[#7;R][#6;R]([N])=[#7;R]",
                2,
                vec![6, 7],
                "[#7]~[#6](~[N])~[#7]",
            ),
            SubstructTerm::new(
                "aci-nitro",
                "[#6]=[N+]([O-])[OH]",
                -4,
                vec![6, 7, 8],
                "[#6]~[#7](~[#8])~[#8]",
            ),
        ]
    })
}

/// `mol` with symmetrized SSSR ring information (upstream computes it on
/// the const molecule).
fn with_symm_rings(mol: &ROMol) -> Result<Cow<'_, ROMol>> {
    if mol.ring_info().is_symm_sssr() {
        return Ok(Cow::Borrowed(mol));
    }
    let mut m = mol.clone();
    symmetrize_sssr(&mut m, SymmetrizeSSSRAlgorithm::Default, true, false, false)?;
    Ok(Cow::Owned(m))
}

fn score_rings_impl(mol: &ROMol) -> i32 {
    let nb = mol.num_bonds() as usize;
    let mut is_arom = vec![false; nb];
    let mut both_carbon = vec![false; nb];
    for (i, b) in mol.bonds().iter().enumerate() {
        if b.is_aromatic() {
            is_arom[i] = true;
            if mol.atom(b.begin_atom_idx()).atomic_num() == 6
                && mol.atom(b.end_atom_idx()).atomic_num() == 6
            {
                both_carbon[i] = true;
            }
        }
    }
    let mut score = 0;
    for bring in mol.ring_info().bond_rings() {
        let mut all_c = true;
        let mut all_aromatic = true;
        for &bidx in bring {
            if !is_arom[bidx as usize] {
                all_aromatic = false;
                break;
            }
            if !both_carbon[bidx as usize] {
                all_c = false;
            }
        }
        if all_aromatic {
            score += 100;
            if all_c {
                score += 150;
            }
        }
    }
    score
}

/// `scoreRings`: 100 per aromatic ring, 250 per all-carbon aromatic ring.
pub fn score_rings(mol: &ROMol) -> Result<i32> {
    Ok(score_rings_impl(&*with_symm_rings(mol)?))
}

/// `scoreSubstructs`
pub fn score_substructs(mol: &ROMol, terms: &[SubstructTerm]) -> Result<i32> {
    let mut score = 0;
    let ps = SubstructMatchParameters::default();
    for term in terms {
        if term.matcher.num_atoms() == 0 {
            rd_log_error!(" matcher for term {} is invalid, ignoring it.", term.name);
            continue;
        }
        score += substruct_match_count(mol, &term.matcher, &ps)? as i32 * term.score;
    }
    Ok(score)
}

/// `getRelevantSubstructTermIndices`: the terms whose required elements
/// and connectivity are present in `mol` (tautomerization only moves Hs and
/// changes bond orders).
pub fn relevant_substruct_term_indices(mol: &ROMol, terms: &[SubstructTerm]) -> Result<Vec<usize>> {
    let present: BTreeSet<u32> = mol.atoms().iter().map(|a| a.atomic_num()).collect();
    let ps = SubstructMatchParameters {
        max_matches: 1,
        ..Default::default()
    };
    let mut res = Vec::with_capacity(terms.len());
    for (i, term) in terms.iter().enumerate() {
        if !term.required_elements.iter().all(|e| present.contains(e)) {
            continue;
        }
        if term.connectivity_matcher.num_atoms() > 0
            && substruct_match(mol, &term.connectivity_matcher, &ps)?.is_empty()
        {
            continue;
        }
        res.push(i);
    }
    Ok(res)
}

/// `[#elem1]=,:[#elem2]`
fn count_double_or_aromatic_bonds(mol: &ROMol, elem1: u32, elem2: u32) -> u32 {
    mol.bonds()
        .iter()
        .filter(|b| {
            matches!(b.bond_type(), BondType::Double | BondType::Aromatic) && {
                let a1 = mol.atom(b.begin_atom_idx()).atomic_num();
                let a2 = mol.atom(b.end_atom_idx()).atomic_num();
                (a1 == elem1 && a2 == elem2) || (a1 == elem2 && a2 == elem1)
            }
        })
        .count() as u32
}

/// `[CX4H3]`
fn count_methyls(mol: &ROMol) -> Result<u32> {
    let mut count = 0;
    for i in 0..mol.num_atoms() {
        if mol.atom(i).atomic_num() == 6
            && mol.total_degree(i)? == 4
            && mol.total_num_hs(i, false)? == 3
        {
            count += 1;
        }
    }
    Ok(count)
}

/// `[C]=[!#1;!#6]`
fn count_carbon_double_hetero(mol: &ROMol) -> u32 {
    let is_c_to_het = |c: u32, het: u32| {
        let (c, het) = (mol.atom(c), mol.atom(het));
        c.atomic_num() == 6 && !c.is_aromatic() && het.atomic_num() != 1 && het.atomic_num() != 6
    };
    mol.bonds()
        .iter()
        .filter(|b| {
            b.bond_type() == BondType::Double
                && (is_c_to_het(b.begin_atom_idx(), b.end_atom_idx())
                    || is_c_to_het(b.end_atom_idx(), b.begin_atom_idx()))
        })
        .count() as u32
}

/// `[c]=!@[N]`
fn count_aromatic_carbon_exocyclic_n(mol: &ROMol) -> Result<u32> {
    let is_c_to_n = |c: u32, n: u32| {
        let (c, n) = (mol.atom(c), mol.atom(n));
        c.atomic_num() == 6 && c.is_aromatic() && n.atomic_num() == 7
    };
    let mut count = 0;
    for (i, b) in mol.bonds().iter().enumerate() {
        if b.bond_type() != BondType::Double || mol.ring_info().num_bond_rings(i as u32)? > 0 {
            continue;
        }
        if is_c_to_n(b.begin_atom_idx(), b.end_atom_idx())
            || is_c_to_n(b.end_atom_idx(), b.begin_atom_idx())
        {
            count += 1;
        }
    }
    Ok(count)
}

/// `scoreSubstructsFiltered`: the terms at `relevant` (indices into the
/// default terms get the dedicated matchers).
pub fn score_substructs_filtered(
    mol: &ROMol,
    terms: &[SubstructTerm],
    relevant: &[usize],
) -> Result<i32> {
    let ps = SubstructMatchParameters::default();
    let mut score = 0;
    for &idx in relevant {
        let term = &terms[idx];
        if term.matcher.num_atoms() == 0 {
            continue;
        }
        let n = match idx {
            2 => count_double_or_aromatic_bonds(mol, 6, 8),
            3 => count_double_or_aromatic_bonds(mol, 7, 8),
            4 => count_double_or_aromatic_bonds(mol, 15, 8),
            5 => count_carbon_double_hetero(mol),
            7 => count_aromatic_carbon_exocyclic_n(mol)?,
            8 => count_methyls(mol)?,
            _ => substruct_match_count(mol, &term.matcher, &ps)?,
        };
        score += n as i32 * term.score;
    }
    Ok(score)
}

/// `scoreHeteroHs`: -1 per H on P, S, Se or Te.
pub fn score_hetero_hs(mol: &ROMol) -> Result<i32> {
    let mut score = 0;
    for i in 0..mol.num_atoms() {
        if matches!(mol.atom(i).atomic_num(), 15 | 16 | 34 | 52) {
            score -= mol.total_num_hs(i, false)? as i32;
        }
    }
    Ok(score)
}

/// `scoreTautomer`
pub fn score_tautomer(mol: &ROMol) -> Result<i32> {
    let mol = with_symm_rings(mol)?;
    Ok(score_rings_impl(&mol)
        + score_substructs(&mol, default_tautomer_score_substructs())?
        + score_hetero_hs(&mol)?)
}

/// `makeOptimizedScorer`: `scoreTautomer` with the substructure terms
/// pre-filtered for `mol`.
pub fn make_optimized_scorer(mol: &ROMol) -> Result<impl Fn(&ROMol) -> Result<i32> + Sync + use<>> {
    let relevant = relevant_substruct_term_indices(mol, default_tautomer_score_substructs())?;
    Ok(move |taut: &ROMol| {
        let taut = with_symm_rings(taut)?;
        Ok(score_rings_impl(&taut)
            + score_substructs_filtered(&taut, default_tautomer_score_substructs(), &relevant)?
            + score_hetero_hs(&taut)?)
    })
}

/// `detail::countSpecifiedStereo`: tie breaker of the canonical tautomer.
fn count_specified_stereo(mol: &ROMol) -> u32 {
    let atoms = mol
        .atoms()
        .iter()
        .filter(|a| a.chiral_tag() != ChiralType::Unspecified)
        .count();
    let bonds = mol
        .bonds()
        .iter()
        .filter(|b| (b.stereo() as u8) > (BondStereo::StereoAny as u8))
        .count();
    (atoms + bonds) as u32
}

// ---------------------------------------------------------------- enumerator

/// `TautomerEnumeratorStatus`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TautomerEnumeratorStatus {
    #[default]
    Completed = 0,
    MaxTautomersReached = 1,
    MaxTransformsReached = 2,
    Canceled = 3,
}

impl TautomerEnumeratorStatus {
    fn message(self) -> &'static str {
        match self {
            Self::Completed => "completed",
            Self::MaxTautomersReached => "max tautomers reached",
            Self::MaxTransformsReached => "max transforms reached",
            Self::Canceled => "canceled",
        }
    }
}

/// `Tautomer`: a tautomer and (lazily) its kekulized form.
#[derive(Debug, Clone)]
pub struct Tautomer {
    pub tautomer: Arc<ROMol>,
    kekulized: OnceLock<Arc<ROMol>>,
    num_modified_atoms: usize,
    num_modified_bonds: usize,
    done: bool,
}

impl Tautomer {
    fn new(tautomer: ROMol, num_modified_atoms: usize, num_modified_bonds: usize) -> Self {
        Tautomer {
            tautomer: Arc::new(tautomer),
            kekulized: OnceLock::new(),
            num_modified_atoms,
            num_modified_bonds,
            done: false,
        }
    }

    /// `getKekulized`: `Kekulize(t, markAtomsBonds = false, canonical = true)`.
    pub fn kekulized(&self) -> Result<Arc<ROMol>> {
        if let Some(k) = self.kekulized.get() {
            return Ok(k.clone());
        }
        let mut k = (*self.tautomer).clone();
        kekulize(&mut k, false, true, 100)?;
        let k = Arc::new(k);
        let _ = self.kekulized.set(k.clone());
        Ok(k)
    }
}

/// `TautomerEnumeratorResult`: tautomers keyed (and ordered) by SMILES.
#[derive(Debug, Clone, Default)]
pub struct TautomerEnumeratorResult {
    pub tautomers: BTreeMap<String, Tautomer>,
    pub status: TautomerEnumeratorStatus,
    pub modified_atoms: Vec<bool>,
    pub modified_bonds: Vec<bool>,
}

impl TautomerEnumeratorResult {
    pub fn len(&self) -> usize {
        self.tautomers.len()
    }

    pub fn is_empty(&self) -> bool {
        self.tautomers.is_empty()
    }

    /// `tautomers()`
    pub fn tautomer_mols(&self) -> Vec<Arc<ROMol>> {
        self.tautomers
            .values()
            .map(|t| t.tautomer.clone())
            .collect()
    }

    /// `smiles()`
    pub fn smiles(&self) -> Vec<String> {
        self.tautomers.keys().cloned().collect()
    }

    /// `at(pos)`
    pub fn at(&self, pos: usize) -> Option<&Arc<ROMol>> {
        self.tautomers.values().nth(pos).map(|t| &t.tautomer)
    }

    /// Indices of the atoms modified by the transforms.
    pub fn modified_atom_indices(&self) -> Vec<u32> {
        indices(&self.modified_atoms)
    }

    /// Indices of the bonds modified by the transforms.
    pub fn modified_bond_indices(&self) -> Vec<u32> {
        indices(&self.modified_bonds)
    }
}

fn indices(bits: &[bool]) -> Vec<u32> {
    bits.iter()
        .enumerate()
        .filter(|(_, b)| **b)
        .map(|(i, _)| i as u32)
        .collect()
}

/// `TautomerEnumeratorCallback`: return `false` to stop the enumeration.
pub type TautomerEnumeratorCallback =
    dyn Fn(&ROMol, &TautomerEnumeratorResult) -> Result<bool> + Send + Sync;

/// A tautomer scoring function.
pub type ScoreFunc<'a> = &'a (dyn Fn(&ROMol) -> Result<i32> + Sync);

/// `TautomerEnumerator`
#[derive(Clone)]
pub struct TautomerEnumerator {
    params: Arc<TautomerCatalogParams>,
    callback: Option<Arc<TautomerEnumeratorCallback>>,
    pub max_tautomers: u32,
    pub max_transforms: u32,
    pub remove_sp3_stereo: bool,
    pub remove_bond_stereo: bool,
    pub remove_isotopic_hs: bool,
    pub reassign_stereo: bool,
}

impl std::fmt::Debug for TautomerEnumerator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TautomerEnumerator")
            .field("num_transforms", &self.params.transforms().len())
            .field("has_callback", &self.callback.is_some())
            .field("max_tautomers", &self.max_tautomers)
            .field("max_transforms", &self.max_transforms)
            .finish()
    }
}

impl Default for TautomerEnumerator {
    /// `TautomerEnumerator()`: the default transforms and parameters.
    fn default() -> Self {
        Self::with_params(TautomerCatalogParams::defaults().clone())
    }
}

/// The bond stereo of a cleared tautomeric double bond: `STEREOANY` unless
/// the bond is in a ring or joins two ring atoms.
fn cleared_tautomer_bond_stereo(mol: &ROMol, bond: u32) -> Result<BondStereo> {
    let b = mol.bond(bond);
    if b.bond_type() != BondType::Double {
        return Ok(BondStereo::StereoNone);
    }
    let ri = mol.ring_info();
    let at_ring = ri.num_bond_rings(bond)? != 0
        || (ri.num_atom_rings(b.begin_atom_idx())? != 0
            && ri.num_atom_rings(b.end_atom_idx())? != 0);
    Ok(if at_ring {
        BondStereo::StereoNone
    } else {
        BondStereo::StereoAny
    })
}

fn has_preservable_double_bond_stereo(mol: &ROMol, bond: u32) -> bool {
    let b = mol.bond(bond);
    let stereo = b.stereo() as u8;
    b.bond_type() == BondType::Double
        && (stereo == BondStereo::StereoAny as u8
            || ((BondStereo::StereoZ as u8..=BondStereo::StereoTrans as u8).contains(&stereo)
                && b.stereo_atoms().len() == 2))
}

fn ensure_fast_rings(mol: &mut ROMol) -> Result<()> {
    if !mol.ring_info().is_find_fast_or_better() {
        fast_find_rings(mol)?;
    }
    Ok(())
}

impl TautomerEnumerator {
    pub fn with_params(params: TautomerCatalogParams) -> Self {
        TautomerEnumerator {
            params: Arc::new(params),
            callback: None,
            max_tautomers: 1000,
            max_transforms: 1000,
            remove_sp3_stereo: true,
            remove_bond_stereo: true,
            remove_isotopic_hs: true,
            reassign_stereo: true,
        }
    }

    /// `TautomerEnumerator(const CleanupParameters &)`
    pub fn from_params(params: &CleanupParameters) -> Result<Self> {
        let catalog = if !params.tautomer_transform_data.is_empty() {
            TautomerCatalogParams::from_data(&params.tautomer_transform_data)?
        } else if params.tautomer_transforms.is_empty() {
            TautomerCatalogParams::defaults().clone()
        } else {
            TautomerCatalogParams::from_file(&params.tautomer_transforms)?
        };
        let mut te = Self::with_params(catalog);
        te.max_tautomers = params.max_tautomers as u32;
        te.max_transforms = params.max_transforms as u32;
        te.remove_sp3_stereo = params.tautomer_remove_sp3_stereo;
        te.remove_bond_stereo = params.tautomer_remove_bond_stereo;
        te.remove_isotopic_hs = params.tautomer_remove_isotopic_hs;
        te.reassign_stereo = params.tautomer_reassign_stereo;
        Ok(te)
    }

    /// `getV1TautomerEnumerator`
    pub fn v1() -> Self {
        Self::with_params(TautomerCatalogParams::defaults_v1().clone())
    }

    pub fn params(&self) -> &TautomerCatalogParams {
        &self.params
    }

    pub fn set_callback(&mut self, callback: Option<Arc<TautomerEnumeratorCallback>>) {
        self.callback = callback;
    }

    pub fn callback(&self) -> Option<&Arc<TautomerEnumeratorCallback>> {
        self.callback.as_ref()
    }

    /// `setTautomerStereoAndIsoHs`: resets the stereo of the atoms and bonds
    /// involved in tautomerism; whether `taut` changed.
    fn set_tautomer_stereo_and_iso_hs(
        &self,
        mol: &ROMol,
        taut: &mut ROMol,
        modified_atoms: &[bool],
        modified_bonds: &[bool],
    ) -> Result<bool> {
        let mut modified = false;
        for atom_idx in indices(modified_atoms) {
            let atom = mol.atom(atom_idx);
            let total_hs = taut.total_num_hs(atom_idx, false)?;
            let taut_atom = taut.atom_mut(atom_idx);
            // clear the chiral tag of sp2 atoms (and sp3 ones if requested)
            if taut_atom.hybridization() == HybridizationType::SP2 || self.remove_sp3_stereo {
                modified |= taut_atom.chiral_tag() != ChiralType::Unspecified;
                taut_atom.set_chiral_tag(ChiralType::Unspecified);
                if taut_atom.has_prop(common_properties::_CIPCode) {
                    taut_atom.clear_prop(common_properties::_CIPCode);
                }
            } else {
                modified |= taut_atom.chiral_tag() != atom.chiral_tag();
                taut_atom.set_chiral_tag(atom.chiral_tag());
                if let Some(cip) = atom.get_prop(common_properties::_CIPCode) {
                    taut_atom
                        .dict_mut()
                        .set_val(common_properties::_CIPCode, cip.to_rd_string());
                }
            }
            // remove isotopic Hs
            if taut_atom.has_prop(common_properties::_isotopicHs)
                && (self.remove_isotopic_hs || total_hs == 0)
            {
                taut_atom.clear_prop(common_properties::_isotopicHs);
            }
        }
        // the stereo of bonds that are part of a tautomeric path
        for bond_idx in indices(modified_bonds) {
            let bond = mol.bond(bond_idx);
            let mut bonds_to_clear_dirs = Vec::new();
            if bond.bond_type() == BondType::Double
                && (bond.stereo() as u8) > (BondStereo::StereoAny as u8)
            {
                // bonds with a direction around the begin and end atoms
                for atom in [bond.begin_atom_idx(), bond.end_atom_idx()] {
                    for obnd in mol.atom_bonds(atom) {
                        if matches!(obnd.bond_dir(), BondDir::EndDownRight | BondDir::EndUpRight) {
                            bonds_to_clear_dirs.push(obnd.idx());
                        }
                    }
                }
            }
            let taut_is_double = taut.bond(bond_idx).bond_type() == BondType::Double;
            if !taut_is_double
                || self.remove_bond_stereo
                || !has_preservable_double_bond_stereo(mol, bond_idx)
            {
                // STEREOANY (not STEREONONE) for acyclic double bonds keeps
                // downstream tools (InChI) from inferring E/Z from 2D
                // coordinates after the bond orders changed
                if taut_is_double {
                    ensure_fast_rings(taut)?;
                }
                let target = cleared_tautomer_bond_stereo(taut, bond_idx)?;
                let tb = taut.bond_mut(bond_idx);
                modified |= tb.stereo() != target;
                tb.set_stereo(target)?;
                tb.stereo_atoms_mut().clear();
                for &bi in &bonds_to_clear_dirs {
                    taut.bond_mut(bi).set_bond_dir(BondDir::None);
                }
            } else {
                let sa = bond.stereo_atoms().to_vec();
                let tb = taut.bond(bond_idx);
                modified |= tb.stereo() != bond.stereo() || sa.len() != tb.stereo_atoms().len();
                if sa.len() == 2 {
                    taut.set_bond_stereo_atoms(bond_idx, sa[0] as u32, sa[1] as u32)?;
                }
                taut.bond_mut(bond_idx).set_stereo(bond.stereo())?;
                for &bi in &bonds_to_clear_dirs {
                    taut.bond_mut(bi).set_bond_dir(mol.bond(bi).bond_dir());
                }
            }
        }
        if self.reassign_stereo {
            assign_stereochemistry(taut, true, true, false)?;
            // assignStereochemistry() may overwrite the explicit "undefined"
            // bond stereo set above: re-apply it
            if self.remove_bond_stereo {
                ensure_fast_rings(taut)?;
                for bond_idx in indices(modified_bonds) {
                    let target = if taut.bond(bond_idx).bond_type() != BondType::Double {
                        BondStereo::StereoNone
                    } else {
                        cleared_tautomer_bond_stereo(taut, bond_idx)?
                    };
                    let tb = taut.bond_mut(bond_idx);
                    tb.set_stereo(target)?;
                    tb.stereo_atoms_mut().clear();
                }
            }
        } else {
            taut.dict_mut()
                .set_val(common_properties::_StereochemDone, 1i32);
        }
        Ok(modified)
    }

    /// `enumerate`: all the tautomers of `mol`. Atoms with the `_protected`
    /// property are not considered.
    pub fn enumerate(&self, mol: &ROMol) -> Result<TautomerEnumeratorResult> {
        let transforms = self.params.transforms();
        let protected: Vec<bool> = mol
            .atoms()
            .iter()
            .map(|a| a.has_prop(common_properties::_protected))
            .collect();

        let mut taut = mol.clone();
        let smi = smiles(&mut taut)?;
        if taut.needs_update_property_cache() {
            taut.update_property_cache(false)?;
        }
        if !taut.ring_info().is_symm_sssr() {
            symmetrize_sssr(
                &mut taut,
                SymmetrizeSSSRAlgorithm::Default,
                true,
                false,
                false,
            )?;
        }
        let mut res = TautomerEnumeratorResult {
            modified_atoms: vec![false; mol.num_atoms() as usize],
            modified_bonds: vec![false; mol.num_bonds() as usize],
            ..Default::default()
        };
        res.tautomers.insert(smi, Tautomer::new(taut, 0, 0));
        let mut num_modified_atoms = 0usize;
        let mut num_modified_bonds = 0usize;
        let mut completed = false;
        let mut bail_out = false;
        let mut n_transforms = 0u32;

        while !completed && !bail_out {
            // std::map iteration: tautomers inserted after the current key
            // are visited in the same pass
            let mut cursor = res.tautomers.keys().next().cloned();
            while let Some(key) = cursor {
                if !res.tautomers[&key].done {
                    for transform in transforms {
                        if bail_out {
                            break;
                        }
                        let kmol = res.tautomers[&key].kekulized()?;
                        let matches: Vec<MatchVect> = substruct_match(
                            &kmol,
                            &transform.mol,
                            &SubstructMatchParameters::default(),
                        )?
                        .into_iter()
                        .filter(|m| !m.iter().any(|&(_, a)| protected[a as usize]))
                        .collect();
                        if matches.is_empty() {
                            continue;
                        }
                        n_transforms += 1;
                        for m in &matches {
                            if n_transforms >= self.max_transforms {
                                res.status = TautomerEnumeratorStatus::MaxTransformsReached;
                                bail_out = true;
                            } else if res.tautomers.len() >= self.max_tautomers as usize {
                                res.status = TautomerEnumeratorStatus::MaxTautomersReached;
                                bail_out = true;
                            } else if let Some(cb) = &self.callback {
                                if !cb(mol, &res)? {
                                    res.status = TautomerEnumeratorStatus::Canceled;
                                    bail_out = true;
                                }
                            }
                            if bail_out {
                                break;
                            }
                            // modify a copy of the kekulized form
                            let mut product = ROMol::copy_from(&kmol, true, -1);
                            // move a hydrogen from the first matched atom to the last
                            let first_idx = m[0].1 as u32;
                            let last_idx = m[m.len() - 1].1 as u32;
                            for idx in [first_idx, last_idx] {
                                if !res.modified_atoms[idx as usize] {
                                    res.modified_atoms[idx as usize] = true;
                                    num_modified_atoms += 1;
                                }
                            }
                            let first_hs = product.total_num_hs(first_idx, false)? as i32;
                            product
                                .atom_mut(first_idx)
                                .set_num_explicit_hs((first_hs - 1).max(0) as u32);
                            let last_hs = product.total_num_hs(last_idx, false)?;
                            product.atom_mut(last_idx).set_num_explicit_hs(last_hs + 1);
                            product.atom_mut(first_idx).set_no_implicit(true);
                            product.atom_mut(last_idx).set_no_implicit(true);
                            // adjust the bond orders
                            let mut bi = 0;
                            for tbond in transform.mol.bonds() {
                                let b1 = m[tbond.begin_atom_idx() as usize].1 as u32;
                                let b2 = m[tbond.end_atom_idx() as usize].1 as u32;
                                let bidx = product.bond_idx_between(b1, b2).ok_or_else(|| {
                                    RDKitError::invariant("required bond not found")
                                })?;
                                let bond = product.bond_mut(bidx);
                                if !transform.bond_types.is_empty() {
                                    bond.set_bond_type(transform.bond_types[bi]);
                                    bi += 1;
                                } else {
                                    let bt = bond.bond_type();
                                    if bt == BondType::Single {
                                        bond.set_bond_type(BondType::Double);
                                    }
                                    if bt == BondType::Double {
                                        bond.set_bond_type(BondType::Single);
                                    }
                                }
                                if !res.modified_bonds[bidx as usize] {
                                    res.modified_bonds[bidx as usize] = true;
                                    num_modified_bonds += 1;
                                }
                            }
                            if !transform.charges.is_empty() {
                                for (ci, &(_, aidx)) in m.iter().enumerate() {
                                    let atom = product.atom_mut(aidx as u32);
                                    atom.set_formal_charge(
                                        atom.formal_charge() + transform.charges[ci],
                                    );
                                }
                            }
                            // the graph (and the rings) did not change: only
                            // the sanitization steps that are needed
                            match sanitize_product(&mut product) {
                                Ok(()) => {}
                                Err(RDKitError::Sanitize(SanitizeError::Kekulize { .. })) => {
                                    continue;
                                }
                                Err(e) => return Err(e),
                            }
                            self.set_tautomer_stereo_and_iso_hs(
                                mol,
                                &mut product,
                                &res.modified_atoms,
                                &res.modified_bonds,
                            )?;
                            let tsmiles = smiles(&mut product)?;
                            if res.tautomers.contains_key(&tsmiles) {
                                continue;
                            }
                            // sanitization may modify bonds too, e.g. Cc1nc2ccccc2[nH]1
                            for (i, (mb, pb)) in mol.bonds().iter().zip(product.bonds()).enumerate()
                            {
                                if mb.bond_type() != pb.bond_type() && !res.modified_bonds[i] {
                                    res.modified_bonds[i] = true;
                                    num_modified_bonds += 1;
                                }
                            }
                            res.tautomers.insert(
                                tsmiles,
                                Tautomer::new(product, num_modified_atoms, num_modified_bonds),
                            );
                        }
                    }
                    if let Some(t) = res.tautomers.get_mut(&key) {
                        t.done = true;
                    }
                }
                cursor = next_key(&res.tautomers, &key);
            }

            completed = true;
            let max_atoms = num_modified_atoms;
            let max_bonds = num_modified_bonds;
            let mut cursor = res.tautomers.keys().next().cloned();
            while let Some(key) = cursor {
                let t = &res.tautomers[&key];
                if !t.done {
                    completed = false;
                }
                let stale = t.num_modified_atoms < max_atoms || t.num_modified_bonds < max_bonds;
                let mut next = None;
                if stale {
                    let mut t = res.tautomers.remove(&key).expect("present");
                    let changed = self.set_tautomer_stereo_and_iso_hs(
                        mol,
                        Arc::make_mut(&mut t.tautomer),
                        &res.modified_atoms,
                        &res.modified_bonds,
                    )?;
                    if changed {
                        // erase + insert under the new SMILES; a duplicate is
                        // dropped and iteration continues after the old key,
                        // otherwise at the new one
                        t.num_modified_atoms = max_atoms;
                        t.num_modified_bonds = max_bonds;
                        let new_key = smiles(Arc::make_mut(&mut t.tautomer))?;
                        next = Some(next_key(&res.tautomers, &key));
                        if !res.tautomers.contains_key(&new_key) {
                            res.tautomers.insert(new_key.clone(), t);
                            next = Some(Some(new_key));
                        }
                    } else {
                        res.tautomers.insert(key.clone(), t);
                    }
                }
                cursor = match next {
                    Some(n) => n,
                    None => next_key(&res.tautomers, &key),
                };
            }
            if bail_out
                && res.tautomers.len() < self.max_tautomers as usize
                && res.status == TautomerEnumeratorStatus::MaxTautomersReached
            {
                res.status = TautomerEnumeratorStatus::Completed;
                bail_out = false;
            }
        }
        if !completed {
            rd_log_warning!(
                "Tautomer enumeration stopped at {} tautomers: {}",
                res.tautomers.len(),
                res.status.message()
            );
        }
        Ok(res)
    }

    /// `pickCanonical(const TautomerEnumeratorResult &, scoreFunc)`
    pub fn pick_canonical(
        &self,
        res: &TautomerEnumeratorResult,
        score: Option<ScoreFunc<'_>>,
    ) -> Result<ROMol> {
        let score: ScoreFunc<'_> = score.unwrap_or(&score_tautomer);
        let best: &ROMol = if res.tautomers.len() == 1 {
            &res.tautomers
                .values()
                .next()
                .expect("one tautomer")
                .tautomer
        } else {
            let mut best_score = i32::MIN;
            let mut best_stereo = 0u32;
            let mut best_smiles = "";
            let mut best_mol: Option<&ROMol> = None;
            for (smi, t) in &res.tautomers {
                let s = score(&t.tautomer)?;
                if s < best_score {
                    continue;
                }
                let n_stereo = count_specified_stereo(&t.tautomer);
                let better = s > best_score
                    || n_stereo > best_stereo
                    || (n_stereo == best_stereo && smi.as_str() < best_smiles);
                if better {
                    best_score = s;
                    best_stereo = n_stereo;
                    best_smiles = smi;
                    best_mol = Some(&t.tautomer);
                }
            }
            best_mol.ok_or_else(|| RDKitError::invariant("no tautomers"))?
        };
        let mut res = best.clone();
        assign_stereochemistry(&mut res, true, true, false)?;
        Ok(res)
    }

    /// `pickCanonical(iterable, scoreFunc)`: the canonical tautomer among
    /// `tautomers`.
    pub fn pick_canonical_mols(
        &self,
        tautomers: &[ROMol],
        score: Option<ScoreFunc<'_>>,
    ) -> Result<ROMol> {
        let score: ScoreFunc<'_> = score.unwrap_or(&score_tautomer);
        let best: &ROMol = if tautomers.len() == 1 {
            &tautomers[0]
        } else {
            let mut best_score = i32::MIN;
            let mut best_stereo = 0u32;
            let mut best_smiles = String::new();
            let mut best_mol: Option<&ROMol> = None;
            for t in tautomers {
                let s = score(t)?;
                if s < best_score {
                    continue;
                }
                let n_stereo = count_specified_stereo(t);
                let mut smi = None;
                let better = if s > best_score || n_stereo > best_stereo {
                    true
                } else if n_stereo == best_stereo {
                    let ts = smiles(&mut t.clone())?;
                    let better = ts < best_smiles;
                    smi = Some(ts);
                    better
                } else {
                    false
                };
                if better {
                    best_smiles = match smi {
                        Some(s) => s,
                        None => smiles(&mut t.clone())?,
                    };
                    best_score = s;
                    best_stereo = n_stereo;
                    best_mol = Some(t);
                }
            }
            best_mol.ok_or_else(|| RDKitError::invariant("no tautomers"))?
        };
        let mut res = best.clone();
        assign_stereochemistry(&mut res, true, true, false)?;
        Ok(res)
    }

    /// `canonicalize`: the canonical tautomer of `mol`.
    pub fn canonicalize(&self, mol: &ROMol, score: Option<ScoreFunc<'_>>) -> Result<ROMol> {
        let mut this = self.clone();
        this.reassign_stereo = false;
        let res = this.enumerate(mol)?;
        if res.is_empty() {
            rd_log_warning!("no tautomers found, returning input molecule");
            return Ok(mol.clone());
        }
        let optimized;
        let score: ScoreFunc<'_> = match score {
            Some(s) => s,
            None => {
                optimized = make_optimized_scorer(mol)?;
                &optimized
            }
        };
        let mut canonical = self.pick_canonical(&res, Some(score))?;
        // the enumeration's quick copies have no properties or conformers
        canonical.update_props(mol.dict(), false);
        for conf in mol.conformers() {
            canonical.add_conformer(conf.clone(), true)?;
        }
        Ok(canonical)
    }

    /// `canonicalizeInPlace`
    pub fn canonicalize_in_place(
        &self,
        mol: &mut RWMol,
        score: Option<ScoreFunc<'_>>,
    ) -> Result<()> {
        let mut this = self.clone();
        this.reassign_stereo = false;
        let res = this.enumerate(mol)?;
        if res.is_empty() {
            rd_log_warning!("no tautomers found, molecule unchanged");
            return Ok(());
        }
        let optimized;
        let score: ScoreFunc<'_> = match score {
            Some(s) => s,
            None => {
                optimized = make_optimized_scorer(mol)?;
                &optimized
            }
        };
        let tmp = self.pick_canonical(&res, Some(score))?;
        let check = |ok: bool| {
            if ok {
                Ok(())
            } else {
                Err(RDKitError::invariant(
                    "canonical tautomer does not match the input",
                ))
            }
        };
        check(tmp.num_atoms() == mol.num_atoms())?;
        check(tmp.num_bonds() == mol.num_bonds())?;
        for i in 0..tmp.num_atoms() {
            let ta = tmp.atom(i);
            check(ta.atomic_num() == mol.atom(i).atomic_num())?;
            let a = mol.atom_mut(i);
            a.set_formal_charge(ta.formal_charge());
            a.set_no_implicit(ta.no_implicit());
            a.set_is_aromatic(ta.is_aromatic());
            a.set_num_explicit_hs(ta.num_explicit_hs());
            a.set_num_radical_electrons(ta.num_radical_electrons());
            a.set_chiral_tag(ta.chiral_tag());
        }
        for i in 0..tmp.num_bonds() {
            let tb = tmp.bond(i);
            check(tb.begin_atom_idx() == mol.bond(i).begin_atom_idx())?;
            check(tb.end_atom_idx() == mol.bond(i).end_atom_idx())?;
            let b = mol.bond_mut(i);
            b.set_bond_type(tb.bond_type());
            b.set_bond_dir(tb.bond_dir());
            b.set_is_aromatic(tb.is_aromatic());
            b.set_is_conjugated(tb.is_conjugated());
            if tb.stereo_atoms().len() == 2 {
                let sa = tb.stereo_atoms().to_vec();
                mol.set_bond_stereo_atoms(i, sa[0] as u32, sa[1] as u32)?;
            }
            mol.bond_mut(i).set_stereo(tb.stereo())?;
        }
        mol.update_property_cache(false)
    }
}

/// The smallest key after `key` (`++it` on a `std::map`).
fn next_key(map: &BTreeMap<String, Tautomer>, key: &str) -> Option<String> {
    map.range::<str, _>((Bound::Excluded(key), Bound::Unbounded))
        .next()
        .map(|(k, _)| k.clone())
}

/// The sanitization steps of a new tautomer (its rings are unchanged).
fn sanitize_product(product: &mut RWMol) -> Result<()> {
    product.clear_computed_props(false);
    product.update_property_cache(false)?;
    kekulize(product, true, true, 100)?;
    set_aromaticity(product, AromaticityModel::Default)?;
    set_conjugation(product);
    set_hybridization(product);
    adjust_hs(product)
}

/// `tautomerEnumeratorFromParams`
pub fn tautomer_enumerator_from_params(params: &CleanupParameters) -> Result<TautomerEnumerator> {
    TautomerEnumerator::from_params(params)
}
