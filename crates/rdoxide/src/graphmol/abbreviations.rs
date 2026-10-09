//! Molecular abbreviations: finding groups such as `CF3` or `OMe` and
//! condensing them into a single labelled pseudo-atom (or marking them with
//! "SUP" substance groups).
//!
//! Ported from RDKit `Code/GraphMol/Abbreviations/{Abbreviations,AbbreviationsUtils}.{h,cpp}`
//! (2026.09.1). The default abbreviations are adapted from the OpenBabel
//! superatom table (GPL-2.0 data originally from the NCI/CADD OSRA project,
//! as carried by RDKit).
//!
//! Copyright (C) 2020 Greg Landrum and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::sync::{Arc, OnceLock};

use super::adjust_query::{AdjustQueryParameters, adjust_flags, adjust_query_properties_in_place};
use super::atom::HybridizationType;
use super::bond::BondType;
use super::find_rings::{SymmetrizeSSSRAlgorithm, fast_find_rings, symmetrize_sssr};
use super::queries::{CompositeQueryType, make_atom_num_query};
use super::query_ops::atom_expand_query;
use super::romol::ROMol;
use super::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use super::substance_group::{SubstanceGroup, add_substance_group};
use super::substruct::{MatchVect, SubstructMatchParameters, substruct_match};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_error, rd_log_warning};

/// Number of dummy atoms of an abbreviation query.
pub const NUM_DUMMIES: &str = "_numDummies";
/// Original atom indices, kept up to date when abbreviations are applied.
pub const ORIG_ATOM_MAPPING: &str = "_origAtomMapping";
/// Original bond indices, kept up to date when abbreviations are applied.
pub const ORIG_BOND_MAPPING: &str = "_origBondMapping";

/// `AbbreviationDefinition`
#[derive(Debug, Clone)]
pub struct AbbreviationDefinition {
    pub label: String,
    pub display_label: String,
    pub display_label_w: String,
    pub smarts: String,
    /// the query molecule (optional)
    pub mol: Option<Arc<ROMol>>,
    pub extra_attach_atoms: Vec<u32>,
    /// whether the definition includes the bonds to the rest of the molecule
    /// (then atom 0 of the query is a dummy)
    pub includes_x_bonds: bool,
}

impl Default for AbbreviationDefinition {
    fn default() -> Self {
        AbbreviationDefinition {
            label: String::new(),
            display_label: String::new(),
            display_label_w: String::new(),
            smarts: String::new(),
            mol: None,
            extra_attach_atoms: Vec::new(),
            includes_x_bonds: true,
        }
    }
}

impl PartialEq for AbbreviationDefinition {
    fn eq(&self, other: &Self) -> bool {
        self.label == other.label
            && self.display_label == other.display_label
            && self.display_label_w == other.display_label_w
            && self.smarts == other.smarts
    }
}

/// `AbbreviationMatch`
#[derive(Debug, Clone, Default, PartialEq)]
pub struct AbbreviationMatch {
    pub match_: MatchVect,
    pub abbrev: AbbreviationDefinition,
}

// ------------------------------------------------------------------ definitions

/// Translations of superatom labels to SMILES: `label SMARTS displayLabel displayLabelW`.
/// The first atom of the SMARTS is the one connected to the rest of the molecule.
const DEFAULT_ABBREVIATIONS: &str = r"CO2Et    C(=O)OCC  CO<sub>2</sub>Et    EtO<sub>2</sub>C
COOEt    C(=O)OCC  CO<sub>2</sub>Et    EtO<sub>2</sub>C
OiBu     OCC(C)C   OiBu     iBuO
nDec     CCCCCCCCCC  nDec
nNon     CCCCCCCCC   nNon
nOct     CCCCCCCC    nOct
nHept    CCCCCCC     nHept
nHex     CCCCCC      nHex
nPent    CCCCC       nPent
iPent    C(C)CCC     iPent
tBu      C(C)(C)C    tBu
iBu      C(C)CC      iBu
nBu      CCCC        nBu
iPr      C(C)C       iPr
nPr      CCC         nPr
Et       CC          Et
NCF3     NC(F)(F)F NCF<sub>3</sub>     F<sub>3</sub>CN
CF3      C(F)(F)F  CF<sub>3</sub>      F<sub>3</sub>C
CCl3     C(Cl)(Cl)Cl CCl<sub>3</sub>     Cl<sub>3</sub>C
CN       C#N       CN       NC
NC       [N+]#[C-] NC       CN 
N(OH)CH3 N([OH])C    N(OH)CH<sub>3</sub> CH<sub>3</sub>(OH)N
NO2      [N+](=O)[O-]  NO<sub>2</sub>      O<sub>2</sub>N
NO       N=O   NO       ON
SO3H     S(=O)(=O)[OH] SO<sub>3</sub>H     HO<sub>3</sub>S
CO2H     C(=O)[OH] CO<sub>2</sub>H     HO<sub>2</sub>C
COOH     C(=O)[OH] COOH     HOOC
OEt      OCC   OEt      EtO
OAc      OC(=O)C   OAc      AcO
NHAc     NC(=O)C   NHAc     AcNH
Ac       C(=O)C    Ac
CHO      C=O   CHO      OHC
NMe      NC    NMe      MeN
SMe      SC    SMe      MeS
OMe      OC    OMe      MeO
CO2-     C(=O)[O-]   COO<sup>-</sup>     <sup>-</sup>OOC
COO-     C(=O)[O-]   COO<sup>-</sup>     <sup>-</sup>OOC";

/// Translations of linker labels to SMILES: the first atom is a dummy
/// connected to the rest of the molecule, the other dummies show the other
/// attachments.
const DEFAULT_LINKERS: &str = r"PEG6  *OCCOCCOCCOCCOCCOCC* PEG6
PEG5  *OCCOCCOCCOCCOCC* PEG5
PEG4  *OCCOCCOCCOCC* PEG4
PEG3  *OCCOCCOCC* PEG3
Dec   *CCCCCCCCCC*
Non   *CCCCCCCCC*
Oct   *CCCCCCCC*
Hept  *CCCCCCC*";

/// `Utils::detail::createAbbreviationMol`
fn create_abbreviation_mol(
    txt: &str,
    remove_extra_dummies: bool,
    allow_connection_to_dummies: bool,
) -> Result<Option<ROMol>> {
    let smarts = if txt.starts_with('*') {
        txt.to_string()
    } else {
        format!("*{txt}")
    };
    let Some(mut q) = mol_from_smarts(&smarts, &SmartsParserParams::default())? else {
        return Ok(None);
    };
    if q.num_atoms() < 2 {
        rd_log_error!("abbreviation with <2 atoms ignored\n");
        return Ok(None);
    }
    let ps = AdjustQueryParameters {
        adjust_degree: true,
        adjust_degree_flags: adjust_flags::ADJUST_IGNOREDUMMIES,
        adjust_ring_count: true,
        adjust_ring_count_flags: adjust_flags::ADJUST_IGNOREDUMMIES,
        ..AdjustQueryParameters::default()
    };
    adjust_query_properties_in_place(&mut q, &ps)?;
    if !allow_connection_to_dummies {
        let mut qry = make_atom_num_query(0);
        qry.negate = true;
        atom_expand_query(q.atom_mut(0), qry, CompositeQueryType::And, true)?;
    }
    let mut n_dummies = smarts.matches('*').count() as u32;
    if remove_extra_dummies {
        q.begin_batch_edit()?;
        for idx in 1..q.num_atoms() {
            let is_null = q
                .atom(idx)
                .query()
                .is_some_and(|query| query.description() == "AtomNull");
            if is_null {
                q.remove_atom(idx)?;
                n_dummies -= 1;
            }
        }
        q.commit_batch_edit()?;
    }
    q.set_prop(NUM_DUMMIES, n_dummies, false)?;
    Ok(Some(q))
}

/// `Utils::parseAbbreviations`: one definition per line
/// (`label SMARTS [displayLabel [displayLabelW]]`). Lines that cannot be
/// turned into a query are skipped.
pub fn parse_abbreviations(
    text: &str,
    remove_extra_dummies: bool,
    allow_connection_to_dummies: bool,
) -> Result<Vec<AbbreviationDefinition>> {
    let mut res = Vec::new();
    for line in text.split('\n').filter(|l| !l.is_empty()) {
        let mut fields = line.split([' ', '\t']).filter(|f| !f.is_empty());
        let Some(label) = fields.next() else {
            continue;
        };
        // a line without a SMARTS yields an empty (and then rejected) query
        let smarts = fields.next().unwrap_or_default();
        let mut defn = AbbreviationDefinition {
            label: label.to_string(),
            smarts: smarts.to_string(),
            display_label: fields.next().unwrap_or_default().to_string(),
            display_label_w: fields.next().unwrap_or_default().to_string(),
            ..AbbreviationDefinition::default()
        };
        defn.mol =
            create_abbreviation_mol(smarts, remove_extra_dummies, allow_connection_to_dummies)?
                .map(Arc::new);
        if defn.mol.is_some() {
            res.push(defn);
        }
    }
    Ok(res)
}

/// `Utils::parseLinkers`
pub fn parse_linkers(text: &str) -> Result<Vec<AbbreviationDefinition>> {
    parse_abbreviations(text, true, true)
}

/// `Utils::getDefaultAbbreviations`
pub fn default_abbreviations() -> Vec<AbbreviationDefinition> {
    static DEFS: OnceLock<Vec<AbbreviationDefinition>> = OnceLock::new();
    DEFS.get_or_init(|| {
        parse_abbreviations(DEFAULT_ABBREVIATIONS, false, false)
            .expect("the built-in abbreviations are valid")
    })
    .clone()
}

/// `Utils::getDefaultLinkers`
pub fn default_linkers() -> Vec<AbbreviationDefinition> {
    static DEFS: OnceLock<Vec<AbbreviationDefinition>> = OnceLock::new();
    DEFS.get_or_init(|| {
        parse_abbreviations(DEFAULT_LINKERS, true, true).expect("the built-in linkers are valid")
    })
    .clone()
}

// ------------------------------------------------------------------ applying

fn uint_vec_prop(mol: &ROMol, key: &str) -> Option<Vec<u32>> {
    match mol.get_prop(key)? {
        RDValue::VecUnsignedInt(v) => Some(v.clone()),
        _ => None,
    }
}

fn atom_bond_indices(mol: &ROMol, atom: u32) -> Vec<u32> {
    mol.atom_bonds(atom).map(|b| b.idx()).collect()
}

/// `applyMatches`: replaces every matched group by one labelled atom.
pub fn apply_matches(mol: &mut ROMol, matches: &[AbbreviationMatch]) -> Result<()> {
    let mut atoms_to_remove = vec![false; mol.num_atoms() as usize];
    let mut bonds_to_remove = vec![false; mol.num_bonds() as usize];
    let prev = uint_vec_prop(mol, ORIG_ATOM_MAPPING).zip(uint_vec_prop(mol, ORIG_BOND_MAPPING));
    let orig_bond = |idx: u32| prev.as_ref().map_or(idx, |(_, b)| b[idx as usize]);
    let orig_atom = |idx: u32| prev.as_ref().map_or(idx, |(a, _)| a[idx as usize]);
    let mut added_bonds: Vec<u32> = Vec::new();

    for amatch in matches {
        // with bonds, atom 0 of the match is the dummy and atom 1 becomes the
        // abbreviation; without, atom 0 is the abbreviation
        let which_atom = usize::from(amatch.abbrev.includes_x_bonds);
        let connect_idx = amatch.match_[which_atom].1 as u32;
        let connecting = mol.atom_mut(connect_idx);
        connecting.set_prop(cp::atomLabel, amatch.abbrev.label.as_str(), false)?;
        if !amatch.abbrev.display_label.is_empty() {
            connecting.set_prop(
                cp::_displayLabel,
                amatch.abbrev.display_label.as_str(),
                false,
            )?;
        }
        if !amatch.abbrev.display_label_w.is_empty() {
            connecting.set_prop(
                cp::_displayLabelW,
                amatch.abbrev.display_label_w.as_str(),
                false,
            )?;
        }
        connecting.set_formal_charge(0);
        connecting.set_atomic_num(0);
        connecting.set_isotope(0);
        connecting.set_is_aromatic(false);
        // so that these are drawn linearly
        connecting.set_hybridization(HybridizationType::SP);

        for &(query_idx, mol_idx) in &amatch.match_[which_atom + 1..] {
            let mol_idx = mol_idx as u32;
            if atoms_to_remove[mol_idx as usize] {
                return Err(RDKitError::invariant("overlapping matches"));
            }
            atoms_to_remove[mol_idx as usize] = true;
            let bonds = atom_bond_indices(mol, mol_idx);
            for &b in &bonds {
                bonds_to_remove[b as usize] = true;
            }
            // with a query molecule, check whether bonds to atoms outside the
            // match must be re-attached to the abbreviation
            let Some(query) = &amatch.abbrev.mol else {
                continue;
            };
            let query_degree = query.atom_bonds(query_idx as u32).count();
            if bonds.len() <= query_degree {
                continue;
            }
            for &b in &bonds {
                let nbr = mol.bond(b).other_atom_idx(mol_idx)?;
                if !amatch.match_.iter().any(|&(_, m)| m == nbr as i32) {
                    mol.add_bond(nbr, connect_idx, BondType::Single)?;
                    added_bonds.push(orig_bond(b));
                }
            }
        }
        // connect the extra attachment atoms to the abbreviation
        for &oaidx in &amatch.abbrev.extra_attach_atoms {
            let mut bond_idx = None;
            for b in atom_bond_indices(mol, oaidx) {
                if bonds_to_remove.get(b as usize).copied().unwrap_or(false) {
                    if bond_idx.is_some() {
                        return Err(RDKitError::invariant("bondIdx must be unique"));
                    }
                    bond_idx = Some(b);
                }
            }
            let bond_idx = bond_idx.ok_or_else(|| RDKitError::invariant("bondIdx not found"))?;
            mol.add_bond(oaidx, connect_idx, BondType::Single)?;
            added_bonds.push(orig_bond(bond_idx));
        }
    }

    let mut atom_mapping = Vec::with_capacity(mol.num_atoms() as usize);
    let mut bond_mapping = Vec::with_capacity(mol.num_bonds() as usize);
    mol.begin_batch_edit()?;
    for (i, &remove) in atoms_to_remove.iter().enumerate() {
        if remove {
            mol.remove_atom(i as u32)?;
        } else {
            atom_mapping.push(orig_atom(i as u32));
        }
    }
    for (i, &remove) in bonds_to_remove.iter().enumerate() {
        if !remove {
            bond_mapping.push(orig_bond(i as u32));
        }
    }
    mol.commit_batch_edit()?;
    bond_mapping.extend(added_bonds);
    if atom_mapping.len() != mol.num_atoms() as usize {
        return Err(RDKitError::invariant("atomMapping mismatch"));
    }
    if bond_mapping.len() != mol.num_bonds() as usize {
        return Err(RDKitError::invariant("bondMapping mismatch"));
    }
    mol.set_prop(ORIG_ATOM_MAPPING, atom_mapping, false)?;
    mol.set_prop(ORIG_BOND_MAPPING, bond_mapping, false)?;
    Ok(())
}

/// `labelMatches`: adds a "SUP" substance group per match.
pub fn label_matches(mol: &mut ROMol, matches: &[AbbreviationMatch]) -> Result<()> {
    for amatch in matches {
        // atom 0 of the match is the dummy
        let mut sg = SubstanceGroup::new("SUP");
        sg.set_prop("LABEL", amatch.abbrev.label.as_str(), false)?;
        for &(_, mol_idx) in &amatch.match_[1..] {
            sg.add_atom_with_idx(mol, mol_idx as u32)?;
        }
        let bond = mol
            .get_bond_between_atoms(amatch.match_[0].1 as u32, amatch.match_[1].1 as u32)?
            .map(|b| b.idx())
            .ok_or_else(|| RDKitError::invariant("bond to attachment point not found"))?;
        sg.add_bond_with_idx(mol, bond)?;
        sg.add_attach_point(amatch.match_[1].1 as u32, amatch.match_[0].1, "1");
        add_substance_group(mol, sg);
    }
    Ok(())
}

/// `findApplicableAbbreviationMatches`: the non-overlapping matches of
/// `abbrevs` (used in order). Abbreviations covering `max_coverage` or more
/// of the atoms (dummies not counted) are ignored.
pub fn find_applicable_abbreviation_matches(
    mol: &ROMol,
    abbrevs: &[AbbreviationDefinition],
    max_coverage: f64,
) -> Result<Vec<AbbreviationMatch>> {
    let n_atoms = mol.num_atoms();
    if n_atoms == 0 || abbrevs.is_empty() {
        return Ok(Vec::new());
    }
    // ring information is needed by the queries
    let mut with_rings;
    let mol = if mol.ring_info().is_find_fast_or_better() {
        mol
    } else {
        with_rings = mol.clone();
        fast_find_rings(&mut with_rings)?;
        &with_rings
    };

    let n = n_atoms as usize;
    let mut first_ats = vec![false; n];
    let mut covered = vec![false; n];
    let mut tentative = Vec::new();
    for abbrev in abbrevs {
        let query = abbrev
            .mol
            .as_ref()
            .ok_or_else(|| RDKitError::invariant("molecule is null"))?;
        if max_coverage > 0.0 {
            let n_dummies = query
                .get_prop(NUM_DUMMIES)
                .and_then(RDValue::as_uint)
                .ok_or_else(|| RDKitError::Key(NUM_DUMMIES.into()))?;
            if f64::from(query.num_atoms() - n_dummies) / f64::from(n_atoms) >= max_coverage {
                continue;
            }
        }
        for m in substruct_match(mol, query, &SubstructMatchParameters::default())? {
            if m.len() < 2 {
                return Err(RDKitError::invariant("bad match size"));
            }
            let first = m[1].1 as usize;
            // skip matches whose first non-dummy atom is already covered or
            // was used as a first atom, and those overlapping covered atoms
            if first_ats[first] || covered[first] || m[2..].iter().any(|p| covered[p.1 as usize]) {
                continue;
            }
            for p in &m[1..] {
                covered[p.1 as usize] = true;
            }
            first_ats[first] = true;
            if !first_ats[m[0].1 as usize] {
                tentative.push(AbbreviationMatch {
                    match_: m,
                    abbrev: abbrev.clone(),
                });
            }
        }
    }
    // the dummy of a match must not be the first atom of another match
    Ok(tentative
        .into_iter()
        .filter(|itm| !first_ats[itm.match_[0].1 as usize])
        .collect())
}

/// `condenseMolAbbreviations`
pub fn condense_mol_abbreviations(
    mol: &mut ROMol,
    abbrevs: &[AbbreviationDefinition],
    max_coverage: f64,
    sanitize: bool,
) -> Result<()> {
    let applicable = find_applicable_abbreviation_matches(mol, abbrevs, max_coverage)?;
    apply_matches(mol, &applicable)?;
    if sanitize && !mol.ring_info().is_symm_sssr() {
        symmetrize_sssr(mol, SymmetrizeSSSRAlgorithm::Default, true, false, false)?;
    }
    Ok(())
}

/// `labelMolAbbreviations`
pub fn label_mol_abbreviations(
    mol: &mut ROMol,
    abbrevs: &[AbbreviationDefinition],
    max_coverage: f64,
) -> Result<()> {
    let applicable = find_applicable_abbreviation_matches(mol, abbrevs, max_coverage)?;
    label_matches(mol, &applicable)
}

/// `condenseAbbreviationSubstanceGroups`: replaces the "SUP" groups by labelled atoms.
pub fn condense_abbreviation_substance_groups(mol: &mut ROMol) -> Result<()> {
    let mut matches = Vec::new();
    for sg in mol.substance_groups() {
        if sg.get_prop("TYPE").map(RDValue::to_rd_string).as_deref() != Some("SUP") {
            continue;
        }
        let label = sg
            .get_prop("LABEL")
            .map_or_else(|| "abbrev".to_string(), RDValue::to_rd_string);
        let mut abbrev = AbbreviationDefinition {
            label,
            ..AbbreviationDefinition::default()
        };
        let mut atoms = sg.atoms().to_vec();
        let bonds = sg.bonds();
        if bonds.is_empty() {
            rd_log_warning!("SUP group without any bonds\n");
            abbrev.includes_x_bonds = false;
        } else {
            abbrev.includes_x_bonds = true;
            let mut first_attach_found = false;
            for &b in bonds {
                let bond = mol.bond(b);
                let (begin, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
                // m_at: the group atom of the bond; o_at: the other one
                let (m_at, o_at) = if atoms.contains(&begin) {
                    (begin, end)
                } else if atoms.contains(&end) {
                    (end, begin)
                } else {
                    rd_log_warning!(
                        "SUP group includes bond not connected to any of the abbreviation atoms\n"
                    );
                    continue;
                };
                if first_attach_found {
                    abbrev.extra_attach_atoms.push(o_at);
                } else {
                    // the atom at the first attachment point comes first in the match
                    if atoms.first() != Some(&m_at) {
                        atoms.retain(|&a| a != m_at);
                        atoms.insert(0, m_at);
                    }
                    atoms.insert(0, o_at);
                    first_attach_found = true;
                }
            }
        }
        matches.push(AbbreviationMatch {
            match_: atoms
                .iter()
                .enumerate()
                .map(|(i, &a)| (i as i32, a as i32))
                .collect(),
            abbrev,
        });
    }
    if matches.is_empty() {
        rd_log_warning!("no suitable SubstanceGroups found\n");
        Ok(())
    } else {
        apply_matches(mol, &matches)
    }
}
