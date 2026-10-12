//! Charge handling: the `Reionizer` (strongest acids ionized first) and the
//! `Uncharger` (neutralizes by adding/removing hydrogens).
//!
//! Ported from RDKit `Code/GraphMol/MolStandardize/Charge.{h,cpp}`.

use std::collections::BTreeSet;
use std::sync::{Arc, OnceLock};

use super::CleanupParameters;
use super::catalogs::AcidBaseCatalogParams;
use crate::error::{RDKitError, Result};
use crate::graphmol::atom::is_early_atom;
use crate::graphmol::mol_ops::structure::get_formal_charge;
use crate::graphmol::new_canon::{RankOptions, rank_mol_atoms};
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::smarts::{SmartsParserParams, mol_from_smarts};
use crate::graphmol::substruct::{MatchVect, SubstructMatchParameters, substruct_match};
use crate::rd_log_info;
use crate::rdgeneral::common_properties;
use crate::rdgeneral::dict::RDProps;

/// `ChargeCorrection`: atoms matching `smarts` get the formal charge `charge`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChargeCorrection {
    pub name: String,
    pub smarts: String,
    pub charge: i32,
}

impl ChargeCorrection {
    pub fn new(name: &str, smarts: &str, charge: i32) -> Self {
        ChargeCorrection {
            name: name.to_string(),
            smarts: smarts.to_string(),
            charge,
        }
    }
}

/// `CHARGE_CORRECTIONS`: the default charge corrections.
pub fn charge_corrections() -> Vec<ChargeCorrection> {
    vec![
        ChargeCorrection::new("[Li,Na,K]", "[Li,Na,K;X0+0]", 1),
        ChargeCorrection::new("[Mg,Ca]", "[Mg,Ca;X0+0]", 2),
        ChargeCorrection::new("[Cl]", "[Cl;X0+0]", -1),
    ]
}

fn smarts(text: &str) -> Result<ROMol> {
    mol_from_smarts(text, &SmartsParserParams::default())?
        .ok_or_else(|| RDKitError::Value(format!("Failed parsing SMARTS: {text}")))
}

/// `SubstructMatch(mol, query, std::vector<MatchVectType> &)`
fn all_matches(mol: &ROMol, query: &ROMol) -> Result<Vec<MatchVect>> {
    substruct_match(mol, query, &SubstructMatchParameters::default())
}

/// `SubstructMatch(mol, query, MatchVectType &)`: the first match.
fn first_match(mol: &ROMol, query: &ROMol) -> Result<Option<MatchVect>> {
    let ps = SubstructMatchParameters {
        max_matches: 1,
        ..Default::default()
    };
    Ok(substruct_match(mol, query, &ps)?.into_iter().next())
}

fn mol_name(mol: &ROMol) -> String {
    mol.get_prop(common_properties::_Name)
        .map(|v| v.to_rd_string())
        .unwrap_or_default()
}

// ------------------------------------------------------------------ Reionizer

/// `Reionizer`: ionizes the strongest acids first and moves protons from
/// stronger to weaker acid groups.
#[derive(Debug, Clone)]
pub struct Reionizer {
    params: Arc<AcidBaseCatalogParams>,
    ccs: Vec<ChargeCorrection>,
}

impl Default for Reionizer {
    /// `Reionizer()`: the default acid/base pairs and charge corrections.
    fn default() -> Self {
        Reionizer {
            params: Arc::new(AcidBaseCatalogParams::defaults().clone()),
            ccs: charge_corrections(),
        }
    }
}

/// (position in the acid/base list, atom indices of the match)
type AcidMatch = (usize, Vec<u32>);

impl Reionizer {
    /// `Reionizer(acidbaseFile[, ccs])`; an empty file name selects the
    /// default pairs.
    pub fn from_file(file_name: &str, ccs: Option<Vec<ChargeCorrection>>) -> Result<Self> {
        let params = if file_name.is_empty() {
            AcidBaseCatalogParams::defaults().clone()
        } else {
            AcidBaseCatalogParams::from_file(file_name)?
        };
        Ok(Self::with_params(params, ccs))
    }

    /// `Reionizer(std::istream &, ccs)`
    pub fn from_text(text: &str, ccs: Option<Vec<ChargeCorrection>>) -> Result<Self> {
        Ok(Self::with_params(
            AcidBaseCatalogParams::from_text(text)?,
            ccs,
        ))
    }

    /// `Reionizer(data[, ccs])`
    pub fn from_data(
        data: &[(String, String, String)],
        ccs: Option<Vec<ChargeCorrection>>,
    ) -> Result<Self> {
        Ok(Self::with_params(
            AcidBaseCatalogParams::from_data(data)?,
            ccs,
        ))
    }

    pub fn with_params(params: AcidBaseCatalogParams, ccs: Option<Vec<ChargeCorrection>>) -> Self {
        Reionizer {
            params: Arc::new(params),
            ccs: ccs.unwrap_or_else(charge_corrections),
        }
    }

    pub fn params(&self) -> &AcidBaseCatalogParams {
        &self.params
    }

    pub fn charge_corrections(&self) -> &[ChargeCorrection] {
        &self.ccs
    }

    /// `reionize`
    pub fn reionize(&self, mol: &ROMol) -> Result<ROMol> {
        let mut res = mol.clone();
        self.reionize_in_place(&mut res)?;
        Ok(res)
    }

    /// `reionizeInPlace`
    pub fn reionize_in_place(&self, mol: &mut RWMol) -> Result<()> {
        let abpairs = self.params.pairs();
        if mol.needs_update_property_cache() {
            mol.update_property_cache(false)?;
        }
        let start_charge = get_formal_charge(mol);
        for cc in &self.ccs {
            let ccmol = smarts(&cc.smarts)?;
            for m in all_matches(mol, &ccmol)? {
                for &(_, idx) in &m {
                    let symbol = mol.atom(idx as u32).symbol();
                    rd_log_info!(
                        "Applying charge correction {} {symbol} {}",
                        cc.name,
                        cc.charge
                    );
                    mol.atom_mut(idx as u32).set_formal_charge(cc.charge);
                }
            }
        }
        let current_charge = get_formal_charge(mol);
        let mut charge_diff = current_charge - start_charge;
        // if the molecule is now neutral, assume everything is fixed; if the
        // charge has become more positive, ionize more protonated acid groups
        if current_charge != 0 {
            while charge_diff > 0 {
                let Some((ppos, poccur)) = Self::strongest_protonated(mol, abpairs)? else {
                    break;
                };
                rd_log_info!(
                    "Ionizing {} to balance previous charge corrections",
                    mol_name(&abpairs[ppos].0)
                );
                let pidx = *poccur.last().expect("non-empty match");
                let patom = mol.atom_mut(pidx);
                patom.set_formal_charge(patom.formal_charge() - 1);
                if patom.num_explicit_hs() > 0 {
                    patom.set_num_explicit_hs(patom.num_explicit_hs() - 1);
                }
                mol.update_atom_property_cache(pidx, true)?;
                charge_diff -= 1;
            }
        }

        let mut already_moved: BTreeSet<(u32, u32)> = BTreeSet::new();
        loop {
            let sp = Self::strongest_protonated(mol, abpairs)?;
            let wi = Self::weakest_ionized(mol, abpairs)?;
            let (Some((ppos, poccur)), Some((ipos, ioccur))) = (sp, wi) else {
                break;
            };
            if ppos >= ipos {
                break;
            }
            let pidx = *poccur.last().expect("non-empty match");
            let iidx = *ioccur.last().expect("non-empty match");
            if pidx == iidx {
                // the H would not be moved, resulting in an infinite loop
                rd_log_info!("Aborted reionization due to unexpected situation");
                break;
            }
            if !already_moved.insert((pidx.min(iidx), pidx.max(iidx))) {
                rd_log_info!(
                    "Aborting reionization to avoid infinite loop due \t\t\t\t\t\t\t\tto it being ambiguous where to put a Hydrogen"
                );
                break;
            }
            rd_log_info!(
                "Moved proton from {} to {}",
                mol_name(&abpairs[ppos].0),
                mol_name(&abpairs[ipos].0)
            );
            // remove a hydrogen from the strongest protonated group
            let patom = mol.atom_mut(pidx);
            patom.set_formal_charge(patom.formal_charge() - 1);
            // if there are no implicit Hs to remove automatically, remove an
            // explicit one
            if patom.num_implicit_hs()? == 0 && patom.num_explicit_hs() > 0 {
                patom.set_num_explicit_hs(patom.num_explicit_hs() - 1);
            }
            mol.update_atom_property_cache(pidx, true)?;

            // add a hydrogen to the weakest ionized group
            let iatom = mol.atom_mut(iidx);
            iatom.set_formal_charge(iatom.formal_charge() + 1);
            let valence_list = PeriodicTable::get().get_valence_list(iatom.atomic_num())?;
            let found = valence_list.contains(&(iatom.total_valence()? as i32));
            let patom = mol.atom(pidx);
            let p_aromatic_n_or_p = matches!(patom.atomic_num(), 7 | 15) && patom.is_aromatic();
            let iatom = mol.atom_mut(iidx);
            if iatom.no_implicit() || p_aromatic_n_or_p || !found {
                iatom.set_num_explicit_hs(iatom.num_explicit_hs() + 1);
            }
            mol.update_atom_property_cache(iidx, true)?;
        }
        Ok(())
    }

    /// `strongestProtonated`: the first acid (in list order) that matches.
    fn strongest_protonated(mol: &ROMol, abpairs: &[(ROMol, ROMol)]) -> Result<Option<AcidMatch>> {
        for (position, (acid, _)) in abpairs.iter().enumerate() {
            if let Some(m) = first_match(mol, acid)? {
                return Ok(Some((position, m.iter().map(|&(_, a)| a as u32).collect())));
            }
        }
        Ok(None)
    }

    /// `weakestIonized`: the last base (in list order) that matches.
    fn weakest_ionized(mol: &ROMol, abpairs: &[(ROMol, ROMol)]) -> Result<Option<AcidMatch>> {
        for (position, (_, base)) in abpairs.iter().enumerate().rev() {
            if let Some(m) = first_match(mol, base)? {
                return Ok(Some((position, m.iter().map(|&(_, a)| a as u32).collect())));
            }
        }
        Ok(None)
    }
}

/// `reionizerFromParams`
pub fn reionizer_from_params(params: &CleanupParameters) -> Result<Reionizer> {
    if params.acidbase_data.is_empty() {
        Reionizer::from_file(&params.acidbase_file, None)
    } else {
        Reionizer::from_data(&params.acidbase_data, None)
    }
}

// ------------------------------------------------------------------- Uncharger

struct UnchargerQueries {
    pos_h: ROMol,
    pos_noh: ROMol,
    neg: ROMol,
    neg_acid: ROMol,
}

fn uncharger_queries() -> &'static UnchargerQueries {
    static QUERIES: OnceLock<UnchargerQueries> = OnceLock::new();
    QUERIES.get_or_init(|| UnchargerQueries {
        pos_h: smarts("[+,+2,+3,+4;!h0;!$(*~[-]),$(*(~[-])~[-])]").expect("pos_h"),
        pos_noh: smarts("[+,+2,+3,+4;h0;!$(*~[-])]").expect("pos_noh"),
        neg: smarts("[-!$(*~[+,+2,+3,+4])]").expect("neg"),
        neg_acid: smarts(concat!(
            // carboxylate, carbonate, sulfi(a)te, and their thio-analogues
            // (among other less likely structures)
            "[$([O,S;-][C,S;+0]=[O,S]),",
            // phosphi(a)te, nitrate and their thio-analogues
            "$([O,S;-][N,P;+](=[O,S])[O,S;-]),",
            // hali(a)te, perhalate
            "$([O-][Cl,Br,I;+,+2,+3][O-]),",
            // tetrazole
            "$([n-]1nnnc1),$([n-]1ncnn1)]"
        ))
        .expect("neg_acid"),
    })
}

/// `Uncharger`: neutralizes molecules by adding/removing hydrogens; charges
/// that cannot be removed are balanced unless `force` is set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Uncharger {
    pub canonical_ordering: bool,
    pub force: bool,
    pub protonation_only: bool,
}

impl Default for Uncharger {
    fn default() -> Self {
        Uncharger {
            canonical_ordering: true,
            force: false,
            protonation_only: false,
        }
    }
}

/// `removeCharge`
fn remove_charge(mol: &mut RWMol, idx: u32, charge: i32, h_delta: i32) -> Result<()> {
    let total_hs = mol.total_num_hs(idx, false)? as i32;
    let atom = mol.atom_mut(idx);
    atom.set_num_explicit_hs((total_hs + h_delta) as u32);
    atom.set_no_implicit(true);
    atom.set_formal_charge(atom.formal_charge() - charge);
    rd_log_info!(
        "Removed {} charge.",
        if charge > 0 { "positive" } else { "negative" }
    );
    // the H count changed, so the other valence parameters need updating
    mol.update_atom_property_cache(idx, false)
}

fn h_delta_removing_neg(mol: &ROMol, idx: u32, protonation_only: bool) -> Result<i32> {
    let early = is_early_atom(mol.atom(idx).atomic_num() as i32);
    let has_hs = mol.total_num_hs(idx, false)? > 0;
    if early && (!has_hs || protonation_only) {
        return Ok(0);
    }
    Ok(if early { -1 } else { 1 })
}

fn h_delta_removing_pos(mol: &ROMol, idx: u32, protonation_only: bool) -> Result<i32> {
    let anum = mol.atom(idx).atomic_num();
    // the special case for C here was github #2792
    let carbon_or_early = anum == 6 || is_early_atom(anum as i32);
    if carbon_or_early && protonation_only {
        return Ok(0);
    }
    let has_hs = mol.total_num_hs(idx, false)? > 0;
    if !carbon_or_early && !has_hs {
        return Ok(0);
    }
    Ok(if carbon_or_early { 1 } else { -1 })
}

impl Uncharger {
    /// `Uncharger(canonicalOrdering, force, protonationOnly)`
    pub fn new(canonical_ordering: bool, force: bool, protonation_only: bool) -> Self {
        Uncharger {
            canonical_ordering,
            force,
            protonation_only,
        }
    }

    /// `uncharge`
    pub fn uncharge(&self, mol: &ROMol) -> Result<ROMol> {
        let mut res = mol.clone();
        self.uncharge_in_place(&mut res)?;
        Ok(res)
    }

    /// `unchargeInPlace`
    pub fn uncharge_in_place(&self, mol: &mut RWMol) -> Result<()> {
        rd_log_info!("Running Uncharger");
        if mol.needs_update_property_cache() {
            mol.update_property_cache(false)?;
        }
        let q = uncharger_queries();
        let p_matches = all_matches(mol, &q.pos_h)?;
        let q_matches = all_matches(mol, &q.pos_noh)?;
        let n_matches = all_matches(mol, &q.neg)?;
        let a_matches = all_matches(mol, &q.neg_acid)?;

        // the amount of positive charge that cannot be removed (an unsigned
        // int upstream)
        let mut q_matched: u32 = 0;
        for m in &q_matches {
            q_matched = q_matched.wrapping_add(mol.atom(m[0].1 as u32).formal_charge() as u32);
        }
        for m in &p_matches {
            let idx = m[0].1 as u32;
            if h_delta_removing_pos(mol, idx, self.protonation_only)? == 0 {
                q_matched = q_matched.wrapping_add(mol.atom(idx).formal_charge() as u32);
            }
        }

        let needs_neutralization =
            q_matched > 0 && (!n_matches.is_empty() || !a_matches.is_empty());
        let atom_ranks: Vec<u32> = if self.canonical_ordering && needs_neutralization {
            rank_mol_atoms(mol, RankOptions::default())?
        } else {
            (0..mol.num_atoms()).collect()
        };
        let rank_idx = |mv: &MatchVect| {
            let aidx = mv[0].1;
            (atom_ranks[aidx as usize] as i32, aidx)
        };
        let mut n_atoms: Vec<(i32, i32)> = n_matches.iter().map(rank_idx).collect();
        let mut a_atoms: Vec<(i32, i32)> = a_matches.iter().map(rank_idx).collect();
        if self.canonical_ordering {
            n_atoms.sort();
            a_atoms.sort();
        }

        // the negatively charged sites to neutralize, in order: first the
        // ones that are not acid anions, then the acid anions
        let mut non_acid = vec![true; mol.num_atoms() as usize];
        for &(_, idx) in &a_atoms {
            non_acid[idx as usize] = false;
        }
        let mut neg_atoms: Vec<(i32, i32)> = n_atoms
            .iter()
            .copied()
            .filter(|&(_, idx)| non_acid[idx as usize])
            .collect();
        // the anions of monoprotic acids must not be protonated repeatedly
        let mut skip_charge_sep = vec![0i32; mol.num_atoms() as usize];
        for &(_, idx) in &a_atoms {
            let idx = idx as u32;
            let nbrs: Vec<u32> = mol.atom_neighbors(idx).collect();
            for nbr_idx in nbrs {
                // if the neighbor has a positive charge, neutralize only the
                // negative charges that are not already balanced within the
                // functional group (normally, at most once e.g., NO3-)
                let nbr_charge = mol.atom(nbr_idx).formal_charge();
                if nbr_charge > 0 {
                    if skip_charge_sep[nbr_idx as usize] < nbr_charge {
                        skip_charge_sep[nbr_idx as usize] += 1;
                        skip_charge_sep[idx as usize] = 1;
                    }
                    break;
                }
            }
        }
        neg_atoms.extend(
            a_atoms
                .iter()
                .copied()
                .filter(|&(_, idx)| skip_charge_sep[idx as usize] == 0),
        );

        // surplus negative charges: all the neutralizable ones, minus (unless
        // forced) the positive charge that must stay balanced
        let mut neg_surplus = neg_atoms.len() as i32;
        if !self.force {
            neg_surplus = neg_surplus.wrapping_sub(q_matched as i32);
        }
        if neg_surplus != 0 {
            for &(_, idx) in &neg_atoms {
                let idx = idx as u32;
                let h_delta = h_delta_removing_neg(mol, idx, self.protonation_only)?;
                if h_delta != 0 {
                    remove_charge(mol, idx, -1, h_delta)?;
                    neg_surplus -= 1;
                    if neg_surplus == 0 {
                        break;
                    }
                }
            }
        }

        // the net charge after neutralizing the negative sites
        let mut net_charge: i32 = mol.atoms().iter().map(|a| a.formal_charge()).sum();

        // neutralize the protonated sites until there is no net charge left
        // (or all of them, if forced)
        if net_charge > 0 || self.force {
            let p_idx_matches: Vec<u32> = p_matches
                .iter()
                .flat_map(|m| m.iter().map(|&(_, a)| a as u32))
                .collect();
            for idx in p_idx_matches {
                while mol.atom(idx).formal_charge() > 0 && (net_charge > 0 || self.force) {
                    let h_delta = h_delta_removing_pos(mol, idx, self.protonation_only)?;
                    if h_delta == 0 {
                        break;
                    }
                    remove_charge(mol, idx, 1, h_delta)?;
                    net_charge -= 1;
                }
                if net_charge == 0 && !self.force {
                    break;
                }
            }
        }
        Ok(())
    }
}
