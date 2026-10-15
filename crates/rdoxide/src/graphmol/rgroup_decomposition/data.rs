//! The R groups found for a molecule and the matches they belong to
//! (upstream `RGroupDecomposition/RGroupData.{h,cpp}` and `RGroupMatch.h`).

use std::collections::{BTreeMap, BTreeSet};
use std::sync::{Arc, Mutex, OnceLock, RwLock};

use super::utils::{RLABEL, SIDECHAIN_RLABELS};
use crate::datastructs::{BitVect, ExplicitBitVect};
use crate::geometry::Point3D;
use crate::graphmol::chem_transforms::combine_mols;
use crate::graphmol::fingerprints::legacy::{MorganKind, MorganResult, morgan_fingerprint};
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, remove_hs};
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_debug};

/// Number of bits of the fingerprints used for scoring.
pub const FINGERPRINT_SIZE: usize = 512;

/// Morgan fingerprint of an R group, used by the fingerprint variance score.
#[derive(Debug, Clone)]
pub struct RGroupFingerprint {
    pub bits: ExplicitBitVect,
    pub on_bits: Vec<i32>,
}

/// `addFingerprintToRGroupData`: the fingerprint of an R group molecule, with
/// the attachment points turned into boron atoms.
pub fn rgroup_fingerprint(mol: &RWMol, smiles: &str) -> Result<RGroupFingerprint> {
    let mut mol = mol.clone();
    for idx in 0..mol.num_atoms() {
        let atom = mol.atom_mut(idx);
        if atom.atomic_num() == 0 {
            atom.set_atomic_num(5);
            if atom.isotope() > 0 {
                atom.set_isotope(0);
            }
        }
    }
    if sanitize_mol(&mut mol, SANITIZE_ALL).is_err() {
        rd_log_debug!("Failed to sanitize RGroup fingerprint mol for {smiles}");
    }
    let MorganResult::Bits(bits) = morgan_fingerprint(
        &mol,
        2,
        MorganKind::BitVect(FINGERPRINT_SIZE as u32),
        None,
        None,
        false,
        true,
        None,
        false,
    )?
    else {
        return Err(RDKitError::invariant("expected a bit vector"));
    };
    let on_bits = bits.on_bits();
    Ok(RGroupFingerprint { bits, on_bits })
}

/// `std::regex_replace(smiles, "\\[\\d*\\*\\]", "*")`: drops the isotope
/// labels of dummy atoms.
fn remove_dummy_isotopes(smiles: &str) -> String {
    let bytes = smiles.as_bytes();
    let mut out = String::with_capacity(smiles.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'[' {
            let mut j = i + 1;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            if j + 1 < bytes.len() && bytes[j] == b'*' && bytes[j + 1] == b']' {
                out.push('*');
                i = j + 2;
                continue;
            }
        }
        // copy one UTF-8 character
        let ch_len = smiles[i..].chars().next().map_or(1, char::len_utf8);
        out.push_str(&smiles[i..i + ch_len]);
        i += ch_len;
    }
    out
}

/// `RGroupData`: all the side chains attached to one R-label.
#[derive(Debug)]
pub struct RGroupData {
    /// all the side chains as one molecule
    pub combined_mol: RWMol,
    /// all the molecules in the R group
    pub mols: Vec<Arc<ROMol>>,
    /// used for R group equivalence
    pub smiles_vect: Vec<String>,
    /// smiles of all the molecules in the R group (with attachments)
    pub smiles: String,
    /// core attachment points
    pub attachments: BTreeSet<i32>,
    fingerprint: OnceLock<RGroupFingerprint>,
    pub is_hydrogen: bool,
    pub single_fragment: bool,
    pub is_linker: bool,
    pub labelled: bool,
}

/// An R group shared between the matches of a decomposition.
pub type RData = Arc<RwLock<RGroupData>>;

impl Default for RGroupData {
    fn default() -> Self {
        RGroupData {
            combined_mol: RWMol::new(),
            mols: Vec::new(),
            smiles_vect: Vec::new(),
            smiles: String::new(),
            attachments: BTreeSet::new(),
            fingerprint: OnceLock::new(),
            is_hydrogen: false,
            single_fragment: true,
            is_linker: false,
            labelled: false,
        }
    }
}

impl RGroupData {
    pub fn new() -> Self {
        Self::default()
    }

    /// `getRGroupLabel`
    pub fn rgroup_label(rlabel: i32) -> String {
        format!("R{rlabel}")
    }

    pub const CORE_LABEL: &'static str = "Core";
    pub const MOL_LABEL: &'static str = "Mol";

    /// `isMolHydrogen`: only hydrogens and labelled attachment points.
    pub fn is_mol_hydrogen(mol: &ROMol) -> bool {
        mol.atoms()
            .iter()
            .all(|a| a.atomic_num() == 1 || (a.atomic_num() == 0 && a.has_prop(SIDECHAIN_RLABELS)))
    }

    fn smiles_of_parts(&self) -> String {
        self.smiles_vect.join(".")
    }

    fn merge_into_combined_mol(&mut self, mol: &Arc<ROMol>, has_combined: bool) -> Result<()> {
        if !has_combined {
            self.combined_mol = (**mol).clone();
        } else {
            self.combined_mol = combine_mols(&self.combined_mol, mol, Point3D::default())?;
            self.single_fragment = false;
        }
        self.smiles = self.smiles_of_parts();
        let smiles = self.smiles.clone();
        self.combined_mol
            .dict_mut()
            .set_val(cp::internalRgroupSmiles, smiles);
        // upstream also merges the target atom and bond indices of the
        // molecules here, but never stores the result: the combined molecule
        // keeps those of the first side chain
        Ok(())
    }

    /// `RGroupData::add`: adds a side chain attached at `rlabel_attachments`.
    pub fn add(&mut self, new_mol: &Arc<ROMol>, rlabel_attachments: &[i32]) -> Result<()> {
        if self.mols.iter().any(|m| Arc::ptr_eq(m, new_mol)) {
            return Ok(());
        }
        if !self.mols.is_empty() {
            if Self::is_mol_hydrogen(new_mol) {
                return Ok(());
            }
            if self.is_hydrogen {
                self.combined_mol = RWMol::new();
                self.smiles_vect.clear();
                self.attachments.clear();
                self.mols.clear();
                self.fingerprint = OnceLock::new();
            }
        }
        self.labelled = false;
        self.attachments.extend(rlabel_attachments.iter().copied());
        let has_combined = !self.mols.is_empty();
        self.mols.push(new_mol.clone());
        let mut copy = (**new_mol).clone();
        let smiles = mol_to_smiles(&mut copy, &SmilesWriteParams::default())?;
        self.smiles_vect.push(remove_dummy_isotopes(&smiles));
        self.merge_into_combined_mol(new_mol, has_combined)?;
        self.is_hydrogen = self.mols.iter().all(|m| Self::is_mol_hydrogen(m));
        self.is_linker = self.single_fragment && self.attachments.len() > 1;
        Ok(())
    }

    /// `getNumBondsToRlabels`
    pub fn num_bonds_to_rlabels(&self) -> BTreeMap<i32, i32> {
        let mut used = BTreeMap::new();
        for atom in self.combined_mol.atoms() {
            if let Some(rlabel) = atom.get_prop(RLABEL).and_then(RDValue::as_int) {
                *used.entry(rlabel).or_insert(0) += 1;
            }
        }
        used
    }

    /// The fingerprint used for scoring, computed on first use.
    pub fn fingerprint(&self) -> Result<&RGroupFingerprint> {
        if let Some(fp) = self.fingerprint.get() {
            return Ok(fp);
        }
        let fp = rgroup_fingerprint(&self.combined_mol, &self.smiles)?;
        Ok(self.fingerprint.get_or_init(|| fp))
    }

    /// `toString`
    pub fn to_string_repr(&self) -> String {
        let attachments = self
            .attachments
            .iter()
            .map(i32::to_string)
            .collect::<Vec<_>>()
            .join(",");
        format!("RG {attachments} {}", self.smiles_of_parts())
    }
}

/// The target molecule of a match, kept for highlighting the R groups.
#[derive(Debug)]
struct HighlightTarget {
    mol: RWMol,
    hs_trimmed: bool,
}

/// `RGroupMatch`: one way of decomposing a molecule on a core.
#[derive(Debug, Clone)]
pub struct RGroupMatch {
    /// index of the matching core
    pub core_idx: i32,
    pub number_missing_user_rgroups: usize,
    /// R-label to R group
    pub rgroups: BTreeMap<i32, RData>,
    /// the core with dummy or query atoms and bonds matched
    pub matched_core: Arc<Mutex<RWMol>>,
    target_mol_for_highlights: Option<Arc<Mutex<HighlightTarget>>>,
}

impl RGroupMatch {
    pub fn new(
        core_idx: i32,
        number_missing_user_rgroups: usize,
        rgroups: BTreeMap<i32, RData>,
        matched_core: RWMol,
    ) -> Self {
        RGroupMatch {
            core_idx,
            number_missing_user_rgroups,
            rgroups,
            matched_core: Arc::new(Mutex::new(matched_core)),
            target_mol_for_highlights: None,
        }
    }

    /// `toString`
    pub fn to_string_repr(&self) -> String {
        let mut groups = String::new();
        for (label, rgroup) in &self.rgroups {
            let rg = rgroup.read().expect("rgroup lock");
            groups.push_str(&format!("\n\t({label}:{})", rg.to_string_repr()));
        }
        format!(
            "Match coreIdx {} missing count {} {groups}",
            self.core_idx, self.number_missing_user_rgroups
        )
    }

    /// `setTargetMoleculeForHighlights`
    pub fn set_target_molecule_for_highlights(&mut self, target: RWMol) {
        self.target_mol_for_highlights = Some(Arc::new(Mutex::new(HighlightTarget {
            mol: target,
            hs_trimmed: false,
        })));
    }

    /// `getTargetMoleculeForHighlights`: with `trim_hs` the hydrogens of the
    /// target are removed, keeping the atom and bond indices of the R groups
    /// consistent.
    pub fn target_molecule_for_highlights(&self, trim_hs: bool) -> Result<Option<RWMol>> {
        let Some(shared) = &self.target_mol_for_highlights else {
            return Ok(None);
        };
        let mut target = shared.lock().expect("highlight target lock");
        if trim_hs && !target.hs_trimmed {
            self.trim_target_hs(&mut target.mol)?;
            target.hs_trimmed = true;
        }
        Ok(Some(target.mol.clone()))
    }

    fn trim_target_hs(&self, target: &mut RWMol) -> Result<()> {
        let num_atoms = target.num_atoms() as usize;
        let mut stored_map_nums = vec![0; num_atoms];
        for (idx, stored) in stored_map_nums.iter_mut().enumerate() {
            let atom = target.atom_mut(idx as u32);
            *stored = atom.atom_map_num();
            atom.set_atom_map_num(0, false)?;
        }
        let mut old_bond_ends = vec![(-1i32, -1i32); target.num_bonds() as usize];
        for rgroup in self.rgroups.values() {
            let rg = rgroup.read().expect("rgroup lock");
            let Some(RDValue::VecInt(bond_indices)) =
                rg.combined_mol.get_prop(cp::_rgroupTargetBonds)
            else {
                continue;
            };
            for &bond_idx in bond_indices {
                let bond = target.bond(bond_idx as u32);
                let (begin, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
                old_bond_ends[bond_idx as usize] = (begin as i32, end as i32);
                target
                    .atom_mut(begin)
                    .set_atom_map_num(begin as i32 + 1, false)?;
                target
                    .atom_mut(end)
                    .set_atom_map_num(end as i32 + 1, false)?;
            }
        }
        let mut old_to_new = vec![-1i32; num_atoms];
        let rhps = RemoveHsParameters {
            remove_mapped: false,
            ..Default::default()
        };
        remove_hs(target, &rhps, true)?;
        for idx in 0..target.num_atoms() {
            let atom = target.atom_mut(idx);
            let map_num = atom.atom_map_num();
            if map_num != 0 {
                let old = (map_num - 1) as usize;
                old_to_new[old] = idx as i32;
                atom.set_atom_map_num(stored_map_nums[old], false)?;
            }
        }
        for rgroup in self.rgroups.values() {
            let mut rg = rgroup.write().expect("rgroup lock");
            let combined = &mut rg.combined_mol;
            let mut atom_indices: Vec<i32> = match combined.get_prop(cp::_rgroupTargetAtoms) {
                Some(RDValue::VecInt(v)) => v.clone(),
                _ => Vec::new(),
            };
            for idx in &mut atom_indices {
                let new_idx = old_to_new[*idx as usize];
                if new_idx == -1 {
                    return Err(RDKitError::invariant("newAtomIdx must be >=0"));
                }
                *idx = new_idx;
            }
            combined
                .dict_mut()
                .set_val(cp::_rgroupTargetAtoms, atom_indices);
            let bond_indices = match combined.get_prop(cp::_rgroupTargetBonds) {
                Some(RDValue::VecInt(v)) => Some(v.clone()),
                _ => None,
            };
            let mut bond_indices = bond_indices.unwrap_or_default();
            for idx in &mut bond_indices {
                let (old_begin, old_end) = old_bond_ends[*idx as usize];
                if old_begin == -1 || old_end == -1 {
                    return Err(RDKitError::invariant("oldPair members must be >=0"));
                }
                let (new_begin, new_end) =
                    (old_to_new[old_begin as usize], old_to_new[old_end as usize]);
                if new_begin == -1 || new_end == -1 {
                    return Err(RDKitError::invariant(
                        "newBeginAtomIdx and newEndAtomIdx must be >=0",
                    ));
                }
                let bond = target
                    .bond_idx_between(new_begin as u32, new_end as u32)
                    .ok_or_else(|| RDKitError::invariant("bond must not be null"))?;
                *idx = bond as i32;
            }
            combined
                .dict_mut()
                .set_val(cp::_rgroupTargetBonds, bond_indices);
        }
        Ok(())
    }
}
