//! Fragment removal and the largest-fragment chooser.
//!
//! Ported from RDKit `Code/GraphMol/MolStandardize/Fragment.{h,cpp}`.

use std::sync::Arc;

use super::CleanupParameters;
use super::catalogs::FragmentCatalogParams;
use crate::error::Result;
use crate::graphmol::mol_ops::frags::{get_mol_frags_as_mols, get_mol_frags_atoms};
use crate::graphmol::mol_ops::sanitize_flags::SANITIZE_ALL;
use crate::graphmol::mol_ops::sanitize_mol;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_fragment_to_smiles};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rd_log_info;
use crate::rdgeneral::common_properties;
use crate::rdgeneral::dict::RDProps;

/// `FragmentRemover`: removes the fragments that exactly match one of the
/// catalog's patterns.
#[derive(Debug, Clone)]
pub struct FragmentRemover {
    params: Arc<FragmentCatalogParams>,
    /// leave at least one fragment, even if all of them match
    pub leave_last: bool,
    /// return the molecule unchanged if every fragment matches
    pub skip_if_all_match: bool,
}

impl Default for FragmentRemover {
    /// `FragmentRemover()`: the default fragments, `leave_last`.
    fn default() -> Self {
        rd_log_info!("Initializing FragmentRemover");
        FragmentRemover {
            params: Arc::new(FragmentCatalogParams::defaults().clone()),
            leave_last: true,
            skip_if_all_match: false,
        }
    }
}

impl FragmentRemover {
    /// `FragmentRemover(fragmentFile, leave_last, skip_if_all_match)`; an
    /// empty file name selects the default fragments.
    pub fn from_file(file_name: &str, leave_last: bool, skip_if_all_match: bool) -> Result<Self> {
        Ok(Self::with_params(
            FragmentCatalogParams::from_file(file_name)?,
            leave_last,
            skip_if_all_match,
        ))
    }

    /// `FragmentRemover(std::istream &, leave_last, skip_if_all_match)`
    pub fn from_text(text: &str, leave_last: bool, skip_if_all_match: bool) -> Result<Self> {
        Ok(Self::with_params(
            FragmentCatalogParams::from_text(text)?,
            leave_last,
            skip_if_all_match,
        ))
    }

    /// `FragmentRemover(data, leave_last, skip_if_all_match)`
    pub fn from_data(
        data: &[(String, String)],
        leave_last: bool,
        skip_if_all_match: bool,
    ) -> Result<Self> {
        Ok(Self::with_params(
            FragmentCatalogParams::from_data(data)?,
            leave_last,
            skip_if_all_match,
        ))
    }

    pub fn with_params(
        params: FragmentCatalogParams,
        leave_last: bool,
        skip_if_all_match: bool,
    ) -> Self {
        FragmentRemover {
            params: Arc::new(params),
            leave_last,
            skip_if_all_match,
        }
    }

    pub fn params(&self) -> &FragmentCatalogParams {
        &self.params
    }

    /// `remove`
    pub fn remove(&self, mol: &ROMol) -> Result<ROMol> {
        let mut res = mol.clone();
        self.remove_in_place(&mut res)?;
        Ok(res)
    }

    /// `removeInPlace`
    pub fn remove_in_place(&self, mol: &mut RWMol) -> Result<()> {
        rd_log_info!("Running FragmentRemover");
        let (frag_mols, atom_frag_mapping) = get_mol_frags_as_mols(mol, false, false)?;
        // (fragment, original fragment index)
        let mut frags: Vec<(ROMol, usize)> = frag_mols.into_iter().zip(0..).collect();
        let ps = SubstructMatchParameters {
            max_matches: 1,
            ..Default::default()
        };
        for fgci in self.params.func_groups() {
            let o_count = frags.len();
            if o_count == 0 {
                break;
            }
            let mut kept = Vec::with_capacity(o_count);
            for frag in &frags {
                let matches = fgci.num_atoms() == frag.0.num_atoms()
                    && fgci.num_bonds() == frag.0.num_bonds()
                    && !substruct_match(&frag.0, fgci, &ps)?.is_empty();
                kept.push(!matches);
            }
            if self.leave_last && !self.skip_if_all_match && !kept.contains(&true) {
                // all the remaining fragments match this pattern: keep them all
                break;
            }
            let mut keep = kept.into_iter();
            frags.retain(|_| keep.next().unwrap_or(true));
            if frags.len() != o_count {
                let name = fgci
                    .get_prop(common_properties::_Name)
                    .map(|v| v.to_rd_string())
                    .unwrap_or_default();
                rd_log_info!("Removed fragment: {name}");
            }
        }
        if frags.is_empty() {
            if self.skip_if_all_match {
                rd_log_info!("All fragments matched; original molecule returned.");
            } else {
                mol.begin_batch_edit()?;
                for i in 0..mol.num_atoms() {
                    mol.remove_atom(i)?;
                }
                mol.commit_batch_edit()?;
            }
            return Ok(());
        }
        let mut remove = vec![true; mol.num_atoms() as usize];
        for (_, frag_idx) in &frags {
            for &atom_idx in &atom_frag_mapping[*frag_idx] {
                remove[atom_idx as usize] = false;
            }
        }
        mol.begin_batch_edit()?;
        for (i, _) in remove.iter().enumerate().filter(|(_, r)| **r) {
            mol.remove_atom(i as u32)?;
        }
        mol.commit_batch_edit()
    }
}

/// `fragmentRemoverFromParams(params, leave_last, skip_if_all_match)`
pub fn fragment_remover_from_params(
    params: &CleanupParameters,
    leave_last: bool,
    skip_if_all_match: bool,
) -> Result<FragmentRemover> {
    if params.fragment_data.is_empty() {
        if params.fragment_file.is_empty() {
            return Ok(FragmentRemover::with_params(
                FragmentCatalogParams::defaults().clone(),
                leave_last,
                skip_if_all_match,
            ));
        }
        FragmentRemover::from_file(&params.fragment_file, leave_last, skip_if_all_match)
    } else {
        FragmentRemover::from_data(&params.fragment_data, leave_last, skip_if_all_match)
    }
}

/// `LargestFragmentChooser`: keeps the largest fragment (by atom count
/// including Hs, then a rough weight, then SMILES).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LargestFragmentChooser {
    pub prefer_organic: bool,
    pub use_atom_count: bool,
    pub count_heavy_atoms_only: bool,
}

impl Default for LargestFragmentChooser {
    fn default() -> Self {
        LargestFragmentChooser::new(false)
    }
}

impl LargestFragmentChooser {
    /// `LargestFragmentChooser(preferOrganic)`
    pub fn new(prefer_organic: bool) -> Self {
        LargestFragmentChooser {
            prefer_organic,
            use_atom_count: true,
            count_heavy_atoms_only: false,
        }
    }

    /// `LargestFragmentChooser(const CleanupParameters &)`
    pub fn from_params(params: &CleanupParameters) -> Self {
        LargestFragmentChooser {
            prefer_organic: params.prefer_organic,
            use_atom_count: params.largest_fragment_chooser_use_atom_count,
            count_heavy_atoms_only: params.largest_fragment_chooser_count_heavy_atoms_only,
        }
    }

    /// `choose`: the largest fragment, sanitized.
    pub fn choose(&self, mol: &ROMol) -> Result<ROMol> {
        let mut res = mol.clone();
        self.choose_in_place(&mut res)?;
        sanitize_mol(&mut res, SANITIZE_ALL).map_err(|(_, e)| e)?;
        Ok(res)
    }

    /// `chooseInPlace`
    pub fn choose_in_place(&self, mol: &mut RWMol) -> Result<()> {
        rd_log_info!("Running LargestFragmentChooser");
        if mol.num_atoms() == 0 {
            return Ok(());
        }
        let frags = get_mol_frags_atoms(mol);
        if frags.len() == 1 {
            return Ok(());
        }
        let ps = SmilesWriteParams::default();
        let mut best: Option<usize> = None;
        let (mut l_smiles, mut l_num_atoms, mut l_weight, mut l_organic) =
            (String::new(), 0u32, 0.0f64, false);
        for (fidx, frag) in frags.iter().enumerate() {
            let smiles = mol_fragment_to_smiles(mol, &ps, frag, None, None, None)?;
            rd_log_info!("Fragment: {smiles}");
            let organic = frag.iter().any(|&i| mol.atom(i as u32).atomic_num() == 6);
            if self.prefer_organic {
                // skip non-organic fragments once an organic one was found
                if best.is_some() && l_organic && !organic {
                    continue;
                }
                // an organic fragment beats any non-organic one
                if best.is_some() && organic && !l_organic {
                    best = None;
                }
            }
            let mut num_atoms = 0u32;
            if self.use_atom_count {
                for &idx in frag {
                    num_atoms += 1;
                    if !self.count_heavy_atoms_only {
                        num_atoms += mol.total_num_hs(idx as u32, false)?;
                    }
                }
                if best.is_some() && num_atoms < l_num_atoms {
                    continue;
                }
            }
            let mut weight = 0.0f64;
            for &idx in frag {
                let atom = mol.atom(idx as u32);
                // it's not important to be perfect here
                weight += f64::from(100 * atom.atomic_num() + atom.isotope())
                    - f64::from(atom.formal_charge()) * 0.1;
                if !self.count_heavy_atoms_only {
                    weight += f64::from(mol.total_num_hs(idx as u32, false)?);
                }
            }
            let same_size = !self.use_atom_count || num_atoms == l_num_atoms;
            if best.is_some() && same_size && weight < l_weight {
                continue;
            }
            if best.is_some() && same_size && weight == l_weight && smiles > l_smiles {
                continue;
            }
            rd_log_info!("New largest fragment: {smiles} ({num_atoms})");
            best = Some(fidx);
            l_smiles = smiles;
            l_num_atoms = num_atoms;
            l_weight = weight;
            l_organic = organic;
        }
        mol.begin_batch_edit()?;
        for (fi, frag) in frags.iter().enumerate() {
            if Some(fi) == best {
                continue;
            }
            for &i in frag {
                mol.remove_atom(i as u32)?;
            }
        }
        mol.commit_batch_edit()
    }
}
