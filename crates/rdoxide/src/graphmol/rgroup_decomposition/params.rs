//! Parameters of the R-group decomposition (upstream
//! `RGroupDecomposition/RGroupDecompParams.{h,cpp}` and
//! `RGroupDecompJSONParsers.cpp`).

use std::collections::BTreeSet;
use std::sync::Arc;

use super::utils::{
    Labelling, RLABEL, UNLABELED_CORE_ATTACHMENT, clear_input_labels, int_prop, set_int_prop,
    set_label,
};
use crate::graphmol::adjust_query::{AdjustQueryParameters, adjust_query_properties_in_place};
use crate::graphmol::atom::{Atom, ChiralType};
use crate::graphmol::bond::{Bond, BondStereo, BondType};
use crate::graphmol::fmcs::{
    AtomTyper, MCSAtomCompareParameters, MCSParameters, find_mcs, mcs_atom_compare_elements,
};
use crate::graphmol::mol_ops::add_hs::set_terminal_atom_coords;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::queries::make_bond_null_query;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::ptree::PTree;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_warning};

/// `RGroupLabels`: where the R-labels of the cores are encoded.
pub mod rgroup_labels {
    pub const ISOTOPE_LABELS: u32 = 0x01;
    pub const ATOM_MAP_LABELS: u32 = 0x02;
    pub const ATOM_INDEX_LABELS: u32 = 0x04;
    pub const RELABEL_DUPLICATE_LABELS: u32 = 0x08;
    pub const MDL_RGROUP_LABELS: u32 = 0x10;
    /// R groups, but they will get relabelled
    pub const DUMMY_ATOM_LABELS: u32 = 0x20;
    pub const AUTO_DETECT: u32 = 0xFF;

    /// (name, value) pairs, for the Python enums and the JSON parser
    pub const NAMES: [(&str, u32); 7] = [
        ("IsotopeLabels", ISOTOPE_LABELS),
        ("AtomMapLabels", ATOM_MAP_LABELS),
        ("AtomIndexLabels", ATOM_INDEX_LABELS),
        ("RelabelDuplicateLabels", RELABEL_DUPLICATE_LABELS),
        ("MDLRGroupLabels", MDL_RGROUP_LABELS),
        ("DummyAtomLabels", DUMMY_ATOM_LABELS),
        ("AutoDetect", AUTO_DETECT),
    ];
}

/// `RGroupMatching`: how the R groups of the molecules are assigned.
pub mod rgroup_matching {
    pub const GREEDY: u32 = 0x01;
    pub const GREEDY_CHUNKS: u32 = 0x02;
    /// not really useful for large sets
    pub const EXHAUSTIVE: u32 = 0x04;
    pub const NO_SYMMETRIZATION: u32 = 0x08;
    pub const GA: u32 = 0x10;

    /// (name, value) pairs, for the Python enums and the JSON parser
    pub const NAMES: [(&str, u32); 5] = [
        ("Greedy", GREEDY),
        ("GreedyChunks", GREEDY_CHUNKS),
        ("Exhaustive", EXHAUSTIVE),
        ("NoSymmetrization", NO_SYMMETRIZATION),
        ("GA", GA),
    ];
}

/// `RGroupLabelling`: where the R-labels are stored on the output.
pub mod rgroup_labelling {
    pub const ATOM_MAP: u32 = 0x01;
    pub const ISOTOPE: u32 = 0x02;
    pub const MDL_RGROUP: u32 = 0x04;

    /// (name, value) pairs, for the Python enums and the JSON parser
    pub const NAMES: [(&str, u32); 3] = [
        ("AtomMap", ATOM_MAP),
        ("Isotope", ISOTOPE),
        ("MDLRGroup", MDL_RGROUP),
    ];
}

/// `RGroupCoreAlignment`
pub mod rgroup_core_alignment {
    pub const NO_ALIGNMENT: u32 = 0x0;
    pub const MCS: u32 = 0x01;

    /// (name, value) pairs, for the Python enums and the JSON parser
    pub const NAMES: [(&str, u32); 2] = [("NoAlignment", NO_ALIGNMENT), ("MCS", MCS)];
}

/// `RGroupScore`
pub mod rgroup_score {
    pub const MATCH: u32 = 0x1;
    pub const FINGERPRINT_VARIANCE: u32 = 0x4;

    /// (name, value) pairs, for the Python enums and the JSON parser
    pub const NAMES: [(&str, u32); 2] = [
        ("Match", MATCH),
        ("FingerprintVariance", FINGERPRINT_VARIANCE),
    ];
}

use rgroup_labels::{
    ATOM_INDEX_LABELS, ATOM_MAP_LABELS, AUTO_DETECT, DUMMY_ATOM_LABELS, ISOTOPE_LABELS,
    MDL_RGROUP_LABELS, RELABEL_DUPLICATE_LABELS,
};

/// `RGroupDecompositionParameters`
#[derive(Clone, Debug)]
pub struct RGroupDecompositionParameters {
    pub labels: u32,
    pub matching_strategy: u32,
    pub score_method: u32,
    pub rgroup_labelling: u32,
    pub alignment: u32,
    pub chunk_size: u32,
    pub only_match_at_rgroups: bool,
    pub remove_all_hydrogen_rgroups: bool,
    pub remove_all_hydrogen_rgroups_and_labels: bool,
    pub remove_hydrogens_post_match: bool,
    pub allow_non_terminal_rgroups: bool,
    pub allow_multiple_rgroups_on_unlabelled: bool,
    pub allow_multiple_cores_in_same_mol: bool,
    pub do_tautomers: bool,
    pub do_enumeration: bool,
    pub include_target_mol_in_results: bool,
    /// timeout in seconds, <= 0 means no timeout
    pub timeout: f64,
    pub ga_population_size: i32,
    pub ga_maximum_operations: i32,
    pub ga_number_operations_without_improvement: i32,
    pub ga_random_seed: i32,
    pub ga_number_runs: i32,
    pub ga_parallel_runs: bool,
    pub substructmatch_params: SubstructMatchParameters,
    /// Internal state: the labels of core atoms without R-labels count down
    /// from here, over all the cores prepared with these parameters.
    pub index_offset: i32,
}

impl Default for RGroupDecompositionParameters {
    fn default() -> Self {
        RGroupDecompositionParameters {
            labels: AUTO_DETECT,
            matching_strategy: rgroup_matching::GREEDY_CHUNKS,
            score_method: rgroup_score::MATCH,
            rgroup_labelling: rgroup_labelling::ATOM_MAP | rgroup_labelling::MDL_RGROUP,
            alignment: rgroup_core_alignment::MCS,
            chunk_size: 5,
            only_match_at_rgroups: false,
            remove_all_hydrogen_rgroups: true,
            remove_all_hydrogen_rgroups_and_labels: true,
            remove_hydrogens_post_match: true,
            allow_non_terminal_rgroups: false,
            allow_multiple_rgroups_on_unlabelled: false,
            allow_multiple_cores_in_same_mol: false,
            do_tautomers: false,
            do_enumeration: false,
            include_target_mol_in_results: false,
            timeout: -1.0,
            ga_population_size: -1,
            ga_maximum_operations: -1,
            ga_number_operations_without_improvement: -1,
            ga_random_seed: -1,
            ga_number_runs: 1,
            ga_parallel_runs: false,
            substructmatch_params: SubstructMatchParameters {
                use_chirality: true,
                ..Default::default()
            },
            index_offset: -1,
        }
    }
}

fn has_label(mol: &ROMol, idx: u32, auto_labels: u32) -> bool {
    let atom = mol.atom(idx);
    let mut has = false;
    if auto_labels & MDL_RGROUP_LABELS != 0 {
        has |= atom.has_prop(cp::_MolFileRLabel);
    }
    if auto_labels & ISOTOPE_LABELS != 0 {
        has |= atom.isotope() > 0;
    }
    if auto_labels & ATOM_MAP_LABELS != 0 {
        has |= atom.atom_map_num() > 0;
    }
    if auto_labels & DUMMY_ATOM_LABELS != 0 {
        has |= atom.atomic_num() == 0 && mol.atom_degree(idx) == 1;
    }
    if let Some(label) = int_prop(atom, RLABEL) {
        has |= label > 0;
    }
    has
}

/// When comparing atoms for the MCS overlay, a core atom with a labelled
/// R group must not be overlaid on one without. The label may be on a
/// terminal dummy atom bonded to the core atom.
fn has_attached_labels(mol: &ROMol, idx: u32, auto_labels: u32) -> bool {
    mol.atom_neighbors(idx).any(|n| {
        mol.atom(n).atomic_num() == 0 && mol.atom_degree(n) == 1 && has_label(mol, n, auto_labels)
    })
}

/// `rgdAtomCompare`: the MCS atom comparator used to align cores.
fn rgd_atom_compare(
    p: &MCSAtomCompareParameters,
    mol1: &ROMol,
    atom1: u32,
    mol2: &ROMol,
    atom2: u32,
    auto_labels: u32,
) -> Result<bool> {
    if !mcs_atom_compare_elements(p, mol1, atom1, mol2, atom2)? {
        return Ok(false);
    }
    let l1 = has_label(mol1, atom1, auto_labels) || has_attached_labels(mol1, atom1, auto_labels);
    let l2 = has_label(mol2, atom2, auto_labels) || has_attached_labels(mol2, atom2, auto_labels);
    Ok(l1 == l2)
}

impl RGroupDecompositionParameters {
    /// `autoGetLabels`: which kind of labels `core` carries.
    pub fn auto_get_labels(&self, core: &RWMol) -> u32 {
        let mut auto_labels = 0;
        if !self.only_match_at_rgroups {
            auto_labels = ATOM_INDEX_LABELS;
        }
        let (mut mdl, mut atom_map, mut isotopes, mut dummies) = (false, false, false, false);
        for (idx, atom) in core.atoms().iter().enumerate() {
            isotopes |= atom.isotope() != 0;
            atom_map |= atom.atom_map_num() != 0;
            mdl |= atom.has_prop(cp::_MolFileRLabel);
            dummies |= atom.atomic_num() == 0 && core.atom_degree(idx as u32) == 1;
        }
        if mdl {
            auto_labels | MDL_RGROUP_LABELS
        } else if atom_map {
            auto_labels | ATOM_MAP_LABELS
        } else if isotopes {
            auto_labels | ISOTOPE_LABELS
        } else if dummies {
            auto_labels | DUMMY_ATOM_LABELS
        } else {
            auto_labels
        }
    }

    /// `checkNonTerminal`
    fn check_non_terminal(&self, degree: u32) -> Result<()> {
        if self.allow_non_terminal_rgroups || degree == 1 {
            return Ok(());
        }
        rd_log_warning!(
            "Non terminal R group defined.  To allow set allowNonTerminalRGroups in RGroupDecompositionParameters"
        );
        Err(RDKitError::Value("Non terminal R group defined.".into()))
    }

    /// Copies the R-labels of `align_core` onto the atoms of `core` that the
    /// MCS of both overlays. Returns the next free label.
    fn align_core_labels(
        &self,
        core: &mut RWMol,
        align_core: &RWMol,
        auto_labels: u32,
    ) -> Result<i32> {
        let mut max_label = 1;
        let mut mcs_params = MCSParameters::default();
        if auto_labels != AUTO_DETECT {
            mcs_params.atom_typer = AtomTyper::Custom(Arc::new(move |p, m1, a1, m2, a2| {
                rgd_atom_compare(p, m1, a1, m2, a2, auto_labels)
            }));
        }
        let res = find_mcs(&[&*core, align_core], Some(&mcs_params))?;
        if !res.is_completed() {
            return Ok(max_label);
        }
        let Some(m) = res.query_mol else {
            return Ok(max_label);
        };
        let first_match = |mol: &RWMol| -> Result<Option<Vec<(i32, i32)>>> {
            let ps = SubstructMatchParameters {
                max_matches: 1,
                ..Default::default()
            };
            Ok(substruct_match(mol, &m, &ps)?.into_iter().next())
        };
        // the matches are (query atom, target atom) pairs
        let match1 = first_match(core)?;
        let match2 = first_match(align_core)?;
        let (Some(match1), Some(match2)) = (match1, match2) else {
            return Ok(max_label);
        };
        if match1.len() != match2.len() {
            return Err(RDKitError::invariant(
                "Matches should be the same size in prepareCore",
            ));
        }
        for (&(q1, core_idx), &(q2, align_idx)) in match1.iter().zip(&match2) {
            if q1 != q2 {
                return Err(RDKitError::invariant("query atoms aren't the same"));
            }
            let align_label = int_prop(align_core.atom(align_idx as u32), RLABEL);
            let core_atom = core.atom_mut(core_idx as u32);
            let _ = core_atom.set_atom_map_num(0, false);
            if core_atom.has_prop(cp::_MolFileRLabel) {
                core_atom.clear_prop(cp::_MolFileRLabel);
                core_atom.set_isotope(0);
            }
            if let Some(rlabel) = align_label {
                max_label = max_label.max(rlabel + 1);
                set_int_prop(core_atom, RLABEL, rlabel);
            }
        }
        Ok(max_label)
    }

    /// `prepareCore`: finds and normalizes the R-labels of `core`, aligning
    /// them to those of `align_core` when requested. `false` when the core
    /// has no R groups to find.
    pub fn prepare_core(&mut self, core: &mut RWMol, align_core: Option<&RWMol>) -> Result<bool> {
        let relabel = self.labels & RELABEL_DUPLICATE_LABELS != 0;
        let mut auto_labels = self.labels;
        if self.labels == AUTO_DETECT {
            auto_labels = self.auto_get_labels(core);
            if auto_labels == 0 {
                rd_log_warning!(
                    "RGroupDecomposition auto detect found no rgroups and onlyMatchAtRgroups is set to true"
                );
                return Ok(false);
            }
        } else if !self.only_match_at_rgroups {
            auto_labels |= ATOM_INDEX_LABELS;
        }

        if !self.substructmatch_params.use_chirality {
            for idx in 0..core.num_atoms() {
                core.atom_mut(idx).set_chiral_tag(ChiralType::Unspecified);
            }
            for idx in 0..core.num_bonds() {
                core.bond_mut(idx).set_stereo(BondStereo::StereoNone)?;
            }
        }
        if !self.substructmatch_params.use_chirality
            || !self.substructmatch_params.use_enhanced_stereo
        {
            core.set_stereo_groups(Vec::new())?;
        }

        let mut max_label = 1;
        if let Some(align_core) = align_core {
            if !self.only_match_at_rgroups
                && self.alignment & super::params::rgroup_core_alignment::MCS != 0
            {
                max_label = self.align_core_labels(core, align_core, auto_labels)?;
            }
        }

        let mut found_labels = BTreeSet::new();
        let mut atom_to_label = Vec::new();
        for idx in 0..core.num_atoms() {
            let degree = core.atom_degree(idx);
            let atom = core.atom_mut(idx);
            let mut found = false;
            if let Some(label) = int_prop(atom, RLABEL) {
                if set_label(
                    atom,
                    label,
                    &mut found_labels,
                    &mut max_label,
                    relabel,
                    Labelling::InternalLabels,
                )? {
                    found = true;
                }
            }
            if atom.atomic_num() == 0 {
                if !found && auto_labels & MDL_RGROUP_LABELS != 0 {
                    if let Some(rgroup) =
                        atom.get_prop(cp::_MolFileRLabel).and_then(RDValue::as_uint)
                    {
                        if set_label(
                            atom,
                            rgroup as i32,
                            &mut found_labels,
                            &mut max_label,
                            relabel,
                            Labelling::RGroupLabels,
                        )? {
                            found = true;
                            self.check_non_terminal(degree)?;
                        }
                    }
                }
                if !found && auto_labels & ISOTOPE_LABELS != 0 && atom.isotope() > 0 {
                    let isotope = atom.isotope() as i32;
                    if set_label(
                        atom,
                        isotope,
                        &mut found_labels,
                        &mut max_label,
                        relabel,
                        Labelling::IsotopeLabels,
                    )? {
                        self.check_non_terminal(degree)?;
                        found = true;
                    }
                }
                if !found && auto_labels & ATOM_MAP_LABELS != 0 && atom.atom_map_num() > 0 {
                    let map_num = atom.atom_map_num();
                    if set_label(
                        atom,
                        map_num,
                        &mut found_labels,
                        &mut max_label,
                        relabel,
                        Labelling::AtomMapLabels,
                    )? {
                        self.check_non_terminal(degree)?;
                        found = true;
                    }
                }
                if !found
                    && auto_labels & DUMMY_ATOM_LABELS != 0
                    && degree == 1
                    && !atom.has_prop(UNLABELED_CORE_ATTACHMENT)
                {
                    let start = max_label;
                    if set_label(
                        atom,
                        start,
                        &mut found_labels,
                        &mut max_label,
                        true,
                        Labelling::DummyLabels,
                    )? {
                        found = true;
                    }
                }
            }
            if !found && auto_labels & ATOM_INDEX_LABELS != 0 {
                set_label(
                    atom,
                    self.index_offset - idx as i32,
                    &mut found_labels,
                    &mut max_label,
                    relabel,
                    Labelling::IndexLabels,
                )?;
            }
            clear_input_labels(atom);
            if let Some(rlabel) = int_prop(atom, RLABEL) {
                atom_to_label.push((idx, rlabel));
            }
        }
        self.index_offset -= core.num_atoms() as i32;

        let adjust = AdjustQueryParameters {
            make_dummies_queries: true,
            adjust_degree: false,
            ..Default::default()
        };
        adjust_query_properties_in_place(core, &adjust)?;
        for (idx, label) in atom_to_label {
            set_int_prop(core.atom_mut(idx), RLABEL, label);
        }
        Ok(true)
    }

    /// `addDummyAtomsToUnlabelledCoreAtoms`: when several R groups are
    /// allowed on unlabelled core atoms, every free valence of such an atom
    /// gets a dummy atom to hold one.
    pub fn add_dummy_atoms_to_unlabelled_core_atoms(&self, core: &mut RWMol) -> Result<()> {
        if !self.allow_multiple_rgroups_on_unlabelled {
            return Ok(());
        }
        let unlabelled: Vec<u32> = (0..core.num_atoms())
            .filter(|&i| core.atom(i).atomic_num() != 1 && !has_label(core, i, self.labels))
            .collect();
        for atom_idx in unlabelled {
            core.calc_implicit_valence(atom_idx, false)?;
            let max_dummies = 4 - core.atom_degree(atom_idx) as i32;
            let dummies_to_add = if core.atom(atom_idx).atomic_num() == 0 {
                max_dummies
            } else {
                let mut bond_order = 0.0;
                for bond in core.atom_bonds(atom_idx) {
                    let mut contrib = bond.valence_contrib(core, atom_idx);
                    if contrib == 0.0 && bond.has_query() {
                        contrib = 1.0;
                    }
                    bond_order += contrib;
                }
                let valences =
                    PeriodicTable::get().get_valence_list(core.atom(atom_idx).atomic_num())?;
                let valence = valences.iter().copied().max().unwrap_or(0);
                (valence - (bond_order + 0.51) as i32).min(max_dummies)
            };
            let mut new_indices = Vec::new();
            for _ in 0..dummies_to_add {
                let mut dummy = Atom::new(0);
                dummy.dict_mut().set_val(UNLABELED_CORE_ATTACHMENT, true);
                let new_idx = core.add_atom(dummy, false);
                new_indices.push(new_idx);
                let mut qb = Bond::new(BondType::Unspecified);
                qb.set_query(Some(make_bond_null_query()));
                qb.set_begin_atom_idx_unchecked(atom_idx);
                qb.set_end_atom_idx_unchecked(new_idx);
                core.add_bond_object(qb)?;
                core.update_atom_property_cache(new_idx, true)?;
            }
            core.update_atom_property_cache(atom_idx, false)?;
            if core.num_conformers() > 0 {
                for new_idx in new_indices {
                    set_terminal_atom_coords(core, new_idx, atom_idx)?;
                }
            }
        }
        Ok(())
    }

    /// `updateRGroupDecompositionParametersFromJSON`
    pub fn update_from_json(&mut self, details_json: &str) -> Result<()> {
        if details_json.is_empty() {
            return Ok(());
        }
        let pt = PTree::from_json(details_json)?;
        let named = |key: &str, names: &[(&str, u32)]| -> Option<u32> {
            let value = pt.get_str(key)?;
            names.iter().find(|(n, _)| *n == value).map(|&(_, v)| v)
        };
        if let Some(v) = named("labels", &rgroup_labels::NAMES) {
            self.labels = v;
        }
        if let Some(v) = named("matchingStrategy", &rgroup_matching::NAMES) {
            self.matching_strategy = v;
        }
        if let Some(v) = named("scoreMethod", &rgroup_score::NAMES) {
            self.score_method = v;
        }
        if let Some(v) = named("rgroupLabelling", &rgroup_labelling::NAMES) {
            self.rgroup_labelling = v;
        }
        if let Some(v) = named("alignment", &rgroup_core_alignment::NAMES) {
            self.alignment = v;
        }
        self.chunk_size = pt.get_u32("chunkSize", self.chunk_size);
        self.only_match_at_rgroups = pt.get_bool("onlyMatchAtRGroups", self.only_match_at_rgroups);
        self.remove_all_hydrogen_rgroups =
            pt.get_bool("removeAllHydrogenRGroups", self.remove_all_hydrogen_rgroups);
        self.remove_all_hydrogen_rgroups_and_labels = pt.get_bool(
            "removeAllHydrogenRGroupsAndLabels",
            self.remove_all_hydrogen_rgroups_and_labels,
        );
        self.remove_hydrogens_post_match =
            pt.get_bool("removeHydrogensPostMatch", self.remove_hydrogens_post_match);
        self.allow_non_terminal_rgroups =
            pt.get_bool("allowNonTerminalRGroups", self.allow_non_terminal_rgroups);
        self.allow_multiple_rgroups_on_unlabelled = pt.get_bool(
            "allowMultipleRGroupsOnUnlabelled",
            self.allow_multiple_rgroups_on_unlabelled,
        );
        self.do_tautomers = pt.get_bool("doTautomers", self.do_tautomers);
        self.do_enumeration = pt.get_bool("doEnumeration", self.do_enumeration);
        self.include_target_mol_in_results = pt.get_bool(
            "includeTargetMolInResults",
            self.include_target_mol_in_results,
        );
        if let Some(t) = pt.get_str("timeout") {
            self.timeout = t
                .parse()
                .map_err(|_| RDKitError::Value(format!("bad timeout value '{t}'")))?;
        }
        Ok(())
    }
}
