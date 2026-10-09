//! `molzip`: joins molecules (or the pieces of one molecule) at labelled
//! dummy atoms, keeping atom and bond stereochemistry.
//!
//! Ported from RDKit 2026.09.1 `GraphMol/ChemTransforms/MolZip.cpp`.
//!
//! Copyright (C) 2020-2025 Brian Kelley and other RDKit contributors.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).

use std::collections::{BTreeMap, BTreeSet};

use super::combine_mols;
use crate::geometry::{Point3D, Transform3D};
use crate::graphmol::atom::ChiralType;
use crate::graphmol::bond::{BondDir, BondStereo, BondType};
use crate::graphmol::conformer::Conformer;
use crate::graphmol::depictor::{Compute2DCoordParameters, compute_2d_coords};
use crate::graphmol::mol_ops::frags::get_mol_frags_mapping;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles};
use crate::rdgeneral::common_properties as cp;
use crate::rdgeneral::{RDProps, RDValue};
use crate::{RDKitError, Result, rd_log_error, rd_log_warning};

/// How the dummy atoms that are zipped together are labelled.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum MolzipLabel {
    /// atom map numbers on dummy atoms
    #[default]
    AtomMapNumber,
    /// isotope labels on dummy atoms
    Isotope,
    /// the output of `fragmentOnBonds` (the labels are atom indices)
    FragmentOnBonds,
    /// the atoms whose symbols are in [`MolzipParams::atom_symbols`]
    AtomType,
    /// an unsigned integer atom property
    AtomProperty,
}

impl MolzipLabel {
    pub fn from_u8(v: u8) -> Option<Self> {
        Some(match v {
            0 => Self::AtomMapNumber,
            1 => Self::Isotope,
            2 => Self::FragmentOnBonds,
            3 => Self::AtomType,
            4 => Self::AtomProperty,
            _ => return None,
        })
    }
}

/// `MolzipParams`
#[derive(Clone, Debug)]
pub struct MolzipParams {
    pub label: MolzipLabel,
    pub atom_symbols: Vec<String>,
    pub atom_property: String,
    pub enforce_valence_rules: bool,
    /// generate 2D coordinates for the zipped molecule and the input
    /// fragments ([`molzip_fragments`] only)
    pub generate_coordinates: bool,
    /// align fragments with coordinates along the connection vectors
    pub align_coordinates: bool,
}

impl Default for MolzipParams {
    fn default() -> Self {
        MolzipParams {
            label: MolzipLabel::AtomMapNumber,
            atom_symbols: Vec::new(),
            atom_property: String::new(),
            enforce_valence_rules: true,
            generate_coordinates: false,
            align_coordinates: false,
        }
    }
}

const NOLABEL: u32 = u32::MAX;
const SOURCE_FILE: &str = "Code/GraphMol/ChemTransforms/MolZip.cpp";

/// A failed `CHECK_INVARIANT` of `MolZip.cpp` (`line` is the last line of the
/// statement).
fn invariant(line: u32, expr: &str, message: impl std::fmt::Display) -> RDKitError {
    RDKitError::invariant_at(message, expr, SOURCE_FILE, line)
}

/// A failed `PRECONDITION` of `MolZip.cpp`.
fn precondition(line: u32, expr: &str, message: &str) -> RDKitError {
    RDKitError::precondition_at(message, expr, SOURCE_FILE, line)
}
/// Atom property holding the position of an atom in a set of fragments to be
/// zipped with generated coordinates.
const INDEX_PROP: &str = "__zipIndex";
const MOLZIP_PREFIX: &str = "__molzip";

/// The label of an atom, or [`NOLABEL`].
fn get_label(mol: &ROMol, atom: u32, params: &MolzipParams) -> Result<u32> {
    let a = mol.atom(atom);
    Ok(match params.label {
        MolzipLabel::AtomMapNumber if a.atomic_num() == 0 => {
            let mapno = a.atom_map_num();
            if mapno != 0 { mapno as u32 } else { NOLABEL }
        }
        MolzipLabel::Isotope if a.atomic_num() == 0 => {
            let iso = a.isotope();
            if iso != 0 { iso } else { NOLABEL }
        }
        MolzipLabel::AtomMapNumber | MolzipLabel::Isotope => NOLABEL,
        MolzipLabel::AtomType => params
            .atom_symbols
            .iter()
            .position(|s| *s == a.symbol())
            .map_or(NOLABEL, |i| i as u32),
        MolzipLabel::FragmentOnBonds => {
            return Err(invariant(
                60,
                "0",
                "FragmentOnBonds is not an atom label, it is an atom index",
            ));
        }
        MolzipLabel::AtomProperty => a
            .get_prop(&params.atom_property)
            .and_then(RDValue::as_uint)
            .unwrap_or(NOLABEL),
    })
}

/// The only neighbor of `atom`, if it has exactly one.
fn get_other_atom(mol: &ROMol, atom: u32) -> Option<u32> {
    if mol.atom_degree(atom) != 1 {
        return None;
    }
    mol.atom_neighbors(atom).next()
}

/// `num_swaps_to_interconvert`: the number of swaps between the identity and
/// the permutation `orders`.
fn num_swaps_to_interconvert(orders: &[u32]) -> Result<u32> {
    let mut nswaps = 0;
    let mut seen = vec![false; orders.len()];
    for i in 0..orders.len() {
        if !seen[i] {
            let mut j = i;
            while orders[j] as usize != i {
                j = orders[j] as usize;
                if j >= orders.len() {
                    return Err(invariant(
                        95,
                        "j < orders.size()",
                        "molzip: bond index outside of number of bonds for atom",
                    ));
                }
                seen[j] = true;
                nswaps += 1;
            }
        }
    }
    Ok(nswaps)
}

/// Bookkeeping for one bond to make between two atoms (`ZipBond`).
#[derive(Clone, Debug, Default)]
struct ZipBond {
    /// atom being bonded
    a: Option<u32>,
    /// labelled atom, i.e. `[*:1]-C` bonds the C to something
    a_dummy: Option<u32>,
    b: Option<u32>,
    b_dummy: Option<u32>,
    /// link bonds have six atoms, `[*:a][*b].[*:a]C.[*:b]D`: four dummies and
    /// two atoms
    a_link: Option<u32>,
    b_link: Option<u32>,
    is_linker: bool,
    linker_bond_type: BondType,
}

impl ZipBond {
    /// The four atoms of a complete bond.
    fn atoms(&self) -> Option<(u32, u32, u32, u32)> {
        Some((self.a?, self.b?, self.a_dummy?, self.b_dummy?))
    }
}

/// Chirality information noted before bonding and used to restore it after.
#[derive(Default)]
struct StereoMarks {
    /// chiral atom -> neighbor -> position of the neighbor around the atom
    /// (a new neighbor takes the position of the dummy it replaces)
    nbr_order: BTreeMap<u32, BTreeMap<u32, i32>>,
    /// stereo bond -> (stereo atoms, stereo), with a dummy replaced
    bond_stereo: BTreeMap<u32, (Vec<u32>, BondStereo)>,
}

impl StereoMarks {
    /// Back up the original chirality (`ZipBond::mark`): the neighbors' order
    /// around `chiral_atom`, which `new_atom` will take over from
    /// `dummy_atom`, and the stereo atoms of its bonds.
    fn mark(
        &mut self,
        mol: &ROMol,
        chiral_atom: u32,
        dummy_atom: u32,
        new_atom: u32,
    ) -> Result<()> {
        if mol.atom(chiral_atom).chiral_tag() != ChiralType::Unspecified {
            let order: Vec<(u32, i32)> = mol
                .atom_neighbors(chiral_atom)
                .enumerate()
                .map(|(i, n)| (n, i as i32))
                .collect();
            let dummy_order = order
                .iter()
                .find(|(n, _)| *n == dummy_atom)
                .map(|&(_, o)| o)
                .ok_or_else(|| RDKitError::Key(format!("__molzip_mark_{chiral_atom}")))?;
            // atoms can be bonded more than once, so the marks accumulate
            let marks = self.nbr_order.entry(chiral_atom).or_default();
            marks.extend(order);
            marks.insert(new_atom, dummy_order);
        }
        for bond in mol.atom_bonds(chiral_atom) {
            if bond.stereo() == BondStereo::StereoNone {
                continue;
            }
            // a bond cut at both ends is marked twice: build on the first mark
            // so that both dummies end up replaced
            let current: Vec<u32> = match self.bond_stereo.get(&bond.idx()) {
                Some((atoms, _)) => atoms.clone(),
                None => bond.stereo_atoms().iter().map(|&a| a as u32).collect(),
            };
            if current.contains(&dummy_atom) {
                let atoms = current
                    .iter()
                    .map(|&a| if a == dummy_atom { new_atom } else { a })
                    .collect();
                self.bond_stereo.insert(bond.idx(), (atoms, bond.stereo()));
            }
        }
        Ok(())
    }
}

impl ZipBond {
    /// `mark_chirality`; must be called before the bond is made.
    fn mark_chirality(&self, mol: &ROMol, marks: &mut StereoMarks) -> Result<()> {
        let a = self
            .a
            .ok_or_else(|| precondition(122, "a", "Must have a begin atom to bond"))?;
        let b = self
            .b
            .ok_or_else(|| precondition(123, "b", "Must have an end atom to bond"))?;
        let a_dummy = self
            .a_dummy
            .ok_or_else(|| precondition(124, "a_dummy", "Must have a begin dummy atom"))?;
        let b_dummy = self
            .b_dummy
            .ok_or_else(|| precondition(125, "b_dummy", "Must have an end dummy atom"))?;
        marks.mark(mol, a, a_dummy, b)?;
        marks.mark(mol, b, b_dummy, a)
    }

    /// Bonds a and b (`ZipBond::bond`); `false` when the labelling is
    /// incomplete or `true` once the bond exists. Records the dummies to
    /// remove in `used`.
    fn bond(
        &self,
        mol: &mut RWMol,
        params: &MolzipParams,
        used: &mut BTreeSet<u32>,
    ) -> Result<bool> {
        let (Some(a), Some(b), Some(a_dummy), Some(b_dummy)) =
            (self.a, self.b, self.a_dummy, self.b_dummy)
        else {
            rd_log_warning!("Incomplete atom labelling, cannot make bond\n");
            return Ok(false);
        };
        // fragment on bonds allows multiple links to the same atom, i.e.
        // C.[1C].[1C]; otherwise it's an error
        let existing = mol.bond_idx_between(a, b);
        if params.label != MolzipLabel::FragmentOnBonds && existing.is_some() {
            return Err(invariant(
                148,
                "params.label == MolzipLabel::FragmentOnBonds || !a->getOwningMol().getBondBetweenAtoms(a->getIdx(), b->getIdx())",
                "molzip: zipped Bond already exists, perhaps labels are duplicated",
            ));
        }
        if existing.is_none() {
            if self.is_linker {
                // the easy bit: just link a and b, and schedule the dummies
                // for deletion
                let (Some(a_link), Some(b_link)) = (self.a_link, self.b_link) else {
                    return Err(invariant(
                        158,
                        "a && b && a_dummy && b_dummy && a_link && b_link",
                        "molzip: Link Bond is missing one or more labelled atoms",
                    ));
                };
                mol.add_bond(a, b, self.linker_bond_type)?;
                used.insert(a_link);
                used.insert(b_link);
            } else {
                self.bond_with_directions(mol, (a, a_dummy), (b, b_dummy))?;
            }
        }
        used.insert(a_dummy);
        used.insert(b_dummy);
        Ok(true)
    }

    /// Bonds `a` and `b`, taking the bond type and wedging of the bonds to
    /// the dummies replaced.
    fn bond_with_directions(
        &self,
        mol: &mut RWMol,
        (a, a_dummy): (u32, u32),
        (b, b_dummy): (u32, u32),
    ) -> Result<()> {
        let dummy_bond = |atom: u32, dummy: u32, line: u32, message: &str| -> Result<_> {
            let idx = mol
                .bond_idx_between(atom, dummy)
                .ok_or_else(|| invariant(line, "bnd != nullptr", message))?;
            let bond = mol.bond(idx);
            Ok((
                bond.bond_type(),
                bond.bond_dir(),
                bond.begin_atom_idx() == atom,
            ))
        };
        let (bond_type_a, mut bond_dir_a, a_is_start) = dummy_bond(
            a,
            a_dummy,
            167,
            "molzip: begin atom and specified dummy atom connection are not bonded.",
        )?;
        let (bond_type_b, mut bond_dir_b, b_is_start) = dummy_bond(
            b,
            b_dummy,
            175,
            "molzip: end atom and specified dummy connection atom are not bonded.",
        )?;

        // Fusion bond-dir logic table, < = wedge:
        //   a-* b-* => a-b
        //   a<* b-* => a<b
        //   a>* b-* => a>b
        //   a-* b>* => a<b
        //   a-* b<* => a>b
        let mut bond_dir = BondDir::None;
        let (mut start, mut end) = (a, b);
        if bond_dir_a != BondDir::None && bond_dir_b != BondDir::None {
            // consistent between the two bond orders? For fragment on bonds
            // a<* and b>* or a>* and b<*, where < is either a hash or wedge
            // bond but not both
            let consistent = bond_dir_a == bond_dir_b && a_is_start != b_is_start;
            if consistent {
                bond_dir_b = BondDir::None;
            } else {
                rd_log_warning!(
                    "inconsistent bond directions when merging fragments, ignoring...\n"
                );
                bond_dir_a = BondDir::None;
                bond_dir_b = BondDir::None;
            }
        }
        if bond_dir_a != BondDir::None {
            if !a_is_start {
                (start, end) = (b, a);
            }
            bond_dir = bond_dir_a;
        } else if bond_dir_b != BondDir::None {
            if b_is_start {
                (start, end) = (b, a);
            }
            bond_dir = bond_dir_b;
        }

        let bond_type = if bond_type_a != BondType::Single {
            bond_type_a
        } else if bond_type_b != BondType::Single {
            bond_type_b
        } else {
            BondType::Single
        };
        let n_bonds = mol.add_bond(start, end, bond_type)?;
        mol.bond_mut(n_bonds - 1).set_bond_dir(bond_dir);
        Ok(())
    }
}

/// Rotates and translates the fragment of `b` so that its connection vector
/// (`b` to `b_dummy`) lines up with the one of the molecule (`a` to `a_dummy`).
fn rotate_fragment_to_bond_vector(
    mol: &mut ROMol,
    (a, b, a_dummy, b_dummy): (u32, u32, u32, u32),
    fragment_for_atom: &[i32],
    atoms_in_fragment: &BTreeMap<i32, Vec<i32>>,
) -> Result<()> {
    if mol.num_conformers() == 0 {
        return Ok(());
    }
    let conf = mol.conformer_mut(-1)?;
    // Notation: Pmc molecule connection point (the atom that will be removed
    // from the molecule), Pma molecule attachment point (the atom to which
    // we'll form the bond), Psc/Psa the same for the sidechain; Um and Us are
    // the unit attachment vectors.
    let p_ma = conf.atom_pos(a)?;
    let p_sa = conf.atom_pos(b)?;
    let p_mc = conf.atom_pos(a_dummy)?;
    let p_sc = conf.atom_pos(b_dummy)?;

    let u_m = p_ma.direction_vector(&p_mc)?;
    // note the opposite direction here:
    let u_s = p_sc.direction_vector(&p_sa)?;

    let mut template = Transform3D::identity();
    template.set_translation(&p_ma);

    let cos_t = u_s.dot_product(&u_m);
    if cos_t.abs() < 1.0 {
        let sin_t = (1.0 - cos_t * cos_t).sqrt();
        let mut axis = u_s.cross_product(&u_m);
        axis.normalize()?;
        template.set_rotation_cs(cos_t, sin_t, &axis);
    } else {
        let mut normal = Point3D::new(1.0, 0.0, 0.0);
        if u_s.dot_product(&normal).abs() == 1.0 {
            normal = Point3D::new(0.0, 1.0, 0.0);
        }
        let axis = u_s.cross_product(&normal);
        template.set_rotation_cs(-1.0, 0.0, &axis);
    }

    // we use the second attachment vector to set the bond distance
    let mut translate = Transform3D::identity();
    translate.set_translation(&(p_sc * -1.0));
    template.mul_assign(&translate);

    // transform the atomic positions in the sidechain
    let frag_id = fragment_for_atom[b as usize];
    for &atom_idx in &atoms_in_fragment[&frag_id] {
        let pos = conf.atom_pos_mut(atom_idx as u32)?;
        *pos = template.transform_point(pos);
    }
    Ok(())
}

/// The zip bonds found in a molecule, by label, and what they refer to.
#[derive(Default)]
struct ZipPlan {
    /// label -> bond to make
    mappings: BTreeMap<u32, ZipBond>,
    /// atom -> labels of the bonds that involve it
    mappings_by_atom: BTreeMap<u32, Vec<u32>>,
    /// atoms to remove if they were used
    deletions: Vec<u32>,
    /// dummy index -> index of the atom that replaces it in the zipped
    /// molecule (atoms carrying [`INDEX_PROP`] only)
    attachment_mapping: BTreeMap<i32, i32>,
}

fn index_prop(mol: &ROMol, atom: u32) -> Option<i32> {
    mol.atom(atom)
        .get_prop(INDEX_PROP)
        .and_then(RDValue::as_int)
}

impl ZipPlan {
    fn record_attachment(&mut self, mol: &ROMol, dummy: Option<u32>, other: Option<u32>) {
        if let (Some(dummy_index), Some(other_index)) = (
            dummy.and_then(|d| index_prop(mol, d)),
            other.and_then(|o| index_prop(mol, o)),
        ) {
            self.attachment_mapping.insert(dummy_index, other_index);
        }
    }

    /// `MolzipLabel::FragmentOnBonds`: each dummy's isotope is the index of
    /// the atom it was cut from. `None` if the atom to bond cannot be found.
    fn from_fragment_on_bonds(mol: &ROMol) -> Result<Option<Self>> {
        let mut plan = ZipPlan::default();
        for atom in 0..mol.num_atoms() {
            if mol.atom(atom).atomic_num() != 0 {
                continue;
            }
            let molno = mol.atom(atom).isotope();
            let attached = get_other_atom(mol, atom).ok_or_else(|| {
                RDKitError::invariant("molzip: a dummy atom must have exactly one neighbor")
            })?;
            if molno >= mol.num_atoms() {
                // newmol->getAtomWithIdx(molno)
                return Err(RDKitError::range_at(
                    "idx",
                    molno,
                    mol.num_atoms(),
                    "Code/GraphMol/ROMol.cpp",
                    202,
                ));
            }
            let target = molno;
            let b_dummy = mol
                .atom_neighbors(target)
                .find(|&n| mol.atom(n).atomic_num() == 0 && mol.atom(n).isotope() == attached);
            let Some(b_dummy) = b_dummy else {
                rd_log_error!("Cannot find atom to bond using FragmentOnBond labelling\n");
                return Ok(None);
            };
            let bond = plan.mappings.entry(molno).or_default();
            bond.a = Some(attached);
            bond.a_dummy = Some(atom);
            bond.b = Some(target);
            bond.b_dummy = Some(b_dummy);
            plan.mappings_by_atom.entry(atom).or_default().push(molno);
            plan.deletions.push(atom);
            plan.record_attachment(mol, Some(atom), Some(target));
        }
        Ok(Some(plan))
    }

    /// Any other labelling: the dummies (or atoms) with the same label are
    /// joined, and `[*:1][*:2].[*:1]C.[*:2]S` style link bonds join their two
    /// attachments.
    fn from_labels(mol: &ROMol, params: &MolzipParams) -> Result<Self> {
        let mut plan = ZipPlan::default();
        // link bonds resolve to the same zip bond by using the lowest label,
        // so [*:1][*:2] sets linker_bonds[2] to the bond of label 1
        let mut linker_bonds: BTreeMap<u32, u32> = BTreeMap::new();
        for start_atom in 0..mol.num_atoms() {
            let mut atom = start_atom;
            let mut molno = get_label(mol, atom, params)?;
            if molno == NOLABEL {
                continue;
            }
            let attached_atom = get_other_atom(mol, atom);
            let mut attached_molno = match attached_atom {
                Some(a) => get_label(mol, a, params)?,
                None => NOLABEL,
            };
            if let Some(attached) = attached_atom.filter(|_| attached_molno != NOLABEL) {
                // we have a link bond; its bond comes first. Use the minimum
                // label to assign the bonds to link.
                let mut attached = attached;
                if molno > attached_molno {
                    std::mem::swap(&mut molno, &mut attached_molno);
                    std::mem::swap(&mut atom, &mut attached);
                }
                let link_bond = mol.bond_idx_between(atom, attached).ok_or_else(|| {
                    invariant(
                        497,
                        "link_bond",
                        format!(
                            "molzip: link bond with labels: {molno},{attached_molno} is missing"
                        ),
                    )
                })?;
                let bond_type = mol.bond(link_bond).bond_type();
                if let Some(bond) = plan.mappings.get(&molno) {
                    // we find this bond twice, so make sure it is set up
                    // correctly
                    if bond_type != bond.linker_bond_type {
                        return Err(invariant(
                            520,
                            "bondType == bond.linkerBondType",
                            format!(
                                "molzip: LINKER bond with labels: {molno},{attached_molno} has inconsistent bond types"
                            ),
                        ));
                    }
                    if !bond.is_linker {
                        return Err(invariant(
                            525,
                            "bond.isLinker",
                            format!(
                                "molzip: LINKER bond with labels: {molno},{attached_molno} found but not encountered first in the molecules to be zipped."
                            ),
                        ));
                    }
                    let same = (bond.a_link == Some(atom) && bond.b_link == Some(attached))
                        || (bond.b_link == Some(atom) && bond.a_link == Some(attached));
                    if !same {
                        return Err(invariant(
                            531,
                            "(bond.a_link == atom && bond.b_link == attached_atom) || (bond.b_link == atom && bond.a_link == attached_atom)",
                            format!(
                                "molzip: Linker Bond with labels {molno},{attached_molno} not setup correctly"
                            ),
                        ));
                    }
                } else {
                    if linker_bonds.contains_key(&attached_molno) {
                        return Err(invariant(
                            504,
                            "linkerBonds.find(attached_molno) == linkerBonds.end()",
                            format!(
                                "molzip: Linker attachment point with label: {attached_molno} found before linker bond"
                            ),
                        ));
                    }
                    linker_bonds.insert(attached_molno, molno);
                    plan.mappings.insert(
                        molno,
                        ZipBond {
                            is_linker: true,
                            linker_bond_type: bond_type,
                            a_link: Some(atom),
                            b_link: Some(attached),
                            ..ZipBond::default()
                        },
                    );
                    plan.deletions.push(atom);
                    plan.deletions.push(attached);
                }
            } else if !plan.mappings.contains_key(&molno) && !linker_bonds.contains_key(&molno) {
                // normal linkage C[*:1].S[*:1] links C and S
                let bond = plan.mappings.entry(molno).or_default();
                if bond.a.is_some() {
                    return Err(invariant(
                        541,
                        "!bond.a",
                        format!("molzip: bond info already setup for bgn atom with label:{molno}"),
                    ));
                }
                bond.a = attached_atom;
                bond.a_dummy = Some(atom);
            } else {
                let key = linker_bonds.get(&molno).copied().unwrap_or(molno);
                let bond = plan.mappings.entry(key).or_default();
                if bond.is_linker {
                    if bond.a.is_some() && bond.b.is_some() {
                        return Err(invariant(
                            552,
                            "!bond.a || !bond.b",
                            format!(
                                "molzip: Linker bond has multiple attachments for label: {molno}"
                            ),
                        ));
                    }
                    if bond.a.is_some() {
                        bond.b = attached_atom;
                        bond.b_dummy = Some(atom);
                    } else {
                        bond.a = attached_atom;
                        bond.a_dummy = Some(atom);
                    }
                    plan.deletions.push(atom);
                } else {
                    if bond.a.is_none() {
                        return Err(invariant(
                            566,
                            "bond.a",
                            format!(
                                "molzip: bond info not properly setup for bgn atom with label:{molno}"
                            ),
                        ));
                    }
                    if bond.b.is_some() {
                        return Err(invariant(
                            570,
                            "!bond.b",
                            format!(
                                "molzip: bond info already exists for end atom with label:{molno}"
                            ),
                        ));
                    }
                    bond.b = attached_atom;
                    bond.b_dummy = Some(atom);
                }
                let bond = bond.clone();
                if let Some(a) = bond.a {
                    plan.mappings_by_atom.entry(a).or_default().push(key);
                }
                plan.record_attachment(mol, bond.a_dummy, bond.b);
                plan.record_attachment(mol, bond.b_dummy, bond.a);
            }
            plan.deletions.push(atom);
        }
        Ok(plan)
    }
}

/// How the atoms and bonds of a molecule are renumbered by deleting atoms.
struct Renumbering {
    atom: Vec<Option<u32>>,
    bond: Vec<Option<u32>>,
}

impl Renumbering {
    /// For `mol` (before the deletion) and the atoms to `remove`.
    fn new(mol: &ROMol, remove: &BTreeSet<u32>) -> Self {
        let mut kept = 0;
        let atom = (0..mol.num_atoms())
            .map(|a| {
                (!remove.contains(&a)).then(|| {
                    kept += 1;
                    kept - 1
                })
            })
            .collect();
        let mut kept = 0;
        let bond =
            mol.bonds()
                .iter()
                .map(|b| {
                    (!remove.contains(&b.begin_atom_idx()) && !remove.contains(&b.end_atom_idx()))
                        .then(|| {
                            kept += 1;
                            kept - 1
                        })
                })
                .collect();
        Renumbering { atom, bond }
    }
}

/// Restores the chirality of `atom` after bonding (`ZipBond::restore`).
fn restore_atom_chirality(
    mol: &mut RWMol,
    atom: u32,
    original: u32,
    marks: &StereoMarks,
    inverse: &BTreeMap<u32, u32>,
) -> Result<()> {
    if mol.atom(atom).chiral_tag() == ChiralType::Unspecified {
        return Ok(());
    }
    let key_error = || RDKitError::Key(format!("__molzip_mark_{original}"));
    let order = marks.nbr_order.get(&original).ok_or_else(key_error)?;
    let orders = mol
        .atom_neighbors(atom)
        .map(|n| {
            let nbr_original = inverse.get(&n).ok_or_else(key_error)?;
            order
                .get(nbr_original)
                .map(|&o| o as u32)
                .ok_or_else(key_error)
        })
        .collect::<Result<Vec<_>>>()?;
    if num_swaps_to_interconvert(&orders)? % 2 == 1 {
        mol.atom_mut(atom).invert_chirality();
    }
    Ok(())
}

/// `molzip(a, b, params)` with the mapping of dummies to the atoms that
/// replaced them (for the atoms carrying `__zipIndex`). `None` when the
/// molecules cannot be zipped (a `FragmentOnBonds` partner is missing).
fn molzip_with_mapping(
    a: &ROMol,
    b: &ROMol,
    params: &MolzipParams,
) -> Result<Option<(RWMol, BTreeMap<i32, i32>)>> {
    let mut new_mol = if b.num_atoms() > 0 {
        combine_mols(a, b, Point3D::default())?
    } else {
        a.clone()
    };

    // aligning coordinates is quicker if we know which atoms are in which
    // fragments
    let mut fragment_for_atom = Vec::new();
    let mut atoms_in_fragment: BTreeMap<i32, Vec<i32>> = BTreeMap::new();
    if params.align_coordinates {
        fragment_for_atom = get_mol_frags_mapping(&new_mol).1;
        for (i, &frag) in fragment_for_atom.iter().enumerate() {
            atoms_in_fragment.entry(frag).or_default().push(i as i32);
        }
    }

    let plan = if params.label == MolzipLabel::FragmentOnBonds {
        match ZipPlan::from_fragment_on_bonds(&new_mol)? {
            Some(plan) => plan,
            None => return Ok(None),
        }
    } else {
        ZipPlan::from_labels(&new_mol, params)?
    };

    // mark the existing chirality so we can try and restore it later
    let mut marks = StereoMarks::default();
    for labels in plan.mappings_by_atom.values() {
        for label in labels {
            plan.mappings[label].mark_chirality(&new_mol, &mut marks)?;
        }
    }

    // make all the bonds
    let mut used = BTreeSet::new();
    let mut aligned = BTreeSet::new();
    for zip_bond in plan.mappings.values() {
        if zip_bond.bond(&mut new_mol, params, &mut used)? && params.align_coordinates {
            let atoms = zip_bond.atoms().expect("a made bond has all its atoms");
            // make sure we haven't aligned this fragment already
            if aligned.insert(fragment_for_atom[atoms.1 as usize]) {
                rotate_fragment_to_bond_vector(
                    &mut new_mol,
                    atoms,
                    &fragment_for_atom,
                    &atoms_in_fragment,
                )?;
            }
        }
    }

    // remove the used atoms
    let removed: BTreeSet<u32> = plan
        .deletions
        .iter()
        .copied()
        .filter(|d| used.contains(d))
        .collect();
    let renumbering = Renumbering::new(&new_mol, &removed);
    new_mol.begin_batch_edit()?;
    for &atom in &removed {
        new_mol.remove_atom(atom)?;
    }
    new_mol.commit_batch_edit()?;

    // try and restore the chirality now that we have new bonds
    let inverse: BTreeMap<u32, u32> = renumbering
        .atom
        .iter()
        .enumerate()
        .filter_map(|(old, new)| new.map(|n| (n, old as u32)))
        .collect();
    let mut already_checked = BTreeSet::new();
    for labels in plan.mappings_by_atom.values() {
        for label in labels {
            let (a, b, _, _) = plan.mappings[label]
                .atoms()
                .expect("marked bonds are complete");
            for atom in [a, b] {
                if already_checked.insert(atom) {
                    if let Some(new) = renumbering.atom[atom as usize] {
                        restore_atom_chirality(&mut new_mol, new, atom, &marks, &inverse)?;
                    }
                }
            }
        }
    }
    if !plan.mappings_by_atom.is_empty() {
        for (&bond, (atoms, stereo)) in &marks.bond_stereo {
            let new_atoms: Option<Vec<i32>> = atoms
                .iter()
                .map(|&a| renumbering.atom[a as usize].map(|n| n as i32))
                .collect();
            if let (Some(new_bond), Some(new_atoms)) = (renumbering.bond[bond as usize], new_atoms)
            {
                let b = new_mol.bond_mut(new_bond);
                *b.stereo_atoms_mut() = new_atoms;
                b.set_stereo(*stereo)?;
            }
        }
    }

    // remove all molzip tags
    for atom in 0..new_mol.num_atoms() {
        for prop in new_mol.atom(atom).get_prop_list(true, true) {
            if prop.starts_with(MOLZIP_PREFIX) {
                new_mol.atom_mut(atom).clear_prop(&prop);
            }
        }
    }
    for bond in 0..new_mol.num_bonds() {
        for prop in new_mol.bond(bond).get_prop_list(true, true) {
            if prop.starts_with(MOLZIP_PREFIX) {
                new_mol.bond_mut(bond).clear_prop(&prop);
            }
        }
    }
    new_mol.update_property_cache(params.enforce_valence_rules)?;
    new_mol.set_prop(cp::_StereochemDone, true, false)?;
    Ok(Some((new_mol, plan.attachment_mapping)))
}

/// `molzip(a, b, params)`: joins `a` and `b` at their labelled atoms.
pub fn molzip(a: &ROMol, b: &ROMol, params: &MolzipParams) -> Result<Option<RWMol>> {
    Ok(molzip_with_mapping(a, b, params)?.map(|(mol, _)| mol))
}

/// `molzip(a, params)`: joins the labelled atoms within one molecule.
pub fn molzip_single(a: &ROMol, params: &MolzipParams) -> Result<Option<RWMol>> {
    molzip(a, &ROMol::new(), params)
}

/// `molzip(decomposition, params)`: zips the pieces of an R group
/// decomposition (the core first). With `generate_coordinates` the zipped
/// molecule gets 2D coordinates, which are also assigned to the pieces.
pub fn molzip_fragments(
    decomposition: &mut [RWMol],
    params: &MolzipParams,
) -> Result<Option<RWMol>> {
    if params.generate_coordinates {
        let mut index = 0;
        for mol in decomposition.iter_mut() {
            for atom in 0..mol.num_atoms() {
                mol.atom_mut(atom).set_prop(INDEX_PROP, index, false)?;
                index += 1;
            }
        }
    }
    if decomposition.is_empty() {
        return Ok(None);
    }

    // when the rgroup decomposition splits a ring, it puts it in both
    // rgroups, so remove these
    let mut unique: Vec<usize> = Vec::new();
    if params.label != MolzipLabel::FragmentOnBonds && decomposition.len() > 1 {
        let mut existing = Vec::new();
        for idx in 1..decomposition.len() {
            let smiles = mol_to_smiles(&mut decomposition[idx], &SmilesWriteParams::default())?;
            if !existing.contains(&smiles) {
                existing.push(smiles);
                unique.push(idx);
            }
        }
    }
    let mut combined = decomposition[0].clone();
    for &idx in &unique {
        combined = combine_mols(&combined, &decomposition[idx], Point3D::default())?;
    }

    let Some((mut zipped, attachment_mapping)) =
        molzip_with_mapping(&combined, &ROMol::new(), params)?
    else {
        return Ok(None);
    };

    if params.generate_coordinates && zipped.num_atoms() > 0 {
        let conf_id = compute_2d_coords(&mut zipped, &Compute2DCoordParameters::default())?;
        let zipped_conf = zipped.conformer(conf_id as i32)?.clone();
        for mol in decomposition.iter_mut() {
            let mut new_conf = Conformer::new(mol.num_atoms() as usize);
            new_conf.set_3d(false);
            for atom in 0..mol.num_atoms() {
                let mut zipped_index = mol
                    .atom(atom)
                    .get_prop(INDEX_PROP)
                    .and_then(RDValue::as_int)
                    .ok_or_else(|| RDKitError::Key(INDEX_PROP.into()))?;
                mol.atom_mut(atom).clear_prop(INDEX_PROP);
                if let Some(&replacement) = attachment_mapping.get(&zipped_index) {
                    zipped_index = replacement;
                }
                let zipped_atom = (0..zipped.num_atoms())
                    .find(|&z| index_prop(&zipped, z) == Some(zipped_index))
                    .ok_or_else(|| {
                        invariant(
                            735,
                            "zippedAtom != zipppedAtoms.end()",
                            "molzip: cannot map atom into zipped molecule",
                        )
                    })?;
                new_conf.set_atom_pos(atom, zipped_conf.atom_pos(zipped_atom)?)?;
            }
            mol.add_conformer(new_conf, true)?;
        }
        for atom in 0..zipped.num_atoms() {
            zipped.atom_mut(atom).clear_prop(INDEX_PROP);
        }
    }
    Ok(Some(zipped))
}

/// `molzip(row, params)`: zips an R group row (`Core` plus the R groups, in
/// key order) back into a molecule.
pub fn molzip_row(
    row: &mut BTreeMap<String, RWMol>,
    params: &MolzipParams,
) -> Result<Option<RWMol>> {
    let core = row.remove("Core").ok_or_else(|| {
        precondition(
            753,
            "core != row.end()",
            "RGroup has no Core, cannot molzip",
        )
    })?;
    let mut mols = vec![core];
    mols.extend(std::mem::take(row).into_values());
    molzip_fragments(&mut mols, params)
}
