//! Molecule validation: the MolVS checks (no atoms, fragments, charge,
//! isotopes), allowed/disallowed atoms, radicals, features, 2D layout and
//! wedged-bond stereo checks.
//!
//! Ported from RDKit `Code/GraphMol/MolStandardize/Validate.{h,cpp}`.

use std::f64::consts::PI;

use super::catalogs::FragmentCatalogParams;
use crate::error::{RDKitError, Result};
use crate::graphmol::atom::Atom;
use crate::graphmol::bond::{BondDir, BondType};
use crate::graphmol::canon::can_have_direction;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::find_rings::fast_find_rings;
use crate::graphmol::mol_ops::frags::get_mol_frags_as_mols;
use crate::graphmol::mol_ops::structure::get_formal_charge;
use crate::graphmol::periodic_table::PeriodicTable;
use crate::graphmol::queries::is_atom_dummy;
use crate::graphmol::romol::ROMol;
use crate::graphmol::smiles_parse::{SmilesParserParams, mol_from_smiles};
use crate::graphmol::substruct::{SubstructMatchParameters, substruct_match};
use crate::rdgeneral::common_properties;
use crate::rdgeneral::cxx_sort::libcxx_sort_by;
use crate::rdgeneral::dict::RDProps;

/// `FeaturesValidation`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct FeaturesValidation {
    pub allow_enhanced_stereo: bool,
    pub allow_aromatic_bond_type: bool,
    pub allow_dative_bond_type: bool,
    pub allow_queries: bool,
    pub allow_dummies: bool,
    pub allow_atom_aliases: bool,
}

/// `Layout2DValidation`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Layout2DValidation {
    pub clash_limit: f64,
    pub bond_length_limit: f64,
    pub allow_long_bonds_in_rings: bool,
    pub allow_atom_bond_clash_exemption: bool,
    pub min_median_bond_length: f64,
}

impl Default for Layout2DValidation {
    fn default() -> Self {
        Layout2DValidation {
            clash_limit: 0.15,
            bond_length_limit: 25.0,
            allow_long_bonds_in_rings: true,
            allow_atom_bond_clash_exemption: true,
            min_median_bond_length: 1e-3,
        }
    }
}

/// A `ValidationMethod`.
#[derive(Debug, Clone)]
pub enum ValidationMethod {
    /// `RDKitValidation`: valence errors
    RDKit { allow_empty_molecules: bool },
    /// `NoAtomValidation`
    NoAtom,
    /// `FragmentValidation`: fragments from the default fragment list
    Fragment,
    /// `NeutralValidation`
    Neutral,
    /// `IsotopeValidation`
    Isotope { strict: bool },
    /// `CompositeValidation` / `MolVSValidation`
    Composite(Vec<ValidationMethod>),
    /// `AllowedAtomsValidation`
    AllowedAtoms(Vec<Atom>),
    /// `DisallowedAtomsValidation`
    DisallowedAtoms(Vec<Atom>),
    /// `DisallowedRadicalValidation`
    DisallowedRadical,
    /// `FeaturesValidation`
    Features(FeaturesValidation),
    /// `Is2DValidation`
    Is2D { threshold: f64 },
    /// `Layout2DValidation`
    Layout2D(Layout2DValidation),
    /// `StereoValidation`
    Stereo,
}

/// `MolVSValidation()`: no atoms, fragments, neutral, isotopes.
pub fn mol_vs_validation() -> ValidationMethod {
    ValidationMethod::Composite(vec![
        ValidationMethod::NoAtom,
        ValidationMethod::Fragment,
        ValidationMethod::Neutral,
        ValidationMethod::Isotope { strict: false },
    ])
}

/// `Atom::Match` / `QueryAtom::Match` of a template atom against `mol`'s atom.
fn template_matches(tmpl: &Atom, mol: &ROMol, idx: u32) -> Result<bool> {
    match tmpl.query() {
        Some(q) => q.matches_with(mol, idx, None),
        None => Ok(tmpl.matches(mol.atom(idx))),
    }
}

impl ValidationMethod {
    /// `validate(mol, reportAllFailures)`: the error/info messages.
    pub fn validate(&self, mol: &ROMol, report_all_failures: bool) -> Result<Vec<String>> {
        match self {
            ValidationMethod::RDKit {
                allow_empty_molecules,
            } => rdkit_validation(mol, *allow_empty_molecules, report_all_failures),
            ValidationMethod::NoAtom => Ok(if mol.num_atoms() == 0 {
                vec!["ERROR: [NoAtomValidation] Molecule has no atoms".into()]
            } else {
                Vec::new()
            }),
            ValidationMethod::Fragment => fragment_validation(mol),
            ValidationMethod::Neutral => {
                let charge = get_formal_charge(mol);
                Ok(if charge != 0 {
                    let charge_str = if charge > 0 {
                        format!("+{charge}")
                    } else {
                        charge.to_string()
                    };
                    vec![format!(
                        "INFO: [NeutralValidation] Not an overall neutral system ({charge_str})"
                    )]
                } else {
                    Vec::new()
                })
            }
            ValidationMethod::Isotope { strict } => {
                isotope_validation(mol, *strict, report_all_failures)
            }
            ValidationMethod::Composite(validations) => {
                let mut errors = Vec::new();
                for v in validations {
                    let partial = v.validate(mol, report_all_failures)?;
                    if !partial.is_empty() {
                        errors.extend(partial);
                        if !report_all_failures {
                            break;
                        }
                    }
                }
                Ok(errors)
            }
            ValidationMethod::AllowedAtoms(allowed) => {
                let mut errors = Vec::new();
                for i in 0..mol.num_atoms() {
                    if !report_all_failures && !errors.is_empty() {
                        break;
                    }
                    let mut matched = false;
                    for a in allowed {
                        if template_matches(a, mol, i)? {
                            matched = true;
                            break;
                        }
                    }
                    if !matched {
                        errors.push(format!(
                            "INFO: [AllowedAtomsValidation] Atom {} is not in allowedAtoms list",
                            mol.atom(i).symbol()
                        ));
                    }
                }
                Ok(errors)
            }
            ValidationMethod::DisallowedAtoms(disallowed) => {
                let mut errors = Vec::new();
                for i in 0..mol.num_atoms() {
                    if !report_all_failures && !errors.is_empty() {
                        break;
                    }
                    let mut matched = false;
                    for a in disallowed {
                        if template_matches(a, mol, i)? {
                            matched = true;
                        }
                    }
                    if matched {
                        errors.push(format!(
                            "INFO: [DisallowedAtomsValidation] Atom {} is in disallowedAtoms list",
                            mol.atom(i).symbol()
                        ));
                    }
                }
                Ok(errors)
            }
            ValidationMethod::DisallowedRadical => {
                disallowed_radical_validation(mol, report_all_failures)
            }
            ValidationMethod::Features(f) => features_validation(mol, f, report_all_failures),
            ValidationMethod::Is2D { threshold } => {
                is_2d_validation(mol, *threshold, report_all_failures)
            }
            ValidationMethod::Layout2D(l) => layout_2d_validation(mol, l, report_all_failures),
            ValidationMethod::Stereo => {
                let mut errors = Vec::new();
                for i in 0..mol.num_atoms() {
                    if has_stereo_bond(mol, i) {
                        check_stereo(mol, i, report_all_failures, &mut errors)?;
                    }
                    if !errors.is_empty() && !report_all_failures {
                        break;
                    }
                }
                Ok(errors)
            }
        }
    }
}

fn rdkit_validation(
    mol: &ROMol,
    allow_empty_molecules: bool,
    report_all_failures: bool,
) -> Result<Vec<String>> {
    let mut mol_copy = mol.clone();
    let mut errors = Vec::new();
    let na = mol.num_atoms();
    if na == 0 && !allow_empty_molecules {
        errors.push("ERROR: [NoAtomValidation] Molecule has no atoms".to_string());
    }
    for i in 0..na {
        if !report_all_failures && !errors.is_empty() {
            break;
        }
        match mol_copy.calc_explicit_valence(i, true) {
            Ok(_) => {}
            Err(e @ RDKitError::Sanitize(_)) => {
                errors.push(format!("INFO: [ValenceValidation] {e}"));
            }
            Err(e) => return Err(e),
        }
    }
    Ok(errors)
}

fn fragment_validation(mol: &ROMol) -> Result<Vec<String>> {
    let mut errors = Vec::new();
    let fparams = FragmentCatalogParams::defaults();
    let (frags, mut atom_mapping) = get_mol_frags_as_mols(mol, true, true)?;
    for fgrp in fparams.func_groups() {
        let fname = fgrp
            .get_prop(common_properties::_Name)
            .map(|v| v.to_rd_string())
            .unwrap_or_default();
        let res = substruct_match(mol, fgrp, &SubstructMatchParameters::default())?;
        if res.is_empty() || frags.is_empty() {
            continue;
        }
        let mut substructmap: Vec<Vec<i32>> = res
            .iter()
            .map(|m| m.iter().map(|&(_, a)| a).collect())
            .collect();
        // report each fragment once, even if present several times
        let mut fpresent = false;
        for molfragidx in &mut atom_mapping {
            molfragidx.sort_unstable();
            for substructidx in &mut substructmap {
                substructidx.sort_unstable();
                if molfragidx == substructidx && !fpresent {
                    errors.push(format!("INFO: [FragmentValidation] {fname} is present"));
                    fpresent = true;
                }
            }
        }
    }
    Ok(errors)
}

fn isotope_validation(mol: &ROMol, strict: bool, report_all_failures: bool) -> Result<Vec<String>> {
    let mut errors = Vec::new();
    for atom in mol.atoms() {
        let isotope = atom.isotope();
        if isotope == 0 {
            continue;
        }
        let symbol = atom.symbol();
        let atomic_num = atom.atomic_num();
        if atomic_num != 0 && strict {
            let mass = PeriodicTable::get().get_mass_for_isotope(atomic_num, isotope)?;
            if mass == 0.0 {
                errors.push(format!(
                    "ERROR: [IsotopeValidation] The molecule contains an unknown isotope: {isotope}{symbol}"
                ));
            }
        } else {
            errors.push(format!(
                "INFO: [IsotopeValidation] Molecule contains isotope {isotope}{symbol}"
            ));
        }
        if !errors.is_empty() && !report_all_failures {
            break;
        }
    }
    Ok(errors)
}

fn disallowed_radical_validation(mol: &ROMol, report_all_failures: bool) -> Result<Vec<String>> {
    let mut errors = Vec::new();
    for (idx, atom) in mol.atoms().iter().enumerate() {
        let n_rad = atom.num_radical_electrons();
        if n_rad == 0 {
            continue;
        }
        let atomic_num = atom.atomic_num();
        let degree = mol.atom_degree(idx as u32);
        if (atomic_num == 7 || atomic_num == 8) && n_rad == 1 && degree == 1 {
            let (mut nbr_num, mut bond_type) = (0, BondType::Unspecified);
            for b in mol.atom_bonds(idx as u32) {
                // a single iteration, the degree is 1
                nbr_num = mol.atom(b.other_atom_idx(idx as u32)?).atomic_num();
                bond_type = b.bond_type();
            }
            // nitric oxide
            if atomic_num == 7 && nbr_num == 8 && bond_type == BondType::Double {
                continue;
            }
            // aminoxyl
            if atomic_num == 8 && nbr_num == 7 && bond_type == BondType::Single {
                continue;
            }
        }
        errors.push(format!(
            "ERROR: [DisallowedRadicalValidation] The radical at atom {idx} is not allowed"
        ));
        if !report_all_failures {
            break;
        }
    }
    Ok(errors)
}

fn features_validation(
    mol: &ROMol,
    f: &FeaturesValidation,
    report_all_failures: bool,
) -> Result<Vec<String>> {
    let mut errors = Vec::new();
    // query and dummy atoms, aliases
    for (idx, atom) in mol.atoms().iter().enumerate() {
        if !f.allow_queries && atom.has_query() {
            errors.push(format!(
                "ERROR: [FeaturesValidation] Query atom {idx} is not allowed"
            ));
            if !report_all_failures {
                return Ok(errors);
            }
        } else if !f.allow_dummies && is_atom_dummy(mol, idx as u32) {
            errors.push(format!(
                "ERROR: [FeaturesValidation] Dummy atom {idx} is not allowed"
            ));
            if !report_all_failures {
                return Ok(errors);
            }
        }
        if !f.allow_atom_aliases {
            if let Some(alias) = atom.get_prop(common_properties::molFileAlias) {
                errors.push(format!(
                    "ERROR: [FeaturesValidation] Atom {idx} with alias '{}' is not allowed",
                    alias.to_rd_string()
                ));
                if !report_all_failures {
                    return Ok(errors);
                }
            }
        }
    }
    // query, aromatic or dative bonds
    for (idx, bond) in mol.bonds().iter().enumerate() {
        if !f.allow_queries && bond.has_query() {
            errors.push(format!(
                "ERROR: [FeaturesValidation] Query bond {idx} is not allowed"
            ));
            if !report_all_failures {
                return Ok(errors);
            }
        }
        if !f.allow_aromatic_bond_type && bond.bond_type() == BondType::Aromatic {
            errors.push(format!(
                "ERROR: [FeaturesValidation] Bond {idx} of aromatic type is not allowed"
            ));
            if !report_all_failures {
                return Ok(errors);
            }
        }
        if !f.allow_dative_bond_type && bond.bond_type() == BondType::Dative {
            errors.push(format!(
                "ERROR: [FeaturesValidation] Bond {idx} of dative type is not allowed"
            ));
            if !report_all_failures {
                return Ok(errors);
            }
        }
    }
    // enhanced stereochemistry
    if !f.allow_enhanced_stereo && !mol.stereo_groups().is_empty() {
        errors.push(
            "ERROR: [FeaturesValidation] Enhanced stereochemistry features are not allowed".into(),
        );
    }
    Ok(errors)
}

fn is_2d_validation(mol: &ROMol, threshold: f64, report_all_failures: bool) -> Result<Vec<String>> {
    let mut errors = Vec::new();
    if mol.num_conformers() == 0 {
        errors.push("ERROR: [Is2DValidation] The molecule has no coordinates".into());
        return Ok(errors);
    }
    let conf = mol.conformer(-1)?;
    if conf.is_3d() {
        errors.push("ERROR: [Is2DValidation] The molecule includes non-null Z coordinates".into());
        return Ok(errors);
    }
    // conf.is_3d() comes from the parser; also check molecules built otherwise
    let max_absz = conf
        .positions()
        .iter()
        .fold(0.0f64, |m, p| m.max(p.z.abs()));
    if max_absz > threshold {
        errors.push("ERROR: [Is2DValidation] The molecule includes non-null Z coordinates".into());
        if !report_all_failures {
            return Ok(errors);
        }
    }
    if conf.num_atoms() < 2 {
        return Ok(errors);
    }
    // are all the atoms at the same position (often: no coordinates
    // assigned, all atoms at the origin)? (upstream starts the maxima at
    // `numeric_limits<double>::min()`, the smallest positive value)
    let (mut min_x, mut max_x) = (f64::MAX, f64::MIN_POSITIVE);
    let (mut min_y, mut max_y) = (f64::MAX, f64::MIN_POSITIVE);
    for p in conf.positions() {
        min_x = p.x.min(min_x);
        max_x = p.x.max(max_x);
        min_y = p.y.min(min_y);
        max_y = p.y.max(max_y);
    }
    let max_delta = (max_x - min_x).max(max_y - min_y);
    if max_delta < threshold {
        errors.push("ERROR: [Is2DValidation] All atoms have the same (x,y) coordinates".into());
        if !report_all_failures {
            return Ok(errors);
        }
    }
    Ok(errors)
}

/// `Layout2DValidation::squaredMedianBondLength`: the squared median bond
/// length, bonds of null length excluded.
pub fn squared_median_bond_length(mol: &ROMol, conf: &Conformer) -> Result<f64> {
    let mut values = Vec::with_capacity(mol.num_bonds() as usize);
    for b in mol.bonds() {
        let v = (conf.atom_pos(b.begin_atom_idx())? - conf.atom_pos(b.end_atom_idx())?).length_sq();
        if v > 0.0 {
            values.push(v);
        }
    }
    if values.is_empty() {
        return Ok(0.0);
    }
    values.sort_by(f64::total_cmp);
    let n = values.len();
    Ok(if n % 2 == 1 {
        values[n / 2]
    } else {
        0.5 * (values[n / 2 - 1] + values[n / 2])
    })
}

fn layout_2d_validation(
    mol: &ROMol,
    l: &Layout2DValidation,
    report_all_failures: bool,
) -> Result<Vec<String>> {
    let mut errors = Vec::new();
    if mol.num_conformers() == 0 {
        errors.push("ERROR: [Layout2DValidation] The molecule has no coordinates".into());
        return Ok(errors);
    }
    let conf = mol.conformer(-1)?;
    let natoms = conf.num_atoms();
    if natoms < 2 {
        return Ok(errors);
    }
    // thresholds relative to the median squared bond length
    let reference = squared_median_bond_length(mol, conf)?;
    if reference < l.min_median_bond_length * l.min_median_bond_length {
        errors.push(
            "ERROR: [Layout2DValidation] The median bond length is smaller than the configured limit"
                .into(),
        );
        if !report_all_failures {
            return Ok(errors);
        }
    }
    // atoms clashing with other atoms
    let atom_clash_threshold = l.clash_limit * l.clash_limit * reference;
    for i in 0..natoms - 1 {
        let pi = conf.atom_pos(i)?;
        for j in i + 1..natoms {
            let d2 = (pi - conf.atom_pos(j)?).length_sq();
            if d2 < atom_clash_threshold {
                errors.push(format!(
                    "ERROR: [Layout2DValidation] Atom {i} is too close to atom {j}"
                ));
                if !report_all_failures {
                    return Ok(errors);
                }
            }
        }
    }
    // ring info is needed for the exemptions
    let mut with_rings;
    let mol = if (l.allow_long_bonds_in_rings || l.allow_atom_bond_clash_exemption)
        && !mol.ring_info().is_initialized()
    {
        with_rings = mol.clone();
        fast_find_rings(&mut with_rings)?;
        &with_rings
    } else {
        mol
    };
    let bond_length_threshold = l.bond_length_limit * l.bond_length_limit * reference;
    for bond in mol.bonds() {
        let (i, j) = (bond.begin_atom_idx(), bond.end_atom_idx());
        let (pi, pj) = (conf.atom_pos(i)?, conf.atom_pos(j)?);
        let ll = (pi - pj).length_sq();
        let bond_rings = || mol.ring_info().num_bond_rings(bond.idx());
        // exceedingly long bonds
        if (!l.allow_long_bonds_in_rings || bond_rings()? == 0) && ll > bond_length_threshold {
            errors.push(format!(
                "ERROR: [Layout2DValidation] The length of bond {} between atoms {i} and {j} exceeds a configured limit",
                bond.idx()
            ));
            if !report_all_failures {
                return Ok(errors);
            }
        }
        // is this bond exempted from atom-bond collision detection?
        if l.allow_atom_bond_clash_exemption && ll > 5.0 * 5.0 * reference && bond_rings()? != 0 {
            continue;
        }
        // atoms clashing with this bond
        for k in 0..natoms {
            if k == i || k == j {
                continue;
            }
            let pk = conf.atom_pos(k)?;
            let vik = pk - pi;
            let vij = pj - pi;
            let rr = vik.length_sq();
            let bb = vij.length_sq();
            let rb = vik.dot_product(&vij);
            // prevents dividing by zero in extreme cases
            const EPS: f64 = 1.0e-7;
            let kb = (rr * bb - rb * rb) / (bb + EPS);
            if rb >= 0.0 && rb <= bb && kb < atom_clash_threshold {
                errors.push(format!(
                    "ERROR: [Layout2DValidation] Atom {k} too close to bond {}",
                    bond.idx()
                ));
                if !report_all_failures {
                    return Ok(errors);
                }
            }
        }
    }
    Ok(errors)
}

// ------------------------------------------------------------------- stereo

fn has_stereo_bond(mol: &ROMol, atom: u32) -> bool {
    mol.atom_bonds(atom).any(|b| {
        b.begin_atom_idx() == atom
            && matches!(
                b.bond_dir(),
                BondDir::BeginDash | BondDir::BeginWedge | BondDir::Unknown
            )
    })
}

#[derive(Debug, Clone)]
struct BondInfo {
    bond: u32,
    other: u32,
    bond_dir: BondDir,
    angle: f64,
}

#[derive(Debug, Default)]
struct BondDirCount {
    wedge: u32,
    dash: u32,
    unknown: u32,
    other: u32,
}

struct NeighborsInfo {
    bonds: Vec<BondInfo>,
    dir_count: BondDirCount,
}

impl NeighborsInfo {
    fn new(mol: &ROMol, atom: u32) -> Result<Self> {
        let mut bonds = Vec::new();
        for b in mol.atom_bonds(atom) {
            // the directions of bonds that begin at the neighbors don't count
            let bond_dir = if b.begin_atom_idx() == atom {
                b.bond_dir()
            } else {
                BondDir::None
            };
            bonds.push(BondInfo {
                bond: b.idx(),
                other: b.other_atom_idx(atom)?,
                bond_dir,
                angle: 0.0,
            });
        }
        let mut dir_count = BondDirCount::default();
        for info in &bonds {
            match info.bond_dir {
                BondDir::BeginDash => dir_count.dash += 1,
                BondDir::BeginWedge => dir_count.wedge += 1,
                BondDir::Unknown => dir_count.unknown += 1,
                // unspecified, or describing an adjacent double bond
                BondDir::None | BondDir::EndUpRight | BondDir::EndDownRight => {}
                _ => dir_count.other += 1,
            }
        }
        // sort the neighbors by the angle they form with the first one
        let conf = mol.conformer(-1)?;
        let p = conf.atom_pos(atom)?;
        let v0 = conf.atom_pos(bonds[0].other)? - p;
        for info in bonds.iter_mut().skip(1) {
            let vn = conf.atom_pos(info.other)? - p;
            info.angle = v0.signed_angle_to(&vn);
        }
        libcxx_sort_by(&mut bonds[1..], |a, b| a.angle < b.angle);
        Ok(NeighborsInfo { bonds, dir_count })
    }
}

fn is_wedge_or_dash(dir: BondDir) -> bool {
    matches!(dir, BondDir::BeginDash | BondDir::BeginWedge)
}

fn check_3_coordinated_stereo(
    mol: &ROMol,
    atom: u32,
    info: &NeighborsInfo,
    errors: &mut Vec<String>,
) -> Result<()> {
    let num_stereo_bonds = info.dir_count.dash + info.dir_count.wedge;
    if num_stereo_bonds == 1 {
        // the stereo bond
        let i = (0..3)
            .find(|&i| is_wedge_or_dash(info.bonds[i].bond_dir))
            .unwrap_or(3);
        // are the stereocenter and the other two ligands colinear?
        let conf = mol.conformer(-1)?;
        let p = conf.atom_pos(atom)?;
        let va = conf.atom_pos(info.bonds[(i + 1) % 3].other)? - p;
        let vb = conf.atom_pos(info.bonds[(i + 2) % 3].other)? - p;
        let angle = va.angle_to(&vb);
        const ANGLE_EPSILON: f64 = PI * 5.0 / 180.0; // 5 degrees
        if angle < ANGLE_EPSILON || (PI - angle) < ANGLE_EPSILON {
            errors.push(format!(
                "ERROR: [StereoValidation] Colinearity of non-stereo bonds at atom {atom}"
            ));
        }
    } else {
        // multiple stereo bonds on a center with 3 explicit ligands are not
        // allowed (as in the AvalonTools' struchk)
        errors.push(format!(
            "ERROR: [StereoValidation] Atom {atom} has 3 explicit substituents and multiple stereo bonds"
        ));
    }
    Ok(())
}

fn check_4_coordinated_stereo(
    mol: &ROMol,
    atom: u32,
    info: &NeighborsInfo,
    report_all_failures: bool,
    errors: &mut Vec<String>,
) -> Result<()> {
    let dirs: Vec<BondDir> = info.bonds.iter().map(|b| b.bond_dir).collect();
    if info.dir_count.dash > 2 || info.dir_count.wedge > 2 {
        errors.push(format!(
            "ERROR: [StereoValidation] Atom {atom} has too many stereo bonds with like orientation"
        ));
        if !report_all_failures {
            return Ok(());
        }
    }
    for i in 0..2 {
        if (dirs[i] == BondDir::BeginDash && dirs[i + 2] == BondDir::BeginWedge)
            || (dirs[i] == BondDir::BeginWedge && dirs[i + 2] == BondDir::BeginDash)
        {
            errors.push(format!(
                "ERROR: [StereoValidation] Atom {atom} has opposing stereo bonds with different up/down orientation"
            ));
            if !report_all_failures {
                return Ok(());
            }
        }
    }
    for i in 0..4 {
        let next = dirs[(i + 1) % 4];
        if (dirs[i] == BondDir::BeginDash && next == BondDir::BeginDash)
            || (dirs[i] == BondDir::BeginWedge && next == BondDir::BeginWedge)
        {
            errors.push(format!(
                "ERROR: [StereoValidation] Atom {atom} has adjacent stereo bonds with like orientation"
            ));
            if !report_all_failures {
                return Ok(());
            }
            // reported once per atom
            break;
        }
    }
    if info.dir_count.dash + info.dir_count.wedge == 1 {
        // a single wedged/dashed bond: check for 'umbrellas' and other
        // geometric violations
        let conf = mol.conformer(-1)?;
        let p = conf.atom_pos(atom)?;
        if let Some(i) = (0..4).find(|&i| is_wedge_or_dash(dirs[i])) {
            // how many of the other bonds lie on the opposite half-plane
            let vi = conf.atom_pos(info.bonds[i].other)? - p;
            let mut opposed = 0;
            for j in (0..4).filter(|&j| j != i) {
                let vj = conf.atom_pos(info.bonds[j].other)? - p;
                if vi.angle_to(&vj) > 95.0 * PI / 180.0 {
                    opposed += 1;
                }
            }
            if opposed == 3 {
                errors.push(format!(
                    "ERROR: [StereoValidation] Atom {atom} has a potentially ambiguous representation: all non-stereo bonds opposite to the only stereo bond"
                ));
            }
            if !report_all_failures {
                return Ok(());
            }
        }
        // colinearity, or a badly placed (too short) middle non-stereo bond
        if let Some(i) = (0..4).find(|&i| is_wedge_or_dash(dirs[i])) {
            let pj = conf.atom_pos(info.bonds[(i + 1) % 4].other)?;
            let pk = conf.atom_pos(info.bonds[(i + 2) % 4].other)?;
            let pl = conf.atom_pos(info.bonds[(i + 3) % 4].other)?;
            let angle = (pj - pk).signed_angle_to(&(pl - pk));
            if angle < 185.0 * PI / 180.0 {
                errors.push(format!(
                    "ERROR: [StereoValidation] Colinearity or triangle rule violation of non-stereo bonds at atom {atom}"
                ));
                if !report_all_failures {
                    return Ok(());
                }
            }
        }
    }
    Ok(())
}

fn check_stereo(
    mol: &ROMol,
    atom: u32,
    report_all_failures: bool,
    errors: &mut Vec<String>,
) -> Result<()> {
    let info = NeighborsInfo::new(mol, atom)?;
    if info.dir_count.other > 0 {
        // unlikely; no point in looking further
        errors.push(format!(
            "ERROR: [StereoValidation] one or more bonds incident to atom {atom} have unexpected direction settings"
        ));
        return Ok(());
    }
    if info.dir_count.unknown > 0 {
        if info.dir_count.dash > 0 || info.dir_count.wedge > 0 {
            errors.push(format!(
                "ERROR: [StereoValidation] Atom {atom} has both unknown and wedged/dashed stereo bonds."
            ));
        }
        return Ok(());
    }
    for b in &info.bonds {
        let is_stereo = matches!(
            b.bond_dir,
            BondDir::BeginDash | BondDir::BeginWedge | BondDir::Unknown
        );
        if is_stereo && !can_have_direction(mol.bond(b.bond).bond_type()) {
            errors.push(format!(
                "ERROR: [StereoValidation] Bond {} has assigned stereo type, but unexpected bond order.",
                b.bond
            ));
            if !report_all_failures {
                return Ok(());
            }
        }
    }
    // only some categories of stereocenters are validated
    let (mut multiple_bond_found, mut possible_allene) = (false, false);
    for b in mol.atom_bonds(atom) {
        if b.bond_type() != BondType::Single {
            multiple_bond_found = true;
            let other = b.other_atom_idx(atom)?;
            if mol.atom_degree(other) == 2 {
                let double_bonds = mol
                    .atom_bonds(other)
                    .filter(|ob| ob.bond_type() == BondType::Double)
                    .count();
                if double_bonds == 2 {
                    possible_allene = true;
                }
            }
        }
    }
    let atomic_num = mol.atom(atom).atomic_num();
    // allenes and P compounds are not validated
    if possible_allene || (multiple_bond_found && atomic_num == 15) {
        return Ok(());
    }
    // stereo bonds at unsaturated atoms have valid uses (wavy bonds,
    // atropisomers) and are not validated
    if multiple_bond_found && atomic_num != 16 {
        return Ok(());
    }
    match mol.atom_degree(atom) {
        1 | 2 => errors.push(format!(
            "ERROR: [StereoValidation] Atom {atom} has stereo bonds, but less than 3 explicit substituents."
        )),
        3 => check_3_coordinated_stereo(mol, atom, &info, errors)?,
        4 => check_4_coordinated_stereo(mol, atom, &info, report_all_failures, errors)?,
        _ => {}
    }
    Ok(())
}

/// `validateSmiles`: the MolVS validation of a SMILES.
pub fn validate_smiles(smiles: &str) -> Result<Vec<String>> {
    let mol = mol_from_smiles(smiles, &SmilesParserParams::default())?.ok_or_else(|| {
        RDKitError::Value(format!(
            "SMILES Parse Error: syntax error for input: {smiles}"
        ))
    })?;
    mol_vs_validation().validate(&mol, true)
}
