//! Molecule hashes (upstream `GraphMol/MolHash`, the NextMove "nmmolhash"
//! functions): canonical SMILES of progressively abstracted molecules.
//!
//! With `use_cx_smiles` the hashes are CXSMILES: the extension block follows
//! the SMILES (except for its coordinates and the fields in
//! `cx_flags_to_skip`).

mod formula;
mod tautomer_v2;

use crate::error::Result;
use crate::graphmol::bond::{Bond, BondStereo, BondType};
use crate::graphmol::chirality::legacy::assign_stereochemistry;
use crate::graphmol::find_rings::fast_find_rings;
use crate::graphmol::mol_ops::add_hs::{RemoveHsParameters, remove_hs};
use crate::graphmol::mol_ops::assign_radicals;
use crate::graphmol::mol_ops::frags::{get_mol_frags_as_mols, get_mol_frags_mapping};
use crate::graphmol::romol::{ROMol, RWMol};
use crate::graphmol::smiles_parse::cx_write::{cx_fields, get_cx_extensions};
use crate::graphmol::smiles_parse::write::{SmilesWriteParams, mol_to_smiles_detail};
use crate::graphmol::{Atom, ChiralType};

/// Atoms or bonds with this property set (any value) are left out of the
/// tautomeric system of the `HetAtomTautomerv2` / `HetAtomProtomerv2` hashes.
pub const EXCLUDE_FROM_TAUTOMERISM_PROP: &str = "_MolHashExcludeFromTautomerism";

/// `MolHash::HashFunction`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HashFunction {
    AnonymousGraph = 1,
    ElementGraph = 2,
    CanonicalSmiles = 3,
    MurckoScaffold = 4,
    ExtendedMurcko = 5,
    MolFormula = 6,
    AtomBondCounts = 7,
    DegreeVector = 8,
    Mesomer = 9,
    HetAtomTautomer = 10,
    HetAtomProtomer = 11,
    RedoxPair = 12,
    Regioisomer = 13,
    NetCharge = 14,
    SmallWorldIndexBR = 15,
    SmallWorldIndexBRL = 16,
    ArthorSubstructureOrder = 17,
    HetAtomTautomerv2 = 18,
    HetAtomProtomerv2 = 19,
}

/// `MolHash::StripType` bit flags for [`strip`].
pub mod strip_type {
    pub const ATOM_STEREO: u32 = 1;
    pub const BOND_STEREO: u32 = 2;
    pub const ISOTOPE: u32 = 4;
    pub const ATOM_MAP: u32 = 8;
    pub const HYDROGEN: u32 = 16;
}

/// Whether the hash is a CXSMILES and which of its fields to leave out
/// (`useCXSmiles`, `cxFlagsToSkip`).
#[derive(Debug, Clone, Copy)]
pub(crate) struct Cx {
    enabled: bool,
    skip: u32,
}

/// Writes the canonical SMILES of `mol` (`convertToSmilesWithCXFlags`);
/// `all_explicit` is what the v2 tautomer hash needs.
fn smiles(mol: &mut RWMol, all_explicit: bool, cx: Cx) -> Result<String> {
    let params = SmilesWriteParams {
        all_bonds_explicit: all_explicit,
        all_hs_explicit: all_explicit,
        ..Default::default()
    };
    let skip_stereo_groups = cx.skip & cx_fields::ENHANCEDSTEREO != 0;
    mol_to_smiles_detail(mol, &params, cx.enabled, !skip_stereo_groups)
}

/// `addCXExtensions`: appends the extension block, without coordinates and
/// the skipped fields (`additional_skips` come on top).
fn add_cx_extensions(
    mol: &mut RWMol,
    result: &mut String,
    cx: Cx,
    additional_skips: u32,
) -> Result<()> {
    if !cx.enabled {
        return Ok(());
    }
    let skip = cx.skip | additional_skips | cx_fields::COORDS;
    let ext = get_cx_extensions(mol, cx_fields::ALL ^ skip)?;
    if !ext.is_empty() {
        *result += " ";
        *result += &ext;
    }
    Ok(())
}

/// Radicals and stereochemistry after the structure was abstracted: stereo
/// elements may just have been destroyed.
fn clean_up(mol: &mut RWMol) -> Result<()> {
    assign_radicals(mol);
    assign_stereochemistry(mol, true, true, false)
}

/// `NMRDKitBondGetOrder`
fn bond_order(bond: &Bond) -> u32 {
    match bond.bond_type() {
        BondType::Aromatic | BondType::Single => 1,
        BondType::Double => 2,
        BondType::Triple => 3,
        BondType::Quadruple => 4,
        BondType::Quintuple => 5,
        BondType::Hextuple => 6,
        _ => 0,
    }
}

/// `NMRDKitSanitizeHydrogens`: all implicit Hs become explicit.
fn sanitize_hydrogens(mol: &mut RWMol) -> Result<()> {
    for idx in 0..mol.num_atoms() {
        let hs = mol.total_num_hs(idx, false)?;
        let atom = mol.atom_mut(idx);
        atom.set_no_implicit(true);
        atom.set_num_explicit_hs(hs);
        mol.update_atom_property_cache(idx, false)?;
    }
    Ok(())
}

fn in_ring(mol: &ROMol, idx: u32) -> bool {
    mol.ring_info().num_atom_rings(idx).is_ok_and(|n| n > 0)
}

fn bond_in_ring(mol: &ROMol, idx: u32) -> bool {
    mol.ring_info().num_bond_rings(idx).is_ok_and(|n| n > 0)
}

fn ensure_rings(mol: &mut ROMol) -> Result<()> {
    if !mol.ring_info().is_find_fast_or_better() {
        fast_find_rings(mol)?;
    }
    Ok(())
}

/// `NormalizeHCount`: the H count that fills the usual valence.
fn normalize_h_count(mol: &mut RWMol, idx: u32) {
    let degree = mol.atom_degree(idx);
    let hs = match mol.atom(idx).atomic_num() {
        9 | 17 | 35 | 53 => 1u32.saturating_sub(degree),
        8 | 16 => 2u32.saturating_sub(degree),
        5 => 3u32.saturating_sub(degree),
        7 | 15 => match degree {
            0..=2 => 3 - degree,
            4 => 1,
            _ => 0,
        },
        6 => 4u32.saturating_sub(degree),
        _ => 0,
    };
    let atom = mol.atom_mut(idx);
    atom.set_no_implicit(true);
    atom.set_num_explicit_hs(hs);
}

/// `AnonymousGraph`: the element graph (`elem`) or the bare skeleton.
fn anonymous_graph(mol: &mut RWMol, elem: bool, cx: Cx) -> Result<String> {
    for idx in 0..mol.num_atoms() {
        let atom = mol.atom_mut(idx);
        atom.set_is_aromatic(false);
        atom.set_formal_charge(0);
        if elem {
            normalize_h_count(mol, idx);
        } else {
            atom.set_num_explicit_hs(0);
            atom.set_no_implicit(true);
            atom.set_atomic_num(0);
            atom.set_isotope(0);
        }
    }
    for idx in 0..mol.num_bonds() {
        let bond = mol.bond_mut(idx);
        bond.set_bond_type(BondType::Single);
        bond.set_is_aromatic(false);
    }
    clean_up(mol)?;
    let mut result = smiles(mol, false, cx)?;
    add_cx_extensions(mol, &mut result, cx, cx_fields::RADICALS)?;
    Ok(result)
}

/// `MesomerHash` (`netq`) and the redox-pair hash.
fn mesomer_hash(mol: &mut RWMol, netq: bool, cx: Cx) -> Result<String> {
    let mut charge = 0;
    for idx in 0..mol.num_atoms() {
        let atom = mol.atom_mut(idx);
        charge += atom.formal_charge();
        atom.set_is_aromatic(false);
        atom.set_formal_charge(0);
    }
    for idx in 0..mol.num_bonds() {
        mol.bond_mut(idx).set_bond_type(BondType::Single);
    }
    clean_up(mol)?;
    let mut result = smiles(mol, false, cx)?;
    if netq {
        result += &format!("_{charge}");
    }
    add_cx_extensions(mol, &mut result, cx, cx_fields::RADICALS)?;
    Ok(result)
}

/// `TautomerHash` (`proto`: the protomer hash).
fn tautomer_hash(mol: &mut RWMol, proto: bool, cx: Cx) -> Result<String> {
    let (mut hcount, mut charge) = (0i32, 0i32);
    for idx in 0..mol.num_atoms() {
        let hs = mol.total_num_hs(idx, false)?;
        let atom = mol.atom_mut(idx);
        charge += atom.formal_charge();
        atom.set_is_aromatic(false);
        atom.set_formal_charge(0);
        if atom.atomic_num() != 6 {
            hcount += hs as i32;
            atom.set_no_implicit(true);
            atom.set_num_explicit_hs(0);
        }
    }
    for idx in 0..mol.num_bonds() {
        let b = mol.bond(idx);
        let carbons = mol.atom(b.begin_atom_idx()).atomic_num() == 6
            && mol.atom(b.end_atom_idx()).atomic_num() == 6;
        if b.bond_type() != BondType::Single && (b.is_conjugated() || !carbons) {
            let b = mol.bond_mut(idx);
            b.set_is_aromatic(false);
            b.set_bond_type(BondType::Single);
            b.set_stereo(BondStereo::StereoNone)?;
        }
    }
    clean_up(mol)?;
    let mut result = smiles(mol, false, cx)?;
    result += &if proto {
        format!("_{}", hcount - charge)
    } else {
        format!("_{hcount}_{charge}")
    };
    add_cx_extensions(mol, &mut result, cx, cx_fields::RADICALS)?;
    Ok(result)
}

/// Removes `atoms` in one batch edit.
fn remove_atoms(mol: &mut RWMol, atoms: &[u32]) -> Result<()> {
    mol.begin_batch_edit()?;
    for &a in atoms {
        mol.remove_atom(a)?;
    }
    mol.commit_batch_edit()
}

/// `TraverseForRing` / `DepthFirstSearchForRing`: can a ring be reached from
/// `nbor` without passing through `root`?
fn reaches_ring(mol: &ROMol, root: u32, nbor: u32) -> bool {
    fn walk(mol: &ROMol, atom: u32, visit: &mut [bool]) -> bool {
        visit[atom as usize] = true;
        for n in mol.atom_neighbors(atom).collect::<Vec<_>>() {
            if !visit[n as usize] && (in_ring(mol, n) || walk(mol, n, visit)) {
                return true;
            }
        }
        false
    }
    let mut visit = vec![false; mol.num_atoms() as usize];
    visit[root as usize] = true;
    walk(mol, nbor, &mut visit)
}

/// `IsInScaffold`: ring atoms and atoms linking two parts that lead to rings.
fn is_in_scaffold(mol: &ROMol, idx: u32) -> bool {
    in_ring(mol, idx)
        || mol
            .atom_neighbors(idx)
            .filter(|&n| reaches_ring(mol, idx, n))
            .count()
            > 1
}

/// `ExtendedMurckoScaffold`: the scaffold plus dummies for its substituents.
fn extended_murcko(mol: &mut RWMol, cx: Cx) -> Result<String> {
    ensure_rings(mol)?;
    let in_scaffold: Vec<bool> = (0..mol.num_atoms())
        .map(|i| is_in_scaffold(mol, i))
        .collect();
    let mut deletion = Vec::new();
    for idx in 0..mol.num_atoms() {
        if in_scaffold[idx as usize] {
            continue;
        }
        if mol.atom_neighbors(idx).any(|n| in_scaffold[n as usize]) {
            let atom = mol.atom_mut(idx);
            atom.set_atomic_num(0);
            atom.set_formal_charge(0);
            atom.set_no_implicit(true);
            atom.set_num_explicit_hs(0);
            atom.set_isotope(0);
        } else {
            deletion.push(idx);
        }
    }
    remove_atoms(mol, &deletion)?;
    clean_up(mol)?;
    let mut result = smiles(mol, false, cx)?;
    add_cx_extensions(mol, &mut result, cx, cx_fields::RADICALS)?;
    Ok(result)
}

/// `MurckoScaffoldHash`: ring systems and their linkers.
fn murcko_scaffold(mol: &mut RWMol, cx: Cx) -> Result<String> {
    loop {
        let mut deletion = Vec::new();
        for idx in 0..mol.num_atoms() {
            let degree = mol.atom_degree(idx);
            if degree >= 2 {
                continue;
            }
            if degree == 1 {
                let bond = mol.atom_bonds(idx).next().expect("degree 1");
                let (order, nbr) = (bond_order(bond), bond.other_atom_idx(idx)?);
                let hs = mol.total_num_hs(nbr, false)?;
                let nbr = mol.atom_mut(nbr);
                nbr.set_num_explicit_hs(hs + order);
                nbr.set_no_implicit(true);
            }
            deletion.push(idx);
        }
        remove_atoms(mol, &deletion)?;
        if deletion.is_empty() {
            break;
        }
    }
    clean_up(mol)?;
    let mut result = smiles(mol, false, cx)?;
    add_cx_extensions(mol, &mut result, cx, cx_fields::RADICALS)?;
    Ok(result)
}

fn net_charge(mol: &ROMol) -> String {
    mol.atoms()
        .iter()
        .map(Atom::formal_charge)
        .sum::<i32>()
        .to_string()
}

/// `SmallWorldHash`: bond, ring and optionally linker counts.
fn small_world(mol: &ROMol, brl: bool) -> String {
    let (atoms, bonds) = (mol.num_atoms(), mol.num_bonds());
    // unsigned arithmetic upstream: disconnected molecules wrap around
    let rings = (bonds + 1).wrapping_sub(atoms);
    if brl {
        let linkers = (0..atoms).filter(|&i| mol.atom_degree(i) == 2).count();
        format!("B{bonds}R{rings}L{linkers}")
    } else {
        format!("B{bonds}R{rings}")
    }
}

fn degree_vector(mol: &ROMol) -> String {
    let mut v = [0u32; 4];
    for idx in 0..mol.num_atoms() {
        if let d @ 1..=4 = mol.atom_degree(idx) {
            v[4 - d as usize] += 1;
        }
    }
    format!("{},{},{},{}", v[0], v[1], v[2], v[3])
}

fn has_double_bond(mol: &ROMol, idx: u32) -> bool {
    mol.atom_bonds(idx).any(|b| bond_order(b) == 2)
}

/// `RegioisomerBond`: how to break a bond: `None` keeps it, otherwise bit 0
/// puts an asterisk on the begin atom and bit 1 on the end atom (else an H).
fn regioisomer_bond(mol: &ROMol, idx: u32) -> Option<u8> {
    let bond = mol.bond(idx);
    if bond_order(bond) != 1 || bond_in_ring(mol, idx) {
        return None;
    }
    let (beg, end) = (bond.begin_atom_idx(), bond.end_atom_idx());
    let (beg_elem, end_elem) = (mol.atom(beg).atomic_num(), mol.atom(end).atomic_num());
    if beg_elem == 0 || end_elem == 0 {
        return None;
    }
    match (in_ring(mol, beg), in_ring(mol, end)) {
        (true, true) => return Some(0),
        (true, false) => return Some(2),
        (false, true) => return Some(1),
        _ => {}
    }
    if beg_elem != 6 && end_elem == 6 && !has_double_bond(mol, end) {
        Some(1)
    } else if beg_elem == 6 && end_elem != 6 && !has_double_bond(mol, beg) {
        Some(2)
    } else {
        None
    }
}

fn clear_ez_stereo(mol: &mut RWMol, atom: u32) -> Result<()> {
    for b in mol.atom_bond_indices(atom).collect::<Vec<_>>() {
        if mol.bond(b).stereo() as u8 > BondStereo::StereoAny as u8 {
            mol.bond_mut(b).set_stereo(BondStereo::StereoAny)?;
        }
    }
    Ok(())
}

/// One end of a broken bond: an asterisk atom or an extra H.
fn cap_end(mol: &mut RWMol, atom: u32, asterisk: bool) -> Result<()> {
    if asterisk {
        let mut star = Atom::new(0);
        star.set_no_implicit(true);
        let star = mol.add_atom(star, true);
        mol.add_bond(atom, star, BondType::Single)?;
    } else {
        let hs = mol.total_num_hs(atom, false)?;
        let a = mol.atom_mut(atom);
        a.set_num_explicit_hs(hs + 1);
        a.set_no_implicit(true);
    }
    Ok(())
}

/// `RegioisomerHash`: the molecule broken into ring and chain fragments.
fn regioisomer(mol: &mut RWMol, cx: Cx) -> Result<String> {
    // decisions are made on a copy, so bonds can be removed from `mol`
    let mut copy = mol.clone();
    ensure_rings(&mut copy)?;
    for idx in (0..copy.num_bonds()).rev() {
        let Some(split) = regioisomer_bond(&copy, idx) else {
            continue;
        };
        let (beg, end) = (
            copy.bond(idx).begin_atom_idx(),
            copy.bond(idx).end_atom_idx(),
        );
        mol.remove_bond(beg, end)?;
        clear_ez_stereo(mol, beg)?;
        clear_ez_stereo(mol, end)?;
        cap_end(mol, beg, split & 1 != 0)?;
        cap_end(mol, end, split & 2 != 0)?;
    }
    assign_stereochemistry(mol, true, true, false)?;
    let mut result = smiles(mol, false, cx)?;
    add_cx_extensions(mol, &mut result, cx, 0)?;
    Ok(result)
}

/// `ArthorSubOrderHash`: a fixed-width hex summary of the composition.
fn arthor_sub_order(mol: &ROMol) -> Result<String> {
    let (atoms, bonds) = (mol.num_atoms(), mol.num_bonds());
    let parts = if atoms == 0 {
        0
    } else {
        get_mol_frags_mapping(mol).0
    };
    let (mut c, mut o, mut z, mut isotopes, mut charged, mut radicals) =
        (0u32, 0u32, 0u32, 0u32, 0u32, 0u32);
    for atom in mol.atoms() {
        let elem = atom.atomic_num();
        let charge = atom.formal_charge();
        let valence = atom.total_valence()?;
        match elem {
            6 => {
                c += 1;
                radicals += u32::from(charge == 0 && valence != 4);
            }
            7 | 15 => {
                o += 1;
                radicals += u32::from(charge == 0 && valence != 3 && valence != 5);
            }
            8 => {
                o += 1;
                radicals += u32::from(charge != 0 && valence != 2);
            }
            9 => {
                o += 1;
                radicals += u32::from(charge != 0 && valence != 1);
            }
            17 | 35 | 53 => {
                o += 1;
                radicals += u32::from(charge == 0 && ![1, 3, 5, 7].contains(&valence));
            }
            16 => {
                o += 1;
                radicals += u32::from(charge == 0 && ![2, 4, 6].contains(&valence));
            }
            _ => {}
        }
        z += elem;
        isotopes += u32::from(atom.isotope() != 0);
        charged += u32::from(charge != 0);
    }
    let (a, b) = (atoms.min(0xffff), bonds.min(0xffff));
    Ok(format!(
        "{:04x}{:04x}{:02x}{:04x}{:04x}{:06x}{:02x}{:02x}{:02x}",
        a,
        b,
        parts.min(0xff),
        c.min(0xffff),
        o.min(0xffff),
        z.min(0xff_ffff),
        radicals.min(0xff),
        charged.min(0xff),
        isotopes.min(0xff)
    ))
}

/// `MolHash::MolHash`: the hash `func` of `mol`, which is modified in place.
pub fn mol_hash(
    mol: &mut RWMol,
    func: HashFunction,
    use_cx_smiles: bool,
    cx_flags_to_skip: u32,
) -> Result<String> {
    let cx = Cx {
        enabled: use_cx_smiles,
        skip: cx_flags_to_skip,
    };
    sanitize_hydrogens(mol)?;
    match func {
        HashFunction::AnonymousGraph => anonymous_graph(mol, false, cx),
        HashFunction::ElementGraph => anonymous_graph(mol, true, cx),
        HashFunction::CanonicalSmiles => {
            let mut result = smiles(mol, false, cx)?;
            add_cx_extensions(mol, &mut result, cx, 0)?;
            Ok(result)
        }
        HashFunction::MurckoScaffold => murcko_scaffold(mol, cx),
        HashFunction::ExtendedMurcko => extended_murcko(mol, cx),
        HashFunction::Mesomer => mesomer_hash(mol, true, cx),
        HashFunction::RedoxPair => mesomer_hash(mol, false, cx),
        HashFunction::HetAtomTautomer => tautomer_hash(mol, false, cx),
        HashFunction::HetAtomTautomerv2 => tautomer_v2::tautomer_hash_v2(mol, false, cx),
        HashFunction::HetAtomProtomer => tautomer_hash(mol, true, cx),
        HashFunction::HetAtomProtomerv2 => tautomer_v2::tautomer_hash_v2(mol, true, cx),
        HashFunction::MolFormula => formula::molecular_formula(mol),
        HashFunction::AtomBondCounts => Ok(format!("{},{}", mol.num_atoms(), mol.num_bonds())),
        HashFunction::NetCharge => Ok(net_charge(mol)),
        HashFunction::SmallWorldIndexBR => Ok(small_world(mol, false)),
        HashFunction::SmallWorldIndexBRL => Ok(small_world(mol, true)),
        HashFunction::DegreeVector => Ok(degree_vector(mol)),
        HashFunction::ArthorSubstructureOrder => arthor_sub_order(mol),
        HashFunction::Regioisomer => regioisomer(mol, cx),
    }
}

/// `MolHash::Strip`: removes the features selected by [`strip_type`] flags.
pub fn strip(mol: &mut RWMol, strip: u32) -> Result<()> {
    if strip & strip_type::ATOM_STEREO != 0 {
        for idx in 0..mol.num_atoms() {
            mol.atom_mut(idx).set_chiral_tag(ChiralType::Unspecified);
        }
        if !mol.stereo_groups().is_empty() {
            mol.set_stereo_groups(Vec::new())?;
        }
    }
    if strip & strip_type::BOND_STEREO != 0 {
        for idx in 0..mol.num_bonds() {
            if mol.bond(idx).stereo() as u8 > BondStereo::StereoAny as u8 {
                mol.bond_mut(idx).set_stereo(BondStereo::StereoAny)?;
            }
        }
    }
    if strip & strip_type::ISOTOPE != 0 {
        for idx in 0..mol.num_atoms() {
            mol.atom_mut(idx).set_isotope(0);
        }
    }
    if strip & strip_type::ATOM_MAP != 0 {
        for idx in 0..mol.num_atoms() {
            mol.atom_mut(idx).set_atom_map_num(0, false)?;
        }
    }
    if strip & strip_type::HYDROGEN != 0 {
        remove_hs(mol, &RemoveHsParameters::default(), true)?;
    }
    Ok(())
}

/// `MolHash::SplitMolecule`: a copy of every fragment.
pub fn split_molecule(mol: &RWMol) -> Result<Vec<RWMol>> {
    Ok(get_mol_frags_as_mols(mol, true, true)?.0)
}
