//! `operator<<` for Atom/Bond and `ROMol::debugMol` (upstream `Atom.cpp`,
//! `Bond.cpp`, `ROMol.cpp`). BSD-3-Clause, RDKit contributors.

use std::fmt::Write;

use super::atom::{ChiralType, HybridizationType};
use super::bond::{BondDir, BondStereo, BondType};
use super::romol::ROMol;
use crate::rdgeneral::RDValue;
use crate::rdgeneral::common_properties as cp;

fn hybridization_str(h: HybridizationType) -> &'static str {
    match h {
        HybridizationType::Unspecified => "",
        HybridizationType::S => "S",
        HybridizationType::SP => "SP",
        HybridizationType::SP2 => "SP2",
        HybridizationType::SP3 => "SP3",
        HybridizationType::SP3D => "SP3D",
        HybridizationType::SP2D => "SP2D",
        HybridizationType::SP3D2 => "SP3D2",
        HybridizationType::Other => "OTHER",
    }
}

fn chirality_str(c: ChiralType) -> &'static str {
    match c {
        ChiralType::Unspecified => "Unspecified",
        ChiralType::TetrahedralCW => "CW",
        ChiralType::TetrahedralCCW => "CCW",
        ChiralType::Other => "Other",
        ChiralType::Tetrahedral => "Td",
        ChiralType::Allene => "Allene",
        ChiralType::SquarePlanar => "SqP",
        ChiralType::TrigonalBipyramidal => "Tbp",
        ChiralType::Octahedral => "Oh",
    }
}

fn bond_type_str(t: BondType) -> &'static str {
    use BondType::*;
    match t {
        Unspecified => "?",
        Single => "1",
        Double => "2",
        Triple => "3",
        Quadruple => "4",
        Quintuple => "5",
        Hextuple => "6",
        OneAndAHalf => "1.5",
        TwoAndAHalf => "2.5",
        ThreeAndAHalf => "3.5",
        FourAndAHalf => "4.5",
        FiveAndAHalf => "5.5",
        Aromatic => "a",
        Ionic => "I",
        Hydrogen => "H",
        ThreeCenter => "3C",
        DativeOne => "D1",
        Dative => "D",
        DativeL => "DL",
        DativeR => "DR",
        Other => "Other",
        Zero => "0",
    }
}

fn bond_dir_str(d: BondDir) -> &'static str {
    match d {
        BondDir::None => "NONE",
        BondDir::BeginWedge => "wedge",
        BondDir::BeginDash => "dash",
        BondDir::EndDownRight => "\\",
        BondDir::EndUpRight => "/",
        BondDir::EitherDouble => "x",
        BondDir::Unknown => "?",
    }
}

fn bond_stereo_str(s: BondStereo) -> &'static str {
    match s {
        BondStereo::StereoNone => "NONE",
        BondStereo::StereoAny => "ANY",
        BondStereo::StereoZ => "Z",
        BondStereo::StereoE => "E",
        BondStereo::StereoCis => "CIS",
        BondStereo::StereoTrans => "TRANS",
        BondStereo::StereoAtropCW => "CW",
        BondStereo::StereoAtropCCW => "CCW",
    }
}

/// `operator<<(std::ostream&, const Atom&)` for an atom of `mol`.
pub fn atom_to_string(mol: &ROMol, idx: u32) -> String {
    let at = mol.atom(idx);
    let mut s = format!("{} {} {}", at.idx(), at.atomic_num(), at.symbol());
    let _ = write!(s, " chg: {}", at.formal_charge());
    let _ = write!(s, "  deg: {}", mol.atom_degree(idx));
    let ev = at.cached_explicit_valence();
    let _ = write!(
        s,
        " exp: {}",
        if ev >= 0 {
            ev.to_string()
        } else {
            "N/A".into()
        }
    );
    s.push_str(" imp: ");
    if at.no_implicit() {
        s.push('0');
    } else {
        let iv = at.cached_implicit_valence();
        s.push_str(&if iv >= 0 {
            iv.to_string()
        } else {
            "N/A".into()
        });
    }
    let _ = write!(s, " hyb: {}", hybridization_str(at.hybridization()));
    if at.is_aromatic() {
        s.push_str(" arom?: 1");
    }
    if at.chiral_tag() != ChiralType::Unspecified {
        let _ = write!(s, " chi: {}", chirality_str(at.chiral_tag()));
        if let Some(perm) = crate::rdgeneral::RDProps::get_prop(at, cp::_chiralPermutation)
            .and_then(RDValue::as_int)
        {
            let _ = write!(s, "({perm})");
        }
        let nbrs: Vec<String> = mol.atom_neighbors(idx).map(|n| n.to_string()).collect();
        let _ = write!(s, " nbrs:[{}]", nbrs.join(" "));
    }
    if at.num_radical_electrons() != 0 {
        let _ = write!(s, " rad: {}", at.num_radical_electrons());
    }
    if at.isotope() != 0 {
        let _ = write!(s, " iso: {}", at.isotope());
    }
    if at.atom_map_num() != 0 {
        let _ = write!(s, " mapno: {}", at.atom_map_num());
    }
    s
}

/// `operator<<(std::ostream&, const Bond&)`
pub fn bond_to_string(mol: &ROMol, idx: u32) -> String {
    let b = mol.bond(idx);
    let mut s = format!(
        "{} {}->{} order: {}",
        b.idx(),
        b.begin_atom_idx(),
        b.end_atom_idx(),
        bond_type_str(b.bond_type())
    );
    if b.bond_dir() != BondDir::None {
        let _ = write!(s, " dir: {}", bond_dir_str(b.bond_dir()));
    }
    if b.stereo() != BondStereo::StereoNone {
        let _ = write!(s, " stereo: {}", bond_stereo_str(b.stereo()));
        if b.stereo_atoms().len() == 2 {
            let _ = write!(s, " ats: ({} {})", b.stereo_atoms()[0], b.stereo_atoms()[1]);
        }
    }
    if b.is_conjugated() {
        s.push_str(" conj?: 1");
    }
    if b.is_aromatic() {
        s.push_str(" aromatic?: 1");
    }
    s
}

/// `ROMol::debugMol`
pub fn debug_mol(mol: &ROMol) -> String {
    let mut s = String::from("Atoms:\n");
    for i in 0..mol.num_atoms() {
        let _ = writeln!(s, "\t{}", atom_to_string(mol, i));
    }
    s.push_str("Bonds:\n");
    for i in 0..mol.num_bonds() {
        let _ = writeln!(s, "\t{}", bond_to_string(mol, i));
    }
    let groups = mol.stereo_groups();
    if !groups.is_empty() {
        s.push_str("Stereo Groups:\n");
        for (i, g) in groups.iter().enumerate() {
            let _ = writeln!(s, "\t{i} {}", stereo_group_to_string(g));
        }
    }
    s
}

/// `operator<<(std::ostream&, const StereoGroup&)`
pub fn stereo_group_to_string(g: &super::StereoGroup) -> String {
    let kind = match g.group_type() {
        super::StereoGroupType::Absolute => "ABS",
        super::StereoGroupType::Or => "OR ",
        super::StereoGroupType::And => "AND",
    };
    let mut s = format!("{kind} rId: {} wId: {}", g.read_id(), g.write_id());
    if !g.atoms().is_empty() {
        s.push_str(" atoms: { ");
        for a in g.atoms() {
            let _ = write!(s, "{a} ");
        }
        s.push('}');
    }
    if !g.bonds().is_empty() {
        s.push_str(" bonds: { ");
        for b in g.bonds() {
            let _ = write!(s, "{b} ");
        }
        s.push('}');
    }
    s
}
