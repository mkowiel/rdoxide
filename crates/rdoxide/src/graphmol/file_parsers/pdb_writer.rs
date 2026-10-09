//! PDB file writer.
//!
//! Ported from RDKit `Code/GraphMol/FileParsers/PDBWriter.cpp` (2026.09.1).
//!
//! Copyright (C) 2013 Greg Landrum and NextMove Software.
//! Licensed under the BSD 3-Clause license (see LICENSES/RDKit-BSD-3-Clause.txt).
//!
//! Flavors (bit flags):
//! * 1: write MODEL/ENDMDL lines around each record
//! * 2: don't write any CONECT records
//! * 4: write CONECT records in both directions
//! * 8: don't use multiple CONECTs to encode bond order
//! * 16: write MASTER record
//! * 32: write TER record

use std::collections::HashMap;

use crate::graphmol::atom::Atom;
use crate::graphmol::bond::BondType;
use crate::graphmol::conformer::Conformer;
use crate::graphmol::kekulize::kekulize;
use crate::graphmol::monomer_info::AtomMonomerType;
use crate::graphmol::romol::{ROMol, RWMol};
use crate::rdgeneral::RDProps;
use crate::rdgeneral::common_properties as cp;
use crate::{RDKitError, Result};

/// `MolOps::Kekulize` default of `maxBackTracks`
const MAX_BACKTRACKS: u32 = 100;

/// The first `n` characters of `s`, right-aligned in `n` columns (`setw(n) << s.substr(0, n)`).
fn right_aligned(s: &str, n: usize) -> String {
    let t: String = s.chars().take(n).collect();
    format!("{t:>n$}")
}

/// The first character of `s`, or a blank.
fn first_or_blank(s: &str) -> char {
    s.chars().next().unwrap_or(' ')
}

/// `GetDefaultAtomNumber`: the next two-character serial of an element
/// (`1`, `2`, ... `9`, `10`, ... `99`, `A0` ... `Z9`, `AA` ... `ZZ`).
fn default_atom_number(atom: &Atom, counts: &mut HashMap<u32, u32>) -> String {
    let mut ret = [b' '; 2];
    let n = counts.entry(atom.atomic_num()).or_insert(0);
    *n += 1;
    let n = *n;
    match n {
        1 => ret[0] = b'1',
        2..=9 => ret[0] = b'0' + n as u8,
        10..=99 => {
            ret[0] = b'0' + (n / 10) as u8;
            ret[1] = b'0' + (n % 10) as u8;
        }
        100..=359 => {
            ret[0] = b'A' + ((n - 100) / 10) as u8;
            ret[1] = b'0' + ((n - 100) % 10) as u8;
        }
        360..=1035 => {
            ret[0] = b'A' + ((n - 360) / 26) as u8;
            ret[1] = b'A' + ((n - 360) % 26) as u8;
        }
        _ => {}
    }
    String::from_utf8_lossy(&ret).into_owned()
}

/// `GetPDBAtomLine`
fn atom_line(
    atom: &Atom,
    idx: u32,
    conf: Option<&Conformer>,
    counts: &mut HashMap<u32, u32>,
) -> Result<String> {
    let symbol = atom.symbol();
    let mut chars = symbol.chars();
    let (at1, at2) = match (chars.next(), chars.next()) {
        (None, _) => (' ', 'X'),
        (Some(c), None) => (' ', c),
        (Some(c1), Some(c2)) => (c1, c2.to_ascii_uppercase()),
    };

    let info = atom
        .monomer_info()
        .filter(|i| i.monomer_type == AtomMonomerType::PdbResidue);
    let pdb = info.map(|i| i.pdb.clone().unwrap_or_default());
    let mut line = String::new();
    match (info, &pdb) {
        (Some(info), Some(pdb)) => {
            line.push_str(if pdb.is_hetero_atom {
                "HETATM"
            } else {
                "ATOM  "
            });
            line.push_str(&format!("{:>5} ", idx + 1));
            if info.name.is_empty() {
                line.push(at1);
                line.push(at2);
                line.push_str(&default_atom_number(atom, counts));
            } else {
                line.push_str(&right_aligned(&info.name, 4));
            }
            line.push(first_or_blank(&pdb.alt_loc));
            line.push_str(&right_aligned(&info.residue_name, 3));
            line.push(' ');
            line.push(first_or_blank(&info.chain_id));
            line.push_str(&format!("{:>4}", info.residue_number));
            line.push(first_or_blank(&pdb.insertion_code));
            line.push_str("   ");
        }
        _ => {
            line.push_str(&format!("HETATM{:>5} ", idx + 1));
            line.push(at1);
            line.push(at2);
            line.push_str(&default_atom_number(atom, counts));
            line.push_str(" UNL     1    ");
        }
    }

    match conf {
        Some(conf) => {
            let p = conf.atom_pos(idx)?;
            line.push_str(&format!("{:8.3}{:8.3}{:8.3}", p.x, p.y, p.z));
        }
        None => line.push_str("   0.000   0.000   0.000"),
    }

    match &pdb {
        Some(pdb) => {
            line.push_str(&format!("{:6.2}{:6.2}", pdb.occupancy, pdb.temp_factor));
            line.push_str("          ");
        }
        None => line.push_str("  1.00  0.00          "),
    }

    line.push(at1);
    line.push(at2);
    match atom.formal_charge() {
        c @ 1..=9 => {
            line.push((b'0' + c as u8) as char);
            line.push('+');
        }
        c @ -9..=-1 => {
            line.push((b'0' + (-c) as u8) as char);
            line.push('-');
        }
        _ => line.push_str("  "),
    }
    Ok(line)
}

/// `GetPDBBondLines`: the CONECT records of one atom.
fn bond_lines(
    mol: &ROMol,
    idx: u32,
    all: bool,
    both: bool,
    mult: bool,
    conect_count: &mut u32,
) -> String {
    let src = idx + 1;
    let mut targets: Vec<u32> = Vec::new();
    for bond in mol.atom_bonds(idx) {
        let Ok(other) = bond.other_atom_idx(idx) else {
            continue;
        };
        let dst = other + 1;
        if dst < src && !both {
            continue;
        }
        let bond_type = if mult {
            bond.bond_type()
        } else {
            BondType::Single
        };
        let repeats = match bond_type {
            BondType::Quadruple => 3,
            BondType::Triple => 2,
            BondType::Double => 1,
            _ => {
                if all {
                    targets.push(dst);
                }
                continue;
            }
        };
        // a double bond is written as two entries, and so on
        targets.extend(std::iter::repeat_n(dst, repeats + 1));
    }
    if targets.is_empty() {
        return String::new();
    }
    targets.sort_unstable();
    let mut out = String::new();
    for (i, dst) in targets.iter().enumerate() {
        if i % 4 == 0 {
            if i != 0 {
                out.push('\n');
            }
            out.push_str(&format!("CONECT{src:>5}"));
            *conect_count += 1;
        }
        out.push_str(&format!("{dst:>5}"));
    }
    out.push('\n');
    out
}

/// Running totals for the MASTER record.
#[derive(Default)]
struct Counts {
    atoms: u32,
    ter: u32,
    conect: u32,
}

/// `MolToPDBBody`
fn mol_to_pdb_body(
    mol: &ROMol,
    conf: Option<&Conformer>,
    flavor: u32,
    counts: &mut Counts,
) -> Result<String> {
    let mut res = String::new();
    let mut last = String::new();
    let mut elements: HashMap<u32, u32> = HashMap::new();
    for (idx, atom) in mol.atoms().iter().enumerate() {
        last = atom_line(atom, idx as u32, conf, &mut elements)?;
        res.push_str(&last);
        res.push('\n');
        counts.atoms += 1;
    }

    if counts.ter == 0 && counts.atoms > 0 && flavor & 32 != 0 {
        res.push_str(&format!("TER   {:>5}", counts.atoms + 1));
        if last.len() >= 27 {
            res.push_str("      ");
            res.push_str(&last[17..27]);
        }
        res.push('\n');
        counts.ter = 1;
    }

    let all = flavor & 2 == 0;
    let both = flavor & 4 != 0;
    let mult = flavor & 8 == 0;
    if all || mult {
        for idx in 0..mol.num_atoms() {
            res.push_str(&bond_lines(mol, idx, all, both, mult, &mut counts.conect));
        }
    }
    Ok(res)
}

/// `MolToPDBBlock`
pub fn mol_to_pdb_block(mol: &ROMol, conf_id: i32, flavor: u32) -> Result<String> {
    let mut rwmol = RWMol::from(mol.clone());
    kekulize(&mut rwmol, true, false, MAX_BACKTRACKS)?;

    let mut res = String::new();
    if let Some(name) = rwmol.get_prop(cp::_Name) {
        let name = name.to_rd_string();
        if !name.is_empty() {
            res.push_str(&format!("COMPND    {name}\n"));
        }
    }

    let mut counts = Counts::default();
    if conf_id < 0 && rwmol.num_conformers() > 1 {
        for id in 0..rwmol.num_conformers() {
            res.push_str(&format!("MODEL     {:>4}\n", id + 1));
            let conf = rwmol.conformer(id as i32)?;
            res.push_str(&mol_to_pdb_body(&rwmol, Some(conf), flavor, &mut counts)?);
            res.push_str("ENDMDL\n");
        }
    } else {
        let conf = if conf_id < 0 && rwmol.num_conformers() == 0 {
            None
        } else {
            Some(rwmol.conformer(conf_id)?)
        };
        res.push_str(&mol_to_pdb_body(&rwmol, conf, flavor, &mut counts)?);
    }

    if flavor & 16 != 0 {
        res.push_str(&format!(
            "MASTER        0    0    0    0    0    0    0    0{:>5}{:>5}{:>5}    0\n",
            counts.atoms, counts.ter, counts.conect
        ));
    }
    res.push_str("END\n");
    Ok(res)
}

/// `MolToPDBFile`
pub fn mol_to_pdb_file(mol: &ROMol, file_name: &str, conf_id: i32, flavor: u32) -> Result<()> {
    let block = pdb_writer_record(mol, conf_id, flavor, 1)?;
    std::fs::write(file_name, block)
        .map_err(|_| RDKitError::BadFile(format!("Bad output file {file_name}")))
}

/// What `PDBWriter::write` emits for the `count`-th molecule: the PDB block,
/// wrapped in MODEL/ENDMDL when flavor bit 1 is set.
pub fn pdb_writer_record(mol: &ROMol, conf_id: i32, flavor: u32, count: u32) -> Result<String> {
    let mut out = String::new();
    if flavor & 1 != 0 {
        out.push_str(&format!("MODEL     {count:>4}\n"));
    }
    out.push_str(&mol_to_pdb_block(mol, conf_id, flavor)?);
    if flavor & 1 != 0 {
        out.push_str("ENDMDL\n");
    }
    Ok(out)
}
