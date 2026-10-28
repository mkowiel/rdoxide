//! Port of the SMILES reading helpers of Avalon Toolkit's `smi2mol.c`: the bracket atom
//! reader (`SmilesToMDLAtom`), the valence preprocessing for denormalization
//! (`ComputeSmilesValenceAndHCount`) and the syntax pre-check (`QuickCheck`).

use super::denormal::cct_denormalize;
use super::graph::ring_list;
use super::layout::{DONT_FLIP_BOND, KEEP_POSITION, layout_molecule, recolor_molecule};
use super::perceive::{AROMATIC, DOUBLE, SINGLE, TRIPLE, setup_neighbourhood};
use super::reaccs::{Atom, Bond, CIS_TRANS_EITHER, Molecule, NONE};
use super::tables::PTABLE;
use super::utilities::{DOWN, UP, atom_symbol_match, compute_implicit_h};

/// flags of `smi2mol.h`
pub const DO_LAYOUT: i32 = 0x0001;
pub const EXPECT_SMARTS: i32 = 0x0002;
pub const DROP_TRIVIAL_HYDROGENS: i32 = 0x0008;
pub const TRUE_DB_STEREO: i32 = 0x0080;
const TRANS_MASK: i32 = 0x10;
const CIS_MASK: i32 = 0x20;
const AT_LEAST_1: i32 = 2;
const ANY_BOND: i32 = 8;
const MAXRINGS: usize = 70;
const MDL_MAXLINE: usize = 80;
const DOUBLET: i32 = 2;
const TRIPLET: i32 = 3;
const SINGLE_AROMATIC: i32 = 6;

/// the MDL attributes of a SMILES atom in square brackets
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SmilesAtom {
    pub symbol: String,
    pub charge: i32,
    pub radical: i32,
    pub isotope: i32,
    pub hydrogen: i32,
    pub stereo: i32,
    pub mapping: i32,
    pub aromatic: bool,
}

/// `SmilesToMDLAtom`: reads the atom at `smiles[pos..]`, which must be at the `[`. Returns the
/// position of the terminal `]` (or where the scan stopped) and the atom; `symbol` is empty
/// on errors. If the text is not at a `[`, `pos` is returned unchanged with an empty atom.
pub fn smiles_to_mdl_atom(smiles: &[u8], pos: usize) -> (usize, SmilesAtom) {
    let at = |i: usize| smiles.get(i).copied().unwrap_or(0);
    let mut a = SmilesAtom::default();
    let mut p = pos;
    if at(p) == b'[' {
        p += 1;
    } else {
        return (p, a);
    }
    while at(p).is_ascii_digit() {
        a.isotope = a
            .isotope
            .wrapping_mul(10)
            .wrapping_add(i32::from(at(p) - b'0'));
        p += 1;
    }
    let mut symbol: Vec<u8> = Vec::new();
    if at(p).is_ascii_uppercase() || at(p) == b'*' {
        symbol.push(at(p));
    } else if at(p).is_ascii_lowercase() {
        a.aromatic = true;
        symbol.push(at(p) - b'a' + b'A');
    } else {
        return (p, a);
    }
    p += 1;
    while at(p).is_ascii_lowercase() {
        symbol.push(at(p));
        p += 1;
    }
    while at(p) == b'@' {
        a.stereo += 1;
        p += 1;
    }
    if a.stereo == 1 && at(p).is_ascii_digit() {
        a.stereo = 0;
        while at(p).is_ascii_digit() {
            a.stereo = a
                .stereo
                .wrapping_mul(10)
                .wrapping_add(i32::from(at(p) - b'0'));
            p += 1;
        }
    }
    let mut sym = String::from_utf8_lossy(&symbol).into_owned();
    if a.isotope > 0 && symbol[0].is_ascii_uppercase() {
        let Some(e) = PTABLE.iter().find(|e| e.0 == sym) else {
            return (p, a);
        };
        if sym == "H" {
            if a.isotope == 2 {
                sym = "D".into();
            } else if a.isotope == 3 {
                sym = "T".into();
            }
            a.isotope = 0;
        } else {
            a.isotope -= (f64::from(e.2) + 0.5) as i32;
        }
    }
    if at(p) == b'H' {
        p += 1;
        if at(p).is_ascii_digit() {
            a.hydrogen = i32::from(at(p) - b'0');
            p += 1;
        } else {
            a.hydrogen = 1;
        }
    }
    a.hydrogen += 1; // hydrogen is always explicit
    while at(p) == b'+' || at(p) == b'-' {
        if a.charge >= 0 && at(p) == b'+' {
            a.charge += 1;
        } else if a.charge <= 0 && at(p) == b'-' {
            a.charge -= 1;
        } else {
            break;
        }
        p += 1;
    }
    if (b'1'..=b'9').contains(&at(p)) && (a.charge == -1 || a.charge == 1) {
        a.charge *= i32::from(at(p) - b'0');
        p += 1;
    }
    if at(p).is_ascii_digit() && a.charge != 0 {
        if a.charge < 0 {
            a.charge = a.charge * 10 - i32::from(at(p) - b'0');
        } else {
            a.charge = a.charge * 10 + i32::from(at(p) - b'0');
        }
        p += 1;
    }
    if at(p) == b':' {
        p += 1;
        while at(p).is_ascii_digit() {
            a.mapping = a
                .mapping
                .wrapping_mul(10)
                .wrapping_add(i32::from(at(p) - b'0'));
            p += 1;
        }
        a.mapping += 1; // MDL mappings start with 1
    }
    if at(p) != b']' {
        sym.clear();
    }
    if sym == "*" || sym == "R" {
        if a.isotope > 0 && a.isotope <= 99 {
            sym = "R#".into();
            a.charge = a.isotope;
            a.isotope = 0;
        } else if a.charge > 0 && a.charge <= 99 {
            sym = "R#".into();
        } else if a.charge == -1 {
            sym = "Frm".into();
            a.charge = 0;
        } else if a.charge == -2 {
            sym = "To".into();
            a.charge = 0;
        } else {
            sym = "R".into();
        }
    }
    a.symbol = sym;
    (p, a)
}

/// `ComputeSmilesValenceAndHCount`: puts the SMILES valence and hydrogen count into the
/// `dummy1` and `dummy2` fields of the atoms
pub fn compute_smiles_valence_and_hcount(mp: &mut Molecule) -> Result<(), &'static str> {
    let nbp = setup_neighbourhood(mp).ok_or("too many neighbours")?;
    for i in 0..mp.atoms.len() {
        let mdl_valence = PTABLE
            .iter()
            .find(|e| e.0 == mp.atoms[i].atom_symbol)
            .map_or(0, |e| e.1);
        let (mut nsingle, mut ndouble, mut ntriple, mut naromatic) = (0, 0, 0, 0);
        for &b in &nbp[i].bonds {
            match mp.bonds[b].bond_type {
                SINGLE => nsingle += 1,
                DOUBLE => ndouble += 1,
                TRIPLE => ntriple += 1,
                _ => naromatic += 1,
            }
        }
        let ap = &mut mp.atoms[i];
        let mut valence = (nsingle * 2 + ndouble * 4 + ntriple * 6 + naromatic * 3) / 2;
        valence += ap.charge.abs();
        if ap.query_h_count > 0 {
            valence += ap.query_h_count - 1;
        }
        if valence == mdl_valence {
            ap.dummy2 = if ap.query_h_count > 0 {
                ap.query_h_count - 1
            } else {
                0
            };
            ap.dummy1 = mdl_valence;
        } else if valence < mdl_valence {
            if ap.query_h_count > 0 {
                ap.dummy2 = ap.query_h_count - 1;
                ap.dummy1 = if naromatic == 2 && nsingle == 0 && ndouble == 0 && mdl_valence == 4 {
                    2
                } else {
                    valence
                };
            } else {
                ap.dummy2 = mdl_valence - valence;
                ap.dummy1 = mdl_valence;
            }
        } else if naromatic > 0 {
            if ap.query_h_count > 0 {
                ap.dummy2 = ap.query_h_count - 1;
                ap.dummy1 = if valence > mdl_valence {
                    mdl_valence + 2 * ((valence - mdl_valence) / 2)
                } else {
                    valence
                };
            } else if valence - mdl_valence == 3 {
                ap.dummy2 = 0;
                ap.dummy1 = mdl_valence + 2;
            } else if valence - mdl_valence > 1 {
                ap.dummy2 = 0;
                ap.dummy1 = valence;
            } else {
                ap.dummy2 = 0;
                ap.dummy1 = mdl_valence;
            }
        } else if ap.query_h_count > 0 {
            ap.dummy2 = ap.query_h_count - 1;
            ap.dummy1 = valence;
        } else {
            ap.dummy2 = (valence - mdl_valence) % 2;
            ap.dummy1 = valence + ap.dummy2;
        }
    }
    Ok(())
}

const SMILES_DIGIT: i32 = 16;
const SMILES_BOND: i32 = 32;
const SMILES_UC: i32 = 1;
const SMILES_CLOSE: i32 = 4;
const SMILES_AROMATIC: i32 = 64;
const SMILES_ARO_BOND: i32 = 128;
const SMILES_BEFORE_RING: i32 =
    SMILES_UC | SMILES_AROMATIC | SMILES_ARO_BOND | SMILES_BOND | SMILES_CLOSE | SMILES_DIGIT;

/// `smiles_class[]` of `InitSMILESClasses` (bytes above 127 have no class)
fn smiles_class(c: u8) -> i32 {
    let mut r = 0;
    if c.is_ascii_uppercase() {
        r |= SMILES_UC;
    }
    if matches!(c, b')' | b']' | b'}') {
        r |= SMILES_CLOSE;
    }
    if c.is_ascii_digit() {
        r |= SMILES_DIGIT;
    }
    if matches!(c, b'-' | b'=' | b'#' | b'\\' | b'/' | b'.') {
        r |= SMILES_BOND;
    }
    if c == b':' {
        r |= SMILES_ARO_BOND;
    }
    if matches!(c, b'c' | b'n' | b'o' | b's' | b'p' | b'*') {
        r |= SMILES_AROMATIC;
    }
    r
}

/// `QuickCheck`: a quick test if `smiles` can be a legal SMILES (matching brackets and
/// ring closures, legal characters before ring digits)
pub fn quick_check(smiles: &str) -> bool {
    let b: Vec<u8> = smiles.bytes().take_while(|&c| c != b'|').collect();
    if smiles.is_empty() {
        return false;
    }
    for (open, close) in [(b'[', b']'), (b'(', b')'), (b'{', b'}')] {
        let (mut bra, mut ket) = (0, 0);
        for &c in &b {
            if c == open {
                bra += 1;
            } else if c == close {
                ket += 1;
            }
            if bra < ket {
                return false;
            }
        }
        if bra != ket {
            return false;
        }
    }
    if b.first().is_some_and(|c| c.is_ascii_digit() || *c == b'%') {
        return false;
    }
    let at = |i: usize| b.get(i).copied().unwrap_or(0);
    let mut closure_count = [0i32; 100];
    let (mut in_atom, mut in_shortcut) = (0, 0);
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        if c == b'{' {
            in_shortcut += 1;
        } else if c == b'}' {
            in_shortcut -= 1;
        }
        if c == b'[' {
            in_atom += 1;
        } else if c == b']' {
            in_atom -= 1;
        }
        if in_atom > 0 || in_shortcut > 0 {
            i += 1;
            continue;
        }
        if smiles_class(c) & SMILES_DIGIT != 0 {
            if i == 0 || smiles_class(b[i - 1]) & SMILES_BEFORE_RING == 0 {
                return false;
            }
            closure_count[usize::from(c - b'0')] += 1;
        } else if c == b'%' {
            if i == 0 || smiles_class(b[i - 1]) & SMILES_BEFORE_RING == 0 {
                return false;
            }
            if smiles_class(at(i + 1)) & SMILES_DIGIT == 0
                || smiles_class(at(i + 2)) & SMILES_DIGIT == 0
            {
                return false;
            }
            closure_count[usize::from(10 * (at(i + 1) - b'0') + (at(i + 2) - b'0'))] += 1;
            i += 2;
        }
        i += 1;
    }
    closure_count.iter().all(|n| n % 2 == 0)
}

/// `FixBisaryl`: changes non-ring aromatic bonds to single/aromatic
fn fix_bisaryl(mp: &mut Molecule) {
    let graph: Vec<[u32; 2]> = mp
        .bonds
        .iter()
        .map(|b| [b.atoms[0] as u32, b.atoms[1] as u32])
        .collect();
    let rings = ring_list(&graph);
    let mut is_ring_bond = vec![false; mp.bonds.len()];
    for r in &rings {
        for (i, flag) in is_ring_bond.iter_mut().enumerate() {
            if r.bond_set.is_member(i) {
                *flag = true;
            }
        }
    }
    for (i, b) in mp.bonds.iter_mut().enumerate() {
        if b.bond_type == AROMATIC && !is_ring_bond[i] {
            b.bond_type = SINGLE_AROMATIC;
        }
    }
}

/// bond to the new atom `to` from `from`, followed by its reversed copy (SMILES stereo is
/// bond based, MOL file stereo is numbering based)
fn push_bond_pair(bonds: &mut Vec<Bond>, from: i32, to: i32, bond_type: i32, bond_dir: i32) {
    let value = if bond_dir == NONE {
        0.0
    } else if bond_dir == UP {
        1.0
    } else {
        -1.0
    };
    let b = Bond {
        atoms: [from, to],
        bond_type,
        value,
        stereo_symbol: NONE,
        ..Default::default()
    };
    let mut r = b.clone();
    r.value = -b.value;
    r.atoms = [to, from];
    bonds.push(b);
    bonds.push(r);
}

/// the `%lf` conversion of `sscanf` at the start of `s` (leading white space skipped)
fn scan_double(s: &[u8]) -> Option<f64> {
    let mut i = 0;
    while i < s.len() && s[i].is_ascii_whitespace() {
        i += 1;
    }
    let st = i;
    if i < s.len() && (s[i] == b'-' || s[i] == b'+') {
        i += 1;
    }
    let mut digits = 0;
    while i < s.len() && s[i].is_ascii_digit() {
        i += 1;
        digits += 1;
    }
    if i < s.len() && s[i] == b'.' {
        i += 1;
        while i < s.len() && s[i].is_ascii_digit() {
            i += 1;
            digits += 1;
        }
    }
    if digits == 0 {
        return None;
    }
    if i < s.len() && (s[i] == b'e' || s[i] == b'E') {
        let mut j = i + 1;
        if j < s.len() && (s[j] == b'-' || s[j] == b'+') {
            j += 1;
        }
        let ed = j;
        while j < s.len() && s[j].is_ascii_digit() {
            j += 1;
        }
        if j > ed {
            i = j;
        }
    }
    std::str::from_utf8(&s[st..i]).ok()?.parse().ok()
}

/// the coordinate string of `SMIToMOL`: `x,y,x,y,...` for the atoms with the given ordinal numbers
fn apply_coordinate_string(mp: &mut Molecule, coords: &str, atnos: &[i32]) {
    let cs = coords.as_bytes();
    let mut cp = 0usize;
    let mut atno = 0;
    let skip_field = |cp: &mut usize| {
        while *cp < cs.len() && cs[*cp] != b',' {
            *cp += 1;
        }
        if *cp < cs.len() {
            *cp += 1;
        }
    };
    loop {
        atno += 1;
        let mut novalue = false;
        let x = scan_double(&cs[cp..]).unwrap_or_else(|| {
            novalue = true;
            0.0
        });
        skip_field(&mut cp);
        let y = scan_double(&cs[cp..]).unwrap_or_else(|| {
            novalue = true;
            0.0
        });
        skip_field(&mut cp);
        if !novalue {
            if let Some(i) = atnos.iter().position(|&n| n == atno) {
                mp.atoms[i].x = x as f32;
                mp.atoms[i].y = y as f32;
                mp.atoms[i].color = KEEP_POSITION;
            }
        }
        if cp >= cs.len() {
            break;
        }
        if atno > mp.atoms.len() as i32 {
            break;
        }
    }
    for b in &mut mp.bonds {
        if mp.atoms[(b.atoms[0] - 1) as usize].color & KEEP_POSITION != 0
            && mp.atoms[(b.atoms[1] - 1) as usize].color & KEEP_POSITION != 0
        {
            b.bond_type |= DONT_FLIP_BOND;
        }
    }
}

/// `SMIToMOL`: converts a SMILES to a molecule; `Ok(None)` if the SMILES is illegal. The
/// atom colours hold the ordinal numbers of the explicit atoms. With `DO_LAYOUT` the molecule gets
/// 2D coordinates; a coordinate string after a `|` fixes the positions of the listed atoms.
pub fn smi_to_mol(smiles: &str, flags: i32) -> Result<Option<Molecule>, &'static str> {
    if !quick_check(smiles) {
        return Ok(None);
    }
    let (smiles, coordp) = match smiles.split_once('|') {
        Some((s, c)) => (s, Some(c)),
        None => (smiles, None),
    };
    let smiles = smiles.as_bytes();
    let at = |i: usize| smiles.get(i).copied().unwrap_or(0);
    let mut mp = Molecule {
        dimensionality: "2D".into(),
        scale1: 1,
        scale2: 1.0,
        version: "V2000".into(),
        ..Default::default()
    };
    let mut open_rings = [0i32; MAXRINGS];
    let mut open_aromatic = [false; MAXRINGS];
    let mut open_bonds = [0usize; MAXRINGS];
    let mut branch_stack: Vec<(i32, bool)> = Vec::new();
    let mut atno = 1;
    let (mut prev_atom, mut prev_aromatic) = (0i32, false);
    let (mut bond_type, mut bond_dir) = (NONE, NONE);
    let mut isomeric_smiles = flags & TRUE_DB_STEREO != 0;
    let mut chiral = false;
    let mut p = 0usize;
    macro_rules! natoms {
        () => {
            mp.atoms.len() as i32
        };
    }
    while p < smiles.len() {
        let c = smiles[p];
        match c {
            b'C' | b'B' | b'N' | b'O' | b'P' | b'S' | b'F' | b'I' | b'*' | b'{' => {
                let mut atext = String::new();
                let mut symbol;
                if c == b'{' {
                    p += 1;
                    let st = p;
                    let mut i = 0;
                    while i < MDL_MAXLINE && at(p) != b'}' && at(p) != 0 {
                        p += 1;
                        i += 1;
                    }
                    if at(p) != b'}' || i == 0 {
                        return Ok(None);
                    }
                    atext = String::from_utf8_lossy(&smiles[st..p]).into_owned();
                    symbol = "R".to_string();
                } else if c == b'C' && at(p + 1) == b'l' {
                    symbol = "Cl".into();
                    p += 1;
                } else if c == b'B' && at(p + 1) == b'r' {
                    symbol = "Br".into();
                    p += 1;
                } else {
                    symbol = (c as char).to_string();
                }
                if symbol == "*" {
                    symbol = "R".into();
                }
                mp.atoms.push(Atom {
                    atom_symbol: symbol,
                    atext,
                    color: atno,
                    ..Default::default()
                });
                atno += 1;
                if prev_atom > 0 {
                    let bt = if bond_type == NONE { SINGLE } else { bond_type };
                    push_bond_pair(&mut mp.bonds, prev_atom, natoms!(), bt, bond_dir);
                }
                bond_type = NONE;
                bond_dir = NONE;
                prev_atom = natoms!();
                prev_aromatic = false;
            }
            b'c' | b'b' | b'n' | b'o' | b'p' | b's' => {
                let symbol = ((c - b'a' + b'A') as char).to_string();
                mp.atoms.push(Atom {
                    atom_symbol: symbol,
                    color: atno,
                    ..Default::default()
                });
                atno += 1;
                if prev_atom > 0 {
                    let bt = if prev_aromatic && bond_type == NONE {
                        AROMATIC
                    } else if bond_type == NONE {
                        SINGLE
                    } else {
                        bond_type
                    };
                    push_bond_pair(&mut mp.bonds, prev_atom, natoms!(), bt, bond_dir);
                }
                bond_type = NONE;
                bond_dir = NONE;
                prev_atom = natoms!();
                prev_aromatic = true;
            }
            b'-' => {
                bond_type = SINGLE;
                bond_dir = NONE;
            }
            b'/' => {
                isomeric_smiles = true;
                bond_type = SINGLE;
                bond_dir = UP;
            }
            b'\\' => {
                isomeric_smiles = true;
                bond_type = SINGLE;
                bond_dir = DOWN;
            }
            b'=' => {
                bond_type = DOUBLE;
                bond_dir = NONE;
            }
            b'#' => {
                bond_type = TRIPLE;
                bond_dir = NONE;
            }
            b':' => {
                bond_type = AROMATIC;
                bond_dir = NONE;
            }
            b'%' | b'1'..=b'9' => {
                let iring: usize = if c == b'%' {
                    p += 1;
                    let mut v = 0usize;
                    let mut idigit = 0;
                    while idigit < 2 && at(p).is_ascii_digit() {
                        v = 10 * v + usize::from(at(p) - b'0');
                        p += 1;
                        idigit += 1;
                    }
                    p -= 1;
                    v
                } else {
                    usize::from(c - b'0')
                };
                if iring >= MAXRINGS {
                    return Ok(None);
                }
                if open_rings[iring] > 0 {
                    let bi = open_bonds[iring];
                    if open_aromatic[iring] && prev_aromatic && bond_type == NONE {
                        mp.bonds[bi].bond_type = AROMATIC;
                    } else if bond_type == NONE {
                        if mp.bonds[bi].bond_type <= SINGLE {
                            mp.bonds[bi].bond_type = SINGLE;
                        }
                    } else if mp.bonds[bi].bond_type <= bond_type {
                        mp.bonds[bi].bond_type = bond_type;
                    }
                    mp.bonds[bi].atoms = [open_rings[iring], prev_atom];
                    open_rings[iring] = 0;
                    let mut rev = mp.bonds[bi].clone();
                    rev.atoms = [mp.bonds[bi].atoms[1], mp.bonds[bi].atoms[0]];
                    if mp.bonds[bi].value == 0.0 {
                        if bond_dir == UP {
                            mp.bonds[bi].value = -1.0;
                        } else if bond_dir == DOWN {
                            mp.bonds[bi].value = 1.0;
                        }
                    }
                    rev.value = -mp.bonds[bi].value;
                    mp.bonds.push(rev);
                } else {
                    open_rings[iring] = prev_atom;
                    open_aromatic[iring] = prev_aromatic;
                    open_bonds[iring] = mp.bonds.len();
                    let value = if bond_dir == NONE {
                        0.0
                    } else if bond_dir == UP {
                        1.0
                    } else {
                        -1.0
                    };
                    mp.bonds.push(Bond {
                        bond_type,
                        stereo_symbol: NONE,
                        value,
                        ..Default::default()
                    });
                }
                bond_type = NONE;
                bond_dir = NONE;
            }
            b'[' => {
                let (np, a) = smiles_to_mdl_atom(smiles, p);
                p = np;
                if a.symbol.is_empty() {
                    return Ok(None);
                }
                let mut ap = Atom {
                    atom_symbol: a.symbol.clone(),
                    ..Default::default()
                };
                if a.symbol == "R#" {
                    // AddNumProperty
                    mp.prop_lines
                        .insert(0, format!("M  RGP  1 {:3} {:3}", natoms!() + 1, a.charge));
                    mp.n_props += 1;
                    ap.charge = 0;
                } else {
                    ap.charge = a.charge;
                }
                ap.radical = a.radical;
                ap.mass_difference = a.isotope;
                if a.isotope != 0 {
                    isomeric_smiles = true;
                }
                ap.query_h_count = a.hydrogen;
                ap.mapping = a.mapping;
                ap.stereo_parity = a.stereo;
                if a.stereo != NONE {
                    isomeric_smiles = true;
                    chiral = true;
                }
                ap.color = atno;
                atno += 1;
                mp.atoms.push(ap);
                if prev_atom > 0 {
                    let bt = if a.aromatic && prev_aromatic && bond_type == NONE {
                        AROMATIC
                    } else if bond_type == NONE {
                        SINGLE
                    } else {
                        bond_type
                    };
                    push_bond_pair(&mut mp.bonds, prev_atom, natoms!(), bt, bond_dir);
                }
                bond_type = NONE;
                bond_dir = NONE;
                prev_atom = natoms!();
                prev_aromatic = a.aromatic;
                if a.stereo != NONE && a.hydrogen == AT_LEAST_1 {
                    push_bond_pair(&mut mp.bonds, prev_atom, prev_atom + 1, SINGLE, NONE);
                    mp.atoms.push(Atom {
                        atom_symbol: "H".into(),
                        ..Default::default()
                    });
                    let last = mp.atoms.len() - 2;
                    mp.atoms[last].query_h_count = NONE;
                }
            }
            b'(' => branch_stack.push((prev_atom, prev_aromatic)),
            b')' => {
                let Some((pa, pr)) = branch_stack.pop() else {
                    return Ok(None);
                };
                prev_atom = pa;
                prev_aromatic = pr;
            }
            b'.' => {
                if at(p + 1).is_ascii_digit() || at(p + 1) == b'%' {
                    bond_type = ANY_BOND;
                    bond_dir = NONE;
                } else {
                    prev_atom = 0;
                    prev_aromatic = false;
                    bond_type = NONE;
                    bond_dir = NONE;
                }
            }
            _ => return Ok(None),
        }
        p += 1;
    }
    mp.chiral_flag = i32::from(chiral);
    if mp.atoms.is_empty() || open_rings.iter().any(|&r| r != 0) {
        return Ok(None);
    }
    let na = mp.atoms.len();
    let nb = mp.bonds.len();
    if mp.bonds.iter().any(|b| {
        b.atoms[0] < 1 || b.atoms[1] < 1 || b.atoms[0] as usize > na || b.atoms[1] as usize > na
    }) {
        return Err("bond to a missing atom");
    }

    // convert SMILES atom stereo parities to atom number parities
    for i in 0..na {
        if mp.atoms[i].stereo_parity == NONE {
            continue;
        }
        let mut ligands: Vec<i32> = Vec::new();
        let mut j = 0;
        while j < nb {
            if mp.bonds[j].atoms[0] == i as i32 + 1 {
                ligands.push(mp.bonds[j].atoms[1]);
                if ligands.len() > 4 {
                    break;
                }
            }
            j += 1;
        }
        if j < nb || ligands.len() < 4 {
            mp.atoms[i].stereo_parity = NONE;
            continue;
        }
        let mut parity = 1;
        for j in 1..ligands.len() {
            for k in (0..j).rev() {
                if ligands[k] > ligands[k + 1] {
                    ligands.swap(k, k + 1);
                    parity *= -1;
                }
            }
        }
        if parity < 0 {
            mp.atoms[i].stereo_parity = 3 - mp.atoms[i].stereo_parity;
        }
    }

    // convert SMILES bond stereo parities to atom numbering base
    for i in 0..nb {
        if mp.bonds[i].bond_type != DOUBLE {
            continue;
        }
        let (at1, at2) = (mp.bonds[i].atoms[0], mp.bonds[i].atoms[1]);
        if at1 > at2 {
            continue;
        }
        let (mut lig1, mut lig2) = (Vec::new(), Vec::new());
        let mut j = 0;
        while j < nb {
            let bh = &mp.bonds[j];
            if bh.atoms[0] == at1 && bh.atoms[1] != at2 {
                if bh.bond_type == AROMATIC {
                    mp.bonds[i].stereo_symbol = CIS_TRANS_EITHER;
                    break;
                }
                if bh.bond_type != SINGLE {
                    mp.bonds[i].stereo_symbol = NONE;
                    break;
                }
                if lig1.len() >= 2 {
                    break;
                }
                lig1.push(j);
            }
            let bh = &mp.bonds[j];
            if bh.atoms[0] == at2 && bh.atoms[1] != at1 {
                if bh.bond_type == AROMATIC {
                    mp.bonds[i].stereo_symbol = CIS_TRANS_EITHER;
                    break;
                }
                if bh.bond_type != SINGLE {
                    mp.bonds[i].stereo_symbol = NONE;
                    break;
                }
                if lig2.len() >= 2 {
                    break;
                }
                lig2.push(j);
            }
            j += 1;
        }
        if j < nb {
            continue;
        }
        if lig1.is_empty() || lig2.is_empty() {
            continue;
        }
        let v = |l: usize| mp.bonds[l].value;
        if (lig1.len() == 1 && v(lig1[0]) == 0.0)
            || (lig2.len() == 1 && v(lig2[0]) == 0.0)
            || (lig1.len() == 2 && v(lig1[0]) == 0.0 && v(lig1[1]) == 0.0)
            || (lig2.len() == 2 && v(lig2[0]) == 0.0 && v(lig2[1]) == 0.0)
        {
            mp.bonds[i].stereo_symbol = CIS_TRANS_EITHER;
            continue;
        }
        let asum = |l: usize| mp.bonds[l].atoms[0] + mp.bonds[l].atoms[1];
        if lig1.len() == 2 && asum(lig1[0]) > asum(lig1[1]) {
            lig1.swap(0, 1);
        }
        if lig2.len() == 2 && asum(lig2[0]) > asum(lig2[1]) {
            lig2.swap(0, 1);
        }
        let v = |l: usize| mp.bonds[l].value;
        if (lig1.len() == 2 && v(lig1[0]) == v(lig1[1]))
            || (lig2.len() == 2 && v(lig2[0]) == v(lig2[1]))
        {
            mp.bonds[i].stereo_symbol = CIS_TRANS_EITHER;
            continue;
        }
        let mut value1 = f64::from(v(lig1[0]));
        if lig1.len() == 2 {
            value1 -= f64::from(v(lig1[1]));
        }
        let mut value2 = f64::from(v(lig2[0]));
        if lig2.len() == 2 {
            value2 -= f64::from(v(lig2[1]));
        }
        if value1 * value2 < 0.0 {
            mp.bonds[i].stereo_symbol |= TRANS_MASK;
        } else {
            mp.bonds[i].stereo_symbol |= CIS_MASK;
        }
    }

    // remove duplicated bonds
    mp.bonds.retain(|b| b.atoms[0] < b.atoms[1]);
    if !isomeric_smiles {
        for b in &mut mp.bonds {
            if b.bond_type == DOUBLE {
                b.stereo_symbol = NONE;
            }
        }
    }
    let atnos: Vec<i32> = mp.atoms.iter().map(|a| a.color).collect();

    // CheckNeighbourhood
    let mut nbcounts = vec![0; na];
    for b in &mp.bonds {
        nbcounts[(b.atoms[0] - 1) as usize] += 1;
        nbcounts[(b.atoms[1] - 1) as usize] += 1;
    }
    if nbcounts.iter().any(|&n| n > 20) {
        return Ok(None);
    }

    if flags & EXPECT_SMARTS == 0 {
        compute_smiles_valence_and_hcount(&mut mp)?;
    } else {
        fix_bisaryl(&mut mp);
    }
    if flags & EXPECT_SMARTS == 0 {
        cct_denormalize(&mut mp);
    }
    for a in &mut mp.atoms {
        a.dummy1 = 0;
        a.dummy2 = 0;
    }

    // fix hydrogen count display
    let mut old_h = vec![0i32; na + 1];
    let mut new_h = vec![0i32; na + 1];
    for (i, a) in mp.atoms.iter_mut().enumerate() {
        old_h[i + 1] = a.query_h_count;
        a.query_h_count = 0;
    }
    compute_implicit_h(&mp, &mut new_h);
    for i in 0..na {
        if old_h[i + 1] != 0 && new_h[i + 1] + 1 != old_h[i + 1] {
            let sym = mp.atoms[i].atom_symbol.clone();
            if atom_symbol_match(&sym, "C,N,O") && new_h[i + 1] == old_h[i + 1] {
                mp.atoms[i].radical = DOUBLET;
            } else if atom_symbol_match(&sym, "C,N") && new_h[i + 1] == old_h[i + 1] + 1 {
                mp.atoms[i].radical = TRIPLET;
            } else {
                mp.atoms[i].query_h_count = old_h[i + 1];
            }
        }
    }
    if flags & DO_LAYOUT != 0 {
        recolor_molecule(&mut mp);
        if let Some(coords) = coordp {
            apply_coordinate_string(&mut mp, coords, &atnos);
        }
        mp = layout_molecule(&mp);
    }
    for (a, c) in mp.atoms.iter_mut().zip(atnos) {
        a.color = c;
    }
    mp.program_name = "NSMI2MOL".into();
    Ok(Some(mp))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::avalon::reaccsio::mol_str_to_mol;

    #[test]
    fn atoms_match_c() {
        let text = include_str!("../../../../../tests/data/avalon_smiles_atoms.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let (pos, a) = smiles_to_mdl_atom(input.as_bytes(), 0);
            let got = if !input.starts_with('[') {
                format!("{pos} [] -")
            } else {
                format!(
                    "{pos} [{}] {} {} {} {} {} {} {}",
                    a.symbol,
                    a.charge,
                    a.radical,
                    a.isotope,
                    a.hydrogen,
                    a.stereo,
                    a.mapping,
                    i32::from(a.aromatic)
                )
            };
            let got = format!("{got}|{}", i32::from(quick_check(input)));
            if got != want {
                bad.push(format!("{input:?}\nwant {want}\ngot  {got}"));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }

    #[test]
    fn valence_matches_c() {
        let text = include_str!("../../../../../tests/data/avalon_smiles_valence.jsonl");
        let mut bad = Vec::new();
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            let want = v["output"].as_str().unwrap();
            let mut mp = mol_str_to_mol(input).expect("molecule");
            compute_smiles_valence_and_hcount(&mut mp).unwrap();
            let mut got = String::new();
            for a in &mp.atoms {
                got += &format!("{} {}\n", a.dummy1, a.dummy2);
            }
            if got != want {
                bad.push(format!("{input}\nwant\n{want}\ngot\n{got}"));
            }
        }
        assert!(bad.is_empty(), "{} differ; first:\n{}", bad.len(), bad[0]);
    }

    #[test]
    fn smi_to_mol_matches_c() {
        use crate::external::avalon::reaccsio::mol_to_mol_str;
        let text = include_str!("../../../../../tests/data/avalon_smi2mol.jsonl");
        let mut bad = Vec::new();
        let mut n = 0;
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().unwrap();
            for (k, flags) in [0, TRUE_DB_STEREO, EXPECT_SMARTS].into_iter().enumerate() {
                let want = v["outputs"][k].as_str().unwrap();
                let got = match smi_to_mol(input, flags).expect("supported") {
                    None => "NULL\n".to_string(),
                    Some(m) => mol_to_mol_str(&m),
                };
                n += 1;
                if got != want {
                    bad.push(format!(
                        "{input:?} flags {flags}\nwant:\n{want}\ngot:\n{got}"
                    ));
                }
            }
        }
        assert!(
            bad.is_empty(),
            "{} of {n} differ; first:\n{}",
            bad.len(),
            bad[0]
        );
    }
}
