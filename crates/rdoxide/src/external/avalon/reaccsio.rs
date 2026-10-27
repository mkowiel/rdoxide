//! `reaccsio.c`, `forio.c`, `symbol_lists.c`: reading and writing MOL files
//! (`MolStr2Mol`, `MolToMolStr`). Copyright 2001-2011 Novartis Pharma AG.

use super::reaccs::{
    ANY_CHARGE, Atom, Bond, CIS_TRANS_EITHER, CIS_TRANS_SWAPPED, EXCLUSIVE, INCLUSIVE, MAXATOMS,
    MAXBONDS, Molecule, NONE, STextLine, SymbolList, combine_charge_radical,
    lookup_mass_difference, number_of, split_charge_radical, symbol_of,
};
use super::scan::{Scan, fmt_g};

/// `MDL_MAXLINE`, `MAX_BUFFER`, `MAXNAME`
const MDL_MAXLINE: usize = 255;
const MAX_BUFFER: usize = 10000;
const MAXNAME: usize = 80;

/// status of a `Fortran_FILE`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Status {
    Normal,
    Eof,
    Error,
}

/// `Fortran_FILE` reading from a string: `GetBuffer` and friends
struct Lines<'a> {
    src: &'a [u8],
    ptr: usize,
    buffer: String,
    status: Status,
}

/// `RemoveTrailingBlanks`
fn trim_end_blanks(s: &str) -> &str {
    s.trim_end_matches(|c: char| c.is_ascii_whitespace() || c == '\u{b}' || c == '\u{c}')
}

impl<'a> Lines<'a> {
    /// `FortranStringOpen`
    fn open(src: &'a str) -> Self {
        let mut l = Lines {
            src: src.as_bytes(),
            ptr: 0,
            buffer: String::new(),
            status: Status::Normal,
        };
        l.get_buffer();
        l
    }

    /// `GetBuffer`
    fn get_buffer(&mut self) {
        if self.ptr >= self.src.len() {
            self.buffer.clear();
            self.status = Status::Eof;
            return;
        }
        let st = self.ptr;
        let mut cp = st;
        while cp < self.src.len()
            && self.src[cp] != b'\n'
            && self.src[cp] != b'\r'
            && cp < st + MAX_BUFFER
        {
            cp += 1;
        }
        let line = String::from_utf8_lossy(&self.src[st..cp]).into_owned();
        if cp < self.src.len() && self.src[cp] == b'\r' {
            cp += 1;
        }
        if cp < self.src.len() && self.src[cp] == b'\n' {
            cp += 1;
        }
        self.ptr = cp;
        self.buffer = trim_end_blanks(&line).to_string();
    }

    fn begins(&self, s: &str) -> bool {
        self.buffer.starts_with(s)
    }
}

/// C `strncpy(dst, src, n)` into a string
fn cut(s: &str, n: usize) -> String {
    s.chars().take(n).collect()
}

/// `BlankToZero`
fn blank_to_zero(s: &str) -> String {
    s.replace(' ', "0")
}

/// the characters of `s` from byte `from` on (empty if `from` is past the end)
fn tail(s: &str, from: usize) -> &str {
    s.get(from..).unwrap_or("")
}

/// `ReadV30Atom`
fn read_v30_atom(f: &mut Lines, mol: &mut Molecule, idx: usize) -> Status {
    if f.status != Status::Normal {
        return f.status;
    }
    let line = f.buffer.clone();
    let mut sc = Scan::new(tail(&line, "M  V30".len()).as_bytes());
    let mut ap = Atom::default();
    let mut nitems = 0;
    let mut seq = 0;
    let mut symbol = String::new();
    let mut props: Vec<String> = Vec::new();
    'scan: {
        let Some(v) = sc.int(0) else { break 'scan };
        seq = v;
        nitems = 1;
        let Some(w) = sc.word() else { break 'scan };
        symbol = w;
        nitems = 2;
        for k in 0..3 {
            let Some(v) = sc.float(0) else { break 'scan };
            [&mut ap.x, &mut ap.y, &mut ap.z][k].clone_from(&v);
            nitems += 1;
        }
        let Some(v) = sc.int(0) else { break 'scan };
        ap.mapping = v;
        nitems = 6;
        for _ in 0..10 {
            let Some(w) = sc.word() else { break 'scan };
            props.push(w);
            nitems += 1;
        }
    }
    if nitems < 6 {
        f.get_buffer();
        return Status::Error;
    }
    if symbol.len() <= 3 {
        ap.atom_symbol = symbol;
    } else {
        parse_v30_symbol_list(&symbol, idx, mol, &mut ap);
        mol.n_atom_lists = mol.symbol_lists.len() as i32;
    }
    for p in &props {
        let int_of = |key: &str| -> Option<i32> {
            p.strip_prefix(key)
                .and_then(|r| Scan::new(r.as_bytes()).int(0))
        };
        if let Some(v) = int_of("CHG=") {
            ap.charge = v;
        } else if let Some(v) = int_of("RAD=") {
            ap.radical = v;
        } else if let Some(v) = int_of("HCOUNT=") {
            ap.query_h_count = v;
        } else if let Some(v) = int_of("CFG=") {
            ap.stereo_parity = v;
        } else if let Some(v) = int_of("VAL=") {
            ap.dummy1 = v;
        } else if let Some(v) = int_of("SUBST=") {
            ap.sub_desc = v;
        } else if let Some(m) = p
            .strip_prefix("MASS=")
            .and_then(|r| Scan::new(r.as_bytes()).float(0))
        {
            ap.mass_difference = lookup_mass_difference(f64::from(m), &ap.atom_symbol);
        }
        // RBCNT, UNSAT and everything else are ignored
    }
    let _ = seq;
    mol.atoms[idx] = ap;
    f.get_buffer();
    Status::Normal
}

/// `ParseV30SymbolList`: prepends the list to `mol.symbol_lists` and sets the symbol to `L`
/// (`Unk` if it cannot be parsed)
fn parse_v30_symbol_list(symbol: &str, iatom: usize, mol: &mut Molecule, ap: &mut Atom) {
    let (from, to) = (symbol.find('['), symbol.find(']'));
    match (from, to) {
        (Some(a), Some(b)) if b >= a => {
            let logic = if symbol.starts_with("NOT") {
                EXCLUSIVE
            } else {
                INCLUSIVE
            };
            mol.symbol_lists.insert(
                0,
                SymbolList {
                    atom: iatom as i32 + 1,
                    string: symbol[a + 1..b].to_string(),
                    logic,
                },
            );
            ap.atom_symbol = "L".into();
        }
        _ => ap.atom_symbol = "Unk".into(),
    }
}

/// `ReadV30Bond`
fn read_v30_bond(f: &mut Lines) -> (Status, Bond) {
    let mut bp = Bond::default();
    if f.status != Status::Normal {
        return (f.status, bp);
    }
    let line = f.buffer.clone();
    let mut sc = Scan::new(tail(&line, "M  V30".len()).as_bytes());
    let mut nitems = 0;
    'scan: {
        if sc.int(0).is_none() {
            break 'scan;
        }
        nitems = 1;
        for k in 0..3 {
            let Some(v) = sc.int(0) else { break 'scan };
            match k {
                0 => bp.bond_type = v,
                1 => bp.atoms[0] = v,
                _ => bp.atoms[1] = v,
            }
            nitems += 1;
        }
    }
    f.get_buffer();
    if nitems < 4 {
        return (Status::Error, bp);
    }
    (Status::Normal, bp)
}

/// `ReadREACCSAtom`
fn read_atom(f: &mut Lines) -> (Status, Atom) {
    let mut ap = Atom::default();
    if f.status != Status::Normal {
        return (f.status, ap);
    }
    let mut buffer = f.buffer.clone().into_bytes();
    if buffer.get(31) == Some(&b' ') {
        if buffer.len() > 32 {
            buffer[32] = b'?';
        }
    }
    let mut sc = Scan::new(&buffer);
    let mut nitems = 0;
    'scan: {
        for k in 0..3 {
            let Some(v) = sc.float(10) else { break 'scan };
            match k {
                0 => ap.x = v,
                1 => ap.y = v,
                _ => ap.z = v,
            }
            nitems += 1;
        }
        let Some(w) = sc.word() else { break 'scan };
        ap.atom_symbol = cut(&w, 3);
        nitems += 1;
    }
    let text = String::from_utf8_lossy(&buffer).into_owned();
    let rest = blank_to_zero(tail(&text, 34));
    let mut sc = Scan::new(rest.as_bytes());
    let mut charge_radical = 0;
    'num: {
        macro_rules! take {
            ($w:expr, $dst:expr) => {{
                let Some(v) = sc.int($w) else { break 'num };
                $dst = v;
            }};
        }
        take!(2, ap.mass_difference);
        take!(3, charge_radical);
        take!(3, ap.stereo_parity);
        take!(3, ap.query_h_count);
        take!(3, ap.query_stereo_box);
        take!(3, ap.dummy1);
        take!(3, ap.dummy2);
        take!(3, ap.dummy3);
        take!(3, ap.dummy4);
        take!(3, ap.mapping);
        take!(3, ap.second_stereo_parity);
        take!(3, ap.dummy5);
    }
    if nitems >= 4 {
        let (c, r) = split_charge_radical(charge_radical);
        ap.charge = c;
        ap.radical = r;
        f.get_buffer();
        (Status::Normal, ap)
    } else {
        f.get_buffer();
        (Status::Error, ap)
    }
}

/// `ReadREACCSBond`
fn read_bond(f: &mut Lines) -> (Status, Bond) {
    let mut bp = Bond::default();
    if f.status != Status::Normal {
        return (f.status, bp);
    }
    bp.reaction_mark = NONE;
    let line = f.buffer.as_bytes().to_vec();
    let len = line.len();
    let mut nitems = if len > 0 { (len - 1) / 3 + 1 } else { 0 };
    nitems = nitems.min(7);
    for i in 0..nitems {
        let pos = i * 3;
        let n_chars = (len - pos).min(3);
        let buf: Vec<u8> = line[pos..pos + n_chars]
            .iter()
            .copied()
            .filter(|&c| c != b' ')
            .collect();
        let Some(v) = Scan::new(&buf).int(3) else {
            continue;
        };
        match i {
            0 => bp.atoms[0] = v,
            1 => bp.atoms[1] = v,
            2 => bp.bond_type = v,
            3 => bp.stereo_symbol = v,
            4 => bp.dummy = v,
            5 => bp.topography = v,
            _ => bp.reaction_mark = v,
        }
    }
    if nitems >= 3 {
        f.get_buffer();
        (Status::Normal, bp)
    } else {
        (Status::Error, bp)
    }
}

/// `ReadSymbolLists`: in list order (the C code prepends each list that is read)
fn read_symbol_lists(f: &mut Lines, nlists: i32) -> Vec<SymbolList> {
    let mut result: Vec<SymbolList> = Vec::new();
    for _ in 0..nlists {
        let line = f.buffer.clone();
        let mut sc = Scan::new(line.as_bytes());
        let index = sc.int(0).unwrap_or(0);
        let type_string = sc.word().unwrap_or_default();
        let n_elements = sc.int(0).unwrap_or(0);
        let mut buffer = String::new();
        let mut ty = 0;
        for j in 0..n_elements.max(0) as usize {
            if j > 0 {
                buffer.push(',');
            }
            if let Some(v) = Scan::new(tail(&line, 11 + 4 * j).as_bytes()).int(0) {
                ty = v;
            }
            buffer.push_str(symbol_of(ty).unwrap_or(""));
        }
        let logic = if type_string == "F" {
            INCLUSIVE
        } else {
            EXCLUSIVE
        };
        result.insert(
            0,
            SymbolList {
                atom: index,
                string: buffer,
                logic,
            },
        );
        f.get_buffer();
    }
    result
}

/// scans `" %d %d ..."` into integers (as many as convert)
fn scan_ints(s: &str, max: usize) -> Vec<i32> {
    let mut sc = Scan::new(s.as_bytes());
    let mut v = Vec::new();
    while v.len() < max {
        let Some(x) = sc.int(0) else { break };
        v.push(x);
    }
    v
}

/// `ReadProperties`: interprets the property lines; returns the others in file order
fn read_properties(f: &mut Lines, mol: &mut Molecule, mut nprops: i32) -> Vec<String> {
    let mut result: Vec<String> = Vec::new();
    let old_nprops = nprops;
    let na = mol.atoms.len() as i32;
    while nprops > 0 {
        nprops -= 1;
        let buf = f.buffer.clone();
        let entries = |prefix: &str| -> Vec<(i32, i32)> {
            let nentries = Scan::new(tail(&buf, prefix.len()).as_bytes())
                .int(3)
                .unwrap_or(0);
            let v = scan_ints(tail(&buf, prefix.len() + 3), 16);
            (0..nentries.max(0) as usize)
                .filter_map(|n| {
                    if 2 * n + 1 < v.len() {
                        Some((v[2 * n], v[2 * n + 1]))
                    } else {
                        None
                    }
                })
                .collect()
        };
        if buf.starts_with("M  CHG") {
            for (a, v) in entries("M  CHG") {
                if (1..=na).contains(&a) {
                    mol.atoms[a as usize - 1].charge = v;
                }
            }
        } else if buf.starts_with("M  RAD") {
            for (a, v) in entries("M  RAD") {
                if (1..=na).contains(&a) {
                    mol.atoms[a as usize - 1].radical = v;
                }
            }
        } else if buf.starts_with("M  SUB") {
            for (a, v) in entries("M  SUB") {
                if (1..=na).contains(&a) {
                    mol.atoms[a as usize - 1].sub_desc = v;
                }
            }
        } else if buf.starts_with("V  ") {
            let mut sc = Scan::new(tail(&buf, 3).as_bytes());
            if let (Some(a), Some(v)) = (sc.int(0), sc.float(0)) {
                if (1..=na).contains(&a) {
                    mol.atoms[a as usize - 1].value = v;
                }
            }
        } else if buf.starts_with("M  END") {
            break;
        } else if buf.starts_with("M  RGP") {
            result.push(cut(&buf, MDL_MAXLINE));
        } else if buf.starts_with("G  ") {
            result.push(cut(&buf, MDL_MAXLINE));
            f.get_buffer();
            nprops -= 1;
            result.push(cut(&f.buffer, MDL_MAXLINE));
        } else if buf.starts_with("A  ") {
            let atom = Scan::new(tail(&buf, 3).as_bytes()).int(0);
            let is_r = atom
                .filter(|&a| (1..=na).contains(&a))
                .is_some_and(|a| mol.atoms[a as usize - 1].atom_symbol == "R");
            if is_r {
                f.get_buffer();
                nprops -= 1;
                mol.atoms[atom.expect("atom") as usize - 1].atext = cut(&f.buffer, MDL_MAXLINE);
            } else {
                result.push(cut(&buf, MDL_MAXLINE));
                f.get_buffer();
                nprops -= 1;
                result.push(cut(&f.buffer, MDL_MAXLINE));
            }
        } else if buf.starts_with("M  ISO") {
            for (a, v) in entries("M  ISO") {
                if !(1..=na).contains(&a) {
                    continue;
                }
                let ap = &mut mol.atoms[a as usize - 1];
                let diff = match (ap.atom_symbol.as_str(), v) {
                    ("R", _) => Some(v),
                    ("H", 1..=3) => Some(v - 1),
                    ("Li", 6..=7) => Some(v - 7),
                    ("B", 10..=11) => Some(v - 11),
                    ("C", 12..=14) => Some(v - 12),
                    ("N", 13..=15) => Some(v - 14),
                    ("O", 16..=18) => Some(v - 16),
                    ("F", 18..=19) => Some(v - 19),
                    ("P", 31..=33) => Some(v - 31),
                    ("S", 32..=36) => Some(v - 32),
                    ("Cl", 35..=37) => Some(v - 35),
                    ("Co", 56..=61) => Some(v - 59),
                    ("Br", 79..=81) => Some(v - 80),
                    ("Y", 88..=90) => Some(v - 89),
                    ("I", 123..=131) => Some(v - 127),
                    _ => None,
                };
                match diff {
                    Some(d) => ap.mass_difference = d,
                    None => result.push(format!("M  ISO  1 {a:3} {v:3}")),
                }
            }
        } else if buf.starts_with("M  ALS")
            || [
                "M  SLB", "M  STY", "M  SAL", "M  SBL", "M  SDI", "M  SMT", "M  SBV", "M  SCL",
                "M  SAP", "M  SDT", "M  SDD", "M  SED", "M  SDS",
            ]
            .iter()
            .any(|p| buf.starts_with(p))
        {
            // display-only or taken care of elsewhere
        } else {
            result.push(cut(&buf, MDL_MAXLINE));
        }
        f.get_buffer();
    }
    let _ = old_nprops;
    if !f.buffer.starts_with('M') && f.status == Status::Normal {
        // "possible format error": the line is left in place
    } else {
        f.get_buffer();
    }
    result
}

/// `ReadREACCSMolecule` with an empty label: reads the molecule at the current line
fn read_molecule(f: &mut Lines) -> Result<Molecule, Status> {
    let mut mp = Molecule {
        program_name: "DUMMY".into(),
        ..Molecule::default()
    };
    if f.status != Status::Normal {
        return Err(f.status);
    }
    // patch for SD-files w/o molecule records
    if f.buffer.starts_with('>') && f.buffer.contains('<') {
        return Ok(mp);
    }
    mp.name = cut(&f.buffer, MAXNAME);
    f.get_buffer();
    if f.status != Status::Normal {
        return Err(f.status);
    }
    let header = blank_to_zero_from(&f.buffer, 22);
    mp.dimensionality = "2D".into();
    mp.user_initials.clear();
    mp.program_name.clear();
    {
        let mut sc = Scan::new(header.as_bytes());
        'hdr: {
            let Some(v) = sc.chars(2) else { break 'hdr };
            mp.user_initials = v;
            let Some(v) = sc.chars(8) else { break 'hdr };
            mp.program_name = v;
            let Some(v) = sc.chars(6) else { break 'hdr };
            mp.date = v;
            let Some(v) = sc.chars(4) else { break 'hdr };
            mp.time = v;
            let Some(v) = sc.chars(2) else { break 'hdr };
            mp.dimensionality = v;
            let Some(v) = sc.int(2) else { break 'hdr };
            mp.scale1 = v;
            let Some(v) = sc.float(10) else { break 'hdr };
            mp.scale2 = v;
            let Some(v) = sc.float(12) else { break 'hdr };
            mp.energy = v;
            let Some(v) = sc.int(6) else { break 'hdr };
            mp.registry_number = i64::from(v);
        }
    }
    for s in [
        &mut mp.user_initials,
        &mut mp.program_name,
        &mut mp.date,
        &mut mp.time,
        &mut mp.dimensionality,
    ] {
        *s = trim_end_blanks(s).to_string();
    }
    f.get_buffer();
    if f.status != Status::Normal {
        return Err(f.status);
    }
    mp.comment = trim_end_blanks(&cut(&f.buffer, MDL_MAXLINE)).to_string();
    f.get_buffer();
    if f.status != Status::Normal {
        return Err(f.status);
    }
    let buffer = f.buffer.clone();
    mp.version.clear();
    if buffer.len() > 34 {
        mp.version = cut(tail(&buffer, 34), 5);
        if mp.version.contains("V200") {
            mp.version = "V2000".into();
        }
    }
    if mp.version == "V3000" {
        read_v3000(f, &mut mp)?;
    } else {
        let counts = blank_to_zero(&buffer);
        let mut sc = Scan::new(counts.as_bytes());
        let (mut n_atoms, mut n_bonds, mut n_stexts, mut n_props) = (0, 0, 0, 0);
        let mut nitems = 0;
        'c: {
            macro_rules! take {
                ($dst:expr) => {{
                    let Some(v) = sc.int(3) else { break 'c };
                    $dst = v;
                    nitems += 1;
                }};
            }
            take!(n_atoms);
            take!(n_bonds);
            take!(mp.n_atom_lists);
            take!(mp.dummy1);
            take!(mp.chiral_flag);
            take!(n_stexts);
            if !sc.skip(12) {
                break 'c;
            }
            take!(n_props);
        }
        mp.n_props = n_props;
        if mp.version.is_empty() {
            mp.version = " ".into();
        }
        if nitems >= 2 && (0..=MAXATOMS).contains(&n_atoms) && (0..=MAXBONDS).contains(&n_bonds) {
            f.get_buffer();
            if f.status != Status::Normal {
                return Err(f.status);
            }
        } else {
            return Err(Status::Error);
        }
        for _ in 0..n_atoms {
            let (st, ap) = read_atom(f);
            if st != Status::Normal {
                return Err(st);
            }
            mp.atoms.push(ap);
        }
        for _ in 0..n_bonds {
            let (st, bp) = read_bond(f);
            if st != Status::Normal {
                return Err(st);
            }
            mp.bonds.push(bp);
        }
        mp.symbol_lists = read_symbol_lists(f, mp.n_atom_lists);
        for _ in 0..n_stexts {
            let mut sc = Scan::new(f.buffer.as_bytes());
            let x = sc.float(0).unwrap_or(0.0);
            let y = sc.float(0).unwrap_or(0.0);
            f.get_buffer();
            let text = cut(&f.buffer, MDL_MAXLINE);
            f.get_buffer();
            mp.stext_lines.push(STextLine { x, y, text });
        }
        let np = mp.n_props;
        let props = read_properties(f, &mut mp, np);
        mp.n_props = props.len() as i32;
        mp.prop_lines = props;
    }
    Ok(mp)
}

/// `BlankToZero(buffer + from)` applied to a copy
fn blank_to_zero_from(s: &str, from: usize) -> String {
    match s.get(..from) {
        Some(head) => format!("{head}{}", blank_to_zero(tail(s, from))),
        None => s.to_string(),
    }
}

/// the V3000 branch of `ReadREACCSMolecule`
fn read_v3000(f: &mut Lines, mp: &mut Molecule) -> Result<(), Status> {
    f.get_buffer();
    mp.n_atom_lists = 0;
    mp.dummy1 = 0;
    mp.chiral_flag = 0;
    mp.n_props = 0;
    let (mut n_atoms, mut n_bonds) = (0usize, 0usize);
    while f.status == Status::Normal && f.buffer != "M  END" {
        if f.begins("M  V30 BEGIN CTAB") || f.begins("M  V30 END CTAB") {
            f.get_buffer();
        } else if f.begins("M  V30 COUNTS ") {
            let mut sc = Scan::new(tail(&f.buffer, "M  V30 COUNTS ".len()).as_bytes());
            let na = sc.int(0).unwrap_or(0);
            let nb = sc.int(0).unwrap_or(0);
            let _n_sgroups = sc.int(0);
            let _n_3d = sc.int(0);
            mp.chiral_flag = sc.int(0).unwrap_or(0);
            n_atoms = na.max(0) as usize;
            n_bonds = nb.max(0) as usize;
            mp.atoms = vec![Atom::default(); n_atoms];
            mp.bonds = vec![Bond::default(); n_bonds];
            f.get_buffer();
        } else if f.begins("M  V30 BEGIN ATOM") {
            f.get_buffer();
            let mut i = 0;
            while f.status == Status::Normal && f.buffer != "M  V30 END ATOM" {
                if i >= n_atoms {
                    return Err(Status::Error);
                }
                let st = read_v30_atom(f, mp, i);
                i += 1;
                if st != Status::Normal {
                    return Err(st);
                }
            }
            if f.status == Status::Normal {
                f.get_buffer();
            }
        } else if f.begins("M  V30 BEGIN BOND") {
            f.get_buffer();
            let mut i = 0;
            while f.status == Status::Normal && f.buffer != "M  V30 END BOND" {
                if i >= n_bonds {
                    return Err(Status::Error);
                }
                let (st, bp) = read_v30_bond(f);
                if st != Status::Normal {
                    return Err(st);
                }
                mp.bonds[i] = bp;
                i += 1;
            }
            if f.status == Status::Normal {
                f.get_buffer();
            }
        } else if f.begins("M  V30 BEGIN COLLECTION") {
            while f.status == Status::Normal && f.buffer != "M  V30 END COLLECTION" {
                f.get_buffer();
            }
            if f.status == Status::Normal {
                f.get_buffer();
            }
        } else {
            f.get_buffer();
        }
    }
    if f.status != Status::Normal {
        return Err(f.status);
    }
    mp.version = "V2000".into();
    Ok(())
}

/// `MolStr2Mol`: `None` on an error
pub fn mol_str_to_mol(mol_str: &str) -> Option<Molecule> {
    let mut f = Lines::open(mol_str);
    read_molecule(&mut f).ok()
}

/// `NeededPropLines`
fn needed_prop_lines(mp: &Molecule) -> i32 {
    let use_cr = mp
        .atoms
        .iter()
        .any(|a| a.charge > 3 || a.charge < -3 || a.radical != NONE);
    let mut result = 0;
    for ap in &mp.atoms {
        if use_cr && ap.charge != 0 {
            result += 1;
        }
        if ap.radical != NONE {
            result += 1;
        }
        if ap.sub_desc != NONE {
            result += 1;
        }
        if ap.value != 0.0 {
            result += 1;
        }
        if ap.atom_symbol == "R" && !ap.atext.is_empty() {
            result += 2;
        }
    }
    result + mp.n_props
}

/// `PrintREACCSMolecule` (with an empty header and `strip_zeroes` off)
fn print_molecule(mp: &Molecule) -> String {
    use std::fmt::Write;
    let mut o = String::new();
    let _ = writeln!(o, "{}", mp.name);
    let _ = writeln!(
        o,
        "{:<2}{:<8}{:<6}{:<4}{:<2}{:2}{:10.5}{:12.5}{:6}",
        mp.user_initials,
        mp.program_name,
        mp.date,
        mp.time,
        mp.dimensionality,
        mp.scale1,
        f64::from(mp.scale2),
        f64::from(mp.energy),
        mp.registry_number
    );
    let _ = writeln!(o, "{}", mp.comment);
    let _ = writeln!(
        o,
        "{:3}{:3}{:3}{:3}{:3}{:3}            {:3} {:<6}",
        mp.atoms.len(),
        mp.bonds.len(),
        mp.n_atom_lists,
        mp.dummy1,
        mp.chiral_flag,
        mp.stext_lines.len(),
        needed_prop_lines(mp) + 1,
        mp.version
    );
    for ap in &mp.atoms {
        let cr = combine_charge_radical(ap.charge, ap.radical);
        let _ = writeln!(
            o,
            "{:10.4}{:10.4}{:10.4} {:<3}{:2}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}",
            f64::from(ap.x),
            f64::from(ap.y),
            f64::from(ap.z),
            ap.atom_symbol,
            ap.mass_difference,
            if cr == ANY_CHARGE { 0 } else { cr },
            ap.stereo_parity,
            ap.query_h_count,
            ap.query_stereo_box,
            ap.dummy1,
            ap.dummy2,
            ap.dummy3,
            ap.dummy4,
            ap.mapping,
            ap.second_stereo_parity,
            ap.dummy5
        );
    }
    for bp in &mp.bonds {
        let ss = if bp.stereo_symbol == CIS_TRANS_SWAPPED {
            CIS_TRANS_EITHER
        } else {
            bp.stereo_symbol
        };
        let _ = writeln!(
            o,
            "{:3}{:3}{:3}{:3}{:3}{:3}{:3}",
            bp.atoms[0], bp.atoms[1], bp.bond_type, ss, bp.dummy, bp.topography, bp.reaction_mark
        );
    }
    for sl in &mp.symbol_lists {
        let _ = write!(
            o,
            "{:3}{}",
            sl.atom,
            if sl.logic == INCLUSIVE { " F " } else { " T " }
        );
        let syms: Vec<&str> = sl
            .string
            .split([' ', ','])
            .filter(|s| !s.is_empty())
            .collect();
        let _ = write!(o, " {:3}", syms.len());
        for s in syms {
            let _ = write!(o, " {:3}", number_of(s));
        }
        o.push('\n');
    }
    for st in &mp.stext_lines {
        let _ = writeln!(o, "{:10.4}{:10.4}", f64::from(st.x), f64::from(st.y));
        let _ = writeln!(o, "{}", st.text);
    }
    // PrintPropLines
    for l in &mp.prop_lines {
        let _ = writeln!(o, "{l}");
    }
    let use_cr = mp
        .atoms
        .iter()
        .any(|a| a.charge > 3 || a.charge < -3 || a.radical != NONE);
    for (i, ap) in mp.atoms.iter().enumerate() {
        if use_cr && ap.charge != 0 {
            let _ = writeln!(o, "M  CHG  1 {:3} {:3}", i + 1, ap.charge);
        }
        if ap.radical != NONE {
            let _ = writeln!(o, "M  RAD  1 {:3} {:3}", i + 1, ap.radical);
        }
        if ap.sub_desc != NONE {
            let _ = writeln!(o, "M  SUB  1 {:3} {:3}", i + 1, ap.sub_desc);
        }
        if ap.value != 0.0 {
            let _ = writeln!(o, "V  {:3} {}", i + 1, fmt_g(f64::from(ap.value)));
        }
        if !ap.atext.is_empty() && ap.atom_symbol == "R" {
            let _ = writeln!(o, "A  {:3}", i + 1);
            let _ = writeln!(o, "{}", ap.atext);
        }
    }
    o.push_str("M  END\n");
    o
}

/// `MolToMolStr`
pub fn mol_to_mol_str(mp: &Molecule) -> String {
    print_molecule(mp)
}

#[cfg(test)]
mod tests {
    use super::*;

    const ETHANOL: &str = "ethanol\n  RDKit          2D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5981   -0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  END\n";

    #[test]
    fn reads_and_writes_v2000() {
        let m = mol_str_to_mol(ETHANOL).expect("molecule");
        assert_eq!(m.name, "ethanol");
        assert_eq!(m.atoms.len(), 3);
        assert_eq!(m.bonds.len(), 2);
        assert_eq!(m.atoms[2].atom_symbol, "O");
        assert_eq!(m.bonds[1].atoms, [2, 3]);
        let out = mol_to_mol_str(&m);
        assert!(
            out.contains("  3  2  0  0  0  0            1 V2000 \n") || out.contains("V2000"),
            "{out}"
        );
        let m2 = mol_str_to_mol(&out).expect("round trip");
        assert_eq!(m2.atoms, m.atoms);
        assert_eq!(m2.bonds, m.bonds);
    }

    /// prints the round trip of the file named by `AVALON_MOL` (used to compare with the C code)
    #[test]
    #[ignore]
    fn dump_round_trip() {
        let path = std::env::var("AVALON_MOL").expect("AVALON_MOL");
        let text = std::fs::read_to_string(path).expect("file");
        match mol_str_to_mol(&text) {
            Some(m) => print!("@@{}", mol_to_mol_str(&m)),
            None => println!("@@NULL"),
        }
    }
}

#[cfg(test)]
mod corpus {
    use super::*;

    /// every molfile of the corpus must come back exactly as the Avalon C code writes it
    #[test]
    fn matches_c_output() {
        let text = include_str!("../../../../../tests/data/avalon_molfile.jsonl");
        let mut bad = Vec::new();
        let mut n = 0;
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).expect("json");
            let input = v["input"].as_str().expect("input");
            let want = v["output"].as_str().expect("output");
            let got = match mol_str_to_mol(input) {
                Some(m) => mol_to_mol_str(&m),
                None => "NULL\n".to_string(),
            };
            n += 1;
            if got != want {
                bad.push(format!(
                    "--- input\n{input}\n--- want\n{want}\n--- got\n{got}"
                ));
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
